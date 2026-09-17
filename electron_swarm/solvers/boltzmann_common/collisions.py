"""Neutral-gas collision data and finite-volume collision assembly."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.constants import (
    AMU_KG,
    BOLTZMANN_J_K,
    ELECTRON_MASS_KG,
    EV_TO_J,
)
from electron_swarm.core.cross_sections import (
    CrossSectionProcess,
    CrossSectionSet,
    ProcessType,
    gas_mass_amu,
    mixture_fraction,
)
from electron_swarm.core.scattering import (
    resolve_species_scattering,
    transport_momentum_processes,
)
from electron_swarm.core.solver_configs import (
    InelasticCollisionConfig,
    MomentumCollisionConfig,
)
from electron_swarm.solvers.boltzmann_common.grid import electron_speed_m_s


@dataclass(slots=True)
class EffectiveCollisionData:
    nu_m: np.ndarray
    nu_m_over_N: np.ndarray
    sigma_m: np.ndarray
    sigma_total_like: np.ndarray
    inelastic_loss_frequency_s_inv: np.ndarray
    elastic_A_eV_s: np.ndarray
    elastic_D_eV2_s: np.ndarray
    processes: list[CrossSectionProcess]
    effective_momentum_model: str
    sigma_total_like_model: str
    suppressed_momentum_processes: tuple[str, ...]


def build_effective_collision_data(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    energy: np.ndarray,
    gas_number_density_m3: float,
    solver_config: MomentumCollisionConfig,
) -> EffectiveCollisionData:
    cfg = solver_config
    speed = electron_speed_m_s(energy)
    sigma_m = np.zeros_like(energy)
    sigma_total_like = np.zeros_like(energy)
    nu_m = np.zeros_like(energy)
    nu_m_over_N = np.zeros_like(energy)
    inelastic_loss = np.zeros_like(energy)
    elastic_A = np.zeros_like(energy)
    elastic_D = np.zeros_like(energy)
    kT_eV = config.conditions.gas_temperature_K * BOLTZMANN_J_K / EV_TO_J
    if not cross_sections.by_type(
        ProcessType.MOMENTUM, ProcessType.EFFECTIVE, ProcessType.ELASTIC
    ):
        raise ValueError(
            "At least one momentum/effective/elastic cross section is required"
        )
    transport_inelastic = {
        ProcessType.EXCITATION,
        ProcessType.IONIZATION,
        ProcessType.ATTACHMENT,
        ProcessType.SUPERELASTIC,
    }
    momentum_ids: set[int] = set()
    total_like_ids: set[int] = set()
    recoil_ids: set[int] = set()
    effective_species: set[str] = set()
    selection_models: list[str] = []
    suppressed: list[str] = []
    for species in cross_sections.species:
        selected, includes_inelastic, model = transport_momentum_processes(
            cross_sections,
            species,
        )
        if not selected:
            continue
        selection_models.append(f"{species}:{model}")
        momentum_ids.update(id(process) for process in selected)
        if includes_inelastic:
            effective_species.add(species)
        try:
            particle = resolve_species_scattering(
                cross_sections,
                species,
                angular_model="isotropic",
            )
        except NotImplementedError:
            # Effective-only data can drive momentum relaxation without
            # identifying a separate physical recoil event.
            total_like_ids.update(id(process) for process in selected)
        else:
            total_like_ids.update(id(process) for process in particle.collision_total)
            recoil_ids.update(id(process) for process in particle.collision_total)

    for proc in cross_sections.processes:
        frac = mixture_fraction(config.conditions, proc.species)
        if frac <= 0.0:
            continue
        sigma = np.asarray(proc.sigma(energy), dtype=float)
        if sigma.shape != energy.shape or not np.all(np.isfinite(sigma)):
            raise ValueError(
                f"Cross section {proc.species}:{proc.process} is non-finite"
            )
        if np.any(sigma < 0.0):
            raise ValueError(
                f"Cross section {proc.species}:{proc.process} is negative"
            )
        include_transport = id(proc) in momentum_ids
        if proc.process_type in transport_inelastic:
            include_transport = proc.species not in effective_species
            inelastic_loss += gas_number_density_m3 * frac * sigma * speed
        if include_transport:
            nuN = frac * sigma * speed
            nu = gas_number_density_m3 * nuN
            sigma_m += frac * sigma
            nu_m_over_N += nuN
            nu_m += nu
        if id(proc) in total_like_ids:
            sigma_total_like += frac * sigma
        if (
            proc.process_type
            in {ProcessType.MOMENTUM, ProcessType.EFFECTIVE, ProcessType.ELASTIC}
            and id(proc) not in momentum_ids
            and id(proc) not in total_like_ids
        ):
            suppressed.append(
                f"{proc.species}:{proc.scattering_role.value if proc.scattering_role else 'unknown'}:{proc.process}"
            )
        if id(proc) not in recoil_ids:
            continue
        mass_amu = proc.mass_amu or gas_mass_amu(config.conditions, proc.species)
        mass_kg = mass_amu * AMU_KG
        nu = gas_number_density_m3 * frac * sigma * speed
        ratio = ELECTRON_MASS_KG / mass_kg
        eps = np.maximum(energy, 1.0e-12)
        elastic_D += 2.0 * ratio * nu * eps * kT_eV
        elastic_A += 2.0 * ratio * nu * (0.5 * kT_eV - eps)
    nu_m = np.maximum(nu_m, 1.0e-60)
    nu_m_over_N = np.maximum(nu_m_over_N, 1.0e-80)
    sigma_m = np.maximum(sigma_m, cfg.min_momentum_cross_section_m2)
    sigma_total_like = np.maximum(
        sigma_total_like, cfg.min_momentum_cross_section_m2
    )
    model = ";".join(selection_models) if selection_models else "missing"
    return EffectiveCollisionData(
        nu_m=nu_m,
        nu_m_over_N=nu_m_over_N,
        sigma_m=sigma_m,
        sigma_total_like=sigma_total_like,
        inelastic_loss_frequency_s_inv=inelastic_loss,
        elastic_A_eV_s=elastic_A,
        elastic_D_eV2_s=np.maximum(elastic_D, 0.0),
        processes=cross_sections.processes,
        effective_momentum_model=model,
        sigma_total_like_model=model,
        suppressed_momentum_processes=tuple(suppressed),
    )


def _deposit_energy(
    mat: sparse.lil_matrix,
    energy: np.ndarray,
    widths: np.ndarray,
    source_i: int,
    target_e: float,
    rate: float,
) -> None:
    n = len(energy)
    if target_e <= energy[0]:
        mat[0, source_i] += rate * widths[source_i] / widths[0]
    elif target_e >= energy[-1]:
        mat[-1, source_i] += rate * widths[source_i] / widths[-1]
    else:
        j = int(np.searchsorted(energy, target_e) - 1)
        j = max(0, min(j, n - 2))
        denom = energy[j + 1] - energy[j]
        wR = (target_e - energy[j]) / denom
        wL = 1.0 - wR
        mat[j, source_i] += rate * widths[source_i] * wL / widths[j]
        mat[j + 1, source_i] += rate * widths[source_i] * wR / widths[j + 1]


def _add_energy_shift_transition(
    mat: sparse.lil_matrix,
    energy: np.ndarray,
    widths: np.ndarray,
    nu: np.ndarray,
    threshold: float,
    *,
    multiplicity: float,
) -> None:
    for i, rate in enumerate(np.asarray(nu, dtype=float)):
        if rate <= 0.0:
            continue
        mat[i, i] += -rate
        target_e = energy[i] - threshold
        _deposit_energy(mat, energy, widths, i, target_e, multiplicity * rate)


def _add_equal_sharing_ionization(
    mat: sparse.lil_matrix,
    energy: np.ndarray,
    widths: np.ndarray,
    nu: np.ndarray,
    threshold: float,
) -> None:
    for i, rate in enumerate(np.asarray(nu, dtype=float)):
        if rate <= 0.0:
            continue
        mat[i, i] += -rate
        target_e = max(0.0, 0.5 * (energy[i] - threshold))
        _deposit_energy(mat, energy, widths, i, target_e, 2.0 * rate)


def _add_primary_secondary_ionization(
    mat: sparse.lil_matrix,
    energy: np.ndarray,
    widths: np.ndarray,
    nu: np.ndarray,
    threshold: float,
    secondary_electron_energy_eV: float,
) -> None:
    secondary_eV = max(float(secondary_electron_energy_eV), 0.0)
    for i, rate in enumerate(np.asarray(nu, dtype=float)):
        if rate <= 0.0:
            continue
        mat[i, i] += -rate
        primary_e = max(0.0, energy[i] - threshold - secondary_eV)
        _deposit_energy(mat, energy, widths, i, primary_e, rate)
        _deposit_energy(mat, energy, widths, i, secondary_eV, rate)


def assemble_collision_operator(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    energy: np.ndarray,
    widths: np.ndarray,
    gas_number_density_m3: float,
    solver_config: InelasticCollisionConfig,
) -> sparse.csr_matrix:
    n = len(energy)
    mat = sparse.lil_matrix((n, n), dtype=float)
    speed = electron_speed_m_s(energy)
    cfg = solver_config
    for proc in cross_sections.processes:
        frac = mixture_fraction(config.conditions, proc.species)
        if frac <= 0.0:
            continue
        if proc.process_type not in {
            ProcessType.EXCITATION,
            ProcessType.SUPERELASTIC,
            ProcessType.IONIZATION,
            ProcessType.ATTACHMENT,
        }:
            continue
        nu = gas_number_density_m3 * frac * proc.sigma(energy) * speed
        if proc.process_type == ProcessType.ATTACHMENT:
            if cfg.nonconservative_model == "ignore":
                continue
            for i, val in enumerate(nu):
                mat[i, i] += -float(val)
            continue
        threshold = float(proc.threshold_eV)
        if proc.process_type == ProcessType.SUPERELASTIC and threshold > 0.0:
            threshold = -threshold
        if proc.process_type in {ProcessType.EXCITATION, ProcessType.SUPERELASTIC}:
            _add_energy_shift_transition(
                mat, energy, widths, nu, threshold, multiplicity=1.0
            )
        elif proc.process_type == ProcessType.IONIZATION:
            if (
                cfg.nonconservative_model == "ignore"
                or cfg.ionization_energy_sharing == "loss_only"
            ):
                _add_energy_shift_transition(
                    mat, energy, widths, nu, threshold, multiplicity=1.0
                )
            elif cfg.ionization_energy_sharing == "primary_secondary":
                _add_primary_secondary_ionization(
                    mat,
                    energy,
                    widths,
                    nu,
                    threshold,
                    cfg.secondary_electron_energy_eV,
                )
            else:
                _add_equal_sharing_ionization(mat, energy, widths, nu, threshold)
    return mat.tocsr()


__all__ = [
    "EffectiveCollisionData",
    "assemble_collision_operator",
    "build_effective_collision_data",
]
