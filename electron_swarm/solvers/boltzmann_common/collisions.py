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
    CrossSectionSet,
    ProcessType,
    ScatteringRole,
    gas_mass_amu,
    mixture_fraction,
)
from electron_swarm.core.scattering import transport_momentum_processes
from electron_swarm.core.solver_configs import (
    InelasticCollisionConfig,
    MomentumCollisionConfig,
)
from electron_swarm.physics.angular_scattering import build_angular_model
from electron_swarm.physics.electron_neutral import ionization_daughter_energy_arrays
from electron_swarm.solvers.boltzmann_common.grid import electron_speed_m_s


@dataclass(slots=True)
class EffectiveCollisionData:
    nu_m: np.ndarray
    sigma_m: np.ndarray
    inelastic_loss_frequency_s_inv: np.ndarray
    elastic_A_eV_s: np.ndarray
    elastic_D_eV2_s: np.ndarray


@dataclass(slots=True)
class InelasticCollisionData:
    matrix: sparse.csr_matrix
    net_creation_frequency_s_inv: np.ndarray


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
    recoil_integrals: dict[int, np.ndarray] = {}
    effective_species: set[str] = set()
    for species in cross_sections.species:
        if mixture_fraction(config.conditions, species) <= 0.0:
            continue
        selected, includes_inelastic, _model = transport_momentum_processes(
            cross_sections,
            species,
        )
        if not selected:
            continue
        momentum_ids.update(id(process) for process in selected)
        if includes_inelastic:
            effective_species.add(species)
        species_processes = cross_sections.by_species(species)
        elastic_momentum = [
            process
            for process in species_processes
            if process.scattering_role == ScatteringRole.ELASTIC_MOMENTUM_TRANSFER
        ]
        elastic_total = [
            process
            for process in species_processes
            if process.scattering_role == ScatteringRole.ELASTIC_TOTAL
        ]
        angular_name = config.physics.angular_scattering.model
        if angular_name == "moment_table" and elastic_total:
            # The PN angular damping uses this same total and first moment.
            # Neither an effective integral nor an unrelated transport total
            # can own the neutral recoil frequency.
            first_moment = build_angular_model(config).moments(energy, 1)[1]
            if elastic_momentum:
                supplied = sum(process.sigma(energy) for process in elastic_momentum)
                reconstructed = sum(
                    process.sigma(energy) for process in elastic_total
                ) * (1.0 - first_moment)
                if not np.allclose(supplied, reconstructed, rtol=1.0e-6, atol=1.0e-30):
                    raise ValueError(
                        f"{species}: elastic momentum-transfer cross sections "
                        "disagree with elastic total and supplied first angular "
                        "moment (rtol=1e-6, atol=1e-30 m^2)"
                    )
                for process in elastic_momentum:
                    recoil_integrals[id(process)] = process.sigma(energy)
            else:
                for process in elastic_total:
                    recoil_integrals[id(process)] = process.sigma(energy) * (
                        1.0 - first_moment
                    )
        elif elastic_momentum:
            for process in elastic_momentum:
                recoil_integrals[id(process)] = process.sigma(energy)
        elif angular_name == "isotropic" and elastic_total:
            for process in elastic_total:
                recoil_integrals[id(process)] = process.sigma(energy)
        else:
            raise NotImplementedError(
                f"{species}: elastic energy relaxation requires an elastic "
                "momentum-transfer cross section, or an elastic total with "
                "an isotropic closure or supplied first angular moment; "
                "effective momentum integrals include inelastic loss and "
                "cannot define elastic recoil"
            )

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
            raise ValueError(f"Cross section {proc.species}:{proc.process} is negative")
        include_transport = id(proc) in momentum_ids
        if proc.process_type in transport_inelastic:
            include_transport = proc.species not in effective_species
            inelastic_loss += gas_number_density_m3 * frac * sigma * speed
        if include_transport:
            sigma_m += frac * sigma
        recoil_sigma = recoil_integrals.get(id(proc))
        if recoil_sigma is None:
            continue
        mass_amu = proc.mass_amu or gas_mass_amu(config.conditions, proc.species)
        mass_kg = mass_amu * AMU_KG
        nu = gas_number_density_m3 * frac * recoil_sigma * speed
        ratio = ELECTRON_MASS_KG / mass_kg
        eps = np.maximum(energy, 1.0e-12)
        elastic_D += 2.0 * ratio * nu * eps * kT_eV
        elastic_A += 2.0 * ratio * nu * (0.5 * kT_eV - eps)
    sigma_m = np.maximum(sigma_m, cfg.min_momentum_cross_section_m2)
    nu_m = np.maximum(gas_number_density_m3 * speed * sigma_m, 1.0e-60)
    return EffectiveCollisionData(
        nu_m=nu_m,
        sigma_m=sigma_m,
        inelastic_loss_frequency_s_inv=inelastic_loss,
        elastic_A_eV_s=elastic_A,
        elastic_D_eV2_s=np.maximum(elastic_D, 0.0),
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


def _add_energy_transfer(
    mat: sparse.lil_matrix,
    energy: np.ndarray,
    widths: np.ndarray,
    nu: np.ndarray,
    daughter_energies: tuple[np.ndarray, ...],
) -> None:
    for i, rate in enumerate(np.asarray(nu, dtype=float)):
        if rate <= 0.0:
            continue
        mat[i, i] += -rate
        for daughter in daughter_energies:
            _deposit_energy(mat, energy, widths, i, float(daughter[i]), rate)


def build_inelastic_collision_data(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    energy: np.ndarray,
    widths: np.ndarray,
    gas_number_density_m3: float,
    solver_config: InelasticCollisionConfig,
) -> InelasticCollisionData:
    n = len(energy)
    mat = sparse.lil_matrix((n, n), dtype=float)
    net_creation = np.zeros(n, dtype=float)
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
            net_creation -= nu
            continue
        threshold = float(proc.threshold_eV)
        if proc.process_type == ProcessType.SUPERELASTIC and threshold > 0.0:
            threshold = -threshold
        if proc.process_type in {ProcessType.EXCITATION, ProcessType.SUPERELASTIC}:
            _add_energy_transfer(mat, energy, widths, nu, (energy - threshold,))
        elif proc.process_type == ProcessType.IONIZATION:
            model = (
                "loss_only"
                if cfg.nonconservative_model == "ignore"
                else cfg.ionization_energy_sharing
            )
            daughters = ionization_daughter_energy_arrays(
                energy,
                threshold,
                model=model,
                secondary_electron_energy_eV=cfg.secondary_electron_energy_eV,
            )
            _add_energy_transfer(mat, energy, widths, nu, daughters)
            net_creation += (len(daughters) - 1) * nu
    return InelasticCollisionData(mat.tocsr(), net_creation)


__all__ = [
    "EffectiveCollisionData",
    "InelasticCollisionData",
    "build_effective_collision_data",
    "build_inelastic_collision_data",
]
