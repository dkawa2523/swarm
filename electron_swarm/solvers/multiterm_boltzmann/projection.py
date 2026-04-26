"""Cross-section normalization and energy-grid projection."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from electron_swarm.core.constants import AMU_KG, ELECTRON_MASS_KG
from electron_swarm.core.cross_sections import (
    CrossSectionProcess,
    ProcessType,
    gas_mass_amu,
    mixture_fraction,
)
from electron_swarm.core.results import RateResult

from .models import MultiTermCase, RateSet

MOMENTUM_TYPES = {ProcessType.MOMENTUM, ProcessType.EFFECTIVE, ProcessType.ELASTIC}
INELASTIC_TYPES = {
    ProcessType.EXCITATION,
    ProcessType.IONIZATION,
    ProcessType.ATTACHMENT,
    ProcessType.SUPERELASTIC,
}


@dataclass(frozen=True, slots=True)
class ProcessNormalization:
    effective_species: tuple[str, ...]
    suppressed_momentum_processes: tuple[str, ...]
    warnings: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ProjectedProcess:
    process: CrossSectionProcess
    fraction: float
    sigma_m2: np.ndarray
    sigma_v_m3_s: np.ndarray
    frequency_s_inv: np.ndarray
    energy_loss_eV: float
    contributes_to_momentum: bool
    contributes_to_elastic_energy: bool


@dataclass(frozen=True, slots=True)
class ProjectedCollisionData:
    processes: tuple[ProjectedProcess, ...]
    momentum_frequency_s_inv: np.ndarray
    elastic_energy_frequency_s_inv: np.ndarray
    total_frequency_s_inv: np.ndarray
    normalization: ProcessNormalization


def _process_label(proc: CrossSectionProcess) -> str:
    return f"{proc.species}:{proc.process_type.value}:{proc.process}"


def _effective_species(case: MultiTermCase) -> set[str]:
    species: set[str] = set()
    for proc in case.cross_sections.processes:
        if proc.process_type != ProcessType.EFFECTIVE:
            continue
        if mixture_fraction(case.config.conditions, proc.species) > 0.0:
            species.add(proc.species)
    return species


def threshold_masked_sigma(proc: CrossSectionProcess, energy: np.ndarray) -> np.ndarray:
    sigma = proc.sigma(energy)
    if (
        proc.process_type != ProcessType.SUPERELASTIC
        and proc.threshold_eV is not None
        and proc.threshold_eV > 0.0
    ):
        sigma = np.where(energy >= proc.threshold_eV, sigma, 0.0)
    return np.asarray(sigma, dtype=float)


def energy_loss_eV(proc: CrossSectionProcess) -> float:
    if proc.process_type in {ProcessType.EXCITATION, ProcessType.IONIZATION}:
        return float(proc.threshold_eV or 0.0)
    if proc.process_type == ProcessType.SUPERELASTIC:
        return -abs(float(proc.threshold_eV or 0.0))
    return 0.0


def project_collision_data(case: MultiTermCase) -> ProjectedCollisionData:
    projected: list[ProjectedProcess] = []
    grid = case.grid
    momentum_frequency = np.zeros_like(grid.centers_eV)
    elastic_energy_frequency = np.zeros_like(grid.centers_eV)
    total_frequency = np.zeros_like(grid.centers_eV)
    effective_species = _effective_species(case)
    suppressed: list[str] = []

    for proc in case.cross_sections.processes:
        frac = mixture_fraction(case.config.conditions, proc.species)
        if frac <= 0.0:
            continue
        sigma = threshold_masked_sigma(proc, grid.centers_eV)
        sigma_v = sigma * grid.speeds_m_s
        frequency = case.gas_number_density_m3 * frac * sigma_v
        total_frequency += frequency

        contributes_to_momentum = False
        contributes_to_elastic_energy = False
        if proc.process_type in MOMENTUM_TYPES:
            if (
                proc.species in effective_species
                and proc.process_type != ProcessType.EFFECTIVE
            ):
                suppressed.append(_process_label(proc))
            else:
                contributes_to_momentum = True
                momentum_frequency += frequency
                if proc.process_type in {ProcessType.MOMENTUM, ProcessType.ELASTIC}:
                    contributes_to_elastic_energy = True
                    mass_amu = proc.mass_amu or gas_mass_amu(
                        case.config.conditions, proc.species
                    )
                    ratio = ELECTRON_MASS_KG / (mass_amu * AMU_KG)
                    elastic_energy_frequency += 2.0 * ratio * frequency
        elif proc.process_type in INELASTIC_TYPES:
            # Moment-closure drift estimate only. The operator uses explicit
            # isotropic source/sink projection from the per-process frequencies.
            contributes_to_momentum = True
            momentum_frequency += 0.2 * frequency

        projected.append(
            ProjectedProcess(
                process=proc,
                fraction=float(frac),
                sigma_m2=sigma,
                sigma_v_m3_s=sigma_v,
                frequency_s_inv=frequency,
                energy_loss_eV=energy_loss_eV(proc),
                contributes_to_momentum=contributes_to_momentum,
                contributes_to_elastic_energy=contributes_to_elastic_energy,
            )
        )

    warnings = ()
    if suppressed:
        warnings = (
            "Suppressed elastic/momentum transport contribution because an "
            "effective cross section exists for the same species.",
        )
    normalization = ProcessNormalization(
        effective_species=tuple(sorted(effective_species)),
        suppressed_momentum_processes=tuple(suppressed),
        warnings=warnings,
    )
    return ProjectedCollisionData(
        processes=tuple(projected),
        momentum_frequency_s_inv=momentum_frequency,
        elastic_energy_frequency_s_inv=elastic_energy_frequency,
        total_frequency_s_inv=total_frequency,
        normalization=normalization,
    )


def compute_rate_set(
    case: MultiTermCase,
    collisions: ProjectedCollisionData,
    energy_pdf_eV_inv: np.ndarray,
    case_id: str,
    solver_name: str,
) -> RateSet:
    grid = case.grid
    F = grid.normalize_energy_pdf(energy_pdf_eV_inv)
    out = []
    ion_freq = 0.0
    att_freq = 0.0
    for item in collisions.processes:
        proc = item.process
        with np.errstate(under="ignore"):
            k = float(np.sum(item.sigma_v_m3_s * F * grid.widths_eV))
        mixture_weighted = item.fraction * k
        freq = case.gas_number_density_m3 * mixture_weighted
        if proc.process_type == ProcessType.IONIZATION:
            ion_freq += freq
        elif proc.process_type == ProcessType.ATTACHMENT:
            att_freq += freq
        out.append(
            RateResult(
                solver=solver_name,
                case_id=case_id,
                e_over_n_Td=case.e_over_n_Td,
                species=proc.species,
                process=proc.process,
                process_type=proc.process_type.value,
                threshold_eV=proc.threshold_eV,
                rate_coefficient_m3_s=k,
                mixture_weighted_rate_m3_s=mixture_weighted,
                frequency_s_inv=freq,
                power_loss_eV_s=freq * item.energy_loss_eV,
            )
        )
    return RateSet(tuple(out), float(ion_freq), float(att_freq))


def total_momentum_frequency(
    case: MultiTermCase, collisions: ProjectedCollisionData, F: np.ndarray
) -> float:
    with np.errstate(under="ignore"):
        return float(
            np.sum(
                np.maximum(collisions.momentum_frequency_s_inv, 1.0)
                * F
                * case.grid.widths_eV
            )
        )


def elastic_power_loss_eV_s(
    case: MultiTermCase,
    collisions: ProjectedCollisionData,
    F: np.ndarray,
    thermal_mean_eV: float,
) -> float:
    with np.errstate(under="ignore"):
        return float(
            np.sum(
                collisions.elastic_energy_frequency_s_inv
                * (case.grid.centers_eV - thermal_mean_eV)
                * F
                * case.grid.widths_eV
            )
        )


def validity_warnings(
    case: MultiTermCase, collisions: ProjectedCollisionData, F: np.ndarray
) -> tuple[str, ...]:
    warnings: list[str] = list(collisions.normalization.warnings)
    tail_mask = case.grid.centers_eV >= 0.9 * case.grid.edges_eV[-1]
    with np.errstate(under="ignore"):
        tail = float(np.sum(F[tail_mask] * case.grid.widths_eV[tail_mask]))
        total_rate = float(
            np.sum(collisions.total_frequency_s_inv * F * case.grid.widths_eV)
        )
        tail_rate = float(
            np.sum(
                collisions.total_frequency_s_inv[tail_mask]
                * F[tail_mask]
                * case.grid.widths_eV[tail_mask]
            )
        )
    tail_rate_fraction = tail_rate / total_rate if total_rate > 0.0 else 0.0
    if tail > 1.0e-5 or tail_rate_fraction > 1.0e-3:
        warnings.append(
            "EEDF tail is non-negligible near max_eV; increase the energy grid."
        )
    extrapolated: list[str] = []
    for proc in case.cross_sections.processes:
        max_tabulated = float(proc.energy_eV[-1])
        if case.grid.edges_eV[-1] <= max_tabulated:
            continue
        above = case.grid.centers_eV > max_tabulated
        if not np.any(above):
            continue
        mass_above = float(np.sum(F[above] * case.grid.widths_eV[above]))
        if mass_above > 1.0e-8:
            policy = str(proc.metadata.get("high_energy_extrapolation", "hold"))
            extrapolated.append(
                f"{_process_label(proc)} above {max_tabulated:.3g} eV "
                f"(mass={mass_above:.3e}, policy={policy})"
            )
    if extrapolated:
        preview = "; ".join(extrapolated[:3])
        suffix = "" if len(extrapolated) <= 3 else f"; +{len(extrapolated) - 3} more"
        warnings.append(
            "EEDF samples energies beyond tabulated cross-section data: "
            f"{preview}{suffix}."
        )
    with np.errstate(under="ignore"):
        mean_e = float(np.sum(case.grid.centers_eV * F * case.grid.widths_eV))
    sig_total = 0.0
    for proc in case.cross_sections.processes:
        frac = mixture_fraction(case.config.conditions, proc.species)
        sig_total += frac * float(threshold_masked_sigma(proc, np.array([mean_e]))[0])
    if sig_total > 0.0:
        mfp = 1.0 / (case.gas_number_density_m3 * sig_total)
        length_scale = case.config.conditions.length_scale_m
        if length_scale is not None and mfp / length_scale > 0.1:
            warnings.append(
                f"Mean free path ratio lambda/L={mfp / length_scale:.3e}; "
                "check local-field and hydrodynamic assumptions."
            )
        elif length_scale is None and mfp > 0.01:
            warnings.append(
                f"Mean free path is large ({mfp:.3e} m); check hydrodynamic locality."
            )
    return tuple(warnings)
