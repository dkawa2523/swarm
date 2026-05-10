"""Optional electron-electron relaxation post-processing.

This is the first, merge-safe e-e collision hook.  It deliberately avoids a full
Fokker-Planck operator and does not touch Monte Carlo cases.  When enabled in
YAML it relaxes Boltzmann-family EEDFs slightly toward a Maxwellian target that
preserves the mean energy unless configured otherwise.
"""

from __future__ import annotations

from typing import Any, Iterable

import numpy as np

from electron_swarm.collisions.config import raw_config_from_source
from electron_swarm.collisions.electron_electron import mean_energy_eV, relaxation_target
from electron_swarm.core.constants import BOLTZMANN_J_K, ELECTRON_MASS_KG, EV_TO_J
from electron_swarm.core.cross_sections import (
    CrossSectionSet,
    ProcessType,
    mixture_fraction,
)
from electron_swarm.core.results import RateResult, SwarmCaseResult
from electron_swarm.diagnostics.common import widths_from_centers


_BOLTZMANN_SOLVERS = {"boltzmann_two_term", "multiterm_boltzmann"}


def _ee_section(raw: dict[str, Any]) -> dict[str, Any]:
    section = raw.get("electron_electron", raw.get("electron_electron_collision", {}))
    return section if isinstance(section, dict) else {}


def _normalise(energy_eV: np.ndarray, eedf: np.ndarray) -> np.ndarray:
    widths = widths_from_centers(energy_eV)
    total = float(np.sum(np.clip(eedf, 0.0, None) * widths))
    if total <= 0.0 or not np.isfinite(total):
        raise ValueError("Cannot normalize e-e relaxed EEDF")
    return np.clip(eedf, 0.0, None) / total


def _update_eepf(case: SwarmCaseResult) -> None:
    energy = np.asarray(case.energy_eV, dtype=float)
    denom = np.sqrt(np.maximum(energy, 1.0e-30))
    try:
        case.eepf = np.asarray(case.eedf, dtype=float) / denom
    except Exception:
        # Some result implementations may not expose mutable eepf.  The EEDF
        # and metadata are still updated; diagnostics will flag the case.
        case.metadata["ee_eepf_recomputed"] = False
    else:
        case.metadata["ee_eepf_recomputed"] = True


def _gas_number_density(config: object) -> float:
    cond = getattr(config, "conditions")
    if cond.gas_number_density_m3 is not None:
        return float(cond.gas_number_density_m3)
    return float(cond.pressure_Pa / (BOLTZMANN_J_K * cond.gas_temperature_K))


def _energy_loss_eV(process_type: ProcessType, threshold_eV: float | None) -> float:
    if process_type in {ProcessType.EXCITATION, ProcessType.IONIZATION}:
        return float(threshold_eV or 0.0)
    if process_type == ProcessType.SUPERELASTIC:
        return -abs(float(threshold_eV or 0.0))
    return 0.0


def _recompute_rates(
    case: SwarmCaseResult,
    config: object,
    cross_sections: CrossSectionSet | None,
) -> bool:
    if cross_sections is None:
        return False
    energy = np.asarray(case.energy_eV, dtype=float)
    eedf = np.asarray(case.eedf, dtype=float)
    widths = widths_from_centers(energy)
    if len(energy) == 0 or len(eedf) != len(energy):
        return False

    speed = np.sqrt(np.maximum(2.0 * EV_TO_J * energy / ELECTRON_MASS_KG, 0.0))
    number_density = _gas_number_density(config)
    rates: list[RateResult] = []
    ion_rate = 0.0
    attach_rate = 0.0
    for proc in cross_sections.processes:
        frac = mixture_fraction(config.conditions, proc.species)
        if frac <= 0.0:
            continue
        k = float(np.sum(proc.sigma(energy) * speed * eedf * widths))
        kmix = frac * k
        freq = number_density * kmix
        if proc.process_type == ProcessType.IONIZATION:
            ion_rate += kmix
        elif proc.process_type == ProcessType.ATTACHMENT:
            attach_rate += kmix
        rates.append(
            RateResult(
                solver=case.solver,
                case_id=case.case_id,
                e_over_n_Td=case.e_over_n_Td,
                species=proc.species,
                process=proc.process,
                process_type=proc.process_type.value,
                threshold_eV=proc.threshold_eV,
                rate_coefficient_m3_s=k,
                mixture_weighted_rate_m3_s=kmix,
                frequency_s_inv=freq,
                power_loss_eV_s=freq
                * _energy_loss_eV(proc.process_type, proc.threshold_eV),
            )
        )
    net_freq = number_density * (ion_rate - attach_rate)
    case.rates = rates
    case.net_ionization_frequency_s = net_freq
    case.effective_townsend_m2 = net_freq / max(
        abs(case.drift_velocity_m_s) * number_density,
        1.0e-300,
    )
    case.metadata.update(
        {
            "convolution_ionization_rate_coefficient_m3_s": ion_rate,
            "convolution_attachment_rate_coefficient_m3_s": attach_rate,
            "convolution_effective_rate_coefficient_m3_s": ion_rate - attach_rate,
            "convolution_net_ionization_frequency_s-1": net_freq,
            "effective_townsend_1_m": net_freq
            / max(abs(case.drift_velocity_m_s), 1.0e-300),
        }
    )
    return True


def _apply_to_case(
    case: SwarmCaseResult,
    section: dict[str, Any],
    config: object,
    cross_sections: CrossSectionSet | None,
) -> SwarmCaseResult:
    if case.solver not in _BOLTZMANN_SOLVERS:
        case.metadata.setdefault("ee_relaxation_skipped", "non_boltzmann_solver")
        return case
    energy = np.asarray(case.energy_eV, dtype=float)
    old = np.asarray(case.eedf, dtype=float)
    widths = widths_from_centers(energy)
    if len(energy) == 0 or len(old) != len(energy):
        case.metadata["ee_relaxation_skipped"] = "invalid_eedf_shape"
        return case

    alpha = section.get("relaxation_fraction", section.get("strength_scale", 0.05))
    alpha = float(np.clip(float(alpha), 0.0, 1.0))
    model = str(section.get("model", "relaxation")).lower()
    conserve = bool(section.get("conserve_mean_energy", True))
    fallback_te = float(
        section.get("fallback_temperature_eV", section.get("temperature_eV", 2.0))
    )
    before = mean_energy_eV(energy, widths, old)
    target = relaxation_target(
        energy,
        widths,
        old,
        conserve_mean_energy=conserve,
        fallback_temperature_eV=fallback_te,
    )
    new = _normalise(energy, (1.0 - alpha) * old + alpha * target)
    after = mean_energy_eV(energy, widths, new)

    case.eedf = new
    case.mean_energy_eV = after
    _update_eepf(case)
    rates_recomputed = _recompute_rates(case, config, cross_sections)
    case.metadata.update(
        {
            "ee_enabled": True,
            "ee_model": model,
            "ee_model_implementation": "relaxation_postprocess",
            "ee_relaxation_fraction": alpha,
            "ee_conserve_mean_energy": conserve,
            "ee_mean_energy_before_eV": before,
            "ee_mean_energy_after_eV": after,
            "ee_rates_recomputed": rates_recomputed,
            "ee_transport_recomputed": False,
            "ee_transport_stale_after_relaxation": True,
        }
    )
    return case


def apply_electron_electron_relaxation_from_config(
    cases: Iterable[SwarmCaseResult],
    config: object,
    cross_sections: CrossSectionSet | None = None,
) -> list[SwarmCaseResult]:
    """Apply optional e-e relaxation to Boltzmann-family cases only."""

    raw = raw_config_from_source(config)
    section = _ee_section(raw)
    if not bool(section.get("enabled", False)):
        return list(cases)
    if str(section.get("model", "relaxation")).lower() not in {
        "relaxation",
        "relaxation_postprocess",
    }:
        return list(cases)
    return [_apply_to_case(case, section, config, cross_sections) for case in cases]
