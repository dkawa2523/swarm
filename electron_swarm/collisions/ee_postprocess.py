"""Electron-electron post-processing for Boltzmann-family product solvers."""

from __future__ import annotations

from typing import Iterable

import numpy as np

from electron_swarm.collisions.ee_fp_energy import apply_fp_energy_operator
from electron_swarm.collisions.electron_electron import mean_energy_eV, relaxation_target
from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.core.numerics import widths_from_centers
from electron_swarm.solvers.kinetic import (
    compute_rates_from_eedf,
    eepf_from_eedf,
    gas_number_density,
)


_BOLTZMANN_SOLVERS = {
    "two_term",
    "multi_term",
}


def _ee_config(config: object) -> object:
    return config.physics.electron_electron


def _normalise(
    energy_eV: np.ndarray,
    eedf: np.ndarray,
    widths_eV: np.ndarray | None = None,
) -> np.ndarray:
    widths = widths_from_centers(energy_eV) if widths_eV is None else widths_eV
    total = float(np.sum(np.clip(eedf, 0.0, None) * widths))
    if total <= 0.0 or not np.isfinite(total):
        raise ValueError("Cannot normalize e-e relaxed EEDF")
    return np.clip(eedf, 0.0, None) / total


def _case_widths(case: SwarmCaseResult) -> np.ndarray:
    if case.energy_widths_eV is not None:
        widths = np.asarray(case.energy_widths_eV, dtype=float)
        if widths.shape == np.asarray(case.energy_eV).shape and np.all(widths > 0.0):
            return widths
    return widths_from_centers(case.energy_eV)


def _update_eepf(case: SwarmCaseResult) -> None:
    case.eepf = eepf_from_eedf(case.energy_eV, case.eedf)


def _recompute_rates(
    case: SwarmCaseResult,
    config: object,
    cross_sections: CrossSectionSet | None,
) -> bool:
    if cross_sections is None:
        return False
    energy = np.asarray(case.energy_eV, dtype=float)
    eedf = np.asarray(case.eedf, dtype=float)
    widths = _case_widths(case)
    if len(energy) == 0 or len(eedf) != len(energy):
        return False

    convolution = compute_rates_from_eedf(
        config,
        cross_sections,
        energy,
        widths,
        eedf,
        case_id=case.case_id,
        e_over_n_Td=case.e_over_n_Td,
        solver_name=case.solver,
    )
    net_freq = convolution.net_ionization_frequency_s
    number_density = gas_number_density(config)
    case.rates = convolution.rates
    case.net_ionization_frequency_s = net_freq
    case.effective_townsend_m2 = net_freq / max(
        abs(case.drift_velocity_m_s) * number_density,
        1.0e-300,
    )
    return True


def _apply_to_case(
    case: SwarmCaseResult,
    ee_config: object,
    config: object,
    cross_sections: CrossSectionSet | None,
) -> SwarmCaseResult:
    if case.solver not in _BOLTZMANN_SOLVERS:
        raise ValueError(
            f"{case.solver} cannot apply electron_electron postprocess; "
            "solver plan should have failed or skipped this case"
        )
    energy = np.asarray(case.energy_eV, dtype=float)
    old = np.asarray(case.eedf, dtype=float)
    widths = _case_widths(case)
    if len(energy) == 0 or len(old) != len(energy):
        raise ValueError("Cannot apply electron_electron relaxation to invalid EEDF")

    alpha = float(ee_config.relaxation_fraction)
    model = str(ee_config.model).lower()
    operator_metadata: dict[str, object] = {}
    if model == "relaxation_postprocess":
        conserve = bool(ee_config.conserve_mean_energy)
        fallback_te = float(ee_config.fallback_temperature_eV)
        target = relaxation_target(
            energy,
            widths,
            old,
            conserve_mean_energy=conserve,
            fallback_temperature_eV=fallback_te,
        )
        new = _normalise(energy, (1.0 - alpha) * old + alpha * target, widths)
    elif model == "fp_energy":
        if ee_config.strength_model == "density_based":
            raise NotImplementedError(
                "electron_electron fp_energy strength_model='density_based' "
                "is not implemented"
            )
        new, operator_metadata = apply_fp_energy_operator(
            energy,
            widths,
            old,
            relaxation_fraction=ee_config.relaxation_fraction,
            conserve_mean_energy=ee_config.conserve_mean_energy,
            fallback_temperature_eV=ee_config.fallback_temperature_eV,
        )
    else:
        raise ValueError("Unsupported electron_electron model for postprocess")
    after = mean_energy_eV(energy, widths, new)

    case.eedf = new
    case.mean_energy_eV = after
    _update_eepf(case)
    _recompute_rates(case, config, cross_sections)
    case.diagnostics["electron_electron"] = operator_metadata
    case.metadata.update(
        {
            "electron_electron_treatment": model,
            "electron_electron_transport_stale": True,
        }
    )
    return case


def apply_electron_electron_relaxation_from_config(
    cases: Iterable[SwarmCaseResult],
    config: object,
    cross_sections: CrossSectionSet | None = None,
) -> list[SwarmCaseResult]:
    """Apply optional e-e relaxation to Boltzmann-family cases only."""

    ee_config = _ee_config(config)
    if not bool(ee_config.enabled):
        return list(cases)
    model = str(ee_config.model).lower()
    if model not in {"relaxation_postprocess", "fp_energy"}:
        raise ValueError("Unsupported electron_electron model for postprocess")
    return [_apply_to_case(case, ee_config, config, cross_sections) for case in cases]
