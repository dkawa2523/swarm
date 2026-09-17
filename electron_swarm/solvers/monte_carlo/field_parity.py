"""Common-random-number field-parity response for Monte Carlo transport."""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np

from electron_swarm.core.results import SwarmCaseResult
import electron_swarm.solvers.monte_carlo.evidence as _mc_evidence


_ODD_RESPONSE_FIELDS = (
    "drift_velocity_m_s",
    "mobility_m2_V_s",
    "energy_mobility_m2_V_s",
)
_STATIONARITY_KEYS = (
    "time_stationarity_window_estimate_early",
    "time_stationarity_window_estimate_late",
)


def _paired_value(primary: float, mirror: float) -> float:
    value = 0.5 * (float(primary) + float(mirror))
    if not np.isfinite(value):
        raise RuntimeError("field-parity transport response is non-finite")
    return value


def _odd_response_values(result: SwarmCaseResult) -> dict[str, float]:
    energy_mobility = result.transport.electron_energy_mobility_m2_V_s
    if energy_mobility is None:
        raise RuntimeError("field-parity response requires energy mobility")
    return {
        "drift_velocity_m_s": float(result.transport.drift_velocity_m_s),
        "mobility_m2_V_s": float(result.transport.mobility_m2_V_s),
        "energy_mobility_m2_V_s": float(energy_mobility),
    }


def _leg_status(diagnostics: dict[str, Any]) -> dict[str, bool]:
    return {
        key: bool(diagnostics.get(key, False))
        for key in (
            "finite",
            "strictly_positive",
            "complete_observation",
            "multiple_time_origins",
            "lineage_qualified",
        )
        if key in diagnostics
    }


def apply_field_parity_response(
    primary: SwarmCaseResult,
    mirror: SwarmCaseResult,
) -> SwarmCaseResult:
    """Replace only odd transport responses by the paired central response.

    Field inversion symmetry makes the arithmetic mean the central response
    without adding algebraic estimator bias.  Each finite-warmup leg can still
    carry transient bias, which remains visible to the ordinary stationarity
    gate.  Common random numbers reduce variance while paths remain correlated;
    collision branching can remove that covariance.  Even observables, EEDF,
    and reaction rates remain owned by the positive-field primary trajectory.
    """

    if (
        primary.solver != "monte_carlo"
        or mirror.solver != "monte_carlo"
        or primary.e_over_n_Td != mirror.e_over_n_Td
        or primary.transport.definition != mirror.transport.definition
        or primary.transport.gas_number_density_m3
        != mirror.transport.gas_number_density_m3
    ):
        raise ValueError("field-parity legs must describe the same Monte Carlo case")

    primary_values = _odd_response_values(primary)
    mirror_values = _odd_response_values(mirror)
    paired = {
        name: _paired_value(primary_values[name], mirror_values[name])
        for name in _ODD_RESPONSE_FIELDS
    }
    density = primary.transport.gas_number_density_m3
    primary.transport = replace(
        primary.transport,
        drift_velocity_m_s=paired["drift_velocity_m_s"],
        reduced_mobility_m2_V_s_m3=paired["mobility_m2_V_s"] * density,
        reduced_electron_energy_mobility_m2_V_s_m3=(
            paired["energy_mobility_m2_V_s"] * density
        ),
    )
    primary.effective_townsend_m2 = primary.net_ionization_frequency_s / max(
        abs(paired["drift_velocity_m_s"]) * density,
        1.0e-300,
    )

    primary_diagnostics = primary.diagnostics["internal_monte_carlo_transport"]
    mirror_diagnostics = mirror.diagnostics["internal_monte_carlo_transport"]
    production = primary_diagnostics["production_estimates"]
    production.update(paired)
    paired_windows: dict[str, dict[str, float]] = {}
    for key in _STATIONARITY_KEYS:
        primary_window = primary_diagnostics.get(key)
        mirror_window = mirror_diagnostics.get(key)
        if not isinstance(primary_window, dict) or not isinstance(mirror_window, dict):
            continue
        window_pair = {
            name: _paired_value(primary_window[name], mirror_window[name])
            for name in _ODD_RESPONSE_FIELDS
            if name in primary_window and name in mirror_window
        }
        primary_window.update(window_pair)
        paired_windows[key] = window_pair

    primary_status = _leg_status(primary_diagnostics)
    mirror_status = _leg_status(mirror_diagnostics)
    candidate_status = "direct_mc_candidate_for_ensemble_qualification"
    primary_candidate = primary_diagnostics.get("energy_transport_status") == (
        candidate_status
    )
    mirror_candidate = mirror_diagnostics.get("energy_transport_status") == (
        candidate_status
    )
    for key in primary_status.keys() & mirror_status.keys():
        primary_diagnostics[key] = primary_status[key] and mirror_status[key]
    primary_diagnostics["field_parity_response"] = {
        "estimator_schema_version": (
            _mc_evidence.FIELD_PARITY_RESPONSE_ESTIMATOR_SCHEMA_VERSION
        ),
        "model": "paired_positive_negative_field_odd_response",
        "rng_coupling": "identical_bit_generator_state_at_leg_initialization",
        "pair_weighting": "equal_arithmetic_mean_of_field_aligned_leg_estimates",
        "response_fields": list(_ODD_RESPONSE_FIELDS),
        "primary_positive_field": primary_values,
        "mirror_negative_field": mirror_values,
        "paired_response": paired,
        "paired_stationarity_windows": paired_windows,
        "primary_leg_status": primary_status,
        "mirror_leg_status": mirror_status,
        "eedf_rates_even_transport_source": "primary_positive_field_leg",
        "primary_sampling_stream": (
            "bitwise_identical_to_single_field_same_seed_and_budget"
        ),
        "mirror_tail_sampling": False,
        "symmetry_bias_scope": (
            "no_additional_algebraic_bias_under_field_inversion_symmetry"
        ),
        "finite_time_bias_control": "existing_per_leg_stationarity_evidence",
        "decorrelation_effect": "can_reduce_covariance_benefit",
    }
    provenance = primary_diagnostics.get("mc_run_provenance")
    if isinstance(provenance, dict):
        provenance.update(
            {
                "transport_estimator": "paired_field_parity",
                "transport_field_legs": 2,
                "transport_budget_interpretation": "configured_budget_per_field_leg",
                "batch_rng_advance": "primary_positive_field_leg_only",
                "primary_sampling_stream": (
                    "bitwise_identical_to_single_field_same_seed_and_budget"
                ),
            }
        )
    primary_diagnostics["energy_transport_status"] = (
        candidate_status
        if primary_candidate
        and mirror_candidate
        and all(value > 0.0 for value in paired.values())
        else "direct_mc_unqualified"
    )
    return primary


__all__ = ["apply_field_parity_response"]
