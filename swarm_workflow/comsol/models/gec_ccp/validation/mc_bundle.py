"""Validate Monte Carlo bundle quality, censoring, and interpolation."""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.interpolate import PchipInterpolator


from swarm_workflow.comsol.models.gec_ccp.closure import (
    _mc_nonnegative_rate_table,
    _reaction_uses_preintegrated_rate,
    _uses_direct_mc_rates,
    _uses_function_eedf,
)
from swarm_workflow.comsol.models.gec_ccp.contracts import (
    GecCcpMapping,
    GecCcpWorkflowError,
)
from swarm_workflow.comsol.models.gec_ccp.data import (
    _lookup,
    _read_csv,
    _required_float,
)
from swarm_workflow.tables.contracts import (
    ZERO_EVENT_CONFIDENCE,
)


REQUIRED_MC_QUALITY_COLUMNS = (
    "aggregate_quality_passed",
    "failure_reasons_json",
    "aggregate_failure_reasons_json",
    "mobility_rse",
    "diffusion_L_rse",
    "diffusion_T_rse",
    "energy_mobility_rse",
    "energy_diffusion_L_rse",
    "energy_diffusion_T_rse",
    "max_major_rate_rse",
    "valid_replicates",
    "uncertainty_available",
    "solver_diagnostics_available",
    "solver_transport_qualified",
    "solver_transport_replicates",
    "solver_origin_stationarity_limiting_field",
    "solver_origin_stationarity_limiting_relative_ci95_bound",
    "solver_mean_energy_stationarity_relative_ci95_bound",
    "solver_mobility_stationarity_absolute_log_drift",
    "solver_lag_convergence_max_relative_ci95_bound",
    "solver_transport_mean_energy_max_relative_ci95_bound",
    "solver_population_growth_gate_mode",
    "solver_population_growth_max_relative_ci95_bound",
    "solver_population_growth_poisson_interval_max_ratio",
    "solver_population_growth_sparse_max_metric",
    "solver_origin_stationarity_limiting_relative_tolerance",
    "solver_lag_convergence_relative_tolerance",
    "solver_transport_mean_energy_relative_tolerance",
    "solver_population_growth_relative_tolerance",
    "quality_source",
)


def _mc_rate_censoring_audit(
    mapping: GecCcpMapping,
    rate_rows: list[dict[str, Any]],
    *,
    required_rate_min_process_peak_fraction: float,
) -> dict[str, Any] | None:
    """Validate MC rates and zero-event tails against raw event evidence.

    Direct-rate closures use the nominal rates as active COMSOL inputs.  A
    Function-EEDF closure keeps the same rates as independent evidence for the
    sampled EEDF tail.  Both uses require the same censoring qualification.
    """

    if mapping.bundle.expected_source != "monte_carlo" or not (
        _uses_direct_mc_rates(mapping) or _uses_function_eedf(mapping.closure)
    ):
        return None
    evidence_only = not _uses_direct_mc_rates(mapping)
    relevance_fraction = float(required_rate_min_process_peak_fraction)
    if (
        not math.isfinite(relevance_fraction)
        or relevance_fraction < 0.0
        or relevance_fraction > 1.0
    ):
        raise GecCcpWorkflowError(
            "Monte Carlo censored-rate relevance fraction must be in [0, 1]"
        )
    evidence_rows = _read_csv(mapping.bundle.path / "rate_evidence.csv")
    process_results: dict[str, Any] = {}
    censored_anchors: list[dict[str, Any]] = []
    failed_censored_anchors: list[dict[str, Any]] = []
    for reaction in mapping.reactions:
        if evidence_only:
            selected = reaction.process_type in {"excitation", "ionization"}
        else:
            selected = _reaction_uses_preintegrated_rate(mapping.closure, reaction)
        if not selected:
            continue
        active, _ = _mc_nonnegative_rate_table(
            rate_rows, process_type=reaction.process_type
        )
        process_peak = max(
            _required_float(row, "rate_coefficient_m3_s") for row in active
        )
        if not math.isfinite(process_peak) or process_peak <= 0.0:
            raise GecCcpWorkflowError(
                "Monte Carlo censored-rate qualification requires a positive "
                f"process peak: {reaction.process_type}"
            )
        upper_limit = relevance_fraction * process_peak
        evidence = [
            row
            for row in evidence_rows
            if row.get("process_type") == reaction.process_type
        ]
        checked = 0
        process_censored: list[dict[str, Any]] = []
        for rate_row in active:
            mean_energy = _required_float(rate_row, "mean_energy_eV")
            rate = _required_float(rate_row, "rate_coefficient_m3_s")
            matches = [
                row
                for row in evidence
                if math.isclose(
                    _required_float(row, "mean_energy_eV"),
                    mean_energy,
                    rel_tol=1.0e-12,
                    abs_tol=1.0e-14,
                )
            ]
            if len(matches) != 1:
                raise GecCcpWorkflowError(
                    "Monte Carlo rate evidence must contain exactly one row "
                    f"for {reaction.process_type} at {mean_energy:.17g} eV"
                )
            evidence_row = matches[0]
            evidence_rate = _required_float(evidence_row, "rate_coefficient_mean_m3_s")
            if not math.isclose(
                evidence_rate,
                rate,
                rel_tol=1.0e-12,
                abs_tol=0.0,
            ):
                raise GecCcpWorkflowError(
                    "Monte Carlo canonical rate differs from raw trajectory "
                    f"evidence for {reaction.process_type} at "
                    f"{mean_energy:.17g} eV"
                )
            checked += 1
            if rate != 0.0:
                continue
            pooled_event_count = _required_float(evidence_row, "pooled_event_count")
            pooled_exposure = _required_float(
                evidence_row, "pooled_target_exposure_s_m3"
            )
            upper = _required_float(evidence_row, "pooled_all_zero_upper_95_m3_s")
            recomputed_upper = (
                -math.log(1.0 - ZERO_EVENT_CONFIDENCE) / pooled_exposure
                if pooled_exposure > 0.0
                else math.nan
            )
            if (
                evidence_row.get("estimate_status") != "censored_all_zero"
                or evidence_row.get("pooled_zero_event_status") != "all_zero_upper_95"
                or pooled_event_count != 0.0
                or pooled_exposure <= 0.0
                or upper <= 0.0
                or not math.isfinite(recomputed_upper)
                or not math.isclose(
                    upper,
                    recomputed_upper,
                    rel_tol=1.0e-12,
                    abs_tol=0.0,
                )
                or str(evidence_row.get("uncertainty_available", "")).lower()
                not in {"0", "false"}
            ):
                raise GecCcpWorkflowError(
                    "zero-valued Monte Carlo rate lacks a self-consistent "
                    f"pooled 95% zero-event upper bound: {reaction.process_type} at "
                    f"{mean_energy:.17g} eV"
                )
            qualified = upper < upper_limit
            item = {
                "process_type": reaction.process_type,
                "mean_energy_eV": mean_energy,
                "nominal_rate_m3_s": 0.0,
                "pooled_event_count": 0,
                "pooled_target_exposure_s_m3": pooled_exposure,
                "pooled_all_zero_upper_95_m3_s": upper,
                "process_peak_m3_s": process_peak,
                "required_rate_min_process_peak_fraction": relevance_fraction,
                "required_upper_limit_m3_s": upper_limit,
                "upper_to_process_peak_fraction": upper / process_peak,
                "estimate_status": "censored_all_zero",
                "passed": qualified,
                "criterion": (
                    "pooled_upper_95_strictly_below_required_fraction_of_process_peak"
                ),
            }
            process_censored.append(item)
            censored_anchors.append(item)
            if not qualified:
                failed_censored_anchors.append(item)
        process_results[reaction.process_type] = {
            "process_peak_m3_s": process_peak,
            "required_rate_min_process_peak_fraction": relevance_fraction,
            "required_upper_limit_m3_s": upper_limit,
            "rate_rows_checked": checked,
            "censored_anchor_count": len(process_censored),
            "qualified_censored_anchor_count": sum(
                int(item["passed"]) for item in process_censored
            ),
            "failed_censored_anchor_count": sum(
                int(not item["passed"]) for item in process_censored
            ),
            "censored_anchors": process_censored,
        }
    physics_passed = not failed_censored_anchors
    return {
        "status": "passed" if physics_passed else "failed",
        "nominal_numerical_binding": {
            "passed": True,
            "estimator": "monte_carlo_trajectory_time_average_sigma_v",
            "role": (
                "independent_function_eedf_tail_evidence"
                if evidence_only
                else "active_preintegrated_rate_input"
            ),
            "active_comsol_rate_input": not evidence_only,
            "zero_estimates_bound_exactly": True,
            "artificial_floor": False,
            "two_term_substitution": False,
        },
        "processes": process_results,
        "censored_anchor_count": len(censored_anchors),
        "failed_censored_anchor_count": len(failed_censored_anchors),
        "physics_qualification": {
            "passed": physics_passed,
            "status": (
                "qualified_no_censored_rate_anchors"
                if not censored_anchors
                else (
                    "qualified_censored_rate_upper_bounds_below_relevance_floor"
                    if physics_passed
                    else "unqualified_censored_rate_upper_bound_not_below_"
                    "relevance_floor"
                )
            ),
            "reason": (
                None
                if physics_passed
                else "at least one pooled 95% zero-event upper bound is not "
                "strictly below the configured fraction of its process peak"
            ),
            "pooled_upper_bound_relevance_gate_performed": True,
            "upper_bound_sensitivity_performed": False,
            "censored_anchors": censored_anchors,
            "failed_censored_anchors": failed_censored_anchors,
        },
    }


def _mc_rate_interpolation_audit(
    mapping: GecCcpMapping,
    rate_rows: list[dict[str, Any]],
) -> dict[str, Any] | None:
    """Independently qualify the nonnegative normalized MC PCHIP tables."""

    if not _uses_direct_mc_rates(mapping):
        return None
    processes: dict[str, Any] = {}
    passed = True
    for reaction in mapping.reactions:
        if not _reaction_uses_preintegrated_rate(mapping.closure, reaction):
            continue
        active, metadata = _mc_nonnegative_rate_table(
            rate_rows, process_type=reaction.process_type
        )
        energy, rate = _lookup(active, "mean_energy_eV", "rate_coefficient_m3_s")
        x = np.log(np.asarray(energy, dtype=float))
        scale = float(metadata["normalization_rate_m3_s"])
        y = np.asarray(rate, dtype=float) / scale
        interpolator = PchipInterpolator(x, y, extrapolate=False)
        anchor_error = float(np.max(np.abs(interpolator(x) - y)))
        minimum = math.inf
        maximum_interval_overshoot = 0.0
        for index in range(len(x) - 1):
            values = np.asarray(
                interpolator(np.linspace(x[index], x[index + 1], 257)),
                dtype=float,
            )
            minimum = min(minimum, float(np.min(values)))
            low = min(y[index], y[index + 1])
            high = max(y[index], y[index + 1])
            maximum_interval_overshoot = max(
                maximum_interval_overshoot,
                float(np.max(np.maximum(low - values, values - high))),
            )
        derivative = interpolator.derivative()
        maximum_derivative_jump = 0.0
        for index in range(1, len(x) - 1):
            delta = 1.0e-7 * min(
                x[index] - x[index - 1],
                x[index + 1] - x[index],
            )
            left = float(derivative(x[index] - delta))
            right = float(derivative(x[index] + delta))
            derivative_scale = max(
                1.0,
                abs(left),
                abs(right),
                abs(float(derivative(x[index]))),
            )
            maximum_derivative_jump = max(
                maximum_derivative_jump,
                abs(left - right) / derivative_scale,
            )
        zero_indices = np.flatnonzero(y == 0.0)
        zero_anchor_error = (
            float(np.max(np.abs(interpolator(x[zero_indices]))))
            if zero_indices.size
            else 0.0
        )
        process_passed = bool(
            np.isfinite(minimum)
            and minimum >= -1.0e-14
            and anchor_error <= 1.0e-12
            and zero_anchor_error <= 1.0e-15
            and maximum_interval_overshoot <= 1.0e-12
            and maximum_derivative_jump <= 1.0e-4
        )
        passed = passed and process_passed
        processes[reaction.process_type] = {
            "passed": process_passed,
            "anchors": len(x),
            "zero_anchors": int(zero_indices.size),
            "minimum_normalized_rate": minimum,
            "maximum_anchor_absolute_error": anchor_error,
            "maximum_zero_anchor_absolute_error": zero_anchor_error,
            "maximum_interval_shape_overshoot": (maximum_interval_overshoot),
            "maximum_normalized_one_sided_derivative_jump": (maximum_derivative_jump),
        }
    result = {
        "status": "passed" if passed else "failed",
        "passed": passed,
        "reference": "scipy_PchipInterpolator",
        "comsol_interpolator": "piecewisecubic",
        "abscissa": "log_mean_energy",
        "ordinate": "raw_rate_over_process_max",
        "contract": "nonnegative_shape_preserving_Hermite_C1",
        "processes": processes,
    }
    if not passed:
        raise GecCcpWorkflowError(
            "Monte Carlo zero-preserving reaction-rate interpolation audit failed"
        )
    return result
