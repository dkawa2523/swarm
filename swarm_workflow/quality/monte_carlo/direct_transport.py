"""Quality evaluation for fixed-population direct Monte Carlo transport."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from typing import Any

import electron_swarm.solvers.monte_carlo.evidence as _mc_evidence

from .contracts import (
    MC_TRANSPORT_COMPONENT_ESTIMATORS,
    MC_TRANSPORT_DIFFUSION_LAG_RELATIVE_TOLERANCE,
    MC_TRANSPORT_LAG_FIELDS,
    MC_TRANSPORT_MIN_ENSEMBLE_REPLICATES,
    MC_TRANSPORT_PRODUCTION_FIELDS,
    MC_TRANSPORT_SAMPLING_FIELDS,
    MC_TRANSPORT_STATE_RELATIVE_TOLERANCE,
    MC_TRANSPORT_STATE_STATIONARITY_FIELDS,
    _normalized_mc_run_settings,
    _paired_log_ratio_absolute_mean,
    _worst_stationarity_evidence,
    paired_log_ratio_ci95_bound,
)


@dataclass(slots=True)
class _DirectTransportEvidence:
    """Normalized replica evidence shared by all direct-MC quality stages."""

    reasons: list[str] = field(default_factory=list)
    stationarity_reference: dict[str, list[float]] = field(
        default_factory=lambda: {
            name: [] for name in MC_TRANSPORT_STATE_STATIONARITY_FIELDS
        }
    )
    stationarity_candidate: dict[str, list[float]] = field(
        default_factory=lambda: {
            name: [] for name in MC_TRANSPORT_STATE_STATIONARITY_FIELDS
        }
    )
    lag_reference: dict[str, list[float]] = field(
        default_factory=lambda: {name: [] for name in MC_TRANSPORT_LAG_FIELDS}
    )
    lag_candidate: dict[str, list[float]] = field(
        default_factory=lambda: {name: [] for name in MC_TRANSPORT_LAG_FIELDS}
    )
    case_mean_energy: list[float] = field(default_factory=list)
    transport_mean_energy: list[float] = field(default_factory=list)
    seeds: list[int] = field(default_factory=list)
    common_lag_s: tuple[float, ...] | None = None
    common_run_settings: str | None = None
    common_sampling_contract: str | None = None
    common_component_contract: str | None = None
    complete: bool = True

    def reject(self, reason: str) -> None:
        if reason not in self.reasons:
            self.reasons.append(reason)

    def incomplete(self, reason: str) -> None:
        self.reject(reason)
        self.complete = False


@dataclass(frozen=True, slots=True)
class _DirectPrecisionEvidence:
    stationarity_bounds: dict[str, float]
    lag_bounds: dict[str, float]
    mean_energy_bound: float | None
    stationarity_worst_field: str | None
    stationarity_max: float | None
    stationarity_tolerance: float | None
    lag_max: float | None


def _collect_direct_run_provenance(
    item: dict[str, Any], evidence: _DirectTransportEvidence
) -> None:
    provenance = item.get("mc_run_provenance")
    try:
        normalized_run_settings = (
            _normalized_mc_run_settings(provenance)
            if isinstance(provenance, dict)
            else None
        )
    except (KeyError, TypeError, ValueError, OverflowError):
        normalized_run_settings = None
    if normalized_run_settings is None:
        evidence.incomplete("mc_transport_run_provenance_invalid")
        return

    seed = provenance.get("seed")
    integers_valid = all(
        isinstance(provenance.get(name), int)
        and not isinstance(provenance.get(name), bool)
        for name in (
            "particles",
            "warmup_collisions",
            "production_collisions",
            "tail_max_collisions",
            "tail_collisions_executed",
        )
    )
    numeric_names = (
        "gas_number_density_m3",
        "electric_field_V_m",
        "trial_collision_frequency_s_inv",
        "tail_rate_rse_trigger",
        "magnetic_B_T",
        "magnetic_angle_EB_deg",
        "max_energy_limit_eV",
    )
    try:
        numeric_values = [float(provenance[name]) for name in numeric_names]
    except (TypeError, ValueError, OverflowError):
        numeric_values = [math.nan]
    strings_valid = all(
        isinstance(provenance.get(name), str) and provenance.get(name)
        for name in (
            "population_model",
            "collision_clock",
            "angular_scattering_model",
            "ionization_source_model",
            "high_energy_extrapolation",
        )
    )
    provenance_valid = bool(
        isinstance(seed, int)
        and not isinstance(seed, bool)
        and integers_valid
        and int(provenance["particles"]) >= 2
        and int(provenance["warmup_collisions"]) >= 0
        and int(provenance["production_collisions"]) > 0
        and int(provenance["tail_max_collisions"]) >= 0
        and int(provenance["tail_collisions_executed"])
        in {0, int(provenance["tail_max_collisions"])}
        and (
            provenance["population_model"] == "weighted_branching"
            or int(provenance["tail_collisions_executed"]) == 0
        )
        and strings_valid
        and isinstance(provenance.get("magnetic_enabled"), bool)
        and all(math.isfinite(value) for value in numeric_values)
        and float(provenance["gas_number_density_m3"]) > 0.0
        and float(provenance["electric_field_V_m"]) > 0.0
        and float(provenance["trial_collision_frequency_s_inv"]) > 0.0
        and 0.0 < float(provenance["tail_rate_rse_trigger"]) <= 1.0
        and float(provenance["magnetic_B_T"]) >= 0.0
        and float(provenance["max_energy_limit_eV"]) > 0.0
    )
    if not provenance_valid:
        evidence.incomplete("mc_transport_run_provenance_invalid")
        return

    evidence.seeds.append(int(seed))
    if evidence.common_run_settings is None:
        evidence.common_run_settings = normalized_run_settings
    elif normalized_run_settings != evidence.common_run_settings:
        evidence.incomplete("mc_transport_replica_settings_mismatch")


def _validate_direct_replica_header(
    item: dict[str, Any], evidence: _DirectTransportEvidence
) -> None:
    if (
        item.get("estimator") != "direct_mc_fixed_lag_block_helfand_energy_flux"
        or item.get("estimator_schema_version")
        != _mc_evidence.DIRECT_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
    ):
        evidence.incomplete("mc_transport_estimator_provenance_mismatch")
    if (
        item.get("energy_transport_status")
        != "direct_mc_candidate_for_ensemble_qualification"
        or item.get("finite") is not True
        or item.get("strictly_positive") is not True
        or item.get("multiple_time_origins") is not True
    ):
        evidence.incomplete("mc_direct_transport_basic_unqualified")
    components = item.get("component_estimators")
    if components != MC_TRANSPORT_COMPONENT_ESTIMATORS:
        evidence.incomplete("mc_transport_component_estimator_contract_mismatch")
    else:
        component_contract = json.dumps(
            components, sort_keys=True, separators=(",", ":")
        )
        if evidence.common_component_contract is None:
            evidence.common_component_contract = component_contract
        elif component_contract != evidence.common_component_contract:
            evidence.incomplete("mc_transport_component_estimator_contract_mismatch")
    _collect_direct_run_provenance(item, evidence)


def _validate_direct_sampling(
    sampling: dict[str, Any],
    lag_times: tuple[float, ...],
    completed_origins: list[Any],
) -> str:
    cadence_s = float(sampling["origin_cadence_s"])
    requested_planes = int(sampling["observation_planes_requested"])
    complete_planes = int(sampling["observation_planes_complete"])
    minimum_planes = int(sampling["minimum_observation_planes"])
    if (
        sampling["lag_grid_source"]
        != "fixed_physical_time_independent_of_production_length"
        or not math.isfinite(cadence_s)
        or cadence_s <= 0.0
        or sampling["origin_cadence_trial_periods"]
        != _mc_evidence.DIRECT_MC_TRANSPORT_CADENCE_TRIAL_PERIODS
        or sampling["rolling_origins"] is not True
        or sampling["production_length_controls_lag"] is not False
        or sampling["warmup_accumulators_reset"] is not True
        or minimum_planes != _mc_evidence.DIRECT_MC_TRANSPORT_OBSERVATION_PLANES
        or requested_planes < minimum_planes
        or complete_planes != requested_planes
        or any(
            not math.isclose(
                value,
                plane * cadence_s,
                rel_tol=1.0e-12,
                abs_tol=0.0,
            )
            for value, plane in zip(
                lag_times,
                _mc_evidence.DIRECT_MC_TRANSPORT_LAG_PLANES,
                strict=True,
            )
        )
    ):
        raise ValueError("invalid fixed-lag sampling contract")
    origin_counts = [int(value) for value in completed_origins]
    if (
        any(
            isinstance(value, bool) or not isinstance(value, int)
            for value in completed_origins
        )
        or any(value < 2 for value in origin_counts)
        or origin_counts
        != [
            requested_planes - lag
            for lag in _mc_evidence.DIRECT_MC_TRANSPORT_LAG_PLANES
        ]
    ):
        raise ValueError("invalid rolling-origin counts")
    return json.dumps(sampling, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _collect_direct_transport_values(
    item: dict[str, Any], evidence: _DirectTransportEvidence
) -> None:
    early = item.get("origin_stationarity_window_mean_early")
    late = item.get("origin_stationarity_window_mean_late")
    production = item.get("production_estimates")
    lag_scan = item.get("lag_scan")
    sampling = item.get("transport_sampling")
    if (
        not isinstance(early, dict)
        or set(early) != set(MC_TRANSPORT_STATE_STATIONARITY_FIELDS)
        or not isinstance(late, dict)
        or set(late) != set(MC_TRANSPORT_STATE_STATIONARITY_FIELDS)
        or not isinstance(production, dict)
        or set(production) != set(MC_TRANSPORT_PRODUCTION_FIELDS)
    ):
        evidence.incomplete("mc_transport_ensemble_evidence_missing")
        return
    if (
        not isinstance(lag_scan, dict)
        or set(lag_scan)
        != {
            "lag_planes",
            "lag_s",
            "estimates",
            "completed_origins",
            "selected_lag_index",
        }
        or not isinstance(sampling, dict)
        or set(sampling) != MC_TRANSPORT_SAMPLING_FIELDS
    ):
        evidence.incomplete("mc_transport_ensemble_evidence_missing")
        return
    lag_planes = lag_scan.get("lag_planes")
    lag_s = lag_scan.get("lag_s")
    estimates = lag_scan.get("estimates")
    completed_origins = lag_scan.get("completed_origins")
    selected_lag_index = lag_scan.get("selected_lag_index")
    if (
        lag_planes != list(_mc_evidence.DIRECT_MC_TRANSPORT_LAG_PLANES)
        or not isinstance(lag_s, list)
        or len(lag_s) != len(_mc_evidence.DIRECT_MC_TRANSPORT_LAG_PLANES)
        or not isinstance(estimates, dict)
        or set(estimates) != set(MC_TRANSPORT_LAG_FIELDS)
        or not isinstance(completed_origins, list)
        or len(completed_origins) != len(_mc_evidence.DIRECT_MC_TRANSPORT_LAG_PLANES)
        or selected_lag_index != len(_mc_evidence.DIRECT_MC_TRANSPORT_LAG_PLANES) - 1
    ):
        evidence.incomplete("mc_transport_lag_grid_mismatch")
        return
    try:
        lag_times = tuple(float(value) for value in lag_s)
        if any(not math.isfinite(value) or value <= 0.0 for value in lag_times) or any(
            right <= left for left, right in zip(lag_times, lag_times[1:])
        ):
            raise ValueError("invalid lag times")
        sampling_contract = _validate_direct_sampling(
            sampling, lag_times, completed_origins
        )
        if evidence.common_sampling_contract is None:
            evidence.common_sampling_contract = sampling_contract
        elif sampling_contract != evidence.common_sampling_contract:
            evidence.incomplete("mc_transport_replica_sampling_mismatch")
        if evidence.common_lag_s is None:
            evidence.common_lag_s = lag_times
        elif lag_times != evidence.common_lag_s:
            evidence.incomplete("mc_transport_replica_lag_times_mismatch")

        estimate_series: dict[str, list[float]] = {}
        for name in MC_TRANSPORT_LAG_FIELDS:
            raw_values = estimates[name]
            if not isinstance(raw_values, list) or len(raw_values) != len(lag_planes):
                raise ValueError("lag scan length mismatch")
            values = [float(value) for value in raw_values]
            if any(not math.isfinite(value) or value <= 0.0 for value in values):
                raise ValueError("nonpositive lag estimate")
            estimate_series[name] = values

        production_values = {
            name: float(production[name]) for name in MC_TRANSPORT_PRODUCTION_FIELDS
        }
        if any(not math.isfinite(value) for value in production_values.values()) or any(
            production_values[name] <= 0.0
            for name in (
                "mobility_m2_V_s",
                "energy_mobility_m2_V_s",
                "mean_energy_eV",
            )
        ):
            raise ValueError("invalid production transport estimate")
        for name in MC_TRANSPORT_STATE_STATIONARITY_FIELDS:
            early_value = float(early[name])
            late_value = float(late[name])
            if (
                not math.isfinite(early_value)
                or not math.isfinite(late_value)
                or early_value <= 0.0
                or late_value <= 0.0
            ):
                raise ValueError("nonpositive stationarity evidence")
            evidence.stationarity_reference[name].append(early_value)
            evidence.stationarity_candidate[name].append(late_value)
        for name in MC_TRANSPORT_LAG_FIELDS:
            values = estimate_series[name]
            evidence.lag_reference[name].append(float(values[-2]))
            evidence.lag_candidate[name].append(float(values[-1]))

        _validate_direct_reported_transport(
            item, production_values, estimate_series, selected_lag_index, evidence
        )
        case_mean = float(item["_case_mean_energy_eV"])
        selected_mean = production_values["mean_energy_eV"]
        if not math.isfinite(case_mean) or case_mean <= 0.0:
            raise ValueError("case mean energy is not positive")
        evidence.case_mean_energy.append(case_mean)
        evidence.transport_mean_energy.append(selected_mean)
    except (KeyError, TypeError, ValueError, OverflowError):
        evidence.incomplete("mc_transport_positive_ensemble_evidence_required")


def _validate_direct_reported_transport(
    item: dict[str, Any],
    production_values: dict[str, float],
    estimate_series: dict[str, list[float]],
    selected_lag_index: Any,
    evidence: _DirectTransportEvidence,
) -> None:
    reported = item.get("_case_reported_transport")
    if not isinstance(reported, dict):
        raise ValueError("case transport values missing")
    for name in (
        "drift_velocity_m_s",
        "mobility_m2_V_s",
        "energy_mobility_m2_V_s",
    ):
        actual = float(reported[name])
        expected = production_values[name]
        if not math.isfinite(actual) or not math.isclose(
            actual, expected, rel_tol=1.0e-12, abs_tol=0.0
        ):
            evidence.incomplete("mc_transport_report_not_production_estimate")
    for name in (
        "diffusion_L_m2_s",
        "diffusion_T_m2_s",
        "energy_diffusion_L_m2_s",
        "energy_diffusion_T_m2_s",
    ):
        actual = float(reported[name])
        expected = estimate_series[name][selected_lag_index]
        if not math.isfinite(actual) or not math.isclose(
            actual, expected, rel_tol=1.0e-12, abs_tol=0.0
        ):
            evidence.incomplete("mc_transport_report_not_selected_lag")


def _collect_direct_evidence(
    diagnostics: list[dict[str, Any]],
) -> _DirectTransportEvidence:
    evidence = _DirectTransportEvidence()
    if len(diagnostics) < MC_TRANSPORT_MIN_ENSEMBLE_REPLICATES:
        evidence.reject("mc_transport_ensemble_replicates_insufficient")
    for item in diagnostics:
        if not isinstance(item, dict) or item.get("_diagnostic_invalid"):
            evidence.incomplete("mc_transport_replica_diagnostic_invalid")
            continue
        _validate_direct_replica_header(item, evidence)
        _collect_direct_transport_values(item, evidence)
    if len(evidence.seeds) != len(diagnostics) or len(set(evidence.seeds)) != len(
        diagnostics
    ):
        evidence.incomplete("mc_transport_replica_seeds_missing_or_not_distinct")
    return evidence


def _evaluate_direct_precision(
    evidence: _DirectTransportEvidence,
) -> _DirectPrecisionEvidence:
    stationarity_bounds: dict[str, float] = {}
    lag_bounds: dict[str, float] = {}
    mean_energy_bound: float | None = None
    if evidence.complete:
        for name in MC_TRANSPORT_STATE_STATIONARITY_FIELDS:
            bound = paired_log_ratio_ci95_bound(
                evidence.stationarity_reference[name],
                evidence.stationarity_candidate[name],
            )
            if bound is None:
                evidence.complete = False
                break
            stationarity_bounds[name] = bound
    if evidence.complete:
        for name in MC_TRANSPORT_LAG_FIELDS:
            bound = paired_log_ratio_ci95_bound(
                evidence.lag_reference[name], evidence.lag_candidate[name]
            )
            if bound is None:
                evidence.complete = False
                break
            lag_bounds[name] = bound
    if evidence.complete:
        mean_energy_bound = paired_log_ratio_ci95_bound(
            evidence.transport_mean_energy,
            evidence.case_mean_energy,
        )
        if mean_energy_bound is None:
            evidence.complete = False
    if not evidence.complete:
        evidence.reject("mc_transport_ensemble_evidence_missing")
    worst_field, stationarity_max, tolerance = _worst_stationarity_evidence(
        stationarity_bounds
    )
    return _DirectPrecisionEvidence(
        stationarity_bounds=stationarity_bounds,
        lag_bounds=lag_bounds,
        mean_energy_bound=mean_energy_bound,
        stationarity_worst_field=worst_field,
        stationarity_max=stationarity_max,
        stationarity_tolerance=tolerance,
        lag_max=max(lag_bounds.values(), default=None),
    )


def _apply_direct_acceptance(
    evidence: _DirectTransportEvidence, precision: _DirectPrecisionEvidence
) -> None:
    if (
        precision.stationarity_max is not None
        and precision.stationarity_tolerance is not None
        and precision.stationarity_max > precision.stationarity_tolerance
    ):
        evidence.reject(
            f"mc_origin_stationarity_not_converged:{precision.stationarity_worst_field}"
        )
    if precision.mean_energy_bound is not None and (
        precision.mean_energy_bound > MC_TRANSPORT_STATE_RELATIVE_TOLERANCE
    ):
        evidence.reject("mc_transport_case_mean_energy_inconsistent")
    if precision.lag_max is not None and (
        precision.lag_max > MC_TRANSPORT_DIFFUSION_LAG_RELATIVE_TOLERANCE
    ):
        worst = max(precision.lag_bounds, key=precision.lag_bounds.__getitem__)
        evidence.reject(f"mc_transport_lag_plateau_not_converged:{worst}")


def _direct_summary(
    diagnostics: list[dict[str, Any]],
    evidence: _DirectTransportEvidence,
    precision: _DirectPrecisionEvidence,
) -> dict[str, Any]:
    qualified = not evidence.reasons
    return {
        "solver_diagnostics_available": 1,
        "solver_converged": None,
        "solver_transport_qualified": int(qualified),
        "solver_iterations_max": None,
        "solver_residual_L1_max": None,
        "solver_residual_tolerance": None,
        "solver_tail_probability_max": None,
        "solver_tail_probability_target": None,
        "solver_edge_to_peak_max": None,
        "solver_edge_to_peak_target": None,
        "solver_grid_max_eV_max": None,
        "solver_grid_max_limit_eV": None,
        "solver_grid_limit_hit": None,
        "solver_diagnostics_passed": int(qualified),
        "solver_transport_replicates": len(diagnostics),
        "solver_origin_stationarity_limiting_field": (
            precision.stationarity_worst_field
        ),
        "solver_origin_stationarity_limiting_relative_ci95_bound": (
            precision.stationarity_max
        ),
        "solver_mean_energy_stationarity_relative_ci95_bound": (
            precision.stationarity_bounds.get("mean_energy_eV")
        ),
        "solver_mobility_stationarity_absolute_log_drift": (
            _paired_log_ratio_absolute_mean(
                evidence.stationarity_reference.get("mobility_m2_V_s", []),
                evidence.stationarity_candidate.get("mobility_m2_V_s", []),
            )
        ),
        "solver_lag_convergence_max_relative_ci95_bound": precision.lag_max,
        "solver_lag_supplemental_replicates": None,
        "solver_lag_supplemental_max_relative_ci95_bound": None,
        "solver_transport_mean_energy_max_relative_ci95_bound": (
            precision.mean_energy_bound
        ),
        "solver_population_growth_gate_mode": None,
        "solver_population_growth_max_relative_ci95_bound": None,
        "solver_population_growth_poisson_interval_max_ratio": None,
        "solver_population_growth_sparse_max_metric": None,
        "solver_population_growth_pooled_event_count": None,
        "solver_population_growth_pooled_exposure_s": None,
        "solver_population_growth_pooled_frequency_ci95_low_s_inv": None,
        "solver_population_growth_pooled_frequency_ci95_high_s_inv": None,
        "solver_population_growth_pooled_direct_within_ci95": None,
        "solver_population_growth_pooled_population_within_ci95": None,
        "solver_origin_stationarity_limiting_relative_tolerance": (
            precision.stationarity_tolerance
        ),
        "solver_lag_convergence_relative_tolerance": (
            MC_TRANSPORT_DIFFUSION_LAG_RELATIVE_TOLERANCE
        ),
        "solver_transport_mean_energy_relative_tolerance": (
            MC_TRANSPORT_STATE_RELATIVE_TOLERANCE
        ),
        "solver_population_growth_relative_tolerance": None,
        "solver_failure_reasons": evidence.reasons,
    }


def summarize_mc_transport_diagnostics(
    diagnostics: list[dict[str, Any]],
) -> dict[str, Any]:
    """Qualify direct MC transport using an independent-replica ensemble.

    Every raw case remains represented in ``diagnostics``. Missing, malformed,
    mixed-version, mixed-setting, incomplete-plane, and duplicate-seed evidence
    fails closed. Overlapping origins stay within a replica; convergence uses
    only paired log ratios across independent replicas.
    """

    evidence = _collect_direct_evidence(diagnostics)
    precision = _evaluate_direct_precision(evidence)
    _apply_direct_acceptance(evidence, precision)
    return _direct_summary(diagnostics, evidence, precision)
