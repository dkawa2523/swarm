"""Pure Monte Carlo transport and closure-quality evaluation.

This module owns numerical qualification of raw Monte Carlo diagnostics.  It
does not read workflow databases or write table artifacts; callers supply the
diagnostic and aggregate-quality records explicitly.
"""

from __future__ import annotations

import json
import math
from typing import TYPE_CHECKING, Any

import numpy as np
from scipy.stats import chi2, t as student_t

import electron_swarm.solvers.monte_carlo.evidence as _mc_evidence

from ...campaign.statistics import summarize_replicates

if TYPE_CHECKING:
    from ..policy import QualityThresholds


MC_TRANSPORT_STATE_RELATIVE_TOLERANCE = 0.1
MC_TRANSPORT_DIFFUSION_LAG_RELATIVE_TOLERANCE = 0.25
WEIGHTED_MC_TRANSPORT_DIFFUSION_IDENTIFICATION_RELATIVE_TOLERANCE = 0.15
MC_EEDF_RATE_ENSEMBLE_RELATIVE_TOLERANCE = 0.1
MC_POPULATION_GROWTH_RELATIVE_TOLERANCE = 0.1
MC_TRANSPORT_MIN_ENSEMBLE_REPLICATES = 3
MC_TRANSPORT_STATE_STATIONARITY_FIELDS = ("mean_energy_eV",)
MC_TRANSPORT_LAG_FIELDS = (
    "diffusion_L_m2_s",
    "diffusion_T_m2_s",
    "energy_diffusion_L_m2_s",
    "energy_diffusion_T_m2_s",
)
MC_TRANSPORT_PRODUCTION_FIELDS = (
    "drift_velocity_m_s",
    "mobility_m2_V_s",
    "energy_mobility_m2_V_s",
    "mean_energy_eV",
)
MC_TRANSPORT_COMPONENT_ESTIMATORS = {
    "drift_mobility": "production_trajectory_displacement_over_residence_time",
    "particle_diffusion": "block_helfand_mean_square_displacement",
    "energy_mobility": "production_residence_energy_current",
    "energy_diffusion": "restricted_density_packet_energy_flux_fixed_lag",
}
WEIGHTED_MC_TRANSPORT_STATIONARITY_FIELDS = (
    "mobility_m2_V_s",
    "diffusion_L_m2_s",
    "diffusion_T_m2_s",
    "energy_mobility_m2_V_s",
    "energy_diffusion_L_m2_s",
    "energy_diffusion_T_m2_s",
    "mean_energy_eV",
)
WEIGHTED_MC_TRANSPORT_COMPONENT_ESTIMATORS = {
    "drift_mobility": "normalized_weight_residence_velocity_moment",
    "particle_diffusion": ("synchronized_block_lag_position_velocity_flux_covariance"),
    "energy_mobility": "normalized_weight_residence_energy_current",
    "energy_diffusion": (
        "synchronized_block_lag_restricted_density_packet_energy_current_covariance"
    ),
}
MC_TRANSPORT_SAMPLING_FIELDS = {
    "lag_grid_source",
    "origin_cadence_s",
    "origin_cadence_trial_periods",
    "rolling_origins",
    "production_length_controls_lag",
    "warmup_accumulators_reset",
    "observation_planes_requested",
    "observation_planes_complete",
    "minimum_observation_planes",
}
MC_TRANSPORT_RUN_PROVENANCE_FIELDS = (
    "seed",
    "particles",
    "warmup_collisions",
    "production_collisions",
    "tail_max_collisions",
    "tail_rate_rse_trigger",
    "tail_collisions_executed",
    "population_model",
    "gas_number_density_m3",
    "electric_field_V_m",
    "trial_collision_frequency_s_inv",
    "collision_clock",
    "magnetic_enabled",
    "magnetic_B_T",
    "magnetic_angle_EB_deg",
    "angular_scattering_model",
    "ionization_source_model",
    "high_energy_extrapolation",
    "max_energy_limit_eV",
)


def _normalized_mc_run_settings(provenance: dict[str, Any]) -> str:
    """Normalize qualification-critical MC settings except the replica seed."""

    missing = set(MC_TRANSPORT_RUN_PROVENANCE_FIELDS) - set(provenance)
    if missing:
        raise ValueError(
            f"missing Monte Carlo run provenance fields: {sorted(missing)}"
        )
    normalized = {
        key: provenance[key]
        for key in MC_TRANSPORT_RUN_PROVENANCE_FIELDS
        if key not in {"seed", "tail_collisions_executed"}
    }
    return json.dumps(
        normalized,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def _worst_stationarity_evidence(
    bounds: dict[str, float],
) -> tuple[str | None, float | None, float | None]:
    """Return the field with the largest fraction of its stationarity limit."""

    if not bounds:
        return None, None, None

    def tolerance(name: str) -> float:
        return (
            MC_TRANSPORT_DIFFUSION_LAG_RELATIVE_TOLERANCE
            if "diffusion" in name
            else MC_TRANSPORT_STATE_RELATIVE_TOLERANCE
        )

    worst = max(bounds, key=lambda name: bounds[name] / tolerance(name))
    return worst, float(bounds[worst]), tolerance(worst)


def quality_failure_reasons(row: dict[str, Any], name: str) -> list[str]:
    try:
        values = json.loads(str(row.get(name, "[]")))
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid quality reasons in {name}") from exc
    if not isinstance(values, list) or any(
        not isinstance(value, str) or not value for value in values
    ):
        raise ValueError(f"invalid quality reasons in {name}")
    return values


def _quality_number(row: dict[str, Any], name: str) -> float:
    value = float(row[name])
    if not math.isfinite(value):
        raise ValueError(f"quality field {name} is not finite")
    return value


def _quality_integer(row: dict[str, Any], name: str) -> int:
    value = row[name]
    if isinstance(value, bool):
        raise ValueError(f"quality field {name} is not an integer")
    number = int(str(value))
    if float(value) != float(number):
        raise ValueError(f"quality field {name} is not an integer")
    return number


def _threshold_number(
    thresholds: QualityThresholds | dict[str, Any], name: str
) -> float:
    raw = (
        thresholds[name] if isinstance(thresholds, dict) else getattr(thresholds, name)
    )
    value = float(raw)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError(f"quality threshold {name} is invalid")
    return value


def _restricted_lmea_inactive_failure(reason: str) -> bool:
    """Identify diagnostics for MC coefficients not consumed by this closure."""

    inactive_aggregate_prefixes = (
        "reduced_diffusion_L_",
        "reduced_diffusion_T_",
        "reduced_energy_mobility_",
        "reduced_energy_diffusion_",
        "reduced_electron_energy_mobility_",
        "reduced_electron_energy_diffusion_",
    )
    if reason.startswith(inactive_aggregate_prefixes):
        return True
    if reason.startswith(
        (
            "mc_transport_lag_plateau_not_converged:",
            "mc_weighted_growth_lag_not_converged:",
        )
    ):
        return True
    if reason in {
        "mc_transport_report_not_selected_lag",
        "mc_weighted_growth_supplemental_moment_evidence_missing",
    }:
        return True
    stationarity_prefix = "mc_weighted_growth_stationarity_not_converged:"
    if reason.startswith(stationarity_prefix):
        field = reason.removeprefix(stationarity_prefix)
        return field.startswith(
            (
                "mobility_",
                "diffusion_",
                "energy_mobility_",
                "energy_diffusion_",
            )
        )
    return False


def mc_function_eedf_restricted_lmea_failure_reasons(
    row: dict[str, Any],
    *,
    thresholds: QualityThresholds | dict[str, Any],
    estimator_schema: str | None,
) -> list[str]:
    """Recompute failures for the MC quantities active in restricted LMEA.

    Full particle/energy diffusion diagnostics remain in ``quality.csv`` but
    do not disqualify a closure that imports only MC mobility, EEDF, and the
    separately audited elastic-loss moment. Weighted-population lineage,
    growth, EEDF/rate-tail, mobility, and mean-energy evidence remain gates.
    """

    failures = [
        reason
        for reason in quality_failure_reasons(row, "failure_reasons_json")
        if not _restricted_lmea_inactive_failure(reason)
    ]

    def reject(reason: str) -> None:
        if reason not in failures:
            failures.append(reason)

    def number(name: str) -> float | None:
        try:
            return _quality_number(row, name)
        except (KeyError, TypeError, ValueError, OverflowError):
            reject(name)
            return None

    def integer(name: str) -> int | None:
        try:
            return _quality_integer(row, name)
        except (KeyError, TypeError, ValueError, OverflowError):
            reject(name)
            return None

    if estimator_schema != _mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION:
        reject("mc_restricted_lmea_requires_current_weighted_growth")
    valid_replicates = integer("valid_replicates")
    if (
        valid_replicates is not None
        and valid_replicates < MC_TRANSPORT_MIN_ENSEMBLE_REPLICATES
    ):
        reject("valid_replicates")
    transport_replicates = integer("solver_transport_replicates")
    if (
        transport_replicates is not None
        and transport_replicates < MC_TRANSPORT_MIN_ENSEMBLE_REPLICATES
    ):
        reject("solver_transport_replicates")
    diagnostics_available = integer("solver_diagnostics_available")
    if diagnostics_available is not None and diagnostics_available != 1:
        reject("solver_diagnostics_available")

    mobility_rse = number("mobility_rse")
    if mobility_rse is not None and (
        mobility_rse < 0.0
        or mobility_rse > _threshold_number(thresholds, "mobility_rse")
    ):
        reject("mobility_rse")
    if mobility_rse is not None and valid_replicates is not None:
        mobility_precision = float(
            student_t.ppf(0.975, valid_replicates - 1) * mobility_rse
        )
        if (
            not math.isfinite(mobility_precision)
            or mobility_precision > MC_TRANSPORT_STATE_RELATIVE_TOLERANCE
        ):
            reject("mobility_relative_ci95_half_width")
    mobility_drift = number("solver_mobility_stationarity_absolute_log_drift")
    if mobility_drift is not None and (
        mobility_drift < 0.0 or mobility_drift > MC_TRANSPORT_STATE_RELATIVE_TOLERANCE
    ):
        reject("solver_mobility_stationarity_absolute_log_drift")
    mean_energy_stationarity = number(
        "solver_mean_energy_stationarity_relative_ci95_bound"
    )
    if mean_energy_stationarity is not None and (
        mean_energy_stationarity < 0.0
        or mean_energy_stationarity > MC_TRANSPORT_STATE_RELATIVE_TOLERANCE
    ):
        reject("solver_mean_energy_stationarity_relative_ci95_bound")
    normalization_error = number("eedf_normalization_error")
    if normalization_error is not None and (
        normalization_error < 0.0
        or normalization_error
        > _threshold_number(thresholds, "eedf_normalization_error")
    ):
        reject("eedf_normalization_error")
    major_rate_text = str(row.get("max_major_rate_rse", "")).strip()
    if major_rate_text:
        major_rate_rse = number("max_major_rate_rse")
        if major_rate_rse is not None and (
            major_rate_rse < 0.0
            or major_rate_rse > _threshold_number(thresholds, "major_rate_rse")
        ):
            reject("max_major_rate_rse")

    mean_energy_bound = number("solver_transport_mean_energy_max_relative_ci95_bound")
    if mean_energy_bound is not None and (
        mean_energy_bound < 0.0
        or mean_energy_bound > MC_TRANSPORT_STATE_RELATIVE_TOLERANCE
    ):
        reject("solver_transport_mean_energy_max_relative_ci95_bound")
    mean_energy_tolerance = number("solver_transport_mean_energy_relative_tolerance")
    if mean_energy_tolerance is not None and not math.isclose(
        mean_energy_tolerance,
        MC_TRANSPORT_STATE_RELATIVE_TOLERANCE,
        rel_tol=1.0e-12,
        abs_tol=0.0,
    ):
        reject("solver_transport_mean_energy_relative_tolerance")

    growth_tolerance = number("solver_population_growth_relative_tolerance")
    if growth_tolerance is not None and not math.isclose(
        growth_tolerance,
        MC_POPULATION_GROWTH_RELATIVE_TOLERANCE,
        rel_tol=1.0e-12,
        abs_tol=0.0,
    ):
        reject("solver_population_growth_relative_tolerance")
    growth_mode = str(row.get("solver_population_growth_gate_mode", "")).strip()
    if growth_mode == "positive_paired_log":
        bound = number("solver_population_growth_max_relative_ci95_bound")
        if bound is not None and (
            bound < 0.0 or bound > MC_POPULATION_GROWTH_RELATIVE_TOLERANCE
        ):
            reject("solver_population_growth_max_relative_ci95_bound")
    elif growth_mode == "sparse_negligible_transport":
        interval_ratio = number("solver_population_growth_poisson_interval_max_ratio")
        sparse_metric = number("solver_population_growth_sparse_max_metric")
        if interval_ratio is not None and (
            interval_ratio < 0.0 or interval_ratio > 1.0
        ):
            reject("solver_population_growth_poisson_interval_max_ratio")
        if sparse_metric is not None and (
            sparse_metric < 0.0
            or sparse_metric > MC_POPULATION_GROWTH_RELATIVE_TOLERANCE
        ):
            reject("solver_population_growth_sparse_max_metric")
    else:
        reject("solver_population_growth_gate_mode")
    if row.get("quality_source") != "monte_carlo_aggregate_quality":
        reject("quality_source")
    return sorted(set(failures))


def paired_log_ratio_ci95_bound(
    reference_values: list[float],
    candidate_values: list[float],
) -> float | None:
    """Return the Student-t 95% bound for paired replica log ratios.

    Replica pairs are independent; time origins within one trajectory are not.
    Positive candidate/reference pairs are converted to
    ``log(candidate/reference)`` before ensemble summarization.  This is
    dimensionless, symmetric for reciprocal changes, and has no random
    denominator shared with its uncertainty estimate.
    """

    if (
        len(reference_values) != len(candidate_values)
        or len(candidate_values) < MC_TRANSPORT_MIN_ENSEMBLE_REPLICATES
    ):
        return None
    reference = np.asarray(reference_values, dtype=float)
    candidate = np.asarray(candidate_values, dtype=float)
    if (
        np.any(~np.isfinite(reference))
        or np.any(~np.isfinite(candidate))
        or np.any(reference <= 0.0)
        or np.any(candidate <= 0.0)
    ):
        return None
    ratio_stats = summarize_replicates(np.log(candidate / reference))
    if (
        ratio_stats.mean is None
        or ratio_stats.standard_error is None
        or ratio_stats.ci95_critical_value is None
    ):
        return None
    return float(
        abs(ratio_stats.mean)
        + ratio_stats.ci95_critical_value * ratio_stats.standard_error
    )


def _paired_log_ratio_absolute_mean(
    reference_values: list[float],
    candidate_values: list[float],
) -> float | None:
    """Return systematic paired drift without folding noise into the drift."""

    if (
        len(reference_values) != len(candidate_values)
        or len(candidate_values) < MC_TRANSPORT_MIN_ENSEMBLE_REPLICATES
    ):
        return None
    reference = np.asarray(reference_values, dtype=float)
    candidate = np.asarray(candidate_values, dtype=float)
    if (
        np.any(~np.isfinite(reference))
        or np.any(~np.isfinite(candidate))
        or np.any(reference <= 0.0)
        or np.any(candidate <= 0.0)
    ):
        return None
    return float(abs(np.mean(np.log(candidate / reference))))


def _poisson_frequency_interval_95(
    event_count: int,
    exposure_s: float,
) -> tuple[float, float]:
    """Return the exact 95% Poisson frequency interval.

    A zero count has the useful one-sided 95% upper limit.  Positive counts
    use the central 95% Garwood interval rather than a central 90% interval.
    """

    count = int(event_count)
    exposure = float(exposure_s)
    if count < 0 or not math.isfinite(exposure) or exposure <= 0.0:
        raise ValueError("Poisson frequency evidence is invalid")
    if count == 0:
        lower = 0.0
        upper = -math.log(0.05) / exposure
    else:
        lower = 0.5 * float(chi2.ppf(0.025, 2 * count)) / exposure
        upper = 0.5 * float(chi2.ppf(0.975, 2 * (count + 1))) / exposure
    if not (math.isfinite(lower) and math.isfinite(upper) and 0.0 <= lower <= upper):
        raise ValueError("Poisson frequency interval is invalid")
    return lower, upper
