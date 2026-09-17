"""Independent quality and anchor-coverage audit for GEC MC bundles."""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from typing import Any

import electron_swarm.solvers.monte_carlo.evidence as _mc_evidence

from swarm_workflow.quality.monte_carlo.contracts import (
    MC_POPULATION_GROWTH_RELATIVE_TOLERANCE,
    MC_TRANSPORT_DIFFUSION_LAG_RELATIVE_TOLERANCE,
    MC_TRANSPORT_STATE_RELATIVE_TOLERANCE,
    MC_TRANSPORT_STATE_STATIONARITY_FIELDS,
    WEIGHTED_MC_TRANSPORT_DIFFUSION_IDENTIFICATION_RELATIVE_TOLERANCE,
    WEIGHTED_MC_TRANSPORT_STATIONARITY_FIELDS,
    mc_function_eedf_restricted_lmea_failure_reasons,
)
from swarm_workflow.tables.contracts import (
    MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA,
)

from ..contracts import GecCcpWorkflowError
from ..data import _float_or_none, _required_float


@dataclass(frozen=True)
class McQualityPolicy:
    """Frozen policy used to recompute each quality-row decision."""

    source: str
    thresholds: dict[str, Any]
    estimator_schema: str | None
    qualification_profile: str | None
    restricted_lmea: bool


def independent_bundle_quality_audit(
    rows: list[dict[str, Any]],
    *,
    source: str,
    thresholds: dict[str, Any],
    mc_sampling_plan: list[dict[str, Any]] | None,
    mc_estimator_schema: str | None,
    source_policy: dict[str, Any] | None = None,
    transport_rows: list[dict[str, Any]] | None = None,
    valid_ranges: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Recompute the COMSOL-entry quality decision from numeric evidence."""

    if not rows:
        raise GecCcpWorkflowError("bundle quality.csv contains no rows")
    qualification_profile = (
        source_policy.get("qualification_profile")
        if source == "monte_carlo" and source_policy is not None
        else None
    )
    policy = McQualityPolicy(
        source=source,
        thresholds=thresholds,
        estimator_schema=mc_estimator_schema,
        qualification_profile=qualification_profile,
        restricted_lmea=(
            qualification_profile == MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
        ),
    )
    row_results = [_quality_row_result(row, policy) for row in rows]
    quality_fields = [float(item["E_over_N_Td"]) for item in row_results]
    if len(quality_fields) != len(set(quality_fields)):
        raise GecCcpWorkflowError("quality.csv contains duplicate E/N anchors")
    anchor_coverage = (
        _mc_anchor_coverage(
            quality_fields,
            mc_sampling_plan=mc_sampling_plan,
            source_policy=source_policy,
            transport_rows=transport_rows,
            valid_ranges=valid_ranges,
            restricted_lmea=policy.restricted_lmea,
        )
        if source == "monte_carlo"
        else None
    )
    return {
        "passed": all(item["passed"] for item in row_results),
        "source": source,
        "decision": "independently_recomputed_from_quality_columns",
        "rows": row_results,
        "thresholds": thresholds,
        "transport_anchor_coverage": anchor_coverage,
    }


def _quality_row_result(
    row: dict[str, Any],
    policy: McQualityPolicy,
) -> dict[str, Any]:
    field = _required_float(row, "E_over_N_Td")
    normalization_error = _required_float(row, "eedf_normalization_error")
    reasons: list[str] = []
    if policy.restricted_lmea:
        try:
            recomputed = mc_function_eedf_restricted_lmea_failure_reasons(
                row,
                thresholds=policy.thresholds,
                estimator_schema=policy.estimator_schema,
            )
            reported_active_reasons = json.loads(
                str(row["active_closure_failure_reasons_json"])
            )
        except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
            raise GecCcpWorkflowError(
                "Monte Carlo restricted-LMEA quality evidence is invalid"
            ) from exc
        if row.get("qualification_profile") != policy.qualification_profile:
            reasons.append("qualification_profile")
        if not _true_value(row, "active_closure_quality_passed"):
            reasons.append("active_closure_quality_passed")
        if reported_active_reasons != recomputed:
            reasons.append("active_closure_failure_reasons_json")
        reasons.extend(recomputed)
    else:
        if not _true_value(row, "passed"):
            reasons.append("reported_passed_is_false")
        eedf_limit = float(policy.thresholds["eedf_normalization_error"])
        if normalization_error < 0.0 or normalization_error > eedf_limit:
            reasons.append("eedf_normalization_error")
    if policy.source == "monte_carlo" and not policy.restricted_lmea:
        _standard_mc_quality_reasons(row, policy, reasons)
    return {
        "E_over_N_Td": field,
        "passed": not reasons,
        "failure_reasons": reasons,
    }


def _standard_mc_quality_reasons(
    row: dict[str, Any],
    policy: McQualityPolicy,
    reasons: list[str],
) -> None:
    for name in ("failure_reasons_json", "aggregate_failure_reasons_json"):
        try:
            reported_reasons = json.loads(str(row[name]))
        except (KeyError, TypeError, ValueError) as exc:
            raise GecCcpWorkflowError(
                f"Monte Carlo quality column {name} is invalid"
            ) from exc
        if reported_reasons != []:
            reasons.append(name)
    for name in (
        "aggregate_quality_passed",
        "uncertainty_available",
        "solver_diagnostics_available",
        "solver_transport_qualified",
    ):
        if not _true_value(row, name):
            reasons.append(name)
    if _integer(row, "valid_replicates") < 2:
        reasons.append("valid_replicates")
    if _integer(row, "solver_transport_replicates") < 3:
        reasons.append("solver_transport_replicates")
    rse_limits = {
        "mobility_rse": float(policy.thresholds["mobility_rse"]),
        "energy_mobility_rse": float(policy.thresholds["mobility_rse"]),
        "diffusion_L_rse": float(policy.thresholds["diffusion_rse"]),
        "diffusion_T_rse": float(policy.thresholds["diffusion_rse"]),
        "energy_diffusion_L_rse": float(policy.thresholds["diffusion_rse"]),
        "energy_diffusion_T_rse": float(policy.thresholds["diffusion_rse"]),
        "max_major_rate_rse": float(policy.thresholds["major_rate_rse"]),
    }
    for name, limit in rse_limits.items():
        value = _required_float(row, name)
        if value < 0.0 or value > limit:
            reasons.append(name)
    _transport_diagnostic_reasons(row, policy.estimator_schema, reasons)
    if (
        policy.estimator_schema
        == _mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
    ):
        _population_growth_reasons(row, reasons)
    if row.get("quality_source") != "monte_carlo_aggregate_quality":
        reasons.append("quality_source")


def _transport_diagnostic_reasons(
    row: dict[str, Any],
    estimator_schema: str | None,
    reasons: list[str],
) -> None:
    bound_pairs = [
        (
            "solver_origin_stationarity_limiting_relative_ci95_bound",
            "solver_origin_stationarity_limiting_relative_tolerance",
        ),
        (
            "solver_transport_mean_energy_max_relative_ci95_bound",
            "solver_transport_mean_energy_relative_tolerance",
        ),
    ]
    if estimator_schema == _mc_evidence.DIRECT_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION:
        lag_tolerance = MC_TRANSPORT_DIFFUSION_LAG_RELATIVE_TOLERANCE
        allowed_stationarity_fields = set(MC_TRANSPORT_STATE_STATIONARITY_FIELDS)
        bound_pairs.insert(
            1,
            (
                "solver_lag_convergence_max_relative_ci95_bound",
                "solver_lag_convergence_relative_tolerance",
            ),
        )
        if any(
            str(row.get(name, "")).strip()
            for name in (
                "solver_population_growth_gate_mode",
                "solver_population_growth_max_relative_ci95_bound",
                "solver_population_growth_poisson_interval_max_ratio",
                "solver_population_growth_sparse_max_metric",
                "solver_population_growth_relative_tolerance",
            )
        ):
            reasons.append("unexpected_population_growth_evidence")
    elif (
        estimator_schema == _mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
    ):
        lag_tolerance = (
            WEIGHTED_MC_TRANSPORT_DIFFUSION_IDENTIFICATION_RELATIVE_TOLERANCE
        )
        allowed_stationarity_fields = set(WEIGHTED_MC_TRANSPORT_STATIONARITY_FIELDS)
        bound_pairs[1:1] = [
            (
                "solver_lag_convergence_max_relative_ci95_bound",
                "solver_lag_convergence_relative_tolerance",
            ),
        ]
    else:
        lag_tolerance = MC_TRANSPORT_DIFFUSION_LAG_RELATIVE_TOLERANCE
        allowed_stationarity_fields = set()
        reasons.append("mc_transport_estimator_schema")
    stationarity_field = str(
        row.get("solver_origin_stationarity_limiting_field", "")
    ).strip()
    if stationarity_field not in allowed_stationarity_fields:
        reasons.append("solver_origin_stationarity_limiting_field")
    stationarity_tolerance = (
        MC_TRANSPORT_DIFFUSION_LAG_RELATIVE_TOLERANCE
        if "diffusion" in stationarity_field
        else MC_TRANSPORT_STATE_RELATIVE_TOLERANCE
    )
    expected_tolerances = {
        "solver_origin_stationarity_limiting_relative_tolerance": (
            stationarity_tolerance
        ),
        "solver_lag_convergence_relative_tolerance": lag_tolerance,
        "solver_transport_mean_energy_relative_tolerance": (
            MC_TRANSPORT_STATE_RELATIVE_TOLERANCE
        ),
        "solver_population_growth_relative_tolerance": (
            MC_POPULATION_GROWTH_RELATIVE_TOLERANCE
        ),
    }
    for bound_name, tolerance_name in bound_pairs:
        bound = _required_float(row, bound_name)
        tolerance = _required_float(row, tolerance_name)
        expected_tolerance = expected_tolerances[tolerance_name]
        if (
            bound < 0.0
            or tolerance <= 0.0
            or not math.isclose(
                tolerance,
                expected_tolerance,
                rel_tol=1.0e-12,
                abs_tol=0.0,
            )
            or bound > expected_tolerance
        ):
            reasons.append(bound_name)


def _population_growth_reasons(
    row: dict[str, Any],
    reasons: list[str],
) -> None:
    growth_mode = str(row.get("solver_population_growth_gate_mode", "")).strip()
    growth_tolerance = _required_float(
        row,
        "solver_population_growth_relative_tolerance",
    )
    if not math.isclose(
        growth_tolerance,
        MC_POPULATION_GROWTH_RELATIVE_TOLERANCE,
        rel_tol=1.0e-12,
        abs_tol=0.0,
    ):
        reasons.append("solver_population_growth_relative_tolerance")
    if growth_mode == "positive_paired_log":
        bound = _required_float(
            row,
            "solver_population_growth_max_relative_ci95_bound",
        )
        interval_ratio = _required_float(
            row,
            "solver_population_growth_poisson_interval_max_ratio",
        )
        sparse_metric = _required_float(
            row,
            "solver_population_growth_sparse_max_metric",
        )
        if bound < 0.0 or bound > MC_POPULATION_GROWTH_RELATIVE_TOLERANCE:
            reasons.append("solver_population_growth_max_relative_ci95_bound")
        if interval_ratio < 0.0:
            reasons.append("solver_population_growth_poisson_interval_max_ratio")
        if sparse_metric < 0.0:
            reasons.append("solver_population_growth_sparse_max_metric")
    elif growth_mode == "sparse_negligible_transport":
        interval_ratio = _required_float(
            row,
            "solver_population_growth_poisson_interval_max_ratio",
        )
        sparse_metric = _required_float(
            row,
            "solver_population_growth_sparse_max_metric",
        )
        if interval_ratio < 0.0 or interval_ratio > 1.0:
            reasons.append("solver_population_growth_poisson_interval_max_ratio")
        if (
            sparse_metric < 0.0
            or sparse_metric > MC_POPULATION_GROWTH_RELATIVE_TOLERANCE
        ):
            reasons.append("solver_population_growth_sparse_max_metric")
        paired_text = str(
            row.get("solver_population_growth_max_relative_ci95_bound", "")
        ).strip()
        if paired_text:
            paired_bound = _required_float(
                row,
                "solver_population_growth_max_relative_ci95_bound",
            )
            if paired_bound < 0.0:
                reasons.append("solver_population_growth_max_relative_ci95_bound")
    else:
        reasons.append("solver_population_growth_gate_mode")


def _mc_anchor_coverage(
    quality_fields: list[float],
    *,
    mc_sampling_plan: list[dict[str, Any]] | None,
    source_policy: dict[str, Any] | None,
    transport_rows: list[dict[str, Any]] | None,
    valid_ranges: dict[str, Any] | None,
    restricted_lmea: bool,
) -> dict[str, Any]:
    if mc_sampling_plan is None:
        raise GecCcpWorkflowError(
            "Monte Carlo quality audit lacks sampling-plan provenance"
        )
    sampling_fields_raw = [float(row["e_over_n_Td"]) for row in mc_sampling_plan]
    if len(sampling_fields_raw) != len(set(sampling_fields_raw)):
        raise GecCcpWorkflowError(
            "Monte Carlo sampling plan contains duplicate E/N anchors"
        )
    sampling_fields = sorted(sampling_fields_raw)
    qualified_fields = sorted(quality_fields)
    sampling_set = set(sampling_fields)
    qualified_set = set(qualified_fields)
    if not qualified_set.issubset(sampling_set):
        raise GecCcpWorkflowError(
            "Monte Carlo quality contains anchors outside the sampling plan"
        )
    first = sampling_fields.index(qualified_fields[0])
    last = sampling_fields.index(qualified_fields[-1])
    expected_qualified = sampling_fields[first : last + 1]
    if qualified_fields != expected_qualified:
        internal_failures = sorted(set(expected_qualified) - qualified_set)
        raise GecCcpWorkflowError(
            "Monte Carlo transport support crosses failed internal planned "
            "E/N anchors: " + ", ".join(f"{value:.17g}" for value in internal_failures)
        )
    excluded_fields = sorted(sampling_set - qualified_set)
    if source_policy is not None:
        _validate_source_policy_anchors(
            source_policy,
            sampling_fields=sampling_fields,
            qualified_fields=qualified_fields,
            excluded_fields=excluded_fields,
            restricted_lmea=restricted_lmea,
        )
        if valid_ranges is None:
            raise GecCcpWorkflowError(
                "Monte Carlo bundle lacks valid E/N transport support"
            )
    if transport_rows is not None:
        _validate_transport_anchors(
            transport_rows,
            qualified_fields=qualified_fields,
            source_policy=source_policy,
            restricted_lmea=restricted_lmea,
        )
    if valid_ranges is not None:
        _validate_valid_range(valid_ranges, qualified_fields)
    return {
        "passed": True,
        "policy": (
            "qualified_active_closure_anchors_are_a_contiguous_"
            "subsequence_of_planned_E_over_N"
            if restricted_lmea
            else (
                "qualified_transport_anchors_are_a_contiguous_"
                "subsequence_of_planned_E_over_N"
            )
        ),
        "planned_E_over_N_Td": sampling_fields,
        "qualified_E_over_N_Td": qualified_fields,
        "exported_support_E_over_N_Td": [
            qualified_fields[0],
            qualified_fields[-1],
        ],
        "failed_outside_exported_support_E_over_N_Td": excluded_fields,
        "failed_internal_E_over_N_Td": [],
        "interpolation_across_failed_planned_anchor": False,
    }


def _validate_source_policy_anchors(
    source_policy: dict[str, Any],
    *,
    sampling_fields: list[float],
    qualified_fields: list[float],
    excluded_fields: list[float],
    restricted_lmea: bool,
) -> None:
    raw_eligible = source_policy.get("transport_eligible_E_over_N_Td")
    raw_excluded = source_policy.get("transport_excluded_E_over_N_Td")
    if not isinstance(raw_eligible, list) or not isinstance(raw_excluded, dict):
        raise GecCcpWorkflowError(
            "Monte Carlo source_policy lacks canonical transport "
            "anchor qualification metadata"
        )
    try:
        policy_eligible = [float(value) for value in raw_eligible]
        policy_excluded = {float(key): reasons for key, reasons in raw_excluded.items()}
    except (TypeError, ValueError) as exc:
        raise GecCcpWorkflowError(
            "Monte Carlo source_policy transport anchors are invalid"
        ) from exc
    if policy_eligible != qualified_fields:
        raise GecCcpWorkflowError(
            "Monte Carlo source_policy eligible anchors disagree with quality.csv"
        )
    if sorted(policy_excluded) != excluded_fields or any(
        not isinstance(reasons, list)
        or not reasons
        or any(not isinstance(reason, str) or not reason.strip() for reason in reasons)
        for reasons in policy_excluded.values()
    ):
        raise GecCcpWorkflowError(
            "Monte Carlo source_policy excluded anchors are incomplete or inconsistent"
        )
    if source_policy.get("transport_eligible_anchor_points") != len(qualified_fields):
        raise GecCcpWorkflowError(
            "Monte Carlo source_policy eligible anchor count is inconsistent"
        )
    if restricted_lmea:
        _validate_restricted_lmea_policy(source_policy, sampling_fields)


def _validate_restricted_lmea_policy(
    source_policy: dict[str, Any],
    sampling_fields: list[float],
) -> None:
    if source_policy.get("qualification_outputs") != [
        "function_eedf",
        "reduced_mobility",
        "direct_elastic_energy_loss",
    ]:
        raise GecCcpWorkflowError(
            "Monte Carlo restricted-LMEA qualification-output metadata is inconsistent"
        )
    if source_policy.get("additional_mc_evidence") != [
        "particle_diffusion_L_T",
        "energy_mobility",
        "restricted_energy_diffusion_L_T",
        "direct_reaction_rates",
    ]:
        raise GecCcpWorkflowError(
            "Monte Carlo restricted-LMEA additional evidence metadata is inconsistent"
        )
    raw_full_eligible = source_policy.get("full_transport_eligible_E_over_N_Td")
    raw_full_excluded = source_policy.get("full_transport_excluded_E_over_N_Td")
    if not isinstance(raw_full_eligible, list) or not isinstance(
        raw_full_excluded, dict
    ):
        raise GecCcpWorkflowError(
            "Monte Carlo restricted-LMEA bundle lacks retained "
            "full-transport qualification evidence"
        )
    sampling_set = set(sampling_fields)
    full_eligible = sorted(float(value) for value in raw_full_eligible)
    full_excluded = sorted(float(value) for value in raw_full_excluded)
    if (
        any(value not in sampling_set for value in full_eligible)
        or any(value not in sampling_set for value in full_excluded)
        or sorted(set(full_eligible).union(full_excluded)) != sampling_fields
        or source_policy.get("full_transport_eligible_anchor_points")
        != len(full_eligible)
    ):
        raise GecCcpWorkflowError(
            "Monte Carlo retained full-transport qualification anchors are inconsistent"
        )


def _validate_transport_anchors(
    transport_rows: list[dict[str, Any]],
    *,
    qualified_fields: list[float],
    source_policy: dict[str, Any] | None,
    restricted_lmea: bool,
) -> None:
    transport_fields = [_required_float(row, "E_over_N_Td") for row in transport_rows]
    if (
        len(transport_fields) != len(set(transport_fields))
        or sorted(transport_fields) != qualified_fields
    ):
        raise GecCcpWorkflowError(
            "Monte Carlo transport table anchors disagree with quality.csv"
        )
    if not restricted_lmea or source_policy is None:
        return
    means = [_required_float(row, "mean_energy_eV") for row in transport_rows]
    ratios = [right / left for left, right in zip(means, means[1:])]
    gaps = [right - left for left, right in zip(means, means[1:])]
    reported_means = source_policy.get("active_mean_energy_anchors_eV")
    try:
        reported_ratio = float(source_policy["maximum_adjacent_mean_energy_ratio"])
        reported_gap = float(source_policy["maximum_adjacent_mean_energy_gap_eV"])
    except (KeyError, TypeError, ValueError) as exc:
        raise GecCcpWorkflowError(
            "Monte Carlo restricted-LMEA mean-energy anchor spacing metadata is invalid"
        ) from exc
    if (
        reported_means != means
        or not math.isclose(
            reported_ratio,
            max(ratios),
            rel_tol=1.0e-12,
            abs_tol=0.0,
        )
        or not math.isclose(
            reported_gap,
            max(gaps),
            rel_tol=1.0e-12,
            abs_tol=0.0,
        )
    ):
        raise GecCcpWorkflowError(
            "Monte Carlo restricted-LMEA mean-energy anchor spacing metadata "
            "is inconsistent"
        )


def _validate_valid_range(
    valid_ranges: dict[str, Any],
    qualified_fields: list[float],
) -> None:
    raw_range = valid_ranges.get("E_over_N_Td")
    if not isinstance(raw_range, list) or len(raw_range) != 2:
        raise GecCcpWorkflowError(
            "Monte Carlo valid_ranges lacks its E/N transport support"
        )
    support = [_float_or_none(value) for value in raw_range]
    expected_support = [qualified_fields[0], qualified_fields[-1]]
    if any(value is None for value in support) or any(
        not math.isclose(
            float(actual),
            expected,
            rel_tol=1.0e-12,
            abs_tol=0.0,
        )
        for actual, expected in zip(support, expected_support, strict=True)
    ):
        raise GecCcpWorkflowError(
            "Monte Carlo valid E/N range disagrees with qualified transport support"
        )


def _integer(row: dict[str, Any], name: str) -> int:
    try:
        return int(str(row[name]))
    except (KeyError, TypeError, ValueError) as exc:
        raise GecCcpWorkflowError(
            f"bundle quality column {name} must contain integers"
        ) from exc


def _true_value(row: dict[str, Any], name: str) -> bool:
    return str(row.get(name, "")).strip().lower() in {"1", "true"}


__all__ = ["independent_bundle_quality_audit"]
