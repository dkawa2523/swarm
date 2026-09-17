"""Bounded Monte Carlo convergence and whole-closure selection policy."""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass, replace
from enum import Enum
import json
import math
from pathlib import Path
from typing import Any, Literal

from electron_swarm.solvers.monte_carlo.source_identity import (
    monte_carlo_source_sha256,
)

from ...selection import (
    MC_QUALIFICATION_TABLE,
    ClosureSelectionError,
    compatible_physics,
    file_sha256,
    read_selection,
    table_evidence,
)


Action = Literal[
    "accept_monte_carlo",
    "extend_time",
    "add_replicas",
    "increase_particles",
    "extend_tail",
    "select_two_term",
    "blocked",
]


class FailureAxis(str, Enum):
    """Finite remediation axes emitted by MC qualification."""

    STRUCTURAL = "structural"
    TIME_HORIZON = "time_horizon"
    REPLICA_PRECISION = "replica_precision"
    POPULATION_LINEAGE = "population_lineage"
    RARE_TAIL = "rare_tail"


class MonteCarloPolicyError(RuntimeError):
    """Raised when solver-selection evidence is incomplete or inconsistent."""


@dataclass(frozen=True, slots=True)
class MonteCarloConvergencePolicy:
    """Small, solver-independent budget for an MC qualification campaign."""

    initial_replicas: int = 4
    replica_increment: int = 4
    maximum_replicas: int = 8
    maximum_attempts: int = 2
    time_extension_factor: int = 4
    particle_increase_factor: int = 2
    tail_extension_factor: int = 4
    maximum_particle_barriers_per_replica: int = 500_000_000
    maximum_total_particle_barriers: int = 10_000_000_000


@dataclass(frozen=True, slots=True)
class SamplingPlanEntry:
    e_over_n_Td: float
    particles: int
    warmup_collisions: int
    max_collisions: int
    tail_max_collisions: int
    replicas: int
    transport_correlation_lag_barriers: int
    transport_estimator: Literal["single_field", "paired_field_parity"]


@dataclass(frozen=True, slots=True)
class PolicyDecisionSummary:
    artifact: Path
    action: Action
    selected_solver: str | None


@dataclass(frozen=True, slots=True)
class _AnchorEvidence:
    e_over_n_Td: float
    passed: bool
    failure_reasons: tuple[str, ...]
    failure_axes: tuple[FailureAxis, ...]
    valid_replicates: int | None


_TIME_HORIZON_FAILURE_CODES = frozenset(
    {
        "mean_energy_stationarity",
        "mc_origin_stationarity_not_converged",
        "mc_transport_case_mean_energy_inconsistent",
        "mc_transport_lag_plateau_not_converged",
        "mc_weighted_growth_lag_not_converged",
        "mc_weighted_growth_population_rate_unresolved",
        "mc_weighted_growth_stationarity_not_converged",
        "solver_lag_convergence_max_relative_ci95_bound",
        "solver_mean_energy_stationarity_relative_ci95_bound",
        "solver_mobility_stationarity_absolute_log_drift",
        "solver_origin_stationarity_limiting_relative_ci95_bound",
        "solver_population_growth_max_relative_ci95_bound",
        "solver_population_growth_poisson_interval_max_ratio",
        "solver_population_growth_sparse_max_metric",
        "solver_transport_mean_energy_max_relative_ci95_bound",
    }
)
_REPLICA_PRECISION_FAILURE_CODES = frozenset(
    {
        "diffusion_L_rse",
        "diffusion_T_rse",
        "energy_diffusion_L_rse",
        "energy_diffusion_T_rse",
        "energy_mobility_rse",
        "mc_transport_ensemble_replicates_insufficient",
        "mobility_relative_ci95_half_width",
        "mobility_rse",
        "reduced_diffusion_L_rse_unavailable_or_exceeds_threshold",
        "reduced_diffusion_T_rse_unavailable_or_exceeds_threshold",
        "reduced_energy_diffusion_L_rse_unavailable_or_exceeds_threshold",
        "reduced_energy_diffusion_T_rse_unavailable_or_exceeds_threshold",
        "reduced_energy_mobility_rse_unavailable_or_exceeds_threshold",
        "reduced_mobility_rse_unavailable_or_exceeds_threshold",
        "solver_transport_replicates",
        "uncertainty_available",
        "valid_replicates",
    }
)
_POPULATION_LINEAGE_FAILURE_CODES = frozenset(
    {
        "mc_weighted_growth_lineage_degenerate",
    }
)
_RARE_TAIL_FAILURE_CODES = frozenset(
    {
        "major_rate_censored_all_zero",
        "major_rate_rse_exceeds_threshold",
        "major_rate_rse_unavailable",
        "max_major_rate_rse",
        "mc_censored_rate_relevance_failed",
        "mc_unresolved_rate_tail",
        "required_rate_censored_all_zero",
        "required_rate_rse_exceeds_threshold",
        "required_rate_rse_unavailable",
        "unresolved_rate_tail",
    }
)


def failure_axes_for_reasons(
    reasons: tuple[str, ...] | list[str],
) -> tuple[FailureAxis, ...]:
    """Route stable failure codes to remediation axes.

    Dynamic details follow the first ``:`` in a reason and never participate
    in policy routing. Unknown codes are structural so new diagnostics cannot
    accidentally authorize more compute.
    """

    axes: set[FailureAxis] = set()
    for reason in reasons:
        code = reason.partition(":")[0]
        if code in _TIME_HORIZON_FAILURE_CODES:
            axes.add(FailureAxis.TIME_HORIZON)
        elif code in _REPLICA_PRECISION_FAILURE_CODES:
            axes.add(FailureAxis.REPLICA_PRECISION)
        elif code in _POPULATION_LINEAGE_FAILURE_CODES:
            axes.add(FailureAxis.POPULATION_LINEAGE)
        elif code in _RARE_TAIL_FAILURE_CODES:
            axes.add(FailureAxis.RARE_TAIL)
        else:
            axes.add(FailureAxis.STRUCTURAL)
    return tuple(sorted(axes, key=lambda axis: axis.value))


def decide_monte_carlo_closure(
    mc_table_directory: str | Path,
    two_term_table_directory: str | Path,
    *,
    attempt: int,
    output: str | Path,
    policy: MonteCarloConvergencePolicy | None = None,
    previous_decision: str | Path | None = None,
) -> PolicyDecisionSummary:
    """Decide the next bounded MC action or one complete closure source.

    The function never edits coefficient tables and never joins MC and
    two-term anchors.  It consumes immutable table evidence and writes one
    auditable decision artifact.
    """

    resolved_policy = policy or MonteCarloConvergencePolicy()
    _validate_policy(resolved_policy)
    if isinstance(attempt, bool) or not isinstance(attempt, int) or attempt < 1:
        raise MonteCarloPolicyError(
            f"attempt must be between 1 and {resolved_policy.maximum_attempts}"
        )

    mc_dir = Path(mc_table_directory).resolve()
    fallback_dir = Path(two_term_table_directory).resolve()
    mc_manifest_path = mc_dir / "manifest.json"
    fallback_manifest_path = fallback_dir / "manifest.json"
    mc_quality_path = mc_dir / MC_QUALIFICATION_TABLE
    if not mc_quality_path.is_file():
        mc_quality_path = mc_dir / "quality.csv"
    fallback_quality_path = fallback_dir / "quality.csv"
    mc_manifest = _read_manifest(mc_manifest_path, expected_source="monte_carlo")
    fallback_manifest = _read_manifest(
        fallback_manifest_path, expected_source="two_term"
    )
    campaign = mc_manifest.get("mc_campaign")
    if policy is None and isinstance(campaign, dict):
        resolved_policy = parse_convergence_policy(campaign.get("policy"))
    if attempt > resolved_policy.maximum_attempts:
        raise MonteCarloPolicyError("MC attempt budget is exhausted")
    try:
        compatible_physics(mc_manifest, fallback_manifest)
        inputs = {
            "monte_carlo": table_evidence(mc_dir, mc_manifest, mc_quality_path),
            "two_term": table_evidence(
                fallback_dir, fallback_manifest, fallback_quality_path
            ),
        }
    except ClosureSelectionError as exc:
        raise MonteCarloPolicyError(str(exc)) from exc

    plan = _sampling_plan(mc_manifest)
    actual_attempt = validate_sampling_budget(
        plan,
        resolved_policy,
        previous_decision=previous_decision,
    )
    if actual_attempt != attempt:
        raise MonteCarloPolicyError(
            "attempt does not follow the recorded MC decision history"
        )
    if previous_decision is not None:
        previous = read_selection(previous_decision)
        if previous.get("physical_context") != mc_manifest.get("physical_context"):
            raise MonteCarloPolicyError("MC physical context changed between attempts")
        if previous.get("mc_run_identity") != _mc_run_identity(mc_manifest):
            raise MonteCarloPolicyError(
                "MC config, estimator, seed, or cross sections changed between attempts"
            )
    campaign = mc_manifest.get("mc_campaign")
    if campaign is not None:
        expected_campaign = campaign_provenance(resolved_policy, previous_decision)
        try:
            campaign_policy = parse_convergence_policy(campaign.get("policy"))
        except (AttributeError, MonteCarloPolicyError) as exc:
            raise MonteCarloPolicyError(
                "table campaign contains an invalid convergence policy"
            ) from exc
        if (
            campaign_policy != resolved_policy
            or campaign.get("previous_decision")
            != expected_campaign["previous_decision"]
        ):
            raise MonteCarloPolicyError(
                "table campaign differs from the requested decision history"
            )
    quality_scope, evidence_by_anchor = _mc_quality_evidence(
        mc_quality_path,
        mc_manifest,
    )
    plan_by_anchor = {entry.e_over_n_Td: entry for entry in plan}
    if set(plan_by_anchor) != set(evidence_by_anchor):
        missing = sorted(set(plan_by_anchor).difference(evidence_by_anchor))
        extra = sorted(set(evidence_by_anchor).difference(plan_by_anchor))
        raise MonteCarloPolicyError(
            "MC sampling-plan and quality anchors disagree; "
            f"missing quality={missing}, unplanned quality={extra}"
        )
    for anchor, evidence in evidence_by_anchor.items():
        planned = plan_by_anchor[anchor].replicas
        if evidence.valid_replicates is None or evidence.valid_replicates > planned:
            raise MonteCarloPolicyError(f"invalid replica evidence at {anchor:g} Td")
        if evidence.passed and evidence.valid_replicates != planned:
            raise MonteCarloPolicyError(
                f"passing MC anchor {anchor:g} Td lacks all planned replicas"
            )

    fallback_qualified, fallback_failures = _fallback_qualification(
        fallback_quality_path,
        required_anchors=tuple(sorted(plan_by_anchor)),
    )
    failed = tuple(
        evidence_by_anchor[anchor]
        for anchor in sorted(evidence_by_anchor)
        if not evidence_by_anchor[anchor].passed
    )
    action, reason, selected_solver, next_plan = _decide(
        failed,
        plan,
        attempt=attempt,
        fallback_qualified=fallback_qualified,
        policy=resolved_policy,
    )
    if (
        action == "accept_monte_carlo"
        and mc_manifest.get("mc_qualification", {}).get("coefficient_tables_available")
        is False
    ):
        action, reason, selected_solver, next_plan = (
            "blocked",
            "mc_closure_tables_unavailable",
            None,
            None,
        )
    status = {
        "accept_monte_carlo": "selected",
        "select_two_term": "selected",
        "extend_time": "requires_mc_rerun",
        "add_replicas": "requires_mc_rerun",
        "increase_particles": "requires_mc_rerun",
        "extend_tail": "requires_mc_rerun",
        "blocked": "blocked",
    }[action]

    payload: dict[str, Any] = {
        "format_version": 2,
        "stage": "mc-solver-selection",
        "status": status,
        "requested_primary_solver": "monte_carlo",
        "selection_scope": "whole_comsol_closure",
        "quality_scope": quality_scope,
        "attempt": attempt,
        "action": action,
        "reason": reason,
        "selected_solver": selected_solver,
        "failed_anchors_Td": [row.e_over_n_Td for row in failed],
        "mc_results_retained_as_validation": selected_solver == "two_term",
        "policy": asdict(resolved_policy),
        "nominal_compute_budget": sampling_budget_provenance(plan),
        "previous_decision": (
            {
                "path": str(Path(previous_decision).resolve()),
                "sha256": file_sha256(previous_decision),
            }
            if previous_decision is not None
            else None
        ),
        "physical_context": mc_manifest["physical_context"],
        "mc_run_identity": _mc_run_identity(mc_manifest),
        "anchor_evidence": [
            {
                "e_over_n_Td": row.e_over_n_Td,
                "passed": row.passed,
                "failure_reasons": list(row.failure_reasons),
                "failure_axes": [axis.value for axis in row.failure_axes],
                "valid_replicates": row.valid_replicates,
            }
            for row in (evidence_by_anchor[key] for key in sorted(evidence_by_anchor))
        ],
        "current_sampling_plan": [asdict(entry) for entry in plan],
        "next_sampling_plan": (
            [asdict(entry) for entry in next_plan] if next_plan is not None else None
        ),
        "reuse_previous_mc_database": action
        in {
            "add_replicas",
            "extend_tail",
            "extend_time",
            "increase_particles",
        },
        "fallback": {
            "solver": "two_term",
            "qualified_for_required_anchors": fallback_qualified,
            "failures": fallback_failures,
        },
        "inputs": inputs,
    }
    output_path = Path(output).resolve()
    encoded = json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if output_path.exists() and output_path.read_text(encoding="utf-8") != encoded:
        raise MonteCarloPolicyError(
            "decision output already contains different evidence; use a new path"
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(encoded, encoding="utf-8")
    return PolicyDecisionSummary(output_path, action, selected_solver)


def _validate_policy(policy: MonteCarloConvergencePolicy) -> None:
    values = asdict(policy)
    if any(
        isinstance(value, bool) or not isinstance(value, int) or value <= 0
        for value in values.values()
    ):
        raise MonteCarloPolicyError(
            "MC convergence-policy values must be positive integers"
        )
    if policy.initial_replicas > policy.maximum_replicas:
        raise MonteCarloPolicyError("initial_replicas cannot exceed maximum_replicas")
    if policy.time_extension_factor < 2:
        raise MonteCarloPolicyError("time_extension_factor must be at least two")
    if policy.particle_increase_factor < 2:
        raise MonteCarloPolicyError("particle_increase_factor must be at least two")
    if policy.tail_extension_factor < 2:
        raise MonteCarloPolicyError("tail_extension_factor must be at least two")


def parse_convergence_policy(value: object) -> MonteCarloConvergencePolicy:
    if not isinstance(value, dict):
        raise MonteCarloPolicyError("MC convergence policy must be a mapping")
    try:
        policy = MonteCarloConvergencePolicy(**value)
    except TypeError as exc:
        raise MonteCarloPolicyError("unknown MC convergence-policy field") from exc
    _validate_policy(policy)
    return policy


def validate_sampling_budget(
    plan: tuple[SamplingPlanEntry, ...],
    policy: MonteCarloConvergencePolicy,
    *,
    previous_decision: str | Path | None = None,
) -> int:
    """Check the executable plan before spending compute, and again at selection."""
    _validate_policy(policy)
    if not plan or any(entry.replicas > policy.maximum_replicas for entry in plan):
        raise MonteCarloPolicyError(
            "MC sampling plan exceeds the replica budget or is empty"
        )
    costs = sampling_budget_provenance(plan)
    if int(costs["maximum_per_replica_particle_barriers"]) > (
        policy.maximum_particle_barriers_per_replica
    ):
        raise MonteCarloPolicyError(
            "MC sampling plan exceeds maximum_particle_barriers_per_replica: "
            f"{costs['maximum_per_replica_particle_barriers']} > "
            f"{policy.maximum_particle_barriers_per_replica}"
        )
    if int(costs["total_particle_barriers"]) > policy.maximum_total_particle_barriers:
        raise MonteCarloPolicyError(
            "MC sampling plan exceeds maximum_total_particle_barriers: "
            f"{costs['total_particle_barriers']} > "
            f"{policy.maximum_total_particle_barriers}"
        )
    if previous_decision is None:
        if any(entry.replicas < policy.initial_replicas for entry in plan):
            raise MonteCarloPolicyError(
                "first MC attempt must use at least initial_replicas at every anchor"
            )
        return 1
    try:
        previous = read_selection(previous_decision)
    except ClosureSelectionError as exc:
        raise MonteCarloPolicyError(str(exc)) from exc
    if previous.get("status") != "requires_mc_rerun" or previous.get("action") not in {
        "add_replicas",
        "extend_tail",
        "extend_time",
        "increase_particles",
    }:
        raise MonteCarloPolicyError(
            "previous decision does not authorize another MC attempt"
        )
    try:
        previous_policy = parse_convergence_policy(previous.get("policy"))
    except MonteCarloPolicyError as exc:
        raise MonteCarloPolicyError(
            "previous decision contains an invalid MC convergence policy"
        ) from exc
    if previous_policy != policy:
        raise MonteCarloPolicyError("MC convergence policy changed between attempts")
    if previous.get("next_sampling_plan") != [
        asdict(entry) for entry in sorted(plan, key=lambda row: row.e_over_n_Td)
    ]:
        raise MonteCarloPolicyError(
            "MC plan differs from the bounded next_sampling_plan"
        )
    attempt = _positive_int(previous.get("attempt")) + 1
    if attempt > policy.maximum_attempts:
        raise MonteCarloPolicyError("MC attempt budget is exhausted")
    return attempt


def sampling_budget_provenance(
    plan: tuple[SamplingPlanEntry, ...],
) -> dict[str, object]:
    """Return the nominal upper budget, including a possible tail phase."""

    rows = []
    for entry in sorted(plan, key=lambda item: item.e_over_n_Td):
        transport_field_legs = (
            2 if entry.transport_estimator == "paired_field_parity" else 1
        )
        per_replica = int(
            entry.particles
            * (
                transport_field_legs * (entry.warmup_collisions + entry.max_collisions)
                + entry.tail_max_collisions
            )
        )
        rows.append(
            {
                "e_over_n_Td": float(entry.e_over_n_Td),
                "per_replica_particle_barriers": per_replica,
                "replicas": int(entry.replicas),
                "total_particle_barriers": per_replica * int(entry.replicas),
                "transport_estimator": entry.transport_estimator,
                "transport_field_legs": transport_field_legs,
                "tail_phase_is_conditional": entry.tail_max_collisions > 0,
            }
        )
    return {
        "unit": "nominal_particle_barriers",
        "includes_conditional_tail_upper_bound": any(
            entry.tail_max_collisions > 0 for entry in plan
        ),
        "maximum_per_replica_particle_barriers": max(
            (int(row["per_replica_particle_barriers"]) for row in rows),
            default=0,
        ),
        "total_particle_barriers": sum(
            int(row["total_particle_barriers"]) for row in rows
        ),
        "entries": rows,
    }


def _sampling_plan_within_compute_budget(
    plan: tuple[SamplingPlanEntry, ...],
    policy: MonteCarloConvergencePolicy,
) -> bool:
    costs = sampling_budget_provenance(plan)
    return bool(
        int(costs["maximum_per_replica_particle_barriers"])
        <= policy.maximum_particle_barriers_per_replica
        and int(costs["total_particle_barriers"])
        <= policy.maximum_total_particle_barriers
    )


def campaign_provenance(
    policy: MonteCarloConvergencePolicy,
    previous_decision: str | Path | None,
) -> dict[str, Any]:
    return {
        "policy": asdict(policy),
        "previous_decision": (
            {
                "path": str(Path(previous_decision).resolve()),
                "sha256": file_sha256(previous_decision),
            }
            if previous_decision is not None
            else None
        ),
    }


def _mc_run_identity(manifest: dict[str, Any]) -> dict[str, Any]:
    hashes = manifest.get("hashes", {})
    if any(
        not hashes.get(key)
        for key in (
            "base_config_sha256",
            "cross_sections_sha256",
            "mc_transport_estimator_schema_version",
            "mc_eedf_estimator_schema_version",
            "mc_solver_source_sha256",
        )
    ) or not isinstance(manifest.get("mc_sampling_plan", {}).get("base_seed"), int):
        raise MonteCarloPolicyError(
            "MC selection requires config, estimator, solver-source, and verified "
            "base-seed provenance"
        )
    if hashes["mc_solver_source_sha256"] != monte_carlo_source_sha256():
        raise MonteCarloPolicyError(
            "Monte Carlo implementation changed after table generation"
        )
    return {
        **{
            key: hashes.get(key)
            for key in (
                "base_config_sha256",
                "cross_sections_sha256",
                "mc_transport_estimator_schema_version",
                "mc_eedf_estimator_schema_version",
                "mc_tail_estimator_schema_version",
                "mc_seed_derivation_schema_version",
                "mc_solver_source_sha256",
            )
        },
        "mixture": manifest.get("mixture"),
        "base_seed": manifest.get("mc_sampling_plan", {}).get("base_seed"),
    }


def _read_manifest(path: Path, *, expected_source: str) -> dict[str, Any]:
    if not path.is_file():
        raise MonteCarloPolicyError(f"table manifest does not exist: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise MonteCarloPolicyError(f"cannot read table manifest: {path}") from exc
    if not isinstance(value, dict):
        raise MonteCarloPolicyError(f"table manifest must be a mapping: {path}")
    if value.get("source") != expected_source:
        raise MonteCarloPolicyError(
            f"table manifest source must be {expected_source}: {path}"
        )
    return value


def _sampling_plan(manifest: dict[str, Any]) -> tuple[SamplingPlanEntry, ...]:
    payload = manifest.get("mc_sampling_plan")
    rows = payload.get("entries") if isinstance(payload, dict) else None
    if not isinstance(rows, list) or not rows:
        raise MonteCarloPolicyError(
            "monte_carlo table manifest lacks a validated mc_sampling_plan"
        )
    entries: list[SamplingPlanEntry] = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict):
            raise MonteCarloPolicyError(
                f"mc_sampling_plan.entries[{index}] must be a mapping"
            )
        try:
            entry = SamplingPlanEntry(
                e_over_n_Td=_positive_float(row["e_over_n_Td"]),
                particles=_positive_int(row["particles"]),
                warmup_collisions=_nonnegative_int(row["warmup_collisions"]),
                max_collisions=_positive_int(row["max_collisions"]),
                tail_max_collisions=_nonnegative_int(row["tail_max_collisions"]),
                replicas=_positive_int(row["replicas"]),
                transport_correlation_lag_barriers=_positive_int(
                    row["transport_correlation_lag_barriers"]
                ),
                transport_estimator=_transport_estimator(row["transport_estimator"]),
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise MonteCarloPolicyError(
                f"mc_sampling_plan.entries[{index}] is invalid"
            ) from exc
        entries.append(entry)
    anchors = [entry.e_over_n_Td for entry in entries]
    if len(anchors) != len(set(anchors)):
        raise MonteCarloPolicyError("mc_sampling_plan contains duplicate anchors")
    return tuple(sorted(entries, key=lambda item: item.e_over_n_Td))


def _mc_quality_evidence(
    path: Path,
    manifest: dict[str, Any],
) -> tuple[str, dict[float, _AnchorEvidence]]:
    rows = _read_csv(path)
    profile = manifest.get("source_policy", {}).get("qualification_profile")
    uses_active_scope = all(
        name in rows[0]
        for name in (
            "active_closure_quality_passed",
            "active_closure_failure_reasons_json",
        )
    )
    if uses_active_scope:
        passed_field = "active_closure_quality_passed"
        reasons_field = "active_closure_failure_reasons_json"
        axes_field = "active_closure_failure_axes_json"
        quality_scope = str(profile or "active_comsol_closure")
    else:
        passed_field = "passed"
        reasons_field = "failure_reasons_json"
        axes_field = "failure_axes_json"
        quality_scope = "full_transport"

    evidence: dict[float, _AnchorEvidence] = {}
    for index, row in enumerate(rows, start=2):
        try:
            anchor = _positive_float(row["E_over_N_Td"])
            passed = _boolean(row[passed_field])
            reasons = _reason_list(row[reasons_field])
            axes = _failure_axis_list(row[axes_field])
            valid_replicates = (
                _nonnegative_int(row["valid_replicates"])
                if row.get("valid_replicates", "") != ""
                else None
            )
        except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
            raise MonteCarloPolicyError(
                f"invalid MC quality evidence at {path}:{index}"
            ) from exc
        if passed and reasons:
            raise MonteCarloPolicyError(
                f"passing MC anchor {anchor:g} Td contains failure reasons"
            )
        if passed and axes:
            raise MonteCarloPolicyError(
                f"passing MC anchor {anchor:g} Td contains failure axes"
            )
        if not passed and not reasons:
            reasons = ("quality_failure_without_reason",)
        if not passed and not axes:
            raise MonteCarloPolicyError(
                f"failed MC anchor {anchor:g} Td lacks typed failure axes"
            )
        if anchor in evidence:
            raise MonteCarloPolicyError(
                f"MC quality table contains duplicate anchor {anchor:g} Td"
            )
        evidence[anchor] = _AnchorEvidence(
            anchor,
            passed,
            reasons,
            axes,
            valid_replicates,
        )
    return quality_scope, evidence


def _fallback_qualification(
    path: Path,
    *,
    required_anchors: tuple[float, ...],
) -> tuple[bool, dict[str, list[str]]]:
    rows = _read_csv(path)
    indexed: dict[float, dict[str, str]] = {}
    for row in rows:
        anchor = _positive_float(row["E_over_N_Td"])
        if anchor in indexed:
            raise MonteCarloPolicyError(
                f"two_term quality table contains duplicate anchor {anchor:g} Td"
            )
        indexed[anchor] = row

    sorted_anchors = sorted(indexed)
    failures: dict[str, list[str]] = {}
    for anchor, row in indexed.items():
        passed = _boolean(row.get("passed", "0"))
        reasons = _reason_list(row.get("failure_reasons_json", "[]"))
        if passed and reasons:
            raise MonteCarloPolicyError(
                "passing two_term anchor contains failure reasons"
            )
        if not passed:
            failures[f"{anchor:.17g}"] = list(
                reasons or ("quality_failure_without_reason",)
            )
    for required in required_anchors:
        key = f"{required:.17g}"
        exact = next(
            (
                anchor
                for anchor in sorted_anchors
                if math.isclose(anchor, required, rel_tol=1.0e-12, abs_tol=1.0e-12)
            ),
            None,
        )
        if exact is not None:
            support = (exact,)
        else:
            lower = [anchor for anchor in sorted_anchors if anchor < required]
            upper = [anchor for anchor in sorted_anchors if anchor > required]
            if not lower or not upper:
                failures[key] = ["required_anchor_outside_two_term_support"]
                continue
            support = (lower[-1], upper[0])
        failed_support: list[str] = []
        for anchor in support:
            row = indexed[anchor]
            if _boolean(row.get("passed", "0")):
                continue
            reasons = _reason_list(row.get("failure_reasons_json", "[]"))
            failed_support.extend(reasons or ("quality_failure_without_reason",))
        if failed_support:
            failures[key] = list(dict.fromkeys(failed_support))
    return not failures, failures


def _decide(
    failed: tuple[_AnchorEvidence, ...],
    plan: tuple[SamplingPlanEntry, ...],
    *,
    attempt: int,
    fallback_qualified: bool,
    policy: MonteCarloConvergencePolicy,
) -> tuple[Action, str, str | None, tuple[SamplingPlanEntry, ...] | None]:
    if not failed:
        return (
            "accept_monte_carlo",
            "mc_all_required_anchors_qualified",
            "monte_carlo",
            None,
        )

    failure_axes = {axis for row in failed for axis in row.failure_axes}
    if FailureAxis.STRUCTURAL in failure_axes:
        return (
            "blocked",
            "mc_failure_requires_model_or_evidence_correction",
            None,
            None,
        )

    if attempt < policy.maximum_attempts:
        candidates: tuple[tuple[FailureAxis, Action, str], ...] = (
            (
                FailureAxis.POPULATION_LINEAGE,
                "increase_particles",
                "mc_population_budget_not_yet_exhausted",
            ),
            (
                FailureAxis.RARE_TAIL,
                "extend_tail",
                "mc_tail_budget_not_yet_exhausted",
            ),
            (
                FailureAxis.TIME_HORIZON,
                "extend_time",
                "mc_time_horizon_budget_not_yet_exhausted",
            ),
            (
                FailureAxis.REPLICA_PRECISION,
                "add_replicas",
                "mc_precision_budget_not_yet_exhausted",
            ),
        )
        for axis, action, reason in candidates:
            affected = {row.e_over_n_Td for row in failed if axis in row.failure_axes}
            if not affected:
                continue
            next_plan = _remediation_plan(
                plan,
                affected_anchors=affected,
                action=action,
                policy=policy,
            )
            if next_plan is not None and _sampling_plan_within_compute_budget(
                next_plan, policy
            ):
                return action, reason, None, next_plan

    if not fallback_qualified:
        return "blocked", "two_term_fallback_not_qualified", None, None
    return "select_two_term", "mc_budget_exhausted", "two_term", None


def _remediation_plan(
    plan: tuple[SamplingPlanEntry, ...],
    *,
    affected_anchors: set[float],
    action: Action,
    policy: MonteCarloConvergencePolicy,
) -> tuple[SamplingPlanEntry, ...] | None:
    affected = [entry for entry in plan if entry.e_over_n_Td in affected_anchors]
    if not affected:
        return None
    if action == "add_replicas" and any(
        entry.replicas >= policy.maximum_replicas for entry in affected
    ):
        return None
    if action == "extend_tail" and any(
        entry.tail_max_collisions == 0 for entry in affected
    ):
        return None
    updated: list[SamplingPlanEntry] = []
    for entry in plan:
        if entry.e_over_n_Td not in affected_anchors:
            updated.append(entry)
        elif action == "extend_time":
            updated.append(
                replace(
                    entry,
                    warmup_collisions=max(
                        policy.time_extension_factor,
                        entry.warmup_collisions * policy.time_extension_factor,
                    ),
                    max_collisions=(
                        entry.max_collisions * policy.time_extension_factor
                    ),
                )
            )
        elif action == "add_replicas":
            updated.append(
                replace(
                    entry,
                    replicas=min(
                        policy.maximum_replicas,
                        entry.replicas + policy.replica_increment,
                    ),
                )
            )
        elif action == "increase_particles":
            updated.append(
                replace(
                    entry,
                    particles=entry.particles * policy.particle_increase_factor,
                )
            )
        elif action == "extend_tail":
            updated.append(
                replace(
                    entry,
                    tail_max_collisions=(
                        entry.tail_max_collisions * policy.tail_extension_factor
                    ),
                )
            )
        else:
            raise MonteCarloPolicyError(f"unsupported MC remediation action: {action}")
    result = tuple(updated)
    return result if result != plan else None


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise MonteCarloPolicyError(f"quality table does not exist: {path}")
    try:
        with path.open("r", encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream))
    except OSError as exc:
        raise MonteCarloPolicyError(f"cannot read quality table: {path}") from exc
    if not rows:
        raise MonteCarloPolicyError(f"quality table contains no rows: {path}")
    return rows


def _reason_list(value: object) -> tuple[str, ...]:
    parsed = json.loads(str(value))
    if not isinstance(parsed, list) or any(
        not isinstance(item, str) for item in parsed
    ):
        raise ValueError("failure reasons must be a JSON string list")
    return tuple(dict.fromkeys(parsed))


def _failure_axis_list(value: object) -> tuple[FailureAxis, ...]:
    parsed = json.loads(str(value))
    if not isinstance(parsed, list) or any(
        not isinstance(item, str) for item in parsed
    ):
        raise ValueError("failure axes must be a JSON string list")
    try:
        axes = {FailureAxis(item) for item in parsed}
    except ValueError as exc:
        raise ValueError("unknown Monte Carlo failure axis") from exc
    return tuple(sorted(axes, key=lambda axis: axis.value))


def _positive_float(value: object) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError("value must be finite and positive")
    return number


def _transport_estimator(
    value: object,
) -> Literal["single_field", "paired_field_parity"]:
    if value not in {"single_field", "paired_field_parity"}:
        raise ValueError("invalid Monte Carlo transport estimator")
    return value


def _positive_int(value: object) -> int:
    number = _nonnegative_int(value)
    if number <= 0:
        raise ValueError("value must be positive")
    return number


def _nonnegative_int(value: object) -> int:
    if isinstance(value, bool):
        raise ValueError("boolean is not an integer")
    number = int(str(value))
    if float(value) != float(number) or number < 0:
        raise ValueError("value must be a nonnegative integer")
    return number


def _boolean(value: object) -> bool:
    normalized = str(value).strip().lower()
    if normalized in {"1", "true"}:
        return True
    if normalized in {"0", "false"}:
        return False
    raise ValueError("value must be boolean")
