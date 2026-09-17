"""Plan a provenance-bound, low-E/N Monte Carlo fallback.

The fallback is deliberately narrower than a whole-closure solver selection:
after at least one bounded MC follow-up attempt, an explicit low-field interval
may use qualified two-term anchors.  Every higher-field anchor remains the raw
MC result, including its unqualified status.  Values are never blended inside
one EEDF and a failed high-field MC result is never disguised as two-term data.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any

from swarm_workflow._io import write_json

from .closure import (
    ClosureSelectionError,
    compatible_physics,
    contained_file,
    file_sha256,
    read_json_object,
    read_selection,
)


ANCHOR_FALLBACK_SCHEMA = "swarm.low_field_anchor_fallback.v2"
ANCHOR_FALLBACK_SOURCE = "composite"
ANCHOR_FALLBACK_STAGE = "mc-anchor-fallback-plan"
COMPOSITE_QUALITY_SCHEMA = "composite_anchor_quality.v2"
COMPOSITE_QUALITY_COLUMNS = (
    "E_over_N_Td",
    "E_over_N_V_m2",
    "passed",
    "effective_source",
    "source_quality_gate",
    "selection_reason",
    "mc_failure_reasons_json",
    "source_manifest_sha256",
)
_MC_SOURCE = "monte_carlo"
_FALLBACK_SOURCE = "two_term"
_RESTRICTED_LMEA_PROFILE = "function_eedf_restricted_lmea"


class AnchorFallbackError(ClosureSelectionError):
    """The requested anchor fallback is not fully qualified or reproducible."""


@dataclass(frozen=True, slots=True)
class AnchorFallbackSummary:
    output_path: Path
    monte_carlo_anchors: int
    two_term_anchors: int
    retained_unqualified_monte_carlo_anchors: int


def build_mc_anchor_fallback_plan(
    mc_table_directory: str | Path,
    two_term_table_directory: str | Path,
    decision_path: str | Path,
    output_path: str | Path,
    *,
    maximum_fallback_e_over_n_Td: float,
) -> AnchorFallbackSummary:
    """Select two-term only for unresolved anchors in an explicit low-E/N band."""

    mc_dir = Path(mc_table_directory).resolve()
    fallback_dir = Path(two_term_table_directory).resolve()
    decision_file = Path(decision_path).resolve()
    output = Path(output_path).resolve()
    if output.exists():
        raise AnchorFallbackError(
            f"anchor-fallback output already exists; choose a new path: {output}"
        )
    if output in {
        decision_file,
        mc_dir / "manifest.json",
        fallback_dir / "manifest.json",
    }:
        raise AnchorFallbackError("anchor-fallback output must be a new artifact")

    cutoff = _positive_finite_cutoff(maximum_fallback_e_over_n_Td)
    decision = read_selection(decision_file)
    previous_decision = _require_bounded_fallback_decision(decision)
    mc_manifest = _verify_bound_source(mc_dir, decision, source=_MC_SOURCE)
    fallback_manifest = _verify_bound_source(
        fallback_dir, decision, source=_FALLBACK_SOURCE
    )
    compatible_physics(mc_manifest, fallback_manifest)
    _require_source_contracts(mc_manifest, fallback_manifest, decision)

    mc_quality = _read_bound_quality(mc_dir, decision, source=_MC_SOURCE)
    fallback_quality = _read_bound_quality(
        fallback_dir, decision, source=_FALLBACK_SOURCE
    )
    planned = _planned_anchors(decision)
    if cutoff >= max(planned):
        raise AnchorFallbackError(
            "the low-E/N fallback ceiling must leave at least one higher-field "
            "Monte Carlo anchor"
        )
    evidence = _decision_anchor_evidence(decision)
    if set(planned) != set(evidence) or set(planned) != set(mc_quality):
        raise AnchorFallbackError(
            "MC sampling plan, decision anchor evidence, and quality anchors disagree"
        )

    mc_manifest_sha = file_sha256(mc_dir / "manifest.json")
    fallback_manifest_sha = file_sha256(fallback_dir / "manifest.json")
    anchors: list[dict[str, Any]] = []
    mc_count = 0
    fallback_count = 0
    retained_unqualified_count = 0
    for field in planned:
        item = evidence[field]
        passed = bool(item["passed"])
        quality_passed, aggregate_passed, quality_reasons = _mc_quality_outcome(
            mc_quality[field], decision["quality_scope"]
        )
        decision_reasons = tuple(item["failure_reasons"])
        if passed != quality_passed or set(decision_reasons) != set(quality_reasons):
            raise AnchorFallbackError(
                f"decision and MC quality disagree at {field:.17g} Td"
            )
        if not aggregate_passed:
            raise AnchorFallbackError(
                "regional fallback cannot retain a structurally unusable raw MC "
                f"anchor at {field:.17g} Td"
            )
        if passed:
            source = _MC_SOURCE
            reason = "mc_anchor_qualified"
            source_manifest_sha = mc_manifest_sha
            mc_count += 1
        elif field <= cutoff:
            fallback_row = fallback_quality.get(field)
            if fallback_row is None or not _accepted_flag(fallback_row.get("passed")):
                raise AnchorFallbackError(
                    "two_term fallback requires an exact qualified anchor at "
                    f"{field:.17g} Td"
                )
            source = _FALLBACK_SOURCE
            reason = "bounded_mc_unresolved_in_low_e_over_n_fallback_band"
            source_manifest_sha = fallback_manifest_sha
            fallback_count += 1
        else:
            source = _MC_SOURCE
            reason = "mc_anchor_unqualified_retained_above_fallback_band"
            source_manifest_sha = mc_manifest_sha
            mc_count += 1
            retained_unqualified_count += 1
        anchors.append(
            {
                "E_over_N_Td": field,
                "effective_source": source,
                "selection_reason": reason,
                "mc_active_closure_quality_passed": passed,
                "mc_failure_reasons": list(decision_reasons),
                "source_manifest_sha256": source_manifest_sha,
            }
        )

    if fallback_count == 0:
        raise AnchorFallbackError(
            "the explicit low-E/N band contains no unresolved MC anchor"
        )
    failed = [field for field, item in evidence.items() if not item["passed"]]
    if failed != [float(value) for value in decision.get("failed_anchors_Td", ())]:
        raise AnchorFallbackError("decision failed-anchor summary is inconsistent")

    payload = {
        "format_version": 1,
        "stage": ANCHOR_FALLBACK_STAGE,
        "status": (
            "evidence_only" if retained_unqualified_count else "selected"
        ),
        "source": ANCHOR_FALLBACK_SOURCE,
        "source_composition": {
            "schema": ANCHOR_FALLBACK_SCHEMA,
            "primary_solver": _MC_SOURCE,
            "fallback_solver": _FALLBACK_SOURCE,
            "scope": "low_e_over_n_anchor_fallback",
            "selection_rule": (
                "after_bounded_mc_attempts_two_term_for_unresolved_anchors_at_or_"
                "below_ceiling_else_raw_monte_carlo"
            ),
            "fallback_coordinate": "E_over_N_Td",
            "maximum_fallback_e_over_n_Td": cutoff,
            "component_mixing_within_anchor": False,
            "whole_closure_replacement": False,
            "high_e_over_n_monte_carlo_preserved": True,
            "postprocess_repair": False,
            "decision_relationship": _decision_relationship(decision),
            "bounded_attempt_evidence": {
                "completed_attempt": decision["attempt"],
                "previous_decision_sha256": file_sha256(previous_decision),
                "previous_action": read_selection(previous_decision)["action"],
            },
            "anchors": anchors,
        },
        "quality_scope": decision["quality_scope"],
        "decision": {
            "file": decision_file.name,
            "sha256": file_sha256(decision_file),
            "attempt": decision["attempt"],
            "action": decision["action"],
            "reason": decision.get("reason"),
        },
        "inputs": {
            _MC_SOURCE: _portable_input_binding(mc_dir, decision, _MC_SOURCE),
            _FALLBACK_SOURCE: _portable_input_binding(
                fallback_dir, decision, _FALLBACK_SOURCE
            ),
        },
        "physical_context": mc_manifest["physical_context"],
        "mixture": mc_manifest["mixture"],
        "cross_sections_sha256": mc_manifest["hashes"]["cross_sections_sha256"],
    }
    write_json(output, payload)
    return AnchorFallbackSummary(
        output,
        mc_count,
        fallback_count,
        retained_unqualified_count,
    )


def _require_bounded_fallback_decision(decision: dict[str, Any]) -> Path:
    if (
        decision.get("status") != "selected"
        or decision.get("action") != "select_two_term"
        or decision.get("selected_solver") != _FALLBACK_SOURCE
        or decision.get("selection_scope") != "whole_comsol_closure"
        or decision.get("quality_scope") != _RESTRICTED_LMEA_PROFILE
    ):
        raise AnchorFallbackError(
            "anchor fallback requires a terminal restricted-LMEA MC decision "
            "with a qualified two_term fallback"
        )
    attempt = decision.get("attempt")
    previous = decision.get("previous_decision")
    if (
        isinstance(attempt, bool)
        or not isinstance(attempt, int)
        or attempt < 2
        or not isinstance(previous, dict)
        or not isinstance(previous.get("path"), str)
        or not isinstance(previous.get("sha256"), str)
    ):
        raise AnchorFallbackError(
            "low-E/N fallback requires at least one recorded bounded MC follow-up"
        )
    previous_path = Path(previous["path"]).resolve()
    if not previous_path.is_file() or file_sha256(previous_path) != previous["sha256"]:
        raise AnchorFallbackError("previous MC-attempt decision evidence changed")
    previous_decision = read_selection(previous_path)
    if (
        previous_decision.get("status") != "requires_mc_rerun"
        or previous_decision.get("action")
        not in {
            "add_replicas",
            "extend_tail",
            "extend_time",
            "increase_particles",
        }
        or previous_decision.get("attempt") != attempt - 1
        or previous_decision.get("next_sampling_plan")
        != decision.get("current_sampling_plan")
        or previous_decision.get("mc_run_identity") != decision.get("mc_run_identity")
    ):
        raise AnchorFallbackError(
            "previous decision does not prove the executed bounded MC follow-up"
        )
    return previous_path


def _decision_relationship(decision: dict[str, Any]) -> dict[str, str]:
    """Record the narrow policy transition authorized for this composite."""

    return {
        "terminal_scope": str(decision["selection_scope"]),
        "terminal_selection": str(decision["selected_solver"]),
        "status": "superseded_for_low_field_composite",
        "replacement_scope": "low_e_over_n_anchor_fallback",
        "trigger": "explicit_per_anchor_fallback_request",
    }


def _positive_finite_cutoff(value: object) -> float:
    if isinstance(value, bool):
        raise AnchorFallbackError("low-E/N fallback ceiling must be positive")
    try:
        cutoff = float(value)
    except (TypeError, ValueError) as exc:
        raise AnchorFallbackError("low-E/N fallback ceiling must be positive") from exc
    if not math.isfinite(cutoff) or cutoff <= 0.0:
        raise AnchorFallbackError("low-E/N fallback ceiling must be positive")
    return cutoff


def _verify_bound_source(
    directory: Path,
    decision: dict[str, Any],
    *,
    source: str,
) -> dict[str, Any]:
    manifest_path = directory / "manifest.json"
    manifest = read_json_object(manifest_path)
    binding = decision.get("inputs", {}).get(source)
    if not isinstance(binding, dict):
        raise AnchorFallbackError(f"decision lacks {source} input binding")
    if manifest.get("source") != source:
        raise AnchorFallbackError(
            f"anchor-fallback input source must be {source}: {directory}"
        )
    if binding.get("manifest_sha256") != file_sha256(manifest_path):
        raise AnchorFallbackError(
            f"{source} table generation changed after the MC decision"
        )
    expected_tables = binding.get("tables")
    if not isinstance(expected_tables, dict) or not expected_tables:
        raise AnchorFallbackError(f"decision lacks {source} table hashes")
    listed_tables = manifest.get("tables")
    if not isinstance(listed_tables, dict):
        raise AnchorFallbackError(f"{source} manifest lacks table inventory")
    for name, expected in expected_tables.items():
        path = contained_file(directory, str(name))
        actual = file_sha256(path)
        entry = listed_tables.get(name)
        if (
            actual != expected
            or not isinstance(entry, dict)
            or entry.get("sha256") != actual
        ):
            raise AnchorFallbackError(f"{source} table hash mismatch: {name}")
    return manifest


def _require_source_contracts(
    mc: dict[str, Any],
    fallback: dict[str, Any],
    decision: dict[str, Any],
) -> None:
    if (
        mc.get("source_policy", {}).get("qualification_profile")
        != _RESTRICTED_LMEA_PROFILE
        or decision.get("fallback", {}).get("qualified_for_required_anchors")
        is not True
        or fallback.get("quality_summary", {}).get("passed") is not True
    ):
        raise AnchorFallbackError(
            "MC restricted-LMEA or deterministic fallback qualification is incomplete"
        )


def _read_bound_quality(
    directory: Path,
    decision: dict[str, Any],
    *,
    source: str,
) -> dict[float, dict[str, str]]:
    binding = decision["inputs"][source]
    quality_name = Path(str(binding.get("quality", ""))).name
    if not quality_name:
        raise AnchorFallbackError(f"decision lacks {source} quality artifact")
    path = contained_file(directory, quality_name)
    if file_sha256(path) != binding.get("quality_sha256"):
        raise AnchorFallbackError(f"{source} quality evidence changed after decision")
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as stream:
            rows = list(csv.DictReader(stream))
    except OSError as exc:
        raise AnchorFallbackError(f"cannot read {source} quality evidence") from exc
    result: dict[float, dict[str, str]] = {}
    for row in rows:
        try:
            field = float(row["E_over_N_Td"])
        except (KeyError, TypeError, ValueError) as exc:
            raise AnchorFallbackError(
                f"{source} quality evidence has invalid E/N"
            ) from exc
        if field in result:
            raise AnchorFallbackError(
                f"{source} quality evidence repeats {field:.17g} Td"
            )
        result[field] = row
    return result


def _planned_anchors(decision: dict[str, Any]) -> list[float]:
    raw = decision.get("current_sampling_plan")
    if not isinstance(raw, list) or not raw:
        raise AnchorFallbackError("decision lacks its MC sampling plan")
    try:
        fields = [float(item["e_over_n_Td"]) for item in raw]
    except (KeyError, TypeError, ValueError) as exc:
        raise AnchorFallbackError("decision MC sampling plan is malformed") from exc
    if len(fields) != len(set(fields)) or fields != sorted(fields):
        raise AnchorFallbackError("decision MC anchors must be unique and increasing")
    return fields


def _decision_anchor_evidence(
    decision: dict[str, Any],
) -> dict[float, dict[str, Any]]:
    raw = decision.get("anchor_evidence")
    if not isinstance(raw, list) or not raw:
        raise AnchorFallbackError("decision lacks anchor evidence")
    result: dict[float, dict[str, Any]] = {}
    for item in raw:
        if not isinstance(item, dict) or not isinstance(item.get("passed"), bool):
            raise AnchorFallbackError("decision anchor evidence is malformed")
        try:
            field = float(item["e_over_n_Td"])
        except (KeyError, TypeError, ValueError) as exc:
            raise AnchorFallbackError("decision anchor E/N is invalid") from exc
        reasons = item.get("failure_reasons")
        if (
            field in result
            or not isinstance(reasons, list)
            or any(not isinstance(reason, str) or not reason for reason in reasons)
            or bool(reasons) == bool(item["passed"])
        ):
            raise AnchorFallbackError("decision anchor evidence is inconsistent")
        result[field] = {
            "passed": item["passed"],
            "failure_reasons": tuple(reasons),
        }
    return result


def _mc_quality_outcome(
    row: dict[str, str], quality_scope: str
) -> tuple[bool, bool, tuple[str, ...]]:
    if quality_scope != _RESTRICTED_LMEA_PROFILE:
        raise AnchorFallbackError("unsupported MC anchor-fallback quality scope")
    try:
        reasons = json.loads(row.get("active_closure_failure_reasons_json", ""))
    except json.JSONDecodeError as exc:
        raise AnchorFallbackError("MC active failure reasons are invalid") from exc
    if not isinstance(reasons, list) or any(
        not isinstance(item, str) for item in reasons
    ):
        raise AnchorFallbackError("MC active failure reasons are invalid")
    passed = _accepted_flag(row.get("active_closure_quality_passed"))
    if passed == bool(reasons):
        raise AnchorFallbackError("MC active quality gate and reasons disagree")
    if "aggregate_quality_passed" not in row:
        raise AnchorFallbackError("MC quality lacks the basic aggregate gate")
    aggregate_passed = _accepted_flag(row["aggregate_quality_passed"])
    return passed, aggregate_passed, tuple(reasons)


def _accepted_flag(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true"}


def _portable_input_binding(
    directory: Path,
    decision: dict[str, Any],
    source: str,
) -> dict[str, Any]:
    binding = decision["inputs"][source]
    quality_name = Path(str(binding["quality"])).name
    return {
        "manifest_file": f"{source}_source_manifest.json",
        "manifest_sha256": file_sha256(directory / "manifest.json"),
        "quality_file": quality_name,
        "quality_sha256": file_sha256(directory / quality_name),
        "tables": dict(sorted(binding["tables"].items())),
    }
