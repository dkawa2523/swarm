"""Materialize an MC table with bounded low-E/N two-term fallback.

The composite is an input-layer closure, not a fifth solver. Every selected
E/N anchor comes wholly from Monte Carlo aggregate data or a qualified
two-term table row.  Only unresolved anchors in the plan's explicit low-field
band may use two-term; higher-field MC is retained with its real quality state.
No coefficient repair, smoothing, or within-anchor mixing is performed here.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
import shutil
import sqlite3
import tempfile
from typing import Any, Iterable

from swarm_workflow._io import write_csv, write_json
from swarm_workflow.campaign.quality import (
    resolve_evaluation_quality_thresholds,
    resolve_quality_thresholds,
    source_quality_thresholds_json,
)
from swarm_workflow.campaign.repository import read_metadata
from swarm_workflow.campaign.store import provenance_hash_manifest
from swarm_workflow.selection import (
    ANCHOR_FALLBACK_SCHEMA,
    ANCHOR_FALLBACK_SOURCE,
    ANCHOR_FALLBACK_STAGE,
    COMPOSITE_QUALITY_COLUMNS,
    COMPOSITE_QUALITY_SCHEMA,
    AnchorFallbackError,
    file_sha256,
    read_json_object,
    read_selection,
)
from swarm_workflow.quality.policy import (
    quality_policy_provenance,
    quality_thresholds_payload,
)

from . import builder as _builder
from . import contracts as _contracts
from . import monte_carlo as _mc_tables
from . import repository as _repository
from .math import finite_range, strictly_monotonic


_PLAN_FILE = "anchor_fallback_plan.json"
_DECISION_FILE = "mc_solver_decision.json"
_MC_MANIFEST_FILE = "monte_carlo_source_manifest.json"
_FALLBACK_MANIFEST_FILE = "two_term_source_manifest.json"
_MC_QUALITY_FILE = "monte_carlo_qualification.csv"
_FALLBACK_QUALITY_FILE = "two_term_quality.csv"


@dataclass(frozen=True, slots=True)
class CompositeTableSummary:
    database_path: Path
    output_directory: Path
    anchors: int
    monte_carlo_anchors: int
    two_term_anchors: int
    retained_unqualified_monte_carlo_anchors: int


def compose_mc_anchor_fallback_tables(
    mc_database_path: str | Path,
    mc_table_directory: str | Path,
    two_term_table_directory: str | Path,
    plan_path: str | Path,
    decision_path: str | Path,
    output_directory: str | Path,
) -> CompositeTableSummary:
    """Create one model-independent composite table directory.

    The MC database is re-assessed against its recorded restricted-LMEA gates
    before any coefficient is read.
    """

    database = Path(mc_database_path).resolve()
    mc_dir = Path(mc_table_directory).resolve()
    fallback_dir = Path(two_term_table_directory).resolve()
    plan_file = Path(plan_path).resolve()
    decision_file = Path(decision_path).resolve()
    output = Path(output_directory).resolve()
    if not database.is_file():
        raise AnchorFallbackError(f"MC database does not exist: {database}")
    if (
        output in {mc_dir, fallback_dir}
        or mc_dir.is_relative_to(output)
        or fallback_dir.is_relative_to(output)
    ):
        raise AnchorFallbackError(
            "composite output must be separate from source tables"
        )

    plan = _validate_plan(
        plan_file,
        decision_file=decision_file,
        mc_directory=mc_dir,
        fallback_directory=fallback_dir,
    )
    selected = _selected_anchors(plan)
    mc_fields = {
        field
        for field, row in selected.items()
        if row["effective_source"] == _contracts.MC_SOLVER
    }
    fallback_fields = set(selected).difference(mc_fields)

    connection = sqlite3.connect(database)
    try:
        connection.row_factory = sqlite3.Row
        metadata = read_metadata(connection)
        if not _repository._aggregate_rows_are_current(connection):
            raise AnchorFallbackError(
                "MC aggregate rows are stale; rebuild evidence before composition"
            )
        sampling_plan = _repository._validate_mc_sampling_plan_against_cases(
            connection, metadata
        )
        mc_manifest = read_json_object(mc_dir / "manifest.json")
        _validate_database_provenance(metadata, mc_manifest, database=database)
        source_thresholds = resolve_quality_thresholds(connection)
        thresholds = resolve_evaluation_quality_thresholds(connection)
        composite_quality_policy = quality_policy_provenance(
            source_thresholds,
            thresholds,
            source_json=source_quality_thresholds_json(connection),
        )
        evidence = _mc_tables._assess_monte_carlo(
            connection,
            int(plan["mixture"]["mixture_id"]),
            quality_thresholds=thresholds,
            qualification_profile=str(plan["quality_scope"]),
        )
        _validate_reassessed_mc_evidence(evidence, selected)
        # The regional policy intentionally keeps high-E/N MC evidence even
        # when a transport/closure gate is unresolved.  Loading only the
        # active-profile eligible subset would silently drop those anchors and
        # recreate the old whole-anchor substitution bug.
        mixture_id = int(plan["mixture"]["mixture_id"])
        mc_cases = _repository._load_cases(
            connection,
            _contracts.MC_SOLVER,
            mixture_id,
            allowed_e=mc_fields,
        )
        _mc_tables._require_case_mean_energy_strictly_increasing(mc_cases)
        mc_rates = _repository._load_rates(
            connection,
            _contracts.MC_SOLVER,
            mixture_id,
            mc_cases,
            allowed_e=mc_fields,
        )
        mc_eedf = _repository._load_eedf(
            connection,
            _contracts.MC_SOLVER,
            mixture_id,
            mc_cases,
            allowed_e=mc_fields,
        )
        _mc_tables._derive_net_townsend_from_rates(mc_cases, mc_rates)
    finally:
        connection.close()

    if {float(row["E_over_N_Td"]) for row in mc_cases} != mc_fields:
        raise AnchorFallbackError(
            "re-assessed raw MC coefficient anchors disagree with the fallback plan"
        )

    fallback_manifest = read_json_object(fallback_dir / "manifest.json")
    fallback_cases = _select_csv_rows(
        fallback_dir / "transport_vs_en.csv", fallback_fields
    )
    fallback_rates = _select_fallback_rates(fallback_dir, fallback_fields)
    fallback_eedf = _select_csv_rows(fallback_dir / "eedf.csv", fallback_fields)

    cases = [*mc_cases, *fallback_cases]
    rates = [*mc_rates, *fallback_rates]
    eedf = [*mc_eedf, *fallback_eedf]
    cases.sort(key=lambda row: float(row["E_over_N_Td"]))
    rates.sort(key=_rate_sort_key)
    eedf.sort(
        key=lambda row: (
            float(row["E_over_N_Td"]),
            float(row["electron_energy_eV"]),
        )
    )
    _validate_composite_rows(cases, rates, eedf, selected)

    mean_by_field = {
        float(row["E_over_N_Td"]): float(row["mean_energy_eV"]) for row in cases
    }
    rates_by_mean = [
        {
            "mean_energy_eV": mean_by_field[float(row["E_over_N_Td"])],
            **row,
        }
        for row in rates
    ]
    rates_by_mean.sort(
        key=lambda row: (
            float(row["mean_energy_eV"]),
            *_rate_sort_key(row)[1:],
        )
    )
    quality = _composite_quality_rows(selected)
    failed_selected_fields = sorted(
        float(row["E_over_N_Td"])
        for row in quality
        if int(row["passed"]) != 1
    )
    selected_sources_qualified = not failed_selected_fields

    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{output.name}.staging-", dir=output.parent)
    )
    try:
        tables = _write_composite_tables(
            staging,
            cases=cases,
            rates=rates,
            rates_by_mean=rates_by_mean,
            eedf=eedf,
            quality=quality,
        )
        kernel_entry = _copy_shared_collision_kernels(
            mc_dir,
            fallback_dir,
            staging,
            mc_manifest,
            fallback_manifest,
        )
        if kernel_entry is not None:
            tables[_contracts.COLLISION_RATE_KERNEL_TABLE] = kernel_entry
        evidence_files = _copy_source_evidence(
            staging,
            plan_file=plan_file,
            decision_file=decision_file,
            mc_directory=mc_dir,
            fallback_directory=fallback_dir,
            plan=plan,
        )
        source_composition = dict(plan["source_composition"])
        database_hash = file_sha256(database)
        source_composition["mc_database_sha256"] = database_hash
        source_composition["evidence"] = evidence_files
        source_composition["quality_policies"] = {
            _contracts.MC_SOLVER: {
                "quality_thresholds": quality_thresholds_payload(thresholds),
                "quality_policy": composite_quality_policy,
            },
            _contracts.TWO_TERM_SOLVER: {
                "quality_thresholds": fallback_manifest["quality_thresholds"],
                "quality_policy": fallback_manifest["quality_policy"],
            },
        }
        hashes = {
            "anchor_fallback_plan_sha256": file_sha256(plan_file),
            "workflow_config_sha256": mc_manifest["hashes"]["workflow_config_sha256"],
            "base_config_sha256": mc_manifest["hashes"]["base_config_sha256"],
            "cross_sections_sha256": mc_manifest["hashes"]["cross_sections_sha256"],
            "monte_carlo_database_sha256": database_hash,
            "monte_carlo_manifest_sha256": plan["inputs"][_contracts.MC_SOLVER][
                "manifest_sha256"
            ],
            "two_term_manifest_sha256": plan["inputs"][_contracts.TWO_TERM_SOLVER][
                "manifest_sha256"
            ],
        }
        source_policy = {
            "source": ANCHOR_FALLBACK_SOURCE,
            "primary_solver": _contracts.MC_SOLVER,
            "fallback_solver": _contracts.TWO_TERM_SOLVER,
            "selection": (
                "bounded_low_e_over_n_two_term_fallback_else_raw_monte_carlo"
            ),
            "postprocess": "none",
            "field_type": "dc",
            "rf_field_treatment": "none",
            "rf_amplitude_definition": "none",
            "rf_frequency_Hz": None,
            "transport_definition": "anchorwise_solver_native",
            "qualification_profile": plan["quality_scope"],
            "qualification_outputs": ["function_eedf", "reduced_mobility"],
            "component_mixing_within_anchor": False,
            "whole_closure_replacement": False,
            "coefficient_repair": False,
            "solver_physics": {
                _contracts.MC_SOLVER: mc_manifest.get("source_policy", {}).get(
                    "solver_physics"
                ),
                _contracts.TWO_TERM_SOLVER: fallback_manifest.get(
                    "source_policy", {}
                ).get("solver_physics"),
            },
        }
        write_json(
            staging / "manifest.json",
            {
                "format_version": _contracts.FORMAT_VERSION,
                "stage": "build-tables",
                "status": "ok" if selected_sources_qualified else "evidence_only",
                "source": ANCHOR_FALLBACK_SOURCE,
                "physical_context": plan["physical_context"],
                "hashes": hashes,
                "mc_sampling_plan": sampling_plan,
                "mixture": plan["mixture"],
                "source_policy": source_policy,
                "source_composition": source_composition,
                "quality_thresholds": quality_thresholds_payload(thresholds),
                "quality_policy": composite_quality_policy,
                "valid_ranges": {
                    "E_over_N_Td": finite_range(
                        float(row["E_over_N_Td"]) for row in cases
                    ),
                    "E_over_N_V_m2": finite_range(
                        float(row["E_over_N_V_m2"]) for row in cases
                    ),
                    "mean_energy_eV": finite_range(
                        float(row["mean_energy_eV"]) for row in cases
                    ),
                },
                "units": _manifest_units(tables),
                "table_argument": {
                    "primary": "E_over_N_Td",
                    "secondary": "mean_energy_eV",
                    "eedf": "electron_energy_eV,E_over_N_Td",
                },
                "monotonicity": {
                    "mean_energy_strictly_monotonic": True,
                    "reason": None,
                },
                "quality_summary": {
                    "passed": selected_sources_qualified,
                    "failed_points": len(failed_selected_fields),
                    "total_points": len(quality),
                    "monte_carlo_points": len(mc_fields),
                    "two_term_fallback_points": len(fallback_fields),
                    "retained_unqualified_monte_carlo_points": len(
                        failed_selected_fields
                    ),
                    "failed_selected_E_over_N_Td": failed_selected_fields,
                },
                "mc_qualification": {
                    "profile": plan["quality_scope"],
                    "all_planned_anchors": len(selected),
                    "qualified_monte_carlo_anchors": sum(
                        1
                        for field in mc_fields
                        if selected[field]["mc_active_closure_quality_passed"]
                    ),
                    "two_term_fallback_anchors": len(fallback_fields),
                    "retained_unqualified_monte_carlo_anchors": len(
                        failed_selected_fields
                    ),
                    "failed_monte_carlo_anchors_Td": sorted(
                        field
                        for field, row in selected.items()
                        if not row["mc_active_closure_quality_passed"]
                    ),
                    "coefficient_tables_available": selected_sources_qualified,
                    "fallback_applied": True,
                },
                "solver_qualification": None,
                "target_qualification": None,
                "tables": tables,
                "evidence": evidence_files,
            },
        )
        _builder._replace_directory(staging, output)
    finally:
        if staging.exists():
            shutil.rmtree(staging)

    return CompositeTableSummary(
        database,
        output,
        len(selected),
        len(mc_fields),
        len(fallback_fields),
        len(failed_selected_fields),
    )


def _validate_plan(
    path: Path,
    *,
    decision_file: Path,
    mc_directory: Path,
    fallback_directory: Path,
) -> dict[str, Any]:
    plan = read_json_object(path)
    composition = plan.get("source_composition")
    terminal_decision = read_selection(decision_file)
    if (
        plan.get("format_version") != 1
        or plan.get("stage") != ANCHOR_FALLBACK_STAGE
        or plan.get("status") not in {"selected", "evidence_only"}
        or plan.get("source") != ANCHOR_FALLBACK_SOURCE
        or not isinstance(composition, dict)
        or composition.get("schema") != ANCHOR_FALLBACK_SCHEMA
        or composition.get("primary_solver") != _contracts.MC_SOLVER
        or composition.get("fallback_solver") != _contracts.TWO_TERM_SOLVER
        or composition.get("scope") != "low_e_over_n_anchor_fallback"
        or composition.get("selection_rule")
        != (
            "after_bounded_mc_attempts_two_term_for_unresolved_anchors_at_or_"
            "below_ceiling_else_raw_monte_carlo"
        )
        or composition.get("fallback_coordinate") != "E_over_N_Td"
        or not _valid_fallback_ceiling(
            composition.get("maximum_fallback_e_over_n_Td")
        )
        or composition.get("component_mixing_within_anchor") is not False
        or composition.get("whole_closure_replacement") is not False
        or composition.get("high_e_over_n_monte_carlo_preserved") is not True
        or composition.get("postprocess_repair") is not False
        or composition.get("decision_relationship")
        != {
            "terminal_scope": terminal_decision.get("selection_scope"),
            "terminal_selection": terminal_decision.get("selected_solver"),
            "status": "superseded_for_low_field_composite",
            "replacement_scope": "low_e_over_n_anchor_fallback",
            "trigger": "explicit_per_anchor_fallback_request",
        }
        or not _valid_attempt_evidence(composition, terminal_decision)
    ):
        raise AnchorFallbackError("invalid or unsupported MC anchor-fallback plan")
    decision = plan.get("decision")
    if not isinstance(decision, dict) or decision.get("sha256") != file_sha256(
        decision_file
    ):
        raise AnchorFallbackError("anchor-fallback decision binding changed")
    for source, directory in (
        (_contracts.MC_SOLVER, mc_directory),
        (_contracts.TWO_TERM_SOLVER, fallback_directory),
    ):
        binding = plan.get("inputs", {}).get(source)
        if not isinstance(binding, dict) or binding.get(
            "manifest_sha256"
        ) != file_sha256(directory / "manifest.json"):
            raise AnchorFallbackError(
                f"anchor-fallback {source} manifest binding changed"
            )
        manifest = read_json_object(directory / "manifest.json")
        for name, expected in binding.get("tables", {}).items():
            source_path = (directory / str(name)).resolve()
            if (
                not source_path.is_relative_to(directory)
                or not source_path.is_file()
                or file_sha256(source_path) != expected
                or manifest.get("tables", {}).get(name, {}).get("sha256") != expected
            ):
                raise AnchorFallbackError(
                    f"anchor-fallback {source} source table changed: {name}"
                )
    return plan


def _valid_fallback_ceiling(value: object) -> bool:
    if isinstance(value, bool):
        return False
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(number) and number > 0.0


def _valid_attempt_evidence(
    composition: dict[str, Any], terminal_decision: dict[str, Any]
) -> bool:
    evidence = composition.get("bounded_attempt_evidence")
    previous = terminal_decision.get("previous_decision")
    return bool(
        isinstance(evidence, dict)
        and isinstance(previous, dict)
        and evidence.get("completed_attempt") == terminal_decision.get("attempt")
        and evidence.get("previous_decision_sha256") == previous.get("sha256")
        and evidence.get("previous_action")
        in {
            "add_replicas",
            "extend_tail",
            "extend_time",
            "increase_particles",
        }
    )


def _selected_anchors(plan: dict[str, Any]) -> dict[float, dict[str, Any]]:
    raw = plan["source_composition"].get("anchors")
    if not isinstance(raw, list) or not raw:
        raise AnchorFallbackError("anchor-fallback plan has no anchors")
    result: dict[float, dict[str, Any]] = {}
    fields: list[float] = []
    for row in raw:
        if not isinstance(row, dict):
            raise AnchorFallbackError("anchor-fallback anchor is malformed")
        try:
            field = float(row["E_over_N_Td"])
        except (KeyError, TypeError, ValueError) as exc:
            raise AnchorFallbackError("anchor-fallback E/N is invalid") from exc
        source = row.get("effective_source")
        quality_passed = row.get("mc_active_closure_quality_passed")
        if (
            source not in {_contracts.MC_SOLVER, _contracts.TWO_TERM_SOLVER}
            or field in result
            or not isinstance(quality_passed, bool)
        ):
            raise AnchorFallbackError("anchor-fallback source selection is invalid")
        result[field] = row
        fields.append(field)
    if fields != sorted(fields):
        raise AnchorFallbackError("anchor-fallback anchors must be increasing")
    cutoff = float(plan["source_composition"]["maximum_fallback_e_over_n_Td"])
    for field, row in result.items():
        expected_source = (
            _contracts.TWO_TERM_SOLVER
            if not row["mc_active_closure_quality_passed"] and field <= cutoff
            else _contracts.MC_SOLVER
        )
        if row["effective_source"] != expected_source:
            raise AnchorFallbackError(
                "anchor-fallback source violates the explicit low-E/N band"
            )
    return result


def _validate_database_provenance(
    metadata: dict[str, str],
    manifest: dict[str, Any],
    *,
    database: Path,
) -> None:
    actual = provenance_hash_manifest(metadata)
    actual["source_database_sha256"] = file_sha256(database)
    expected = manifest.get("hashes")
    if not isinstance(expected, dict):
        raise AnchorFallbackError("MC source manifest lacks provenance hashes")
    for name, digest in expected.items():
        if actual.get(name) != digest:
            raise AnchorFallbackError(
                f"MC database provenance differs from source evidence: {name}"
            )


def _validate_reassessed_mc_evidence(
    evidence: _contracts.MonteCarloEvidence,
    selected: dict[float, dict[str, Any]],
) -> None:
    quality = {float(row["E_over_N_Td"]): row for row in evidence.quality}
    if set(quality) != set(selected):
        raise AnchorFallbackError(
            "re-assessed MC quality anchors differ from the anchor plan"
        )
    for field, plan_row in selected.items():
        actual_passed = field in evidence.transport_eligible
        expected_passed = plan_row["mc_active_closure_quality_passed"]
        reasons = json.loads(
            str(quality[field].get("active_closure_failure_reasons_json", "[]"))
        )
        if actual_passed != expected_passed or set(reasons) != set(
            plan_row.get("mc_failure_reasons", ())
        ):
            raise AnchorFallbackError(
                f"re-assessed MC quality differs at {field:.17g} Td"
            )


def _select_csv_rows(path: Path, fields: set[float]) -> list[dict[str, Any]]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as stream:
            rows = [dict(row) for row in csv.DictReader(stream)]
    except OSError as exc:
        raise AnchorFallbackError(f"cannot read fallback table: {path}") from exc
    selected = [row for row in rows if float(row["E_over_N_Td"]) in fields]
    found = {float(row["E_over_N_Td"]) for row in selected}
    if found != fields:
        missing = sorted(fields.difference(found))
        raise AnchorFallbackError(
            f"two_term table lacks exact fallback anchors: {missing}"
        )
    return selected


def _select_fallback_rates(directory: Path, fields: set[float]) -> list[dict[str, Any]]:
    """Join deterministic rates to their separately materialized loss moments."""

    rates = _select_csv_rows(directory / "rates_vs_en.csv", fields)
    losses = _select_csv_rows(directory / "energy_loss.csv", fields)
    loss_by_process: dict[tuple[Any, ...], dict[str, Any]] = {}
    for row in losses:
        key = _rate_identity_key(row)
        if key in loss_by_process:
            raise AnchorFallbackError(
                "two_term energy-loss table repeats a reaction process key"
            )
        loss_by_process[key] = row
    rate_keys = {_rate_identity_key(row) for row in rates}
    if rate_keys != set(loss_by_process):
        raise AnchorFallbackError(
            "two_term rate and energy-loss process inventories disagree"
        )
    for row in rates:
        loss = loss_by_process[_rate_identity_key(row)]
        for name in (
            "energy_loss_eV",
            "energy_loss_rate_coefficient_eV_m3_s",
        ):
            value = loss.get(name)
            try:
                number = float(value)
            except (TypeError, ValueError) as exc:
                raise AnchorFallbackError(
                    f"two_term energy-loss table lacks finite {name}"
                ) from exc
            if not math.isfinite(number):
                raise AnchorFallbackError(
                    f"two_term energy-loss table lacks finite {name}"
                )
            row[name] = value
    return rates


def _rate_identity_key(row: dict[str, Any]) -> tuple[Any, ...]:
    try:
        threshold = float(row["threshold_eV"] or 0.0)
        target_fraction = float(row["target_species_fraction"])
        field = float(row["E_over_N_Td"])
    except (KeyError, TypeError, ValueError) as exc:
        raise AnchorFallbackError(
            "rate or energy-loss table contains an invalid process key"
        ) from exc
    return (
        field,
        str(row.get("species", "")),
        str(row.get("process", "")),
        str(row.get("process_type", "")),
        threshold,
        target_fraction,
    )


def _validate_composite_rows(
    cases: list[dict[str, Any]],
    rates: list[dict[str, Any]],
    eedf: list[dict[str, Any]],
    selected: dict[float, dict[str, Any]],
) -> None:
    fields = [float(row["E_over_N_Td"]) for row in cases]
    if fields != list(selected) or len(fields) != len(set(fields)):
        raise AnchorFallbackError(
            "composite coefficients do not contain exactly one case per planned anchor"
        )
    means = [float(row["mean_energy_eV"]) for row in cases]
    if not strictly_monotonic(means):
        raise AnchorFallbackError(
            "composite mean-energy support is not strictly increasing; "
            "coefficient repair is forbidden"
        )
    rate_fields = {float(row["E_over_N_Td"]) for row in rates}
    eedf_fields = {float(row["E_over_N_Td"]) for row in eedf}
    if rate_fields != set(selected) or eedf_fields != set(selected):
        raise AnchorFallbackError(
            "each composite anchor requires complete rate and EEDF rows from "
            "its selected solver"
        )
    process_sets: dict[float, set[tuple[str, str, str]]] = {}
    for row in rates:
        field = float(row["E_over_N_Td"])
        process_sets.setdefault(field, set()).add(
            (
                str(row["species"]),
                str(row["process"]),
                str(row["process_type"]),
            )
        )
    reference = next(iter(process_sets.values()))
    if any(processes != reference for processes in process_sets.values()):
        raise AnchorFallbackError(
            "selected solver anchors do not share one reaction-channel inventory"
        )
    dataset = _contracts.TableDataset(
        source=ANCHOR_FALLBACK_SOURCE,
        mixture_id=0,
        mixture=[],
        cases=cases,
        rates=rates,
        eedf=eedf,
        rate_evidence=None,
        elastic_energy_loss=None,
        elastic_energy_loss_metadata=None,
        quality=[],
        source_policy={},
    )
    try:
        _builder._validate_positive_values(dataset)
    except _contracts.TableBuildError as exc:
        raise AnchorFallbackError(str(exc)) from exc


def _composite_quality_rows(
    selected: dict[float, dict[str, Any]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for field, item in selected.items():
        source = str(item["effective_source"])
        mc_passed = bool(item["mc_active_closure_quality_passed"])
        selected_passed = source == _contracts.TWO_TERM_SOLVER or mc_passed
        rows.append(
            {
                "E_over_N_Td": field,
                "E_over_N_V_m2": field * _contracts.TD_TO_V_M2,
                "passed": int(selected_passed),
                "effective_source": source,
                "source_quality_gate": (
                    "active_closure_quality_passed"
                    if source == _contracts.MC_SOLVER
                    else "passed"
                ),
                "selection_reason": item["selection_reason"],
                "mc_failure_reasons_json": json.dumps(
                    item.get("mc_failure_reasons", ()),
                    separators=(",", ":"),
                ),
                "source_manifest_sha256": item["source_manifest_sha256"],
            }
        )
    return rows


def _write_composite_tables(
    directory: Path,
    *,
    cases: list[dict[str, Any]],
    rates: list[dict[str, Any]],
    rates_by_mean: list[dict[str, Any]],
    eedf: list[dict[str, Any]],
    quality: list[dict[str, Any]],
) -> dict[str, dict[str, Any]]:
    tables: dict[str, dict[str, Any]] = {}
    tables["mean_energy_vs_en.csv"] = _write_table(
        directory / "mean_energy_vs_en.csv",
        ("E_over_N_Td", "E_over_N_V_m2", "mean_energy_eV"),
        cases,
        argument="E_over_N_Td",
    )
    tables["transport_vs_en.csv"] = _write_table(
        directory / "transport_vs_en.csv",
        _contracts.CASE_COLUMNS,
        cases,
        argument="E_over_N_Td",
    )
    tables["rates_vs_en.csv"] = _write_table(
        directory / "rates_vs_en.csv",
        _contracts.RATE_COLUMNS,
        rates,
        argument="E_over_N_Td",
    )
    tables["townsend_vs_en.csv"] = _write_table(
        directory / "townsend_vs_en.csv",
        ("E_over_N_Td", "E_over_N_V_m2", "effective_townsend_m2"),
        cases,
        argument="E_over_N_Td",
    )
    tables["energy_loss.csv"] = _write_table(
        directory / "energy_loss.csv",
        _contracts.ENERGY_LOSS_COLUMNS,
        rates,
        argument="E_over_N_Td",
    )
    tables["quality.csv"] = _write_table(
        directory / "quality.csv",
        COMPOSITE_QUALITY_COLUMNS,
        quality,
        argument="E_over_N_Td",
    )
    tables["quality.csv"].update(
        {
            "schema": COMPOSITE_QUALITY_SCHEMA,
            "evidence_kind": "per_anchor_selected_source_qualification",
        }
    )
    tables["eedf.csv"] = _write_table(
        directory / "eedf.csv",
        _contracts.EEDF_COLUMNS,
        eedf,
        argument="electron_energy_eV,E_over_N_Td",
    )
    tables["eedf.csv"]["artifact_role"] = (
        "raw_swarm_eedf_audit_and_canonicalization_source"
    )
    tables["transport_vs_mean_energy.csv"] = _write_table(
        directory / "transport_vs_mean_energy.csv",
        (
            "mean_energy_eV",
            "E_over_N_Td",
            "E_over_N_V_m2",
            *_contracts.CASE_COLUMNS[3:],
        ),
        cases,
        argument="mean_energy_eV",
    )
    tables["rates_vs_mean_energy.csv"] = _write_table(
        directory / "rates_vs_mean_energy.csv",
        ("mean_energy_eV", *_contracts.RATE_COLUMNS),
        rates_by_mean,
        argument="mean_energy_eV",
    )
    return tables


def _write_table(
    path: Path,
    columns: tuple[str, ...],
    rows: Iterable[dict[str, Any]],
    *,
    argument: str,
) -> dict[str, Any]:
    write_csv(path, columns, rows)
    return {
        "columns": list(columns),
        "sha256": file_sha256(path),
        "units": {
            column: _contracts.UNITS[column]
            for column in columns
            if column in _contracts.UNITS
        },
        "argument": argument,
    }


def _copy_shared_collision_kernels(
    mc_directory: Path,
    fallback_directory: Path,
    output_directory: Path,
    mc_manifest: dict[str, Any],
    fallback_manifest: dict[str, Any],
) -> dict[str, Any] | None:
    """Reuse one immutable kernel table when both sources used the same data."""

    name = _contracts.COLLISION_RATE_KERNEL_TABLE
    entries = [
        manifest.get("tables", {}).get(name)
        for manifest in (mc_manifest, fallback_manifest)
    ]
    if entries == [None, None]:
        return None
    if any(not isinstance(entry, dict) for entry in entries):
        raise AnchorFallbackError(
            "composite sources disagree on collision-rate kernel availability"
        )
    mc_entry, fallback_entry = entries
    assert isinstance(mc_entry, dict) and isinstance(fallback_entry, dict)
    if mc_entry.get("sha256") != fallback_entry.get("sha256") or mc_entry.get(
        "cross_sections_sha256"
    ) != fallback_entry.get("cross_sections_sha256"):
        raise AnchorFallbackError(
            "composite sources use different collision-rate kernels"
        )
    source = mc_directory / name
    fallback = fallback_directory / name
    digest = file_sha256(source)
    if digest != mc_entry.get("sha256") or file_sha256(fallback) != digest:
        raise AnchorFallbackError("collision-rate kernel hash mismatch")
    target = output_directory / name
    target.write_bytes(source.read_bytes())
    return dict(mc_entry)


def _copy_source_evidence(
    directory: Path,
    *,
    plan_file: Path,
    decision_file: Path,
    mc_directory: Path,
    fallback_directory: Path,
    plan: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    mc_quality_name = str(plan["inputs"][_contracts.MC_SOLVER]["quality_file"])
    fallback_quality_name = str(
        plan["inputs"][_contracts.TWO_TERM_SOLVER]["quality_file"]
    )
    copies = (
        (_PLAN_FILE, plan_file, "anchor_fallback_plan"),
        (_DECISION_FILE, decision_file, "terminal_mc_policy_decision"),
        (
            _MC_MANIFEST_FILE,
            mc_directory / "manifest.json",
            "monte_carlo_source_manifest",
        ),
        (
            _FALLBACK_MANIFEST_FILE,
            fallback_directory / "manifest.json",
            "two_term_source_manifest",
        ),
        (
            _MC_QUALITY_FILE,
            mc_directory / mc_quality_name,
            "monte_carlo_anchor_qualification",
        ),
        (
            _FALLBACK_QUALITY_FILE,
            fallback_directory / fallback_quality_name,
            "two_term_anchor_qualification",
        ),
    )
    result: dict[str, dict[str, Any]] = {}
    for name, source, role in copies:
        target = directory / name
        target.write_bytes(source.read_bytes())
        result[name] = {"role": role, "sha256": file_sha256(target)}
    return result


def _rate_sort_key(row: dict[str, Any]) -> tuple[Any, ...]:
    return (
        float(row["E_over_N_Td"]),
        str(row["species"]),
        str(row["process"]),
        str(row["process_type"]),
        float(row["threshold_eV"] or 0.0),
    )


def _manifest_units(tables: dict[str, dict[str, Any]]) -> dict[str, str]:
    return {
        column: _contracts.UNITS[column]
        for column in sorted(
            {
                column
                for table in tables.values()
                for column in table.get("columns", ())
                if column in _contracts.UNITS
            }
        )
    }
