from __future__ import annotations

import csv
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path

import pytest

from swarm_workflow.selection import (
    ANCHOR_FALLBACK_SCHEMA,
    ANCHOR_FALLBACK_STAGE,
    COMPOSITE_QUALITY_COLUMNS,
    COMPOSITE_QUALITY_SCHEMA,
)
from swarm_workflow.comsol.models.gec_icp import (
    GecIcpBundleError,
    validate_gec_icp_bundle,
)
from swarm_workflow.comsol.models.gec_icp.run_contracts import GecIcpWorkflowError
from swarm_workflow.comsol.models.gec_icp.run_mapping import (
    load_gec_icp_run_mapping,
)
from swarm_workflow.tables.contracts import (
    MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA,
)
from test_comsol_gec_icp_bundle import (
    _manifest,
    _write_bundle,
    _write_manifest,
)
from test_comsol_gec_icp_run import (
    _run_mapping_payload,
    _write_run_mapping,
)


def test_gec_icp_accepts_provenance_bound_composite_bundle(
    tmp_path: Path,
) -> None:
    bundle = _write_composite_bundle(tmp_path)

    evidence = validate_gec_icp_bundle(
        bundle,
        "composite",
        expected_transport_definition="anchorwise_flux",
        expected_mc_qualification_profile=(
            MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
        ),
    )

    qualification = evidence["provenance"]["mc_qualification"]
    assert qualification["qualified_monte_carlo_anchors"] == 1
    assert qualification["two_term_fallback_anchors"] == 1
    assert qualification["anchor_sources"] == {
        "1": "two_term",
        "2": "monte_carlo",
    }
    assert evidence["artifact_verification"]["verified_files"] == 9
    assert Path(evidence["table_paths"]["evidence:anchor_fallback_plan.json"]).is_file()


def test_gec_icp_composite_rejects_changed_selection_evidence(
    tmp_path: Path,
) -> None:
    bundle = _write_composite_bundle(tmp_path)
    with (bundle / "mc_solver_decision.json").open("a", encoding="utf-8") as stream:
        stream.write("\n")

    with pytest.raises(GecIcpBundleError, match="evidence SHA-256 mismatch"):
        validate_gec_icp_bundle(bundle, "composite")


def test_gec_icp_composite_rejects_quality_source_switch(
    tmp_path: Path,
) -> None:
    bundle = _write_composite_bundle(tmp_path)
    quality = bundle / "quality.csv"
    with quality.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    rows[0]["effective_source"] = "monte_carlo"
    with quality.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=list(COMPOSITE_QUALITY_COLUMNS),
        )
        writer.writeheader()
        writer.writerows(rows)
    manifest = _manifest(bundle)
    manifest["tables"]["quality.csv"]["sha256"] = _digest(quality)
    _write_manifest(bundle, manifest)

    with pytest.raises(GecIcpBundleError, match="quality evidence disagrees"):
        validate_gec_icp_bundle(bundle, "composite")


def test_run_mapping_loads_composite_plasma_spreadsheet_contract(
    tmp_path: Path,
) -> None:
    payload = _run_mapping_payload("monte_carlo")
    payload["bundle"]["expected_source"] = "composite"
    payload["bundle"]["expected_transport_definition"] = "anchorwise_flux"
    payload["closure"]["function_eedf"]["function_tag"] = "sw_icp_eedf_composite"
    payload["closure"]["function_eedf"]["table"] = "eedf_f0_comsol_2d.csv"
    payload["closure"]["function_eedf"]["interpolation"] = (
        "structured_spreadsheet_linear_projection"
    )

    mapping = load_gec_icp_run_mapping(_write_run_mapping(tmp_path, payload))

    assert mapping.bundle.expected_source == "composite"
    assert mapping.closure.function_eedf.table == "eedf_f0_comsol_2d.csv"
    assert mapping.closure.function_eedf.interpolation == (
        "structured_spreadsheet_linear_projection"
    )


def test_run_mapping_requires_composite_mc_profile(tmp_path: Path) -> None:
    payload = _run_mapping_payload("monte_carlo")
    payload["bundle"]["expected_source"] = "composite"
    del payload["bundle"]["expected_mc_qualification_profile"]

    with pytest.raises(GecIcpWorkflowError, match="composite ICP mapping requires"):
        load_gec_icp_run_mapping(_write_run_mapping(tmp_path, payload))


def _write_composite_bundle(tmp_path: Path) -> Path:
    bundle = _write_bundle(tmp_path, "monte_carlo")
    original = _manifest(bundle)
    _write_manifest(bundle, original)
    original = _manifest(bundle)
    mc_manifest_path = bundle / "monte_carlo_source_manifest.json"
    _write_json(mc_manifest_path, original)
    two_term_manifest = {
        "source": "two_term",
        "physical_context": original["physical_context"],
        "mixture": original["mixture"],
        "hashes": {
            "cross_sections_sha256": original["hashes"]["cross_sections_sha256"]
        },
    }
    two_term_manifest_path = bundle / "two_term_source_manifest.json"
    _write_json(two_term_manifest_path, two_term_manifest)

    profile = MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
    mc_quality_path = bundle / "monte_carlo_qualification.csv"
    _write_csv(
        mc_quality_path,
        (
            "E_over_N_Td",
            "qualification_profile",
            "active_closure_quality_passed",
            "active_closure_failure_reasons_json",
        ),
        (
            (1.0, profile, 0, '["mobility_not_converged"]'),
            (2.0, profile, 1, "[]"),
        ),
    )
    two_term_quality_path = bundle / "two_term_quality.csv"
    _write_csv(
        two_term_quality_path,
        ("E_over_N_Td", "passed"),
        ((1.0, 1), (2.0, 1)),
    )
    decision = {
        "status": "selected",
        "action": "select_two_term",
        "selected_solver": "two_term",
        "selection_scope": "whole_comsol_closure",
        "quality_scope": profile,
        "failed_anchors_Td": [1.0],
        "fallback": {"qualified_for_required_anchors": True},
        "mc_run_identity": {
            "mc_solver_source_sha256": original["hashes"]["mc_solver_source_sha256"]
        },
    }
    decision_path = bundle / "mc_solver_decision.json"
    _write_json(decision_path, decision)

    mc_manifest_hash = _digest(mc_manifest_path)
    two_term_manifest_hash = _digest(two_term_manifest_path)
    anchors = [
        {
            "E_over_N_Td": 1.0,
            "effective_source": "two_term",
            "selection_reason": "mc_anchor_unqualified",
            "mc_failure_reasons": ["mobility_not_converged"],
            "source_manifest_sha256": two_term_manifest_hash,
        },
        {
            "E_over_N_Td": 2.0,
            "effective_source": "monte_carlo",
            "selection_reason": "mc_anchor_qualified",
            "mc_failure_reasons": [],
            "source_manifest_sha256": mc_manifest_hash,
        },
    ]
    plan_composition = {
        "schema": ANCHOR_FALLBACK_SCHEMA,
        "primary_solver": "monte_carlo",
        "fallback_solver": "two_term",
        "scope": "complete_anchor_replacement",
        "selection_rule": "qualified_monte_carlo_else_exact_qualified_two_term",
        "component_mixing_within_anchor": False,
        "postprocess_repair": False,
        "decision_relationship": {
            "terminal_scope": "whole_comsol_closure",
            "terminal_selection": "two_term",
            "status": "superseded_for_composite",
            "replacement_scope": "complete_anchor_replacement",
            "trigger": "explicit_per_anchor_fallback_request",
        },
        "anchors": anchors,
    }
    plan = {
        "format_version": 1,
        "stage": ANCHOR_FALLBACK_STAGE,
        "status": "selected",
        "source": "composite",
        "quality_scope": profile,
        "physical_context": original["physical_context"],
        "mixture": original["mixture"],
        "cross_sections_sha256": original["hashes"]["cross_sections_sha256"],
        "source_composition": plan_composition,
        "decision": {"sha256": _digest(decision_path)},
        "inputs": {
            "monte_carlo": {
                "manifest_sha256": mc_manifest_hash,
                "quality_sha256": _digest(mc_quality_path),
            },
            "two_term": {
                "manifest_sha256": two_term_manifest_hash,
                "quality_sha256": _digest(two_term_quality_path),
            },
        },
    }
    plan_path = bundle / "anchor_fallback_plan.json"
    _write_json(plan_path, plan)
    evidence_roles = {
        "anchor_fallback_plan.json": "anchor_fallback_plan",
        "mc_solver_decision.json": "terminal_mc_policy_decision",
        "monte_carlo_qualification.csv": "monte_carlo_anchor_qualification",
        "monte_carlo_source_manifest.json": "monte_carlo_source_manifest",
        "two_term_quality.csv": "two_term_anchor_qualification",
        "two_term_source_manifest.json": "two_term_source_manifest",
    }
    evidence = {
        name: {"role": role, "sha256": _digest(bundle / name)}
        for name, role in evidence_roles.items()
    }
    composition = {
        **plan_composition,
        "mc_database_sha256": "d" * 64,
        "evidence": evidence,
    }

    quality_path = bundle / "quality.csv"
    _write_csv(
        quality_path,
        COMPOSITE_QUALITY_COLUMNS,
        (
            (
                1.0,
                1.0e-21,
                1,
                "two_term",
                "passed",
                "mc_anchor_unqualified",
                '["mobility_not_converged"]',
                two_term_manifest_hash,
            ),
            (
                2.0,
                2.0e-21,
                1,
                "monte_carlo",
                "active_closure_quality_passed",
                "mc_anchor_qualified",
                "[]",
                mc_manifest_hash,
            ),
        ),
    )
    manifest = deepcopy(original)
    manifest["source"] = "composite"
    manifest["source_policy"] = {
        "source": "composite",
        "field_type": "dc",
        "rf_frequency_Hz": None,
        "postprocess": "none",
        "transport_definition": "anchorwise_flux",
        "primary_solver": "monte_carlo",
        "fallback_solver": "two_term",
        "selection": "qualified_mc_else_exact_qualified_two_term_anchor",
        "qualification_profile": profile,
        "qualification_outputs": ["function_eedf", "reduced_mobility"],
        "component_mixing_within_anchor": False,
        "coefficient_repair": False,
    }
    manifest["hashes"].update(
        {
            "anchor_fallback_plan_sha256": _digest(plan_path),
            "monte_carlo_manifest_sha256": mc_manifest_hash,
            "two_term_manifest_sha256": two_term_manifest_hash,
            "monte_carlo_database_sha256": "d" * 64,
        }
    )
    manifest["source_composition"] = composition
    manifest["evidence"] = evidence
    manifest["mc_qualification"] = {
        "profile": profile,
        "all_planned_anchors": 2,
        "qualified_monte_carlo_anchors": 1,
        "two_term_fallback_anchors": 1,
        "failed_monte_carlo_anchors_Td": [1.0],
        "coefficient_tables_available": True,
        "fallback_applied": True,
    }
    manifest["quality_summary"] = {
        "passed": True,
        "failed_points": 0,
        "total_points": 2,
        "monte_carlo_points": 1,
        "two_term_fallback_points": 1,
    }
    manifest["tables"].pop("mc_qualification.csv")
    manifest["tables"]["quality.csv"] = {
        "argument": "E_over_N_Td",
        "artifact_role": "swarm_quality_audit",
        "columns": list(COMPOSITE_QUALITY_COLUMNS),
        "schema": COMPOSITE_QUALITY_SCHEMA,
        "evidence_kind": "per_anchor_selected_source_qualification",
        "sha256": _digest(quality_path),
    }
    _write_manifest(bundle, manifest)
    (bundle / "mc_qualification.csv").unlink()
    return bundle


def _write_csv(
    path: Path,
    columns: tuple[str, ...],
    rows: tuple[tuple[object, ...], ...],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        writer.writerows(rows)


def _write_json(path: Path, value: object) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()
