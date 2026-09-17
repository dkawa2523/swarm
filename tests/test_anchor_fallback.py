from __future__ import annotations

import csv
from hashlib import sha256
import json
from pathlib import Path

import pytest

from swarm_workflow.selection import (
    ANCHOR_FALLBACK_SCHEMA,
    AnchorFallbackError,
    build_mc_anchor_fallback_plan,
)
from swarm_workflow.quality.monte_carlo.policy import (
    MonteCarloConvergencePolicy,
    decide_monte_carlo_closure,
)
from swarm_workflow.tables.anchor_fallback import (
    _composite_quality_rows,
    _select_fallback_rates,
    _validate_reassessed_mc_evidence,
)
from swarm_workflow.tables.contracts import MonteCarloEvidence
from test_mc_policy import _tables


def _write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _fallback_inputs(
    tmp_path: Path,
    *,
    exact_fallback_anchor: bool = True,
    failed_anchors: tuple[float, ...] = (5.0,),
) -> tuple[Path, Path, Path]:
    fallback_anchors = (5.0, 10.0) if exact_fallback_anchor else (1.0, 10.0)
    mc, fallback = _tables(
        tmp_path,
        failures={
            anchor: ["solver_mean_energy_stationarity_relative_ci95_bound"]
            for anchor in failed_anchors
        },
        active_scope=True,
        fallback_anchors=fallback_anchors,
    )
    mc_quality = mc / "quality.csv"
    with mc_quality.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        field = float(row["E_over_N_Td"])
        failed = field in failed_anchors
        row["aggregate_quality_passed"] = "1"
        row["active_closure_quality_passed"] = "0" if failed else "1"
        row["active_closure_failure_reasons_json"] = json.dumps(
            ["solver_mean_energy_stationarity_relative_ci95_bound"] if failed else []
        )
        row["active_closure_failure_axes_json"] = json.dumps(
            ["time_horizon"] if failed else []
        )
    _write_rows(mc_quality, rows)
    mc_manifest_path = mc / "manifest.json"
    mc_manifest = json.loads(mc_manifest_path.read_text(encoding="utf-8"))
    mc_manifest["tables"]["quality.csv"]["sha256"] = sha256(
        mc_quality.read_bytes()
    ).hexdigest()
    mc_manifest_path.write_text(json.dumps(mc_manifest), encoding="utf-8")

    fallback_manifest_path = fallback / "manifest.json"
    fallback_manifest = json.loads(fallback_manifest_path.read_text(encoding="utf-8"))
    fallback_manifest["quality_summary"] = {
        "passed": True,
        "failed_points": 0,
        "total_points": len(fallback_anchors),
    }
    fallback_manifest_path.write_text(json.dumps(fallback_manifest), encoding="utf-8")

    first_decision = tmp_path / "decision-attempt1.json"
    decide_monte_carlo_closure(
        mc,
        fallback,
        attempt=1,
        output=first_decision,
    )
    first = json.loads(first_decision.read_text())
    assert first["action"] == "extend_time"
    mc_manifest = json.loads(mc_manifest_path.read_text(encoding="utf-8"))
    mc_manifest["mc_sampling_plan"]["entries"] = first["next_sampling_plan"]
    mc_manifest_path.write_text(json.dumps(mc_manifest), encoding="utf-8")

    decision = tmp_path / "decision-attempt2.json"
    decide_monte_carlo_closure(
        mc,
        fallback,
        attempt=2,
        output=decision,
        previous_decision=first_decision,
    )
    assert json.loads(decision.read_text())["action"] == "select_two_term"
    return mc, fallback, decision


def test_anchor_fallback_plan_keeps_passed_mc_and_replaces_only_failure(
    tmp_path: Path,
) -> None:
    mc, fallback, decision = _fallback_inputs(tmp_path)
    output = tmp_path / "anchor-fallback.json"

    summary = build_mc_anchor_fallback_plan(
        mc,
        fallback,
        decision,
        output,
        maximum_fallback_e_over_n_Td=5.0,
    )

    assert summary.monte_carlo_anchors == 1
    assert summary.two_term_anchors == 1
    assert summary.retained_unqualified_monte_carlo_anchors == 0
    result = json.loads(output.read_text(encoding="utf-8"))
    assert result["source"] == "composite"
    assert result["source_composition"]["schema"] == ANCHOR_FALLBACK_SCHEMA
    assert result["source_composition"]["component_mixing_within_anchor"] is False
    assert result["source_composition"]["decision_relationship"] == {
        "terminal_scope": "whole_comsol_closure",
        "terminal_selection": "two_term",
        "status": "superseded_for_low_field_composite",
        "replacement_scope": "low_e_over_n_anchor_fallback",
        "trigger": "explicit_per_anchor_fallback_request",
    }
    assert [
        (row["E_over_N_Td"], row["effective_source"])
        for row in result["source_composition"]["anchors"]
    ] == [(5.0, "two_term"), (10.0, "monte_carlo")]
    assert result["source_composition"]["anchors"][0]["mc_failure_reasons"] == [
        "solver_mean_energy_stationarity_relative_ci95_bound"
    ]


def test_anchor_fallback_requires_exact_qualified_two_term_anchor(
    tmp_path: Path,
) -> None:
    mc, fallback, decision = _fallback_inputs(tmp_path, exact_fallback_anchor=False)

    with pytest.raises(AnchorFallbackError, match="exact qualified anchor"):
        build_mc_anchor_fallback_plan(
            mc,
            fallback,
            decision,
            tmp_path / "anchor-fallback.json",
            maximum_fallback_e_over_n_Td=5.0,
        )


def test_anchor_fallback_rejects_source_change_after_decision(
    tmp_path: Path,
) -> None:
    mc, fallback, decision = _fallback_inputs(tmp_path)
    quality = fallback / "quality.csv"
    quality.write_bytes(quality.read_bytes() + b"\n")

    with pytest.raises(AnchorFallbackError, match="table hash mismatch"):
        build_mc_anchor_fallback_plan(
            mc,
            fallback,
            decision,
            tmp_path / "anchor-fallback.json",
            maximum_fallback_e_over_n_Td=5.0,
        )


def test_anchor_fallback_retains_unqualified_high_field_monte_carlo(
    tmp_path: Path,
) -> None:
    mc, fallback, decision = _fallback_inputs(
        tmp_path,
        failed_anchors=(5.0, 10.0),
    )
    output = tmp_path / "anchor-fallback.json"

    summary = build_mc_anchor_fallback_plan(
        mc,
        fallback,
        decision,
        output,
        maximum_fallback_e_over_n_Td=5.0,
    )

    result = json.loads(output.read_text(encoding="utf-8"))
    assert result["status"] == "evidence_only"
    assert summary.retained_unqualified_monte_carlo_anchors == 1
    assert [
        (row["E_over_N_Td"], row["effective_source"], row["selection_reason"])
        for row in result["source_composition"]["anchors"]
    ] == [
        (
            5.0,
            "two_term",
            "bounded_mc_unresolved_in_low_e_over_n_fallback_band",
        ),
        (
            10.0,
            "monte_carlo",
            "mc_anchor_unqualified_retained_above_fallback_band",
        ),
    ]


def test_anchor_fallback_rejects_immediate_first_attempt_substitution(
    tmp_path: Path,
) -> None:
    mc, fallback = _tables(
        tmp_path,
        failures={5.0: ["solver_mean_energy_stationarity_relative_ci95_bound"]},
        active_scope=True,
    )
    quality = mc / "quality.csv"
    with quality.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    rows[0]["active_closure_quality_passed"] = "0"
    rows[0]["active_closure_failure_reasons_json"] = json.dumps(
        ["solver_mean_energy_stationarity_relative_ci95_bound"]
    )
    rows[0]["active_closure_failure_axes_json"] = json.dumps(["time_horizon"])
    for row in rows:
        row["aggregate_quality_passed"] = "1"
    _write_rows(quality, rows)
    mc_manifest_path = mc / "manifest.json"
    manifest = json.loads(mc_manifest_path.read_text(encoding="utf-8"))
    manifest["tables"]["quality.csv"]["sha256"] = sha256(
        quality.read_bytes()
    ).hexdigest()
    mc_manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    fallback_manifest_path = fallback / "manifest.json"
    fallback_manifest = json.loads(fallback_manifest_path.read_text())
    fallback_manifest["quality_summary"] = {"passed": True}
    fallback_manifest_path.write_text(json.dumps(fallback_manifest), encoding="utf-8")
    decision = tmp_path / "decision.json"
    decide_monte_carlo_closure(
        mc,
        fallback,
        attempt=1,
        output=decision,
        policy=MonteCarloConvergencePolicy(
            maximum_particle_barriers_per_replica=30_000,
            maximum_total_particle_barriers=100_000,
        ),
    )

    with pytest.raises(AnchorFallbackError, match="bounded MC follow-up"):
        build_mc_anchor_fallback_plan(
            mc,
            fallback,
            decision,
            tmp_path / "anchor-fallback.json",
            maximum_fallback_e_over_n_Td=5.0,
        )


def test_composite_quality_does_not_mark_retained_high_field_mc_as_qualified() -> None:
    selected = {
        5.0: {
            "E_over_N_Td": 5.0,
            "effective_source": "two_term",
            "selection_reason": (
                "bounded_mc_unresolved_in_low_e_over_n_fallback_band"
            ),
            "mc_active_closure_quality_passed": False,
            "mc_failure_reasons": ["stationarity"],
            "source_manifest_sha256": "a" * 64,
        },
        10.0: {
            "E_over_N_Td": 10.0,
            "effective_source": "monte_carlo",
            "selection_reason": (
                "mc_anchor_unqualified_retained_above_fallback_band"
            ),
            "mc_active_closure_quality_passed": False,
            "mc_failure_reasons": ["lineage"],
            "source_manifest_sha256": "b" * 64,
        },
    }
    evidence = MonteCarloEvidence(
        quality=[
            {
                "E_over_N_Td": 5.0,
                "active_closure_failure_reasons_json": '["stationarity"]',
            },
            {
                "E_over_N_Td": 10.0,
                "active_closure_failure_reasons_json": '["lineage"]',
            },
        ],
        estimator_schema="test",
        aggregate_eligible={5.0, 10.0},
        full_transport_eligible=set(),
        active_failure_reasons={5.0: ["stationarity"], 10.0: ["lineage"]},
        rate_failures={},
        transport_eligible=set(),
        eedf_rate_eligible=set(),
    )

    _validate_reassessed_mc_evidence(evidence, selected)
    rows = _composite_quality_rows(selected)

    assert [(row["effective_source"], row["passed"]) for row in rows] == [
        ("two_term", 1),
        ("monte_carlo", 0),
    ]


def test_fallback_rates_join_complete_energy_loss_by_process_key(
    tmp_path: Path,
) -> None:
    rate = {
        "E_over_N_Td": "5",
        "species": "Ar",
        "process": "e+Ar=>e+Ars",
        "process_type": "excitation",
        "threshold_eV": "11.5",
        "target_species_fraction": "1",
        "rate_coefficient_m3_s": "2e-15",
    }
    loss = {
        **rate,
        "energy_loss_eV": "11.5",
        "energy_loss_rate_coefficient_eV_m3_s": "2.3e-14",
    }
    _write_rows(tmp_path / "rates_vs_en.csv", [rate])
    _write_rows(tmp_path / "energy_loss.csv", [loss])

    result = _select_fallback_rates(tmp_path, {5.0})

    assert result[0]["energy_loss_eV"] == "11.5"
    assert result[0]["energy_loss_rate_coefficient_eV_m3_s"] == "2.3e-14"


def test_fallback_rates_reject_energy_loss_process_mismatch(
    tmp_path: Path,
) -> None:
    rate = {
        "E_over_N_Td": "5",
        "species": "Ar",
        "process": "excitation",
        "process_type": "excitation",
        "threshold_eV": "11.5",
        "target_species_fraction": "1",
    }
    loss = {
        **rate,
        "process": "different channel",
        "energy_loss_eV": "11.5",
        "energy_loss_rate_coefficient_eV_m3_s": "2.3e-14",
    }
    _write_rows(tmp_path / "rates_vs_en.csv", [rate])
    _write_rows(tmp_path / "energy_loss.csv", [loss])

    with pytest.raises(AnchorFallbackError, match="inventories disagree"):
        _select_fallback_rates(tmp_path, {5.0})
