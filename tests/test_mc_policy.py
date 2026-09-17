from __future__ import annotations

import csv
from hashlib import sha256
import json
from pathlib import Path
import pytest

from electron_swarm import load_config
from electron_swarm.solvers.monte_carlo.source_identity import (
    monte_carlo_source_sha256,
)
from swarm_workflow.selection import physical_context

from swarm_workflow.cli import main as workflow_cli_main
from swarm_workflow.quality.monte_carlo.policy import (
    FailureAxis,
    MonteCarloConvergencePolicy,
    MonteCarloPolicyError,
    decide_monte_carlo_closure,
    failure_axes_for_reasons,
)


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _tables(
    tmp_path: Path,
    *,
    failures: dict[float, list[str]] | None = None,
    failure_axes: dict[float, list[str]] | None = None,
    active_scope: bool = False,
    fallback_anchors: tuple[float, ...] = (5.0, 10.0),
) -> tuple[Path, Path]:
    failures = failures or {}
    failure_axes = failure_axes or {}
    mc = tmp_path / "mc"
    two_term = tmp_path / "two_term"
    mc.mkdir(parents=True)
    two_term.mkdir(parents=True)
    plan = [
        {
            "e_over_n_Td": anchor,
            "particles": 128,
            "warmup_collisions": 16,
            "max_collisions": 32,
            "tail_max_collisions": 32,
            "replicas": 4,
            "transport_correlation_lag_barriers": 4,
            "transport_estimator": "single_field",
        }
        for anchor in (5.0, 10.0)
    ]
    mc_manifest = {
        "source": "monte_carlo",
        "mc_sampling_plan": {"entries": plan, "base_seed": 1729},
        "source_policy": {
            "qualification_profile": "function_eedf_restricted_lmea"
            if active_scope
            else "full_transport"
        },
    }
    (mc / "manifest.json").write_text(json.dumps(mc_manifest), encoding="utf-8")
    (two_term / "manifest.json").write_text(
        json.dumps({"source": "two_term"}), encoding="utf-8"
    )
    mc_rows: list[dict[str, object]] = []
    for anchor in (5.0, 10.0):
        reasons = failures.get(anchor, [])
        axes = failure_axes.get(
            anchor,
            [axis.value for axis in failure_axes_for_reasons(reasons)],
        )
        row: dict[str, object] = {
            "E_over_N_Td": anchor,
            "passed": int(not reasons),
            "failure_reasons_json": json.dumps(reasons),
            "failure_axes_json": json.dumps(axes),
            "valid_replicates": 4,
        }
        if active_scope:
            row.update(
                {
                    "active_closure_quality_passed": 1,
                    "active_closure_failure_reasons_json": "[]",
                    "active_closure_failure_axes_json": "[]",
                }
            )
        mc_rows.append(row)
    _write_csv(mc / "quality.csv", mc_rows)
    _write_csv(
        two_term / "quality.csv",
        [
            {
                "E_over_N_Td": anchor,
                "passed": 1,
                "failure_reasons_json": "[]",
            }
            for anchor in fallback_anchors
        ],
    )
    context = physical_context(
        load_config(
            Path(__file__).resolve().parents[1] / "examples/argon_gec_ccp_base.yaml"
        )
    )
    for directory in (mc, two_term):
        path = directory / "manifest.json"
        manifest = json.loads(path.read_text())
        manifest.update(
            {
                "physical_context": context,
                "mixture": {
                    "mixture_id": 0,
                    "species": [{"species": "Ar", "fraction": 1.0, "mass_amu": 39.948}],
                },
                "hashes": {
                    "base_config_sha256": "a" * 64,
                    "cross_sections_sha256": "b" * 64,
                    "mc_transport_estimator_schema_version": "test_estimator.v1",
                    "mc_eedf_estimator_schema_version": "test_eedf.v1",
                    **(
                        {"mc_solver_source_sha256": (monte_carlo_source_sha256())}
                        if directory == mc
                        else {}
                    ),
                },
                "tables": {
                    "quality.csv": {
                        "sha256": sha256(
                            (directory / "quality.csv").read_bytes()
                        ).hexdigest()
                    }
                },
            }
        )
        path.write_text(json.dumps(manifest), encoding="utf-8")
    return mc, two_term


def _decision(
    tmp_path: Path,
    mc: Path,
    two_term: Path,
    *,
    attempt: int,
) -> dict[str, object]:
    output = tmp_path / "decision.json"
    previous = None
    if attempt == 2:
        previous = tmp_path / "first-decision.json"
        decide_monte_carlo_closure(mc, two_term, attempt=1, output=previous)
        first = json.loads(previous.read_text())
        next_mc = tmp_path / "mc-attempt2"
        next_mc.mkdir()
        for path in mc.iterdir():
            (next_mc / path.name).write_bytes(path.read_bytes())
        mc = next_mc
        manifest = json.loads((mc / "manifest.json").read_text())
        manifest["mc_sampling_plan"]["entries"] = first["next_sampling_plan"]
        with (mc / "quality.csv").open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        for row, entry in zip(rows, first["next_sampling_plan"], strict=True):
            row["valid_replicates"] = entry["replicas"]
        _write_csv(mc / "quality.csv", rows)
        manifest["tables"]["quality.csv"]["sha256"] = sha256(
            (mc / "quality.csv").read_bytes()
        ).hexdigest()
        (mc / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    decide_monte_carlo_closure(
        mc,
        two_term,
        attempt=attempt,
        output=output,
        previous_decision=previous,
    )
    return json.loads(output.read_text(encoding="utf-8"))


def test_policy_accepts_mc_only_when_complete_active_closure_passes(
    tmp_path: Path,
) -> None:
    mc, two_term = _tables(
        tmp_path,
        failures={5.0: ["mc_weighted_growth_lag_not_converged:diffusion_L_m2_s"]},
        active_scope=True,
    )

    decision = _decision(tmp_path, mc, two_term, attempt=1)

    assert decision["quality_scope"] == "function_eedf_restricted_lmea"
    assert decision["action"] == "accept_monte_carlo"
    assert decision["selected_solver"] == "monte_carlo"
    assert decision["selection_scope"] == "whole_comsol_closure"


def test_stationarity_failure_extends_time_once_without_adding_replicas(
    tmp_path: Path,
) -> None:
    mc, two_term = _tables(
        tmp_path,
        failures={
            5.0: ["mc_weighted_growth_stationarity_not_converged:mobility_m2_V_s"]
        },
    )

    decision = _decision(tmp_path, mc, two_term, attempt=1)
    current = decision["current_sampling_plan"][0]
    updated = decision["next_sampling_plan"][0]

    assert decision["action"] == "extend_time"
    assert decision["selected_solver"] is None
    assert updated["particles"] == current["particles"]
    assert updated["replicas"] == current["replicas"] == 4
    assert updated["warmup_collisions"] == 4 * current["warmup_collisions"]
    assert updated["max_collisions"] == 4 * current["max_collisions"]
    assert decision["next_sampling_plan"][1] == decision["current_sampling_plan"][1]
    assert decision["reuse_previous_mc_database"] is True


def test_stationarity_extension_does_not_cross_compute_ceiling(
    tmp_path: Path,
) -> None:
    mc, two_term = _tables(
        tmp_path,
        failures={
            5.0: ["mc_weighted_growth_stationarity_not_converged:mobility_m2_V_s"]
        },
    )
    output = tmp_path / "bounded-decision.json"
    policy = MonteCarloConvergencePolicy(
        maximum_particle_barriers_per_replica=30_000,
        maximum_total_particle_barriers=100_000,
    )

    decide_monte_carlo_closure(
        mc,
        two_term,
        attempt=1,
        output=output,
        policy=policy,
    )
    decision = json.loads(output.read_text(encoding="utf-8"))

    assert decision["action"] == "select_two_term"
    assert decision["reason"] == "mc_budget_exhausted"
    assert decision["next_sampling_plan"] is None


def test_precision_only_failure_adds_four_replicas_and_reuses_cases(
    tmp_path: Path,
) -> None:
    mc, two_term = _tables(tmp_path, failures={5.0: ["mobility_rse"]})

    decision = _decision(tmp_path, mc, two_term, attempt=1)

    assert decision["action"] == "add_replicas"
    assert decision["reuse_previous_mc_database"] is True
    assert decision["next_sampling_plan"][0]["replicas"] == 8
    assert decision["next_sampling_plan"][0]["max_collisions"] == 32


def test_lineage_failure_increases_only_the_failed_anchor_population(
    tmp_path: Path,
) -> None:
    mc, two_term = _tables(
        tmp_path,
        failures={5.0: ["mc_weighted_growth_lineage_degenerate"]},
    )

    decision = _decision(tmp_path, mc, two_term, attempt=1)

    assert decision["action"] == "increase_particles"
    assert decision["anchor_evidence"][0]["failure_axes"] == [
        FailureAxis.POPULATION_LINEAGE.value
    ]
    assert decision["next_sampling_plan"][0]["particles"] == 256
    assert decision["next_sampling_plan"][1] == decision["current_sampling_plan"][1]


def test_censored_rate_failure_extends_only_the_tail_budget(tmp_path: Path) -> None:
    mc, two_term = _tables(
        tmp_path,
        failures={5.0: ["required_rate_censored_all_zero:Ar:ionization:ionization"]},
    )

    decision = _decision(tmp_path, mc, two_term, attempt=1)

    assert decision["action"] == "extend_tail"
    current = decision["current_sampling_plan"][0]
    updated = decision["next_sampling_plan"][0]
    assert updated["tail_max_collisions"] == 4 * current["tail_max_collisions"]
    assert updated["warmup_collisions"] == current["warmup_collisions"]
    assert updated["max_collisions"] == current["max_collisions"]


def test_decision_routes_from_typed_axis_not_failure_text(tmp_path: Path) -> None:
    mc, two_term = _tables(
        tmp_path,
        failures={5.0: ["new_quality_gate_without_name_convention"]},
        failure_axes={5.0: [FailureAxis.RARE_TAIL.value]},
    )

    decision = _decision(tmp_path, mc, two_term, attempt=1)

    assert decision["action"] == "extend_tail"


def test_policy_tries_another_matching_action_before_fallback(tmp_path: Path) -> None:
    mc, two_term = _tables(
        tmp_path,
        failures={
            5.0: [
                "solver_mean_energy_stationarity_relative_ci95_bound",
                "mobility_rse",
            ]
        },
    )
    output = tmp_path / "decision.json"

    decide_monte_carlo_closure(
        mc,
        two_term,
        attempt=1,
        output=output,
        policy=MonteCarloConvergencePolicy(
            maximum_particle_barriers_per_replica=30_000,
            maximum_total_particle_barriers=130_000,
        ),
    )
    decision = json.loads(output.read_text(encoding="utf-8"))

    assert decision["action"] == "add_replicas"


def test_final_failed_attempt_selects_one_two_term_closure(tmp_path: Path) -> None:
    mc, two_term = _tables(
        tmp_path,
        failures={5.0: ["solver_mean_energy_stationarity_relative_ci95_bound"]},
    )

    decision = _decision(tmp_path, mc, two_term, attempt=2)

    assert decision["action"] == "select_two_term"
    assert decision["selected_solver"] == "two_term"
    assert decision["reason"] == "mc_budget_exhausted"
    assert decision["failed_anchors_Td"] == [5.0]
    assert decision["next_sampling_plan"] is None
    assert decision["mc_results_retained_as_validation"] is True


def test_policy_does_not_fallback_for_unclassified_physics_failure(
    tmp_path: Path,
) -> None:
    mc, two_term = _tables(tmp_path, failures={5.0: ["negative_eedf"]})

    decision = _decision(tmp_path, mc, two_term, attempt=1)

    assert decision["action"] == "blocked"
    assert decision["selected_solver"] is None
    assert decision["reason"] == "mc_failure_requires_model_or_evidence_correction"


def test_policy_blocks_when_two_term_does_not_cover_required_anchor(
    tmp_path: Path,
) -> None:
    mc, two_term = _tables(
        tmp_path,
        failures={5.0: ["mobility_rse"]},
        fallback_anchors=(10.0,),
    )

    decision = _decision(tmp_path, mc, two_term, attempt=2)

    assert decision["action"] == "blocked"
    assert decision["reason"] == "two_term_fallback_not_qualified"
    assert decision["fallback"]["failures"]["5"] == [
        "required_anchor_outside_two_term_support"
    ]


def test_two_term_bracketing_support_does_not_require_identical_anchor_grid(
    tmp_path: Path,
) -> None:
    mc, two_term = _tables(
        tmp_path,
        failures={5.0: ["mobility_rse"]},
        fallback_anchors=(1.0, 7.5, 20.0),
    )

    decision = _decision(tmp_path, mc, two_term, attempt=2)

    assert decision["action"] == "select_two_term"
    assert decision["fallback"]["qualified_for_required_anchors"] is True


def test_decision_records_immutable_input_hashes(tmp_path: Path) -> None:
    mc, two_term = _tables(tmp_path)

    decision = _decision(tmp_path, mc, two_term, attempt=1)

    quality = mc / "quality.csv"
    assert (
        decision["inputs"]["monte_carlo"]["quality_sha256"]
        == sha256(quality.read_bytes()).hexdigest()
    )


def test_decide_mc_cli_writes_selection_artifact(tmp_path: Path) -> None:
    mc, two_term = _tables(tmp_path)
    output = tmp_path / "cli-decision.json"

    workflow_cli_main(
        [
            "decide-mc",
            "--mc-tables",
            str(mc),
            "--two-term-tables",
            str(two_term),
            "--attempt",
            "1",
            "--output",
            str(output),
        ]
    )

    assert json.loads(output.read_text(encoding="utf-8"))["action"] == (
        "accept_monte_carlo"
    )


@pytest.mark.parametrize("change", ["mixture", "cross_sections", "gas_temperature"])
def test_selection_rejects_different_physics_before_fallback(
    tmp_path: Path, change: str
) -> None:
    mc, fallback = _tables(tmp_path, failures={5.0: ["mobility_rse"]})
    path = fallback / "manifest.json"
    manifest = json.loads(path.read_text())
    if change == "mixture":
        manifest["mixture"]["species"][0]["species"] = "Ne"
    elif change == "cross_sections":
        manifest["hashes"]["cross_sections_sha256"] = "c" * 64
    else:
        manifest["physical_context"]["gas_temperature_K"] = 600
    path.write_text(json.dumps(manifest))
    with pytest.raises(MonteCarloPolicyError, match="different"):
        _decision(tmp_path, mc, fallback, attempt=1)


@pytest.mark.parametrize("replicas", [1, 16])
def test_first_attempt_replica_bounds_are_enforced_at_decision(
    tmp_path: Path, replicas: int
) -> None:
    mc, fallback = _tables(tmp_path)
    path = mc / "manifest.json"
    manifest = json.loads(path.read_text())
    manifest["mc_sampling_plan"]["entries"][0]["replicas"] = replicas
    path.write_text(json.dumps(manifest))
    with pytest.raises(MonteCarloPolicyError, match="replica"):
        _decision(tmp_path, mc, fallback, attempt=1)


def test_first_attempt_accepts_per_anchor_replicas_within_bounds(
    tmp_path: Path,
) -> None:
    mc, fallback = _tables(tmp_path)
    manifest_path = mc / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["mc_sampling_plan"]["entries"][0]["replicas"] = 8

    quality_path = mc / "quality.csv"
    with quality_path.open(newline="") as stream:
        quality_rows = list(csv.DictReader(stream))
    quality_rows[0]["valid_replicates"] = "8"
    _write_csv(quality_path, quality_rows)
    manifest["tables"]["quality.csv"]["sha256"] = sha256(
        quality_path.read_bytes()
    ).hexdigest()
    manifest_path.write_text(json.dumps(manifest))

    decision = _decision(tmp_path, mc, fallback, attempt=1)

    assert decision["action"] == "accept_monte_carlo"
    assert [row["replicas"] for row in decision["current_sampling_plan"]] == [8, 4]


def test_attempt_two_requires_recorded_previous_plan(tmp_path: Path) -> None:
    mc, fallback = _tables(tmp_path)
    with pytest.raises(MonteCarloPolicyError, match="decision history"):
        decide_monte_carlo_closure(
            mc, fallback, attempt=2, output=tmp_path / "decision.json"
        )


def test_selection_rejects_stale_table_quality_hash(tmp_path: Path) -> None:
    mc, fallback = _tables(tmp_path)
    path = mc / "quality.csv"
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(MonteCarloPolicyError, match="hash mismatch"):
        _decision(tmp_path, mc, fallback, attempt=1)


def test_structural_stationarity_error_cannot_spend_more_time(tmp_path: Path) -> None:
    mc, fallback = _tables(
        tmp_path, failures={5.0: ["stationarity_provenance_missing"]}
    )
    assert _decision(tmp_path, mc, fallback, attempt=1)["action"] == "blocked"


def test_policy_rejects_unknown_typed_failure_axis(tmp_path: Path) -> None:
    mc, fallback = _tables(
        tmp_path,
        failures={5.0: ["mobility_rse"]},
        failure_axes={5.0: ["guess_from_text"]},
    )

    with pytest.raises(MonteCarloPolicyError, match="invalid MC quality evidence"):
        _decision(tmp_path, mc, fallback, attempt=1)


def test_quality_only_evidence_cannot_become_accepted_mc(tmp_path: Path) -> None:
    mc, fallback = _tables(tmp_path)
    path = mc / "manifest.json"
    manifest = json.loads(path.read_text())
    manifest["mc_qualification"] = {"coefficient_tables_available": False}
    path.write_text(json.dumps(manifest))
    assert _decision(tmp_path, mc, fallback, attempt=1)["action"] == "blocked"
