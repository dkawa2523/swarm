"""User-visible configuration, planning, execution, and resume use one document."""

from copy import deepcopy
import json
from pathlib import Path
import sqlite3

import pytest
import yaml

from electron_swarm.orchestration.preparation import prepare_run
from swarm_workflow.campaign.config import load_workflow
from swarm_workflow.campaign.document import read_calculation
from swarm_workflow.campaign.sweep import plan_sweep, run_sweep
from swarm_workflow.campaign.store import WorkflowStore
from swarm_workflow.cli import main
from product_helpers import base_product_config, write_config


def _input(tmp_path: Path, solver="two_term") -> tuple[dict, Path]:
    data = base_product_config(tmp_path, [solver])
    data["run"]["e_over_n_Td"] = [10.0, 30.0]
    data["campaign"] = {"database": "results/cases.sqlite"}
    if solver == "monte_carlo":
        data["solvers"]["monte_carlo"] = {
            "seed": 17,
            "particles": 2,
            "production_trial_periods": 8,
        }
        data["physics"]["ionization"]["population_model"] = (
            "fixed_particle_single_daughter"
        )
        data["campaign"]["mc"] = {"replicas": 2}
    if solver == "propagator":
        data["solvers"]["propagator"] = {
            "energy_cells": 64,
            "polar_cells": 8,
            "max_memory_mb": 128,
        }
    return data, write_config(tmp_path, data)


@pytest.mark.parametrize(
    "solver", ("two_term", "multi_term", "monte_carlo", "propagator")
)
def test_single_document_plans_executes_and_resumes_all_solvers(tmp_path, solver):
    data, path = _input(tmp_path, solver)
    plan = plan_sweep(path)
    assert plan["resolved_solvers"]["0"][solver]
    assert not (tmp_path / "results").exists()
    first = run_sweep(path)
    resumed = run_sweep(path)
    assert (
        first.cases_written
        == plan["expected_case_count"]
        == (4 if solver == "monte_carlo" else 2)
    )
    assert resumed.cases_written == 0
    with WorkflowStore(first.database_path) as store:
        key = next(iter(store.existing_case_keys()))
        restored = store.read_case(
            mixture_id=key[0], solver=key[1], e_over_n_Td=key[2], replicate=key[3]
        )
        assert restored.solver == solver
        assert restored.rf_phase is None
        assert (restored.energy_angle_distribution is not None) == (
            solver == "propagator"
        )
        assert (restored.eedf_counts is not None) == (solver == "monte_carlo")
    with sqlite3.connect(first.database_path) as connection:
        assert connection.execute("SELECT DISTINCT solver FROM cases").fetchall() == [
            (solver,)
        ]
        metadata = dict(connection.execute("SELECT key, value FROM metadata"))
        assert metadata["calculation_path"] == str(path.resolve())
        assert (
            metadata["calculation_sha256"]
            == read_calculation(path).calculation_sha256()
        )
        assert "base_config_path" not in metadata
    assert yaml.safe_load(path.read_text()) == data


@pytest.mark.parametrize("solver", ("two_term", "monte_carlo"))
def test_cli_plans_direct_and_campaign_without_creating_artifacts(
    tmp_path, capsys, solver
):
    data, path = _input(tmp_path, solver)
    main(["run", str(path), "--plan", "--json"])
    plan = json.loads(capsys.readouterr().out)
    assert plan["expected_case_count"] == (4 if solver == "monte_carlo" else 2)
    assert plan["database"] == str((tmp_path / "results/cases.sqlite").resolve())
    assert not (tmp_path / "results").exists()
    main(["run", str(path), "--plan"])
    display = capsys.readouterr().out
    assert f"MC cases={4 if solver == 'monte_carlo' else 0}" in display
    del data["campaign"]
    write_config(tmp_path, data)
    main(["run", str(path), "--plan"])
    assert "2 planned solver cases" in capsys.readouterr().out
    assert not (tmp_path / "prod_manifest.json").exists()


def test_seed_is_required_at_planning_for_direct_and_campaign(tmp_path):
    data, path = _input(tmp_path, "monte_carlo")
    del data["solvers"]["monte_carlo"]["seed"]
    write_config(tmp_path, data)
    with pytest.raises(ValueError, match="seed.*required"):
        plan_sweep(path)
    del data["campaign"]
    write_config(tmp_path, data)
    with pytest.raises(ValueError, match="seed.*required"):
        prepare_run(read_calculation(path).template)
    assert not (tmp_path / "results").exists()


@pytest.mark.parametrize("value", (None, [], "campaign", True))
def test_invalid_campaign_cannot_silently_become_a_direct_run(tmp_path, value):
    data, path = _input(tmp_path)
    data["campaign"] = value
    write_config(tmp_path, data)
    with pytest.raises(ValueError, match="campaign must be a mapping"):
        read_calculation(path)


def test_removed_base_reference_fails_with_migration_error(tmp_path):
    path = write_config(tmp_path, {"base_config": "missing.yaml"})
    with pytest.raises(ValueError, match="migration.*inline"):
        load_workflow(path)


def test_every_mixture_column_is_validated_before_execution(tmp_path):
    data, path = _input(tmp_path)
    data["conditions"]["gas_mixture"][0]["fraction"] = [1.0, 0.8]
    write_config(tmp_path, data)
    with pytest.raises(ValueError, match="sum to 1"):
        plan_sweep(path)
    assert not (tmp_path / "results").exists()


def test_pure_mixture_arrays_do_not_hide_duplicate_work(tmp_path):
    data, path = _input(tmp_path)
    data["conditions"]["gas_mixture"][0]["fraction"] = [1.0, 1.0]
    write_config(tmp_path, data)
    with pytest.raises(ValueError, match="duplicate mixtures"):
        load_workflow(path)


def test_calculation_identity_excludes_publication_and_campaign_controls(tmp_path):
    data, path = _input(tmp_path)
    first = read_calculation(path).calculation_sha256()
    changed = deepcopy(data)
    changed["output"]["directory"] = "elsewhere"
    changed["campaign"]["quality"] = {"major_rate_rse": 0.1}
    write_config(tmp_path, changed)
    assert read_calculation(path).calculation_sha256() == first
    changed["conditions"]["gas_temperature_K"] = 400
    write_config(tmp_path, changed)
    assert read_calculation(path).calculation_sha256() != first


def test_campaign_comparison_request_is_not_silently_discarded(tmp_path):
    data, path = _input(tmp_path)
    data["run"]["solvers"].append({"id": "multi_term"})
    data["comparison"]["enabled"] = True
    write_config(tmp_path, data)
    with pytest.raises(NotImplementedError, match="campaign comparison output"):
        plan_sweep(path)


def test_mc_anchor_selection_is_visible_in_the_plan(tmp_path):
    data, path = _input(tmp_path, "monte_carlo")
    data["campaign"]["mc"]["e_over_n_Td"] = [20.0, 40.0, 80.0]
    write_config(tmp_path, data)
    plan = plan_sweep(path)
    assert plan["solver_e_over_n_Td"]["monte_carlo"] == [20.0, 40.0, 80.0]
    assert plan["expected_case_count"] == 6
    assert {case["e_over_n_Td"] for case in plan["physical_cases"]["0"]} >= {
        20.0,
        40.0,
        80.0,
    }


def test_workers_execute_the_prepared_snapshot_when_yaml_is_edited_and_restored(
    tmp_path, monkeypatch
):
    import swarm_workflow.campaign.sweep as sweep

    _data, path = _input(tmp_path)
    original = path.read_text()
    execute = sweep._execute_deterministic_job

    def interleave(job):
        changed = yaml.safe_load(original)
        changed["conditions"]["gas_temperature_K"] = 900
        path.write_text(yaml.safe_dump(changed), encoding="utf-8")
        try:
            return execute(job)
        finally:
            path.write_text(original, encoding="utf-8")

    monkeypatch.setattr(sweep, "_execute_deterministic_job", interleave)
    summary = run_sweep(path)
    with sqlite3.connect(summary.database_path) as connection:
        for (metadata,) in connection.execute("SELECT metadata_json FROM cases"):
            assert json.loads(metadata)["physical_context"]["gas_temperature_K"] == 300


def test_worker_cannot_borrow_parent_provenance_for_changed_data(tmp_path, monkeypatch):
    import swarm_workflow.campaign.sweep as sweep

    data, path = _input(tmp_path)
    source = Path(data["cross_sections"]["files"][0]["path"])
    local = tmp_path / "cross_sections.csv"
    local.write_bytes(source.read_bytes())
    data["cross_sections"]["files"][0]["path"] = str(local)
    write_config(tmp_path, data)
    original = local.read_bytes()
    execute = sweep._execute_deterministic_job

    def interleave(job):
        local.write_bytes(original + b"\n")
        try:
            return execute(job)
        finally:
            local.write_bytes(original)

    monkeypatch.setattr(sweep, "_execute_deterministic_job", interleave)
    with pytest.raises(RuntimeError, match="worker collision/angular inputs"):
        run_sweep(path)
    with sqlite3.connect(tmp_path / "results/cases.sqlite") as connection:
        assert connection.execute("SELECT COUNT(*) FROM cases").fetchone()[0] == 0


def test_campaign_cli_reports_partial_when_physics_skips_a_requested_solver(
    tmp_path, capsys
):
    data, path = _input(tmp_path)
    data["run"]["solvers"].append({"id": "propagator"})
    data["physics"]["electron_electron"].update(
        enabled=True, model="relaxation_postprocess"
    )
    data["feature_policy"]["unsupported"] = "skip_solver"
    write_config(tmp_path, data)
    with pytest.raises(SystemExit) as stopped:
        main(["run", str(path)])
    assert stopped.value.code == 4
    assert "partial" in capsys.readouterr().out
    with sqlite3.connect(tmp_path / "results/cases.sqlite") as connection:
        assert connection.execute("SELECT DISTINCT solver FROM cases").fetchall() == [
            ("two_term",)
        ]
        (raw,) = connection.execute(
            "SELECT value FROM metadata WHERE key='execution_coverage_json'"
        ).fetchone()
        coverage = json.loads(raw)
        assert coverage == {
            "requested_case_count": 4,
            "expected_case_count": 2,
            "skipped_case_count": 2,
            "coverage_status": "partial",
        }
