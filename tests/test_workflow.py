from __future__ import annotations

from product_helpers import campaign_config

import csv
from dataclasses import asdict
import json
import sqlite3
from pathlib import Path

import pytest
import yaml

from swarm_workflow.campaign.document import load_case_template
import swarm_workflow.campaign.sweep as sweep_module
from swarm_workflow.campaign.aggregate import (
    aggregate_database,
)
from swarm_workflow.quality.policy import QualityThresholds, quality_thresholds_json
from swarm_workflow.campaign.config import load_workflow
from swarm_workflow.campaign.store import WorkflowSchemaError, WorkflowStore
from swarm_workflow.campaign.sweep import run_sweep
from swarm_workflow.cli import main as workflow_cli_main
from swarm_workflow.tables import build_tables
from product_helpers import base_product_config, write_config


def _write_mixture_xs(tmp_path: Path) -> Path:
    path = tmp_path / "ar_o2_minimal.csv"
    rows = [
        "species,process,type,threshold_eV,mass_amu,energy_eV,cross_section_m2",
    ]
    for species, mass, scale in [("Ar", 39.948, 1.0), ("O2", 31.998, 1.4)]:
        for energy in [0.0, 1.0, 10.0, 100.0, 1000.0]:
            rows.append(
                f"{species},momentum,momentum,0.0,{mass},{energy},{scale * 1.0e-20}"
            )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def test_changed_solver_implementation_cannot_resume_or_build_tables(
    tmp_path: Path,
) -> None:
    from electron_swarm.core.source_identity import SOLVER_IMPLEMENTATIONS_KEY
    from swarm_workflow.tables.contracts import TableBuildError

    workflow_path = _write_workflow(tmp_path, _write_workflow_base(tmp_path))
    sweep = run_sweep(workflow_path)
    with sqlite3.connect(sweep.database_path) as connection:
        value = connection.execute(
            "SELECT value FROM metadata WHERE key=?", (SOLVER_IMPLEMENTATIONS_KEY,)
        ).fetchone()[0]
        identities = json.loads(value)
        identities["two_term"]["source_sha256"] = "0" * 64
        connection.execute(
            "UPDATE metadata SET value=? WHERE key=?",
            (json.dumps(identities), SOLVER_IMPLEMENTATIONS_KEY),
        )
    with pytest.raises(ValueError, match="implementation or numerical runtime changed"):
        run_sweep(workflow_path)
    with pytest.raises(
        TableBuildError, match="implementation or numerical runtime changed"
    ):
        build_tables(sweep.database_path, tmp_path / "tables", source="two_term")


def _write_workflow_base(tmp_path: Path) -> Path:
    xs_path = _write_mixture_xs(tmp_path)
    data = base_product_config(tmp_path, ["two_term"])
    data["run"]["e_over_n_Td"] = [50.0]
    data["conditions"]["gas_mixture"] = [
        {"species": "Ar", "fraction": 1.0, "mass_amu": 39.948},
        {"species": "O2", "fraction": 0.0, "mass_amu": 31.998},
    ]
    data["cross_sections"]["files"] = [{"path": xs_path.as_posix(), "format": "csv"}]
    return write_config(tmp_path, data, "base.yaml")


def _write_workflow(tmp_path: Path, base_path: Path) -> Path:
    path = tmp_path / "workflow.yaml"
    path.write_text(
        yaml.safe_dump(
            campaign_config(base_path, database='outputs/swarm.sqlite', e_over_n_Td=[10.0, 20.0], mixtures=[{'Ar': 1.0, 'O2': 0.0}, {'Ar': 0.99, 'O2': 0.01}]),
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    return path


def test_workflow_rejects_duplicate_base_solver_ids(tmp_path: Path) -> None:
    base_path = _write_workflow_base(tmp_path)
    raw = yaml.safe_load(base_path.read_text(encoding="utf-8"))
    raw["run"]["solvers"].append({"id": "two_term", "enabled": False})
    base_path.write_text(
        yaml.safe_dump(raw, sort_keys=False),
        encoding="utf-8",
    )
    workflow_path = _write_workflow(tmp_path, base_path)

    with pytest.raises(
        ValueError,
        match=r"run\.solvers must not contain duplicate solver ids",
    ):
        run_sweep(workflow_path)
    assert not (tmp_path / "outputs" / "swarm.sqlite").exists()


def test_two_term_workflow_sweep_resumes_without_duplicate_rows(tmp_path: Path) -> None:
    base_path = _write_workflow_base(tmp_path)
    workflow_path = _write_workflow(tmp_path, base_path)

    first = run_sweep(workflow_path)
    # Provenance hashes validated semantics, not YAML key order or formatting.
    workflow_path.write_text(
        yaml.safe_dump(
            yaml.safe_load(workflow_path.read_text(encoding="utf-8")),
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    with sqlite3.connect(first.database_path) as connection:
        connection.execute("DELETE FROM aggregate_quality")
        connection.commit()
    second = run_sweep(workflow_path)
    workflow_cli_main(["run", str(workflow_path)])

    assert first.cases_written == 4
    assert second.cases_written == 0
    with WorkflowStore(first.database_path) as store:
        assert len(store.existing_case_keys()) == 4
        assert store.has_case(
            mixture_id=0,
            solver="two_term",
            e_over_n_Td=10.0,
            replicate=0,
        )
    with sqlite3.connect(tmp_path / "outputs" / "swarm.sqlite") as conn:
        assert conn.execute("SELECT COUNT(*) FROM mixtures").fetchone()[0] == 2
        assert conn.execute("SELECT COUNT(*) FROM cases").fetchone()[0] == 4
        assert conn.execute("SELECT COUNT(*) FROM metadata").fetchone()[0] >= 3
        assert conn.execute("SELECT COUNT(*) FROM aggregate_quality").fetchone() == (4,)
        assert conn.execute(
            "SELECT LENGTH(value) FROM metadata WHERE key = 'workflow_config_sha256'"
        ).fetchone() == (64,)
        assert conn.execute(
            "SELECT value FROM metadata WHERE key = 'quality_thresholds_json'"
        ).fetchone() == (quality_thresholds_json(QualityThresholds()),)
        assert {
            row[0]
            for row in conn.execute(
                "SELECT DISTINCT solver FROM cases ORDER BY solver"
            ).fetchall()
        } == {"two_term"}
        norms = conn.execute(
            """
            SELECT mixture_id, solver, e_over_n_Td, replicate,
                   SUM(eedf * energy_width_eV)
            FROM eedf_bins
            GROUP BY mixture_id, solver, e_over_n_Td, replicate
            ORDER BY mixture_id, e_over_n_Td
            """
        ).fetchall()

    assert len(norms) == 4
    for *_case_key, norm in norms:
        assert norm == pytest.approx(1.0, rel=1.0e-8, abs=1.0e-10)


def test_bounded_two_term_workflow_preserves_eedf_support_through_tables(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["two_term"])
    data["physics"]["energy_grid_policy"]["max_eV_limit"] = 20.0
    base_path = write_config(tmp_path, data, "base.yaml")
    workflow_path = tmp_path / "workflow.yaml"
    workflow_path.write_text(
        yaml.safe_dump(
            campaign_config(base_path, database='outputs/swarm.sqlite', e_over_n_Td=[0.1, 0.2], mixtures=[{'Ar': 1.0}]),
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    sweep = run_sweep(workflow_path)
    resumed = run_sweep(workflow_path)
    aggregate_database(sweep.database_path)
    tables = tmp_path / "tables"
    build_tables(sweep.database_path, tables, source="two_term")

    assert sweep.cases_written == 2
    assert resumed.cases_written == 0
    with sqlite3.connect(sweep.database_path) as connection:
        raw = connection.execute(
            "SELECT e_over_n_Td, SUM(energy_width_eV), "
            "SUM(eedf * energy_width_eV) FROM eedf_bins "
            "GROUP BY e_over_n_Td ORDER BY e_over_n_Td"
        ).fetchall()
        aggregate = connection.execute(
            "SELECT e_over_n_Td, MIN(energy_left_eV), MAX(energy_right_eV), "
            "SUM(mean_probability_mass) FROM aggregate_eedf_bins "
            "GROUP BY e_over_n_Td ORDER BY e_over_n_Td"
        ).fetchall()
    for raw_row, aggregate_row in zip(raw, aggregate, strict=True):
        assert raw_row[1] == pytest.approx(20.0)
        assert raw_row[2] == pytest.approx(1.0)
        assert aggregate_row == pytest.approx((raw_row[0], 0.0, 20.0, 1.0))
    with (tables / "mixture_0000" / "eedf.csv").open(
        encoding="utf-8", newline=""
    ) as stream:
        rows = list(csv.DictReader(stream))
    assert rows
    for field in (0.1, 0.2):
        selected = [row for row in rows if float(row["E_over_N_Td"]) == field]
        assert sum(float(row["energy_width_eV"]) for row in selected) == pytest.approx(
            20.0
        )
        assert sum(
            float(row["eedf"]) * float(row["energy_width_eV"]) for row in selected
        ) == pytest.approx(1.0)
    with (tables / "mixture_0000" / "quality.csv").open(
        encoding="utf-8", newline=""
    ) as stream:
        quality = list(csv.DictReader(stream))
    assert len(quality) == 2
    for row in quality:
        assert row["aggregate_quality_passed"] == "1"
        assert row["solver_converged"] == "1"
        assert row["passed"] == "0"
        assert row["solver_grid_limit_hit"] == "1"
        assert json.loads(row["failure_reasons_json"]) == ["energy_grid_limit_hit"]
    manifest = json.loads(
        (tables / "mixture_0000" / "manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["quality_summary"]["passed"] is False
    assert manifest["quality_summary"]["failed_points"] == 2


def test_two_term_workflow_resumes_only_missing_cases_after_interruption(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    base_path = _write_workflow_base(tmp_path)
    workflow_path = _write_workflow(tmp_path, base_path)
    original_write_case = WorkflowStore.write_case
    writes = 0

    def interrupt_after_first_commit(
        store: WorkflowStore,
        **kwargs: object,
    ) -> None:
        nonlocal writes
        original_write_case(store, **kwargs)
        writes += 1
        if writes == 1:
            raise RuntimeError("simulated interruption")

    monkeypatch.setattr(WorkflowStore, "write_case", interrupt_after_first_commit)
    with pytest.raises(RuntimeError, match="simulated interruption"):
        run_sweep(workflow_path)

    database = tmp_path / "outputs" / "swarm.sqlite"
    with sqlite3.connect(database) as connection:
        assert connection.execute("SELECT COUNT(*) FROM cases").fetchone() == (1,)
        assert connection.execute(
            "SELECT COUNT(DISTINCT mixture_id || ':' || solver || ':' || "
            "e_over_n_Td || ':' || replicate) FROM eedf_bins"
        ).fetchone() == (1,)

    monkeypatch.setattr(WorkflowStore, "write_case", original_write_case)
    original_run = sweep_module.run
    resumed_calls: list[tuple[str, tuple[float, ...], tuple[str, ...]]] = []

    def record_run(config: object, **kwargs: object) -> object:
        run_config = config
        resumed_calls.append(
            (
                run_config.run.case_prefix,
                tuple(run_config.run.e_over_n_Td),
                tuple(item.id for item in run_config.run.solvers),
            )
        )
        return original_run(run_config, **kwargs)

    monkeypatch.setattr(sweep_module, "run", record_run)
    resumed = run_sweep(workflow_path)

    assert resumed.cases_written == 3
    assert resumed_calls == [
        ("Prod_m0000_r0000", (20.0,), ("two_term",)),
        ("Prod_m0001_r0000", (10.0, 20.0), ("two_term",)),
    ]
    with sqlite3.connect(database) as connection:
        assert connection.execute("SELECT COUNT(*) FROM cases").fetchone() == (4,)
        assert connection.execute(
            "SELECT COUNT(*) FROM aggregate_quality"
        ).fetchone() == (4,)


def test_write_case_parent_and_children_rollback_as_one_checkpoint(
    tmp_path: Path,
) -> None:
    base_path = _write_workflow_base(tmp_path)
    case = sweep_module.run(load_case_template(base_path), write=False).cases[0]
    database = tmp_path / "atomic.sqlite"

    with WorkflowStore(database) as store:
        store.write_mixture(
            0,
            [("Ar", 1.0, 39.948), ("O2", 0.0, 31.998)],
        )
        store.connection.execute(
            """
            CREATE TRIGGER reject_eedf_insert
            BEFORE INSERT ON eedf_bins
            BEGIN
                SELECT RAISE(ABORT, 'simulated child failure');
            END
            """
        )
        store.connection.commit()

        with pytest.raises(sqlite3.IntegrityError, match="simulated child failure"):
            store.write_case(mixture_id=0, replicate=0, case=case)

        assert not store.has_case(
            mixture_id=0,
            solver="two_term",
            e_over_n_Td=50.0,
            replicate=0,
        )
        for table in ("cases", "rates", "eedf_bins"):
            assert store.connection.execute(
                f"SELECT COUNT(*) FROM {table}"
            ).fetchone() == (0,)


def test_workflow_rejects_changed_config_when_reusing_database(
    tmp_path: Path,
) -> None:
    base_path = _write_workflow_base(tmp_path)
    workflow_path = _write_workflow(tmp_path, base_path)
    run_sweep(workflow_path)

    raw = yaml.safe_load(workflow_path.read_text(encoding="utf-8"))
    raw["run"]["e_over_n_Td"].append(30.0)
    workflow_path.write_text(
        yaml.safe_dump(raw, sort_keys=False),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="different .*; use a separate database"):
        run_sweep(workflow_path)

    with sqlite3.connect(tmp_path / "outputs" / "swarm.sqlite") as connection:
        assert connection.execute("SELECT COUNT(*) FROM cases").fetchone() == (4,)
        assert connection.execute(
            "SELECT COUNT(*) FROM cases WHERE e_over_n_Td = 30.0"
        ).fetchone() == (0,)


def test_workflow_rejects_unused_mc_plan(
    tmp_path: Path,
) -> None:
    base_path = _write_workflow_base(tmp_path)
    workflow_path = _write_workflow(tmp_path, base_path)
    raw = yaml.safe_load(workflow_path.read_text(encoding="utf-8"))
    raw["campaign"]["mc"] = {"replicas": 3}
    workflow_path.write_text(yaml.safe_dump(raw), encoding="utf-8")
    with pytest.raises(ValueError, match="mc settings require an enabled monte_carlo"):
        run_sweep(workflow_path)
    assert not (tmp_path / "outputs" / "swarm.sqlite").exists()


def test_quality_only_workflow_change_reaggregates_without_solver_rerun(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    base_path = _write_workflow_base(tmp_path)
    workflow_path = _write_workflow(tmp_path, base_path)
    first = run_sweep(workflow_path)
    source_quality = QualityThresholds()
    evaluation_quality = QualityThresholds(
        required_rate_min_process_peak_fraction=1.0e-5
    )

    with sqlite3.connect(first.database_path) as connection:
        source_metadata = dict(
            connection.execute("SELECT key, value FROM metadata").fetchall()
        )

    raw = yaml.safe_load(workflow_path.read_text(encoding="utf-8"))
    raw["campaign"]["quality"] = {
        "required_rate_min_process_peak_fraction": 1.0e-5,
    }
    workflow_path.write_text(
        yaml.safe_dump(raw, sort_keys=False),
        encoding="utf-8",
    )

    def reject_solver_call(*_args: object, **_kwargs: object) -> object:
        raise AssertionError("quality-only reevaluation must not run a solver")

    monkeypatch.setattr(sweep_module, "run", reject_solver_call)
    resumed = run_sweep(workflow_path)

    assert resumed.cases_written == 0
    assert resumed.quality_policy_reevaluated is True
    with sqlite3.connect(resumed.database_path) as connection:
        assert (
            dict(connection.execute("SELECT key, value FROM metadata").fetchall())
            == source_metadata
        )
        assert source_metadata["quality_thresholds_json"] == quality_thresholds_json(
            source_quality
        )
        policy_rows = connection.execute(
            """
            SELECT DISTINCT source_thresholds_json, thresholds_json,
                            quality_policy_reevaluated
            FROM aggregate_quality
            """
        ).fetchall()
        assert policy_rows == [
            (
                quality_thresholds_json(source_quality),
                quality_thresholds_json(evaluation_quality),
                1,
            )
        ]

    # A later standalone aggregation restores the last evaluation policy instead
    # of silently reverting to the immutable source policy.
    aggregate_database(resumed.database_path)
    with sqlite3.connect(resumed.database_path) as connection:
        assert connection.execute(
            """
            SELECT DISTINCT source_thresholds_json, thresholds_json,
                            quality_policy_reevaluated
            FROM aggregate_quality
            """
        ).fetchall() == [
            (
                quality_thresholds_json(source_quality),
                quality_thresholds_json(evaluation_quality),
                1,
            )
        ]
        assert (
            dict(connection.execute("SELECT key, value FROM metadata").fetchall())
            == source_metadata
        )


def test_workflow_rejects_pre_rate_relevance_quality_provenance(
    tmp_path: Path,
) -> None:
    base_path = _write_workflow_base(tmp_path)
    workflow_path = _write_workflow(tmp_path, base_path)
    first = run_sweep(workflow_path)
    source_quality = QualityThresholds()
    legacy_quality_payload = asdict(source_quality)
    legacy_quality_payload.pop("required_rate_min_process_peak_fraction")
    legacy_quality_json = json.dumps(
        legacy_quality_payload,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )

    with sqlite3.connect(first.database_path) as connection:
        connection.execute(
            "UPDATE metadata SET value = ? WHERE key = 'quality_thresholds_json'",
            (legacy_quality_json,),
        )
        connection.commit()

    with pytest.raises(
        WorkflowSchemaError,
        match="quality_thresholds_json is not canonical",
    ):
        aggregate_database(first.database_path)


def test_workflow_rejects_legacy_database_without_workflow_hash(
    tmp_path: Path,
) -> None:
    base_path = _write_workflow_base(tmp_path)
    workflow_path = _write_workflow(tmp_path, base_path)
    summary = run_sweep(workflow_path)
    with sqlite3.connect(summary.database_path) as connection:
        connection.execute("DELETE FROM metadata WHERE key = 'workflow_config_sha256'")
        connection.commit()

    with pytest.raises(
        WorkflowSchemaError,
        match="missing workflow_config_sha256",
    ):
        run_sweep(workflow_path)


def test_workflow_rejects_legacy_database_without_quality_policy(
    tmp_path: Path,
) -> None:
    base_path = _write_workflow_base(tmp_path)
    workflow_path = _write_workflow(tmp_path, base_path)
    summary = run_sweep(workflow_path)
    with sqlite3.connect(summary.database_path) as connection:
        connection.execute("DELETE FROM metadata WHERE key = 'quality_thresholds_json'")
        connection.commit()

    with pytest.raises(
        WorkflowSchemaError,
        match="missing quality_thresholds_json",
    ):
        run_sweep(workflow_path)


def test_workflow_loads_independent_monte_carlo_anchor_grid(tmp_path: Path) -> None:
    base_path = _write_workflow_base(tmp_path)
    base = yaml.safe_load(base_path.read_text(encoding="utf-8"))
    base["run"]["solvers"].append({"id": "monte_carlo", "enabled": True})
    base["solvers"]["monte_carlo"] = {"seed": 12345}
    base_path.write_text(yaml.safe_dump(base), encoding="utf-8")
    workflow_path = _write_workflow(tmp_path, base_path)
    raw = yaml.safe_load(workflow_path.read_text(encoding="utf-8"))
    raw["campaign"]["mc"] = {"e_over_n_Td": [20.0, 100.0]}
    workflow_path.write_text(
        yaml.safe_dump(raw, sort_keys=False),
        encoding="utf-8",
    )

    workflow = load_workflow(workflow_path)

    assert workflow.e_over_n_Td == (10.0, 20.0)
    assert workflow.mc_e_over_n_Td == (20.0, 100.0)
