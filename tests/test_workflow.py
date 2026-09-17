from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
import json
import sqlite3
from pathlib import Path

import pytest
import yaml

from electron_swarm import load_config
import swarm_workflow.campaign.provenance as provenance_module
import swarm_workflow.campaign.sweep as sweep_module
from swarm_workflow.campaign.aggregate import (
    aggregate_database,
)
from swarm_workflow.quality.policy import QualityThresholds, quality_thresholds_json
from swarm_workflow.campaign.config import load_workflow
from swarm_workflow.campaign.store import WorkflowSchemaError, WorkflowStore
from swarm_workflow.campaign.sweep import run_sweep
from swarm_workflow.cli import main as workflow_cli_main
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
            {
                "base_config": base_path.name,
                "database": "outputs/swarm.sqlite",
                "e_over_n_Td": [10.0, 20.0],
                "mixtures": [
                    {"Ar": 1.0, "O2": 0.0},
                    {"Ar": 0.99, "O2": 0.01},
                ],
                "mc": {"replicas": 1, "base_seed": 12345},
            },
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
    workflow_cli_main(["sweep", str(workflow_path)])

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
    case = sweep_module.run(load_config(base_path), write=False).cases[0]
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
    raw["e_over_n_Td"].append(30.0)
    workflow_path.write_text(
        yaml.safe_dump(raw, sort_keys=False),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="workflow_config_sha256"):
        run_sweep(workflow_path)

    with sqlite3.connect(tmp_path / "outputs" / "swarm.sqlite") as connection:
        assert connection.execute("SELECT COUNT(*) FROM cases").fetchone() == (4,)
        assert connection.execute(
            "SELECT COUNT(*) FROM cases WHERE e_over_n_Td = 30.0"
        ).fetchone() == (0,)


def test_workflow_rejects_obsolete_hash_with_unused_mc_plan(
    tmp_path: Path,
) -> None:
    base_path = _write_workflow_base(tmp_path)
    workflow_path = _write_workflow(tmp_path, base_path)
    summary = run_sweep(workflow_path)
    workflow = load_workflow(workflow_path)
    assert workflow.mc_enabled is False
    obsolete_payload = {
        "base_config_path": str(workflow.base_config_path),
        "database_path": str(workflow.database_path),
        "e_over_n_Td": list(workflow.e_over_n_Td),
        "mixtures": [
            {
                "mixture_id": mixture.mixture_id,
                "fractions": dict(sorted(mixture.fractions.items())),
            }
            for mixture in workflow.mixtures
        ],
        "mc": {
            "e_over_n_Td": None,
            "replicas": workflow.mc_replicas,
            "base_seed": workflow.mc_base_seed,
            "workers": workflow.mc_workers,
            "sampling_plan": [asdict(entry) for entry in workflow.mc_sampling_plan],
        },
        "quality": asdict(workflow.quality),
    }
    obsolete_hash = sha256(
        json.dumps(
            obsolete_payload,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()
    assert obsolete_hash != provenance_module.workflow_config_sha256(workflow)
    with sqlite3.connect(summary.database_path) as connection:
        connection.execute(
            "UPDATE metadata SET value = ? WHERE key = 'workflow_config_sha256'",
            (obsolete_hash,),
        )
        connection.commit()

    with pytest.raises(ValueError, match="workflow_config_sha256"):
        run_sweep(workflow_path)


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
    raw["quality"] = {
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


def test_workflow_rejects_species_not_present_in_base_config(tmp_path: Path) -> None:
    base_path = _write_workflow_base(tmp_path)
    workflow_path = tmp_path / "bad_workflow.yaml"
    workflow_path.write_text(
        yaml.safe_dump(
            {
                "base_config": base_path.name,
                "database": "outputs/swarm.sqlite",
                "e_over_n_Td": [10.0],
                "mixtures": [{"Ar": 1.0, "N2": 0.1}],
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="unknown species"):
        load_workflow(workflow_path)


def test_workflow_loads_independent_monte_carlo_anchor_grid(tmp_path: Path) -> None:
    base_path = _write_workflow_base(tmp_path)
    workflow_path = _write_workflow(tmp_path, base_path)
    raw = yaml.safe_load(workflow_path.read_text(encoding="utf-8"))
    raw["mc"]["e_over_n_Td"] = [20.0, 100.0]
    workflow_path.write_text(
        yaml.safe_dump(raw, sort_keys=False),
        encoding="utf-8",
    )

    workflow = load_workflow(workflow_path)

    assert workflow.e_over_n_Td == (10.0, 20.0)
    assert workflow.mc_e_over_n_Td == (20.0, 100.0)


def test_electron_swarm_does_not_import_workflow_or_sqlite() -> None:
    root = Path(__file__).resolve().parents[1]
    for path in (root / "electron_swarm").rglob("*.py"):
        text = path.read_text(encoding="utf-8")
        assert "swarm_workflow" not in text
        assert "sqlite3" not in text
