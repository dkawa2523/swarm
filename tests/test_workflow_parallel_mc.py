from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
import sqlite3
from types import SimpleNamespace

import pytest
import yaml

from product_helpers import base_product_config, write_config
from electron_swarm.solvers.monte_carlo.evidence import (
    DIRECT_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
    MC_CASE_SEED_DERIVATION,
    MC_EEDF_ESTIMATOR_SCHEMA_VERSION,
    MC_SEED_DERIVATION_SCHEMA_VERSION,
    WEIGHTED_MC_TAIL_ESTIMATOR_SCHEMA_VERSION,
    WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
)
from electron_swarm.solvers.monte_carlo.source_identity import (
    monte_carlo_source_sha256,
)
import swarm_workflow.campaign.sweep as sweep_module
from swarm_workflow.campaign.sweep import (
    MC_SAMPLING_PLAN_METADATA_KEY,
    run_sweep,
)
from swarm_workflow.campaign.store import (
    WorkflowSchemaError,
    canonical_mc_sampling_plan_json,
    provenance_hash_manifest,
)
from swarm_workflow.quality.propagator_source import (
    PROPAGATOR_SOLVER_SOURCE_METADATA_KEY,
)
from swarm_workflow.campaign.config import load_workflow


def test_provenance_manifest_scopes_solver_specific_hashes() -> None:
    metadata = {
        "workflow_config_sha256": "a" * 64,
        "base_config_sha256": "b" * 64,
        "cross_sections_sha256": "c" * 64,
        "mc_solver_source_sha256": "d" * 64,
        PROPAGATOR_SOLVER_SOURCE_METADATA_KEY: "e" * 64,
    }

    mc_hashes = provenance_hash_manifest(metadata, source="monte_carlo")
    propagator_hashes = provenance_hash_manifest(metadata, source="propagator")

    assert "mc_solver_source_sha256" in mc_hashes
    assert PROPAGATOR_SOLVER_SOURCE_METADATA_KEY not in mc_hashes
    assert "mc_solver_source_sha256" not in propagator_hashes
    assert PROPAGATOR_SOLVER_SOURCE_METADATA_KEY in propagator_hashes


def _parallel_mc_workflow(
    tmp_path: Path,
    *,
    workers: int | None,
    population_model: str = "fixed_particle_single_daughter",
) -> Path:
    tmp_path.mkdir(parents=True, exist_ok=True)
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["run"]["e_over_n_Td"] = [20.0]
    data["solvers"]["monte_carlo"] = {
        "population_model": population_model,
        "particles": 2,
        "max_collisions": 10,
        "seed": 7,
    }
    base = write_config(tmp_path, data, "mc_base.yaml")
    mc = {
        "replicas": 2,
        "base_seed": 321,
        "e_over_n_Td": [20.0, 40.0],
    }
    if workers is not None:
        mc["workers"] = workers
    workflow = tmp_path / "workflow.yaml"
    workflow.write_text(
        yaml.safe_dump(
            {
                "base_config": base.name,
                "database": "outputs/swarm.sqlite",
                "e_over_n_Td": [20.0],
                "mixtures": [{"Ar": 1.0}],
                "mc": mc,
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    return workflow


def test_workflow_mc_workers_defaults_to_one_and_validates(tmp_path: Path) -> None:
    default_path = _parallel_mc_workflow(tmp_path / "default", workers=None)
    assert load_workflow(default_path).mc_workers == 1

    bad_path = _parallel_mc_workflow(tmp_path / "bad", workers=0)
    with pytest.raises(ValueError, match=r"mc\.workers"):
        load_workflow(bad_path)


def test_workflow_resolves_optional_mc_reuse_database(tmp_path: Path) -> None:
    workflow_path = _parallel_mc_workflow(tmp_path, workers=1)
    raw = yaml.safe_load(workflow_path.read_text(encoding="utf-8"))
    raw["mc"]["reuse_database"] = "cache/source.sqlite"
    workflow_path.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")

    workflow = load_workflow(workflow_path)

    assert (
        workflow.mc_reuse_database == (tmp_path / "cache" / "source.sqlite").resolve()
    )


def test_mc_workflow_requires_explicit_integer_base_seed(tmp_path: Path) -> None:
    workflow_path = _parallel_mc_workflow(tmp_path, workers=1)
    raw = yaml.safe_load(workflow_path.read_text(encoding="utf-8"))
    del raw["mc"]["base_seed"]
    workflow_path.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
    with pytest.raises(ValueError, match=r"mc\.base_seed.*required.*integer"):
        load_workflow(workflow_path)

    for invalid in (None, True, 1.5, "321"):
        raw["mc"]["base_seed"] = invalid
        workflow_path.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
        with pytest.raises(ValueError, match=r"mc\.base_seed.*integer"):
            load_workflow(workflow_path)


def test_parallel_mc_results_are_yielded_in_completion_order(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    class FakeFuture:
        def __init__(self, value: object) -> None:
            self.value = value

        def result(self) -> object:
            return self.value

    class FakeExecutor:
        def __init__(self, **_kwargs: object) -> None:
            pass

        def __enter__(self) -> "FakeExecutor":
            return self

        def __exit__(self, *_args: object) -> None:
            return None

        def submit(self, _function: object, job: object) -> FakeFuture:
            return FakeFuture(job)

    monkeypatch.setattr(sweep_module, "ProcessPoolExecutor", FakeExecutor)
    monkeypatch.setattr(
        sweep_module,
        "as_completed",
        lambda futures: reversed(futures),
    )

    assert list(sweep_module._iter_monte_carlo_results([1, 2, 3], workers=2)) == [
        3,
        2,
        1,
    ]


def _write_sampling_plan(workflow_path: Path, rows: list[dict[str, object]]) -> None:
    raw = yaml.safe_load(workflow_path.read_text(encoding="utf-8"))
    raw["mc"]["sampling_plan"] = rows
    workflow_path.write_text(
        yaml.safe_dump(raw, sort_keys=False),
        encoding="utf-8",
    )


def _sampling_row(e_over_n_Td: float, **updates: object) -> dict[str, object]:
    row: dict[str, object] = {
        "e_over_n_Td": e_over_n_Td,
        "particles": 3,
        "warmup_collisions": 4,
        "max_collisions": 20,
        "replicas": 3,
    }
    row.update(updates)
    return row


def test_mc_sampling_plan_resolves_sparse_overrides_against_base_defaults(
    tmp_path: Path,
) -> None:
    workflow_path = _parallel_mc_workflow(tmp_path, workers=1)
    _write_sampling_plan(workflow_path, [_sampling_row(40.0)])

    workflow = load_workflow(workflow_path)

    assert [asdict(entry) for entry in workflow.mc_sampling_plan] == [
        {
            "e_over_n_Td": 20.0,
            "particles": 2,
            "warmup_collisions": 0,
            "max_collisions": 10,
            "tail_max_collisions": 0,
            "replicas": 2,
            "transport_correlation_lag_barriers": 64,
            "transport_estimator": "single_field",
        },
        {
            "e_over_n_Td": 40.0,
            "particles": 3,
            "warmup_collisions": 4,
            "max_collisions": 20,
            "tail_max_collisions": 0,
            "replicas": 3,
            "transport_correlation_lag_barriers": 64,
            "transport_estimator": "single_field",
        },
    ]


def test_mc_sampling_plan_row_inherits_omitted_controls(tmp_path: Path) -> None:
    workflow_path = _parallel_mc_workflow(tmp_path, workers=1)
    _write_sampling_plan(
        workflow_path,
        [{"e_over_n_Td": 40.0, "warmup_collisions": 8}],
    )

    workflow = load_workflow(workflow_path)

    assert asdict(workflow.mc_sampling_plan[1]) == {
        "e_over_n_Td": 40.0,
        "particles": 2,
        "warmup_collisions": 8,
        "max_collisions": 10,
        "tail_max_collisions": 0,
        "replicas": 2,
        "transport_correlation_lag_barriers": 64,
        "transport_estimator": "single_field",
    }


def test_mc_sampling_plan_explicit_tail_opt_in_and_null_disable(
    tmp_path: Path,
) -> None:
    workflow_path = _parallel_mc_workflow(
        tmp_path,
        workers=1,
        population_model="weighted_branching",
    )
    base_path = tmp_path / "mc_base.yaml"
    base = yaml.safe_load(base_path.read_text(encoding="utf-8"))
    base["solvers"]["monte_carlo"]["tail_max_collisions"] = 12
    base_path.write_text(yaml.safe_dump(base, sort_keys=False), encoding="utf-8")
    _write_sampling_plan(
        workflow_path,
        [
            {"e_over_n_Td": 20.0, "tail_max_collisions": None},
            {"e_over_n_Td": 40.0, "tail_max_collisions": 7},
        ],
    )

    workflow = load_workflow(workflow_path)

    assert [entry.tail_max_collisions for entry in workflow.mc_sampling_plan] == [
        0,
        7,
    ]


@pytest.mark.parametrize(
    ("rows", "match"),
    [
        ([_sampling_row(20.0, unsupported=1)], "unsupported fields"),
        ([_sampling_row(20.0), _sampling_row(20.0)], "duplicate E/N anchor"),
        ([_sampling_row(30.0)], "not an MC anchor"),
        ([_sampling_row(20.0, particles=0)], "particles must be a positive"),
        (
            [_sampling_row(20.0, warmup_collisions=-1)],
            "warmup_collisions must be a nonnegative",
        ),
        ([_sampling_row(20.0, max_collisions=-1)], "max_collisions must be a positive"),
        (
            [_sampling_row(20.0, tail_max_collisions=0)],
            "tail_max_collisions must be a positive",
        ),
        ([_sampling_row(20.0, replicas=0)], "replicas must be a positive"),
        (
            [_sampling_row(20.0, transport_correlation_lag_barriers=6)],
            "must be a power of two",
        ),
        (
            [_sampling_row(20.0, transport_correlation_lag_barriers=True)],
            "must be a positive integer",
        ),
        (
            [_sampling_row(20.0, transport_estimator="automatic")],
            "transport_estimator",
        ),
    ],
)
def test_mc_sampling_plan_rejects_invalid_rows(
    tmp_path: Path,
    rows: list[dict[str, object]],
    match: str,
) -> None:
    workflow_path = _parallel_mc_workflow(tmp_path, workers=1)
    _write_sampling_plan(workflow_path, rows)

    with pytest.raises(ValueError, match=match):
        load_workflow(workflow_path)


def test_mc_sampling_plan_rejects_duplicate_anchor_grid(tmp_path: Path) -> None:
    workflow_path = _parallel_mc_workflow(tmp_path, workers=1)
    raw = yaml.safe_load(workflow_path.read_text(encoding="utf-8"))
    raw["mc"]["e_over_n_Td"] = [20.0, 20.0]
    workflow_path.write_text(
        yaml.safe_dump(raw, sort_keys=False),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="MC E/N anchors must be unique"):
        load_workflow(workflow_path)


def test_mc_compute_budget_rejects_oversized_plan_before_execution(
    tmp_path: Path,
) -> None:
    workflow_path = _parallel_mc_workflow(tmp_path, workers=1)
    raw = yaml.safe_load(workflow_path.read_text(encoding="utf-8"))
    raw["mc"]["convergence"] = {
        "initial_replicas": 2,
        "replica_increment": 1,
        "maximum_replicas": 3,
        "maximum_attempts": 2,
        "time_extension_factor": 2,
        "maximum_particle_barriers_per_replica": 10,
        "maximum_total_particle_barriers": 100,
    }
    workflow_path.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")

    with pytest.raises(
        RuntimeError,
        match="maximum_particle_barriers_per_replica",
    ):
        load_workflow(workflow_path)


def test_v5_sampling_plan_provenance_is_not_promoted_to_v6() -> None:
    legacy_plan = [_sampling_row(20.0)]

    with pytest.raises(WorkflowSchemaError, match="must contain exactly"):
        canonical_mc_sampling_plan_json(legacy_plan)


def test_mc_sampling_plan_drives_jobs_and_is_immutable_database_provenance(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    workflow_path = _parallel_mc_workflow(tmp_path, workers=1)
    _write_sampling_plan(
        workflow_path,
        [
            _sampling_row(
                40.0,
                transport_correlation_lag_barriers=256,
                transport_estimator="paired_field_parity",
            )
        ],
    )
    calls: list[dict[str, object]] = []

    def record_config(config: object, **_kwargs: object) -> object:
        mc = config.solvers.monte_carlo
        calls.append(
            {
                "e_over_n_Td": config.run.e_over_n_Td[0],
                "case_prefix": config.run.case_prefix,
                "particles": mc.particles,
                "warmup_collisions": mc.warmup_collisions,
                "max_collisions": mc.max_collisions,
                "tail_max_collisions": mc.tail_max_collisions,
                "seed": mc.seed,
                "transport_correlation_lag_barriers": (
                    mc.transport_correlation_lag_barriers
                ),
                "transport_estimator": mc.transport_estimator,
            }
        )
        return SimpleNamespace(cases=())

    monkeypatch.setattr(sweep_module, "run", record_config)
    monkeypatch.setattr(
        sweep_module,
        "aggregate_workflow_results",
        lambda *_args, **_kwargs: None,
    )

    run_sweep(workflow_path)

    assert len(calls) == 5
    assert [call["e_over_n_Td"] for call in calls] == [
        20.0,
        40.0,
        20.0,
        40.0,
        40.0,
    ]
    by_anchor = {
        e_over_n: [call for call in calls if call["e_over_n_Td"] == e_over_n]
        for e_over_n in (20.0, 40.0)
    }
    assert len(by_anchor[20.0]) == 2
    assert {
        (call["particles"], call["warmup_collisions"], call["max_collisions"])
        for call in by_anchor[20.0]
    } == {(2, 0, 10)}
    assert {call["transport_correlation_lag_barriers"] for call in by_anchor[20.0]} == {
        64
    }
    assert {call["transport_estimator"] for call in by_anchor[20.0]} == {"single_field"}
    assert {call["tail_max_collisions"] for call in by_anchor[20.0]} == {None}
    assert len(by_anchor[40.0]) == 3
    assert {
        (call["particles"], call["warmup_collisions"], call["max_collisions"])
        for call in by_anchor[40.0]
    } == {(3, 4, 20)}
    assert {call["transport_correlation_lag_barriers"] for call in by_anchor[40.0]} == {
        256
    }
    assert {call["transport_estimator"] for call in by_anchor[40.0]} == {
        "paired_field_parity"
    }
    assert {call["tail_max_collisions"] for call in by_anchor[40.0]} == {None}
    assert len({call["seed"] for call in calls}) == 5

    database = tmp_path / "outputs" / "swarm.sqlite"
    with sqlite3.connect(database) as connection:
        metadata = dict(connection.execute("SELECT key, value FROM metadata"))
    assert json.loads(metadata[MC_SAMPLING_PLAN_METADATA_KEY]) == [
        {
            "e_over_n_Td": 20.0,
            "particles": 2,
            "warmup_collisions": 0,
            "max_collisions": 10,
            "tail_max_collisions": 0,
            "replicas": 2,
            "transport_correlation_lag_barriers": 64,
            "transport_estimator": "single_field",
        },
        {
            "e_over_n_Td": 40.0,
            "particles": 3,
            "warmup_collisions": 4,
            "max_collisions": 20,
            "tail_max_collisions": 0,
            "replicas": 3,
            "transport_correlation_lag_barriers": 256,
            "transport_estimator": "paired_field_parity",
        },
    ]
    compute_budget = json.loads(metadata["mc_nominal_compute_budget_json"])
    assert compute_budget["maximum_per_replica_particle_barriers"] == 144
    assert compute_budget["total_particle_barriers"] == 472
    assert compute_budget["includes_conditional_tail_upper_bound"] is False
    assert all(
        entry["tail_phase_is_conditional"] is False
        for entry in compute_budget["entries"]
    )

    _write_sampling_plan(workflow_path, [_sampling_row(40.0, particles=4)])
    with pytest.raises(ValueError, match="workflow_config_sha256"):
        run_sweep(workflow_path)


@pytest.mark.mc
def test_parallel_mc_sweep_writes_results_in_parent_database(tmp_path: Path) -> None:
    workflow_path = _parallel_mc_workflow(tmp_path, workers=2)
    workflow = load_workflow(workflow_path)
    assert workflow.mc_workers == 2

    summary = run_sweep(workflow_path)

    assert summary.cases_written == 4
    with sqlite3.connect(summary.database_path) as connection:
        metadata = dict(connection.execute("SELECT key, value FROM metadata"))
        rows = connection.execute(
            """
            SELECT mixture_id, e_over_n_Td, replicate, diagnostics_json
            FROM cases
            WHERE solver = 'monte_carlo'
            ORDER BY e_over_n_Td, replicate
            """
        ).fetchall()
    assert metadata["mc_transport_estimator_schema_version"] == (
        DIRECT_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
    )
    assert metadata["mc_eedf_estimator_schema_version"] == (
        MC_EEDF_ESTIMATOR_SCHEMA_VERSION
    )
    assert metadata[sweep_module.MC_SEED_DERIVATION_METADATA_KEY] == (
        MC_SEED_DERIVATION_SCHEMA_VERSION
    )
    assert metadata[sweep_module.MC_SOLVER_SOURCE_METADATA_KEY] == (
        monte_carlo_source_sha256()
    )
    assert provenance_hash_manifest(metadata)[
        sweep_module.MC_SOLVER_SOURCE_METADATA_KEY
    ] == monte_carlo_source_sha256()
    assert provenance_hash_manifest(metadata)[
        sweep_module.MC_SEED_DERIVATION_METADATA_KEY
    ] == MC_SEED_DERIVATION_SCHEMA_VERSION
    assert json.loads(metadata[MC_SAMPLING_PLAN_METADATA_KEY]) == [
        {
            "e_over_n_Td": 20.0,
            "particles": 2,
            "warmup_collisions": 0,
            "max_collisions": 10,
            "tail_max_collisions": 0,
            "replicas": 2,
            "transport_correlation_lag_barriers": 64,
            "transport_estimator": "single_field",
        },
        {
            "e_over_n_Td": 40.0,
            "particles": 2,
            "warmup_collisions": 0,
            "max_collisions": 10,
            "tail_max_collisions": 0,
            "replicas": 2,
            "transport_correlation_lag_barriers": 64,
            "transport_estimator": "single_field",
        },
    ]
    compute_budget = json.loads(metadata["mc_nominal_compute_budget_json"])
    assert compute_budget["maximum_per_replica_particle_barriers"] == 20
    assert compute_budget["total_particle_barriers"] == 80
    assert compute_budget["includes_conditional_tail_upper_bound"] is False
    assert "mc_tail_estimator_schema_version" not in metadata
    assert [(row[0], row[1], row[2]) for row in rows] == [
        (0, 20.0, 0),
        (0, 20.0, 1),
        (0, 40.0, 0),
        (0, 40.0, 1),
    ]
    for *_key, diagnostics_json in rows:
        diagnostics = json.loads(diagnostics_json)
        transport = diagnostics["internal_monte_carlo_transport"]
        assert transport["estimator_schema_version"] == (
            DIRECT_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
        )
        assert set(transport["mc_run_provenance"]) >= {
            "seed",
            "particles",
            "warmup_collisions",
            "production_collisions",
            "eedf_estimator_schema_version",
        }
        assert transport["mc_run_provenance"][
            "eedf_estimator_schema_version"
        ] == MC_EEDF_ESTIMATOR_SCHEMA_VERSION
        assert transport["mc_run_provenance"]["case_seed_derivation"] == (
            MC_CASE_SEED_DERIVATION
        )
        assert transport["mc_run_provenance"]["solver_source_sha256"] == (
            metadata[sweep_module.MC_SOLVER_SOURCE_METADATA_KEY]
        )
        rate_diagnostics = diagnostics["internal_monte_carlo_reaction_rates"]
        assert rate_diagnostics["estimator"] == "trajectory_time_average_sigma_v"
        assert all(
            row["batch_means_blocks_requested"] >= 20
            for row in rate_diagnostics["rates"]
        )


@pytest.mark.mc
def test_parallel_mc_resume_executes_only_missing_job(tmp_path: Path) -> None:
    workflow_path = _parallel_mc_workflow(tmp_path, workers=2)
    first = run_sweep(workflow_path)
    assert first.cases_written == 4

    with sqlite3.connect(first.database_path) as connection:
        key = (0, "monte_carlo", 40.0, 1)
        for table in ("rates", "eedf_bins"):
            connection.execute(
                f"""
                DELETE FROM {table}
                WHERE mixture_id = ? AND solver = ? AND e_over_n_Td = ?
                  AND replicate = ?
                """,
                key,
            )
        connection.execute(
            """
            DELETE FROM cases
            WHERE mixture_id = ? AND solver = ? AND e_over_n_Td = ?
              AND replicate = ?
            """,
            key,
        )
        connection.execute(
            """
            UPDATE cases SET case_id = 'preserved_existing_case'
            WHERE mixture_id = 0 AND solver = 'monte_carlo'
              AND e_over_n_Td = 20.0 AND replicate = 0
            """
        )
        connection.commit()

    resumed = run_sweep(workflow_path)

    assert resumed.cases_written == 1
    with sqlite3.connect(first.database_path) as connection:
        assert connection.execute("SELECT COUNT(*) FROM cases").fetchone() == (4,)
        assert connection.execute(
            """
            SELECT case_id FROM cases
            WHERE mixture_id = 0 AND solver = 'monte_carlo'
              AND e_over_n_Td = 20.0 AND replicate = 0
            """
        ).fetchone() == ("preserved_existing_case",)
        assert connection.execute(
            "SELECT COUNT(*) FROM aggregate_quality"
        ).fetchone() == (2,)

    with sqlite3.connect(first.database_path) as connection:
        connection.execute(
            "DELETE FROM metadata WHERE key = 'mc_transport_estimator_schema_version'"
        )
        connection.commit()
    with pytest.raises(
        WorkflowSchemaError, match="missing mc_transport_estimator_schema_version"
    ):
        run_sweep(workflow_path)

    with sqlite3.connect(first.database_path) as connection:
        connection.execute(
            "INSERT OR REPLACE INTO metadata(key, value) VALUES (?, ?)",
            (
                "mc_transport_estimator_schema_version",
                DIRECT_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
            ),
        )
        connection.execute(
            "DELETE FROM metadata WHERE key = 'mc_eedf_estimator_schema_version'"
        )
        connection.commit()
    with pytest.raises(
        WorkflowSchemaError, match="missing mc_eedf_estimator_schema_version"
    ):
        run_sweep(workflow_path)

    with sqlite3.connect(first.database_path) as connection:
        connection.execute(
            "INSERT OR REPLACE INTO metadata(key, value) VALUES (?, ?)",
            (
                "mc_eedf_estimator_schema_version",
                MC_EEDF_ESTIMATOR_SCHEMA_VERSION,
            ),
        )
        connection.execute(
            "DELETE FROM metadata WHERE key = ?",
            (sweep_module.MC_SEED_DERIVATION_METADATA_KEY,),
        )
        connection.commit()
    with pytest.raises(
        WorkflowSchemaError, match="missing mc_seed_derivation_schema_version"
    ):
        run_sweep(workflow_path)

    with sqlite3.connect(first.database_path) as connection:
        connection.execute(
            "INSERT OR REPLACE INTO metadata(key, value) VALUES (?, ?)",
            (
                sweep_module.MC_SEED_DERIVATION_METADATA_KEY,
                MC_SEED_DERIVATION_SCHEMA_VERSION,
            ),
        )
        connection.execute(
            "DELETE FROM metadata WHERE key = ?",
            (MC_SAMPLING_PLAN_METADATA_KEY,),
        )
        connection.commit()
    with pytest.raises(WorkflowSchemaError, match="missing mc_sampling_plan_json"):
        run_sweep(workflow_path)


@pytest.mark.mc
def test_mc_reuse_database_imports_exact_raw_jobs_without_recalculation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source_workflow = _parallel_mc_workflow(tmp_path / "source", workers=1)
    source_summary = run_sweep(source_workflow)
    assert source_summary.cases_written == 4

    target_workflow = _parallel_mc_workflow(tmp_path / "target", workers=2)
    raw = yaml.safe_load(target_workflow.read_text(encoding="utf-8"))
    raw["base_config"] = source_workflow.with_name("mc_base.yaml").as_posix()
    raw["mc"]["reuse_database"] = source_summary.database_path.as_posix()
    target_workflow.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")

    def fail_if_recalculated(_job: object) -> object:
        raise AssertionError("an exact reusable MC case was recalculated")

    monkeypatch.setattr(sweep_module, "_execute_monte_carlo_job", fail_if_recalculated)
    target_summary = run_sweep(target_workflow)

    assert target_summary.cases_written == 4
    with sqlite3.connect(source_summary.database_path) as source_connection:
        source_rows = source_connection.execute(
            """
            SELECT e_over_n_Td, replicate, mean_energy_eV, diagnostics_json
            FROM cases ORDER BY e_over_n_Td, replicate
            """
        ).fetchall()
        source_rate_count = source_connection.execute(
            "SELECT COUNT(*) FROM rates"
        ).fetchone()
    with sqlite3.connect(target_summary.database_path) as target_connection:
        target_rows = target_connection.execute(
            """
            SELECT e_over_n_Td, replicate, mean_energy_eV, diagnostics_json
            FROM cases ORDER BY e_over_n_Td, replicate
            """
        ).fetchall()
        assert (
            target_connection.execute("SELECT COUNT(*) FROM rates").fetchone()
            == source_rate_count
        )
        assert (
            target_connection.execute("SELECT COUNT(*) FROM eedf_bins").fetchone()[0]
            > 0
        )
    assert target_rows == source_rows

    resumed = run_sweep(target_workflow)
    assert resumed.cases_written == 0


@pytest.mark.mc
def test_mc_reuse_database_requires_matching_seed_derivation(
    tmp_path: Path,
) -> None:
    source_workflow = _parallel_mc_workflow(tmp_path / "source_seed", workers=1)
    source_summary = run_sweep(source_workflow)

    with sqlite3.connect(source_summary.database_path) as connection:
        connection.execute(
            "UPDATE metadata SET value = 'obsolete' WHERE key = ?",
            (sweep_module.MC_SEED_DERIVATION_METADATA_KEY,),
        )
        connection.commit()

    target_workflow = _parallel_mc_workflow(tmp_path / "target_seed", workers=1)
    raw = yaml.safe_load(target_workflow.read_text(encoding="utf-8"))
    raw["base_config"] = source_workflow.with_name("mc_base.yaml").as_posix()
    raw["mc"]["reuse_database"] = source_summary.database_path.as_posix()
    target_workflow.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")

    with pytest.raises(
        ValueError, match="different mc_seed_derivation_schema_version"
    ):
        run_sweep(target_workflow)


@pytest.mark.mc
def test_tail_metadata_upgrade_preserves_complete_current_case_evidence(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source_workflow = _parallel_mc_workflow(
        tmp_path / "source_tail_upgrade",
        workers=1,
        population_model="weighted_branching",
    )
    source_base = source_workflow.with_name("mc_base.yaml")
    base_raw = yaml.safe_load(source_base.read_text(encoding="utf-8"))
    base_raw["physics"]["energy_grid_policy"]["threshold_refinement"] = False
    base_raw["solvers"]["monte_carlo"]["tail_max_collisions"] = 10
    base_raw["solvers"]["monte_carlo"]["tail_rate_rse_trigger"] = 0.2
    source_base.write_text(
        yaml.safe_dump(base_raw, sort_keys=False),
        encoding="utf-8",
    )
    source_summary = run_sweep(source_workflow)

    with sqlite3.connect(source_summary.database_path) as connection:
        connection.execute(
            "DELETE FROM metadata WHERE key = ?",
            (sweep_module.MC_TAIL_ESTIMATOR_METADATA_KEY,),
        )
        source_diagnostics = connection.execute(
            "SELECT diagnostics_json FROM cases ORDER BY e_over_n_Td, replicate"
        ).fetchall()
        connection.commit()

    target_workflow = _parallel_mc_workflow(
        tmp_path / "target_tail_upgrade",
        workers=1,
        population_model="weighted_branching",
    )
    target_raw = yaml.safe_load(target_workflow.read_text(encoding="utf-8"))
    target_raw["base_config"] = source_base.as_posix()
    target_raw["mc"]["reuse_database"] = source_summary.database_path.as_posix()
    target_workflow.write_text(
        yaml.safe_dump(target_raw, sort_keys=False),
        encoding="utf-8",
    )

    monkeypatch.setattr(
        sweep_module,
        "_execute_monte_carlo_job",
        lambda _job: (_ for _ in ()).throw(
            AssertionError("a complete current-tail case was recalculated")
        ),
    )
    target_summary = run_sweep(target_workflow)

    assert target_summary.cases_written == 4
    with sqlite3.connect(target_summary.database_path) as connection:
        metadata = dict(connection.execute("SELECT key, value FROM metadata"))
        diagnostics = [
            json.loads(row[0])
            for row in connection.execute(
                "SELECT diagnostics_json FROM cases ORDER BY e_over_n_Td, replicate"
            )
        ]
    assert metadata[sweep_module.MC_TAIL_ESTIMATOR_METADATA_KEY] == (
        WEIGHTED_MC_TAIL_ESTIMATOR_SCHEMA_VERSION
    )
    assert [json.dumps(item, sort_keys=True) for item in diagnostics] == [
        row[0] for row in source_diagnostics
    ]


@pytest.mark.mc
def test_tail_metadata_upgrade_recomputes_incomplete_case_evidence(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source_workflow = _parallel_mc_workflow(
        tmp_path / "source_incomplete_tail",
        workers=1,
        population_model="weighted_branching",
    )
    source_base = source_workflow.with_name("mc_base.yaml")
    base_raw = yaml.safe_load(source_base.read_text(encoding="utf-8"))
    base_raw["physics"]["energy_grid_policy"]["threshold_refinement"] = False
    base_raw["solvers"]["monte_carlo"]["tail_max_collisions"] = 10
    base_raw["solvers"]["monte_carlo"]["tail_rate_rse_trigger"] = 0.2
    source_base.write_text(
        yaml.safe_dump(base_raw, sort_keys=False),
        encoding="utf-8",
    )
    source_summary = run_sweep(source_workflow)

    missing_fields = (
        ("configured_max_collisions",),
        ("executed_collisions",),
        ("activation_rate_rse_trigger",),
        (
            "reported_eedf_includes_main_production",
            "reported_rates_include_main_production",
        ),
    )
    with sqlite3.connect(source_summary.database_path) as connection:
        connection.execute(
            "DELETE FROM metadata WHERE key = ?",
            (sweep_module.MC_TAIL_ESTIMATOR_METADATA_KEY,),
        )
        rows = connection.execute(
            """
            SELECT mixture_id, solver, e_over_n_Td, replicate, diagnostics_json
            FROM cases ORDER BY e_over_n_Td, replicate
            """
        ).fetchall()
        assert len(rows) == len(missing_fields)
        for row, fields in zip(rows, missing_fields, strict=True):
            mixture_id, solver, field, replicate, encoded = row
            diagnostics = json.loads(encoded)
            tail = diagnostics["internal_monte_carlo_transport"]["tail_sampling"]
            for field_name in fields:
                tail.pop(field_name)
            connection.execute(
                """
                UPDATE cases SET diagnostics_json = ?
                WHERE mixture_id = ? AND solver = ? AND e_over_n_Td = ?
                  AND replicate = ?
                """,
                (
                    json.dumps(diagnostics, sort_keys=True),
                    mixture_id,
                    solver,
                    field,
                    replicate,
                ),
            )
        connection.commit()

    target_workflow = _parallel_mc_workflow(
        tmp_path / "target_incomplete_tail",
        workers=1,
        population_model="weighted_branching",
    )
    target_raw = yaml.safe_load(target_workflow.read_text(encoding="utf-8"))
    target_raw["base_config"] = source_base.as_posix()
    target_raw["mc"]["reuse_database"] = source_summary.database_path.as_posix()
    target_workflow.write_text(
        yaml.safe_dump(target_raw, sort_keys=False),
        encoding="utf-8",
    )

    original_execute = sweep_module._execute_monte_carlo_job
    recalculated: list[tuple[float, int]] = []

    def record_recalculation(job: object) -> object:
        recalculated.append((job.e_over_n_Td, job.replicate))
        return original_execute(job)

    monkeypatch.setattr(
        sweep_module,
        "_execute_monte_carlo_job",
        record_recalculation,
    )
    target_summary = run_sweep(target_workflow)

    assert target_summary.cases_written == 4
    assert sorted(recalculated) == [
        (20.0, 0),
        (20.0, 1),
        (40.0, 0),
        (40.0, 1),
    ]
    with sqlite3.connect(target_summary.database_path) as connection:
        repaired = [
            json.loads(row[0])["internal_monte_carlo_transport"]["tail_sampling"]
            for row in connection.execute(
                "SELECT diagnostics_json FROM cases ORDER BY e_over_n_Td, replicate"
            )
        ]
    assert all(
        tail["configured_max_collisions"] == 10
        and tail["executed_collisions"] in {0, 10}
        and tail["activation_rate_rse_trigger"] == pytest.approx(0.2)
        and isinstance(tail["reported_eedf_includes_main_production"], bool)
        and isinstance(tail["reported_rates_include_main_production"], bool)
        for tail in repaired
    )


@pytest.mark.mc
def test_weighted_branching_sweep_freezes_distinct_transport_schema(
    tmp_path: Path,
) -> None:
    workflow_path = _parallel_mc_workflow(
        tmp_path,
        workers=1,
        population_model="weighted_branching",
    )

    summary = run_sweep(workflow_path)

    with sqlite3.connect(summary.database_path) as connection:
        metadata = dict(connection.execute("SELECT key, value FROM metadata"))
        diagnostics = [
            json.loads(row[0])
            for row in connection.execute(
                "SELECT diagnostics_json FROM cases WHERE solver = 'monte_carlo'"
            )
        ]
    assert metadata["mc_transport_estimator_schema_version"] == (
        WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
    )
    assert all(
        payload["internal_monte_carlo_transport"]["estimator_schema_version"]
        == WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
        for payload in diagnostics
    )
