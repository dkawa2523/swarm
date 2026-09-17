from __future__ import annotations

import json
from pathlib import Path
import sqlite3
from types import SimpleNamespace

import pytest
import yaml

from electron_swarm import load_config
from product_helpers import base_product_config, write_config
import swarm_workflow.campaign.sweep as sweep_module
from swarm_workflow.campaign.config import MixtureSpec, load_workflow
from swarm_workflow.campaign.sweep import run_sweep


def _deterministic_workflow(
    tmp_path: Path,
    *,
    solver: str = "two_term",
    fields: tuple[float, ...] = (10.0, 20.0, 30.0, 40.0),
    workers: int | None = None,
    memory_budget_mb: int | None = None,
    propagator_memory_mb: int = 1024,
) -> Path:
    tmp_path.mkdir(parents=True, exist_ok=True)
    data = base_product_config(tmp_path, [solver])
    data["run"]["e_over_n_Td"] = [fields[0]]
    if solver == "propagator":
        data["solvers"]["propagator"] = {
            "energy_cells": 64,
            "polar_cells": 8,
            "max_iterations": 20,
            "convergence_tolerance": 1.0e-7,
            "max_memory_mb": propagator_memory_mb,
        }
    base = write_config(tmp_path, data, "base.yaml")
    raw: dict[str, object] = {
        "base_config": base.name,
        "database": "outputs/swarm.sqlite",
        "e_over_n_Td": list(fields),
        "mixtures": [{"Ar": 1.0}],
    }
    if workers is not None or memory_budget_mb is not None:
        deterministic: dict[str, int] = {}
        if workers is not None:
            deterministic["workers"] = workers
        if memory_budget_mb is not None:
            deterministic["global_memory_budget_mb"] = memory_budget_mb
        raw["execution"] = {"deterministic": deterministic}
    path = tmp_path / "workflow.yaml"
    path.write_text(
        yaml.safe_dump(raw, sort_keys=False),
        encoding="utf-8",
    )
    return path


def test_deterministic_execution_defaults_and_validation(tmp_path: Path) -> None:
    default_path = _deterministic_workflow(tmp_path / "default")
    default = load_workflow(default_path).deterministic_execution
    assert default.workers == 1
    assert default.global_memory_budget_mb is None
    assert default.configured is False

    missing_budget = _deterministic_workflow(
        tmp_path / "missing-budget",
        workers=2,
    )
    with pytest.raises(ValueError, match="global_memory_budget_mb.*required"):
        load_workflow(missing_budget)

    invalid_workers = _deterministic_workflow(
        tmp_path / "invalid-workers",
        workers=0,
        memory_budget_mb=2048,
    )
    with pytest.raises(ValueError, match=r"deterministic\.workers"):
        load_workflow(invalid_workers)

    invalid_memory = _deterministic_workflow(
        tmp_path / "invalid-memory",
        workers=1,
        memory_budget_mb=0,
    )
    with pytest.raises(ValueError, match="global_memory_budget_mb"):
        load_workflow(invalid_memory)


def test_worker_cap_uses_cpu_global_budget_and_propagator_memory_ceiling(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    workflow_path = _deterministic_workflow(
        tmp_path,
        solver="propagator",
        fields=(10.0, 20.0, 30.0, 40.0, 50.0, 60.0),
        workers=8,
        memory_budget_mb=5000,
        propagator_memory_mb=2048,
    )
    workflow = load_workflow(workflow_path)
    base = load_config(workflow.base_config_path)
    monkeypatch.setattr(sweep_module, "_available_cpu_count", lambda: 6)

    plan = sweep_module._resolve_deterministic_execution(
        workflow,
        base,
        ["propagator"],
    )

    assert plan.solver_memory_reservations_mb == (("propagator", 2304),)
    assert plan.worker_memory_reservation_mb == 2304
    assert plan.effective_worker_cap == 2

    raw = yaml.safe_load(workflow_path.read_text(encoding="utf-8"))
    raw["execution"]["deterministic"]["global_memory_budget_mb"] = 2000
    workflow_path.write_text(
        yaml.safe_dump(raw, sort_keys=False),
        encoding="utf-8",
    )
    undersized = load_workflow(workflow_path)
    with pytest.raises(ValueError, match="below.*per-worker reservation"):
        sweep_module._resolve_deterministic_execution(
            undersized,
            load_config(undersized.base_config_path),
            ["propagator"],
        )


def test_small_propagator_ceiling_is_not_replaced_by_generic_reservation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    workflow_path = _deterministic_workflow(
        tmp_path,
        solver="propagator",
        fields=tuple(float(index) for index in range(1, 9)),
        workers=8,
        memory_budget_mb=4096,
        propagator_memory_mb=256,
    )
    workflow = load_workflow(workflow_path)
    base = load_config(workflow.base_config_path)
    monkeypatch.setattr(sweep_module, "_available_cpu_count", lambda: 20)

    plan = sweep_module._resolve_deterministic_execution(
        workflow,
        base,
        ["propagator"],
    )

    assert plan.solver_memory_reservations_mb == (("propagator", 512),)
    assert plan.worker_memory_reservation_mb == 512
    assert plan.effective_worker_cap == 8


def test_parallel_job_builder_keeps_chunks_contiguous_and_ordered(
    tmp_path: Path,
) -> None:
    mixture = MixtureSpec(mixture_id=0, fractions={"Ar": 1.0})
    group = sweep_module._DeterministicSweepGroup(
        base_config_path=tmp_path / "base.yaml",
        mixture=mixture,
        solver_id="propagator",
        indexed_e_over_n_Td=(
            (0, 10.0),
            (1, 20.0),
            (3, 40.0),
            (4, 50.0),
        ),
        propagator_source_sha256="a" * 64,
    )

    jobs = sweep_module._build_deterministic_jobs([group], workers=3)

    assert [job.case_indices for job in jobs] == [(0,), (1,), (3, 4)]
    assert [job.e_over_n_Td for job in jobs] == [
        (10.0,),
        (20.0,),
        (40.0, 50.0),
    ]
    assert all(job.canonicalize_case_ids for job in jobs)

    sequential = sweep_module._build_deterministic_jobs([group], workers=1)
    assert [job.e_over_n_Td for job in sequential] == [
        (10.0, 20.0, 40.0, 50.0)
    ]
    assert sequential[0].canonicalize_case_ids is False


def test_deterministic_worker_runs_one_solver_chunk_and_canonicalizes_ids(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    workflow_path = _deterministic_workflow(tmp_path)
    workflow = load_workflow(workflow_path)
    mixture = workflow.mixtures[0]
    observed: list[tuple[tuple[float, ...], tuple[str, ...]]] = []

    def fake_run(config: object, **_kwargs: object) -> object:
        observed.append(
            (
                tuple(config.run.e_over_n_Td),
                tuple(item.id for item in config.run.solvers),
            )
        )
        cases = []
        for index, value in enumerate(config.run.e_over_n_Td):
            cases.append(
                SimpleNamespace(
                    solver="two_term",
                    e_over_n_Td=value,
                    case_id=f"local_{index}",
                    rates=[SimpleNamespace(case_id=f"local_{index}")],
                    diagnostics={"two_term": {"run_label": f"local_{index}"}},
                )
            )
        return SimpleNamespace(cases=cases)

    monkeypatch.setattr(sweep_module, "run", fake_run)
    result = sweep_module._execute_deterministic_job(
        sweep_module._DeterministicSweepJob(
            base_config_path=workflow.base_config_path,
            mixture=mixture,
            solver_id="two_term",
            case_indices=(2, 3),
            e_over_n_Td=(30.0, 40.0),
            canonicalize_case_ids=True,
            propagator_source_sha256=None,
        )
    )

    assert observed == [((30.0, 40.0), ("two_term",))]
    assert [case.case_id for case in result.cases] == [
        "Prod_m0000_r0000_0002",
        "Prod_m0000_r0000_0003",
    ]
    assert [case.rates[0].case_id for case in result.cases] == [
        "Prod_m0000_r0000_0002",
        "Prod_m0000_r0000_0003",
    ]


def test_parallel_results_preserve_submission_order_and_propagate_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cancelled: list[object] = []

    class FakeFuture:
        def __init__(self, value: object) -> None:
            self.value = value

        def result(self) -> object:
            if isinstance(self.value, BaseException):
                raise self.value
            return self.value

        def cancel(self) -> None:
            cancelled.append(self.value)

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
    assert list(
        sweep_module._iter_deterministic_results([1, 2, 3], workers=2)
    ) == [1, 2, 3]

    class FailingExecutor(FakeExecutor):
        def submit(self, _function: object, job: object) -> FakeFuture:
            if job == 2:
                return FakeFuture(RuntimeError("worker failed"))
            return FakeFuture(job)

    monkeypatch.setattr(sweep_module, "ProcessPoolExecutor", FailingExecutor)
    with pytest.raises(RuntimeError, match="worker failed"):
        list(sweep_module._iter_deterministic_results([1, 2, 3], workers=2))
    assert cancelled[-3] == 1
    assert isinstance(cancelled[-2], RuntimeError)
    assert str(cancelled[-2]) == "worker failed"
    assert cancelled[-1] == 3


def test_parallel_sweep_writes_canonical_order_and_execution_provenance(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    workflow_path = _deterministic_workflow(
        tmp_path,
        workers=2,
        memory_budget_mb=2048,
    )
    calls: list[tuple[float, ...]] = []
    original_run = sweep_module.run

    def record_run(config: object, **kwargs: object) -> object:
        calls.append(tuple(config.run.e_over_n_Td))
        return original_run(config, **kwargs)

    class ImmediateFuture:
        def __init__(self, function: object, job: object) -> None:
            try:
                self.value = function(job)
                self.error: BaseException | None = None
            except BaseException as exc:  # pragma: no cover - defensive parity
                self.value = None
                self.error = exc

        def result(self) -> object:
            if self.error is not None:
                raise self.error
            return self.value

        def cancel(self) -> None:
            return None

    class ImmediateExecutor:
        def __init__(self, **_kwargs: object) -> None:
            pass

        def __enter__(self) -> "ImmediateExecutor":
            return self

        def __exit__(self, *_args: object) -> None:
            return None

        def submit(self, function: object, job: object) -> ImmediateFuture:
            assert not hasattr(job, "database_path")
            return ImmediateFuture(function, job)

    monkeypatch.setattr(sweep_module, "_available_cpu_count", lambda: 4)
    monkeypatch.setattr(sweep_module, "ProcessPoolExecutor", ImmediateExecutor)
    monkeypatch.setattr(sweep_module, "run", record_run)

    summary = run_sweep(workflow_path)

    assert calls == [(10.0, 20.0), (30.0, 40.0)]
    assert summary.cases_written == 4
    with sqlite3.connect(summary.database_path) as connection:
        rows = connection.execute(
            "SELECT e_over_n_Td, case_id FROM cases ORDER BY rowid"
        ).fetchall()
        metadata = dict(connection.execute("SELECT key, value FROM metadata"))
    assert rows == [
        (10.0, "Prod_m0000_r0000_0000"),
        (20.0, "Prod_m0000_r0000_0001"),
        (30.0, "Prod_m0000_r0000_0002"),
        (40.0, "Prod_m0000_r0000_0003"),
    ]
    execution = json.loads(metadata["deterministic_execution_json"])
    assert execution == {
        "schema": "swarm.deterministic_execution.v1",
        "requested_workers": 2,
        "effective_worker_cap": 2,
        "cpu_limit": 4,
        "global_memory_budget_mb": 2048,
        "worker_memory_reservation_mb": 1024,
        "solver_memory_reservations_mb": {"two_term": 1024},
        "chunk_policy": "contiguous_workflow_order_warm_start",
        "result_order": "workflow_mixture_solver_field_order",
        "database_writer": "parent_process_only",
    }

    resumed = run_sweep(workflow_path)
    assert resumed.cases_written == 0
    assert calls == [(10.0, 20.0), (30.0, 40.0)]


def test_spawned_deterministic_workers_return_results_to_parent_database(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    workflow_path = _deterministic_workflow(
        tmp_path,
        fields=(20.0, 40.0),
        workers=2,
        memory_budget_mb=2048,
    )
    monkeypatch.setattr(sweep_module, "_available_cpu_count", lambda: 2)

    summary = run_sweep(workflow_path)

    assert summary.cases_written == 2
    with sqlite3.connect(summary.database_path) as connection:
        assert connection.execute("SELECT COUNT(*) FROM cases").fetchone() == (2,)
        assert connection.execute(
            "SELECT COUNT(*) FROM eedf_bins"
        ).fetchone()[0] > 0
