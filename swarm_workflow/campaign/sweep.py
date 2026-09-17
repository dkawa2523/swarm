"""External E/N and gas-mixture workflow sweeps."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from copy import deepcopy
from dataclasses import asdict, dataclass
import json
import math
import multiprocessing
import os
from pathlib import Path
import sqlite3
import tempfile
from typing import Iterator

import yaml

from electron_swarm import (
    GasComponent,
    RequestedSolverConfig,
    SolverId,
    SwarmConfig,
    SwarmCaseResult,
    load_config,
    run,
)
import electron_swarm.solvers.monte_carlo.evidence as _mc_evidence
from electron_swarm.solvers.monte_carlo.source_identity import (
    monte_carlo_source_sha256,
)
from electron_swarm.orchestration.plan import validate_unique_solver_ids
from swarm_workflow.quality.propagator_source import (
    propagator_qualification_source_fingerprint,
)

from .aggregate import (
    aggregate_workflow_results,
    stable_mc_seed,
)
from .mean_energy_support import (
    continue_deterministic_mean_energy_support,
    propose_upper_mean_energy_extension,
)
from .store import (
    MC_SAMPLING_PLAN_METADATA_KEY,
    WorkflowStore,
    mc_sampling_plan_provenance,
    validate_workflow_schema,
)
from ..selection import (
    PHYSICAL_CONTEXT_KEY,
    file_sha256,
    physical_context,
    read_selection,
)
from ..quality.monte_carlo.policy import (
    MonteCarloPolicyError,
    SamplingPlanEntry,
    campaign_provenance,
    parse_convergence_policy,
    sampling_budget_provenance,
    validate_sampling_budget,
)
from ..quality.policy import (
    QUALITY_THRESHOLDS_METADATA_KEY,
    quality_thresholds_json,
)
from . import config as _workflow_config
from .provenance import (
    DETERMINISTIC_EXECUTION_METADATA_KEY,
    MC_EEDF_ESTIMATOR_METADATA_KEY,
    MC_SEED_DERIVATION_METADATA_KEY,
    MC_SOLVER_SOURCE_METADATA_KEY,
    MC_TAIL_ESTIMATOR_METADATA_KEY,
    MC_TRANSPORT_ESTIMATOR_METADATA_KEY,
    PROPAGATOR_SOLVER_SOURCE_METADATA_KEY,
    mc_sampling_plan_json,
    resume_provenance,
    workflow_config_sha256,
)
from .repository import read_metadata


# A worker imports NumPy/SciPy and constructs solver-local sparse operators in
# addition to the arrays counted by a solver.  The propagator exposes an
# explicit core-memory ceiling, so reserve that ceiling plus process overhead.
# The other deterministic solvers do not expose a product memory knob; use a
# deliberately conservative 1 GiB process reservation, with an lmax-dependent
# increase for unusually large multi-term systems.
_WORKER_PROCESS_OVERHEAD_MB = 256
_DETERMINISTIC_FALLBACK_RESERVATION_MB = 1024
_REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True, slots=True)
class SweepSummary:
    database_path: Path
    mixtures: int
    e_over_n_values: int
    cases_written: int
    quality_policy_reevaluated: bool


@dataclass(frozen=True, slots=True)
class _MonteCarloSweepJob:
    base_config_path: Path
    mixture: _workflow_config.MixtureSpec
    e_over_n_Td: float
    replicate: int
    seed: int | None
    sampling: SamplingPlanEntry


@dataclass(frozen=True, slots=True)
class _MonteCarloSweepResult:
    mixture_id: int
    replicate: int
    cases: tuple[SwarmCaseResult, ...]


@dataclass(frozen=True, slots=True)
class _DeterministicExecutionPlan:
    requested_workers: int
    effective_worker_cap: int
    cpu_limit: int
    global_memory_budget_mb: int | None
    worker_memory_reservation_mb: int
    solver_memory_reservations_mb: tuple[tuple[SolverId, int], ...]


@dataclass(frozen=True, slots=True)
class _DeterministicSweepGroup:
    base_config_path: Path
    mixture: _workflow_config.MixtureSpec
    solver_id: SolverId
    indexed_e_over_n_Td: tuple[tuple[int, float], ...]
    propagator_source_sha256: str | None


@dataclass(frozen=True, slots=True)
class _DeterministicSweepJob:
    base_config_path: Path
    mixture: _workflow_config.MixtureSpec
    solver_id: SolverId
    case_indices: tuple[int, ...]
    e_over_n_Td: tuple[float, ...]
    canonicalize_case_ids: bool
    propagator_source_sha256: str | None


@dataclass(frozen=True, slots=True)
class _DeterministicSweepResult:
    mixture_id: int
    solver_id: SolverId
    cases: tuple[SwarmCaseResult, ...]


def _enabled_solver_ids(config: SwarmConfig) -> list[SolverId]:
    return [item.id for item in config.run.solvers if item.enabled]


def _require_current_propagator_source(expected_sha256: str | None) -> None:
    if expected_sha256 is None:
        return
    current = str(
        propagator_qualification_source_fingerprint(
            _REPOSITORY_ROOT
        )["sha256"]
    )
    if current != expected_sha256:
        raise RuntimeError("Propagator source changed during deterministic sweep")


def _available_cpu_count() -> int:
    process_cpu_count = getattr(os, "process_cpu_count", None)
    count = process_cpu_count() if process_cpu_count is not None else os.cpu_count()
    return max(int(count or 1), 1)


def _deterministic_solver_memory_reservation_mb(
    config: SwarmConfig,
    solver_id: SolverId,
) -> int:
    if solver_id == "propagator":
        return (
            int(config.solvers.propagator.max_memory_mb) + _WORKER_PROCESS_OVERHEAD_MB
        )
    if solver_id == "multi_term":
        angular_workspace_mb = 4 * (int(config.solvers.multi_term.lmax) + 1) ** 2
        return max(
            _DETERMINISTIC_FALLBACK_RESERVATION_MB,
            _WORKER_PROCESS_OVERHEAD_MB + angular_workspace_mb,
        )
    if solver_id == "two_term":
        return _DETERMINISTIC_FALLBACK_RESERVATION_MB
    raise ValueError(f"{solver_id!r} is not a deterministic solver")


def _resolve_deterministic_execution(
    workflow: _workflow_config.WorkflowConfig,
    base_config: SwarmConfig,
    solver_ids: list[SolverId],
) -> _DeterministicExecutionPlan:
    deterministic_ids = tuple(
        solver_id for solver_id in solver_ids if solver_id != "monte_carlo"
    )
    reservations = tuple(
        (
            solver_id,
            _deterministic_solver_memory_reservation_mb(base_config, solver_id),
        )
        for solver_id in deterministic_ids
    )
    cpu_limit = _available_cpu_count()
    requested = workflow.deterministic_execution.workers
    if not reservations:
        return _DeterministicExecutionPlan(
            requested_workers=requested,
            effective_worker_cap=0,
            cpu_limit=cpu_limit,
            global_memory_budget_mb=(
                workflow.deterministic_execution.global_memory_budget_mb
            ),
            worker_memory_reservation_mb=0,
            solver_memory_reservations_mb=(),
        )

    worker_reservation = max(value for _solver, value in reservations)
    memory_budget = workflow.deterministic_execution.global_memory_budget_mb
    if memory_budget is None:
        memory_limit = 1
    else:
        memory_limit = memory_budget // worker_reservation
        if memory_limit < 1:
            raise ValueError(
                "workflow execution.deterministic.global_memory_budget_mb="
                f"{memory_budget} is below the conservative per-worker "
                f"reservation of {worker_reservation} MiB"
            )
    total_cases = (
        len(workflow.mixtures) * len(workflow.e_over_n_Td) * len(deterministic_ids)
    )
    effective = min(requested, cpu_limit, memory_limit, max(total_cases, 1))
    return _DeterministicExecutionPlan(
        requested_workers=requested,
        effective_worker_cap=max(effective, 1),
        cpu_limit=cpu_limit,
        global_memory_budget_mb=memory_budget,
        worker_memory_reservation_mb=worker_reservation,
        solver_memory_reservations_mb=reservations,
    )


def _deterministic_execution_payload(
    plan: _DeterministicExecutionPlan,
) -> dict[str, object]:
    return {
        "schema": "swarm.deterministic_execution.v1",
        "requested_workers": plan.requested_workers,
        "effective_worker_cap": plan.effective_worker_cap,
        "cpu_limit": plan.cpu_limit,
        "global_memory_budget_mb": plan.global_memory_budget_mb,
        "worker_memory_reservation_mb": plan.worker_memory_reservation_mb,
        "solver_memory_reservations_mb": dict(plan.solver_memory_reservations_mb),
        "chunk_policy": "contiguous_workflow_order_warm_start",
        "result_order": "workflow_mixture_solver_field_order",
        "database_writer": "parent_process_only",
    }


def _case_prefix(
    base: SwarmConfig,
    mixture: _workflow_config.MixtureSpec,
    replicate: int,
) -> str:
    return f"{base.run.case_prefix}_m{mixture.mixture_id:04d}_r{replicate:04d}"


def _clone_config(
    base: SwarmConfig,
    *,
    mixture: _workflow_config.MixtureSpec,
    e_over_n_values: tuple[float, ...],
    solver_ids: list[SolverId],
    replicate: int,
    seed: int | None = None,
) -> SwarmConfig:
    config = deepcopy(base)
    by_species = {gas.species: gas for gas in base.conditions.gas_mixture}
    config.conditions.gas_mixture = [
        GasComponent(
            species=species,
            fraction=mixture.fractions[species],
            mass_amu=by_species[species].mass_amu,
        )
        for species in by_species
    ]
    config.run.e_over_n_Td = list(e_over_n_values)
    config.run.solvers = [
        RequestedSolverConfig(id=solver_id, enabled=True) for solver_id in solver_ids
    ]
    config.run.case_prefix = _case_prefix(base, mixture, replicate)
    if "monte_carlo" in solver_ids and seed is not None:
        config.solvers.monte_carlo.seed = seed
    return config


def _mixture_species_rows(
    base: SwarmConfig,
    mixture: _workflow_config.MixtureSpec,
) -> list[tuple[str, float, float]]:
    mass_by_species = {gas.species: gas.mass_amu for gas in base.conditions.gas_mixture}
    return [
        (species, mixture.fractions[species], mass_by_species[species])
        for species in mass_by_species
    ]


def _extend_workflow_mean_energy_support(
    store: WorkflowStore,
    *,
    workflow: _workflow_config.WorkflowConfig,
    base_config: SwarmConfig,
    solver_id: SolverId,
) -> tuple[SwarmCaseResult, ...]:
    """Append only the bounded deterministic anchors missing from this database."""

    support = workflow.mean_energy_support
    if support is None:
        return ()
    mixture = workflow.mixtures[0]
    rows = store.connection.execute(
        "SELECT e_over_n_Td, mean_energy_eV FROM cases "
        "WHERE mixture_id=? AND solver=? AND replicate=0 "
        "ORDER BY e_over_n_Td",
        (mixture.mixture_id, solver_id),
    ).fetchall()
    fields = [float(row[0]) for row in rows]
    means = [float(row[1]) for row in rows]
    request = propose_upper_mean_energy_extension(
        fields,
        means,
        support.required_max_mean_energy_eV,
        relative_guard=support.relative_guard,
        maximum_field_step_factor=support.maximum_field_step_factor,
    )
    if request is None:
        return ()

    configured_fields = set(workflow.e_over_n_Td)
    existing_extensions = sum(field not in configured_fields for field in fields)
    remaining_steps = support.maximum_steps - existing_extensions
    if remaining_steps <= 0:
        raise ValueError(
            "deterministic mean-energy support remains below its guarded target "
            "after the configured maximum continuation steps"
        )
    if (
        support.maximum_e_over_n_Td is not None
        and support.maximum_e_over_n_Td <= fields[-1]
    ):
        raise ValueError(
            "deterministic mean-energy support remains below its guarded target "
            "at the configured maximum E/N"
        )

    config = _clone_config(
        base_config,
        mixture=mixture,
        e_over_n_values=tuple(fields),
        solver_ids=[solver_id],
        replicate=0,
    )
    extension_cases: list[SwarmCaseResult] = []
    for _step in range(remaining_steps):
        continuation = continue_deterministic_mean_energy_support(
            config,
            solver_id,
            fields,
            means,
            support.required_max_mean_energy_eV,
            relative_guard=support.relative_guard,
            maximum_field_step_factor=support.maximum_field_step_factor,
            maximum_steps=1,
            maximum_e_over_n_Td=support.maximum_e_over_n_Td,
        )
        if not continuation.extension_cases:
            break
        case = continuation.extension_cases[0]
        store.write_case(
            mixture_id=mixture.mixture_id,
            replicate=0,
            case=case,
        )
        extension_cases.append(case)
        fields.append(float(case.e_over_n_Td))
        means.append(float(case.mean_energy_eV))
        if continuation.reached_target:
            return tuple(extension_cases)
        if (
            support.maximum_e_over_n_Td is not None
            and support.maximum_e_over_n_Td <= fields[-1]
        ):
            break
    raise ValueError(
        "deterministic mean-energy support remains below its guarded target "
        "after bounded continuation"
    )


def _contiguous_index_runs(
    items: tuple[tuple[int, float], ...],
) -> list[tuple[tuple[int, float], ...]]:
    if not items:
        return []
    runs: list[list[tuple[int, float]]] = [[items[0]]]
    for item in items[1:]:
        if item[0] == runs[-1][-1][0] + 1:
            runs[-1].append(item)
        else:
            runs.append([item])
    return [tuple(run_items) for run_items in runs]


def _balanced_contiguous_chunks(
    items: tuple[tuple[int, float], ...],
    chunks: int,
) -> list[tuple[tuple[int, float], ...]]:
    chunks = min(max(int(chunks), 1), len(items))
    base_size, remainder = divmod(len(items), chunks)
    result: list[tuple[tuple[int, float], ...]] = []
    start = 0
    for index in range(chunks):
        size = base_size + (1 if index < remainder else 0)
        result.append(items[start : start + size])
        start += size
    return result


def _build_deterministic_jobs(
    groups: list[_DeterministicSweepGroup],
    *,
    workers: int,
) -> list[_DeterministicSweepJob]:
    if workers <= 1:
        return [
            _DeterministicSweepJob(
                base_config_path=group.base_config_path,
                mixture=group.mixture,
                solver_id=group.solver_id,
                case_indices=tuple(
                    index for index, _value in group.indexed_e_over_n_Td
                ),
                e_over_n_Td=tuple(value for _index, value in group.indexed_e_over_n_Td),
                canonicalize_case_ids=False,
                propagator_source_sha256=group.propagator_source_sha256,
            )
            for group in groups
            if group.indexed_e_over_n_Td
        ]

    segments: list[tuple[_DeterministicSweepGroup, tuple[tuple[int, float], ...]]] = []
    for group in groups:
        segments.extend(
            (group, run_items)
            for run_items in _contiguous_index_runs(group.indexed_e_over_n_Td)
        )
    if not segments:
        return []

    allocations = [1] * len(segments)
    target_jobs = min(
        max(int(workers), len(segments)),
        sum(len(items) for _group, items in segments),
    )
    while sum(allocations) < target_jobs:
        candidates = [
            index
            for index, (_group, items) in enumerate(segments)
            if allocations[index] < len(items)
        ]
        if not candidates:
            break
        split_index = max(
            candidates,
            key=lambda index: (
                math.ceil(len(segments[index][1]) / allocations[index]),
                len(segments[index][1]),
                -index,
            ),
        )
        allocations[split_index] += 1

    jobs: list[_DeterministicSweepJob] = []
    for (group, items), chunks in zip(segments, allocations, strict=True):
        for chunk in _balanced_contiguous_chunks(items, chunks):
            jobs.append(
                _DeterministicSweepJob(
                    base_config_path=group.base_config_path,
                    mixture=group.mixture,
                    solver_id=group.solver_id,
                    case_indices=tuple(index for index, _value in chunk),
                    e_over_n_Td=tuple(value for _index, value in chunk),
                    canonicalize_case_ids=True,
                    propagator_source_sha256=group.propagator_source_sha256,
                )
            )
    return jobs


def _execute_deterministic_job(
    job: _DeterministicSweepJob,
) -> _DeterministicSweepResult:
    """Run one same-solver E/N chunk without opening the workflow database."""

    if not job.e_over_n_Td or len(job.e_over_n_Td) != len(job.case_indices):
        raise ValueError("deterministic sweep job has inconsistent E/N indices")
    if job.solver_id == "propagator":
        _require_current_propagator_source(job.propagator_source_sha256)
    base_config = load_config(job.base_config_path)
    config = _clone_config(
        base_config,
        mixture=job.mixture,
        e_over_n_values=job.e_over_n_Td,
        solver_ids=[job.solver_id],
        replicate=0,
    )
    cases = tuple(run(config, write=False).cases)
    if job.solver_id == "propagator":
        _require_current_propagator_source(job.propagator_source_sha256)
    if len(cases) != len(job.e_over_n_Td):
        raise RuntimeError(
            f"deterministic {job.solver_id} chunk returned {len(cases)} cases "
            f"for {len(job.e_over_n_Td)} E/N values"
        )
    prefix = _case_prefix(base_config, job.mixture, 0)
    for case, expected_e_over_n, case_index in zip(
        cases,
        job.e_over_n_Td,
        job.case_indices,
        strict=True,
    ):
        if case.solver != job.solver_id or float(case.e_over_n_Td) != float(
            expected_e_over_n
        ):
            raise RuntimeError(
                "deterministic sweep worker returned cases out of requested order"
            )
        if job.canonicalize_case_ids:
            local_case_id = case.case_id
            case.case_id = f"{prefix}_{case_index:04d}"
            for rate in case.rates:
                rate.case_id = case.case_id
            for solver_diagnostics in case.diagnostics.values():
                if (
                    isinstance(solver_diagnostics, dict)
                    and solver_diagnostics.get("run_label") == local_case_id
                ):
                    solver_diagnostics["run_label"] = case.case_id
    return _DeterministicSweepResult(
        mixture_id=job.mixture.mixture_id,
        solver_id=job.solver_id,
        cases=cases,
    )


def _iter_deterministic_results(
    jobs: list[_DeterministicSweepJob],
    *,
    workers: int,
) -> Iterator[_DeterministicSweepResult]:
    if workers <= 1:
        for job in jobs:
            yield _execute_deterministic_job(job)
        return
    if not jobs:
        return

    context = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=min(int(workers), len(jobs)),
        mp_context=context,
    ) as executor:
        futures = [executor.submit(_execute_deterministic_job, job) for job in jobs]
        try:
            # Read futures in submission order. Execution remains parallel,
            # while database checkpoints and exported rows stay deterministic.
            for future in futures:
                yield future.result()
        except BaseException:
            for future in futures:
                future.cancel()
            raise


def _execute_monte_carlo_job(job: _MonteCarloSweepJob) -> _MonteCarloSweepResult:
    """Run one anchor/replica without opening the workflow database."""

    base_config = load_config(job.base_config_path)
    config = _clone_config(
        base_config,
        mixture=job.mixture,
        e_over_n_values=(job.e_over_n_Td,),
        solver_ids=["monte_carlo"],
        replicate=job.replicate,
        seed=job.seed,
    )
    config.solvers.monte_carlo.particles = job.sampling.particles
    config.solvers.monte_carlo.warmup_collisions = job.sampling.warmup_collisions
    config.solvers.monte_carlo.max_collisions = job.sampling.max_collisions
    config.solvers.monte_carlo.tail_max_collisions = (
        job.sampling.tail_max_collisions
        if job.sampling.tail_max_collisions > 0
        else None
    )
    config.solvers.monte_carlo.transport_correlation_lag_barriers = (
        job.sampling.transport_correlation_lag_barriers
    )
    config.solvers.monte_carlo.transport_estimator = job.sampling.transport_estimator
    return _MonteCarloSweepResult(
        mixture_id=job.mixture.mixture_id,
        replicate=job.replicate,
        cases=tuple(run(config, write=False).cases),
    )


def _iter_monte_carlo_results(
    jobs: list[_MonteCarloSweepJob],
    *,
    workers: int,
) -> Iterator[_MonteCarloSweepResult]:
    if workers <= 1:
        for job in jobs:
            yield _execute_monte_carlo_job(job)
        return
    if not jobs:
        return

    # Explicit spawn avoids inheriting an open SQLite connection on POSIX and
    # uses the only safe process-start model on Windows. Workers return plain
    # result objects; all WorkflowStore writes remain in the parent process.
    context = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=min(int(workers), len(jobs)),
        mp_context=context,
    ) as executor:
        futures = [executor.submit(_execute_monte_carlo_job, job) for job in jobs]
        for future in as_completed(futures):
            yield future.result()


def _mc_reuse_keys(
    source: sqlite3.Connection,
    *,
    workflow: _workflow_config.WorkflowConfig,
    target_provenance: dict[str, str],
    tail_rate_rse_trigger: float,
) -> list[tuple[int, str, float, int]]:
    """Select source MC cases representing exactly the requested raw jobs."""

    source_metadata = read_metadata(source)
    source_tail_schema = source_metadata.get(MC_TAIL_ESTIMATOR_METADATA_KEY)
    target_tail_schema = target_provenance.get(MC_TAIL_ESTIMATOR_METADATA_KEY)
    tail_schema_metadata_upgrade = bool(
        source_tail_schema is None
        and target_tail_schema == _mc_evidence.WEIGHTED_MC_TAIL_ESTIMATOR_SCHEMA_VERSION
    )
    if source_tail_schema != target_tail_schema and not tail_schema_metadata_upgrade:
        raise ValueError(
            f"workflow mc.reuse_database has different {MC_TAIL_ESTIMATOR_METADATA_KEY}"
        )
    for name in (
        "base_config_sha256",
        "cross_sections_sha256",
        MC_TRANSPORT_ESTIMATOR_METADATA_KEY,
        MC_EEDF_ESTIMATOR_METADATA_KEY,
        MC_SEED_DERIVATION_METADATA_KEY,
        MC_SOLVER_SOURCE_METADATA_KEY,
    ):
        if source_metadata.get(name) != target_provenance.get(name):
            raise ValueError(f"workflow mc.reuse_database has different {name}")

    source_plan_provenance = mc_sampling_plan_provenance(
        source_metadata,
        required=True,
    )
    assert source_plan_provenance is not None
    source_plan = {
        float(row["e_over_n_Td"]): row for row in source_plan_provenance["entries"]
    }
    target_plan = {row.e_over_n_Td: row for row in workflow.mc_sampling_plan}
    reusable_controls = (
        "particles",
        "warmup_collisions",
        "max_collisions",
        "tail_max_collisions",
        "transport_correlation_lag_barriers",
        "transport_estimator",
    )
    keys: list[tuple[int, str, float, int]] = []
    for mixture in workflow.mixtures:
        expected_fractions = json.dumps(
            dict(sorted(mixture.fractions.items())),
            sort_keys=True,
        )
        source_mixture = source.execute(
            "SELECT fractions_json FROM mixtures WHERE mixture_id = ?",
            (mixture.mixture_id,),
        ).fetchone()
        if source_mixture is None or source_mixture[0] != expected_fractions:
            continue

        for anchor, target in target_plan.items():
            source_sampling = source_plan.get(anchor)
            if source_sampling is None or any(
                source_sampling[name] != getattr(target, name)
                for name in reusable_controls
            ):
                continue
            available_replicas = min(
                int(source_sampling["replicas"]),
                target.replicas,
            )
            for replicate in range(available_replicas):
                row = source.execute(
                    """
                    SELECT diagnostics_json FROM cases
                    WHERE mixture_id = ? AND solver = 'monte_carlo'
                      AND e_over_n_Td = ? AND replicate = ?
                    """,
                    (mixture.mixture_id, anchor, replicate),
                ).fetchone()
                if row is None:
                    continue
                diagnostics = json.loads(row[0])
                transport_diagnostics = diagnostics.get(
                    "internal_monte_carlo_transport", {}
                )
                run_provenance = transport_diagnostics.get("mc_run_provenance", {})
                expected_seed = stable_mc_seed(
                    base_seed=workflow.mc_base_seed,
                    mixture_id=mixture.mixture_id,
                    e_over_n_Td=anchor,
                    replicate=replicate,
                    solver="monte_carlo",
                )
                expected_run = {
                    "seed": expected_seed,
                    "particles": target.particles,
                    "warmup_collisions": target.warmup_collisions,
                    "production_collisions": target.max_collisions,
                    "tail_max_collisions": target.tail_max_collisions,
                    "transport_estimator": target.transport_estimator,
                    "solver_source_sha256": target_provenance[
                        MC_SOLVER_SOURCE_METADATA_KEY
                    ],
                }
                if any(
                    run_provenance.get(name) != value
                    for name, value in expected_run.items()
                ):
                    continue
                if (
                    run_provenance.get("eedf_estimator_schema_version")
                    != (target_provenance[MC_EEDF_ESTIMATOR_METADATA_KEY])
                ):
                    continue
                if target_tail_schema is not None and not _has_current_tail_evidence(
                    transport_diagnostics,
                    configured_collisions=target.tail_max_collisions,
                    rate_rse_trigger=tail_rate_rse_trigger,
                    estimator_schema=target_tail_schema,
                ):
                    continue
                if target.transport_estimator == "paired_field_parity" and (
                    transport_diagnostics.get("field_parity_response", {}).get(
                        "estimator_schema_version"
                    )
                    != _mc_evidence.FIELD_PARITY_RESPONSE_ESTIMATOR_SCHEMA_VERSION
                ):
                    continue
                keys.append((mixture.mixture_id, "monte_carlo", anchor, replicate))
    return keys


def _has_current_tail_evidence(
    transport_diagnostics: object,
    *,
    configured_collisions: int,
    rate_rse_trigger: float,
    estimator_schema: str,
) -> bool:
    """Return whether one raw case carries the complete current tail contract."""

    if not isinstance(transport_diagnostics, dict):
        return False
    run = transport_diagnostics.get("mc_run_provenance")
    tail = transport_diagnostics.get("tail_sampling")
    if not isinstance(run, dict) or not isinstance(tail, dict):
        return False

    configured = tail.get("configured_max_collisions")
    executed = tail.get("executed_collisions")
    if (
        isinstance(configured, bool)
        or not isinstance(configured, int)
        or configured != configured_collisions
        or isinstance(executed, bool)
        or not isinstance(executed, int)
        or executed not in {0, configured_collisions}
    ):
        return False
    expected_case_schema = estimator_schema if configured_collisions > 0 else None
    if (
        tail.get("estimator_schema_version") != expected_case_schema
        or run.get("tail_max_collisions") != configured_collisions
        or run.get("tail_collisions_executed") != executed
    ):
        return False

    activation_trigger = tail.get("activation_rate_rse_trigger")
    run_trigger = run.get("tail_rate_rse_trigger")
    if (
        isinstance(activation_trigger, bool)
        or not isinstance(activation_trigger, (int, float))
        or not math.isfinite(float(activation_trigger))
        or float(activation_trigger) != float(rate_rse_trigger)
        or isinstance(run_trigger, bool)
        or not isinstance(run_trigger, (int, float))
        or not math.isfinite(float(run_trigger))
        or float(run_trigger) != float(rate_rse_trigger)
    ):
        return False

    tail_executed = executed > 0
    expected_model = (
        "reaction_kernel_weighted_ensemble"
        if tail_executed
        else (
            "ordinary_trajectory_sampling_resolved"
            if configured_collisions > 0
            else "ordinary_trajectory_sampling"
        )
    )
    strata = tail.get("strata_edges_eV")
    strata_are_current = (
        isinstance(strata, list) and len(strata) >= 2
        if configured_collisions > 0
        else strata == []
    )
    return bool(
        tail.get("model") == expected_model
        and strata_are_current
        and tail.get("reported_eedf_includes_main_production") is True
        and tail.get("reported_rates_include_main_production") is True
        and tail.get("transport_ensemble_resampled_for_tail") is False
    )


def _reuse_monte_carlo_cases(
    store: WorkflowStore,
    *,
    workflow: _workflow_config.WorkflowConfig,
    target_provenance: dict[str, str],
    tail_rate_rse_trigger: float,
) -> int:
    source_path = workflow.mc_reuse_database
    if source_path is None:
        return 0
    if not source_path.is_file():
        raise FileNotFoundError(
            f"workflow mc.reuse_database does not exist: {source_path}"
        )
    if source_path.resolve() == store.path.resolve():
        return 0

    source_uri = source_path.resolve().as_uri() + "?mode=ro"
    with sqlite3.connect(source_uri, uri=True) as source:
        validate_workflow_schema(source)
        keys = _mc_reuse_keys(
            source,
            workflow=workflow,
            target_provenance=target_provenance,
            tail_rate_rse_trigger=tail_rate_rse_trigger,
        )
        missing = [key for key in keys if key not in store.existing_case_keys()]
        imported = store.import_cases(source, missing)
        return imported


def advance_monte_carlo_workflow(
    workflow_path: str | Path,
    decision_path: str | Path,
    *,
    output_path: str | Path,
    database_path: str | Path,
) -> Path:
    """Materialize exactly the one follow-up authorized by a pending decision."""
    workflow = _workflow_config.load_workflow(workflow_path)
    decision_file = Path(decision_path).resolve()
    decision = read_selection(decision_file)
    policy = parse_convergence_policy(decision.get("policy"))
    entries = decision.get("next_sampling_plan")
    if not isinstance(entries, list):
        raise MonteCarloPolicyError("this decision has no next MC sampling plan")
    plan = tuple(SamplingPlanEntry(**row) for row in entries)
    validate_sampling_budget(plan, policy, previous_decision=decision_file)
    if [asdict(entry) for entry in workflow.mc_sampling_plan] != decision.get(
        "current_sampling_plan"
    ):
        raise MonteCarloPolicyError(
            "source workflow does not match the completed MC plan"
        )
    output = Path(output_path).resolve()
    database = Path(database_path).resolve()
    if (
        output == workflow.path
        or database == workflow.database_path
        or database.exists()
    ):
        raise MonteCarloPolicyError(
            "the next attempt requires separate workflow and database paths"
        )
    raw, _ = _workflow_config._load_workflow_raw(workflow.path)
    raw["base_config"] = str(workflow.base_config_path)
    raw["database"] = str(database)
    raw["mc"]["convergence"] = asdict(policy)
    raw["mc"]["previous_decision"] = str(decision_file)
    raw["mc"]["previous_decision_sha256"] = file_sha256(decision_file)
    raw["mc"]["sampling_plan"] = entries
    if decision["reuse_previous_mc_database"]:
        raw["mc"]["reuse_database"] = str(workflow.database_path)
    else:
        raw["mc"].pop("reuse_database", None)
    encoded = yaml.safe_dump(raw, sort_keys=False, allow_unicode=True)
    if output.exists() and output.read_text(encoding="utf-8") != encoded:
        raise MonteCarloPolicyError(
            "next workflow already exists with different settings"
        )
    # Validate through the same loader before publishing the user-visible YAML.
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".mc-plan-", dir=output.parent) as staging:
        staged = Path(staging) / "workflow.yaml"
        staged.write_text(encoded, encoding="utf-8")
        _workflow_config.load_workflow(staged)
    output.write_text(encoded, encoding="utf-8")
    return output


def run_sweep(path: str | Path) -> SweepSummary:
    workflow = _workflow_config.load_workflow(path)
    base_config = load_config(workflow.base_config_path)
    validate_unique_solver_ids(base_config)
    solver_ids = _enabled_solver_ids(base_config)
    if not solver_ids:
        raise ValueError("base_config run.solvers has no enabled solvers")
    deterministic_solvers = [
        solver_id for solver_id in solver_ids if solver_id != "monte_carlo"
    ]
    deterministic_execution = _resolve_deterministic_execution(
        workflow,
        base_config,
        solver_ids,
    )
    deterministic_execution_payload = (
        _deterministic_execution_payload(deterministic_execution)
        if workflow.deterministic_execution.configured
        else None
    )

    xs_hashes = _workflow_config._cross_section_file_hashes(base_config)
    mc_estimator_version = (
        (
            _mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
            if base_config.solvers.monte_carlo.population_model == "weighted_branching"
            else _mc_evidence.DIRECT_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
        )
        if "monte_carlo" in solver_ids
        else None
    )
    mc_eedf_estimator_version = (
        (
            _mc_evidence.MC_MAGNETIC_EEDF_ESTIMATOR_SCHEMA_VERSION
            if (
                base_config.physics.field.magnetic_field.enabled
                and base_config.physics.field.magnetic_field.B_T > 0.0
            )
            else _mc_evidence.MC_EEDF_ESTIMATOR_SCHEMA_VERSION
        )
        if "monte_carlo" in solver_ids
        else None
    )
    mc_tail_estimator_version = (
        _mc_evidence.WEIGHTED_MC_TAIL_ESTIMATOR_SCHEMA_VERSION
        if "monte_carlo" in solver_ids
        and base_config.solvers.monte_carlo.population_model == "weighted_branching"
        and any(entry.tail_max_collisions > 0 for entry in workflow.mc_sampling_plan)
        else None
    )
    mc_seed_derivation_version = (
        _mc_evidence.MC_SEED_DERIVATION_SCHEMA_VERSION
        if "monte_carlo" in solver_ids
        else None
    )
    mc_solver_source = (
        monte_carlo_source_sha256() if "monte_carlo" in solver_ids else None
    )
    propagator_solver_source = (
        str(
            propagator_qualification_source_fingerprint(_REPOSITORY_ROOT)["sha256"]
        )
        if "propagator" in solver_ids
        else None
    )
    provenance = {
        "workflow_config_sha256": workflow_config_sha256(
            workflow,
            mc_transport_estimator_schema_version=mc_estimator_version,
            mc_eedf_estimator_schema_version=mc_eedf_estimator_version,
            mc_tail_estimator_schema_version=mc_tail_estimator_version,
            mc_seed_derivation_schema_version=mc_seed_derivation_version,
            mc_solver_source_sha256=mc_solver_source,
            propagator_solver_source_sha256=propagator_solver_source,
            deterministic_execution_payload=deterministic_execution_payload,
        ),
        QUALITY_THRESHOLDS_METADATA_KEY: quality_thresholds_json(workflow.quality),
        "base_config_path": str(workflow.base_config_path),
        "base_config_sha256": _workflow_config._sha256_file(workflow.base_config_path),
        "cross_sections_sha256": _workflow_config._combined_hash(xs_hashes),
        "cross_section_files_json": json.dumps(xs_hashes, sort_keys=True),
        PHYSICAL_CONTEXT_KEY: json.dumps(physical_context(base_config), sort_keys=True),
    }
    if mc_estimator_version is not None:
        provenance[MC_TRANSPORT_ESTIMATOR_METADATA_KEY] = mc_estimator_version
        provenance[MC_SAMPLING_PLAN_METADATA_KEY] = mc_sampling_plan_json(workflow)
        provenance["mc_nominal_compute_budget_json"] = json.dumps(
            sampling_budget_provenance(workflow.mc_sampling_plan),
            sort_keys=True,
        )
        provenance["mc_base_seed"] = str(workflow.mc_base_seed)
    if mc_eedf_estimator_version is not None:
        provenance[MC_EEDF_ESTIMATOR_METADATA_KEY] = mc_eedf_estimator_version
    if mc_seed_derivation_version is not None:
        provenance[MC_SEED_DERIVATION_METADATA_KEY] = mc_seed_derivation_version
    if mc_solver_source is not None:
        provenance[MC_SOLVER_SOURCE_METADATA_KEY] = mc_solver_source
    if propagator_solver_source is not None:
        provenance[PROPAGATOR_SOLVER_SOURCE_METADATA_KEY] = (
            propagator_solver_source
        )
    if workflow.mc_convergence is not None:
        provenance["mc_campaign_json"] = json.dumps(
            campaign_provenance(workflow.mc_convergence, workflow.mc_previous_decision),
            sort_keys=True,
        )
    if mc_tail_estimator_version is not None:
        provenance[MC_TAIL_ESTIMATOR_METADATA_KEY] = mc_tail_estimator_version
    if deterministic_execution_payload is not None:
        provenance[DETERMINISTIC_EXECUTION_METADATA_KEY] = json.dumps(
            deterministic_execution_payload,
            sort_keys=True,
            separators=(",", ":"),
        )

    cases_written = 0
    deterministic_groups: list[_DeterministicSweepGroup] = []
    mc_jobs: list[_MonteCarloSweepJob] = []
    with WorkflowStore(workflow.database_path) as store:
        validation_provenance, quality_policy_reevaluated = resume_provenance(
            store,
            workflow,
            provenance,
            deterministic_execution_payload=deterministic_execution_payload,
        )
        store.set_provenance(validation_provenance)
        existing_case_keys = store.existing_case_keys()
        for mixture in workflow.mixtures:
            store.write_mixture(
                mixture.mixture_id,
                _mixture_species_rows(base_config, mixture),
            )

        cases_written += _reuse_monte_carlo_cases(
            store,
            workflow=workflow,
            target_provenance=provenance,
            tail_rate_rse_trigger=(
                base_config.solvers.monte_carlo.tail_rate_rse_trigger
            ),
        )
        existing_case_keys = store.existing_case_keys()

        for mixture in workflow.mixtures:
            for solver_id in deterministic_solvers:
                indexed_missing_e_over_n = tuple(
                    (index, float(e_over_n))
                    for index, e_over_n in enumerate(workflow.e_over_n_Td)
                    if (
                        mixture.mixture_id,
                        solver_id,
                        float(e_over_n),
                        0,
                    )
                    not in existing_case_keys
                )
                if not indexed_missing_e_over_n:
                    continue
                deterministic_groups.append(
                    _DeterministicSweepGroup(
                        base_config_path=workflow.base_config_path,
                        mixture=mixture,
                        solver_id=solver_id,
                        indexed_e_over_n_Td=indexed_missing_e_over_n,
                        propagator_source_sha256=(
                            propagator_solver_source
                            if solver_id == "propagator"
                            else None
                        ),
                    )
                )

            if "monte_carlo" in solver_ids:
                # Submit one replica across all anchors before the next replica.
                # This preserves every configured job/seed while avoiding a
                # parallel batch made entirely from one potentially slow E/N.
                maximum_replicas = max(
                    sampling.replicas for sampling in workflow.mc_sampling_plan
                )
                for replicate in range(maximum_replicas):
                    for sampling in workflow.mc_sampling_plan:
                        if replicate >= sampling.replicas:
                            continue
                        e_over_n = sampling.e_over_n_Td
                        if (
                            mixture.mixture_id,
                            "monte_carlo",
                            float(e_over_n),
                            replicate,
                        ) in existing_case_keys:
                            continue
                        mc_jobs.append(
                            _MonteCarloSweepJob(
                                base_config_path=workflow.base_config_path,
                                mixture=mixture,
                                e_over_n_Td=float(e_over_n),
                                replicate=replicate,
                                seed=stable_mc_seed(
                                    base_seed=workflow.mc_base_seed,
                                    mixture_id=mixture.mixture_id,
                                    e_over_n_Td=e_over_n,
                                    replicate=replicate,
                                    solver="monte_carlo",
                                ),
                                sampling=sampling,
                            )
                        )

        deterministic_jobs = _build_deterministic_jobs(
            deterministic_groups,
            workers=deterministic_execution.effective_worker_cap,
        )
        for deterministic_result in _iter_deterministic_results(
            deterministic_jobs,
            workers=deterministic_execution.effective_worker_cap,
        ):
            for case in deterministic_result.cases:
                store.write_case(
                    mixture_id=deterministic_result.mixture_id,
                    replicate=0,
                    case=case,
                )
                existing_case_keys.add(
                    (
                        deterministic_result.mixture_id,
                        case.solver,
                        float(case.e_over_n_Td),
                        0,
                    )
                )
                cases_written += 1

        if workflow.mean_energy_support is not None:
            _require_current_propagator_source(propagator_solver_source)
            extension_cases = _extend_workflow_mean_energy_support(
                store,
                workflow=workflow,
                base_config=base_config,
                solver_id=deterministic_solvers[0],
            )
            cases_written += len(extension_cases)
            _require_current_propagator_source(propagator_solver_source)

        for mc_result in _iter_monte_carlo_results(
            mc_jobs,
            workers=workflow.mc_workers,
        ):
            for case in mc_result.cases:
                store.write_case(
                    mixture_id=mc_result.mixture_id,
                    replicate=mc_result.replicate,
                    case=case,
                )
                existing_case_keys.add(
                    (
                        mc_result.mixture_id,
                        case.solver,
                        float(case.e_over_n_Td),
                        mc_result.replicate,
                    )
                )
                cases_written += 1
        _require_current_propagator_source(propagator_solver_source)
        aggregate_workflow_results(
            store.connection,
            workflow.quality,
            _allow_evaluation_policy_override=quality_policy_reevaluated,
        )

    return SweepSummary(
        database_path=workflow.database_path,
        mixtures=len(workflow.mixtures),
        e_over_n_values=len(workflow.e_over_n_Td),
        cases_written=cases_written,
        quality_policy_reevaluated=quality_policy_reevaluated,
    )
