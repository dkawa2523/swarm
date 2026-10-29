"""Command-line and Python entry points for product swarm runs."""

from __future__ import annotations

from pathlib import Path
from electron_swarm import __version__

from electron_swarm.collisions.ee_postprocess import (
    apply_electron_electron_relaxation_from_config,
)
from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.config_parser import load_config
from electron_swarm.core.result_metadata import normalize_product_metadata
from electron_swarm.core.results import SwarmRunResult, case_numerical_status
from electron_swarm.diagnostics import enrich_tail_metrics
from electron_swarm.io.writers import write_outputs
from electron_swarm.orchestration.comparison import build_comparison_summary
from electron_swarm.orchestration.preparation import prepare_run
from electron_swarm.core.results import execution_outcome
from electron_swarm.core.provenance import seal_execution
from electron_swarm.orchestration.executor import execute_solve_plan
from electron_swarm.orchestration.plan import (
    solver_plan_metadata,
)


def run(
    config: SwarmConfig,
    *,
    write: bool = True,
    collect_diagnostics: bool = False,
) -> SwarmRunResult:
    prepared = prepare_run(config, collect_audit=collect_diagnostics)
    config = prepared.config
    active_inputs = prepared.inputs
    plan, internal = prepared.plan, prepared.internal
    cases = execute_solve_plan(config, active_inputs, plan, internal)
    contexts = {case.e_over_n_Td: case for case in prepared.physical_cases}
    for case in cases:
        case.physical_context = contexts[case.e_over_n_Td]
        case.transport_observables = case.transport.observable_definitions
    cases = apply_electron_electron_relaxation_from_config(
        cases,
        config,
        active_inputs,
    )
    cases = enrich_tail_metrics(cases, config, active_inputs)
    for case in cases:
        case.metadata["numerical_status"] = case_numerical_status(case)
    normalize_product_metadata(cases)
    prepared.verify_unchanged()
    result = SwarmRunResult(
        cases=cases,
        metadata={
            "package_version": __version__,
            "source_config": str(config.source_path) if config.source_path else None,
            "schema_version": config.schema_version,
            "solver_plan": solver_plan_metadata(plan),
            "configuration": prepared.configuration,
            "implementations": prepared.implementations,
            "input_files": prepared.input_files,
            "resolved_solvers": prepared.resolved_solvers,
            "input_evidence": prepared.input_evidence,
        },
    )
    build_comparison_summary(result, config.comparison)
    result.metadata["outcome"] = execution_outcome(result)
    seal_execution(result)
    if write:
        paths = write_outputs(result, config.output)
        result.metadata["output_paths"] = {
            key: str(value) for key, value in paths.items()
        }
    return result


def run_from_config(
    path: str | Path,
    *,
    write: bool = True,
    collect_diagnostics: bool = False,
) -> SwarmRunResult:
    return run(load_config(path), write=write, collect_diagnostics=collect_diagnostics)
