"""Command-line and Python entry points for product swarm runs."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from electron_swarm.collisions.ee_postprocess import (
    apply_electron_electron_relaxation_from_config,
)
from electron_swarm.core.config import MIGRATION_ERROR, SwarmConfig, load_config
from electron_swarm.core.cross_sections import load_cross_sections
from electron_swarm.core.result_metadata import normalize_product_metadata
from electron_swarm.core.results import SwarmRunResult
from electron_swarm.core.solver_configs import build_internal_solver_configs
from electron_swarm.diagnostics import enrich_tail_metrics
from electron_swarm.io.writers import write_outputs
from electron_swarm.orchestration.comparison import build_comparison_summary
from electron_swarm.orchestration.executor import execute_solve_plan
from electron_swarm.orchestration.plan import (
    build_solve_plan,
    solver_plan_metadata,
)


def run(
    config: SwarmConfig,
    *,
    write: bool = True,
    collect_diagnostics: bool = False,
) -> SwarmRunResult:
    if config.schema_version != 2:
        raise ValueError(MIGRATION_ERROR)
    cross_sections = load_cross_sections(config.cross_sections, config.conditions)
    internal = build_internal_solver_configs(config.solvers, config.physics)
    internal.monte_carlo.collect_audit = collect_diagnostics
    plan = build_solve_plan(config)
    cases = execute_solve_plan(config, cross_sections, plan, internal)
    cases = apply_electron_electron_relaxation_from_config(
        cases,
        config,
        cross_sections,
    )
    cases = enrich_tail_metrics(cases, config, cross_sections)
    normalize_product_metadata(cases)
    result = SwarmRunResult(
        cases=cases,
        metadata={
            "source_config": str(config.source_path) if config.source_path else None,
            "schema_version": config.schema_version,
            "solver_plan": solver_plan_metadata(plan),
        },
    )
    build_comparison_summary(result, config.comparison)
    if write:
        paths = write_outputs(result, config.output, comparison=config.comparison)
        result.metadata["output_paths"] = {key: str(value) for key, value in paths.items()}
    return result


def run_from_config(
    path: str | Path,
    *,
    write: bool = True,
    collect_diagnostics: bool = False,
) -> SwarmRunResult:
    return run(load_config(path), write=write, collect_diagnostics=collect_diagnostics)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Product electron swarm solver comparison runner"
    )
    parser.add_argument("config", type=Path, help="schema v2 YAML configuration file")
    parser.add_argument(
        "--no-write", action="store_true", help="Run without writing CSV outputs"
    )
    parser.add_argument(
        "--log-level",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        help="Enable Python logging at the selected level",
    )
    args = parser.parse_args(argv)
    if args.log_level:
        logging.basicConfig(level=getattr(logging, args.log_level))
    result = run_from_config(args.config, write=not args.no_write)
    print(f"completed {len(result.cases)} solver cases")
    if result.metadata.get("output_paths"):
        for key, value in result.metadata["output_paths"].items():
            print(f"{key}: {value}")
