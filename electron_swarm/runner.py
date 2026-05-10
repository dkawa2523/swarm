"""Command-line and Python entry points for unified swarm runs."""

from __future__ import annotations

import argparse
from pathlib import Path

from electron_swarm.collisions import augment_cross_sections_from_config
from electron_swarm.collisions.postprocess import apply_case_hooks
from electron_swarm.core.config import SwarmConfig, load_config
from electron_swarm.core.cross_sections import load_cross_sections
from electron_swarm.core.results import SwarmRunResult
from electron_swarm.diagnostics import enrich_run_diagnostics
from electron_swarm.io.writers import write_outputs
from electron_swarm.plotting import write_plots
from electron_swarm.solvers.boltzmann_two_term import BoltzmannTwoTermSolver
from electron_swarm.solvers.monte_carlo_adapter import MonteCarloAdapter
from electron_swarm.solvers.multiterm_boltzmann import MultiTermBoltzmannSolver


def run(config: SwarmConfig, *, write: bool = True) -> SwarmRunResult:
    cross_sections = load_cross_sections(config.cross_sections, config.conditions)
    cross_sections = augment_cross_sections_from_config(cross_sections, config)
    cases = []
    if (
        config.run.mode in {"boltzmann_two_term", "both", "all"}
        and config.boltzmann_two_term.enabled
    ):
        cases.extend(BoltzmannTwoTermSolver(config, cross_sections).solve_all())
    if (
        config.run.mode in {"multiterm_boltzmann", "all"}
        and config.multiterm_boltzmann.enabled
    ):
        cases.extend(MultiTermBoltzmannSolver(config, cross_sections).solve_all())
    if (
        config.run.mode in {"monte_carlo", "both", "all"}
        and config.monte_carlo.enabled
    ):
        cases.extend(MonteCarloAdapter(config, cross_sections).solve_all())
    cases = apply_case_hooks(cases, config, cross_sections)
    result = SwarmRunResult(
        cases=cases,
        metadata={
            "source_config": str(config.source_path) if config.source_path else None
        },
    )
    enrich_run_diagnostics(result)
    if write:
        paths = write_outputs(result, config.output)
        paths.update(write_plots(result, config.output))
        result.metadata["output_paths"] = {k: str(v) for k, v in paths.items()}
    return result


def run_from_config(path: str | Path, *, write: bool = True) -> SwarmRunResult:
    return run(load_config(path), write=write)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Unified particle-MC/Boltzmann two-term electron swarm runner"
    )
    parser.add_argument("config", type=Path, help="YAML configuration file")
    parser.add_argument(
        "--no-write", action="store_true", help="Run without writing CSV or plot outputs"
    )
    args = parser.parse_args(argv)
    result = run_from_config(args.config, write=not args.no_write)
    print(f"completed {len(result.cases)} solver cases")
    if result.metadata.get("output_paths"):
        for key, value in result.metadata["output_paths"].items():
            print(f"{key}: {value}")
