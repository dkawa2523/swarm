"""Command-line adapter; numerical execution stays in runner."""

import argparse
import json
import logging
from pathlib import Path

from electron_swarm import load_config
from electron_swarm.runner import run_from_config
from electron_swarm.orchestration.preparation import prepare_run
from electron_swarm.core.yaml_input import read_yaml_mapping


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Product electron swarm solver comparison runner"
    )
    parser.add_argument("config", type=Path, help="schema v2 YAML configuration file")
    parser.add_argument(
        "--no-write", action="store_true", help="Run without writing CSV outputs"
    )
    parser.add_argument(
        "--plan",
        action="store_true",
        help="Validate inputs and print the solver plan without solving",
    )
    parser.add_argument(
        "--log-level",
        choices=("DEBUG", "INFO", "WARNING", "ERROR"),
        help="Enable Python logging at the selected level",
    )
    args = parser.parse_args(argv)
    if "campaign" in read_yaml_mapping(args.config):
        parser.error("campaign settings require: swarm-workflow run <config.yaml>")
    if args.log_level:
        logging.basicConfig(level=getattr(logging, args.log_level))
    if args.plan:
        prepared = prepare_run(load_config(args.config))
        print(json.dumps(prepared.description(), indent=2))
        return 0
    result = run_from_config(args.config, write=not args.no_write)
    outcome = result.metadata["outcome"]
    print(
        f"{outcome['status']}: {len(result.cases)} solver cases; comparison={outcome['comparison_status']}"
    )
    if result.metadata.get("output_paths"):
        for key, value in result.metadata["output_paths"].items():
            print(f"{key}: {value}")

    return {"complete": 0, "partial": 4, "skipped": 3}[outcome["status"]]
