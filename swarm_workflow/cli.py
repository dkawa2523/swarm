"""Command line entry point for external swarm workflows."""

from __future__ import annotations

import argparse
from pathlib import Path

from .sweep import run_sweep


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="swarm-workflow")
    subparsers = parser.add_subparsers(dest="command", required=True)
    sweep_parser = subparsers.add_parser(
        "sweep",
        help="Run an external E/N and mixture sweep into SQLite",
    )
    sweep_parser.add_argument("workflow", type=Path, help="workflow YAML file")

    args = parser.parse_args(argv)
    if args.command == "sweep":
        summary = run_sweep(args.workflow)
        print(
            "completed "
            f"{summary.cases_written} cases "
            f"({summary.mixtures} mixtures, {summary.e_over_n_values} E/N values)"
        )
        return
    parser.error(f"unsupported command: {args.command}")


if __name__ == "__main__":
    main()
