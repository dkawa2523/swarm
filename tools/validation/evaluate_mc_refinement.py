"""Evaluate raw Monte Carlo population refinement without solver substitution.

This validation-only entry point compares completed Monte Carlo SQLite
campaigns at 512, 1024, and 2048 particles. It never replaces an EEDF with a
two-term result, smooths a histogram, or treats qualification as EEDF distance.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Iterable, Sequence


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from swarm_workflow.plots.eedf_contracts import EedfComparisonError  # noqa: E402
from tools.validation.mc_refinement_contracts import (  # noqa: E402
    DEFAULT_MIXTURES,
    DEFAULT_POPULATIONS,
    DatabaseInput,
    RefinementEvaluationError,
)
from tools.validation.mc_refinement_database import (  # noqa: E402
    project_machine_resolution_cells as _project_machine_resolution_cells,
)
from tools.validation.mc_refinement_report import (  # noqa: E402
    evaluate_mc_refinement,
)


def _parse_database_argument(value: str) -> DatabaseInput:
    parts = value.split("=", 2)
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("database must use MIXTURE=POPULATION=PATH")
    mixture, population_text, path_text = parts
    try:
        population = int(population_text)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("database population must be an integer") from exc
    if not mixture.strip() or population <= 0 or not path_text.strip():
        raise argparse.ArgumentTypeError("database argument contains an empty value")
    return DatabaseInput(mixture.strip(), population, Path(path_text))


def _parse_qualification_argument(value: str) -> tuple[str, Path]:
    parts = value.split("=", 1)
    if len(parts) != 2 or not parts[0].strip() or not parts[1].strip():
        raise argparse.ArgumentTypeError("qualification must use MIXTURE=PATH")
    return parts[0].strip(), Path(parts[1])


def _inputs_from_root(
    root: Path,
    mixtures: Iterable[str],
    populations: Iterable[int],
) -> list[DatabaseInput]:
    return [
        DatabaseInput(
            mixture,
            int(population),
            root / mixture / f"p{int(population)}" / "monte_carlo.sqlite",
        )
        for mixture in mixtures
        for population in populations
    ]


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--root",
        type=Path,
        help="root containing MIXTURE/pPOPULATION/monte_carlo.sqlite",
    )
    source.add_argument(
        "--database",
        action="append",
        type=_parse_database_argument,
        help="explicit input as MIXTURE=POPULATION=PATH; repeat for every cell",
    )
    parser.add_argument("--mixtures", nargs="+", default=list(DEFAULT_MIXTURES))
    parser.add_argument(
        "--populations", nargs="+", type=int, default=list(DEFAULT_POPULATIONS)
    )
    parser.add_argument(
        "--qualification",
        action="append",
        type=_parse_qualification_argument,
        default=[],
        help="optional qualification evidence as MIXTURE=PATH",
    )
    parser.add_argument("--mixture-id", type=int, default=0)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    if args.root is not None:
        inputs = _inputs_from_root(args.root, args.mixtures, args.populations)
    else:
        inputs = list(args.database)
    qualification_paths = dict(args.qualification)
    if len(qualification_paths) != len(args.qualification):
        parser.error("duplicate --qualification mixture")
    try:
        paths = evaluate_mc_refinement(
            inputs,
            output_directory=args.output_directory,
            qualification_paths=qualification_paths,
            mixture_id=args.mixture_id,
            populations=args.populations,
            overwrite=args.force,
        )
    except (RefinementEvaluationError, EedfComparisonError) as exc:
        parser.error(str(exc))
    print(paths["manifest"])
    return 0


__all__ = [
    "DatabaseInput",
    "RefinementEvaluationError",
    "_project_machine_resolution_cells",
    "evaluate_mc_refinement",
    "main",
]


if __name__ == "__main__":
    raise SystemExit(main())
