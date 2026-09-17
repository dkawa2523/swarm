"""Run the repository's canonical raw-solver EEDF comparison."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from swarm_workflow.plots.eedf_comparison import (  # noqa: E402
    generate_eedf_solver_comparison,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--two-term-directory", required=True)
    parser.add_argument("--propagator-directory", required=True)
    parser.add_argument("--mc-database", required=True)
    parser.add_argument("--mc-qualification-csv", required=True)
    parser.add_argument("--output-directory", required=True)
    parser.add_argument("--mixture-id", type=int, default=0)
    parser.add_argument(
        "--qualification-column", default="active_closure_quality_passed"
    )
    parser.add_argument(
        "--representative-fields-td", type=float, nargs="+", required=True
    )
    parser.add_argument("--tail-thresholds-eV", type=float, nargs="+", required=True)
    args = parser.parse_args()

    summary = generate_eedf_solver_comparison(
        two_term_directory=args.two_term_directory,
        propagator_directory=args.propagator_directory,
        mc_database=args.mc_database,
        mc_qualification_csv=args.mc_qualification_csv,
        output_directory=args.output_directory,
        mixture_id=args.mixture_id,
        qualification_column=args.qualification_column,
        representative_fields_td=args.representative_fields_td,
        tail_thresholds_eV=args.tail_thresholds_eV,
    )
    print(summary.manifest)


if __name__ == "__main__":
    main()
