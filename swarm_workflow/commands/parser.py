"""Argument definitions for the ``swarm-workflow`` command."""

from __future__ import annotations

import argparse
from pathlib import Path

from swarm_workflow.commands import handlers
from swarm_workflow.tables.contracts import (
    MC_QUALIFICATION_PROFILES,
    SOURCE_CHOICES,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="swarm-workflow")
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_campaign_commands(subparsers)
    _add_input_commands(subparsers)
    _add_comsol_model_commands(subparsers)
    return parser


def _add_campaign_commands(subparsers: argparse._SubParsersAction) -> None:
    sweep = subparsers.add_parser(
        "sweep",
        help="Run an external E/N and mixture sweep into SQLite",
    )
    sweep.add_argument("workflow", type=Path, help="workflow YAML file")
    sweep.set_defaults(handler=handlers.run_sweep_command)

    aggregate = subparsers.add_parser(
        "aggregate",
        help="Aggregate workflow SQLite results into CSV and manifest files",
    )
    aggregate.add_argument("database", type=Path, help="workflow SQLite database")
    aggregate.add_argument(
        "--output",
        required=True,
        type=Path,
        help="aggregate output directory",
    )
    aggregate.set_defaults(handler=handlers.aggregate_command)

    decide_mc = subparsers.add_parser(
        "decide-mc",
        help="Choose the next bounded MC action or one complete closure source",
    )
    decide_mc.add_argument(
        "--mc-tables",
        required=True,
        type=Path,
        help="monte_carlo mixture table directory",
    )
    decide_mc.add_argument(
        "--two-term-tables",
        required=True,
        type=Path,
        help="qualified two_term mixture table directory",
    )
    decide_mc.add_argument(
        "--attempt",
        required=True,
        type=int,
        help="completed MC attempt number",
    )
    decide_mc.add_argument(
        "--previous-decision",
        type=Path,
        help="immutable decision authorizing this follow-up attempt",
    )
    decide_mc.add_argument(
        "--output",
        required=True,
        type=Path,
        help="solver-selection JSON artifact",
    )
    decide_mc.set_defaults(handler=handlers.decide_mc_command)

    advance = subparsers.add_parser(
        "advance-mc",
        help="Write the exact bounded follow-up workflow from an MC decision",
    )
    advance.add_argument("workflow", type=Path)
    advance.add_argument("--decision", required=True, type=Path)
    advance.add_argument("--output", required=True, type=Path)
    advance.add_argument("--database", required=True, type=Path)
    advance.set_defaults(handler=handlers.advance_mc_command)


def _add_input_commands(subparsers: argparse._SubParsersAction) -> None:
    tables = subparsers.add_parser(
        "build-tables",
        help="Build per-mixture COMSOL-ready workflow table directories",
    )
    tables.add_argument("database", type=Path, help="workflow SQLite database")
    tables.add_argument(
        "--output",
        required=True,
        type=Path,
        help="table output directory",
    )
    tables.add_argument(
        "--source",
        required=True,
        choices=SOURCE_CHOICES,
        help=(
            "canonical solver source; deterministic sources are validated "
            "once per anchor, while monte_carlo keeps only independently "
            "qualified raw MC anchors and never substitutes another solver"
        ),
    )
    tables.add_argument(
        "--allow-unqualified-mc",
        action="store_true",
        help=(
            "retain all MC quality evidence when statistics cannot yet form "
            "coefficient tables"
        ),
    )
    tables.add_argument(
        "--solver-qualification",
        type=Path,
        help=(
            "immutable solver qualification JSON; required for propagator "
            "tables and copied into every downstream bundle"
        ),
    )
    tables.add_argument(
        "--target-qualification",
        type=Path,
        help=(
            "optional target-specific refinement JSON; validated against "
            "the propagator workflow database and copied downstream"
        ),
    )
    tables.add_argument(
        "--mc-qualification-profile",
        choices=MC_QUALIFICATION_PROFILES,
        default="full_transport",
        help=(
            "Monte Carlo evidence gate; restricted LMEA qualifies only the "
            "MC inputs consumed by that COMSOL closure while retaining all "
            "other MC diagnostics"
        ),
    )
    tables.set_defaults(handler=handlers.build_tables_command)

    export = subparsers.add_parser(
        "export-comsol",
        help="Export built workflow tables as COMSOL CSV bundles",
    )
    export.add_argument(
        "table_directory",
        type=Path,
        help="build-tables root directory or one mixture table directory",
    )
    export.add_argument(
        "--selection",
        type=Path,
        help="bind this bundle to a completed whole-closure solver selection",
    )
    export.add_argument(
        "--output",
        required=True,
        type=Path,
        help="output bundle directory",
    )
    export.set_defaults(handler=handlers.export_comsol_command)


def _add_comsol_model_commands(subparsers: argparse._SubParsersAction) -> None:
    positive_column = subparsers.add_parser(
        "run-positive-column",
        help="Apply, verify, run, and export the positive-column COMSOL model",
    )
    positive_column.add_argument(
        "mapping",
        type=Path,
        help="COMSOL mapping YAML file",
    )
    positive_column.add_argument(
        "--dry-run",
        action="store_true",
        help="validate and generate Java without executing COMSOL",
    )
    positive_column.add_argument(
        "--comsol",
        type=Path,
        help="COMSOL executable path used for batch execution",
    )
    positive_column.add_argument(
        "--bundle",
        required=True,
        type=Path,
        help="export-comsol bundle directory to import",
    )
    positive_column.add_argument(
        "--pressure-Pa",
        type=float,
        help="explicit positive-column gas pressure in Pa",
    )
    positive_column.add_argument(
        "--gas-temperature-K",
        type=float,
        help="explicit positive-column gas temperature in K",
    )
    positive_column.add_argument(
        "--mesh-elements",
        type=int,
        choices=(200, 400),
        help="explicit comp1/mesh1/edg1/dis1 element count",
    )
    positive_column.set_defaults(handler=handlers.run_positive_column_command)

    gec = subparsers.add_parser(
        "run-gec-ccp",
        help="Prepare or run an external-Swarm argon GEC CCP target",
    )
    gec.add_argument("mapping", type=Path, help="GEC CCP mapping YAML file")
    gec.add_argument(
        "--dry-run",
        action="store_true",
        help="inspect the MPH, validate tables, and generate Java only",
    )
    gec.add_argument(
        "--comsol",
        type=Path,
        help="COMSOL executable path used for batch execution",
    )
    gec.add_argument(
        "--bundle",
        required=True,
        type=Path,
        help="export-comsol bundle directory to import",
    )
    gec.set_defaults(handler=handlers.run_gec_ccp_command)

    icp = subparsers.add_parser(
        "run-gec-icp",
        help="Prepare or run an external-Swarm argon GEC ICP target",
    )
    icp.add_argument("mapping", type=Path, help="GEC ICP run mapping YAML file")
    icp.add_argument(
        "--dry-run",
        action="store_true",
        help="inspect the MPH, validate the bundle, and generate Java only",
    )
    icp.add_argument(
        "--comsol",
        type=Path,
        help="COMSOL executable path used for batch execution",
    )
    icp.add_argument(
        "--bundle",
        required=True,
        type=Path,
        help="one-mixture export-comsol bundle directory to import",
    )
    icp.set_defaults(handler=handlers.run_gec_icp_command)

    plot_gec = subparsers.add_parser(
        "plot-gec-ccp",
        help="Plot GEC CCP Swarm inputs and available COMSOL comparisons",
    )
    plot_gec.add_argument(
        "--bundle",
        required=True,
        type=Path,
        help="export-comsol bundle directory",
    )
    plot_gec.add_argument(
        "--comsol-results",
        type=Path,
        help="directory containing swarm_tables and any optional built-in reference",
    )
    plot_gec.add_argument(
        "--output",
        required=True,
        type=Path,
        help="plot output directory",
    )
    plot_gec.set_defaults(handler=handlers.plot_gec_ccp_command)

    compare_gec = subparsers.add_parser(
        "plot-gec-ccp-comparison",
        help="Compare accepted two_term and monte_carlo GEC CCP spatial results",
    )
    compare_gec.add_argument("--two-term-bundle", required=True, type=Path)
    compare_gec.add_argument("--two-term-results", required=True, type=Path)
    compare_gec.add_argument("--monte-carlo-bundle", required=True, type=Path)
    compare_gec.add_argument("--monte-carlo-results", required=True, type=Path)
    compare_gec.add_argument(
        "--output",
        required=True,
        type=Path,
        help="comparison output directory",
    )
    compare_gec.set_defaults(handler=handlers.compare_gec_ccp_command)

    compare_positive = subparsers.add_parser(
        "compare-positive-column",
        help="Compare external-table and built-in-Boltzmann result CSV files",
    )
    compare_positive.add_argument(
        "external_csv",
        type=Path,
        help="external Swarm-table result CSV",
    )
    compare_positive.add_argument(
        "reference_csv",
        type=Path,
        help="built-in Boltzmann reference result CSV",
    )
    compare_positive.add_argument(
        "--output",
        required=True,
        type=Path,
        help="comparison output directory",
    )
    compare_positive.add_argument(
        "--external-runtime",
        type=float,
        help="optional external solve time in seconds",
    )
    compare_positive.add_argument(
        "--reference-runtime",
        type=float,
        help="optional built-in reference solve time in seconds",
    )
    compare_positive.add_argument(
        "--plot",
        action="store_true",
        help="write a PNG comparison of the spatial profiles",
    )
    compare_positive.set_defaults(handler=handlers.compare_positive_column_command)
