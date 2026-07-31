"""Command line entry point for external swarm workflows."""

from __future__ import annotations

import argparse
from pathlib import Path

from .aggregate import aggregate_database
from .comsol_adapter import ComsolAdapterError
from .comsol_export import export_comsol_bundle
from .comsol_mapping import ComsolMappingError
from .comsol_positive_column import (
    PositiveColumnWorkflowError,
    execute_positive_column_run,
    format_positive_column_plan,
    format_positive_column_summary,
    prepare_positive_column_run,
)
from .comsol_gec_ccp import (
    GecCcpWorkflowError,
    execute_gec_ccp_run,
    format_gec_ccp_plan,
    format_gec_ccp_summary,
    prepare_gec_ccp_run,
)
from .gec_ccp_plots import GecCcpPlotError, plot_gec_ccp_results
from .comsol_compare import (
    ComsolCompareError,
    compare_comsol_profiles,
    format_comsol_comparison,
)
from .comsol_verify import ComsolVerifyError
from .sweep import run_sweep
from .tables import SOURCE_CHOICES, build_tables


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="swarm-workflow")
    subparsers = parser.add_subparsers(dest="command", required=True)
    sweep_parser = subparsers.add_parser(
        "sweep",
        help="Run an external E/N and mixture sweep into SQLite",
    )
    sweep_parser.add_argument("workflow", type=Path, help="workflow YAML file")
    aggregate_parser = subparsers.add_parser(
        "aggregate",
        help="Aggregate workflow SQLite results into CSV and manifest files",
    )
    aggregate_parser.add_argument("database", type=Path, help="workflow SQLite database")
    aggregate_parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="aggregate output directory",
    )
    tables_parser = subparsers.add_parser(
        "build-tables",
        help="Build per-mixture COMSOL-ready workflow table directories",
    )
    tables_parser.add_argument("database", type=Path, help="workflow SQLite database")
    tables_parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="table output directory",
    )
    tables_parser.add_argument(
        "--source",
        required=True,
        choices=SOURCE_CHOICES,
        help="table source policy",
    )
    export_parser = subparsers.add_parser(
        "export-comsol",
        help="Export built workflow tables as COMSOL CSV bundles",
    )
    export_parser.add_argument(
        "table_directory",
        type=Path,
        help="build-tables root directory or one mixture table directory",
    )
    export_parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="output bundle directory",
    )
    run_parser = subparsers.add_parser(
        "run-comsol",
        help="Apply, verify, run, and export a COMSOL model with Swarm tables",
    )
    run_parser.add_argument("mapping", type=Path, help="COMSOL mapping YAML file")
    run_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="validate and generate Java without executing COMSOL",
    )
    run_parser.add_argument(
        "--comsol",
        type=Path,
        help="COMSOL executable path used for batch execution",
    )
    run_parser.add_argument(
        "--bundle",
        required=True,
        type=Path,
        help="export-comsol bundle directory to import",
    )
    run_parser.add_argument(
        "--pressure-Pa",
        type=float,
        help="explicit positive-column gas pressure in Pa",
    )
    run_parser.add_argument(
        "--gas-temperature-K",
        type=float,
        help="explicit positive-column gas temperature in K",
    )
    run_parser.add_argument(
        "--mesh-elements",
        type=int,
        choices=(200, 400),
        help="explicit comp1/mesh1/edg1/dis1 element count",
    )
    gec_parser = subparsers.add_parser(
        "run-gec-ccp",
        help="Prepare or run the time-periodic argon GEC CCP comparison",
    )
    gec_parser.add_argument("mapping", type=Path, help="GEC CCP mapping YAML file")
    gec_parser.add_argument(
        "--dry-run",
        action="store_true",
        help="inspect the MPH, validate tables, and generate Java only",
    )
    gec_parser.add_argument(
        "--comsol",
        type=Path,
        help="COMSOL executable path used for batch execution",
    )
    gec_parser.add_argument(
        "--bundle",
        required=True,
        type=Path,
        help="export-comsol bundle directory to import",
    )
    gec_parser.add_argument(
        "--reuse-baseline",
        action="store_true",
        help="reuse an already converged direct baseline MPH",
    )
    plot_gec_parser = subparsers.add_parser(
        "plot-gec-ccp",
        help="Plot GEC CCP Swarm inputs and available COMSOL comparisons",
    )
    plot_gec_parser.add_argument(
        "--bundle",
        required=True,
        type=Path,
        help="export-comsol bundle directory",
    )
    plot_gec_parser.add_argument(
        "--comsol-results",
        type=Path,
        help="directory containing builtin_druyvesteyn and swarm_tables results",
    )
    plot_gec_parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="plot output directory",
    )
    benchmark_parser = subparsers.add_parser(
        "compare-comsol",
        help="Compare external-table and built-in-Boltzmann result CSV files",
    )
    benchmark_parser.add_argument(
        "external_csv",
        type=Path,
        help="external Swarm-table result CSV",
    )
    benchmark_parser.add_argument(
        "reference_csv",
        type=Path,
        help="built-in Boltzmann reference result CSV",
    )
    benchmark_parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="comparison output directory",
    )
    benchmark_parser.add_argument(
        "--external-runtime",
        type=float,
        help="optional external solve time in seconds",
    )
    benchmark_parser.add_argument(
        "--reference-runtime",
        type=float,
        help="optional built-in reference solve time in seconds",
    )
    benchmark_parser.add_argument(
        "--plot",
        action="store_true",
        help="write a PNG comparison of the spatial profiles",
    )
    args = parser.parse_args(argv)
    if args.command == "sweep":
        summary = run_sweep(args.workflow)
        print(
            "completed "
            f"{summary.cases_written} cases "
            f"({summary.mixtures} mixtures, {summary.e_over_n_values} E/N values)"
        )
        return
    if args.command == "aggregate":
        summary = aggregate_database(args.database, args.output)
        print(
            "aggregated "
            f"{summary.e_over_n_points} solver/E/N point(s) "
            f"from {summary.database_path}"
        )
        return
    if args.command == "build-tables":
        summary = build_tables(args.database, args.output, source=args.source)
        print(
            "built "
            f"{summary.mixtures} table director"
            f"{'y' if summary.mixtures == 1 else 'ies'} "
            f"from {summary.database_path}"
        )
        return
    if args.command == "export-comsol":
        summary = export_comsol_bundle(
            args.table_directory,
            args.output,
        )
        print(
            "exported "
            f"{summary.bundles_written} COMSOL bundle(s) "
            f"from {summary.table_directory}"
        )
        return
    if args.command == "run-comsol":
        try:
            if args.dry_run:
                print(
                    format_positive_column_plan(
                        prepare_positive_column_run(
                            args.mapping,
                            bundle_path=args.bundle,
                            write_java=True,
                            pressure_Pa=args.pressure_Pa,
                            gas_temperature_K=args.gas_temperature_K,
                            mesh_elements=args.mesh_elements,
                        )
                    )
                )
                return
            summary = execute_positive_column_run(
                args.mapping,
                bundle_path=args.bundle,
                comsol_executable=args.comsol,
                pressure_Pa=args.pressure_Pa,
                gas_temperature_K=args.gas_temperature_K,
                mesh_elements=args.mesh_elements,
            )
            print(format_positive_column_summary(summary))
            return
        except (
            ComsolMappingError,
            ComsolAdapterError,
            ComsolVerifyError,
            PositiveColumnWorkflowError,
        ) as exc:
            parser.error(str(exc))
    if args.command == "run-gec-ccp":
        try:
            if args.dry_run:
                print(
                    format_gec_ccp_plan(
                        prepare_gec_ccp_run(
                            args.mapping,
                            bundle_path=args.bundle,
                            write_java=True,
                        )
                    )
                )
                return
            print(
                format_gec_ccp_summary(
                    execute_gec_ccp_run(
                        args.mapping,
                        bundle_path=args.bundle,
                        comsol_executable=args.comsol,
                        reuse_baseline=args.reuse_baseline,
                    )
                )
            )
            return
        except (GecCcpWorkflowError, ComsolAdapterError) as exc:
            parser.error(str(exc))
    if args.command == "plot-gec-ccp":
        try:
            summary = plot_gec_ccp_results(
                args.bundle,
                output_dir=args.output,
                comsol_results_dir=args.comsol_results,
            )
            print(
                "wrote "
                f"{len(summary.figures)} GEC CCP figure(s); "
                f"manifest: {summary.manifest}"
            )
            return
        except GecCcpPlotError as exc:
            parser.error(str(exc))
    if args.command == "compare-comsol":
        try:
            summary = compare_comsol_profiles(
                args.external_csv,
                args.reference_csv,
                output_dir=args.output,
                external_runtime_s=args.external_runtime,
                reference_runtime_s=args.reference_runtime,
                create_plot=args.plot,
            )
            print(format_comsol_comparison(summary))
            return
        except ComsolCompareError as exc:
            parser.error(str(exc))
    parser.error(f"unsupported command: {args.command}")


if __name__ == "__main__":
    main()
