"""Thin CLI handlers that translate arguments into workflow calls."""

from __future__ import annotations

import argparse
from typing import Any

from swarm_workflow.campaign.aggregate import aggregate_database
from swarm_workflow.campaign.sweep import (
    advance_monte_carlo_workflow,
    run_sweep,
)
from swarm_workflow.comsol.input import ComsolExportError, export_comsol_bundle
from swarm_workflow.selection import ClosureSelectionError
from swarm_workflow.comsol.models.gec_ccp import (
    GecCcpWorkflowError,
    execute_gec_ccp_run,
    format_gec_ccp_plan,
    format_gec_ccp_summary,
    prepare_gec_ccp_run,
)
from swarm_workflow.comsol.models.gec_ccp.plots.workflow import (
    GecCcpPlotError,
    plot_gec_ccp_results,
    plot_gec_ccp_solver_comparison,
)
from swarm_workflow.comsol.models.gec_icp import (
    GecIcpWorkflowError,
    execute_gec_icp_run,
    format_gec_icp_plan,
    format_gec_icp_summary,
    prepare_gec_icp_run,
)
from swarm_workflow.comsol.models.positive_column import (
    ComsolCompareError,
    ComsolMappingError,
    ComsolVerifyError,
    PositiveColumnWorkflowError,
    compare_comsol_profiles,
    execute_positive_column_run,
    format_comsol_comparison,
    format_positive_column_plan,
    format_positive_column_summary,
    prepare_positive_column_run,
)
from swarm_workflow.comsol.runtime import ComsolAdapterError
from swarm_workflow.quality.monte_carlo.policy import (
    MonteCarloPolicyError,
    decide_monte_carlo_closure,
)
from swarm_workflow.tables import build_tables


def run_sweep_command(args: Any, _parser: argparse.ArgumentParser) -> None:
    summary = run_sweep(args.workflow)
    print(
        "completed "
        f"{summary.cases_written} cases "
        f"({summary.mixtures} mixtures, {summary.e_over_n_values} E/N values)"
    )


def aggregate_command(args: Any, _parser: argparse.ArgumentParser) -> None:
    summary = aggregate_database(args.database, args.output)
    print(
        "aggregated "
        f"{summary.e_over_n_points} solver/E/N point(s) "
        f"from {summary.database_path}"
    )


def decide_mc_command(args: Any, parser: argparse.ArgumentParser) -> None:
    try:
        summary = decide_monte_carlo_closure(
            args.mc_tables,
            args.two_term_tables,
            attempt=args.attempt,
            output=args.output,
            previous_decision=args.previous_decision,
        )
    except (MonteCarloPolicyError, ClosureSelectionError) as exc:
        parser.error(str(exc))
    selected = summary.selected_solver or "pending"
    print(
        f"MC policy action={summary.action}; selected_solver={selected}; "
        f"artifact={summary.artifact}"
    )


def advance_mc_command(args: Any, parser: argparse.ArgumentParser) -> None:
    try:
        output = advance_monte_carlo_workflow(
            args.workflow,
            args.decision,
            output_path=args.output,
            database_path=args.database,
        )
    except (ValueError, MonteCarloPolicyError) as exc:
        parser.error(str(exc))
    print(f"wrote bounded MC workflow: {output}")


def build_tables_command(args: Any, _parser: argparse.ArgumentParser) -> None:
    summary = build_tables(
        args.database,
        args.output,
        source=args.source,
        mc_qualification_profile=args.mc_qualification_profile,
        allow_unqualified_mc=args.allow_unqualified_mc,
        solver_qualification_path=args.solver_qualification,
        target_qualification_path=args.target_qualification,
    )
    print(
        "built "
        f"{summary.mixtures} table director"
        f"{'y' if summary.mixtures == 1 else 'ies'} "
        f"from {summary.database_path}"
    )
    if summary.unqualified_mixtures:
        print(
            "MC evidence only, without coefficient tables: "
            f"{summary.unqualified_mixtures} mixture(s)"
        )


def export_comsol_command(args: Any, parser: argparse.ArgumentParser) -> None:
    try:
        summary = export_comsol_bundle(
            args.table_directory,
            args.output,
            selection_path=args.selection,
        )
    except ComsolExportError as exc:
        parser.error(str(exc))
    print(
        "exported "
        f"{summary.bundles_written} COMSOL bundle(s) "
        f"from {summary.table_directory}"
    )


def run_positive_column_command(
    args: Any,
    parser: argparse.ArgumentParser,
) -> None:
    try:
        if args.dry_run:
            plan = prepare_positive_column_run(
                args.mapping,
                bundle_path=args.bundle,
                write_java=True,
                pressure_Pa=args.pressure_Pa,
                gas_temperature_K=args.gas_temperature_K,
                mesh_elements=args.mesh_elements,
            )
            print(format_positive_column_plan(plan))
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
    except (
        ComsolMappingError,
        ComsolAdapterError,
        ComsolVerifyError,
        PositiveColumnWorkflowError,
    ) as exc:
        parser.error(str(exc))


def run_gec_ccp_command(args: Any, parser: argparse.ArgumentParser) -> None:
    try:
        if args.dry_run:
            plan = prepare_gec_ccp_run(
                args.mapping,
                bundle_path=args.bundle,
                write_java=True,
            )
            print(format_gec_ccp_plan(plan))
            return
        summary = execute_gec_ccp_run(
            args.mapping,
            bundle_path=args.bundle,
            comsol_executable=args.comsol,
        )
        print(format_gec_ccp_summary(summary))
    except (GecCcpWorkflowError, ComsolAdapterError) as exc:
        parser.error(str(exc))


def run_gec_icp_command(args: Any, parser: argparse.ArgumentParser) -> None:
    try:
        if args.dry_run:
            plan = prepare_gec_icp_run(
                args.mapping,
                bundle_path=args.bundle,
                write_java=True,
            )
            print(format_gec_icp_plan(plan))
            return
        summary = execute_gec_icp_run(
            args.mapping,
            bundle_path=args.bundle,
            comsol_executable=args.comsol,
        )
        print(format_gec_icp_summary(summary))
    except (GecIcpWorkflowError, ComsolAdapterError) as exc:
        parser.error(str(exc))


def plot_gec_ccp_command(args: Any, parser: argparse.ArgumentParser) -> None:
    try:
        summary = plot_gec_ccp_results(
            args.bundle,
            output_dir=args.output,
            comsol_results_dir=args.comsol_results,
        )
    except GecCcpPlotError as exc:
        parser.error(str(exc))
    print(
        f"wrote {len(summary.figures)} GEC CCP figure(s); "
        f"manifest: {summary.manifest}"
    )


def compare_gec_ccp_command(args: Any, parser: argparse.ArgumentParser) -> None:
    try:
        summary = plot_gec_ccp_solver_comparison(
            args.two_term_bundle,
            args.monte_carlo_bundle,
            two_term_results_dir=args.two_term_results,
            monte_carlo_results_dir=args.monte_carlo_results,
            output_dir=args.output,
        )
    except GecCcpPlotError as exc:
        parser.error(str(exc))
    print(
        f"wrote {len(summary.figures)} GEC CCP solver comparison figure(s); "
        f"manifest: {summary.manifest}"
    )


def compare_positive_column_command(
    args: Any,
    parser: argparse.ArgumentParser,
) -> None:
    try:
        summary = compare_comsol_profiles(
            args.external_csv,
            args.reference_csv,
            output_dir=args.output,
            external_runtime_s=args.external_runtime,
            reference_runtime_s=args.reference_runtime,
            create_plot=args.plot,
        )
    except ComsolCompareError as exc:
        parser.error(str(exc))
    print(format_comsol_comparison(summary))
