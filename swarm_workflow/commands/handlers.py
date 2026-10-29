"""Thin CLI handlers that translate arguments into workflow calls."""

from __future__ import annotations

import argparse
import json
import sqlite3

from electron_swarm.runner import run
from electron_swarm.orchestration.preparation import prepare_run
from swarm_workflow.campaign.document import read_calculation
from typing import Any

from swarm_workflow.campaign.aggregate import aggregate_database
from swarm_workflow.campaign.sweep import (
    advance_monte_carlo_workflow,
    run_sweep,
    plan_sweep,
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
from electron_swarm.core.output import OutputConfig
from swarm_workflow.reaction_rates import export_reaction_rates
from swarm_workflow.figures import plot_results
from swarm_workflow.campaign.store import WorkflowSchemaError


def run_command(args: Any, parser: argparse.ArgumentParser) -> None:
    if args.json and not args.plan:
        parser.error("--json requires --plan")
    try:
        document = read_calculation(args.config)
        if document.campaign is not None:
            if args.plan:
                description = plan_sweep(args.config)
            else:
                summary = run_sweep(args.config)
                partial = summary.status == "partial"
                status = "partial" if partial else "completed"
                print(
                    f"{status} {summary.cases_written} new cases; database={summary.database_path}"
                )
                if document.template.output.export is not None:
                    print(
                        f"reaction_rates: {document.template.output.directory / (document.template.output.base_name + '_manifest.json')}"
                    )
                if document.template.output.figures:
                    print(
                        f"figures: {document.template.output.directory / (document.template.output.base_name + '_manifest.json')}"
                    )
                if partial:
                    raise SystemExit(4)
                return
        else:
            if len(document.mixtures) != 1:
                raise ValueError(
                    "multiple mixtures require campaign.database for explicit persistence and resume"
                )
            if args.plan:
                description = prepare_run(document.template).description()
            else:
                result = run(document.template)
                outcome = result.metadata["outcome"]
                print(f"{outcome['status']}: {len(result.cases)} solver cases")
                if document.template.output.export is not None:
                    print(
                        f"reaction_rates: {result.metadata['output_paths']['manifest_json']}"
                    )
                if document.template.output.figures:
                    print(
                        f"figures: {result.metadata['output_paths']['manifest_json']}"
                    )
                if outcome["status"] != "complete":
                    raise SystemExit({"partial": 4, "skipped": 3}[outcome["status"]])
                return
    except (ValueError, OSError, NotImplementedError) as exc:
        parser.error(str(exc))
    if args.json:
        print(json.dumps(description, indent=2))
        return
    print(
        f"{description['status']}: {description['expected_case_count']} planned solver cases (before resume)"
    )
    if description.get("execution") == "campaign":
        print(
            f"mixtures={len(description['mixtures'])}; E/N points={len(description['e_over_n_Td'])}; MC cases={description['monte_carlo_cases']}"
        )
        print(f"database: {description['database']}")
        items = [item for plan in description["solver_plan"].values() for item in plan]
    else:
        items = description["solver_plan"]
    for item in items:
        print(
            f"  {item['solver']}: {item['state']}"
            + (f"; {item['warnings']}" if item.get("warnings") else "")
        )
    if description.get("export"):
        export = description["export"]
        print(
            f"export: reaction_rates from {export['solver']}; {export['directory']}/{export['csv']}"
        )
    if description.get("figures"):
        print(
            f"figures: {', '.join(description['figures']['requested'])}; PNG/SVG and numerical CSV"
        )
    print(
        "Use --plan --json for resolved settings, physical cases, and resource details."
    )


def plot_command(args: Any, parser: argparse.ArgumentParser) -> None:
    try:
        paths = plot_results(
            args.source,
            OutputConfig(args.output, args.base_name, figures=tuple(args.figures)),
        )
    except (
        ValueError,
        OSError,
        NotImplementedError,
        sqlite3.Error,
        WorkflowSchemaError,
    ) as exc:
        parser.error(str(exc))
    print(f"figures: {paths['manifest_json']}")


def export_reaction_rates_command(args: Any, parser: argparse.ArgumentParser) -> None:
    try:
        paths = export_reaction_rates(
            args.source, OutputConfig(args.output, args.base_name), solver=args.solver
        )
    except (
        ValueError,
        OSError,
        NotImplementedError,
        sqlite3.Error,
        WorkflowSchemaError,
    ) as exc:
        parser.error(str(exc))
    print(f"reaction_rates: {paths['manifest_json']}")


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
        f"wrote {len(summary.figures)} GEC CCP figure(s); manifest: {summary.manifest}"
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
            args.external_results,
            args.reference_results,
            external_solution=args.external_solution,
            reference_solution=args.reference_solution,
            output_dir=args.output,
            external_runtime_s=args.external_runtime,
            reference_runtime_s=args.reference_runtime,
            create_plot=args.plot,
        )
    except ComsolCompareError as exc:
        parser.error(str(exc))
    print(format_comsol_comparison(summary))
