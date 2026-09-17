"""Execute and accept one prepared GEC-ICP COMSOL target."""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any

from swarm_workflow._io import write_json
from swarm_workflow.comsol.input.eedf_audit import (
    analyze_comsol_eedf_audit,
    extract_comsol_eedf_audit_log,
)
from swarm_workflow.comsol.runtime import (
    ComsolAdapterError,
    ComsolExecutionSummary,
    ComsolRunContext,
    execute_generated_comsol_java,
)

from .closure_readback import audit_gec_icp_closure_readback
from .convergence import assess_gec_icp_convergence
from .prepare import prepare_gec_icp_run
from .run_contracts import (
    GecIcpPlan,
    GecIcpQualityError,
    GecIcpRunSummary,
    GecIcpWorkflowError,
)


def execute_gec_icp_run(
    mapping_path: str | Path,
    *,
    bundle_path: str | Path | None = None,
    comsol_executable: str | Path | None = None,
) -> GecIcpRunSummary:
    """Apply one solver bundle, cold-solve one ICP case, and enforce quality."""

    plan = prepare_gec_icp_run(mapping_path, bundle_path=bundle_path)
    status_path = plan.output_directory / "run_status.json"
    _clear_result_files(plan)
    executions: list[ComsolExecutionSummary] = []
    closure_readback: dict[str, Any] | None = None
    native_eedf: dict[str, Any] | None = None
    convergence: dict[str, Any] | None = None
    solve_status = "not_started"
    write_json(
        status_path,
        _status_payload(plan, status="running", solve_status=solve_status),
    )
    try:
        closure_result = execute_generated_comsol_java(
            _context(plan),
            plan.closure_java,
            operation="gec_icp_closure",
            comsol_executable=comsol_executable,
            study="not_applicable",
            require_fresh_output=True,
            support_java_paths=plan.closure_support_java,
        )
        executions.append(closure_result)
        closure_readback = audit_gec_icp_closure_readback(
            plan.closure_readback,
            closure_result.stdout_paths[-1],
        )
        if not closure_readback["passed"]:
            raise GecIcpQualityError(
                "saved COMSOL closure settings failed post-apply readback"
            )
        extract_comsol_eedf_audit_log(
            plan.native_eedf_audit, closure_result.stdout_paths[-1]
        )
        native_eedf = analyze_comsol_eedf_audit(
            plan.native_eedf_audit,
            cross_section_csv=(
                _repo_root(plan.mapping.path)
                / "examples"
                / "cross_sections"
                / "argon_application_library.csv"
            ),
        )
        native_eedf["status"] = "passed" if native_eedf["passed"] else "failed"
        native_path = plan.output_directory / "comsol_eedf_audit.json"
        write_json(native_path, native_eedf)
        if not native_eedf["passed"]:
            raise GecIcpQualityError(
                "saved COMSOL Function-EEDF binding failed native evaluation"
            )

        solve_status = "running"
        solve_result = execute_generated_comsol_java(
            _context(plan),
            plan.solve_java,
            operation="gec_icp_solve",
            comsol_executable=comsol_executable,
            study=plan.mapping.model_mapping.model.study,
            require_fresh_output=True,
        )
        executions.append(solve_result)
        solve_status = "completed"
        convergence = assess_gec_icp_convergence(plan)
        accepted = bool(
            closure_readback["passed"]
            and native_eedf["passed"]
            and convergence["passed"]
        )
        final = _status_payload(
            plan,
            status="passed" if accepted else "failed",
            executions=executions,
            closure_readback=closure_readback,
            native_eedf=native_eedf,
            convergence=convergence,
            solve_status=solve_status,
            quality_accepted_for_declared_closure=accepted,
            physical_target_accepted=accepted,
        )
        write_json(status_path, final)
        if not accepted:
            failed = [
                name
                for name, gate in convergence["gates"].items()
                if not gate["passed"]
            ]
            raise GecIcpQualityError(
                "GEC ICP solve completed but convergence failed: " + ", ".join(failed)
            )
    except GecIcpQualityError as exc:
        write_json(
            status_path,
            _status_payload(
                plan,
                status="failed",
                executions=executions,
                closure_readback=closure_readback,
                native_eedf=native_eedf,
                convergence=convergence,
                solve_status=solve_status,
                reason=str(exc),
                quality_accepted_for_declared_closure=False,
                physical_target_accepted=False,
            ),
        )
        raise
    except (ComsolAdapterError, GecIcpWorkflowError, OSError, ValueError) as exc:
        write_json(
            status_path,
            _status_payload(
                plan,
                status="blocked",
                executions=executions,
                closure_readback=closure_readback,
                native_eedf=native_eedf,
                convergence=convergence,
                solve_status=solve_status,
                reason=str(exc),
                quality_accepted_for_declared_closure=False,
                physical_target_accepted=False,
            ),
        )
        raise
    return GecIcpRunSummary(plan, tuple(executions), status_path)


def _context(plan: GecIcpPlan) -> ComsolRunContext:
    mapping = plan.mapping
    return ComsolRunContext(
        root=mapping.root,
        mapping_path=mapping.path,
        input_mph=mapping.model_mapping.model.input_mph,
        output_mph=mapping.run.output_mph,
        log_path=mapping.log_path,
        bundle_path=mapping.bundle.path,
    )


def _clear_result_files(plan: GecIcpPlan) -> None:
    """Prevent a successful process from inheriting stale CSV evidence."""

    for path in plan.expected_result_files:
        if path.is_file():
            path.unlink()
    for path in (
        plan.native_eedf_audit.values_path,
        plan.native_eedf_audit.contract_path,
    ):
        if path.is_file():
            path.unlink()


def _status_payload(
    plan: GecIcpPlan,
    *,
    status: str,
    executions: list[ComsolExecutionSummary] | None = None,
    closure_readback: dict[str, Any] | None = None,
    native_eedf: dict[str, Any] | None = None,
    convergence: dict[str, Any] | None = None,
    solve_status: str = "not_started",
    quality_accepted_for_declared_closure: bool | None = None,
    physical_target_accepted: bool = False,
    reason: str | None = None,
) -> dict[str, Any]:
    source = plan.mapping.bundle.expected_source
    external_input_handling = (
        "external_swarm_complete_anchor_selection"
        if source == "composite"
        else "external_swarm_direct"
    )
    return {
        "schema": "swarm.gec_icp_comsol_run.v1",
        "stage": "run-gec-icp",
        "status": status,
        "reason": reason,
        "model_adapter": "gec_icp",
        "source": source,
        "result_role": "physical_target",
        "quality_acceptance_scope": "declared_restricted_lmea_closure",
        "quality_accepted_for_declared_closure": (
            quality_accepted_for_declared_closure
        ),
        "physical_target_accepted": physical_target_accepted,
        "solve_status": solve_status,
        "plan": str(plan.plan_json),
        "output_mph": str(plan.mapping.run.output_mph),
        "coefficient_sweep": False,
        "closure": asdict(plan.mapping.closure),
        "physics_feature_handling": {
            "function_eedf": external_input_handling,
            "reduced_mobility": external_input_handling,
            "anchor_fallback": (
                "qualified_monte_carlo_else_exact_qualified_two_term_complete_anchor"
                if source == "composite"
                else "not_requested"
            ),
            "particle_diffusion": "comsol_SpecifyMueOnly_Einstein_closure",
            "electron_energy_transport": "comsol_local_mean_energy",
            "ground_state_reactions": "comsol_integral_of_external_function_eedf",
            "metastable_superelastic_and_stepwise_ionization": (
                "comsol_embedded_cross_sections_with_external_function_eedf"
            ),
            "rf_and_magnetic_field": "comsol_frequency_transient",
            "homogeneous_solver_field": "steady_dc_local_mean_energy_closure",
        },
        "closure_readback": closure_readback,
        "native_function_eedf": native_eedf,
        "convergence": convergence,
        "expected_artifacts": {
            "closure_readback": str(plan.closure_readback.report_path),
            "native_function_eedf": str(
                plan.output_directory / "comsol_eedf_audit.json"
            ),
            "convergence": str(plan.output_directory / "convergence.json"),
        },
        "executions": [
            {
                "operation": item.operation,
                "total_time_s": item.total_time_s,
                "log_dir": str(item.log_dir),
                "result": str(item.result_json),
                "provenance": str(item.provenance_json),
                "comsol_version": item.comsol_version,
                "comsol_build": item.comsol_build,
            }
            for item in (executions or [])
        ],
    }


def _repo_root(start: Path) -> Path:
    for candidate in (start.parent, *start.parents):
        if (candidate / "pyproject.toml").is_file():
            return candidate
    raise GecIcpWorkflowError("could not locate repository root")


__all__ = ["execute_gec_icp_run"]
