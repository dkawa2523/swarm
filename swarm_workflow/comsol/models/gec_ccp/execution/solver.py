"""Run the ordered COMSOL operations for a prepared GEC-CCP plan."""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

from swarm_workflow._io import write_json
from swarm_workflow._paths import discover_repo_root
from swarm_workflow.comsol.input.eedf_audit import (
    ComsolEedfAuditError,
    analyze_comsol_eedf_audit,
    extract_comsol_eedf_audit_log,
)
from swarm_workflow.comsol.runtime import (
    ComsolAdapterError,
    ComsolExecutionSummary,
    ComsolRunContext,
    execute_generated_comsol_java,
)

from ..contracts import GecCcpMapping, GecCcpWorkflowError
from .context import PreparedExecution, SolverExecution
from .inputs import validate_gec_plan_inputs
from .status import status_contract_metadata


def _comsol_run_context(mapping: GecCcpMapping) -> ComsolRunContext:
    return ComsolRunContext(
        root=mapping.root,
        mapping_path=mapping.path,
        input_mph=mapping.model.input_mph,
        output_mph=mapping.model.output_mph,
        log_path=mapping.logs.path,
        bundle_path=mapping.bundle.path,
    )


def run_solver_steps(
    prepared: PreparedExecution,
    *,
    comsol_executable: str | Path | None,
) -> SolverExecution:
    """Apply inputs, run/export target cases, and retain native EEDF evidence."""

    plan = prepared.plan
    # Validate, solve, and export the physical target first. The optional
    # built-in reference is intentionally outside the normal production path.
    steps = [("apply", plan.apply_java, "not_applicable")]
    if plan.native_eedf_audit is not None:
        steps.append(
            (
                "gec_native_eedf_audit",
                plan.native_eedf_audit.java_path,
                "not_applicable",
            )
        )
    steps.extend(
        (
            (
                "gec_external_run",
                plan.external_run_java,
                plan.mapping.model.external_study,
            ),
            (
                "gec_external_export",
                plan.external_export_java,
                plan.mapping.model.external_study,
            ),
        )
    )
    if plan.mapping.run.include_builtin_reference:
        if plan.baseline_run_java is None or plan.baseline_export_java is None:
            raise GecCcpWorkflowError("built-in reference Java was not prepared")
        steps.extend(
            (
                (
                    "gec_baseline_run",
                    plan.baseline_run_java,
                    plan.mapping.model.study,
                ),
                (
                    "gec_baseline_export",
                    plan.baseline_export_java,
                    plan.mapping.model.study,
                ),
            )
        )
    executions: list[ComsolExecutionSummary] = []
    native_eedf_audit_path: Path | None = None
    native_eedf_status = "not_applicable"
    try:
        for operation, java_file, executed_study in steps:
            if operation in {"apply", "gec_external_run"}:
                try:
                    validate_gec_plan_inputs(plan)
                except GecCcpWorkflowError as exc:
                    write_json(
                        prepared.status_path,
                        {
                            "stage": "run-gec-ccp",
                            "status": "blocked",
                            "blocker_code": "plan_input_provenance_mismatch",
                            "reason": str(exc),
                            "completed_operations": [
                                item.operation for item in executions
                            ],
                            "solve_status": "external_not_started",
                            "quality_status": "not_evaluated",
                            "plan": str(plan.plan_json),
                            "plan_sha256": prepared.plan_sha256,
                            **status_contract_metadata(
                                plan.mapping,
                                gradient_response=prepared.gradient_response,
                            ),
                            **prepared.input_provenance,
                        },
                    )
                    raise
            execution = execute_generated_comsol_java(
                _comsol_run_context(plan.mapping),
                java_file,
                operation=operation,
                comsol_executable=comsol_executable,
                study=executed_study,
                require_fresh_output=operation == "apply",
            )
            executions.append(execution)
            if operation == "gec_native_eedf_audit":
                try:
                    extract_comsol_eedf_audit_log(
                        plan.native_eedf_audit,
                        execution.stdout_paths[-1],
                    )
                    native_result = analyze_comsol_eedf_audit(
                        plan.native_eedf_audit,
                        cross_section_csv=(
                            discover_repo_root(__file__)
                            / "examples"
                            / "cross_sections"
                            / "argon_application_library.csv"
                        ),
                    )
                    native_eedf_status = (
                        "passed" if native_result["passed"] else "failed"
                    )
                    native_result["status"] = native_eedf_status
                    native_eedf_audit_path = (
                        plan.native_eedf_audit.values_path.parent
                        / "comsol_eedf_audit.json"
                    )
                    write_json(native_eedf_audit_path, native_result)
                except (ComsolEedfAuditError, OSError) as exc:
                    write_json(
                        prepared.status_path,
                        {
                            "stage": "run-gec-ccp",
                            "status": "blocked",
                            "blocker_code": (
                                "native_function_eedf_binding_audit_incomplete"
                            ),
                            "reason": str(exc),
                            "solve_status": "external_not_started",
                            "quality_status": "not_evaluated",
                            "completed_operations": [
                                item.operation for item in executions
                            ],
                            "plan": str(plan.plan_json),
                            "plan_sha256": prepared.plan_sha256,
                            **status_contract_metadata(
                                plan.mapping,
                                gradient_response=prepared.gradient_response,
                            ),
                            **prepared.input_provenance,
                        },
                    )
                    raise GecCcpWorkflowError(
                        "native COMSOL Function-EEDF binding audit is incomplete"
                    ) from exc
    except ComsolAdapterError as exc:
        write_json(
            prepared.status_path,
            {
                "stage": "run-gec-ccp",
                "status": "blocked",
                "blocker_code": _comsol_blocker_code(exc),
                "reason": str(exc),
                "completed_operations": [item.operation for item in executions],
                "solve_status": "blocked",
                "quality_status": "not_evaluated",
                "native_function_eedf_quality_status": native_eedf_status,
                "native_function_eedf_audit": (
                    str(native_eedf_audit_path)
                    if native_eedf_audit_path is not None
                    else None
                ),
                "coefficient_sweep": False,
                "model_contract": asdict(plan.contract),
                "plan": str(plan.plan_json),
                "plan_sha256": prepared.plan_sha256,
                **status_contract_metadata(
                    plan.mapping,
                    gradient_response=prepared.gradient_response,
                ),
                **prepared.input_provenance,
            },
        )
        raise
    return SolverExecution(
        executions=tuple(executions),
        native_eedf_audit_path=native_eedf_audit_path,
        native_eedf_status=native_eedf_status,
    )


def _comsol_blocker_code(exc: ComsolAdapterError) -> str:
    result = exc.step_result or {}
    text_parts: list[str] = []
    for field in ("stdout", "stderr"):
        value = result.get(field)
        if not isinstance(value, str):
            continue
        path = Path(value)
        if path.exists():
            text_parts.append(
                path.read_text(encoding="utf-8", errors="replace").lower()
            )
    text = "\n".join(text_parts)
    if "license has expired" in text or "product has expired" in text:
        return "comsol_license_expired"
    return "comsol_execution_failed"


__all__ = ["run_solver_steps"]
