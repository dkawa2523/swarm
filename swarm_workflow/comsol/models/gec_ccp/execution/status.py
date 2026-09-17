"""Canonical status serialization and final acceptance enforcement."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Iterable

from swarm_workflow._io import write_json

from ..closure import _restricted_gradient_response_metadata
from ..contracts import GecCcpMapping, GecCcpQualityError, GecCcpWorkflowError
from .context import PostsolveAudits, PreparedExecution, RunQuality, SolverExecution


def status_contract_metadata(
    mapping: GecCcpMapping,
    *,
    gradient_response: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Return the result-role and closure scope required on every status."""

    response = (
        gradient_response
        if gradient_response is not None
        else _restricted_gradient_response_metadata(
            mapping,
            scalar_diffusion_audit=None,
        )
    )
    full_gradient = bool(
        isinstance(response, dict)
        and response.get("full_gradient_response_identified") is True
    )
    full_physical = bool(
        isinstance(response, dict)
        and response.get("full_physical_eligible") is True
    )
    return {
        "result_role": mapping.results.role,
        "external_swarm_source": mapping.bundle.expected_source,
        "mc_qualification_profile": (
            mapping.bundle.expected_mc_qualification_profile
        ),
        "electron_transport": mapping.closure.electron_transport,
        "quality_acceptance_scope": "declared_restricted_closure",
        "quality_accepted_for_declared_closure": None,
        "physical_target_accepted": False,
        "gradient_response": response,
        "full_gradient_response_identified": full_gradient,
        "full_physical_eligible": full_physical,
    }


def canonical_result_artifacts(
    output_directory: Path,
    expected_result_files: Iterable[Path],
) -> list[dict[str, Any]]:
    """Describe the exact COMSOL CSV set accepted by the run."""

    root = output_directory.resolve()
    artifacts: list[dict[str, Any]] = []
    for path in expected_result_files:
        resolved = path.resolve()
        try:
            relative = resolved.relative_to(root).as_posix()
        except ValueError as exc:
            raise GecCcpWorkflowError(
                f"expected result is outside the output directory: {path}"
            ) from exc
        data = resolved.read_bytes()
        artifacts.append(
            {
                "path": relative,
                "size_bytes": len(data),
                "sha256": hashlib.sha256(data).hexdigest(),
            }
        )
    return sorted(artifacts, key=lambda item: item["path"])


def write_final_status(
    prepared: PreparedExecution,
    solver: SolverExecution,
    audits: PostsolveAudits,
    quality: RunQuality,
    comsol_runtime: dict[str, Any],
) -> None:
    """Write the stable accepted/rejected GEC result artifact schema."""

    plan = prepared.plan
    result_artifacts = canonical_result_artifacts(
        plan.output_directory,
        plan.expected_result_files,
    )
    write_json(
        prepared.status_path,
        {
            "stage": "run-gec-ccp",
            "status": "completed" if quality.accepted else "rejected",
            "solve_status": "completed",
            "quality_status": quality.physics_status,
            "quality_accepted": quality.accepted,
            **status_contract_metadata(
                plan.mapping,
                gradient_response=quality.gradient_response,
            ),
            "quality_accepted_for_declared_closure": quality.accepted,
            "physical_target_accepted": bool(
                quality.accepted
                and plan.mapping.results.role == "physical_target"
            ),
            "full_gradient_response_identified": (
                quality.full_gradient_response_identified
            ),
            "full_physical_eligible": quality.full_physical_eligible,
            "gradient_response": quality.gradient_response,
            "numerical_quality_status": quality.numerical_status,
            "closure_numerical_quality_status": quality.numerical_status,
            "physics_quality_status": quality.physics_status,
            "mc_rate_censoring_quality_status": quality.mc_censoring_status,
            "mc_statistical_physics_qualification_status": (
                quality.mc_statistical_status
            ),
            "mc_rate_censoring": quality.mc_censoring,
            "elastic_energy_loss_model": (
                plan.mapping.closure.elastic_energy_loss_model
            ),
            "elastic_energy_loss_statistical_evidence": (
                quality.elastic_statistics
            ),
            "mc_statistical_physics_qualification": {
                "status": quality.mc_statistical_status,
                "upstream_active_closure_quality": (
                    quality.upstream_mc_quality_status
                ),
                "transport_runtime_audit": quality.transport_status,
                "rate_censoring": quality.mc_censoring_status,
                "elastic_energy_loss_uncertainty": (
                    quality.elastic_statistics_status
                ),
            },
            "rejection_code": (
                None
                if quality.accepted
                else (
                    "mc_censored_rate_physics_unqualified"
                    if quality.numerical_status == "passed"
                    and quality.mc_censoring_status == "failed"
                    else "physics_quality_failed"
                )
            ),
            "plan_sha256": prepared.plan_sha256,
            **prepared.input_provenance,
            "comsol_runtime": comsol_runtime,
            "output_mph": {
                "path": str(plan.mapping.model.output_mph),
                "sha256": hashlib.sha256(
                    plan.mapping.model.output_mph.read_bytes()
                ).hexdigest(),
                "size_bytes": plan.mapping.model.output_mph.stat().st_size,
            },
            "operations": [item.operation for item in solver.executions],
            "builtin_reference_performed": (
                plan.mapping.run.include_builtin_reference
            ),
            "power_W": plan.mapping.run.power_W,
            "coefficient_sweep": False,
            "conservation_audit": str(audits.conservation_path),
            "conservation_quality_status": quality.conservation_status,
            "builtin_reference_conservation_status": (
                quality.builtin_reference_conservation_status
            ),
            "transport_audit": (
                str(audits.transport_path)
                if audits.transport_path is not None
                else None
            ),
            "transport_quality_status": quality.transport_status,
            "transport_support_qualification": quality.transport_support,
            "closure_support_audit": (
                str(audits.closure_support_path)
                if audits.closure_support_path is not None
                else None
            ),
            "closure_support_quality_status": quality.closure_support_status,
            "external_closure_support_qualification": (
                quality.closure_support_qualification
            ),
            "function_eedf_audit": (
                str(audits.function_eedf_path)
                if audits.function_eedf_path is not None
                else None
            ),
            "function_eedf_quality_status": quality.function_eedf_status,
            "native_function_eedf_audit": (
                str(solver.native_eedf_audit_path)
                if solver.native_eedf_audit_path is not None
                else None
            ),
            "native_function_eedf_quality_status": (
                solver.native_eedf_status
            ),
            "results": result_artifacts,
        },
    )


def enforce_quality_acceptance(
    prepared: PreparedExecution,
    solver: SolverExecution,
    audits: PostsolveAudits,
    quality: RunQuality,
) -> None:
    """Raise the quality error after the rejected status is safely written."""

    if quality.accepted:
        return
    failed_audits = [
        name
        for name, audit_status in (
            ("conservation", quality.conservation_status),
            ("transport", quality.transport_status),
            ("closure_support", quality.closure_support_status),
            ("function_eedf", quality.function_eedf_status),
            ("native_function_eedf", solver.native_eedf_status),
            ("mc_rate_censoring", quality.mc_censoring_status),
        )
        if audit_status == "failed"
    ]
    error = GecCcpQualityError(
        "COMSOL solve completed, but physics-quality acceptance failed "
        f"({', '.join(failed_audits)}); see {prepared.status_path}"
    )
    if audits.transport_error is not None:
        raise error from audits.transport_error
    if audits.closure_support_error is not None:
        raise error from audits.closure_support_error
    if audits.function_eedf_error is not None:
        raise error from audits.function_eedf_error
    raise error


__all__ = [
    "canonical_result_artifacts",
    "enforce_quality_acceptance",
    "status_contract_metadata",
    "write_final_status",
]
