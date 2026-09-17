"""Run independent GEC audits and derive the fail-closed quality decision."""

from __future__ import annotations

import json
from typing import Any

from ..audits.conservation import audit_gec_ccp_conservation_run
from ..audits.function_eedf import audit_gec_ccp_function_eedf_run
from ..audits.transport_run import (
    audit_gec_ccp_closure_support_run,
    audit_gec_ccp_transport_run,
)
from ..closure import (
    _uses_external_elastic_energy_loss,
    _uses_function_eedf,
    _uses_transport,
)
from ..contracts import GecCcpPlan, GecCcpWorkflowError
from .context import PostsolveAudits, PreparedExecution, RunQuality, SolverExecution


def run_postsolve_audits(plan: GecCcpPlan) -> PostsolveAudits:
    """Run every audit applicable to the plan while retaining failure causes."""

    conservation_audit = audit_gec_ccp_conservation_run(plan)
    transport_audit = None
    transport_error = None
    if _uses_transport(plan.mapping.closure):
        try:
            transport_audit = audit_gec_ccp_transport_run(plan)
        except GecCcpWorkflowError as exc:
            transport_error = exc
            candidate = plan.output_directory / "transport_audit.json"
            if candidate.exists():
                transport_audit = candidate
    closure_support_audit = None
    closure_support_error = None
    if (
        not _uses_transport(plan.mapping.closure)
        and plan.mapping.closure.reaction_model == "external_rates"
    ):
        try:
            closure_support_audit = audit_gec_ccp_closure_support_run(plan)
        except GecCcpWorkflowError as exc:
            closure_support_error = exc
            candidate = plan.output_directory / "closure_support_audit.json"
            if candidate.exists():
                closure_support_audit = candidate
    function_eedf_audit = None
    function_eedf_error = None
    if _uses_function_eedf(plan.mapping.closure):
        try:
            function_eedf_audit = audit_gec_ccp_function_eedf_run(plan)
        except GecCcpWorkflowError as exc:
            function_eedf_error = exc
            candidate = plan.output_directory / "function_eedf_audit.json"
            if candidate.exists():
                function_eedf_audit = candidate
    return PostsolveAudits(
        conservation_path=conservation_audit,
        transport_path=transport_audit,
        closure_support_path=closure_support_audit,
        function_eedf_path=function_eedf_audit,
        transport_error=transport_error,
        closure_support_error=closure_support_error,
        function_eedf_error=function_eedf_error,
    )


def evaluate_run_quality(
    prepared: PreparedExecution,
    solver: SolverExecution,
    audits: PostsolveAudits,
) -> RunQuality:
    """Combine run audits with frozen upstream physics qualification."""

    plan = prepared.plan
    conservation_payload = json.loads(
        audits.conservation_path.read_text(encoding="utf-8")
    )
    conservation_status = conservation_payload["status"]
    builtin_reference_conservation_status = (
        conservation_payload.get("cases", {})
        .get("baseline", {})
        .get("status", "not_requested")
    )
    transport_status = (
        json.loads(audits.transport_path.read_text(encoding="utf-8"))["status"]
        if audits.transport_path is not None
        else "not_applicable"
    )
    closure_support_status = (
        json.loads(audits.closure_support_path.read_text(encoding="utf-8"))[
            "status"
        ]
        if audits.closure_support_path is not None
        else (
            "failed"
            if audits.closure_support_error is not None
            else "not_applicable"
        )
    )
    closure_support_qualification: dict[str, Any] | None = None
    if audits.closure_support_path is not None:
        closure_support_payload = json.loads(
            audits.closure_support_path.read_text(encoding="utf-8")
        ).get("support")
        if isinstance(closure_support_payload, dict):
            closure_support_qualification = closure_support_payload
    transport_support: dict[str, Any] | None = None
    if audits.transport_path is not None:
        transport_payload = json.loads(
            audits.transport_path.read_text(encoding="utf-8")
        )
        candidate_support = transport_payload.get("transport", {}).get(
            "operating_mean_energy_range"
        )
        if isinstance(candidate_support, dict):
            guard = candidate_support.get("low_energy_guard_influence")
            transport_support = {
                "status": (
                    "passed"
                    if candidate_support.get("passed") is True
                    else "failed"
                ),
                "policy": candidate_support.get("support_policy", "strict"),
                "acceptance_basis": candidate_support.get(
                    "support_acceptance_basis", "strict_primary_support"
                ),
                "strict_primary_support_passed": candidate_support.get(
                    "strict_primary_support_passed",
                    candidate_support.get("passed") is True,
                ),
                "constant_extrapolation_accepted": candidate_support.get(
                    "constant_extrapolation_accepted", False
                ),
                "phase_local_mean_energy_eV": candidate_support.get(
                    "phase_local_mean_energy_eV"
                ),
                "primary_support_mean_energy_eV": candidate_support.get(
                    "common_intersection_mean_energy_eV"
                ),
                "guard_maximum_relative_influence": (
                    guard.get("maximum_relative_influence")
                    if isinstance(guard, dict)
                    else None
                ),
                "guard_relative_influence_limit": (
                    guard.get("relative_influence_limit")
                    if isinstance(guard, dict)
                    else None
                ),
                "scope_limit": (
                    guard.get("scope_limit")
                    if isinstance(guard, dict)
                    else None
                ),
            }
    function_status = (
        json.loads(audits.function_eedf_path.read_text(encoding="utf-8"))[
            "status"
        ]
        if audits.function_eedf_path is not None
        else (
            "failed"
            if audits.function_eedf_error is not None
            else "not_applicable"
        )
    )
    numerical_quality_status = (
        "passed"
        if conservation_status == "passed"
        and transport_status in {"passed", "not_applicable"}
        and closure_support_status in {"passed", "not_applicable"}
        and function_status in {"passed", "not_applicable"}
        and solver.native_eedf_status in {"passed", "not_applicable"}
        else "failed"
    )
    gradient_response = prepared.gradient_response
    full_gradient_response_identified = bool(
        isinstance(gradient_response, dict)
        and gradient_response.get("full_gradient_response_identified") is True
    )
    full_physical_eligible = bool(
        isinstance(gradient_response, dict)
        and gradient_response.get("full_physical_eligible") is True
    )
    elastic_energy_loss = prepared.plan_payload.get("closure", {}).get(
        "elastic_energy_loss"
    )
    mc_censoring = prepared.plan_payload.get("closure", {}).get(
        "mc_rate_censoring"
    )
    mc_physics = (
        mc_censoring.get("physics_qualification")
        if isinstance(mc_censoring, dict)
        else None
    )
    mc_censoring_status = (
        "passed"
        if isinstance(mc_physics, dict) and mc_physics.get("passed") is True
        else "failed" if isinstance(mc_physics, dict) else "not_applicable"
    )
    elastic_statistics = (
        elastic_energy_loss.get("statistical_evidence")
        if isinstance(elastic_energy_loss, dict)
        else None
    )
    elastic_statistics_status = _elastic_statistics_status(
        plan,
        elastic_statistics,
    )
    upstream_mc_quality = prepared.plan_payload.get("bundle", {}).get(
        "independent_quality_audit"
    )
    upstream_mc_quality_status = (
        (
            "passed"
            if isinstance(upstream_mc_quality, dict)
            and upstream_mc_quality.get("passed") is True
            else "failed"
        )
        if plan.mapping.bundle.expected_source == "monte_carlo"
        else "not_applicable"
    )
    mc_statistical_status = (
        (
            "passed"
            if upstream_mc_quality_status == "passed"
            and mc_censoring_status in {"passed", "not_applicable"}
            and elastic_statistics_status in {"passed", "not_applicable"}
            else "failed"
        )
        if plan.mapping.bundle.expected_source == "monte_carlo"
        else "not_applicable"
    )
    physics_quality_status = (
        "passed"
        if numerical_quality_status == "passed"
        and mc_statistical_status in {"passed", "not_applicable"}
        else "failed"
    )
    return RunQuality(
        conservation_status=conservation_status,
        builtin_reference_conservation_status=(
            builtin_reference_conservation_status
        ),
        transport_status=transport_status,
        closure_support_status=closure_support_status,
        function_eedf_status=function_status,
        numerical_status=numerical_quality_status,
        physics_status=physics_quality_status,
        mc_censoring_status=mc_censoring_status,
        mc_statistical_status=mc_statistical_status,
        upstream_mc_quality_status=upstream_mc_quality_status,
        elastic_statistics_status=elastic_statistics_status,
        accepted=physics_quality_status == "passed",
        gradient_response=gradient_response,
        full_gradient_response_identified=full_gradient_response_identified,
        full_physical_eligible=full_physical_eligible,
        transport_support=transport_support,
        closure_support_qualification=closure_support_qualification,
        mc_censoring=mc_censoring,
        elastic_statistics=elastic_statistics,
    )


def _elastic_statistics_status(
    plan: GecCcpPlan,
    elastic_statistics: dict[str, Any] | None,
) -> str:
    external_mc_elastic = (
        plan.mapping.bundle.expected_source == "monte_carlo"
        and _uses_external_elastic_energy_loss(plan.mapping.closure)
    )
    if not external_mc_elastic:
        return "not_applicable"
    if (
        isinstance(elastic_statistics, dict)
        and elastic_statistics.get("uncertainty_available_at_every_anchor")
        is True
        and elastic_statistics.get("relative_standard_error_qualified") is True
    ):
        return "passed"
    return "failed"


__all__ = ["evaluate_run_quality", "run_postsolve_audits"]
