"""Typed state passed between the GEC-CCP execution pipeline stages."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from swarm_workflow.comsol.runtime import ComsolExecutionSummary

from ..contracts import GecCcpPlan, GecCcpWorkflowError


@dataclass(frozen=True)
class PreparedExecution:
    """Immutable preflight result consumed by the COMSOL runner."""

    plan: GecCcpPlan
    status_path: Path
    plan_payload: dict[str, Any]
    gradient_response: dict[str, Any] | None
    plan_sha256: str
    input_provenance: dict[str, Any]


@dataclass(frozen=True)
class SolverExecution:
    """COMSOL operations and the optional native EEDF audit outcome."""

    executions: tuple[ComsolExecutionSummary, ...]
    native_eedf_audit_path: Path | None
    native_eedf_status: str


@dataclass(frozen=True)
class PostsolveAudits:
    """Paths and retained exceptions from the independent run audits."""

    conservation_path: Path
    transport_path: Path | None
    closure_support_path: Path | None
    function_eedf_path: Path | None
    transport_error: GecCcpWorkflowError | None
    closure_support_error: GecCcpWorkflowError | None
    function_eedf_error: GecCcpWorkflowError | None


@dataclass(frozen=True)
class RunQuality:
    """All quality decisions needed to serialize and enforce acceptance."""

    conservation_status: str
    builtin_reference_conservation_status: str
    transport_status: str
    closure_support_status: str
    function_eedf_status: str
    numerical_status: str
    physics_status: str
    mc_censoring_status: str
    mc_statistical_status: str
    upstream_mc_quality_status: str
    elastic_statistics_status: str
    accepted: bool
    gradient_response: dict[str, Any] | None
    full_gradient_response_identified: bool
    full_physical_eligible: bool
    transport_support: dict[str, Any] | None
    closure_support_qualification: dict[str, Any] | None
    mc_censoring: dict[str, Any] | None
    elastic_statistics: dict[str, Any] | None


__all__ = [
    "PostsolveAudits",
    "PreparedExecution",
    "RunQuality",
    "SolverExecution",
]
