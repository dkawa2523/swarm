"""Prepare, run, and accept canonical external-Swarm GEC-CCP models.

Both production targets import a Function EEDF, reduced mobility, and elastic
energy loss.  COMSOL retains particle diffusion, electron-energy transport,
RF fields, species, and boundaries under its restricted local-mean-energy
formulation.  The built-in Druyvesteyn solve is an optional comparison and is
never repeated by a normal physical-target run.
"""

from __future__ import annotations

import json as _json

from swarm_workflow.comsol.models.gec_ccp.contracts import (
    GecCcpPlan,
    GecCcpQualityError as GecCcpQualityError,
    GecCcpRunSummary,
    GecCcpWorkflowError as GecCcpWorkflowError,
)
from swarm_workflow.comsol.models.gec_ccp.execution import execute_gec_ccp_run as execute_gec_ccp_run
from swarm_workflow.comsol.models.gec_ccp.prepare import prepare_gec_ccp_run as prepare_gec_ccp_run

__all__ = (
    "GecCcpPlan",
    "GecCcpQualityError",
    "GecCcpRunSummary",
    "GecCcpWorkflowError",
    "execute_gec_ccp_run",
    "format_gec_ccp_plan",
    "format_gec_ccp_summary",
    "prepare_gec_ccp_run",
)


def format_gec_ccp_plan(plan: GecCcpPlan) -> str:
    closure_description = (
        "COMSOL-solved mean energy; electron_transport="
        f"{plan.mapping.closure.electron_transport}; reaction_model="
        f"{plan.mapping.closure.reaction_model}; elastic_energy_loss_model="
        f"{plan.mapping.closure.elastic_energy_loss_model}"
    )
    return "\n".join(
        [
            "GEC CCP COMSOL dry-run ready",
            f"model: {plan.mapping.model.input_mph}",
            f"detected EEDF: {plan.contract.original_eedf}",
            f"bundle: {plan.mapping.bundle.path}",
            f"result role: {plan.mapping.results.role}",
            f"closure: {closure_description}",
            f"plan: {plan.plan_json}",
            f"generated Java: {plan.output_directory}",
            "COMSOL command: not executed",
        ]
    )


def format_gec_ccp_summary(summary: GecCcpRunSummary) -> str:
    status = _json.loads(summary.status_json.read_text(encoding="utf-8"))
    return "\n".join(
        [
            "GEC CCP COMSOL workflow completed",
            f"result role: {status.get('result_role')}",
            "physical target accepted: "
            f"{status.get('physical_target_accepted')}",
            f"status: {summary.status_json}",
            f"results: {summary.plan.output_directory}",
        ]
    )
