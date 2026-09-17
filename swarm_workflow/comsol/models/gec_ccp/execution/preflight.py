"""Prepare and freeze a GEC-CCP execution before COMSOL is launched."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from swarm_workflow._io import write_json

from ..contracts import GecCcpWorkflowError
from ..mapping import load_gec_ccp_mapping
from ..prepare import prepare_gec_ccp_run
from .context import PreparedExecution
from .inputs import validate_gec_plan_inputs
from .status import status_contract_metadata


def prepare_execution(
    mapping_path: str | Path,
    *,
    bundle_path: str | Path | None,
) -> PreparedExecution:
    """Validate, materialize, freeze, and clean one execution plan."""

    mapping = load_gec_ccp_mapping(mapping_path, bundle_path=bundle_path)
    status_path = mapping.results.output_directory / "gec_ccp_run_status.json"
    status_path.parent.mkdir(parents=True, exist_ok=True)
    status_path.unlink(missing_ok=True)
    write_json(
        status_path,
        {
            "stage": "preflight-gec-ccp",
            "status": "running",
            "solve_status": "not_started",
            "quality_status": "not_evaluated",
            **status_contract_metadata(mapping),
        },
    )
    try:
        plan = prepare_gec_ccp_run(
            mapping_path,
            bundle_path=bundle_path,
            write_java=True,
        )
    except GecCcpWorkflowError as exc:
        write_json(
            status_path,
            {
                "stage": "preflight-gec-ccp",
                "status": "blocked",
                "blocker_code": "preflight_validation_failed",
                "reason": str(exc),
                "solve_status": "not_started",
                "quality_status": "not_evaluated",
                **status_contract_metadata(mapping),
            },
        )
        raise
    plan_payload = json.loads(plan.plan_json.read_text(encoding="utf-8"))
    plan_gradient_response = plan_payload.get("closure", {}).get(
        "gradient_response"
    )
    plan_sha256 = hashlib.sha256(plan.plan_json.read_bytes()).hexdigest()
    try:
        input_provenance = validate_gec_plan_inputs(plan)
    except GecCcpWorkflowError as exc:
        write_json(
            status_path,
            {
                "stage": "run-gec-ccp",
                "status": "blocked",
                "blocker_code": "plan_input_provenance_mismatch",
                "reason": str(exc),
                "solve_status": "not_started",
                "quality_status": "not_evaluated",
                "plan": str(plan.plan_json),
                "plan_sha256": plan_sha256,
                **status_contract_metadata(
                    plan.mapping,
                    gradient_response=plan_gradient_response,
                ),
            },
        )
        raise
    stale_outputs = [
        *plan.expected_result_files,
        plan.output_directory / "conservation_audit.json",
        plan.output_directory / "transport_audit.json",
        plan.output_directory / "function_eedf_audit.json",
        plan.output_directory / "function_eedf_rate_audit.csv",
    ]
    if plan.native_eedf_audit is not None:
        stale_outputs.extend(
            (
                plan.native_eedf_audit.values_path,
                plan.native_eedf_audit.contract_path,
                plan.native_eedf_audit.values_path.parent
                / "comsol_eedf_audit.json",
            )
        )
    for stale_output in stale_outputs:
        stale_output.unlink(missing_ok=True)
    return PreparedExecution(
        plan=plan,
        status_path=status_path,
        plan_payload=plan_payload,
        gradient_response=plan_gradient_response,
        plan_sha256=plan_sha256,
        input_provenance=input_provenance,
    )


__all__ = ["prepare_execution"]
