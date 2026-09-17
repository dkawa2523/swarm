"""Compose GEC-CCP closure-support and saved-transport run audits."""

from __future__ import annotations

import hashlib
from pathlib import Path
from zipfile import ZipFile

from swarm_workflow._io import write_json
from swarm_workflow.comsol.models.gec_ccp.closure import _uses_transport
from swarm_workflow.comsol.models.gec_ccp.audits.conservation import _read_comsol_numeric_csv
from swarm_workflow.comsol.models.gec_ccp.contracts import GecCcpPlan, GecCcpWorkflowError
from swarm_workflow.comsol.models.gec_ccp.validation.joint_consistency import _active_closure_mean_energy_support
from swarm_workflow.comsol.models.gec_ccp.audits.saved_model import _audit_transport_binding
from swarm_workflow.comsol.models.gec_ccp.audits.transport_support import _operating_mean_energy_support_audit


def audit_gec_ccp_closure_support_run(plan: GecCcpPlan) -> Path:
    """Audit converged support for external rates without external transport."""

    mapping = plan.mapping
    domain_path = (
        mapping.results.output_directory
        / "swarm_tables"
        / "domain_phase_closure.csv"
    )
    headers, values = _read_comsol_numeric_csv(domain_path)
    support = _operating_mean_energy_support_audit(
        mapping,
        headers,
        values,
        _active_closure_mean_energy_support(mapping),
    )
    audit_path = plan.output_directory / "closure_support_audit.json"
    write_json(
        audit_path,
        {
            "stage": "audit-gec-ccp-closure-support",
            "status": "passed" if support.get("passed") else "failed",
            "support": support,
        },
    )
    if not support.get("passed"):
        raise GecCcpWorkflowError(
            "converged COMSOL result left the qualified external-table "
            f"support; see {audit_path}"
        )
    return audit_path


def audit_gec_ccp_transport_run(plan: GecCcpPlan) -> Path:
    """Audit saved bindings and converged values for external transport."""

    mapping = plan.mapping
    if not _uses_transport(mapping.closure):
        raise GecCcpWorkflowError(
            "transport audit requires an external electron-transport closure"
        )
    try:
        with ZipFile(mapping.model.output_mph) as archive:
            model_xml = archive.read("dmodel.xml").decode(
                "utf-8", errors="replace"
            )
    except (OSError, KeyError) as exc:
        raise GecCcpWorkflowError(
            f"cannot audit saved external-transport MPH: {mapping.model.output_mph}"
        ) from exc
    transport = _audit_transport_binding(mapping, model_xml)
    audit_path = plan.output_directory / "transport_audit.json"
    write_json(
        audit_path,
        {
            "stage": "audit-gec-ccp-transport",
            "status": "passed" if transport["passed"] else "failed",
            "model": {
                "path": str(mapping.model.output_mph),
                "sha256": hashlib.sha256(
                    mapping.model.output_mph.read_bytes()
                ).hexdigest(),
            },
            "closure": mapping.closure.electron_transport,
            "transport": transport,
        },
    )
    if not transport["passed"]:
        raise GecCcpWorkflowError(
            "converged COMSOL result failed the external-transport audit; "
            f"see {audit_path}"
        )
    return audit_path
