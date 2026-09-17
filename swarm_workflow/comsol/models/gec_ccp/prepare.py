"""Prepare a validated, provenance-bound GEC-CCP COMSOL run plan."""

from __future__ import annotations

from pathlib import Path

import swarm_workflow.comsol.models.gec_ccp.prepare_stages as prepare_stages
from swarm_workflow.comsol.models.gec_ccp.contracts import (
    GecCcpPlan,
)
from swarm_workflow.comsol.models.gec_ccp.mapping import load_gec_ccp_mapping
from swarm_workflow.comsol.models.gec_ccp.mph import (
    _validate_contract,
    inspect_gec_ccp_mph,
)
from swarm_workflow.comsol.models.gec_ccp.planning.evidence import (
    build_closure_evidence,
)


def prepare_gec_ccp_run(
    mapping_path: str | Path,
    *,
    bundle_path: str | Path | None = None,
    write_java: bool = True,
) -> GecCcpPlan:
    mapping = load_gec_ccp_mapping(mapping_path, bundle_path=bundle_path)
    contract = inspect_gec_ccp_mph(mapping.model.input_mph)
    _validate_contract(mapping, contract)
    evidence = build_closure_evidence(mapping)
    layout = prepare_stages.prepare_layout(
        mapping,
        contract,
        write_java=write_java,
    )
    prepare_stages.materialize_java(
        layout,
        preintegrated_rate_rows=evidence.preintegrated_rate_rows,
        write_java=write_java,
    )
    prepare_stages.write_plan_manifest(layout, evidence, write_java=write_java)
    return layout.plan
