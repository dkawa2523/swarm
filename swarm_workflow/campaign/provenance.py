"""Canonical workflow identity and safe campaign-resume decisions."""

from __future__ import annotations

from dataclasses import asdict, replace
from hashlib import sha256
import json

from ..quality.monte_carlo.policy import campaign_provenance
from ..quality.policy import (
    QUALITY_THRESHOLDS_METADATA_KEY,
    quality_thresholds_json,
)
from . import config as _workflow_config
from .quality import resolve_quality_thresholds, source_quality_thresholds_json
from .store import (
    PROPAGATOR_SOLVER_SOURCE_METADATA_KEY,
    WorkflowStore,
    canonical_mc_sampling_plan_json,
)


MC_TRANSPORT_ESTIMATOR_METADATA_KEY = "mc_transport_estimator_schema_version"
MC_EEDF_ESTIMATOR_METADATA_KEY = "mc_eedf_estimator_schema_version"
MC_TAIL_ESTIMATOR_METADATA_KEY = "mc_tail_estimator_schema_version"
MC_SEED_DERIVATION_METADATA_KEY = "mc_seed_derivation_schema_version"
MC_SOLVER_SOURCE_METADATA_KEY = "mc_solver_source_sha256"
DETERMINISTIC_EXECUTION_METADATA_KEY = "deterministic_execution_json"


def mc_sampling_plan_payload(
    workflow: _workflow_config.WorkflowConfig,
) -> list[dict[str, int | float | str | None]]:
    """Return the complete canonical MC sampling plan represented by a workflow."""

    return [asdict(entry) for entry in workflow.mc_sampling_plan]


def mc_sampling_plan_json(workflow: _workflow_config.WorkflowConfig) -> str:
    """Serialize the workflow MC plan with the database provenance contract."""

    return canonical_mc_sampling_plan_json(mc_sampling_plan_payload(workflow))


def workflow_config_sha256(
    workflow: _workflow_config.WorkflowConfig,
    *,
    quality_payload: dict[str, object] | None = None,
    mc_transport_estimator_schema_version: str | None = None,
    mc_eedf_estimator_schema_version: str | None = None,
    mc_tail_estimator_schema_version: str | None = None,
    mc_seed_derivation_schema_version: str | None = None,
    mc_solver_source_sha256: str | None = None,
    propagator_solver_source_sha256: str | None = None,
    deterministic_execution_payload: dict[str, object] | None = None,
) -> str:
    """Hash validated workflow semantics with every filesystem path resolved."""

    payload: dict[str, object] = {
        "base_config_path": str(workflow.base_config_path),
        "database_path": str(workflow.database_path),
        "e_over_n_Td": list(workflow.e_over_n_Td),
        "mixtures": [
            {
                "mixture_id": mixture.mixture_id,
                "fractions": dict(sorted(mixture.fractions.items())),
            }
            for mixture in workflow.mixtures
        ],
        "quality": (
            asdict(workflow.quality) if quality_payload is None else quality_payload
        ),
    }
    if workflow.mc_enabled:
        payload["mc"] = {
            "e_over_n_Td": (
                list(workflow.mc_e_over_n_Td)
                if workflow.mc_e_over_n_Td is not None
                else None
            ),
            "replicas": workflow.mc_replicas,
            "base_seed": workflow.mc_base_seed,
            "workers": workflow.mc_workers,
            "sampling_plan": mc_sampling_plan_payload(workflow),
        }
    if deterministic_execution_payload is not None:
        payload["execution"] = {"deterministic": deterministic_execution_payload}
    if workflow.mean_energy_support is not None:
        payload["mean_energy_support"] = asdict(workflow.mean_energy_support)
    if mc_transport_estimator_schema_version is not None:
        payload[MC_TRANSPORT_ESTIMATOR_METADATA_KEY] = (
            mc_transport_estimator_schema_version
        )
    if mc_eedf_estimator_schema_version is not None:
        payload[MC_EEDF_ESTIMATOR_METADATA_KEY] = mc_eedf_estimator_schema_version
    if workflow.mc_convergence is not None:
        payload["mc_campaign"] = campaign_provenance(
            workflow.mc_convergence,
            workflow.mc_previous_decision,
        )
    if mc_tail_estimator_schema_version is not None:
        payload[MC_TAIL_ESTIMATOR_METADATA_KEY] = mc_tail_estimator_schema_version
    if mc_seed_derivation_schema_version is not None:
        payload[MC_SEED_DERIVATION_METADATA_KEY] = mc_seed_derivation_schema_version
    if mc_solver_source_sha256 is not None:
        payload[MC_SOLVER_SOURCE_METADATA_KEY] = mc_solver_source_sha256
    if propagator_solver_source_sha256 is not None:
        payload[PROPAGATOR_SOLVER_SOURCE_METADATA_KEY] = (
            propagator_solver_source_sha256
        )
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return sha256(encoded).hexdigest()


def resume_provenance(
    store: WorkflowStore,
    workflow: _workflow_config.WorkflowConfig,
    current: dict[str, str],
    *,
    deterministic_execution_payload: dict[str, object] | None = None,
) -> tuple[dict[str, str], bool]:
    """Accept an existing database only when the workflow change is quality-only."""

    metadata = store.metadata()
    stored_hash = metadata.get("workflow_config_sha256")
    current_hash = current["workflow_config_sha256"]
    if stored_hash is None or stored_hash == current_hash:
        return current, False

    source_quality = resolve_quality_thresholds(store.connection)
    source_quality_json = source_quality_thresholds_json(store.connection)
    source_workflow = replace(workflow, quality=source_quality)
    source_quality_payload = json.loads(source_quality_json)
    compatible_hash = workflow_config_sha256(
        source_workflow,
        quality_payload=source_quality_payload,
        mc_transport_estimator_schema_version=current.get(
            MC_TRANSPORT_ESTIMATOR_METADATA_KEY
        ),
        mc_eedf_estimator_schema_version=current.get(MC_EEDF_ESTIMATOR_METADATA_KEY),
        mc_tail_estimator_schema_version=current.get(MC_TAIL_ESTIMATOR_METADATA_KEY),
        mc_seed_derivation_schema_version=current.get(MC_SEED_DERIVATION_METADATA_KEY),
        mc_solver_source_sha256=current.get(MC_SOLVER_SOURCE_METADATA_KEY),
        propagator_solver_source_sha256=current.get(
            PROPAGATOR_SOLVER_SOURCE_METADATA_KEY
        ),
        deterministic_execution_payload=deterministic_execution_payload,
    )
    if stored_hash != compatible_hash:
        return current, False

    evaluation_quality_json = quality_thresholds_json(workflow.quality)
    return (
        {
            **current,
            "workflow_config_sha256": stored_hash,
            QUALITY_THRESHOLDS_METADATA_KEY: source_quality_json,
        },
        source_quality_json != evaluation_quality_json,
    )
