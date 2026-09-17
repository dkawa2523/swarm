"""Validate GEC-CCP bundle identity, provenance, and artifact hashes."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import re
from typing import Any

import electron_swarm.solvers.monte_carlo.evidence as _mc_evidence
from electron_swarm.solvers.monte_carlo.source_identity import (
    monte_carlo_source_sha256,
)

from swarm_workflow.campaign.store import (
    WorkflowSchemaError,
    canonical_mc_sampling_plan_json,
)
from swarm_workflow.quality.policy import (
    parse_quality_thresholds,
    quality_thresholds_payload,
)
from swarm_workflow.quality.solver import (
    PropagatorTargetRequirement,
    SolverQualificationError,
    validate_copied_qualification,
    validate_propagator_core_qualification,
    validate_propagator_target_qualification,
)
from swarm_workflow.tables.contracts import (
    MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA,
)

from ..contracts import GecCcpMapping, GecCcpWorkflowError
from .bundle_context import (
    GEC_OPTIONAL_PROVENANCE_HASH_KEYS,
    GEC_PROVENANCE_HASH_KEYS,
    ManifestEvidence,
)


_GEC_CCP_PROPAGATOR_TARGET_REQUIREMENT = PropagatorTargetRequirement(
    target="argon_gec_ccp_restricted_local_mean_energy",
    fields_Td=(3000.0, 3500.0, 4000.0),
    operating_range_bracket_Td=(3000.0, 3500.0),
    table_support_cap_Td=4000.0,
    mixture_id=0,
    medium_grid=(300, 48),
    fine_grid=(600, 72),
    max_scalar_refinement_limit=0.01,
    max_eedf_weighted_L1_limit=0.02,
)


def validate_manifest(
    mapping: GecCcpMapping, *, uses_transport: bool
) -> ManifestEvidence:
    (
        manifest_path,
        manifest,
        source,
        solver_qualification,
        target_qualification,
        source_policy,
        expected_mc_profile,
    ) = _validate_manifest_identity_and_policy(mapping, uses_transport=uses_transport)
    (
        hashes,
        mc_sampling_plan,
        mc_estimator_schema,
        mc_eedf_estimator_schema,
        mc_tail_estimator_schema,
        quality_thresholds,
    ) = _validate_manifest_provenance(mapping, manifest, source=source)
    tables, artifact_digests = _verify_manifest_artifacts(mapping, manifest)
    return ManifestEvidence(
        path=manifest_path,
        manifest=manifest,
        source=source,
        solver_qualification=solver_qualification,
        target_qualification=target_qualification,
        source_policy=source_policy,
        expected_mc_profile=expected_mc_profile,
        hashes=hashes,
        mc_sampling_plan=mc_sampling_plan,
        mc_estimator_schema=mc_estimator_schema,
        mc_eedf_estimator_schema=mc_eedf_estimator_schema,
        mc_tail_estimator_schema=mc_tail_estimator_schema,
        quality_thresholds=quality_thresholds,
        tables=tables,
        artifact_digests=artifact_digests,
    )


def _validate_manifest_identity_and_policy(
    mapping: GecCcpMapping, *, uses_transport: bool
) -> tuple[
    Path,
    dict[str, Any],
    str,
    dict[str, Any] | None,
    dict[str, Any] | None,
    dict[str, Any],
    str | None,
]:
    manifest_path = mapping.bundle.path / "manifest.json"
    if not manifest_path.exists():
        raise GecCcpWorkflowError(f"missing COMSOL bundle manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "ok":
        raise GecCcpWorkflowError("COMSOL bundle manifest status is not ok")
    source = manifest.get("source")
    if source != mapping.bundle.expected_source:
        raise GecCcpWorkflowError(
            "COMSOL bundle source mismatch: expected "
            f"{mapping.bundle.expected_source}, found {source}"
        )
    solver_qualification, target_qualification = _propagator_qualifications(
        mapping, manifest, source=str(source)
    )
    source_policy = manifest.get("source_policy")
    if not isinstance(source_policy, dict):
        raise GecCcpWorkflowError("COMSOL bundle lacks source_policy metadata")
    expected_mc_profile = mapping.bundle.expected_mc_qualification_profile
    if source == "monte_carlo" and expected_mc_profile is not None:
        actual_mc_profile = source_policy.get("qualification_profile")
        if actual_mc_profile != expected_mc_profile:
            raise GecCcpWorkflowError(
                "Monte Carlo bundle qualification profile mismatch: expected "
                f"{expected_mc_profile}, found {actual_mc_profile}"
            )
        if expected_mc_profile == MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA:
            expected_solver_physics = {
                "population_model": "weighted_branching",
                "angular_model": "isotropic",
                "angular_moment_source": "isotropic_closure",
                "exact_dcs_based": False,
                "ordinary_integral_xs_closure": True,
                "ionization_source_model": "equal",
                "ionization_source_treatment": "equal",
                "electron_electron_treatment": "none",
                "electron_electron_transport_stale": False,
                "magnetic_field_treatment": "none",
                "tail_refinement_treatment": "approximate",
                "high_energy_extrapolation": "zero",
            }
            if source_policy.get("solver_physics") != expected_solver_physics:
                raise GecCcpWorkflowError(
                    "Monte Carlo restricted-LMEA solver-physics assumptions "
                    "are missing or inconsistent"
                )
    _validate_source_policy(
        mapping,
        source_policy,
        source=str(source),
        uses_transport=uses_transport,
    )
    return (
        manifest_path,
        manifest,
        str(source),
        solver_qualification,
        target_qualification,
        source_policy,
        expected_mc_profile,
    )


def _propagator_qualifications(
    mapping: GecCcpMapping,
    manifest: dict[str, Any],
    *,
    source: str,
) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
    if source != "propagator":
        if manifest.get("solver_qualification") is not None:
            raise GecCcpWorkflowError(
                "solver qualification evidence is only valid for a propagator bundle"
            )
        if manifest.get("target_qualification") is not None:
            raise GecCcpWorkflowError(
                "target qualification evidence is only valid for a propagator bundle"
            )
        return None, None
    entry = manifest.get("solver_qualification")
    if not isinstance(entry, dict):
        raise GecCcpWorkflowError(
            "Propagator COMSOL bundle lacks bound numerical qualification evidence"
        )
    try:
        qualification_path = validate_copied_qualification(mapping.bundle.path, entry)
        validated = validate_propagator_core_qualification(qualification_path)
    except SolverQualificationError as exc:
        raise GecCcpWorkflowError(str(exc)) from exc
    expected_entry = validated.manifest_entry(file=qualification_path.name)
    if entry != expected_entry:
        raise GecCcpWorkflowError(
            "Propagator qualification manifest metadata is inconsistent"
        )
    target_entry = manifest.get("target_qualification")
    if not isinstance(target_entry, dict):
        raise GecCcpWorkflowError(
            "Propagator GEC-CCP bundle lacks bound target-range qualification evidence"
        )
    try:
        target_path = validate_copied_qualification(mapping.bundle.path, target_entry)
        validated_target = validate_propagator_target_qualification(
            target_path,
            core_qualification=validated,
            provenance_hashes=manifest.get("hashes", {}),
            requirement=_GEC_CCP_PROPAGATOR_TARGET_REQUIREMENT,
        )
    except SolverQualificationError as exc:
        raise GecCcpWorkflowError(str(exc)) from exc
    expected_target = validated_target.manifest_entry(file=target_path.name)
    if target_entry != expected_target:
        raise GecCcpWorkflowError(
            "Propagator target qualification manifest metadata is inconsistent"
        )
    return expected_entry, expected_target


def _validate_source_policy(
    mapping: GecCcpMapping,
    source_policy: dict[str, Any],
    *,
    source: str,
    uses_transport: bool,
) -> None:
    field_type = source_policy.get("field_type")
    if field_type != mapping.bundle.expected_field_type:
        raise GecCcpWorkflowError(
            "COMSOL bundle field_type mismatch: expected "
            f"{mapping.bundle.expected_field_type}, found {field_type}"
        )
    if mapping.closure.source_field == "steady_dc" and field_type != "dc":
        raise GecCcpWorkflowError(
            "closure.source_field=steady_dc requires bundle field_type=dc"
        )
    rf_frequency = source_policy.get("rf_frequency_Hz")
    if field_type == "dc" and rf_frequency not in {None, ""}:
        raise GecCcpWorkflowError(
            "DC bundle source_policy must not specify rf_frequency_Hz"
        )
    transport_definition = source_policy.get("transport_definition")
    if transport_definition != mapping.bundle.expected_transport_definition:
        raise GecCcpWorkflowError(
            "COMSOL bundle transport_definition mismatch: expected "
            f"{mapping.bundle.expected_transport_definition}, found "
            f"{transport_definition}"
        )
    if uses_transport and "flux" not in str(transport_definition).lower():
        raise GecCcpWorkflowError(
            "external electron transport requires a flux transport definition"
        )
    if source_policy.get("source") not in {None, source}:
        raise GecCcpWorkflowError(
            "bundle source_policy.source disagrees with manifest source"
        )
    if source_policy.get("postprocess") != "none":
        raise GecCcpWorkflowError(
            "GEC-CCP requires source_policy.postprocess=none; rebuild a "
            "qualified pure-solver bundle without prior substitution or "
            "coefficient/EEDF postprocessing"
        )


def _validate_manifest_provenance(
    mapping: GecCcpMapping,
    manifest: dict[str, Any],
    *,
    source: str,
) -> tuple[
    dict[str, Any],
    list[dict[str, Any]] | None,
    str | None,
    str | None,
    str | None,
    dict[str, Any],
]:
    hashes = manifest.get("hashes")
    required_hashes = set(GEC_PROVENANCE_HASH_KEYS)
    allowed_hashes = required_hashes | set(GEC_OPTIONAL_PROVENANCE_HASH_KEYS)
    if (
        not isinstance(hashes, dict)
        or not required_hashes.issubset(hashes)
        or not set(hashes).issubset(allowed_hashes)
        or any(
            re.fullmatch(r"[0-9a-f]{64}", str(hashes.get(name, ""))) is None
            for name in GEC_PROVENANCE_HASH_KEYS
        )
        or (
            "source_database_sha256" in hashes
            and re.fullmatch(r"[0-9a-f]{64}", str(hashes["source_database_sha256"]))
            is None
        )
    ):
        raise GecCcpWorkflowError(
            "bundle manifest requires the canonical SHA-256 chain for "
            "workflow config, base config, and cross sections"
        )
    sampling_plan = None
    estimator_schema = None
    eedf_estimator_schema = None
    tail_schema = None
    if source == "monte_carlo":
        (
            sampling_plan,
            estimator_schema,
            eedf_estimator_schema,
            tail_schema,
        ) = _mc_provenance(mapping, hashes)
    elif any(
        name in hashes
        for name in GEC_OPTIONAL_PROVENANCE_HASH_KEYS
        if name.startswith("mc_")
    ):
        raise GecCcpWorkflowError(
            "non-Monte-Carlo bundle must not contain MC sampling provenance"
        )
    raw_thresholds = manifest.get("quality_thresholds")
    if not isinstance(raw_thresholds, dict):
        raise GecCcpWorkflowError("bundle manifest lacks quality_thresholds provenance")
    try:
        thresholds = quality_thresholds_payload(
            parse_quality_thresholds(raw_thresholds)
        )
    except (TypeError, ValueError) as exc:
        raise GecCcpWorkflowError(
            "bundle manifest has invalid quality_thresholds provenance"
        ) from exc
    if raw_thresholds != thresholds:
        raise GecCcpWorkflowError(
            "bundle manifest quality_thresholds is incomplete or noncanonical"
        )
    return (
        hashes,
        sampling_plan,
        estimator_schema,
        eedf_estimator_schema,
        tail_schema,
        thresholds,
    )


def _mc_provenance(
    mapping: GecCcpMapping, hashes: dict[str, Any]
) -> tuple[list[dict[str, Any]], str, str, str | None]:
    sampling_json = hashes.get("mc_sampling_plan_json")
    sampling_sha256 = hashes.get("mc_sampling_plan_sha256")
    estimator_version = hashes.get("mc_transport_estimator_schema_version")
    eedf_estimator_version = hashes.get("mc_eedf_estimator_schema_version")
    tail_version = hashes.get("mc_tail_estimator_schema_version")
    seed_derivation = hashes.get("mc_seed_derivation_schema_version")
    solver_source = hashes.get("mc_solver_source_sha256")
    try:
        decoded = json.loads(str(sampling_json))
        canonical = canonical_mc_sampling_plan_json(decoded)
        decoded = json.loads(canonical)
    except (TypeError, ValueError, WorkflowSchemaError) as exc:
        raise GecCcpWorkflowError(
            "Monte Carlo bundle has invalid mc_sampling_plan_json"
        ) from exc
    for row in decoded:
        field = float(row["e_over_n_Td"])
        particles = int(row["particles"])
        warmup = int(row["warmup_collisions"])
        collisions = int(row["max_collisions"])
        replicas = int(row["replicas"])
        if (
            not math.isfinite(field)
            or field <= 0.0
            or particles < 2
            or warmup < 0
            or collisions <= 0
            or replicas < 2
        ):
            raise GecCcpWorkflowError(
                "Monte Carlo sampling-plan values are outside the qualification domain"
            )
    expected_estimator = mapping.bundle.expected_mc_transport_estimator_schema
    expected_tail = mapping.bundle.expected_mc_tail_estimator_schema
    if (
        sampling_json != canonical
        or hashlib.sha256(canonical.encode("utf-8")).hexdigest() != sampling_sha256
        or estimator_version
        not in {
            _mc_evidence.DIRECT_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
            _mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
        }
        or eedf_estimator_version != _mc_evidence.MC_EEDF_ESTIMATOR_SCHEMA_VERSION
        or (
            seed_derivation is not None
            and seed_derivation != _mc_evidence.MC_SEED_DERIVATION_SCHEMA_VERSION
        )
        or solver_source != monte_carlo_source_sha256()
        or (expected_estimator is not None and estimator_version != expected_estimator)
        or tail_version
        not in {
            None,
            _mc_evidence.WEIGHTED_MC_TAIL_ESTIMATOR_SCHEMA_VERSION,
        }
        or (expected_tail is not None and tail_version != expected_tail)
    ):
        raise GecCcpWorkflowError(
            "Monte Carlo sampling-plan provenance or estimator version is inconsistent"
        )
    return (
        decoded,
        str(estimator_version),
        str(eedf_estimator_version),
        None if tail_version is None else str(tail_version),
    )


def _verify_manifest_artifacts(
    mapping: GecCcpMapping, manifest: dict[str, Any]
) -> tuple[dict[str, dict[str, Any]], dict[str, str]]:
    tables = manifest.get("tables")
    if not isinstance(tables, dict) or not tables:
        raise GecCcpWorkflowError("bundle manifest.tables must be a mapping")
    artifact_digests: dict[str, str] = {}
    bundle_root = mapping.bundle.path.resolve()
    for name, metadata in tables.items():
        relative_name = PurePosixPath(name) if isinstance(name, str) else None
        if (
            not isinstance(name, str)
            or not name
            or "\\" in name
            or relative_name is None
            or relative_name.is_absolute()
            or relative_name.as_posix() != name
            or any(part in {"", ".", ".."} for part in relative_name.parts)
            or not isinstance(metadata, dict)
        ):
            raise GecCcpWorkflowError(
                "bundle manifest contains an invalid table artifact entry"
            )
        artifact_path = (bundle_root / Path(*relative_name.parts)).resolve()
        try:
            artifact_path.relative_to(bundle_root)
        except ValueError as exc:
            raise GecCcpWorkflowError(
                "bundle manifest contains an invalid table artifact entry"
            ) from exc
        expected_digest = metadata.get("sha256")
        role = metadata.get("artifact_role")
        if not artifact_path.is_file():
            raise GecCcpWorkflowError(f"bundle artifact does not exist: {name}")
        if (
            re.fullmatch(r"[0-9a-f]{64}", str(expected_digest or "")) is None
            or not isinstance(role, str)
            or not role.strip()
        ):
            raise GecCcpWorkflowError(
                f"bundle artifact lacks canonical SHA-256 or role: {name}"
            )
        actual_digest = hashlib.sha256(artifact_path.read_bytes()).hexdigest()
        if actual_digest != expected_digest:
            raise GecCcpWorkflowError(f"bundle artifact SHA-256 mismatch: {name}")
        artifact_digests[name] = actual_digest
    return tables, artifact_digests
