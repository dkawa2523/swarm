"""Execution-time verification of every artifact frozen by a GEC plan."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath
import re
from typing import Any

from swarm_workflow.comsol.input.bundle_selection import validate_bundle_selection
from swarm_workflow.selection import ClosureSelectionError

from ..contracts import GecCcpPlan, GecCcpWorkflowError
from ..validation.bundle_guards import _validate_low_energy_guard_bundle


def validate_gec_plan_inputs(plan: GecCcpPlan) -> dict[str, Any]:
    """Fail if mapping, MPH, bundle, guard, or generated Java has changed."""

    payload = json.loads(plan.plan_json.read_text(encoding="utf-8"))
    try:
        selection = validate_bundle_selection(
            plan.mapping.bundle.path,
            required=plan.mapping.bundle.require_solver_selection,
        )
    except ClosureSelectionError as exc:
        raise GecCcpWorkflowError(str(exc)) from exc
    if selection != payload.get("bundle", {}).get("solver_selection"):
        raise GecCcpWorkflowError(
            "solver selection changed after GEC plan generation"
        )
    mapping_metadata = payload.get("mapping")
    actual_mapping_sha256 = hashlib.sha256(
        plan.mapping.path.read_bytes()
    ).hexdigest()
    if (
        not isinstance(mapping_metadata, dict)
        or mapping_metadata.get("sha256") != actual_mapping_sha256
    ):
        raise GecCcpWorkflowError("GEC mapping changed after plan generation")
    model = payload.get("model", {})
    expected_mph = model.get("input_mph", {}) if isinstance(model, dict) else {}
    actual_mph_sha256 = hashlib.sha256(
        plan.mapping.model.input_mph.read_bytes()
    ).hexdigest()
    if expected_mph.get("sha256") != actual_mph_sha256:
        raise GecCcpWorkflowError("input MPH changed after GEC plan generation")
    bundle = payload.get("bundle")
    artifact_verification = (
        bundle.get("artifact_verification")
        if isinstance(bundle, dict)
        else None
    )
    expected_artifacts = (
        artifact_verification.get("sha256")
        if isinstance(artifact_verification, dict)
        else None
    )
    manifest_path = plan.mapping.bundle.path / "manifest.json"
    if (
        not isinstance(bundle, dict)
        or not isinstance(expected_artifacts, dict)
        or not expected_artifacts
        or not manifest_path.is_file()
    ):
        raise GecCcpWorkflowError("GEC plan lacks bundle artifact provenance")
    actual_manifest_sha256 = hashlib.sha256(
        manifest_path.read_bytes()
    ).hexdigest()
    if bundle.get("manifest_sha256") != actual_manifest_sha256:
        raise GecCcpWorkflowError(
            "bundle manifest changed after GEC plan generation"
        )
    try:
        current_manifest = json.loads(
            manifest_path.read_text(encoding="utf-8")
        )
    except (OSError, json.JSONDecodeError) as exc:
        raise GecCcpWorkflowError(
            "bundle manifest is unreadable at execution time"
        ) from exc
    current_tables = current_manifest.get("tables")
    if (
        not isinstance(current_tables, dict)
        or set(current_tables) != set(expected_artifacts)
    ):
        raise GecCcpWorkflowError(
            "bundle artifact set changed after GEC plan generation"
        )
    bundle_root = plan.mapping.bundle.path.resolve()
    actual_artifact_sha256: dict[str, str] = {}
    for name, expected_digest in expected_artifacts.items():
        relative_name = PurePosixPath(name) if isinstance(name, str) else None
        if (
            not isinstance(name, str)
            or not name
            or "\\" in name
            or relative_name is None
            or relative_name.is_absolute()
            or relative_name.as_posix() != name
            or any(part in {"", ".", ".."} for part in relative_name.parts)
            or re.fullmatch(r"[0-9a-f]{64}", str(expected_digest or ""))
            is None
        ):
            raise GecCcpWorkflowError(
                "GEC plan contains invalid bundle artifact provenance"
            )
        artifact_path = (bundle_root / Path(*relative_name.parts)).resolve()
        try:
            artifact_path.relative_to(bundle_root)
        except ValueError as exc:
            raise GecCcpWorkflowError(
                "GEC plan contains invalid bundle artifact provenance"
            ) from exc
        if not artifact_path.is_file():
            raise GecCcpWorkflowError(
                f"bundle artifact is absent at execution time: {name}"
            )
        actual_digest = hashlib.sha256(artifact_path.read_bytes()).hexdigest()
        manifest_digest = current_tables.get(name, {}).get("sha256")
        if (
            actual_digest != expected_digest
            or manifest_digest != expected_digest
        ):
            raise GecCcpWorkflowError(
                f"bundle artifact changed after GEC plan generation: {name}"
            )
        actual_artifact_sha256[name] = actual_digest
    expected_guard = bundle.get("low_energy_guard")
    primary_hashes = current_manifest.get("hashes")
    if not isinstance(primary_hashes, dict):
        raise GecCcpWorkflowError(
            "primary bundle hash provenance is unavailable"
        )
    actual_guard = _validate_low_energy_guard_bundle(
        plan.mapping,
        primary_cross_sections_sha256=str(
            primary_hashes.get("cross_sections_sha256", "")
        ),
    )
    if actual_guard != expected_guard:
        raise GecCcpWorkflowError(
            "low-energy guard bundle changed after GEC plan generation"
        )
    generated = payload.get("generated_java")
    if not isinstance(generated, dict):
        raise GecCcpWorkflowError("GEC plan lacks generated Java provenance")
    paths = {
        "apply": plan.apply_java,
        "external_run": plan.external_run_java,
        "external_export": plan.external_export_java,
    }
    if plan.mapping.run.include_builtin_reference:
        if plan.baseline_run_java is None or plan.baseline_export_java is None:
            raise GecCcpWorkflowError(
                "built-in reference Java is absent from the prepared plan"
            )
        paths.update(
            {
                "baseline_run": plan.baseline_run_java,
                "baseline_export": plan.baseline_export_java,
            }
        )
    if plan.native_eedf_audit is not None:
        paths["native_eedf_audit"] = plan.native_eedf_audit.java_path
    java_hashes: dict[str, str] = {}
    for name, path in paths.items():
        metadata = generated.get(name)
        if (
            not isinstance(metadata, dict)
            or metadata.get("materialized") is not True
            or not path.is_file()
        ):
            raise GecCcpWorkflowError(
                f"GEC plan Java artifact is not materialized: {name}"
            )
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if metadata.get("sha256") != digest:
            raise GecCcpWorkflowError(
                f"generated Java changed after GEC plan generation: {name}"
            )
        java_hashes[name] = digest
    native_provenance: dict[str, Any] = {}
    if plan.native_eedf_audit is not None:
        native_metadata = payload.get("native_function_eedf_audit")
        if not isinstance(native_metadata, dict):
            raise GecCcpWorkflowError(
                "GEC plan lacks native Function-EEDF audit provenance"
            )
        active_table = native_metadata.get("active_table")
        query_csv = native_metadata.get("query_csv")
        if not isinstance(active_table, dict) or not isinstance(query_csv, dict):
            raise GecCcpWorkflowError(
                "GEC plan native Function-EEDF audit provenance is invalid"
            )
        table_sha256 = hashlib.sha256(
            plan.native_eedf_audit.table_path.read_bytes()
        ).hexdigest()
        query_sha256 = hashlib.sha256(
            plan.native_eedf_audit.query_path.read_bytes()
        ).hexdigest()
        if active_table.get("sha256") != table_sha256:
            raise GecCcpWorkflowError(
                "active Function-EEDF table changed after GEC plan generation"
            )
        if query_csv.get("sha256") != query_sha256:
            raise GecCcpWorkflowError(
                "native Function-EEDF audit query changed after plan generation"
            )
        native_provenance = {
            "active_function_eedf_table_sha256": table_sha256,
            "native_function_eedf_query_sha256": query_sha256,
        }
    return {
        "mapping_sha256": actual_mapping_sha256,
        "input_mph_sha256": actual_mph_sha256,
        "generated_java_sha256": java_hashes,
        "bundle_artifacts_verified": True,
        "bundle_manifest_sha256": actual_manifest_sha256,
        "bundle_artifact_sha256": actual_artifact_sha256,
        **native_provenance,
    }


__all__ = ["validate_gec_plan_inputs"]
