"""Manifest construction and source-evidence validation for COMSOL export."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from typing import Any

from swarm_workflow.quality.policy import (
    parse_quality_thresholds,
    quality_policy_provenance,
    quality_thresholds_payload,
)
from swarm_workflow.tables.contracts import (
    RATE_EVIDENCE_TABLE,
    UNITS,
    ZERO_EVENT_CONFIDENCE,
)

from .contracts import (
    COMSOL_FUNCTION_EEDF_TABLE,
    COPIED_TABLE_ROLES,
    FORMAT_VERSION,
    ComsolExportError,
)


def _canonical_eedf_source(
    table_dir: Path,
    output_dir: Path,
    source_manifest: dict[str, Any],
    listed_tables: dict[str, Any],
) -> Path:
    if "eedf.csv" not in listed_tables:
        return output_dir / "eedf.csv"
    policy = source_manifest.get("source_policy")
    postprocess = policy.get("postprocess") if isinstance(policy, dict) else None
    if postprocess not in {None, "none"}:
        raise ComsolExportError(
            "postprocessed or prior-filled EEDF tables are obsolete; rebuild "
            "a qualified pure-solver table with source_policy.postprocess=none"
        )
    source = table_dir / "eedf.csv"
    copied_source = output_dir / "eedf.csv"
    if _sha256_file(copied_source) != _sha256_file(source):
        raise ComsolExportError("copied Function-EEDF source hash mismatch")
    return copied_source


def _validate_source_artifact(
    path: Path,
    metadata: object,
    *,
    expected_roles: set[str],
) -> None:
    if not path.is_file() or not isinstance(metadata, dict):
        raise ComsolExportError(f"invalid source artifact: {path}")
    if metadata.get("artifact_role") not in expected_roles:
        raise ComsolExportError(f"invalid artifact role: {path.name}")
    expected_hash = metadata.get("sha256")
    actual_hash = _sha256_file(path)
    if expected_hash != actual_hash:
        raise ComsolExportError(f"artifact hash mismatch: {path.name}")


def _bundle_manifest(
    source_manifest: dict[str, Any],
    *,
    status: str,
    tables: dict[str, dict[str, Any]],
    missing_required_coefficients: list[str],
) -> dict[str, Any]:
    units = {
        column: unit
        for table in tables.values()
        for column, unit in table.get("units", {}).items()
    }
    raw_eedf = tables.get("eedf.csv")
    function_eedf_name = next(
        (
            name
            for name, metadata in tables.items()
            if metadata.get("artifact_role") == "canonical_comsol_function_eedf_input"
            and metadata.get("canonical_comsol_input") is True
        ),
        None,
    )
    function_eedf = (
        tables.get(function_eedf_name) if function_eedf_name is not None else None
    )
    return {
        "format_version": FORMAT_VERSION,
        "stage": "export-comsol",
        "status": status,
        "source": source_manifest.get("source"),
        "source_composition": source_manifest.get("source_composition"),
        "evidence": source_manifest.get("evidence"),
        "physical_context": source_manifest.get("physical_context"),
        "mc_campaign": source_manifest.get("mc_campaign"),
        "mc_qualification": source_manifest.get("mc_qualification"),
        "solver_qualification": source_manifest.get("solver_qualification"),
        "target_qualification": source_manifest.get("target_qualification"),
        "hashes": source_manifest.get("hashes", {}),
        "mc_sampling_plan": source_manifest.get("mc_sampling_plan"),
        "quality_thresholds": _validated_quality_thresholds(source_manifest),
        "quality_policy": _validated_quality_policy(source_manifest),
        "mixture": source_manifest.get("mixture"),
        "source_policy": source_manifest.get("source_policy", {}),
        "valid_ranges": source_manifest.get("valid_ranges", {}),
        "units": dict(sorted(units.items())),
        "table_argument": dict(source_manifest.get("table_argument", {})),
        "monotonicity": source_manifest.get("monotonicity", {}),
        "quality_summary": source_manifest.get("quality_summary", {}),
        "missing_required_coefficients": missing_required_coefficients,
        "eedf_artifacts": {
            "raw_audit": _eedf_artifact_reference("eedf.csv", raw_eedf),
            "canonicalization_source": _eedf_artifact_reference(
                "eedf.csv",
                raw_eedf,
            ),
            "comsol_input": _eedf_artifact_reference(
                function_eedf_name or COMSOL_FUNCTION_EEDF_TABLE,
                function_eedf,
            ),
        },
        "tables": tables,
    }


def _eedf_artifact_reference(
    name: str,
    metadata: dict[str, Any] | None,
) -> dict[str, str] | None:
    if not isinstance(metadata, dict):
        return None
    role = metadata.get("artifact_role")
    digest = metadata.get("sha256")
    if not isinstance(role, str) or not isinstance(digest, str):
        return None
    return {"file": name, "sha256": digest, "role": role}


def _sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _table_metadata(
    columns: tuple[str, ...],
    source_manifest: dict[str, Any],
    name: str,
) -> dict[str, Any]:
    source_tables = source_manifest.get("tables", {})
    source_table = (
        source_tables.get(name, {}) if isinstance(source_tables, dict) else {}
    )
    metadata = {
        "columns": list(columns),
        "units": {column: UNITS[column] for column in columns if column in UNITS},
        "argument": source_table.get("argument", "E_over_N_Td"),
        "artifact_role": COPIED_TABLE_ROLES[name],
    }
    if name == "quality.csv":
        for field in ("schema", "evidence_kind"):
            if field in source_table:
                metadata[field] = source_table[field]
    if name == RATE_EVIDENCE_TABLE:
        expected_statistics = {
            "replicate_interval": "two_sided_student_t_95",
            "zero_event_confidence": ZERO_EVENT_CONFIDENCE,
            "zero_event_upper_bound": (
                "poisson_zero_count_over_pooled_target_exposure"
            ),
        }
        if source_table.get("statistics") != expected_statistics:
            raise ComsolExportError(
                "rate_evidence.csv lacks the canonical pooled zero-event "
                "statistics contract"
            )
        metadata["statistics"] = expected_statistics
    return metadata


def _is_root_manifest(manifest: dict[str, Any]) -> bool:
    return isinstance(manifest.get("mixtures"), list)


def _read_manifest(directory: Path) -> dict[str, Any]:
    path = directory / "manifest.json"
    if not path.exists():
        raise ComsolExportError(f"missing table manifest: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ComsolExportError(f"invalid table manifest: {path}")
    return data


def _validated_quality_thresholds(manifest: dict[str, Any]) -> dict[str, object]:
    raw = manifest.get("quality_thresholds")
    if not isinstance(raw, dict):
        raise ComsolExportError(
            "source table manifest lacks quality_thresholds provenance"
        )
    try:
        canonical = quality_thresholds_payload(parse_quality_thresholds(raw))
    except (TypeError, ValueError) as exc:
        raise ComsolExportError(
            "source table manifest has invalid quality_thresholds provenance"
        ) from exc
    if raw != canonical:
        raise ComsolExportError(
            "source table manifest quality_thresholds is incomplete or noncanonical"
        )
    return canonical


def _validated_quality_policy(manifest: dict[str, Any]) -> dict[str, object]:
    raw = manifest.get("quality_policy")
    if not isinstance(raw, dict):
        raise ComsolExportError("source table manifest lacks quality_policy")
    source_raw = raw.get("source")
    evaluation_raw = raw.get("evaluation")
    try:
        source = parse_quality_thresholds(source_raw)
        evaluation = parse_quality_thresholds(evaluation_raw)
    except (TypeError, ValueError) as exc:
        raise ComsolExportError(
            "source table manifest has invalid quality_policy"
        ) from exc
    canonical = quality_policy_provenance(
        source,
        evaluation,
        source_json=json.dumps(
            source_raw,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ),
    )
    if raw != canonical:
        raise ComsolExportError(
            "source table manifest quality_policy is incomplete or noncanonical"
        )
    if quality_thresholds_payload(evaluation) != _validated_quality_thresholds(
        manifest
    ):
        raise ComsolExportError(
            "source table evaluation policy differs from quality_thresholds"
        )
    return canonical
