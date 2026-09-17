"""Run directories and provenance for model-independent COMSOL execution."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
from pathlib import Path
import re
from typing import Any

from ..._io import write_json
from .types import ComsolRunContext, ResolvedComsolExecutable


def new_operation_log_dir(log_path: Path, operation: str) -> Path:
    log_path.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    candidate = log_path / f"{operation}_{stamp}"
    suffix = 1
    while candidate.exists():
        suffix += 1
        candidate = log_path / f"{operation}_{stamp}_{suffix}"
    candidate.mkdir(parents=True)
    return candidate


def file_artifact(path: Path | None) -> dict[str, Any] | None:
    if path is None:
        return None
    resolved = path.resolve()
    artifact: dict[str, Any] = {
        "path": str(resolved),
        "exists": resolved.is_file(),
    }
    if not resolved.is_file():
        return artifact
    stat = resolved.stat()
    artifact.update(
        {
            "size_bytes": stat.st_size,
            "modified_at_utc": _timestamp_utc(stat.st_mtime),
            "sha256": sha256_file(resolved),
        }
    )
    return artifact


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_execution_provenance(
    path: Path,
    *,
    status: str,
    operation: str,
    context: ComsolRunContext,
    resolved: ResolvedComsolExecutable,
    java_path: Path | None,
    source_java_path: Path | None,
    support_java_paths: tuple[Path, ...],
    source_support_java_paths: tuple[Path, ...],
    started_at_utc: str,
    results: list[dict[str, Any]],
    total_time_s: float | None,
    preexisting_output_mph: dict[str, Any] | None,
    study: str,
) -> dict[str, Any]:
    identity = extract_comsol_identity(results)
    source_origin = source_java_path if source_java_path is not None else java_path
    support_origins = (
        source_support_java_paths if source_support_java_paths else support_java_paths
    )
    payload = {
        "format_version": 1,
        "status": status,
        "operation": operation,
        "recorded_at_utc": utc_now(),
        "started_at_utc": started_at_utc,
        "total_time_s": total_time_s,
        "comsol": {
            "executable": resolved.path,
            "executable_source": resolved.source,
            "license_file": file_artifact(
                Path(resolved.license_path)
                if resolved.license_path is not None
                else None
            ),
            "license_source": resolved.license_source,
            "version": identity.get("version"),
            "build": identity.get("build"),
            "progress_version": identity.get("progress_version"),
        },
        "inputs": {
            "mapping": file_artifact(context.mapping_path),
            "input_mph": file_artifact(context.input_mph),
            "bundle_manifest": file_artifact(
                context.bundle_path / "manifest.json"
                if context.bundle_path is not None
                else None
            ),
            "preexisting_output_mph": preexisting_output_mph,
        },
        "generated_artifacts": {
            "source_origin": file_artifact(source_origin),
            "staged_java": file_artifact(java_path),
            "compiled_class": file_artifact(
                java_path.with_suffix(".class") if java_path is not None else None
            ),
            "support_classes": [
                {
                    "source_origin": file_artifact(origin),
                    "staged_java": file_artifact(staged),
                    "compiled_class": file_artifact(staged.with_suffix(".class")),
                }
                for origin, staged in zip(support_origins, support_java_paths)
            ],
            "output_mph": file_artifact(context.output_mph),
        },
        "study": study,
        "steps": [
            {
                "name": result.get("name"),
                "status": result.get("status"),
                "return_code": result.get("return_code"),
                "detected_comsol_error": result.get("detected_comsol_error"),
            }
            for result in results
        ],
    }
    write_json(path, payload)
    return payload


def extract_comsol_identity(results: list[dict[str, Any]]) -> dict[str, str]:
    texts: list[str] = []
    for result in results:
        for key in ("stdout", "stderr"):
            value = result.get(key)
            if isinstance(value, str) and Path(value).is_file():
                texts.append(Path(value).read_text(encoding="utf-8", errors="replace"))
    joined = "\n".join(texts)
    identity: dict[str, str] = {}
    banner = re.search(
        r"COMSOL Multiphysics\s+([0-9.]+)\s+\(Build:\s*([^)]+)\)",
        joined,
        flags=re.IGNORECASE,
    )
    if banner is not None:
        identity["version"] = banner.group(1)
        identity["build"] = banner.group(2).strip()
    progress = re.search(
        r"\*+COMSOL\s+([0-9]+(?:\.[0-9]+)+)\s+progress output file",
        joined,
        flags=re.IGNORECASE,
    )
    if progress is not None:
        identity["progress_version"] = progress.group(1)
    return identity


def utc_now() -> str:
    return datetime.now(tz=timezone.utc).isoformat().replace("+00:00", "Z")


def _timestamp_utc(value: float) -> str:
    return (
        datetime.fromtimestamp(value, tz=timezone.utc)
        .isoformat()
        .replace("+00:00", "Z")
    )
