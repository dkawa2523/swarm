"""Resolve a local COMSOL installation and optional license file."""

from __future__ import annotations

import os
from pathlib import Path
import shutil

from .types import ComsolAdapterError, ResolvedComsolExecutable


def resolve_comsol_executable(
    comsol_executable: str | Path | None = None,
) -> ResolvedComsolExecutable:
    if comsol_executable is not None:
        return _resolve_candidate(str(comsol_executable), source="cli")
    for env_name in ("COMSOL_BATCH", "COMSOL_EXECUTABLE"):
        value = os.environ.get(env_name)
        if value:
            return _resolve_candidate(value, source=f"env:{env_name}")
    for name in ("comsol", "comsolbatch"):
        found = shutil.which(name)
        if found:
            return _resolved_candidate(Path(found).resolve(), source=f"path:{name}")
    raise ComsolAdapterError(
        "COMSOL executable not found; pass --comsol, set COMSOL_BATCH or "
        "COMSOL_EXECUTABLE, or put comsol/comsolbatch on PATH"
    )


def _resolve_candidate(value: str, *, source: str) -> ResolvedComsolExecutable:
    path = Path(value)
    if path.exists():
        return _resolved_candidate(path.resolve(), source=source)
    found = shutil.which(value)
    if found:
        return _resolved_candidate(Path(found).resolve(), source=source)
    raise ComsolAdapterError(f"COMSOL executable does not exist: {value}")


def _resolved_candidate(
    path: Path,
    *,
    source: str,
) -> ResolvedComsolExecutable:
    license_path, license_source = _resolve_comsol_license_file(path)
    return ResolvedComsolExecutable(
        path=str(path),
        source=source,
        license_path=str(license_path) if license_path is not None else None,
        license_source=license_source,
    )


def _resolve_comsol_license_file(
    executable: Path,
) -> tuple[Path | None, str | None]:
    configured = os.environ.get("COMSOL_LICENSE_FILE")
    if configured:
        path = Path(configured).expanduser().resolve()
        if not path.is_file():
            raise ComsolAdapterError(
                "COMSOL_LICENSE_FILE does not name a readable license file: "
                f"{path}"
            )
        return path, "env:COMSOL_LICENSE_FILE"
    for parent in executable.parents:
        candidate = parent / "license" / "license.dat"
        if candidate.is_file():
            return candidate.resolve(), "installation:license/license.dat"
    return None, None
