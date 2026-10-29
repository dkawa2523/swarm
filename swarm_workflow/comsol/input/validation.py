"""Model-independent integrity checks for manifest-listed CSV tables."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
import csv
from dataclasses import dataclass
from hashlib import sha256
import io
from pathlib import Path, PurePosixPath, PureWindowsPath
import re
from typing import Any


class ComsolInputArtifactError(ValueError):
    """A manifest table is missing, changed, or outside its input directory."""


@dataclass(frozen=True, slots=True)
class ValidatedTableArtifact:
    name: str
    path: Path
    columns: tuple[str, ...]
    sha256: str


def contained_artifact_path(directory: str | Path, name: object) -> Path:
    """Resolve one canonical relative artifact name without following an escape."""

    root = Path(directory).resolve()
    if not isinstance(name, str) or not name or "\\" in name or "\x00" in name:
        raise ComsolInputArtifactError(f"invalid artifact path: {name!r}")
    relative = PurePosixPath(name)
    if (
        relative.is_absolute()
        or PureWindowsPath(name).drive
        or relative.as_posix() != name
        or any(part in {"", ".", ".."} for part in relative.parts)
    ):
        raise ComsolInputArtifactError(f"invalid artifact path: {name!r}")
    candidate = (root / Path(*relative.parts)).resolve()
    if not candidate.is_relative_to(root) or not candidate.is_file():
        raise ComsolInputArtifactError(f"invalid or missing artifact: {name}")
    return candidate


def validate_manifest_tables(
    directory: str | Path,
    manifest: Mapping[str, Any],
    *,
    required_names: Iterable[str] = (),
    allowed_names: Iterable[str] | None = None,
    require_roles: bool = False,
) -> dict[str, ValidatedTableArtifact]:
    """Verify every listed table's path, SHA256, and exact CSV column contract.

    Unknown listed tables fail when the caller supplies its supported inventory.
    Unlisted files are never consumed. Solver identity, physics, and qualification
    are separate contracts owned by the caller.
    """

    tables = manifest.get("tables")
    if not isinstance(tables, dict) or not tables:
        raise ComsolInputArtifactError("manifest.tables must be a nonempty mapping")
    missing = set(required_names).difference(tables)
    if missing:
        raise ComsolInputArtifactError(
            "manifest omits required tables: " + ", ".join(sorted(missing))
        )
    if allowed_names is not None:
        unknown = set(tables).difference(allowed_names)
        if unknown:
            raise ComsolInputArtifactError(
                "manifest lists unsupported tables: "
                + ", ".join(sorted(str(name) for name in unknown))
            )
    validated: dict[str, ValidatedTableArtifact] = {}
    for name, metadata in tables.items():
        path = contained_artifact_path(directory, name)
        if not isinstance(metadata, dict):
            raise ComsolInputArtifactError(f"invalid table metadata: {name}")
        expected_hash = metadata.get("sha256")
        if (
            not isinstance(expected_hash, str)
            or re.fullmatch(r"[0-9a-f]{64}", expected_hash) is None
        ):
            raise ComsolInputArtifactError(f"table lacks a valid SHA256: {name}")
        columns = metadata.get("columns")
        if (
            not isinstance(columns, list)
            or not columns
            or any(not isinstance(column, str) or not column.strip() for column in columns)
            or len(columns) != len(set(columns))
        ):
            raise ComsolInputArtifactError(f"invalid manifest columns: {name}")
        role = metadata.get("artifact_role")
        if require_roles and (not isinstance(role, str) or not role.strip()):
            raise ComsolInputArtifactError(f"table lacks an artifact role: {name}")
        try:
            content = path.read_bytes()
        except OSError as exc:
            raise ComsolInputArtifactError(f"cannot read table: {name}") from exc
        if sha256(content).hexdigest() != expected_hash:
            raise ComsolInputArtifactError(f"artifact hash mismatch: {name}")
        try:
            reader = csv.reader(io.StringIO(content.decode("utf-8-sig")), strict=True)
            header = next(reader, None)
            if header != columns:
                raise ComsolInputArtifactError(f"CSV columns differ from manifest: {name}")
            if any(len(row) != len(columns) for row in reader):
                raise ComsolInputArtifactError(f"CSV row width differs from columns: {name}")
        except (UnicodeError, csv.Error) as exc:
            raise ComsolInputArtifactError(f"invalid CSV table: {name}") from exc
        validated[name] = ValidatedTableArtifact(name, path, tuple(columns), expected_hash)
    return validated


__all__ = (
    "ComsolInputArtifactError",
    "ValidatedTableArtifact",
    "contained_artifact_path",
    "validate_manifest_tables",
)
