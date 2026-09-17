"""Small model-independent contracts for local COMSOL execution."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


class ComsolAdapterError(RuntimeError):
    """Raised when local COMSOL batch execution cannot be planned or run."""

    def __init__(
        self,
        message: str,
        *,
        step_result: dict[str, Any] | None = None,
    ) -> None:
        super().__init__(message)
        self.step_result = step_result


@dataclass(frozen=True, slots=True)
class ComsolRunContext:
    """Filesystem inputs required by the COMSOL runtime.

    Model mappings remain owned by model packages.  The runtime only needs
    concrete paths for execution and provenance.
    """

    root: Path
    mapping_path: Path
    input_mph: Path
    output_mph: Path
    log_path: Path
    bundle_path: Path | None = None


@dataclass(frozen=True, slots=True)
class ResolvedComsolExecutable:
    path: str
    source: str
    license_path: str | None
    license_source: str | None


@dataclass(frozen=True, slots=True)
class ComsolExecutionSummary:
    operation: str
    log_dir: Path
    command_json: Path
    result_json: Path
    stdout_paths: tuple[Path, ...]
    stderr_paths: tuple[Path, ...]
    return_codes: tuple[int, ...]
    java_path: Path | None
    output_mph: Path
    total_time_s: float
    provenance_json: Path | None = None
    comsol_version: str | None = None
    comsol_build: str | None = None
