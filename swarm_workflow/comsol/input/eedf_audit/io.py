"""Validated readers and scalar conversion for COMSOL Function-EEDF audits."""

from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path

import numpy as np

from ..function_eedf import (
    ComsolFunctionEedfGrid,
    FunctionEedfError,
    read_comsol_function_eedf_grid,
)
from .contracts import ComsolEedfAuditError, ComsolEedfAuditPlan


def _validate_audit_inputs(plan: ComsolEedfAuditPlan) -> None:
    for name, path, expected in (
        ("table", plan.table_path, plan.table_sha256),
        ("query plan", plan.query_path, plan.query_sha256),
        ("saved model", plan.model_path, plan.model_sha256),
    ):
        if not path.is_file():
            raise ComsolEedfAuditError(f"COMSOL audit {name} does not exist: {path}")
        if (
            expected is not None
            and hashlib.sha256(path.read_bytes()).hexdigest() != expected
        ):
            raise ComsolEedfAuditError(f"COMSOL audit {name} changed after preparation")


def read_native_eedf_grid(path: str | Path) -> ComsolFunctionEedfGrid:
    """Read the active grid through the shared COMSOL Function-EEDF parser."""

    try:
        return read_comsol_function_eedf_grid(path)
    except FunctionEedfError as exc:
        raise ComsolEedfAuditError(str(exc)) from exc


def _read_cross_sections(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    grouped: dict[str, list[tuple[float, float]]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        required = {"type", "energy_eV", "cross_section_m2"}
        if not required.issubset(reader.fieldnames or ()):
            raise ComsolEedfAuditError("cross-section CSV lacks required columns")
        for row in reader:
            grouped.setdefault(str(row["type"]), []).append(
                (
                    _finite_float(row["energy_eV"]),
                    _finite_float(row["cross_section_m2"]),
                )
            )
    result: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for process, points in grouped.items():
        points.sort()
        energy = np.asarray([item[0] for item in points])
        sigma = np.asarray([item[1] for item in points])
        if len(energy) < 2 or np.any(np.diff(energy) <= 0.0) or np.any(sigma < 0.0):
            raise ComsolEedfAuditError(f"invalid cross section: {process}")
        result[process] = (energy, sigma)
    return result


def _read_contract(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as stream:
        for line in stream:
            key, separator, value = line.rstrip("\r\n").partition("\t")
            if not separator or key in result:
                raise ComsolEedfAuditError("invalid COMSOL EEDF contract output")
            result[key] = value
    return result


def _read_rows(path: Path, columns: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != columns:
            raise ComsolEedfAuditError(f"unexpected audit CSV columns: {path}")
        return list(reader)


def _finite_float(value: object) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ComsolEedfAuditError(f"nonnumeric audit value: {value!r}") from exc
    if not math.isfinite(parsed):
        raise ComsolEedfAuditError(f"nonfinite audit value: {value!r}")
    return parsed
