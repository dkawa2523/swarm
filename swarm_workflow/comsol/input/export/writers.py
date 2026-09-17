"""Common CSV and artifact writers for COMSOL export."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any

from swarm_workflow._io import write_csv as _write_csv
from swarm_workflow.comsol.input.bundle_selection import (
    SELECTION_FILE,
    SOURCE_MANIFEST_FILE,
)
from swarm_workflow.quality.solver import (
    PROPAGATOR_CORE_QUALIFICATION_FILE,
    PROPAGATOR_TARGET_QUALIFICATION_FILE,
)
from swarm_workflow.tables.contracts import UNITS

from .contracts import (
    BASE_TABLES,
    COMSOL_FUNCTION_TABLES,
    COMSOL_FUNCTION_EEDF_TABLE,
    DERIVED_FUNCTION_ROLE,
    FUNCTION_EEDF_TABLE,
    OPTIONAL_TABLES,
    ComsolExportError,
)
from .manifest import _sha256_file


def _missing_required_coefficients(table_dir: Path) -> list[str]:
    missing: list[str] = []
    rows = _read_csv(table_dir / "transport_vs_en.csv")
    if not rows:
        return ["transport_vs_en.csv"]
    for row in rows:
        e_over_n = row.get("E_over_N_Td", "?")
        if _float_or_none(row.get("reduced_mobility_m2_V_s_m3")) is None:
            missing.append(f"reduced_mobility_m2_V_s_m3@{e_over_n}")
    return sorted(set(missing))


def _write_comsol_function_tables(output_dir: Path) -> dict[str, dict[str, Any]]:
    tables: dict[str, dict[str, Any]] = {}
    for (
        function_name,
        source_name,
        source_column,
        output_column,
    ) in COMSOL_FUNCTION_TABLES:
        source = output_dir / source_name
        if not source.exists():
            continue
        rows = _function_rows(
            _read_csv(source),
            source_column=source_column,
            output_column=output_column,
        )
        if not rows:
            continue
        columns = ("E_over_N_V_m2", output_column)
        target = output_dir / function_name
        _write_csv(target, columns, rows)
        tables[function_name] = {
            "argument": "E_over_N_V_m2",
            "artifact_role": DERIVED_FUNCTION_ROLE,
            "columns": list(columns),
            "derived_for_comsol_interpolation": True,
            "units": {
                "E_over_N_V_m2": "V m^2",
                output_column: UNITS.get(output_column, _function_unit(output_column)),
            },
            "sha256": _sha256_file(target),
        }
    return tables


def _function_rows(
    rows: list[dict[str, str]],
    *,
    source_column: str,
    output_column: str,
) -> list[dict[str, Any]]:
    return [
        {
            "E_over_N_V_m2": _required_float(row, "E_over_N_V_m2"),
            output_column: _required_float(row, source_column),
        }
        for row in rows
    ]


def _function_unit(column: str) -> str:
    if column.endswith("_townsend_m2"):
        return "m^2"
    if column.endswith("_rate_coefficient_m3_s"):
        return "m^3/s"
    if column.endswith("_energy_loss_rate_coefficient_eV_m3_s"):
        return "eV m^3/s"
    return ""


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise ComsolExportError(f"missing required table: {path.name}")
    with path.open("r", encoding="utf-8", newline="") as fp:
        return list(csv.DictReader(fp))


def _csv_header(path: Path) -> tuple[str, ...]:
    with path.open("r", encoding="utf-8", newline="") as fp:
        reader = csv.reader(fp)
        return tuple(next(reader))


def _clear_known_outputs(directory: Path) -> None:
    for name in (
        "manifest.json",
        *BASE_TABLES,
        *OPTIONAL_TABLES,
        "eedf_f0.csv",
        FUNCTION_EEDF_TABLE,
        COMSOL_FUNCTION_EEDF_TABLE,
        # Remove artifacts produced by the retired native-Grid binding when an
        # existing bundle directory is regenerated with the Spreadsheet path.
        "eedf_f0_comsol_grid.txt",
        "eedf_f0_binding_seed.csv",
        SELECTION_FILE,
        SOURCE_MANIFEST_FILE,
        PROPAGATOR_CORE_QUALIFICATION_FILE,
        PROPAGATOR_TARGET_QUALIFICATION_FILE,
        *(name for name, *_ in COMSOL_FUNCTION_TABLES),
    ):
        path = directory / name
        if path.exists() and path.is_file():
            path.unlink()
    function_dir = directory / "comsol_functions"
    if function_dir.exists():
        for pattern in ("*.csv", "*.txt"):
            for path in function_dir.glob(pattern):
                path.unlink()


def _required_float(row: dict[str, Any], name: str) -> float:
    value = _float_or_none(row.get(name))
    if value is None:
        raise ComsolExportError(f"missing required value {name}")
    return value


def _float_or_none(value: object) -> float | None:
    if value is None:
        return None
    if isinstance(value, str) and not value.strip():
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None
