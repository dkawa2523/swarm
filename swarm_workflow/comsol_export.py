"""Export built workflow table directories as COMSOL CSV bundles."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any

from ._io import write_csv as _write_csv
from ._io import write_json as _write_json
from .tables import UNITS


FORMAT_VERSION = 1

BASE_TABLES = (
    "mean_energy_vs_en.csv",
    "transport_vs_en.csv",
    "rates_vs_en.csv",
    "townsend_vs_en.csv",
    "quality.csv",
)
OPTIONAL_TABLES = (
    "rates_vs_mean_energy.csv",
    "transport_vs_mean_energy.csv",
    "eedf.csv",
)
COMSOL_FUNCTION_TABLES = (
    (
        "comsol_functions/sw_meanE.csv",
        "mean_energy_vs_en.csv",
        "mean_energy_eV",
        "mean_energy_eV",
    ),
    (
        "comsol_functions/sw_muN.csv",
        "transport_vs_en.csv",
        "reduced_mobility_m2_V_s_m3",
        "reduced_mobility_m2_V_s_m3",
    ),
)


class ComsolExportError(RuntimeError):
    """Raised when built tables cannot be exported safely for COMSOL."""


@dataclass(frozen=True, slots=True)
class ComsolExportSummary:
    table_directory: Path
    output_directory: Path
    bundles_written: int


def export_comsol_bundle(
    table_directory: str | Path,
    output_directory: str | Path,
) -> ComsolExportSummary:
    table_dir = Path(table_directory)
    output_dir = Path(output_directory)
    manifest = _read_manifest(table_dir)
    if _is_root_manifest(manifest):
        entries = manifest["mixtures"]
        output_dir.mkdir(parents=True, exist_ok=True)
        bundle_entries = []
        for entry in entries:
            source_manifest = table_dir / str(entry["path"])
            source_dir = source_manifest.parent
            target_dir = output_dir / f"mixture_{int(entry['mixture_id']):04d}"
            _export_single(source_dir, target_dir)
            bundle_entries.append(
                {
                    "mixture_id": int(entry["mixture_id"]),
                    "path": f"mixture_{int(entry['mixture_id']):04d}/manifest.json",
                }
            )
        _write_json(
            output_dir / "manifest.json",
            {
                "format_version": FORMAT_VERSION,
                "stage": "export-comsol",
                "source": manifest.get("source"),
                "hashes": manifest.get("hashes", {}),
                "bundles": bundle_entries,
            },
        )
        return ComsolExportSummary(table_dir, output_dir, len(bundle_entries))

    output_dir.mkdir(parents=True, exist_ok=True)
    _export_single(table_dir, output_dir)
    return ComsolExportSummary(table_dir, output_dir, 1)


def _export_single(
    table_dir: Path,
    output_dir: Path,
) -> None:
    source_manifest = _read_manifest(table_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    _clear_known_outputs(output_dir)
    missing = _missing_required_coefficients(table_dir)
    if missing:
        failed = _bundle_manifest(
            source_manifest,
            status="failed",
            tables={},
            missing_required_coefficients=missing,
        )
        _write_json(output_dir / "manifest.json", failed)
        raise ComsolExportError(
            "COMSOL export missing required coefficients: " + ", ".join(missing)
        )

    tables: dict[str, dict[str, Any]] = {}
    for name in (*BASE_TABLES, *OPTIONAL_TABLES):
        source = table_dir / name
        if source.exists():
            rows = _read_csv(source)
            columns = _csv_header(source)
            _write_csv(output_dir / name, columns, rows)
            tables[name] = _table_metadata(columns, source_manifest, name)
    tables.update(_write_comsol_function_tables(output_dir))
    tables.update(_write_eedf_probability_table(output_dir))
    _write_json(
        output_dir / "manifest.json",
        _bundle_manifest(
            source_manifest,
            status="ok",
            tables=tables,
            missing_required_coefficients=[],
        ),
    )


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
    for function_name, source_name, source_column, output_column in COMSOL_FUNCTION_TABLES:
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
        _write_csv(output_dir / function_name, columns, rows)
        tables[function_name] = {
            "argument": "E_over_N_V_m2",
            "columns": list(columns),
            "derived_for_comsol_interpolation": True,
            "units": {
                "E_over_N_V_m2": "V m^2",
                output_column: UNITS.get(output_column, _function_unit(output_column)),
            },
        }
    return tables


def _write_eedf_probability_table(
    output_dir: Path,
) -> dict[str, dict[str, Any]]:
    """Write COMSOL's EEPF convention f0 = EEDF/sqrt(energy).

    This is an audit/visualization artifact for the GEC CCP workflow.  The
    production COMSOL closure imports transport and reaction-rate lookups by
    mean electron energy; it does not try to force a spatially uniform EEDF.
    """

    source = output_dir / "eedf.csv"
    if not source.exists():
        return {}
    rows: list[dict[str, Any]] = []
    for row in _read_csv(source):
        energy = _required_float(row, "electron_energy_eV")
        eedf = _required_float(row, "eedf")
        if energy <= 0.0:
            raise ComsolExportError(
                "eedf.csv electron_energy_eV must be positive for EEPF conversion"
            )
        if eedf < 0.0:
            raise ComsolExportError("eedf.csv contains a negative EEDF value")
        rows.append(
            {
                "electron_energy_eV": energy,
                "mean_energy_eV": _required_float(row, "mean_energy_eV"),
                "E_over_N_Td": _required_float(row, "E_over_N_Td"),
                "eepf_eV_m32": eedf / math.sqrt(energy),
            }
        )
    if not rows:
        return {}
    name = "eedf_f0.csv"
    columns = (
        "electron_energy_eV",
        "mean_energy_eV",
        "E_over_N_Td",
        "eepf_eV_m32",
    )
    _write_csv(output_dir / name, columns, rows)
    return {
        name: {
            "argument": "electron_energy_eV,E_over_N_Td",
            "columns": list(columns),
            "derived_for_comsol_interpolation": True,
            "definition": "eepf_eV_m32 = eedf / sqrt(electron_energy_eV)",
            "units": {
                "electron_energy_eV": "eV",
                "mean_energy_eV": "eV",
                "E_over_N_Td": "Td",
                "eepf_eV_m32": "eV^-3/2",
            },
        }
    }


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
    return {
        "format_version": FORMAT_VERSION,
        "stage": "export-comsol",
        "status": status,
        "source": source_manifest.get("source"),
        "hashes": source_manifest.get("hashes", {}),
        "mixture": source_manifest.get("mixture"),
        "source_policy": source_manifest.get("source_policy", {}),
        "valid_ranges": source_manifest.get("valid_ranges", {}),
        "units": dict(sorted(units.items())),
        "table_argument": dict(source_manifest.get("table_argument", {})),
        "monotonicity": source_manifest.get("monotonicity", {}),
        "quality_summary": source_manifest.get("quality_summary", {}),
        "missing_required_coefficients": missing_required_coefficients,
        "tables": tables,
    }


def _table_metadata(
    columns: tuple[str, ...],
    source_manifest: dict[str, Any],
    name: str,
) -> dict[str, Any]:
    source_tables = source_manifest.get("tables", {})
    source_table = source_tables.get(name, {}) if isinstance(source_tables, dict) else {}
    return {
        "columns": list(columns),
        "units": {column: UNITS[column] for column in columns if column in UNITS},
        "argument": source_table.get("argument", "E_over_N_Td"),
    }


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
        "energy_loss.csv",
        *(name for name, *_ in COMSOL_FUNCTION_TABLES),
    ):
        path = directory / name
        if path.exists() and path.is_file():
            path.unlink()
    function_dir = directory / "comsol_functions"
    if function_dir.exists():
        for path in function_dir.glob("*.csv"):
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
