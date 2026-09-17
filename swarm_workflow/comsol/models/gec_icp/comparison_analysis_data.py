"""Fail-closed loaders and numerical checks for GEC-ICP comparisons."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

from .comparison_analysis_contracts import (
    COIL_COLUMNS,
    GecIcpComparisonAnalysisError,
    GecIcpComparisonCaseData,
    GecIcpIntegratedState,
    METRIC_KEYS,
    SOLUTION_TIME_COLUMNS,
    VOLUME_COLUMNS,
)
from .comparison_export import EXPECTED_EXPORT_NAMES


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_comparison_case(
    case_id: str,
    label: str,
    export_directory: str | Path,
    *,
    common_time_s: float,
    stationarity_limits: Mapping[str, float],
) -> GecIcpComparisonCaseData:
    """Load one exporter directory only after its provenance closes."""

    directory = Path(export_directory).resolve()
    if not directory.is_dir():
        raise GecIcpComparisonAnalysisError(
            f"comparison case directory does not exist: {directory}"
        )
    manifest_path = directory / "comparison_export_manifest.json"
    manifest = _load_and_validate_manifest(case_id, directory, manifest_path)
    source_hashes = {
        name: sha256_file(directory / name) for name in EXPECTED_EXPORT_NAMES
    }

    saved_times = _read_solution_times(directory / "solution_times.csv")
    common = _load_state(
        directory / "volume_common_time.csv",
        directory / "coil_power_common_time.csv",
        expected_rows=1,
        context="common-time",
    )[0]
    terminal = _load_state(
        directory / "volume_terminal.csv",
        directory / "coil_power_terminal.csv",
        expected_rows=1,
        context="terminal",
    )[0]
    series = _load_state(
        directory / "volume_time_series.csv",
        directory / "coil_power_time_series.csv",
        expected_rows=None,
        context="time-series",
    )
    _require_time(common.time_s, common_time_s, "common-time export")
    manifest_time = _finite_float(
        manifest.get("time_selection", {}).get("common", {}).get("time_s"),
        "manifest common time",
    )
    _require_time(manifest_time, common_time_s, "manifest common time")
    _require_matching_time_axes(saved_times, tuple(row.time_s for row in series))
    _require_time(terminal.time_s, saved_times[-1], "terminal saved time")
    _require_state_close(terminal, series[-1], "terminal versus time-series")

    if len(series) < 4:
        raise GecIcpComparisonAnalysisError(
            f"{case_id}: stationarity requires at least four saved states"
        )
    stationarity_changes = {
        key: _maximum_adjacent_relative_change(
            [float(getattr(row, key)) for row in series[-4:]]
        )
        for key in METRIC_KEYS
    }
    stationarity_passed = {
        key: stationarity_changes[key] <= stationarity_limits[key]
        for key in METRIC_KEYS
    }
    input_record = manifest["input_mph"]
    mapping_record = manifest["model_mapping"]
    return GecIcpComparisonCaseData(
        case_id=case_id,
        label=label,
        export_directory=directory,
        input_mph=Path(input_record["path"]).resolve(),
        input_mph_sha256=str(input_record["sha256"]),
        source_manifest=manifest_path,
        source_manifest_sha256=sha256_file(manifest_path),
        model_mapping=Path(mapping_record["path"]).resolve(),
        model_mapping_sha256=str(mapping_record["sha256"]),
        common=common,
        terminal=terminal,
        time_series=tuple(series),
        solution_times_s=saved_times,
        source_hashes=source_hashes,
        stationarity_changes=stationarity_changes,
        stationarity_limits=dict(stationarity_limits),
        stationarity_passed=stationarity_passed,
    )


def _load_and_validate_manifest(
    case_id: str,
    directory: Path,
    manifest_path: Path,
) -> dict[str, object]:
    if not manifest_path.is_file():
        raise GecIcpComparisonAnalysisError(
            f"{case_id}: missing comparison export manifest"
        )
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise GecIcpComparisonAnalysisError(
            f"{case_id}: invalid comparison export manifest"
        ) from exc
    if not isinstance(payload, dict):
        raise GecIcpComparisonAnalysisError(
            f"{case_id}: comparison export manifest must be an object"
        )
    required = {
        "schema_version": 1,
        "stage": "execute-gec-icp-comparison-export",
        "status": "completed",
        "case_id": case_id,
        "read_only": True,
    }
    for key, expected in required.items():
        if payload.get(key) != expected:
            raise GecIcpComparisonAnalysisError(
                f"{case_id}: manifest {key!r} is not {expected!r}"
            )

    input_record = payload.get("input_mph")
    if not isinstance(input_record, dict):
        raise GecIcpComparisonAnalysisError(
            f"{case_id}: manifest lacks input MPH provenance"
        )
    source = Path(str(input_record.get("path", ""))).resolve()
    recorded_hash = input_record.get("sha256")
    if (
        not source.is_file()
        or not _valid_sha256(recorded_hash)
        or input_record.get("post_execution_sha256") != recorded_hash
        or sha256_file(source) != recorded_hash
    ):
        raise GecIcpComparisonAnalysisError(
            f"{case_id}: input MPH provenance is not closed"
        )

    mapping_record = payload.get("model_mapping")
    if not isinstance(mapping_record, dict):
        raise GecIcpComparisonAnalysisError(
            f"{case_id}: manifest lacks model mapping provenance"
        )
    mapping_path = Path(str(mapping_record.get("path", ""))).resolve()
    mapping_hash = mapping_record.get("sha256")
    if (
        not mapping_path.is_file()
        or not _valid_sha256(mapping_hash)
        or sha256_file(mapping_path) != mapping_hash
    ):
        raise GecIcpComparisonAnalysisError(
            f"{case_id}: model mapping provenance is not closed"
        )

    expected_outputs = payload.get("expected_outputs")
    if not isinstance(expected_outputs, list):
        raise GecIcpComparisonAnalysisError(
            f"{case_id}: manifest expected_outputs must be a list"
        )
    expected_names = [Path(str(value)).name for value in expected_outputs]
    if len(expected_names) != len(set(expected_names)) or set(expected_names) != set(
        EXPECTED_EXPORT_NAMES
    ):
        raise GecIcpComparisonAnalysisError(
            f"{case_id}: manifest has an incomplete output contract"
        )

    output_records = payload.get("outputs")
    if not isinstance(output_records, list):
        raise GecIcpComparisonAnalysisError(
            f"{case_id}: manifest lacks completed output hashes"
        )
    by_name: dict[str, dict[str, object]] = {}
    for record in output_records:
        if not isinstance(record, dict):
            raise GecIcpComparisonAnalysisError(
                f"{case_id}: malformed output provenance record"
            )
        name = Path(str(record.get("path", ""))).name
        if name in by_name:
            raise GecIcpComparisonAnalysisError(
                f"{case_id}: duplicate output provenance for {name}"
            )
        by_name[name] = record
    if set(by_name) != set(EXPECTED_EXPORT_NAMES):
        raise GecIcpComparisonAnalysisError(
            f"{case_id}: completed output hashes do not cover the export contract"
        )
    for name in EXPECTED_EXPORT_NAMES:
        path = directory / name
        record = by_name[name]
        recorded_output_hash = record.get("sha256")
        if (
            not path.is_file()
            or path.stat().st_size <= 0
            or record.get("size_bytes") != path.stat().st_size
            or not _valid_sha256(recorded_output_hash)
            or sha256_file(path) != recorded_output_hash
        ):
            raise GecIcpComparisonAnalysisError(
                f"{case_id}: output provenance does not match {name}"
            )
    return payload


def _read_solution_times(path: Path) -> tuple[float, ...]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as stream:
            rows = list(csv.reader(stream))
    except (OSError, UnicodeError, csv.Error) as exc:
        raise GecIcpComparisonAnalysisError(
            f"invalid solution-time export: {path}"
        ) from exc
    if not rows or tuple(value.strip() for value in rows[0]) != SOLUTION_TIME_COLUMNS:
        raise GecIcpComparisonAnalysisError(
            f"solution-time header does not match its contract: {path}"
        )
    if len(rows) < 2:
        raise GecIcpComparisonAnalysisError(f"empty solution-time export: {path}")
    values: list[float] = []
    for expected_index, row in enumerate(rows[1:], start=1):
        if len(row) != 2:
            raise GecIcpComparisonAnalysisError(f"malformed solution-time row: {path}")
        try:
            index_value = int(row[0])
        except ValueError as exc:
            raise GecIcpComparisonAnalysisError(
                f"invalid solution index: {path}"
            ) from exc
        if index_value != expected_index:
            raise GecIcpComparisonAnalysisError(
                f"solution indices are not contiguous: {path}"
            )
        values.append(_finite_float(row[1], f"solution time in {path}"))
    _require_strictly_increasing(values, f"solution times in {path}")
    if values[0] < 0.0:
        raise GecIcpComparisonAnalysisError(
            f"solution time axis starts below zero: {path}"
        )
    return tuple(values)


def _read_comsol_table(
    path: Path,
    expected_columns: tuple[str, ...],
) -> list[list[float]]:
    try:
        lines = path.read_text(encoding="utf-8-sig").splitlines()
    except (OSError, UnicodeError) as exc:
        raise GecIcpComparisonAnalysisError(
            f"cannot read COMSOL table: {path}"
        ) from exc
    header_index: int | None = None
    for index, line in enumerate(lines):
        stripped = line.lstrip()
        if not stripped.startswith("%"):
            continue
        cells = tuple(
            value.strip() for value in next(csv.reader([stripped[1:].lstrip()]))
        )
        if cells == expected_columns:
            if header_index is not None:
                raise GecIcpComparisonAnalysisError(
                    f"duplicate COMSOL table header: {path}"
                )
            header_index = index
    if header_index is None:
        raise GecIcpComparisonAnalysisError(
            f"COMSOL table header does not match its contract: {path}"
        )
    rows: list[list[float]] = []
    for line in lines[header_index + 1 :]:
        if not line.strip():
            continue
        if line.lstrip().startswith("%"):
            raise GecIcpComparisonAnalysisError(
                f"unexpected COMSOL comment after table header: {path}"
            )
        try:
            cells = next(csv.reader([line]))
        except csv.Error as exc:
            raise GecIcpComparisonAnalysisError(
                f"invalid COMSOL CSV row: {path}"
            ) from exc
        if len(cells) != len(expected_columns):
            raise GecIcpComparisonAnalysisError(
                f"COMSOL table row has the wrong width: {path}"
            )
        rows.append(
            [_finite_float(value, f"numeric value in {path}") for value in cells]
        )
    if not rows:
        raise GecIcpComparisonAnalysisError(f"empty COMSOL table: {path}")
    return rows


def _load_state(
    volume_path: Path,
    coil_path: Path,
    *,
    expected_rows: int | None,
    context: str,
) -> list[GecIcpIntegratedState]:
    volume_rows = _read_comsol_table(volume_path, VOLUME_COLUMNS)
    coil_rows = _read_comsol_table(coil_path, COIL_COLUMNS)
    if expected_rows is not None and (
        len(volume_rows) != expected_rows or len(coil_rows) != expected_rows
    ):
        raise GecIcpComparisonAnalysisError(
            f"{context} export must contain exactly {expected_rows} row(s)"
        )
    if len(volume_rows) != len(coil_rows):
        raise GecIcpComparisonAnalysisError(
            f"{context} volume and coil row counts differ"
        )
    times = [row[0] for row in volume_rows]
    if len(times) > 1:
        _require_strictly_increasing(times, f"{context} time axis")
    states: list[GecIcpIntegratedState] = []
    for volume, coil in zip(volume_rows, coil_rows, strict=True):
        _require_time(volume[0], coil[0], f"{context} volume/coil time")
        if volume[1] <= 0.0 or volume[2] <= 0.0:
            raise GecIcpComparisonAnalysisError(
                f"{context} has a nonpositive volume or electron inventory"
            )
        mean_energy = volume[3] / volume[2]
        positive = (mean_energy, volume[5], coil[1])
        nonnegative = (volume[4], volume[6])
        if any(value <= 0.0 for value in positive) or any(
            value < 0.0 for value in nonnegative
        ):
            raise GecIcpComparisonAnalysisError(
                f"{context} contains an unphysical integrated quantity"
            )
        states.append(
            GecIcpIntegratedState(
                time_s=volume[0],
                axisymmetric_volume_m3=volume[1],
                electron_inventory=volume[2],
                electron_weighted_mean_energy_eV=mean_energy,
                metastable_inventory=volume[4],
                ion_inventory=volume[5],
                absorbed_power_W=volume[6],
                coil_power_W=coil[1],
            )
        )
    return states


def _finite_float(value: object, context: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise GecIcpComparisonAnalysisError(f"invalid {context}") from exc
    if not math.isfinite(parsed):
        raise GecIcpComparisonAnalysisError(f"nonfinite {context}")
    return parsed


def _valid_sha256(value: object) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdef" for character in value)
    )


def _require_strictly_increasing(values: Sequence[float], context: str) -> None:
    if any(right <= left for left, right in zip(values, values[1:])):
        raise GecIcpComparisonAnalysisError(f"{context} is not strictly increasing")


def _require_time(observed: float, expected: float, context: str) -> None:
    if not math.isclose(observed, expected, rel_tol=1.0e-10, abs_tol=1.0e-15):
        raise GecIcpComparisonAnalysisError(
            f"{context} differs: observed {observed:.17g}, expected {expected:.17g}"
        )


def _require_matching_time_axes(
    saved: Sequence[float],
    exported: Sequence[float],
) -> None:
    if len(saved) != len(exported):
        raise GecIcpComparisonAnalysisError(
            "saved-solution and integrated time axes have different lengths"
        )
    for left, right in zip(saved, exported, strict=True):
        _require_time(right, left, "saved-solution time axis")


def _require_state_close(
    left: GecIcpIntegratedState,
    right: GecIcpIntegratedState,
    context: str,
) -> None:
    for key in (
        "time_s",
        "axisymmetric_volume_m3",
        *METRIC_KEYS,
    ):
        left_value = float(getattr(left, key))
        right_value = float(getattr(right, key))
        if not math.isclose(left_value, right_value, rel_tol=1.0e-9, abs_tol=1.0e-15):
            raise GecIcpComparisonAnalysisError(f"{context} differs for {key}")


def _relative_change(left: float, right: float) -> float:
    return abs(right - left) / max(abs(left), abs(right), 1.0e-300)


def _maximum_adjacent_relative_change(values: Sequence[float]) -> float:
    return max(_relative_change(left, right) for left, right in zip(values, values[1:]))


__all__ = ("load_comparison_case", "sha256_file")
