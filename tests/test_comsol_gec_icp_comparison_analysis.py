from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

from swarm_workflow.comsol.models.gec_icp.comparison_analysis import (
    generate_gec_icp_saved_solution_comparison,
)
from swarm_workflow.comsol.models.gec_icp.comparison_analysis_contracts import (
    COIL_COLUMNS,
    GecIcpComparisonAnalysisError,
    METRIC_KEYS,
    VOLUME_COLUMNS,
)
from swarm_workflow.comsol.models.gec_icp.comparison_export import (
    EXPECTED_EXPORT_NAMES,
)


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _table(path: Path, columns: tuple[str, ...], rows: list[list[float]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        stream.write("% Model,test.mph\n")
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow([f"% {columns[0]}", *columns[1:]])
        writer.writerows(rows)


def _series_rows(
    times: list[float],
    *,
    factor: float,
    drifting_metastable: bool,
) -> tuple[list[list[float]], list[list[float]]]:
    volume_rows: list[list[float]] = []
    coil_rows: list[list[float]] = []
    metastable = [60.0, 65.0, 70.0, 79.0, 88.0, 93.0]
    for index, time_s in enumerate(times):
        small_change = 1.0 + 5.0e-4 * index
        electron = 100.0 * factor * small_change
        mean_energy = 5.0 * factor * small_change
        metastable_value = (
            metastable[index] if drifting_metastable else 80.0 * factor * small_change
        )
        volume_rows.append(
            [
                time_s,
                0.005,
                electron,
                electron * mean_energy,
                metastable_value,
                101.0 * factor * small_change,
                1200.0 * factor * small_change,
            ]
        )
        coil_rows.append([time_s, 1500.0 * small_change])
    return volume_rows, coil_rows


def _write_case(
    root: Path,
    mapping: Path,
    case_id: str,
    *,
    factor: float,
    original: bool = False,
) -> Path:
    directory = root / case_id
    directory.mkdir(parents=True)
    times = (
        [0.0, 2.0e-4, 4.0e-4, 6.0e-4, 8.0e-4, 1.0e-3]
        if original
        else [0.0, 2.0e-4, 5.0e-4, 1.0e-3, 5.0e-3, 1.0e-2]
    )
    volume_rows, coil_rows = _series_rows(
        times,
        factor=factor,
        drifting_metastable=original,
    )
    common_index = times.index(1.0e-3)
    _table(
        directory / "volume_common_time.csv",
        VOLUME_COLUMNS,
        [volume_rows[common_index]],
    )
    _table(
        directory / "coil_power_common_time.csv",
        COIL_COLUMNS,
        [coil_rows[common_index]],
    )
    _table(
        directory / "volume_terminal.csv",
        VOLUME_COLUMNS,
        [volume_rows[-1]],
    )
    _table(
        directory / "coil_power_terminal.csv",
        COIL_COLUMNS,
        [coil_rows[-1]],
    )
    _table(directory / "volume_time_series.csv", VOLUME_COLUMNS, volume_rows)
    _table(directory / "coil_power_time_series.csv", COIL_COLUMNS, coil_rows)
    for name in ("fields_common_time.csv", "fields_terminal.csv"):
        (directory / name).write_text("% field evidence\n0,1\n", encoding="utf-8")
    with (directory / "solution_times.csv").open(
        "w", encoding="utf-8", newline=""
    ) as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(["solution_index", "time_s"])
        writer.writerows((index, value) for index, value in enumerate(times, start=1))

    mph = root / f"{case_id}.mph"
    mph.write_bytes(f"saved-{case_id}".encode())
    outputs = [
        {
            "path": str(directory / name),
            "size_bytes": (directory / name).stat().st_size,
            "sha256": _hash(directory / name),
        }
        for name in EXPECTED_EXPORT_NAMES
    ]
    source_hash = _hash(mph)
    manifest = {
        "schema_version": 1,
        "stage": "execute-gec-icp-comparison-export",
        "status": "completed",
        "case_id": case_id,
        "read_only": True,
        "input_mph": {
            "path": str(mph),
            "sha256": source_hash,
            "post_execution_sha256": source_hash,
        },
        "model_mapping": {
            "path": str(mapping),
            "sha256": _hash(mapping),
        },
        "time_selection": {
            "common": {
                "mode": "exact_transient_interpolation",
                "time_s": 1.0e-3,
            },
            "terminal": {"mode": "last_saved_solution"},
        },
        "expected_outputs": [str(directory / name) for name in EXPECTED_EXPORT_NAMES],
        "outputs": outputs,
    }
    (directory / "comparison_export_manifest.json").write_text(
        json.dumps(manifest), encoding="utf-8"
    )
    return directory


def _four_cases(tmp_path: Path) -> dict[str, Path]:
    mapping = tmp_path / "mapping.yaml"
    mapping.write_text("schema_version: 2\n", encoding="utf-8")
    return {
        "original": _write_case(
            tmp_path, mapping, "original", factor=1.0, original=True
        ),
        "two_term": _write_case(tmp_path, mapping, "two_term", factor=1.2),
        "propagator": _write_case(tmp_path, mapping, "propagator", factor=0.9),
        "composite": _write_case(tmp_path, mapping, "composite", factor=1.1),
    }


def test_comparison_separates_common_and_terminal_and_detects_drift(
    tmp_path: Path,
) -> None:
    pytest.importorskip("matplotlib")
    cases = _four_cases(tmp_path)

    summary = generate_gec_icp_saved_solution_comparison(
        cases=cases,
        output_directory=tmp_path / "analysis",
        case_labels={"composite": "MC + two-term fallback"},
    )

    assert summary.manifest.is_file()
    assert len(summary.figures) == 4
    assert all(path.is_file() and path.stat().st_size > 0 for path in summary.figures)
    manifest = json.loads(summary.manifest.read_text(encoding="utf-8"))
    assert manifest["comparison_contract"]["common_time"]["like_for_like"] is True
    assert manifest["comparison_contract"]["terminal"]["like_for_like_time"] is False
    assert manifest["assessment"]["baseline_terminal_stationarity_passed"] is False
    assert (
        "metastable_inventory"
        in manifest["assessment"]["baseline_failed_stationarity_metrics"]
    )
    assert len(manifest["outputs"]) == 6
    assert len(manifest["cases"]) == 4
    for case in manifest["cases"]:
        assert set(case["source"]["csv_sha256"]) == set(EXPECTED_EXPORT_NAMES)

    with summary.comparison_csv.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 8
    common = [row for row in rows if row["basis"] == "common_time_exact"]
    terminal = [row for row in rows if row["basis"] == "terminal_model_specific"]
    assert all(row["like_for_like_time"] == "1" for row in common)
    assert all(row["like_for_like_time"] == "0" for row in terminal)
    two_term = next(row for row in common if row["case_id"] == "two_term")
    expected_change = (120.0 * 1.0015 - 100.0 * 1.0025) / (100.0 * 1.0025)
    assert float(
        two_term["electron_inventory_signed_change_from_baseline"]
    ) == pytest.approx(expected_change)

    with summary.stationarity_csv.open(encoding="utf-8", newline="") as stream:
        stationarity = list(csv.DictReader(stream))
    assert len(stationarity) == 4 * len(METRIC_KEYS)
    original_metastable = next(
        row
        for row in stationarity
        if row["case_id"] == "original" and row["metric"] == "metastable_inventory"
    )
    assert original_metastable["passed"] == "0"
    assert float(original_metastable["maximum_adjacent_relative_change"]) > 0.03


def test_comparison_rejects_csv_changed_after_export(tmp_path: Path) -> None:
    cases = _four_cases(tmp_path)
    target = cases["two_term"] / "volume_common_time.csv"
    target.write_text(target.read_text(encoding="utf-8") + "\n", encoding="utf-8")

    with pytest.raises(GecIcpComparisonAnalysisError, match="output provenance"):
        generate_gec_icp_saved_solution_comparison(
            cases={"original": cases["original"], "two_term": cases["two_term"]},
            output_directory=tmp_path / "analysis",
        )


def test_comparison_rejects_wrong_comsol_header_even_with_updated_hash(
    tmp_path: Path,
) -> None:
    cases = _four_cases(tmp_path)
    directory = cases["propagator"]
    target = directory / "volume_terminal.csv"
    target.write_text(
        target.read_text(encoding="utf-8").replace(
            "electron inventory (1)", "unexpected electron column"
        ),
        encoding="utf-8",
    )
    manifest_path = directory / "comparison_export_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    record = next(
        item for item in manifest["outputs"] if Path(item["path"]).name == target.name
    )
    record["sha256"] = _hash(target)
    record["size_bytes"] = target.stat().st_size
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(GecIcpComparisonAnalysisError, match="header"):
        generate_gec_icp_saved_solution_comparison(
            cases={
                "original": cases["original"],
                "propagator": cases["propagator"],
            },
            output_directory=tmp_path / "analysis",
        )


def test_comparison_rejects_more_than_four_cases(tmp_path: Path) -> None:
    cases = {f"case_{index}": tmp_path / str(index) for index in range(5)}

    with pytest.raises(GecIcpComparisonAnalysisError, match="2 to 4"):
        generate_gec_icp_saved_solution_comparison(
            cases=cases,
            output_directory=tmp_path / "analysis",
            baseline_case_id="case_0",
        )
