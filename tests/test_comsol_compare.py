from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from swarm_workflow.cli import main as workflow_cli_main
from swarm_workflow.comsol_compare import ComsolCompareError, compare_comsol_profiles


def test_compare_comsol_profiles_writes_minimal_metrics(tmp_path: Path) -> None:
    external = _write_profile(tmp_path / "external.csv", scale=0.9)
    reference = _write_profile(tmp_path / "reference.csv", scale=1.0)

    summary = compare_comsol_profiles(
        external,
        reference,
        output_dir=tmp_path / "comparison",
        external_runtime_s=10.0,
        reference_runtime_s=40.0,
    )

    assert summary.speedup == pytest.approx(4.0)
    assert summary.quantities == ("electron_density", "E_over_N")
    rows = list(csv.DictReader(summary.metrics_csv.open(encoding="utf-8")))
    assert {row["quantity"] for row in rows} == {"electron_density", "E_over_N"}
    density = next(row for row in rows if row["quantity"] == "electron_density")
    assert float(density["l2_relative_error"]) == pytest.approx(0.1)
    payload = json.loads(summary.summary_json.read_text(encoding="utf-8"))
    assert payload["runtime"]["speedup_reference_over_external"] == 4.0
    assert payload["operating_conditions"] == {
        "applied_voltage": 200.0,
        "gas_pressure": 100.0,
    }
    assert payload["current_conservation"] is None


def test_compare_comsol_profiles_rejects_different_voltage(tmp_path: Path) -> None:
    external = _write_profile(tmp_path / "external.csv", scale=1.0)
    reference = _write_profile(tmp_path / "reference.csv", scale=1.0)
    reference.write_text(
        reference.read_text(encoding="utf-8").replace(",200.0,100.0", ",100.0,100.0"),
        encoding="utf-8",
    )

    with pytest.raises(ComsolCompareError, match="applied_voltage.*differs"):
        compare_comsol_profiles(external, reference, output_dir=tmp_path / "out")


def test_compare_comsol_profiles_rejects_nonfinite_values(tmp_path: Path) -> None:
    external = _write_profile(tmp_path / "external.csv", scale=1.0)
    reference = _write_profile(tmp_path / "reference.csv", scale=1.0)
    external.write_text(
        external.read_text(encoding="utf-8").replace("20.0,", "nan,"),
        encoding="utf-8",
    )

    with pytest.raises(ComsolCompareError, match="non-finite"):
        compare_comsol_profiles(external, reference, output_dir=tmp_path / "out")


def test_compare_comsol_profiles_reports_bulk_current_conservation(
    tmp_path: Path,
) -> None:
    content = (
        "x,total_current_density,applied_voltage,gas_pressure\n"
        "0.0,-10.0,200.0,100.0\n"
        "0.5,-1.0,200.0,100.0\n"
        "1.0,5.0,200.0,100.0\n"
    )
    external = tmp_path / "external.csv"
    reference = tmp_path / "reference.csv"
    external.write_text(content, encoding="utf-8")
    reference.write_text(content, encoding="utf-8")

    summary = compare_comsol_profiles(
        external,
        reference,
        output_dir=tmp_path / "comparison",
    )

    payload = json.loads(summary.summary_json.read_text(encoding="utf-8"))
    current = payload["current_conservation"]["external"]
    assert current["full_domain"]["min_A_per_m2"] == -10.0
    assert current["bulk_central_80_percent"]["mean_A_per_m2"] == -1.0
    assert current["bulk_central_80_percent"]["relative_standard_deviation"] == 0.0


def test_compare_comsol_cli(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    pytest.importorskip("matplotlib")
    external = _write_profile(tmp_path / "external.csv", scale=1.0)
    reference = _write_profile(tmp_path / "reference.csv", scale=1.0)

    workflow_cli_main(
        [
            "compare-comsol",
            str(external),
            str(reference),
            "--output",
            str(tmp_path / "out"),
            "--external-runtime",
            "5",
            "--reference-runtime",
            "20",
            "--plot",
        ]
    )

    output = capsys.readouterr().out
    assert "speedup: 4x" in output
    assert "spatial_profile_comparison.png" in output
    assert (tmp_path / "out" / "spatial_profile_comparison.png").is_file()


def _write_profile(path: Path, *, scale: float) -> Path:
    path.write_text(
        "x,electron_density,E_over_N,applied_voltage,gas_pressure\n"
        f"0.0,{10.0 * scale},{100.0 * scale},200.0,100.0\n"
        f"0.5,{20.0 * scale},{200.0 * scale},200.0,100.0\n"
        f"1.0,{30.0 * scale},{300.0 * scale},200.0,100.0\n",
        encoding="utf-8",
    )
    return path
