"""Unit-bearing saved results, explicit selection and mesh-independent metrics."""

import csv
from hashlib import sha256
import json
import math

import pytest

from swarm_workflow.cli import main
from swarm_workflow.comsol.models.positive_column.comparison import (
    ComsolCompareError,
    compare_comsol_profiles,
)
from swarm_workflow.comsol.models.positive_column.results import (
    RESULT_SCHEMA,
    read_positive_column_results,
)
from swarm_workflow.comsol.models.positive_column.result_units import QUANTITIES


def _artifact(
    tmp_path,
    name,
    *,
    x=(0.0, 0.5, 1.0),
    values=None,
    units=None,
    time=None,
    solutions=(1,),
):
    directory = tmp_path / name
    directory.mkdir()
    values = values or {
        "electron_density": [10.0, 20.0, 30.0],
        "E_over_N": [1e-20, 2e-20, 3e-20],
    }
    columns = ["x", *values, "applied_voltage", "gas_pressure"]
    units = {name: QUANTITIES[name].unit for name in columns} | (units or {})
    rows = []
    for solution in solutions:
        for point, position in enumerate(x):
            rows.append(
                {
                    "dataset": "dset1",
                    "solution_index": solution,
                    "point_index": point,
                    "time_s": (time[solution] if isinstance(time, dict) else time)
                    if time is not None
                    else "",
                    "parameter_name": "t" if time is not None else "",
                    "parameter_unit": "s" if time is not None else "",
                    "parameter_value": (
                        time[solution] if isinstance(time, dict) else time
                    )
                    if time is not None
                    else "",
                    "x": position,
                    "applied_voltage": 200.0,
                    "gas_pressure": 100.0,
                    **{name: values[point] for name, values in values.items()},
                }
            )
    data = directory / "positive_column_results.csv"
    with data.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metadata = {
        "artifact_schema": RESULT_SCHEMA,
        "status": "completed",
        "result_csv": data.name,
        "result_csv_sha256": sha256(data.read_bytes()).hexdigest(),
        "voltage_parameter": "V0",
        "result_selection": {
            "dataset": "dset1",
            "solutions": "all",
            "field_units": units,
        },
        "physical_context": {
            "gas_temperature_K": 293.15,
            "gas_mixture": [{"species": "Ar", "fraction": 1.0, "mass_amu": 39.948}],
        },
    }
    summary = directory / "positive_column_run_summary.json"
    summary.write_text(json.dumps(metadata))
    return summary


def _modify(summary, update):
    metadata = json.loads(summary.read_text())
    path = summary.parent / metadata["result_csv"]
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    update(rows)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metadata["result_csv_sha256"] = sha256(path.read_bytes()).hexdigest()
    summary.write_text(json.dumps(metadata))


@pytest.mark.parametrize(
    "nodes", [(0.0, 0.5, 1.0), (0.0, 0.01, 0.02, 0.1, 1.0), (0.0, 0.9, 0.99, 1.0)]
)
@pytest.mark.parametrize("scale", [1.0, 1e-200, 1e200])
def test_spatial_l2_is_independent_of_sampling_density(tmp_path, nodes, scale):
    external = _artifact(
        tmp_path,
        "external",
        x=nodes,
        values={"electron_density": [(1 + x) * scale for x in nodes]},
    )
    reference = _artifact(
        tmp_path, "reference", x=(0.0, 1.0), values={"electron_density": [scale, scale]}
    )
    summary = compare_comsol_profiles(
        external, reference, output_dir=tmp_path / "comparison"
    )
    [row] = list(csv.DictReader(summary.metrics_csv.open()))
    assert float(row["l2_relative_error"]) == pytest.approx(1 / math.sqrt(3))
    assert float(row["integral_ratio_external_over_reference"]) == pytest.approx(1.5)
    assert row["shape_correlation"] == ""
    assert row["unit"] == "1/m^3"
    payload = json.loads(summary.summary_json.read_text())
    assert payload["common_x_range_m"] == [0.0, 1.0]
    assert payload["sources"]["external"]["selection"]["solution_index"] == 1
    assert summary.profiles_csv.is_file()


def test_units_are_normalized_before_comparison_and_curve_export(tmp_path):
    a = _artifact(
        tmp_path,
        "a",
        x=(0.0, 50.0, 100.0),
        values={"electron_density": [1.0, 2.0, 3.0], "E_over_N": [10.0, 20.0, 30.0]},
        units={
            "x": "cm",
            "electron_density": "1/cm^3",
            "E_over_N": "Td",
            "applied_voltage": "kV",
            "gas_pressure": "kPa",
        },
    )
    _modify(
        a,
        lambda rows: [
            row.update(applied_voltage=0.2, gas_pressure=0.1) for row in rows
        ],
    )
    b = _artifact(
        tmp_path,
        "b",
        values={"electron_density": [1e6, 2e6, 3e6], "E_over_N": [1e-20, 2e-20, 3e-20]},
    )
    summary = compare_comsol_profiles(a, b, output_dir=tmp_path / "comparison")
    metrics = list(csv.DictReader(summary.metrics_csv.open()))
    assert all(float(row["l2_relative_error"]) < 1e-14 for row in metrics)
    profiles = list(csv.DictReader(summary.profiles_csv.open()))
    assert max(float(row["x_m"]) for row in profiles) == 1
    density = [row for row in profiles if row["quantity"] == "electron_density"]
    assert [float(row["external"]) for row in density] == pytest.approx([1e6, 2e6, 3e6])
    artifact = read_positive_column_results(a)
    assert artifact.units["electron_density"] == "1/m^3"
    assert artifact.metadata["source_units"]["electron_density"] == "1/cm^3"


@pytest.mark.parametrize(
    "fault,reason",
    [
        ("voltage", "applied_voltage.*differs"),
        ("time", "time_s.*differs"),
        ("mixed_time", "time_s.*constant"),
        ("dataset", "dataset differs"),
        ("duplicate", "duplicate solution/point"),
        ("nan", "non-finite"),
        ("unit", "unsupported"),
        ("physical", "physical_context differs"),
        ("negative_density", "electron_density.*nonnegative"),
        ("negative_source", "ionization_source.*nonnegative"),
        ("selection", "missing units or saved-solution selection"),
    ],
)
def test_incompatible_or_ambiguous_sources_fail_before_publication(
    tmp_path, fault, reason
):
    a = _artifact(
        tmp_path, "a", time=0.001 if fault in {"time", "mixed_time"} else None
    )
    b = _artifact(
        tmp_path,
        "b",
        time=0.002 if fault == "time" else 0.001 if fault == "mixed_time" else None,
    )
    if fault == "voltage":
        _modify(a, lambda rows: [row.update(applied_voltage=100) for row in rows])
    if fault == "mixed_time":
        _modify(a, lambda rows: rows[1].update(time_s=0.002))
    if fault == "dataset":
        _modify(a, lambda rows: rows[0].update(dataset="other"))
    if fault == "duplicate":
        _modify(a, lambda rows: rows[1].update(point_index=0))
    if fault == "nan":
        _modify(a, lambda rows: rows[0].update(electron_density="nan"))
    if fault == "negative_density":
        _modify(a, lambda rows: rows[0].update(electron_density=-1))
    if fault == "negative_source":
        a = _artifact(tmp_path, "source", values={"ionization_source": [-1, 1, 2]})
    if fault == "selection":
        data = json.loads(a.read_text())
        del data["result_selection"]["field_units"]
        a.write_text(json.dumps(data))
    if fault in {"unit", "physical"}:
        data = json.loads(a.read_text())
        if fault == "unit":
            data["result_selection"]["field_units"]["electron_density"] = "m^3"
        else:
            data["physical_context"]["gas_temperature_K"] = 300.0
        a.write_text(json.dumps(data))
    with pytest.raises(ComsolCompareError, match=reason):
        compare_comsol_profiles(a, b, output_dir=tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_explicit_solution_selection_and_complete_common_interval(tmp_path):
    a = _artifact(
        tmp_path, "a", x=(0.0, 0.4, 1.0), solutions=(2, 3), time={2: 0.001, 3: 0.002}
    )
    b = _artifact(
        tmp_path, "b", x=(0.1, 0.5, 0.9), solutions=(8, 9), time={8: 0.001, 9: 0.002}
    )
    with pytest.raises(ComsolCompareError, match="one saved solution"):
        compare_comsol_profiles(a, b, output_dir=tmp_path / "missing")
    summary = compare_comsol_profiles(
        a, b, external_solution=3, reference_solution=9, output_dir=tmp_path / "out"
    )
    payload = json.loads(summary.summary_json.read_text())
    assert payload["common_x_range_m"] == [0.1, 0.9]
    assert payload["sources"]["external"]["selection"]["solution_index"] == 3
    x = {float(row["x_m"]) for row in csv.DictReader(summary.profiles_csv.open())}
    assert x == {0.1, 0.4, 0.5, 0.9}
    with pytest.raises(ComsolCompareError, match="unavailable"):
        compare_comsol_profiles(
            a,
            b,
            external_solution=99,
            reference_solution=8,
            output_dir=tmp_path / "bad",
        )


def test_current_statistics_integrate_the_spatial_profile(tmp_path):
    a = _artifact(tmp_path, "a", values={"total_current_density": [-10.0, -1.0, 5.0]})
    b = _artifact(tmp_path, "b", values={"total_current_density": [-10.0, -1.0, 5.0]})
    summary = compare_comsol_profiles(a, b, output_dir=tmp_path / "out")
    current = json.loads(summary.summary_json.read_text())[
        "conductive_current_uniformity"
    ]["external"]
    assert current["full_domain"]["mean_A_per_m2"] == pytest.approx(-1.75)
    assert current["central_80_percent"]["mean_A_per_m2"] == pytest.approx(-1.6)
    assert current["central_80_percent"][
        "standard_deviation_A_per_m2"
    ] == pytest.approx(math.sqrt(12.12))


def test_comparison_cli_outputs_native_data_png_and_svg_without_overwriting_sources(
    tmp_path, capsys
):
    pytest.importorskip("matplotlib")
    a, b = (_artifact(tmp_path, name) for name in ("a", "b"))
    before = a.read_bytes()
    main(
        [
            "compare-positive-column",
            str(a),
            str(b),
            "--output",
            str(tmp_path / "out"),
            "--plot",
            "--external-runtime",
            "5",
            "--reference-runtime",
            "20",
        ]
    )
    output = capsys.readouterr().out
    assert "speedup: 4x" in output
    assert (
        (tmp_path / "out/spatial_profile_comparison.png")
        .read_bytes()
        .startswith(b"\x89PNG")
    )
    assert (tmp_path / "out/spatial_profile_comparison.svg").is_file()
    assert a.read_bytes() == before
    with pytest.raises(ComsolCompareError, match="already exists"):
        compare_comsol_profiles(a, b, output_dir=tmp_path / "out")


def test_unitless_or_altered_results_do_not_acquire_success(tmp_path):
    a = _artifact(tmp_path, "a")
    b = _artifact(tmp_path, "b")
    csv_path = a.parent / "positive_column_results.csv"
    with pytest.raises(ComsolCompareError, match="bare CSV"):
        compare_comsol_profiles(csv_path, b, output_dir=tmp_path / "bare")
    csv_path.write_text(csv_path.read_text().replace("10.0", "11.0"))
    with pytest.raises(ComsolCompareError, match="altered"):
        compare_comsol_profiles(a, b, output_dir=tmp_path / "altered")


def test_zero_reference_and_two_point_current_keep_undefined_metrics_missing(tmp_path):
    a = _artifact(tmp_path, "a", x=(0.0, 1.0), values={"total_current_density": [0, 0]})
    b = _artifact(tmp_path, "b", x=(0.0, 1.0), values={"total_current_density": [0, 0]})
    summary = compare_comsol_profiles(a, b, output_dir=tmp_path / "out")
    [row] = list(csv.DictReader(summary.metrics_csv.open()))
    assert row["l2_relative_error"] == row["shape_correlation"] == ""
    assert row["integral_ratio_external_over_reference"] == ""
    current = json.loads(summary.summary_json.read_text())[
        "conductive_current_uniformity"
    ]["external"]
    assert current["central_80_percent"]["mean_A_per_m2"] == 0
    assert current["central_80_percent"]["standard_deviation_A_per_m2"] == 0
    assert current["central_80_percent"]["relative_standard_deviation"] is None


def test_failed_plot_does_not_publish_a_partial_comparison(tmp_path, monkeypatch):
    a, b = (_artifact(tmp_path, name) for name in ("a", "b"))

    def fail_plot(directory, *args):
        (directory / "spatial_profile_comparison.png").write_bytes(b"incomplete")
        raise RuntimeError("plot interrupted")

    monkeypatch.setattr(
        "swarm_workflow.comsol.models.positive_column.comparison._plot_profiles",
        fail_plot,
    )
    with pytest.raises(RuntimeError, match="plot interrupted"):
        compare_comsol_profiles(a, b, output_dir=tmp_path / "out", create_plot=True)
    assert not (tmp_path / "out").exists()
    assert not list(tmp_path.glob(".positive-column-comparison-*"))
    assert read_positive_column_results(a).metadata["status"] == "completed"


@pytest.mark.parametrize("times", [(1.0, None), (-1.0, 2.0), (1e-300, 1e300)])
def test_invalid_runtime_pairs_do_not_publish_speedup(tmp_path, times):
    a, b = (_artifact(tmp_path, name) for name in ("a", "b"))
    with pytest.raises(
        ComsolCompareError,
        match="both solve runtimes|positive and finite|floating-point range",
    ):
        compare_comsol_profiles(
            a,
            b,
            output_dir=tmp_path / "out",
            external_runtime_s=times[0],
            reference_runtime_s=times[1],
        )
    assert not (tmp_path / "out").exists()
