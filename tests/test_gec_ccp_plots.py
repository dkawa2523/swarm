from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path

import numpy as np
import pytest

from swarm_workflow.quality.policy import QualityThresholds, quality_thresholds_payload
from swarm_workflow.comsol.models.gec_ccp.plots.eedf import (
    _function_eedf_slice,
    _weighted_quantiles,
)
from swarm_workflow.comsol.models.gec_ccp.plots.data import (
    _comsol_field_matrix,
    _comsol_radial_midplane_data,
    _comsol_sorted_phase_field,
    _comsol_spatial_data,
    _comsol_waveform_data,
    _spatial_comparison_metrics,
    _validated_plot_provenance,
    _waveform_comparison_metrics,
)
from swarm_workflow.comsol.models.gec_ccp.plots.workflow import (
    GecCcpPlotError,
    plot_gec_ccp_solver_comparison,
)
from swarm_workflow.comsol.models.gec_ccp.plots.spatial import (
    _gec_domain_geometry,
    _gec_domain_polygon,
    _plot_domain_2d_comparison,
)
from swarm_workflow.comsol.input.function_eedf import ComsolFunctionEedfGrid


def test_waveform_metrics_integrate_periodic_endpoints_on_each_native_grid(
    tmp_path: Path,
) -> None:
    paths = [tmp_path / "baseline.csv", tmp_path / "external.csv"]
    for path, samples, offset, amplitude in (
        (paths[0], 31, -70.0, 0.1),
        (paths[1], 41, -50.0, 0.2),
    ):
        phase = np.linspace(0.0, 1.0, samples)
        np.savetxt(
            path,
            np.column_stack((
                phase,
                offset + 100.0 * np.cos(2.0 * np.pi * phase),
                amplitude * np.sin(2.0 * np.pi * phase),
            )),
            delimiter=",",
            header="x1_ptp,ptp.mct1.V (V),ptp.mct1.I (A)",
            comments="% ",
        )
    metrics = _waveform_comparison_metrics(*paths)
    assert metrics["averaging_method"] == "trapezoidal_over_each_exported_rf_period"
    assert metrics["voltage_mean_V"] == pytest.approx({"baseline": -70.0, "external": -50.0})
    assert metrics["current_rms_A"]["baseline"] == pytest.approx(0.1 / np.sqrt(2.0))
    assert metrics["current_rms_A"]["external"] == pytest.approx(0.2 / np.sqrt(2.0))
    assert metrics["current_rms_A"]["external_to_baseline_ratio"] == pytest.approx(2.0)


def test_waveform_metrics_weight_nonuniform_phases_and_preserve_native_peaks(
    tmp_path: Path,
) -> None:
    paths = [tmp_path / "baseline.csv", tmp_path / "external.csv"]
    for path, phase, voltage, current in (
        (paths[0], [1.0, 0.1, 0.0, 0.7, 0.2], [0.0, 5.0, 0.0, 3.75, 10.0], 2.0),
        (paths[1], [0.0, 0.4, 0.8, 1.0], [0.0, 30.0, 10.0, 0.0], 6.0),
    ):
        np.savetxt(
            path,
            np.column_stack((phase, voltage, np.full(len(phase), current))),
            delimiter=",",
            header="x1_ptp,ptp.mct1.V (V),ptp.mct1.I (A)",
            comments="% ",
        )
    metrics = _waveform_comparison_metrics(*paths)
    # The voltage is a triangular pulse: its period mean is half its peak.
    assert metrics["voltage_mean_V"] == pytest.approx({"baseline": 5.0, "external": 15.0})
    assert metrics["voltage_peak_to_peak_V"]["external"] == pytest.approx(30.0)
    assert metrics["current_rms_A"]["external"] == pytest.approx(6.0)


@pytest.mark.parametrize("phase", [[0.0, 0.0], [0.0, 0.5, 0.5, 1.0]])
def test_waveform_metrics_reject_duplicate_phases(tmp_path: Path, phase: list[float]) -> None:
    path = tmp_path / "ambiguous.csv"
    np.savetxt(
        path,
        np.column_stack((phase, np.ones(len(phase)), np.ones(len(phase)))),
        delimiter=",",
        header="x1_ptp,ptp.mct1.V (V),ptp.mct1.I (A)",
        comments="% ",
    )
    with pytest.raises(GecCcpPlotError, match="distinct RF phases"):
        _waveform_comparison_metrics(path, path)


def test_plot_provenance_binds_bundle_hashes_and_closure_axes(
    tmp_path: Path,
) -> None:
    bundle = tmp_path / "bundle"
    results = tmp_path / "results"
    bundle.mkdir()
    (results / "swarm_tables").mkdir(parents=True)
    hashes = {
        "workflow_config_sha256": "c" * 64,
        "base_config_sha256": "a" * 64,
        "cross_sections_sha256": "b" * 64,
    }
    quality_thresholds = quality_thresholds_payload(QualityThresholds())
    bundle_manifest = {
        "source": "two_term",
        "hashes": hashes,
        "quality_thresholds": quality_thresholds,
    }
    bundle_manifest_path = bundle / "manifest.json"
    bundle_manifest_path.write_text(
        json.dumps(bundle_manifest), encoding="utf-8"
    )
    bundle_manifest_sha256 = sha256(
        bundle_manifest_path.read_bytes()
    ).hexdigest()
    input_mph_sha256 = "d" * 64
    java_hashes = {"apply": "e" * 64}
    output_mph = results / "external.mph"
    output_mph.write_bytes(b"saved mph")
    result_file = results / "swarm_tables" / "domain_period_average.csv"
    result_file.write_text("% r,z,ne\n0,0,1\n", encoding="utf-8")
    plan_path = results / "gec_ccp_plan.json"
    plan_path.write_text(
        json.dumps(
            {
                "status": "ready",
                "model": {"input_mph": {"sha256": input_mph_sha256}},
                "generated_java": {
                    "apply": {"path": "apply.java", "sha256": java_hashes["apply"]}
                },
                "bundle": {
                    "path": str(bundle.resolve()),
                    "source": "two_term",
                    "hashes": hashes,
                    "quality_thresholds": quality_thresholds,
                    "manifest_sha256": bundle_manifest_sha256,
                },
                "closure": {
                    "electron_transport": "comsol_specify_all_restricted",
                    "reaction_model": "function_eedf",
                },
                "expected_results": [str(result_file)],
            }
        ),
        encoding="utf-8",
    )
    (results / "gec_ccp_run_status.json").write_text(
        json.dumps(
            {
                "status": "completed",
                "solve_status": "completed",
                "quality_status": "passed",
                "quality_accepted": True,
                "plan_sha256": sha256(plan_path.read_bytes()).hexdigest(),
                "input_mph_sha256": input_mph_sha256,
                "generated_java_sha256": java_hashes,
                "bundle_artifacts_verified": True,
                "comsol_runtime": {
                    "version": "6.4",
                    "build": "429",
                    "executable": "comsolbatch.exe",
                },
                "output_mph": {
                    "path": str(output_mph),
                    "sha256": sha256(output_mph.read_bytes()).hexdigest(),
                },
                "results": _result_metadata(results, [result_file]),
            }
        ),
        encoding="utf-8",
    )

    provenance = _validated_plot_provenance(
        bundle.resolve(), bundle_manifest, results.resolve()
    )

    assert provenance["quality_status"] == "passed"
    assert "two-term Swarm" in provenance["external_label"]
    assert "Function-EEDF" in provenance["external_label"]
    assert provenance["results"][0]["path"] == (
        "swarm_tables/domain_period_average.csv"
    )

    solver_specific = _validated_plot_provenance(
        bundle.resolve(),
        {
            **bundle_manifest,
            "hashes": {
                **hashes,
                "mc_sampling_plan_json": "[]",
                "mc_sampling_plan_sha256": "f" * 64,
                "mc_transport_estimator_schema_version": (
                    "direct_mc_transport.v3"
                ),
            },
        },
        results.resolve(),
    )
    assert solver_specific["quality_status"] == "passed"

    with pytest.raises(GecCcpPlotError, match="hashes differ"):
        _validated_plot_provenance(
            bundle.resolve(),
            {
                "source": "two_term",
                "hashes": {**hashes, "extra": "bad"},
                "quality_thresholds": quality_thresholds,
            },
            results.resolve(),
        )

    status_path = results / "gec_ccp_run_status.json"
    status = json.loads(status_path.read_text(encoding="utf-8"))
    status.update(
        {
            "status": "rejected",
            "quality_status": "failed",
            "quality_accepted": False,
        }
    )
    status_path.write_text(json.dumps(status), encoding="utf-8")
    with pytest.raises(GecCcpPlotError, match="physics-quality acceptance"):
        _validated_plot_provenance(
            bundle.resolve(), bundle_manifest, results.resolve()
        )

    status.update(
        {
            "status": "completed",
            "quality_status": "passed",
            "quality_accepted": True,
        }
    )
    status["plan_sha256"] = "0" * 64
    status_path.write_text(json.dumps(status), encoding="utf-8")
    with pytest.raises(GecCcpPlotError, match="does not belong"):
        _validated_plot_provenance(
            bundle.resolve(), bundle_manifest, results.resolve()
        )


def test_plot_provenance_rejects_tampered_result_csv(tmp_path: Path) -> None:
    bundle, results = _write_comparison_case(
        tmp_path, "two_term", external_scale=1.1
    )
    manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    _validated_plot_provenance(bundle.resolve(), manifest, results.resolve())

    result = results / "swarm_tables" / "domain_period_average.csv"
    result.write_text(
        result.read_text(encoding="utf-8") + "0,0,999,999,999,999\n",
        encoding="utf-8",
    )

    with pytest.raises(GecCcpPlotError, match="changed since the solve"):
        _validated_plot_provenance(bundle.resolve(), manifest, results.resolve())


def test_plot_provenance_rejects_missing_result_csv(tmp_path: Path) -> None:
    bundle, results = _write_comparison_case(
        tmp_path, "two_term", external_scale=1.1
    )
    manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    (results / "swarm_tables" / "domain_period_average.csv").unlink()

    with pytest.raises(GecCcpPlotError, match="missing or unreadable"):
        _validated_plot_provenance(bundle.resolve(), manifest, results.resolve())


def test_plot_provenance_rejects_legacy_path_only_result_status(
    tmp_path: Path,
) -> None:
    bundle, results = _write_comparison_case(
        tmp_path, "two_term", external_scale=1.1
    )
    manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    status_path = results / "gec_ccp_run_status.json"
    status = json.loads(status_path.read_text(encoding="utf-8"))
    status["results"] = [item["path"] for item in status["results"]]
    status_path.write_text(json.dumps(status), encoding="utf-8")

    with pytest.raises(GecCcpPlotError, match="provenance entry is invalid"):
        _validated_plot_provenance(bundle.resolve(), manifest, results.resolve())


def test_comsol_spatial_data_uses_varying_coordinate_after_2d_coordinates(
    tmp_path: Path,
) -> None:
    exported = tmp_path / "axis.csv"
    exported.write_text(
        "\n".join(
            [
                "% Dimension,2",
                "% R,Z,ne,Te",
                "0,0.02,2,4",
                "0,0.00,1,3",
                "0,0.01,1.5,3.5",
            ]
        ),
        encoding="utf-8",
    )

    coordinate, fields = _comsol_spatial_data(exported)

    np.testing.assert_allclose(coordinate, [0.0, 0.01, 0.02])
    np.testing.assert_allclose(
        fields,
        [[1.0, 3.0], [1.5, 3.5], [2.0, 4.0]],
    )

    metrics = _spatial_comparison_metrics(exported, exported)
    density = metrics["fields"]["electron_density_m3"]
    assert density["external_to_baseline_abs_peak_ratio"] == 1.0
    assert density["relative_l2_difference"] == 0.0


def test_comsol_phase_fields_and_density_weighted_quantiles(
    tmp_path: Path,
) -> None:
    exported = tmp_path / "phase.csv"
    exported.write_text(
        "\n".join(
            [
                "% Dimension,2",
                "% R,Z,ptp.ebar (V) @ t=0,ptp.ebar (V) @ t=1",
                "0,0.02,2,4",
                "0,0.00,3,5",
            ]
        ),
        encoding="utf-8",
    )

    fields = _comsol_field_matrix(exported, "ptp.ebar")
    np.testing.assert_allclose(fields, [[2.0, 4.0], [3.0, 5.0]])
    coordinate, sorted_fields = _comsol_sorted_phase_field(
        exported, "ptp.ebar"
    )
    np.testing.assert_allclose(coordinate, [0.0, 0.02])
    np.testing.assert_allclose(sorted_fields, [[3.0, 5.0], [2.0, 4.0]])
    quantiles = _weighted_quantiles(
        np.asarray([2.0, 4.0, 6.0]),
        np.asarray([1.0, 8.0, 1.0]),
        np.asarray([0.5]),
    )
    assert 3.5 < quantiles[0] < 4.5


def test_comsol_waveform_data_accepts_row_wise_and_native_wide_exports(
    tmp_path: Path,
) -> None:
    row_wise = tmp_path / "row_wise.csv"
    row_wise.write_text(
        "\n".join(
            [
                "% Dimension,1",
                "% x1_ptp,ptp.mct1.V (V),ptp.mct1.I (A)",
                "0,10,0.1",
                "0.5,-20,-0.2",
                "1,10,0.1",
            ]
        ),
        encoding="utf-8",
    )
    wide = tmp_path / "wide.csv"
    wide.write_text(
        "\n".join(
            [
                "% Dimension,2",
                "% R,Z,ptp.mct1.V (V) @ t=0,ptp.mct1.I (A) @ t=0,"
                "ptp.mct1.V (V) @ t=5E-9,ptp.mct1.I (A) @ t=5E-9,"
                "ptp.mct1.V (V) @ t=1E-8,ptp.mct1.I (A) @ t=1E-8",
                "0,0,10,0.1,-20,-0.2,10,0.1",
                "0.1,0.2,10,0.1,-20,-0.2,10,0.1",
            ]
        ),
        encoding="utf-8",
    )

    for exported in (row_wise, wide):
        phase, voltage, current = _comsol_waveform_data(exported)
        np.testing.assert_allclose(phase, [0.0, 0.5, 1.0])
        np.testing.assert_allclose(voltage, [10.0, -20.0, 10.0])
        np.testing.assert_allclose(current, [0.1, -0.2, 0.1])


def test_comsol_waveform_data_rejects_non_global_wide_terminal_value(
    tmp_path: Path,
) -> None:
    exported = tmp_path / "wide.csv"
    exported.write_text(
        "\n".join(
            [
                "% Dimension,2",
                "% R,Z,ptp.mct1.V (V) @ t=0,ptp.mct1.I (A) @ t=0,"
                "ptp.mct1.V (V) @ t=1,ptp.mct1.I (A) @ t=1",
                "0,0,10,0.1,-20,-0.2",
                "0.1,0.2,11,0.1,-20,-0.2",
            ]
        ),
        encoding="utf-8",
    )

    with pytest.raises(GecCcpPlotError, match="voltage is not global"):
        _comsol_waveform_data(exported)


def test_operating_eedf_slice_uses_active_physical_grid_and_linear_mean() -> None:
    grid = ComsolFunctionEedfGrid(
        electron_energies_eV=np.asarray([0.0, 1.0, 4.0]),
        mean_energies_eV=np.asarray([1.0, 3.0]),
        values_eV_m32=np.asarray(
            [[4.0, 2.0, 1.0], [8.0, 6.0, 3.0]]
        ),
        normalization_error_max=0.0,
        mean_energy_relative_error_max=0.0,
        minimum_value=1.0,
    )

    plotted = _function_eedf_slice(grid, 2.0)

    np.testing.assert_array_equal(plotted["energy_eV"], [0.0, 1.0, 4.0])
    np.testing.assert_allclose(
        plotted["eedf_eV_inv"],
        np.sqrt([0.0, 1.0, 4.0]) * np.asarray([6.0, 4.0, 2.0]),
    )
    assert plotted["source"] == "active_comsol_function_eedf_physical_grid"


def test_gec_domain_geometry_is_inferred_from_stepped_mesh() -> None:
    coordinates = np.asarray(
        [
            [0.0, 0.0],
            [0.0, 0.0254],
            [0.0508, 0.0],
            [0.0538, -0.0381],
            [0.0538, 0.0635],
            [0.1016, -0.0381],
            [0.1016, 0.0635],
        ]
    )

    geometry = _gec_domain_geometry(coordinates)

    assert geometry["gap_z_min_cm"] == 0.0
    assert geometry["gap_z_max_cm"] == 2.54
    assert geometry["powered_radius_cm"] == 5.08
    assert geometry["throat_radius_cm"] == 5.38
    assert geometry["chamber_radius_cm"] == 10.16
    polygon = _gec_domain_polygon(geometry)
    assert (0.0, 0.0) in polygon
    assert (10.16, 6.35) in polygon


def test_domain_surface_comparison_writes_only_four_png_figures(
    tmp_path: Path,
) -> None:
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg", force=True)
    plt = pytest.importorskip("matplotlib.pyplot")
    coordinates = (
        (0.0000, 0.0000),
        (0.0000, 0.0254),
        (0.0538, 0.0000),
        (0.0538, 0.0254),
        (0.0538, -0.0381),
        (0.0538, 0.0635),
        (0.1016, -0.0381),
        (0.1016, 0.0635),
    )

    def write_domain(path: Path, scale: float) -> None:
        rows = ["% Dimension,2", "% r,z,ne,Te,V,Re"]
        for index, (radius, axial) in enumerate(coordinates, start=1):
            rows.append(
                f"{radius},{axial},{scale * index:.6g},"
                f"{2.0 + 0.1 * index:.6g},{-10.0 + index:.6g},"
                f"{scale * index * 1.0e5:.6g}"
            )
        path.write_text("\n".join(rows) + "\n", encoding="utf-8")

    baseline = tmp_path / "baseline.csv"
    external = tmp_path / "external.csv"
    write_domain(baseline, 1.0)
    write_domain(external, 1.1)

    summary, figures = _plot_domain_2d_comparison(
        baseline,
        external,
        tmp_path / "domain_2d_period_average_comparison",
        plt,
    )

    assert len(figures) == 4
    assert all(path.suffix == ".png" and path.exists() for path in figures)
    assert summary["presentation"]["color_scale"] == "linear"
    assert all(
        field["display_scale"] == "linear"
        for field in summary["fields"].values()
    )
    assert summary["fields"]["electron_density_m3"][
        "shared_display_range"
    ][0] == 0.0
    assert summary["fields"]["ionization_source_m3_s"][
        "shared_display_range"
    ][0] == 0.0
    assert not (tmp_path / "domain_2d_period_average_comparison.png").exists()
    assert not list(tmp_path.glob("comsol_style_*_comparison.svg"))


def test_two_solver_comparison_is_accepted_and_fail_closed(
    tmp_path: Path,
) -> None:
    pytest.importorskip("matplotlib")
    two_term_bundle, two_term_results = _write_comparison_case(
        tmp_path, "two_term", external_scale=1.1
    )
    mc_bundle, mc_results = _write_comparison_case(
        tmp_path, "monte_carlo", external_scale=0.9
    )

    summary = plot_gec_ccp_solver_comparison(
        two_term_bundle,
        mc_bundle,
        two_term_results_dir=two_term_results,
        monte_carlo_results_dir=mc_results,
        output_dir=tmp_path / "comparison",
    )

    manifest = json.loads(summary.manifest.read_text(encoding="utf-8"))
    assert len(summary.figures) == 12
    assert sum(path.suffix == ".png" for path in summary.figures) == 6
    assert sum(path.suffix == ".svg" for path in summary.figures) == 6
    assert all(path.exists() for path in summary.figures)
    assert manifest["stage"] == "plot-gec-ccp-comparison"
    assert manifest["run_provenance"]["two_term"]["quality_status"] == "passed"
    assert manifest["run_provenance"]["monte_carlo"]["quality_status"] == "passed"
    assert manifest["run_provenance"]["two_term"]["closure"][
        "reaction_model"
    ] == "function_eedf"
    assert manifest["run_provenance"]["two_term"]["closure"][
        "electron_transport"
    ] == "swarm_mobility_einstein"
    assert "Function-EEDF" in manifest["comparison_labels"]["two_term"]
    assert "preintegrated Swarm rates" in manifest["comparison_labels"][
        "monte_carlo"
    ]
    assert "COMSOL Einstein diffusion" in manifest["comparison_labels"][
        "two_term"
    ]
    assert manifest["comparison_metrics"]["reference_solver"] == "two_term"
    assert manifest["comparison_metrics"]["comparison_solver"] == (
        "monte_carlo"
    )
    assert manifest["comparison_metrics"]["presentation"][
        "electrical_boundaries_identified"
    ] is True
    radial = manifest["comparison_metrics"]["radial_midplane"]
    assert radial["center_height_cm"] == pytest.approx(1.27)
    assert radial["coordinate_match"] is True
    assert radial["sample_counts"] == {
        "two_term": 5,
        "monte_carlo": 5,
    }
    assert radial["presentation"]["y_scale"] == "linear"
    assert all(
        field["display_scale"] == "linear"
        for field in radial["fields"].values()
    )
    assert {
        path.name for path in summary.figures
    }.issuperset(
        {
            "comsol_style_radial_midplane_electron_density_comparison.png",
            "comsol_style_radial_midplane_electron_temperature_comparison.png",
            "comsol_style_radial_midplane_electric_potential_comparison.png",
        }
    )
    assert not any("ionization_source" in path.name for path in summary.figures)

    status_path = mc_results / "gec_ccp_run_status.json"
    status = json.loads(status_path.read_text(encoding="utf-8"))
    status["physical_target_accepted"] = False
    status_path.write_text(json.dumps(status), encoding="utf-8")
    with pytest.raises(GecCcpPlotError, match="accepted physical-target"):
        plot_gec_ccp_solver_comparison(
            two_term_bundle,
            mc_bundle,
            two_term_results_dir=two_term_results,
            monte_carlo_results_dir=mc_results,
            output_dir=tmp_path / "diagnostic",
        )

    status.update(
        {
            "status": "rejected",
            "quality_status": "failed",
            "quality_accepted": False,
            "physical_target_accepted": False,
        }
    )
    status_path.write_text(json.dumps(status), encoding="utf-8")
    with pytest.raises(GecCcpPlotError, match="failed physics-quality"):
        plot_gec_ccp_solver_comparison(
            two_term_bundle,
            mc_bundle,
            two_term_results_dir=two_term_results,
            monte_carlo_results_dir=mc_results,
            output_dir=tmp_path / "rejected",
        )


@pytest.mark.parametrize(
    "reaction_model",
    ["function_eedf", "external_rates"],
)
def test_solver_comparison_rejects_noncanonical_reaction_model(
    tmp_path: Path,
    reaction_model: str,
) -> None:
    pytest.importorskip("matplotlib")
    two_term_bundle, two_term_results = _write_comparison_case(
        tmp_path, "two_term", external_scale=1.1
    )
    mc_bundle, mc_results = _write_comparison_case(
        tmp_path,
        "monte_carlo",
        external_scale=0.9,
        reaction_model=reaction_model,
    )

    with pytest.raises(GecCcpPlotError, match="current source-specific closure"):
        plot_gec_ccp_solver_comparison(
            two_term_bundle,
            mc_bundle,
            two_term_results_dir=two_term_results,
            monte_carlo_results_dir=mc_results,
            output_dir=tmp_path / "comparison",
        )


def _write_comparison_case(
    root: Path,
    solver: str,
    *,
    external_scale: float,
    reaction_model: str | None = None,
) -> tuple[Path, Path]:
    bundle = root / f"{solver}_bundle"
    results = root / f"{solver}_results"
    bundle.mkdir()
    (results / "swarm_tables").mkdir(parents=True)
    hashes = {
        "workflow_config_sha256": "c" * 64,
        "base_config_sha256": "a" * 64,
        "cross_sections_sha256": "b" * 64,
    }
    quality_thresholds = quality_thresholds_payload(QualityThresholds())
    bundle_manifest = {
        "source": solver,
        "hashes": hashes,
        "quality_thresholds": quality_thresholds,
    }
    manifest_path = bundle / "manifest.json"
    manifest_path.write_text(json.dumps(bundle_manifest), encoding="utf-8")
    plan_path = results / "gec_ccp_plan.json"
    domain_files = [results / "swarm_tables" / "domain_period_average.csv"]
    radial_files = [results / "swarm_tables" / "radial_period_average.csv"]
    result_files = [*domain_files, *radial_files]
    _write_domain_export(domain_files[0], external_scale)
    _write_radial_export(radial_files[0], external_scale)
    if reaction_model is None:
        reaction_model = (
            "function_eedf"
            if solver == "two_term"
            else "function_eedf_preintegrated_inelastic"
        )
    plan_path.write_text(
        json.dumps(
            {
                "status": "ready",
                "model": {"input_mph": {"sha256": "d" * 64}},
                "generated_java": {
                    "apply": {"path": "apply.java", "sha256": "e" * 64}
                },
                "bundle": {
                    "path": str(bundle.resolve()),
                    "source": solver,
                    "hashes": hashes,
                    "quality_thresholds": quality_thresholds,
                    "manifest_sha256": sha256(manifest_path.read_bytes()).hexdigest(),
                },
                "closure": {
                    "electron_transport": "swarm_mobility_einstein",
                    "reaction_model": reaction_model,
                },
                "result_role": "physical_target",
                "expected_results": [str(path) for path in result_files],
            }
        ),
        encoding="utf-8",
    )
    output_mph = results / "external.mph"
    output_mph.write_bytes(b"saved mph")
    (results / "gec_ccp_run_status.json").write_text(
        json.dumps(
            {
                "status": "completed",
                "solve_status": "completed",
                "quality_status": "passed",
                "quality_accepted": True,
                "result_role": "physical_target",
                "physical_target_accepted": True,
                "plan_sha256": sha256(plan_path.read_bytes()).hexdigest(),
                "input_mph_sha256": "d" * 64,
                "generated_java_sha256": {"apply": "e" * 64},
                "bundle_artifacts_verified": True,
                "comsol_runtime": {
                    "version": "6.4",
                    "build": "429",
                    "executable": "comsolbatch.exe",
                },
                "output_mph": {
                    "path": str(output_mph),
                    "sha256": sha256(output_mph.read_bytes()).hexdigest(),
                },
                "results": _result_metadata(results, result_files),
            }
        ),
        encoding="utf-8",
    )
    return bundle, results


def _result_metadata(root: Path, paths: list[Path]) -> list[dict[str, object]]:
    metadata = []
    for path in paths:
        data = path.read_bytes()
        metadata.append(
            {
                "path": path.resolve().relative_to(root.resolve()).as_posix(),
                "size_bytes": len(data),
                "sha256": sha256(data).hexdigest(),
            }
        )
    return sorted(metadata, key=lambda item: str(item["path"]))


def _write_domain_export(path: Path, scale: float) -> None:
    coordinates = (
        (0.0000, 0.0000),
        (0.0000, 0.0254),
        (0.0538, 0.0000),
        (0.0538, 0.0254),
        (0.0538, -0.0381),
        (0.0538, 0.0635),
        (0.1016, -0.0381),
        (0.1016, 0.0635),
    )
    rows = ["% Dimension,2", "% r,z,ne,Te,V,Re"]
    for index, (radius, axial) in enumerate(coordinates, start=1):
        rows.append(
            f"{radius},{axial},{scale * index:.6g},"
            f"{scale * (2.0 + 0.1 * index):.6g},"
            f"{scale * (-10.0 + index):.6g},"
            f"{scale * index * 1.0e5:.6g}"
        )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def _write_radial_export(
    path: Path,
    scale: float,
    *,
    height_cm: float = 1.27,
) -> None:
    radii_m = (0.1016, 0.0254, 0.0, 0.0762, 0.0508)
    rows = ["% Dimension,2", "% R,Z,ne,Te,V,Re"]
    for index, radius_m in enumerate(radii_m, start=1):
        rows.append(
            f"{radius_m},{height_cm / 100.0},{scale * index:.6g},"
            f"{scale * (2.0 + 0.1 * index):.6g},"
            f"{scale * (-2.0 + index):.6g},"
            f"{scale * index * 1.0e5:.6g}"
        )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def test_radial_midplane_data_sorts_radius_and_rejects_wrong_height(
    tmp_path: Path,
) -> None:
    radial = tmp_path / "radial.csv"
    _write_radial_export(radial, 1.0)

    radius_cm, fields, height_cm = _comsol_radial_midplane_data(radial, 1.27)

    assert np.all(np.diff(radius_cm) > 0.0)
    assert fields.shape == (5, 4)
    assert height_cm == pytest.approx(1.27)

    wrong_height = tmp_path / "wrong_height.csv"
    _write_radial_export(wrong_height, 1.0, height_cm=1.30)
    with pytest.raises(GecCcpPlotError, match="not at the GEC gap midpoint"):
        _comsol_radial_midplane_data(wrong_height, 1.27)
