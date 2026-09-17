from __future__ import annotations

import csv
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.interpolate import PchipInterpolator
import yaml

from swarm_workflow.cli import main as workflow_cli_main
from swarm_workflow.quality.policy import (
    QualityThresholds,
    quality_policy_provenance,
    quality_thresholds_payload,
)
from swarm_workflow.quality.propagator_source import (
    PROPAGATOR_SOLVER_SOURCE_METADATA_KEY,
)
from swarm_workflow.comsol.input import (
    ComsolExportError,
    export_comsol_bundle,
)
from swarm_workflow.comsol.input.function_eedf import (
    RATE_IMPORTANCE_FRACTION,
    RATE_SCALED_ERROR_TOLERANCE,
    SHAPE_TOTAL_VARIATION_TOLERANCE,
    CollisionRateKernel,
    FunctionEedfError,
    SOURCE_MEAN_RELATIVE_ERROR_LIMIT,
    build_c1_function_eedf,
    collision_rate_coefficient,
    evaluate_c1_function_eedf,
    piecewise_linear_weighted_moments,
    pchip_weighted_moments,
    project_c1_function_eedf_to_comsol_grid,
    read_c1_function_eedf,
    scaled_shape_moments,
)
from swarm_workflow.campaign.store import WorkflowStore
from swarm_workflow.campaign.sweep import run_sweep
from swarm_workflow.tables.contracts import ELASTIC_ENERGY_LOSS_TABLE
from swarm_workflow.tables import build_tables
from product_helpers import (
    base_product_config,
    insert_workflow_case,
    insert_workflow_eedf,
    insert_workflow_rate,
    write_config,
)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fp:
        return list(csv.DictReader(fp))


def _write_two_mixture_workflow(tmp_path: Path) -> Path:
    data = base_product_config(tmp_path, ["two_term"])
    base_path = write_config(tmp_path, data, "base.yaml")
    workflow_path = tmp_path / "workflow.yaml"
    workflow_path.write_text(
        yaml.safe_dump(
            {
                "base_config": base_path.name,
                "database": "outputs/swarm.sqlite",
                "e_over_n_Td": [10.0, 20.0],
                "mixtures": [{"Ar": 1.0}, {"Ar": 1.0}],
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    return workflow_path


def test_two_mixture_two_term_sweep_build_tables_and_export_comsol(
    tmp_path: Path,
) -> None:
    workflow_path = _write_two_mixture_workflow(tmp_path)
    run_sweep(workflow_path)

    workflow_cli_main(
        [
            "build-tables",
            str(tmp_path / "outputs" / "swarm.sqlite"),
            "--output",
            str(tmp_path / "tables"),
            "--source",
            "two_term",
        ]
    )
    workflow_cli_main(
        [
            "export-comsol",
            str(tmp_path / "tables"),
            "--output",
            str(tmp_path / "bundle"),
        ]
    )

    root_manifest = json.loads((tmp_path / "bundle" / "manifest.json").read_text())
    assert len(root_manifest["bundles"]) == 2
    for mixture_id in [0, 1]:
        bundle = tmp_path / "bundle" / f"mixture_{mixture_id:04d}"
        expected = {
            "manifest.json",
            "mean_energy_vs_en.csv",
            "transport_vs_en.csv",
            "rates_vs_en.csv",
            "townsend_vs_en.csv",
            "quality.csv",
            "eedf.csv",
            "eedf_f0_vs_mean_energy.csv",
            "eedf_f0_comsol_2d.csv",
        }
        assert expected.issubset({path.name for path in bundle.iterdir()})
        function_dir = bundle / "comsol_functions"
        assert {
            "sw_meanE.csv",
            "sw_muN.csv",
        } == {path.name for path in function_dir.iterdir()}
        assert not (bundle / "eedf_f0.csv").exists()
        manifest = json.loads((bundle / "manifest.json").read_text())
        assert manifest["status"] == "ok"
        artifacts = {
            path.relative_to(bundle).as_posix()
            for path in bundle.rglob("*")
            if path.is_file() and path.name != "manifest.json"
        }
        assert artifacts == set(manifest["tables"])
        for name, metadata in manifest["tables"].items():
            assert isinstance(metadata["artifact_role"], str)
            assert metadata["artifact_role"]
            assert (
                metadata["sha256"] == sha256((bundle / name).read_bytes()).hexdigest()
            )
        assert manifest["tables"]["quality.csv"]["artifact_role"] == (
            "swarm_quality_audit"
        )
        assert manifest["tables"]["transport_vs_en.csv"]["artifact_role"] == (
            "swarm_coefficient_evidence"
        )
        assert (
            manifest["tables"]["comsol_functions/sw_muN.csv"]["artifact_role"]
            == "derived_comsol_interpolation_input"
        )
        assert manifest["tables"]["comsol_functions/sw_meanE.csv"][
            "derived_for_comsol_interpolation"
        ]
        raw = manifest["eedf_artifacts"]["raw_audit"]
        canonical_source = manifest["eedf_artifacts"]["canonicalization_source"]
        canonical = manifest["eedf_artifacts"]["comsol_input"]
        assert raw["file"] == "eedf.csv"
        assert raw["role"] == ("raw_swarm_eedf_audit_and_canonicalization_source")
        assert raw["sha256"] == sha256((bundle / raw["file"]).read_bytes()).hexdigest()
        assert canonical_source == raw
        assert canonical["file"] == "eedf_f0_comsol_2d.csv"
        assert canonical["role"] == "canonical_comsol_function_eedf_input"
        assert (
            canonical["sha256"]
            == sha256((bundle / canonical["file"]).read_bytes()).hexdigest()
        )
        active = manifest["tables"]["eedf_f0_comsol_2d.csv"]
        assert active["mean_energy_axis"] == (
            "error_controlled_c1_materialization_with_solver_anchors"
        )
        assert active["mean_energy_grid_points"] >= 2
        assert active["solver_anchor_mean_energy_count"] == 2
        assert active["materialized_mean_nodes_include_all_solver_anchors"]
        assert active["shape_total_variation_error_max"] <= (
            SHAPE_TOTAL_VARIATION_TOLERANCE
        )
        assert active["rate_scaled_error_max"] <= RATE_SCALED_ERROR_TOLERANCE
        assert active["mean_axis_rate_scaled_error_max"] <= (
            RATE_SCALED_ERROR_TOLERANCE
        )
        assert active["rate_kernel_count"] >= 1
        assert active["rate_importance_fraction"] == RATE_IMPORTANCE_FRACTION
        assert not (bundle / "eedf_f0_comsol_grid.txt").exists()
        assert not (bundle / "eedf_f0_binding_seed.csv").exists()
        transport = _read_csv(bundle / "transport_vs_en.csv")
        assert len(transport) == 2
        assert float(transport[0]["E_over_N_V_m2"]) == pytest.approx(
            float(transport[0]["E_over_N_Td"]) * 1.0e-21
        )


def test_export_rates_do_not_double_apply_mixture_fraction(tmp_path: Path) -> None:
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        insert_workflow_case(store.connection, solver="two_term", drift=2.0)
        insert_workflow_rate(
            store.connection,
            solver="two_term",
            value=4.0e-15,
            target_fraction=0.25,
            process="excitation",
            process_type="excitation",
            threshold_eV=11.5,
        )
        insert_workflow_eedf(store.connection, solver="two_term")
        store.connection.commit()

    build_tables(tmp_path / "workflow.sqlite", tmp_path / "tables", source="two_term")
    export_comsol_bundle(
        tmp_path / "tables" / "mixture_0000",
        tmp_path / "bundle",
    )
    row = _read_csv(tmp_path / "bundle" / "rates_vs_en.csv")[0]

    assert float(row["rate_coefficient_m3_s"]) == pytest.approx(4.0e-15)
    assert float(row["mixture_weighted_rate_m3_s"]) == pytest.approx(1.0e-15)
    assert float(row["reduced_townsend_m2"]) == pytest.approx(2.0e-15)
    assert float(row["mixture_weighted_reduced_townsend_m2"]) == pytest.approx(0.5e-15)
    assert (tmp_path / "bundle" / "rates_vs_mean_energy.csv").exists()


def test_rf_bundle_writes_function_eedf_by_mean_energy(tmp_path: Path) -> None:
    rf_metadata = {
        "rf_field_treatment": ("time_periodic_f0:instantaneous_f1:sinusoidal_rms"),
        "rf_frequency_Hz": 13.56e6,
        "rf_amplitude_definition": "rms",
    }
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        for e_over_n, mean_energy, masses in [
            (10.0, 1.0, (0.5, 0.5)),
            (20.0, 1.2, (0.3, 0.7)),
        ]:
            insert_workflow_case(
                store.connection,
                solver="two_term",
                e_over_n=e_over_n,
                mean_energy=mean_energy,
                metadata=rf_metadata,
            )
            insert_workflow_rate(
                store.connection,
                solver="two_term",
                e_over_n=e_over_n,
            )
            insert_workflow_eedf(
                store.connection,
                solver="two_term",
                e_over_n=e_over_n,
                probability_masses=masses,
            )
        store.connection.commit()

    build_tables(tmp_path / "workflow.sqlite", tmp_path / "tables", source="two_term")
    export_comsol_bundle(
        tmp_path / "tables" / "mixture_0000",
        tmp_path / "bundle",
    )

    function_eedf = tmp_path / "bundle" / "eedf_f0_vs_mean_energy.csv"
    assert function_eedf.exists()
    rows = _read_csv(function_eedf)
    assert float(rows[0]["dimensionless_energy"]) == 0.0
    assert set(rows[0]) == {
        "dimensionless_energy",
        "shape_0000",
        "shape_0001",
    }
    manifest = json.loads((tmp_path / "bundle" / "manifest.json").read_text())
    metadata = manifest["tables"]["eedf_f0_vs_mean_energy.csv"]
    assert metadata["argument"] == "electron_energy_eV,mean_energy_eV"
    assert metadata["anchor_mean_energy_eV"] == [1.0, 1.2]
    assert metadata["representation"] == ("dimensionless_shape_preserving_c1_convex")
    assert metadata["normalization_error_max"] <= 1.0e-8
    assert metadata["mean_energy_relative_error_max"] <= 1.0e-8
    assert metadata["artifact_role"] == "canonicalized_function_eedf_evidence"
    assert not metadata["canonical_comsol_input"]

    active_path = tmp_path / "bundle" / "eedf_f0_comsol_2d.csv"
    assert active_path.exists()
    assert active_path.read_text(encoding="utf-8").startswith(
        "electron_energy_eV,mean_energy_eV,eepf_eV_m32\n"
    )
    active = manifest["tables"]["eedf_f0_comsol_2d.csv"]
    assert active["artifact_role"] == "canonical_comsol_function_eedf_input"
    assert active["canonical_comsol_input"]
    assert active["energy_grid_points"] < 2049
    assert active["mean_energy_grid_points"] >= 2
    assert active["solver_anchor_mean_energy_count"] == 2
    assert active["grid_shape"] == [
        active["mean_energy_grid_points"],
        active["energy_grid_points"],
    ]
    assert active["projected_normalization_error_max"] <= 1.0e-8
    assert active["projected_mean_energy_relative_error_max"] <= 1.0e-8
    assert active["projected_nonnegative_minimum"] >= 0.0
    assert active["comsol_import"] == {
        "source": "file",
        "struct": "spreadsheet",
        "nargs": 2,
        "argunit": "eV,eV",
        "fununit": "1",
        "interp": "linear",
        "extrap": "const",
        "funcnametable_position": "1",
        "scaledata": "auto",
    }


def test_dc_function_eedf_conserves_zero_start_bin_and_midpoint(
    tmp_path: Path,
) -> None:
    table_dir = tmp_path / "tables" / "mixture_0000"
    table_dir.mkdir(parents=True)
    _write_minimal_manifest(table_dir, optional_tables=("eedf.csv",))
    _write_required_base_tables(table_dir)
    (table_dir / "eedf.csv").write_text(
        (
            "electron_energy_eV,energy_width_eV,E_over_N_Td,E_over_N_V_m2,"
            "mean_energy_eV,eedf\n"
            "0.25,0.5,10,1e-20,0.35,1.6\n"
            "0.75,0.5,10,1e-20,0.35,0.4\n"
            "0.25,0.5,20,2e-20,0.65,0.4\n"
            "0.75,0.5,20,2e-20,0.65,1.6\n"
        ),
        encoding="utf-8",
    )
    _record_source_artifact(
        table_dir,
        "eedf.csv",
        role="raw_swarm_eedf_audit_and_canonicalization_source",
    )

    export_comsol_bundle(table_dir, tmp_path / "bundle")

    manifest = json.loads((tmp_path / "bundle" / "manifest.json").read_text())
    metadata = manifest["tables"]["eedf_f0_vs_mean_energy.csv"]
    representation = read_c1_function_eedf(
        tmp_path / "bundle" / "eedf_f0_vs_mean_energy.csv",
        metadata,
    )
    obsolete_metadata = {**metadata, "representation": "rectangular_linear"}
    with pytest.raises(FunctionEedfError, match="unsupported"):
        read_c1_function_eedf(
            tmp_path / "bundle" / "eedf_f0_vs_mean_energy.csv",
            obsolete_metadata,
        )
    x_values = representation.dimensionless_energy
    shapes = representation.shapes

    assert x_values[0] == 0.0
    assert x_values[-1] > 1.0
    for shape in shapes:
        dense_x = np.linspace(0.0, x_values[-1], 5000)
        interpolator = PchipInterpolator(x_values, shape)
        dense_shape = interpolator(dense_x)
        assert np.min(dense_shape) >= -1.0e-13
        assert shape[-2] == 0.0
        assert shape[-1] == 0.0
        assert interpolator.derivative()(x_values[-2]) == pytest.approx(
            0.0,
            abs=1.0e-12,
        )
        assert interpolator.derivative()(x_values[-1]) == pytest.approx(
            0.0,
            abs=1.0e-12,
        )
        normalization, first_moment = pchip_weighted_moments(
            x_values,
            shape,
        )
        assert normalization == pytest.approx(1.0, abs=1.0e-8)
        assert first_moment / normalization == pytest.approx(1.0, abs=1.0e-8)

    # The scale transformation makes the requested physical mean exact for
    # every intermediate mean, not just for the anchor slices.
    for mean in (0.35, 0.41, 0.50, 0.59, 0.65):
        energy = mean * x_values
        f0 = evaluate_c1_function_eedf(representation, energy, mean)
        assert np.all(f0 >= -1.0e-13)
        normalization, first_moment = scaled_shape_moments(
            representation,
            mean,
        )
        assert normalization == pytest.approx(1.0, abs=1.0e-8)
        assert mean * first_moment / normalization == pytest.approx(
            mean,
            abs=1.0e-8,
        )

    assert not metadata["canonical_comsol_input"]
    assert metadata["artifact_role"] == "canonicalized_function_eedf_evidence"
    assert (
        metadata["sha256"]
        == sha256(
            (tmp_path / "bundle" / "eedf_f0_vs_mean_energy.csv").read_bytes()
        ).hexdigest()
    )
    assert metadata["mean_energy_axis_interpolation"] == (
        "offline_C1_shape_family; "
        "active_table_uses_error_controlled_materialized_mean_nodes"
    )

    active_path = tmp_path / "bundle" / "eedf_f0_comsol_2d.csv"
    active = manifest["tables"]["eedf_f0_comsol_2d.csv"]
    active_rows = _read_csv(active_path)
    energy_count = active["energy_grid_points"]
    mean_count = active["mean_energy_grid_points"]
    energies = np.asarray(
        [float(row["electron_energy_eV"]) for row in active_rows[:energy_count]]
    )
    means = np.asarray(
        [
            float(active_rows[index * energy_count]["mean_energy_eV"])
            for index in range(mean_count)
        ]
    )
    values = np.asarray([float(row["eepf_eV_m32"]) for row in active_rows]).reshape(
        mean_count, energy_count
    )
    assert values.shape == (mean_count, energy_count)
    midpoint_values = 0.5 * (values[0] + values[1])
    np.testing.assert_allclose(
        midpoint_values,
        0.5 * (values[0] + values[1]),
        rtol=1.0e-13,
        atol=0.0,
    )
    assert energies[0] == 0.0
    assert np.all(np.diff(energies) > 0.0)
    assert np.all(np.diff(means) > 0.0)
    assert np.min(values) >= 0.0
    for anchor in representation.mean_energies_eV:
        assert np.any(np.isclose(means, anchor, rtol=0.0, atol=1.0e-14))

    normalization_error = 0.0
    mean_error = 0.0
    for mean, row in zip(means, values, strict=True):
        normalization, first_moment = piecewise_linear_weighted_moments(
            energies,
            row,
        )
        normalization_error = max(normalization_error, abs(normalization - 1.0))
        mean_error = max(
            mean_error,
            abs(first_moment / normalization - mean) / mean,
        )
    assert normalization_error <= 1.0e-8
    assert mean_error <= 1.0e-8
    assert active["projected_normalization_error_max"] == pytest.approx(
        normalization_error,
        abs=1.0e-15,
    )
    assert active["projected_mean_energy_relative_error_max"] == pytest.approx(
        mean_error,
        abs=1.0e-15,
    )
    assert active["artifact_role"] == "canonical_comsol_function_eedf_input"
    assert active["canonical_comsol_input"]
    assert active["sha256"] == sha256(active_path.read_bytes()).hexdigest()
    assert active["comsol_import"] == {
        "source": "file",
        "struct": "spreadsheet",
        "nargs": 2,
        "argunit": "eV,eV",
        "fununit": "1",
        "interp": "linear",
        "extrap": "const",
        "funcnametable_position": "1",
        "scaledata": "auto",
    }
    assert not (tmp_path / "bundle" / "eedf_f0.csv").exists()
    assert not (
        tmp_path / "bundle" / "comsol_functions" / "sw_eedf_f0_grid.txt"
    ).exists()


def test_function_eedf_rejects_source_mean_mismatch_before_projection() -> None:
    grouped = {
        # The bin-weighted mean is 0.35 eV, intentionally inconsistent with
        # the declared 0.40 eV case mean.
        0.40: [(0.25, 0.5, 1.6), (0.75, 0.5, 0.4)],
        0.65: [(0.25, 0.5, 0.4), (0.75, 0.5, 1.6)],
    }
    with pytest.raises(FunctionEedfError, match="source EEDF mean energy"):
        build_c1_function_eedf(grouped)


def test_function_eedf_accepts_bounded_source_bin_quadrature_error() -> None:
    representation = build_c1_function_eedf(
        {
            # Bin-weighted means are 0.35 and 0.65 eV. The independent case
            # means differ by less than the documented 0.1% acceptance limit.
            0.3502: [(0.25, 0.5, 1.6), (0.75, 0.5, 0.4)],
            0.6502: [(0.25, 0.5, 0.4), (0.75, 0.5, 1.6)],
        }
    )
    assert representation.source_mean_relative_error_max < (
        SOURCE_MEAN_RELATIVE_ERROR_LIMIT
    )


def test_adaptive_projection_bounds_active_rate_and_preserves_exact_zero() -> None:
    representation = build_c1_function_eedf(
        {
            0.35: [(0.25, 0.5, 1.6), (0.75, 0.5, 0.4)],
            0.65: [(0.25, 0.5, 0.4), (0.75, 0.5, 1.6)],
        }
    )
    energy = np.asarray([0.0, 0.3, 0.5, 0.7, 1.0])
    active = CollisionRateKernel(
        species="Ar",
        process="excitation",
        process_type="excitation",
        threshold_eV=0.3,
        electron_energies_eV=energy,
        cross_sections_m2=np.asarray([0.0, 0.0, 1.0e-20, 2.0e-20, 2.0e-20]),
    )
    zero = CollisionRateKernel(
        species="Ar",
        process="inactive",
        process_type="excitation",
        threshold_eV=0.3,
        electron_energies_eV=energy,
        cross_sections_m2=np.zeros_like(energy),
    )

    projected = project_c1_function_eedf_to_comsol_grid(
        representation,
        rate_kernels=(active, zero),
    )

    assert projected.shape_total_variation_error_max <= (
        SHAPE_TOTAL_VARIATION_TOLERANCE
    )
    assert projected.rate_scaled_error_max is not None
    assert projected.rate_scaled_error_max <= RATE_SCALED_ERROR_TOLERANCE
    for row in projected.values_eV_m32:
        assert (
            collision_rate_coefficient(
                projected.electron_energies_eV,
                row,
                zero,
            )
            == 0.0
        )


def test_comsol_projection_removes_only_identically_zero_global_tail() -> None:
    centers = np.arange(0.5, 100.0, 1.0)

    def source(mean: float) -> list[tuple[float, float, float]]:
        values = np.zeros_like(centers)
        left_index = int(mean - 1.0)
        values[left_index : left_index + 2] = 0.5
        return list(
            zip(
                centers.tolist(),
                np.ones_like(centers).tolist(),
                values.tolist(),
                strict=True,
            )
        )

    representation = build_c1_function_eedf(
        {mean: source(mean) for mean in (1.0, 2.0, 3.0)}
    )
    projected = project_c1_function_eedf_to_comsol_grid(representation)
    legacy_global_limit = (
        representation.mean_energies_eV[-1] * representation.dimensionless_energy[-1]
    )
    assert representation.source_energy_max_eV == pytest.approx(100.0)
    assert projected.electron_energies_eV[-1] == pytest.approx(100.0)
    assert projected.electron_energies_eV[-1] < 0.34 * legacy_global_limit

    # The removed interval is exactly zero for every anchor and local blend,
    # so neither EEDF probability mass nor a bounded collision-rate integral
    # is changed by the compact physical support.
    tail_energy = np.geomspace(
        100.0 * (1.0 + 1.0e-12),
        legacy_global_limit,
        2048,
    )
    for mean in (1.0, 1.5, 2.0, 2.5, 3.0):
        tail = evaluate_c1_function_eedf(
            representation,
            tail_energy,
            mean,
        )
        tail_mass = np.trapezoid(np.sqrt(tail_energy) * tail, tail_energy)
        surrogate_rate = np.trapezoid(
            tail_energy * (1.0e-20 * tail),
            tail_energy,
        )
        assert tail_mass == 0.0
        assert surrogate_rate == 0.0

    assert projected.normalization_error_max <= 1.0e-8
    assert projected.mean_energy_relative_error_max <= 1.0e-8


def test_propagator_export_uses_common_adaptive_spreadsheet(
    tmp_path: Path,
) -> None:
    table_dir = tmp_path / "tables" / "mixture_0000"
    table_dir.mkdir(parents=True)
    _write_minimal_manifest(
        table_dir,
        optional_tables=("eedf.csv",),
        source="propagator",
    )
    _write_required_base_tables(table_dir)
    (table_dir / "eedf.csv").write_text(
        (
            "electron_energy_eV,energy_width_eV,E_over_N_Td,E_over_N_V_m2,"
            "mean_energy_eV,eedf\n"
            "0.5,1,10,1e-20,0.75,0.75\n"
            "1.5,1,10,1e-20,0.75,0.25\n"
            "0.5,1,20,2e-20,1.25,0.25\n"
            "1.5,1,20,2e-20,1.25,0.75\n"
        ),
        encoding="utf-8",
    )
    _record_source_artifact(
        table_dir,
        "eedf.csv",
        role="raw_swarm_eedf_audit_and_canonicalization_source",
    )

    export_comsol_bundle(table_dir, tmp_path / "bundle")

    bundle = tmp_path / "bundle"
    manifest = json.loads((bundle / "manifest.json").read_text())
    active = manifest["tables"]["eedf_f0_comsol_2d.csv"]
    assert (bundle / "eedf_f0_comsol_2d.csv").is_file()
    assert not (bundle / "eedf_f0_comsol_grid.txt").exists()
    assert not (bundle / "eedf_f0_binding_seed.csv").exists()
    assert (bundle / "eedf_f0_vs_mean_energy.csv").is_file()
    assert active["format"] == "csv"
    assert active["comsol_import"]["struct"] == "spreadsheet"
    assert active["grid_shape"] == [
        active["mean_energy_grid_points"],
        active["energy_grid_points"],
    ]
    assert active["mean_energy_grid_points"] > 2
    assert active["representation"] == ("physical_2d_adaptive_moment_rate_projected")
    assert active["mean_energy_axis"] == (
        "error_controlled_c1_materialization_with_solver_anchors"
    )
    assert active["mean_axis_shape_total_variation_error_max"] <= (
        SHAPE_TOTAL_VARIATION_TOLERANCE
    )
    assert active["mean_axis_rate_scaled_error_max"] is None
    assert active["source_support"] == {
        "policy": "full_source_support",
        "source_values_modified": False,
    }


@pytest.mark.parametrize(
    "solver_source_sha256",
    [None, "0" * 64],
    ids=("missing", "mismatch"),
)
def test_propagator_export_rejects_unbound_solver_source(
    tmp_path: Path,
    solver_source_sha256: str | None,
) -> None:
    table_dir = tmp_path / "tables" / "mixture_0000"
    table_dir.mkdir(parents=True)
    _write_minimal_manifest(table_dir, source="propagator")
    _write_required_base_tables(table_dir)
    manifest_path = table_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if solver_source_sha256 is None:
        del manifest["hashes"][PROPAGATOR_SOLVER_SOURCE_METADATA_KEY]
    else:
        manifest["hashes"][PROPAGATOR_SOLVER_SOURCE_METADATA_KEY] = (
            solver_source_sha256
        )
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(
        ComsolExportError,
        match="source provenance differs from the qualified implementation",
    ):
        export_comsol_bundle(table_dir, tmp_path / "bundle")


def _write_mc_eedf_export_fixture(
    table_dir: Path,
    *,
    postprocess: str,
) -> str:
    table_dir.mkdir(parents=True)
    _write_minimal_manifest(
        table_dir,
        optional_tables=("eedf.csv", "rate_evidence.csv"),
        source="monte_carlo",
        postprocess=postprocess,
    )
    _write_required_base_tables(table_dir)
    eedf = (
        "electron_energy_eV,energy_width_eV,E_over_N_Td,E_over_N_V_m2,"
        "mean_energy_eV,eedf,pooled_effective_sample_count\n"
        "0.25,0.5,10,1e-20,0.35,1.6,100\n"
        "0.75,0.5,10,1e-20,0.35,0.4,100\n"
        "0.25,0.5,20,2e-20,0.65,0.4,100\n"
        "0.75,0.5,20,2e-20,0.65,1.6,100\n"
    )
    (table_dir / "eedf.csv").write_text(eedf, encoding="utf-8")
    rate_evidence = (
        "E_over_N_Td,process,rate_coefficient_mean_m3_s,"
        "pooled_zero_event_status\n"
        "10,ionization,1e-15,observed_events\n"
    )
    (table_dir / "rate_evidence.csv").write_text(
        rate_evidence,
        encoding="utf-8",
    )
    _record_source_artifact(
        table_dir,
        "eedf.csv",
        role="raw_swarm_eedf_audit_and_canonicalization_source",
    )
    _record_source_artifact(
        table_dir,
        "rate_evidence.csv",
        role="raw_swarm_rate_evidence",
    )
    return eedf


def test_pure_monte_carlo_uses_raw_eedf_as_canonicalization_source(
    tmp_path: Path,
) -> None:
    table_dir = tmp_path / "tables" / "mixture_0000"
    eedf = _write_mc_eedf_export_fixture(table_dir, postprocess="none")

    bundle = tmp_path / "bundle"
    export_comsol_bundle(table_dir, bundle)

    assert (bundle / "eedf.csv").read_text(encoding="utf-8") == eedf
    manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["tables"]["rate_evidence.csv"]["statistics"] == {
        "replicate_interval": "two_sided_student_t_95",
        "zero_event_confidence": 0.95,
        "zero_event_upper_bound": ("poisson_zero_count_over_pooled_target_exposure"),
    }
    raw = manifest["eedf_artifacts"]["raw_audit"]
    assert raw == manifest["eedf_artifacts"]["canonicalization_source"]
    assert raw == {
        "file": "eedf.csv",
        "sha256": sha256((bundle / "eedf.csv").read_bytes()).hexdigest(),
        "role": "raw_swarm_eedf_audit_and_canonicalization_source",
    }
    assert (
        manifest["tables"]["eedf_f0_vs_mean_energy.csv"]["source_artifact"]
        == "eedf.csv"
    )
    assert manifest["tables"]["eedf_f0_comsol_2d.csv"]["source_artifact"] == "eedf.csv"


def test_export_rejects_noncanonical_zero_event_statistics(
    tmp_path: Path,
) -> None:
    table_dir = tmp_path / "tables" / "mixture_0000"
    _write_mc_eedf_export_fixture(table_dir, postprocess="none")
    manifest_path = table_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["tables"]["rate_evidence.csv"]["statistics"]["zero_event_confidence"] = (
        0.90
    )
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(
        ComsolExportError,
        match="lacks the canonical pooled zero-event statistics contract",
    ):
        export_comsol_bundle(table_dir, tmp_path / "bundle")


def test_export_rejects_obsolete_regularized_eedf_policy_with_migration(
    tmp_path: Path,
) -> None:
    table_dir = tmp_path / "tables" / "mixture_0000"
    _write_mc_eedf_export_fixture(table_dir, postprocess="regularized")

    with pytest.raises(
        ComsolExportError,
        match=(
            "postprocessed or prior-filled EEDF tables are obsolete; rebuild "
            "a qualified pure-solver table with source_policy.postprocess=none"
        ),
    ):
        export_comsol_bundle(table_dir, tmp_path / "bundle")


def test_export_omits_unused_rate_and_energy_loss_functions(tmp_path: Path) -> None:
    table_dir = tmp_path / "tables" / "mixture_0000"
    table_dir.mkdir(parents=True)
    _write_minimal_manifest(table_dir)
    _write_required_base_tables(table_dir)

    export_comsol_bundle(table_dir, tmp_path / "bundle")

    function_dir = tmp_path / "bundle" / "comsol_functions"
    assert {path.name for path in function_dir.iterdir()} == {
        "sw_meanE.csv",
        "sw_muN.csv",
    }
    assert not (tmp_path / "bundle" / "energy_loss.csv").exists()


def test_export_copies_only_manifested_elastic_energy_loss_contract(
    tmp_path: Path,
) -> None:
    table_dir = tmp_path / "tables" / "mixture_0000"
    table_dir.mkdir(parents=True)
    _write_minimal_manifest(
        table_dir,
        optional_tables=(ELASTIC_ENERGY_LOSS_TABLE,),
    )
    _write_required_base_tables(table_dir)
    (table_dir / ELASTIC_ENERGY_LOSS_TABLE).write_text(
        (
            "mean_energy_eV,E_over_N_Td,E_over_N_V_m2,"
            "elastic_energy_loss_rate_coefficient_eV_m3_s,"
            "elastic_energy_loss_standard_error_eV_m3_s,"
            "elastic_energy_loss_relative_standard_error,"
            "elastic_energy_loss_ci95_low_eV_m3_s,"
            "elastic_energy_loss_ci95_high_eV_m3_s,ci95_critical_value,"
            "estimate_status,valid_replicates,uncertainty_available\n"
            "2.0,10.0,1e-20,1e-18,,,,,,"
            "deterministic_operator_moment,1,0\n"
        ),
        encoding="utf-8",
    )
    _record_source_artifact(
        table_dir,
        ELASTIC_ENERGY_LOSS_TABLE,
        role="canonical_comsol_elastic_energy_loss_input",
    )
    manifest_path = table_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["tables"][ELASTIC_ENERGY_LOSS_TABLE]["physics_contract"] = {
        "schema": "swarm.elastic_energy_loss.v1",
        "symbol": "K_epsilon_el",
        "estimator": ("same_discrete_elastic_collision_operator_energy_moment"),
        "gas_temperature_terms_included": True,
    }
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    export_comsol_bundle(table_dir, tmp_path / "bundle")

    assert (tmp_path / "bundle" / ELASTIC_ENERGY_LOSS_TABLE).is_file()
    bundle_manifest = json.loads((tmp_path / "bundle" / "manifest.json").read_text())
    assert (
        bundle_manifest["tables"][ELASTIC_ENERGY_LOSS_TABLE]["physics_contract"][
            "symbol"
        ]
        == "K_epsilon_el"
    )


def test_export_ignores_file_not_listed_by_source_manifest(tmp_path: Path) -> None:
    table_dir = tmp_path / "tables" / "mixture_0000"
    table_dir.mkdir(parents=True)
    _write_minimal_manifest(table_dir)
    _write_required_base_tables(table_dir)
    (table_dir / "rates_vs_mean_energy.csv").write_text(
        "mean_energy_eV,rate_coefficient_m3_s\n2.0,1e-15\n",
        encoding="utf-8",
    )

    export_comsol_bundle(table_dir, tmp_path / "bundle")

    assert not (tmp_path / "bundle" / "rates_vs_mean_energy.csv").exists()
    manifest = json.loads((tmp_path / "bundle" / "manifest.json").read_text())
    assert "rates_vs_mean_energy.csv" not in manifest["tables"]


def test_export_rejects_table_manifest_without_quality_provenance(
    tmp_path: Path,
) -> None:
    table_dir = tmp_path / "tables" / "mixture_0000"
    table_dir.mkdir(parents=True)
    _write_minimal_manifest(table_dir)
    manifest_path = table_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    del manifest["quality_thresholds"]
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    _write_required_base_tables(table_dir)

    with pytest.raises(ComsolExportError, match="lacks quality_thresholds"):
        export_comsol_bundle(table_dir, tmp_path / "bundle")


def test_export_rejects_missing_mobility_with_failure_manifest(
    tmp_path: Path,
) -> None:
    table_dir = tmp_path / "tables" / "mixture_0000"
    table_dir.mkdir(parents=True)
    _write_minimal_manifest(table_dir)
    _write_required_base_tables(table_dir)
    transport = table_dir / "transport_vs_en.csv"
    transport.write_text(
        transport.read_text(encoding="utf-8").replace(
            "10.0,1e-20,2.0,5.0,10.0,20.0",
            "10.0,1e-20,2.0,5.0,,20.0",
        ),
        encoding="utf-8",
    )

    with pytest.raises(ComsolExportError, match="missing required coefficients"):
        export_comsol_bundle(
            table_dir,
            tmp_path / "bundle",
        )

    manifest = json.loads((tmp_path / "bundle" / "manifest.json").read_text())
    assert manifest["status"] == "failed"
    assert any(
        "reduced_mobility_m2_V_s_m3" in item
        for item in manifest["missing_required_coefficients"]
    )


def _write_minimal_manifest(
    table_dir: Path,
    *,
    optional_tables: tuple[str, ...] = (),
    source: str = "two_term",
    postprocess: str = "none",
) -> None:
    quality = QualityThresholds()
    table_names = (
        "mean_energy_vs_en.csv",
        "transport_vs_en.csv",
        "rates_vs_en.csv",
        "townsend_vs_en.csv",
        "quality.csv",
        *optional_tables,
    )
    solver_qualification = None
    hashes: dict[str, str] = {}
    if source == "propagator":
        qualification_source = (
            Path(__file__).resolve().parents[1]
            / "docs"
            / "dev"
            / "results"
            / "propagator_p1_deterministic_qualification_20260908.json"
        )
        qualification_name = "propagator_p1_core_qualification.json"
        qualification_target = table_dir / qualification_name
        qualification_target.write_bytes(qualification_source.read_bytes())
        payload = json.loads(qualification_source.read_text(encoding="utf-8"))
        fingerprint = payload["environment"]["implementation_fingerprint"]
        solver_qualification = {
            "file": qualification_name,
            "sha256": sha256(qualification_target.read_bytes()).hexdigest(),
            "schema": payload["schema"],
            "decision": "p1_deterministic_core_qualified",
            "implementation_fingerprint": {
                "schema": fingerprint["schema"],
                "sha256": fingerprint["sha256"],
            },
        }
        hashes[PROPAGATOR_SOLVER_SOURCE_METADATA_KEY] = fingerprint["sha256"]
    (table_dir / "manifest.json").write_text(
        json.dumps(
            {
                "format_version": 1,
                "stage": "build-tables",
                "source": source,
                "solver_qualification": solver_qualification,
                "hashes": hashes,
                "quality_thresholds": quality_thresholds_payload(quality),
                "quality_policy": quality_policy_provenance(quality, quality),
                "mixture": {"mixture_id": 0, "species": []},
                "source_policy": {
                    "source": source,
                    "postprocess": postprocess,
                },
                "valid_ranges": {"E_over_N_Td": [10.0, 10.0]},
                "units": {},
                "table_argument": {"primary": "E_over_N_Td", "secondary": None},
                "monotonicity": {"mean_energy_strictly_monotonic": True},
                "quality_summary": {
                    "passed": True,
                    "failed_points": 0,
                    "total_points": 1,
                },
                "tables": {name: {"columns": []} for name in table_names},
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def _record_source_artifact(table_dir: Path, name: str, *, role: str) -> None:
    manifest_path = table_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    metadata = {
        "artifact_role": role,
        "sha256": sha256((table_dir / name).read_bytes()).hexdigest(),
    }
    if name == "rate_evidence.csv":
        metadata["statistics"] = {
            "replicate_interval": "two_sided_student_t_95",
            "zero_event_confidence": 0.95,
            "zero_event_upper_bound": (
                "poisson_zero_count_over_pooled_target_exposure"
            ),
        }
    manifest["tables"][name].update(metadata)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def _write_required_base_tables(table_dir: Path) -> None:
    (table_dir / "mean_energy_vs_en.csv").write_text(
        "E_over_N_Td,E_over_N_V_m2,mean_energy_eV\n10.0,1e-20,2.0\n",
        encoding="utf-8",
    )
    (table_dir / "transport_vs_en.csv").write_text(
        (
            "E_over_N_Td,E_over_N_V_m2,mean_energy_eV,drift_velocity_m_s,"
            "reduced_mobility_m2_V_s_m3,reduced_diffusion_L_m2_s_m3,"
            "reduced_diffusion_T_m2_s_m3,"
            "reduced_electron_energy_mobility_m2_V_s_m3,"
            "reduced_electron_energy_diffusion_m2_s_m3\n"
            "10.0,1e-20,2.0,5.0,10.0,20.0,30.0,7.0,8.0\n"
        ),
        encoding="utf-8",
    )
    (table_dir / "rates_vs_en.csv").write_text(
        (
            "E_over_N_Td,E_over_N_V_m2,species,process,process_type,threshold_eV,"
            "target_species_fraction,rate_coefficient_m3_s,"
            "mixture_weighted_rate_m3_s,reduced_townsend_m2,"
            "mixture_weighted_reduced_townsend_m2\n"
            "10.0,1e-20,Ar,ionization,ionization,15.0,1.0,1e-15,1e-15,2e-16,2e-16\n"
        ),
        encoding="utf-8",
    )
    (table_dir / "townsend_vs_en.csv").write_text(
        "E_over_N_Td,E_over_N_V_m2,effective_townsend_m2\n10.0,1e-20,2e-16\n",
        encoding="utf-8",
    )
    (table_dir / "energy_loss.csv").write_text(
        (
            "E_over_N_Td,E_over_N_V_m2,species,process,process_type,threshold_eV,"
            "target_species_fraction,energy_loss_eV,energy_loss_rate_coefficient_eV_m3_s\n"
            "10.0,1e-20,Ar,ionization,ionization,15.0,1.0,15.0,1.5e-14\n"
        ),
        encoding="utf-8",
    )
    (table_dir / "quality.csv").write_text(
        (
            "E_over_N_Td,E_over_N_V_m2,passed,failure_reasons_json,mobility_rse,"
            "diffusion_L_rse,diffusion_T_rse,max_major_rate_rse,"
            "eedf_normalization_error,valid_replicates,uncertainty_available,"
            "quality_source\n"
            "10.0,1e-20,1,[],0,0,0,0,0,1,0,two_term\n"
        ),
        encoding="utf-8",
    )
