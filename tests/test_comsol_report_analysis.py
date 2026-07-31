from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sqlite3
import sys

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = (
    ROOT
    / "reports"
    / "comsol_swarm_benchmark_2026"
    / "repro"
    / "benchmark_analysis.py"
)
SPEC = importlib.util.spec_from_file_location("benchmark_analysis_test", MODULE_PATH)
assert SPEC and SPEC.loader
analysis = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = analysis
SPEC.loader.exec_module(analysis)


def test_tail_audit_separates_grid_tail_from_solver_convergence(
    tmp_path: Path,
) -> None:
    database = tmp_path / "swarm.sqlite"
    with sqlite3.connect(database) as connection:
        connection.execute(
            """
            CREATE TABLE cases (
                solver TEXT,
                e_over_n_Td REAL,
                replicate INTEGER,
                diagnostics_json TEXT
            )
            """
        )
        connection.executemany(
            "INSERT INTO cases VALUES (?, ?, ?, ?)",
            [
                (
                    "two_term",
                    1000.0,
                    0,
                    json.dumps(
                        {
                            "two_term": {
                                "converged": True,
                                "iterations": 42,
                                "residual_L1": 1.0e-10,
                                "adaptive_cycles": 4,
                                "grid_max_eV": 400.0,
                                "tail_probability": 1.0e-11,
                                "edge_to_peak": 1.0e-12,
                            }
                        }
                    ),
                ),
                (
                    "two_term",
                    1.0,
                    0,
                    json.dumps(
                        {
                            "two_term": {
                                "converged": False,
                                "iterations": 600,
                                "residual_L1": 2.0e-7,
                                "adaptive_cycles": 1,
                                "grid_max_eV": 100.0,
                                "tail_probability": 1.0e-20,
                                "edge_to_peak": 1.0e-30,
                            }
                        }
                    ),
                ),
            ],
        )
    result = analysis.tail_convergence_audit(
        database,
        tail_probability_target=1.0e-9,
        edge_to_peak_target=1.0e-10,
        max_energy_limit_eV=20000.0,
    )
    assert result["formal_tail_grid_gate_pass"].all()
    assert result["formal_solver_convergence_gate_pass"].tolist() == [
        False,
        True,
    ]
    assert result["formal_swarm_case_gate_pass"].tolist() == [False, True]


def test_weighted_l2_uses_exact_piecewise_linear_products() -> None:
    x = np.array([0.0, 0.1, 1.0])
    reference = np.ones(3)
    external = np.array([2.0, 1.0, 1.0])
    numerator = 0.1 / 3.0
    denominator = 1.0
    expected = np.sqrt(numerator / denominator)
    frame_external = pd.DataFrame(
        {
            "x": x,
            **{
                field: external if field == "electron_density" else reference
                for field, _, _ in analysis.QUANTITIES
            },
        }
    )
    frame_reference = pd.DataFrame(
        {"x": x, **{field: reference for field, _, _ in analysis.QUANTITIES}}
    )
    metrics, _ = analysis.comparison_metrics(
        frame_external, frame_reference, "Td"
    )
    observed = metrics.loc[
        metrics["quantity"] == "electron_density",
        "relative_L2_spatial_weighted",
    ].iloc[0]
    assert observed == pytest.approx(expected)
    assert observed != pytest.approx(np.linalg.norm(external - reference) / np.linalg.norm(reference))
    trapezoidal_sensitivity = metrics.loc[
        metrics["quantity"] == "electron_density",
        "relative_L2_trapezoidal_of_squared_samples_sensitivity",
    ].iloc[0]
    assert trapezoidal_sensitivity == pytest.approx(
        np.sqrt(np.trapezoid((external - reference) ** 2, x))
    )
    assert observed != pytest.approx(trapezoidal_sensitivity)


def test_piecewise_linear_product_integral_is_exact_for_linear_product() -> None:
    x = np.array([0.0, 1.0])
    values = np.array([0.0, 1.0])
    assert analysis.piecewise_linear_product_integral(
        x, values, values
    ) == pytest.approx(1.0 / 3.0)


def test_spatial_rsd_uses_exact_piecewise_linear_variance() -> None:
    x = np.array([0.0, 1.0])
    values = np.array([1.0, 2.0])
    mean, standard_deviation, rsd = analysis.spatial_rsd(x, values)
    assert mean == pytest.approx(1.5)
    assert standard_deviation == pytest.approx(np.sqrt(1.0 / 12.0))
    assert rsd == pytest.approx(np.sqrt(1.0 / 12.0) / 1.5)
    _, trapezoidal_standard_deviation, _ = (
        analysis.spatial_rsd_trapezoidal_sensitivity(x, values)
    )
    assert trapezoidal_standard_deviation == pytest.approx(0.5)


def test_current_rsd_inserts_exact_central_80_percent_boundaries() -> None:
    frame = pd.DataFrame(
        {
            "x": [0.0, 0.25, 0.75, 1.0],
            "total_current_density": [-2.0, -1.0, -1.0, -2.0],
        }
    )
    x, current = analysis.clip_profile_interval(frame, "total_current_density", 0.1, 0.9)
    assert x[0] == pytest.approx(0.1)
    assert x[-1] == pytest.approx(0.9)
    _, _, rsd = analysis.spatial_rsd(x, current)
    assert np.isfinite(rsd)


def test_piecewise_clamp_fraction_is_length_weighted() -> None:
    x = np.array([0.0, 0.1, 1.0])
    field = np.array([0.0, 2.0, 2.0])
    measure = analysis._piecewise_condition_measure(x, field, "below", 1.0)
    assert measure == pytest.approx(0.05)


def _synthetic_profile(path_id: str) -> tuple[analysis.ProfileSource, pd.DataFrame]:
    x = np.array([0.0, 0.25, 0.75, 1.0])
    frame = pd.DataFrame(
        {
            "x": x,
            "electron_density": [1.0, 2.0, 2.0, 1.0],
            "mean_electron_energy": [3.0, 2.0, 2.0, 3.0],
            "electric_potential": [0.0, 1.0, 2.0, 3.0],
            "electron_current_density": [-0.8, -0.8, -0.8, -0.8],
            "ion_current_density": [-0.2, -0.2, -0.2, -0.2],
            "total_current_density": [-1.2, -1.0, -1.0, -1.2],
            "excitation_source": [3.0, 2.0, 2.0, 1.0],
            "ionization_source": [2.0, 1.0, 1.0, 0.5],
            "E_over_N": np.array([600.0, 400.0, 400.0, 600.0]) * 1.0e-21,
            "applied_voltage": 200.0,
            "gas_pressure": 13.3322,
        }
    )
    source = analysis.ProfileSource(
        path=Path(f"{path_id}.csv"),
        label=path_id,
        mesh_elements=200,
        path_id=path_id,
    )
    return source, frame


def test_regional_error_shares_partition_full_domain() -> None:
    _, reference = _synthetic_profile("reference")
    external = reference.copy()
    external["electron_density"] = [2.0, 2.0, 2.0, 1.0]
    regional = analysis.regional_error_attribution(
        external, reference, "V m^2"
    )
    density_share = regional.loc[
        regional["quantity"] == "electron_density",
        "squared_error_share",
    ].sum()
    assert density_share == pytest.approx(1.0)
    assert set(regional["region"]) == {
        "cathode_side_10_percent",
        "central_80_percent",
        "anode_side_10_percent",
    }


def test_current_rsd_sensitivity_covers_both_paths_and_requested_windows() -> None:
    external = _synthetic_profile("external_swarm")
    reference = _synthetic_profile("builtin_reference")
    result = analysis.current_rsd_sensitivity([external, reference])
    assert len(result) == 16
    assert set(result["path"]) == {"external_swarm", "builtin_reference"}
    assert set(result["central_retained_percent"]) == {
        50.0,
        60.0,
        70.0,
        80.0,
        85.0,
        90.0,
        95.0,
        100.0,
    }
    assert np.isfinite(
        result["piecewise_linear_exact_relative_standard_deviation"]
    ).all()


def test_regime_metrics_use_length_and_node_fractions_for_both_paths() -> None:
    external = _synthetic_profile("external_swarm")
    reference = _synthetic_profile("builtin_reference")
    metrics, impact, context = analysis.drift_diffusion_regime_metrics(
        [external, reference],
        "V m^2",
        gas_temperature_K=293.15,
        pressure_Pa=13.3322,
    )
    assert len(metrics) == 2
    assert len(impact) == 2 * len(analysis.QUANTITIES)
    assert np.allclose(metrics["node_fraction"].to_numpy(float), 0.5)
    assert metrics["spatial_fraction"].between(0.0, 1.0).all()
    assert context["neutral_number_density_m3"] == pytest.approx(
        13.3322 / (analysis.BOLTZMANN_CONSTANT_J_K * 293.15)
    )


def test_net_flux_speed_ratio_threshold_measure_uses_continuous_inputs() -> None:
    x = np.array([0.0, 0.2, 1.0])
    density = np.full(3, 2.0e15)
    energy = np.full(3, 3.0)
    energy_speed = np.sqrt(
        2.0
        * analysis.ELEMENTARY_CHARGE_C
        * energy
        / analysis.ELECTRON_MASS_KG
    )
    current = (
        0.2
        * analysis.ELEMENTARY_CHARGE_C
        * density
        * energy_speed
    )
    assert analysis._net_flux_ratio_measure_above(
        x, density, current, energy, 0.1
    ) == pytest.approx(1.0)
    assert analysis._net_flux_ratio_measure_above(
        x, density, current, energy, 1.0
    ) == pytest.approx(0.0)
    assert analysis._net_flux_ratio_piecewise_max(
        x, density, current, energy
    ) == pytest.approx(0.2)


def test_fluid_applicability_diagnostics_report_weak_ionization_proxy() -> None:
    profiles = [
        _synthetic_profile("external_swarm"),
        _synthetic_profile("builtin_reference"),
    ]
    for _, frame in profiles:
        frame["electron_density"] *= 1.0e15
        frame["electron_current_density"] *= 1.0e-4
    result = analysis.fluid_applicability_diagnostics(
        profiles,
        gas_temperature_K=293.15,
        pressure_Pa=13.3322,
    )
    assert len(result) == 2
    assert (result["max_electron_to_neutral_ratio"] < 1.0e-5).all()
    assert result[
        "spatial_fraction_net_flux_speed_ratio_above_guide"
    ].between(0.0, 1.0).all()
    assert np.isfinite(result["max_net_flux_speed_ratio"]).all()


def test_positive_column_source_probes_use_model_neutral_density() -> None:
    text = (
        ROOT / "Model" / "maps" / "positive_column_external.yaml"
    ).read_text(encoding="utf-8")
    assert "alpha(plas.ebar)*plas.Nn*abs(plas.Jelx)/e_const" in text
    assert "p0/(k_B_const*293.15[K])" not in text


def test_historical_activation_audit_is_formulation_and_provenance_aware(
    tmp_path: Path,
) -> None:
    bundle = tmp_path / "mixture_0000"
    bundle.mkdir()
    pd.DataFrame(
        {
            "E_over_N_V_m2": [0.0, 1.0],
            "mean_energy_eV": [1.0, 2.0],
        }
    ).to_csv(bundle / "mean_energy_vs_en.csv", index=False)
    pd.DataFrame(
        {
            "mean_energy_eV": [0.0, 3.0],
            "reduced_mobility_m2_V_s_m3": [10.0, 4.0],
        }
    ).to_csv(bundle / "transport_vs_mean_energy.csv", index=False)
    pd.DataFrame(
        {
            "mean_energy_eV": [0.0, 3.0, 0.0, 3.0],
            "process_type": [
                "excitation",
                "excitation",
                "ionization",
                "ionization",
            ],
            "mixture_weighted_reduced_townsend_m2": [0.0, 3.0, 0.0, 6.0],
        }
    ).to_csv(bundle / "rates_vs_mean_energy.csv", index=False)
    archived = pd.DataFrame(
        {
            "E_over_N": [0.0, 1.0],
            "mean_electron_energy": [0.5, 1.0],
            "reduced_mobility": [9.0, 8.0],
            "excitation_townsend": [0.5, 1.0],
            "ionization_townsend": [1.0, 2.0],
        }
    )
    result = analysis.historical_activation_audit(archived, bundle)
    assert len(result) == 4
    mean_energy = result[
        result["quantity"] == "mean_energy_from_E_over_N"
    ].iloc[0]
    assert mean_energy["formulation_activity"] == "inactive_in_archived_LEA"
    assert "not_activation_verification" in mean_energy["audit_scope"]
    assert "stale_class" in str(mean_energy["status"])
