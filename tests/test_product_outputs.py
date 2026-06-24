from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from electron_swarm import load_config, run
from electron_swarm.core.result_metadata import PRODUCT_CASE_METADATA_KEYS
from electron_swarm.io.writers import (
    SOLVER_PLAN_COLUMNS,
    SUMMARY_COLUMNS,
    SUMMARY_METADATA_COLUMNS,
)

from product_helpers import base_product_config, write_config, write_moment_table


def _expected_summary_columns() -> list[str]:
    return SUMMARY_COLUMNS + [f"meta_{key}" for key in SUMMARY_METADATA_COLUMNS]


def test_canonical_outputs_and_comparison_summary(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["two_term", "multi_term"])
    data["comparison"] = {
        "enabled": True,
        "reference_solver": "two_term",
        "candidate_solvers": ["multi_term"],
        "compare_eedf": True,
        "required": True,
    }
    cfg = load_config(write_config(tmp_path, data))
    result = run(cfg, write=True)

    assert {case.solver for case in result.cases} == {"two_term", "multi_term"}
    for case in result.cases:
        assert set(case.metadata) <= PRODUCT_CASE_METADATA_KEYS
    assert (tmp_path / "prod_summary.csv").exists()
    assert (tmp_path / "prod_rates.csv").exists()
    assert (tmp_path / "prod_eedf.csv").exists()
    assert (tmp_path / "prod_solver_plan.csv").exists()
    assert (tmp_path / "prod_comparison_summary.csv").exists()
    assert set(result.metadata["output_paths"]) == {
        "summary_csv",
        "eedf_csv",
        "rates_csv",
        "solver_plan_csv",
        "comparison_summary_csv",
    }

    assert not (tmp_path / "prod_comparison_rates.csv").exists()
    assert not (tmp_path / "prod_comparison_capabilities.csv").exists()
    assert not list(tmp_path.glob("prod_*.png"))
    assert not (tmp_path / "summary.csv").exists()
    assert not (tmp_path / "summary_boltzmann.csv").exists()
    assert not (tmp_path / "summary_multiterm.csv").exists()

    summary = pd.read_csv(tmp_path / "prod_summary.csv")
    expected_meta_columns = {
        "meta_angular_model",
        "meta_angular_moment_source",
        "meta_moment_table_provenance",
        "meta_exact_dcs_based",
        "meta_ordinary_integral_xs_closure",
        "meta_lmax",
        "meta_ionization_source_treatment",
        "meta_electron_electron_treatment",
        "meta_electron_electron_transport_stale",
        "meta_magnetic_field_treatment",
        "meta_tail_refinement_treatment",
        "meta_transport_definition",
    }
    for column in expected_meta_columns:
        assert column in summary.columns
    assert {c for c in summary.columns if c.startswith("meta_")} == expected_meta_columns
    eedf = pd.read_csv(tmp_path / "prod_eedf.csv")
    assert "energy_width_eV" in eedf.columns
    assert "sample_count" in eedf.columns
    assert "effective_sample_count" in eedf.columns
    assert "relative_standard_error" in eedf.columns
    for _, group in eedf.groupby(["solver", "case_id"], sort=False):
        assert (group["energy_width_eV"] > 0.0).all()
        assert float((group["eedf"] * group["energy_width_eV"]).sum()) == pytest.approx(
            1.0
        )
    assert not any(col.startswith("meta_capability_") for col in summary.columns)
    mt = summary[summary["solver"] == "multi_term"].iloc[0]
    assert mt["meta_angular_moment_source"] == "isotropic_closure"
    assert bool(mt["meta_exact_dcs_based"]) is False
    assert bool(mt["meta_ordinary_integral_xs_closure"]) is True
    assert mt["meta_transport_definition"] == "f0_gradient_reconstruction"
    assert not any(col.startswith("meta_multiterm_") for col in summary.columns)
    assert not any(col.startswith("meta_bulk_") for col in summary.columns)
    assert not any(col.startswith("meta_estimated_bulk_") for col in summary.columns)
    for removed in [
        "meta_pn_residual",
        "meta_negative_mass_fraction",
        "meta_lmax1_regression_target",
        "meta_effective_angular_scattering",
        "meta_hydrodynamic",
        "meta_field_integrator",
        "meta_tail_rate_fraction_max",
        "meta_mc_seed",
        "meta_mc_particles",
        "meta_mc_population_model",
        "meta_mc_histogram_samples",
        "meta_mc_tail_uncertainty_status",
        "meta_mc_tail_comparison_status",
        "meta_mc_energy_balance_status",
    ]:
        assert removed not in summary.columns
    assert "schema_version" not in summary.columns
    assert "reduced_mobility_m2_V_s_m3" not in summary.columns
    assert "reduced_diffusion_L_m2_s_m3" not in summary.columns
    assert "reduced_diffusion_T_m2_s_m3" not in summary.columns

    comparison = pd.read_csv(tmp_path / "prod_comparison_summary.csv")
    assert list(comparison["candidate_solver"]) == ["multi_term"]
    assert "mean_energy_eV_relative_difference" in comparison.columns
    assert "eedf_l1_error" in comparison.columns

    rates = pd.read_csv(tmp_path / "prod_rates.csv")
    assert "tail_fraction" in rates.columns

    plan = pd.read_csv(tmp_path / "prod_solver_plan.csv")
    assert "effective_ionization_source" in plan.columns
    assert "effective_finite_k" in plan.columns
    assert not any(col.startswith("capability_") for col in plan.columns)
    assert "effective_rf_field" in plan.columns


def test_all_solvers_skipped_writes_canonical_empty_summary(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["two_term", "multi_term"])
    data["physics"]["field"]["magnetic_field"] = {
        "enabled": True,
        "B_T": 0.01,
        "angle_EB_deg": 90.0,
    }
    data["feature_policy"]["unsupported"] = "skip_solver"
    cfg = load_config(write_config(tmp_path, data))
    result = run(cfg, write=True)

    assert result.cases == []
    summary = pd.read_csv(tmp_path / "prod_summary.csv")
    assert summary.empty
    assert list(summary.columns) == _expected_summary_columns()

    eedf = pd.read_csv(tmp_path / "prod_eedf.csv")
    rates = pd.read_csv(tmp_path / "prod_rates.csv")
    assert eedf.empty
    assert rates.empty

    plan = pd.read_csv(tmp_path / "prod_solver_plan.csv")
    assert list(plan.columns) == SOLVER_PLAN_COLUMNS
    assert set(plan["solver"]) == {"two_term", "multi_term"}
    assert plan["skipped"].astype(str).str.lower().eq("true").all()


def test_all_solvers_disabled_writes_canonical_solver_plan(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["two_term"])
    data["run"]["solvers"] = [{"id": "two_term", "enabled": False}]
    cfg = load_config(write_config(tmp_path, data))
    result = run(cfg, write=True)

    assert result.cases == []
    plan = pd.read_csv(tmp_path / "prod_solver_plan.csv")
    assert list(plan.columns) == SOLVER_PLAN_COLUMNS
    assert len(plan) == 1
    [row] = plan.to_dict("records")
    assert row["solver"] == "two_term"
    assert bool(row["skipped"]) is True
    assert row["skip_reason"] == "solver disabled"
    for column in SOLVER_PLAN_COLUMNS:
        if column.startswith("effective_"):
            assert pd.isna(row[column])


def test_pn_dcs_moment_table_summary_records_table_moment_path(
    tmp_path: Path,
) -> None:
    table = write_moment_table(tmp_path)
    data = base_product_config(tmp_path, ["multi_term"])
    data["solvers"]["multi_term"]["method"] = "pn_dcs"
    data["physics"]["angular_scattering"] = {
        "model": "moment_table",
        "moment_table": {
            "path": table.as_posix(),
            "format": "normalized_legendre_moments",
            "provenance": "model_derived",
            "extrapolation": "error",
        },
    }
    cfg = load_config(write_config(tmp_path, data))
    run(cfg, write=True)

    summary = pd.read_csv(tmp_path / "prod_summary.csv")
    [row] = summary.to_dict("records")
    assert row["solver_method"] == "pn_dcs"
    assert row["meta_angular_moment_source"] == "moment_table"
    assert row["meta_moment_table_provenance"] == "model_derived"
    assert bool(row["meta_exact_dcs_based"]) is False
    assert bool(row["meta_ordinary_integral_xs_closure"]) is False


def test_fp_energy_summary_output_records_product_treatment(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["two_term"])
    data["physics"]["electron_electron"] = {
        "enabled": True,
        "model": "fp_energy",
        "strength_model": "simple_relaxation",
        "relaxation_fraction": 0.4,
        "conserve_mean_energy": False,
        "fallback_temperature_eV": 0.8,
    }
    cfg = load_config(write_config(tmp_path, data))
    run(cfg, write=True)

    summary = pd.read_csv(tmp_path / "prod_summary.csv")
    [row] = summary.to_dict("records")
    assert row["meta_electron_electron_treatment"] == "fp_energy"
    assert bool(row["meta_electron_electron_transport_stale"]) is True


@pytest.mark.mc
def test_internal_mc_magnetic_summary_output(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "particles": 12,
        "max_collisions": 6,
        "seed": 23,
    }
    data["physics"]["field"]["magnetic_field"] = {
        "enabled": True,
        "B_T": 0.02,
        "angle_EB_deg": 90.0,
    }
    cfg = load_config(write_config(tmp_path, data))
    result = run(cfg, write=True)
    [case] = result.cases

    summary = pd.read_csv(tmp_path / "prod_summary.csv")
    [row] = summary.to_dict("records")
    assert row["meta_magnetic_field_treatment"] == "boris_lorentz_push"
    assert row["meta_transport_definition"] == "mc_flux_particle_tracking_fixed_population"
    assert "meta_mc_seed" not in summary.columns
    assert "meta_field_integrator" not in summary.columns
    assert "meta_mc_orbit_substeps" not in summary.columns
    assert "mc_run" not in case.metadata
    assert "mc_tail_audit" not in case.metadata
    assert "_dev_internal_monte_carlo_audit" not in case.metadata
    assert "internal_monte_carlo_audit" not in case.diagnostics

    eedf = pd.read_csv(tmp_path / "prod_eedf.csv")
    assert "energy_width_eV" in eedf.columns
    assert (eedf["energy_width_eV"] > 0.0).all()
    assert float((eedf["eedf"] * eedf["energy_width_eV"]).sum()) == pytest.approx(
        1.0
    )
    assert eedf["sample_count"].notna().any()
    assert eedf["effective_sample_count"].notna().any()
    positive = eedf["sample_count"] > 0
    assert (eedf.loc[positive, "relative_standard_error"] > 0.0).all()


@pytest.mark.mc
def test_internal_mc_weighted_branching_summary_metadata(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["run"]["e_over_n_Td"] = [600.0]
    data["cross_sections"]["high_energy_extrapolation"] = "hold"
    data["solvers"]["monte_carlo"] = {
        "population_model": "weighted_branching",
        "particles": 8,
        "max_collisions": 30,
        "seed": 41,
    }
    cfg = load_config(write_config(tmp_path, data))
    result = run(cfg, write=True)
    [case] = result.cases

    summary = pd.read_csv(tmp_path / "prod_summary.csv")
    [row] = summary.to_dict("records")
    assert (
        row["meta_transport_definition"]
        == "mc_flux_particle_tracking_weighted_growth_population"
    )
    assert "_dev_internal_monte_carlo_audit" not in case.metadata
    assert "internal_monte_carlo_audit" not in case.diagnostics


def test_matching_angular_pn_mc_comparison_summary(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["multi_term", "monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "particles": 8,
        "max_collisions": 4,
        "seed": 13,
    }
    data["comparison"] = {
        "enabled": True,
        "reference_solver": "multi_term",
        "candidate_solvers": ["monte_carlo"],
        "compare_eedf": True,
        "required": True,
    }
    cfg = load_config(write_config(tmp_path, data))
    run(cfg, write=True)

    comparison = pd.read_csv(tmp_path / "prod_comparison_summary.csv")
    [row] = comparison.to_dict("records")
    assert set(comparison.columns) == {
        "case_id",
        "E_over_N_Td",
        "reference_solver",
        "candidate_solver",
        "angular_model_status",
        "mean_energy_eV_relative_difference",
        "drift_velocity_relative_difference",
        "mobility_relative_difference",
        "diffusion_L_relative_difference",
        "net_ionization_frequency_relative_difference",
        "eedf_l1_error",
    }
    assert row["angular_model_status"] == "match"


def test_monte_carlo_sampler_fallback_is_rejected(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["physics"]["angular_scattering"] = {
        "model": "momentum_power",
        "higher_moment_closure": "power",
    }
    data["solvers"]["monte_carlo"] = {}
    data["feature_policy"] = {
        "unsupported": "approximate",
        "degraded": "record",
        "allow_unsupported_fallback": True,
    }
    with pytest.raises(ValueError, match="allow_unsupported_fallback"):
        load_config(write_config(tmp_path, data))
