from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from electron_swarm import load_config, run

from product_helpers import base_product_config, write_config, write_moment_table


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
    for column in [
        "meta_physics_level",
        "meta_angular_model",
        "meta_angular_moment_source",
        "meta_ordinary_integral_xs_closure",
        "meta_ionization_source_model",
        "meta_ionization_source_treatment",
        "meta_electron_electron_treatment",
        "meta_transport_definition",
        "meta_magnetic_field_treatment",
        "meta_tail_probability",
        "meta_energy_grid_tail_status",
    ]:
        assert column in summary.columns
    assert not any(col.startswith("meta_capability_") for col in summary.columns)
    mt = summary[summary["solver"] == "multi_term"].iloc[0]
    assert mt["meta_angular_model"] == "isotropic"
    assert mt["meta_angular_moment_source"] == "isotropic_closure"
    assert bool(mt["meta_exact_dcs_based"]) is False
    assert bool(mt["meta_ordinary_integral_xs_closure"]) is True
    assert mt["meta_transport_definition"] == "flux"
    assert not any(col.startswith("meta_multiterm_") for col in summary.columns)
    assert not any(col.startswith("meta_bulk_") for col in summary.columns)
    assert not any(col.startswith("meta_estimated_bulk_") for col in summary.columns)
    assert "meta_hydrodynamic" not in summary.columns
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
    assert "capability_angular_scattering" in plan.columns
    assert "capability_ionization_source" in plan.columns
    assert "effective_ionization_source" in plan.columns
    assert "effective_finite_k" in plan.columns
    assert set(plan["capability_bulk_transport"]) == {"unsupported"}
    assert "effective_rf_field" not in plan.columns


def test_pn_dcs_moment_table_summary_output(tmp_path: Path) -> None:
    table = write_moment_table(tmp_path)
    data = base_product_config(tmp_path, ["multi_term"])
    data["solvers"]["multi_term"]["method"] = "pn_dcs"
    data["physics"]["angular_scattering"] = {
        "model": "moment_table",
        "moment_table": {
            "path": table.as_posix(),
            "format": "csv",
            "provenance": "precomputed_moments",
            "extrapolation": "error",
        },
    }
    cfg = load_config(write_config(tmp_path, data))
    run(cfg, write=True)

    summary = pd.read_csv(tmp_path / "prod_summary.csv")
    [row] = summary.to_dict("records")
    assert row["solver"] == "multi_term"
    assert row["solver_method"] == "pn_dcs"
    assert row["meta_physics_level"] == "table_moments"
    assert row["meta_angular_model"] == "moment_table"
    assert row["meta_angular_moment_source"] == "moment_table"
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
    assert bool(row["meta_electron_electron_transport_stale"]) is False


@pytest.mark.mc
def test_internal_mc_magnetic_summary_output(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "backend": "internal",
        "angular_scattering": "same_as_physics",
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
    run(cfg, write=True)

    summary = pd.read_csv(tmp_path / "prod_summary.csv")
    [row] = summary.to_dict("records")
    assert row["meta_magnetic_field_treatment"] == "boris_lorentz_push"
    assert row["meta_field_integrator"] == "boris"
    assert row["meta_magnetic_field_B_T"] == pytest.approx(0.02)
    assert row["meta_transport_definition"] == "mc_particle_tracking"


def test_same_angular_pn_mc_comparison_summary(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["multi_term", "monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "python_api": "product_helpers:fake_mc_matching",
        "angular_scattering": "same_as_physics",
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
    assert row["angular_model_status"] == "match"
    assert bool(row["same_angular_model"]) is True
    assert row["angular_model"] == "isotropic"
    assert pd.isna(row["angular_model_warning"]) or row["angular_model_warning"] == ""
    assert row["reference_angular_model"] == "isotropic"
    assert row["candidate_angular_model"] == "isotropic"
    assert row["reference_angular_moment_source"] == "isotropic_closure"
    assert row["candidate_angular_moment_source"] == "isotropic_closure"


def test_required_pn_mc_comparison_rejects_unknown_angular_model(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["multi_term", "monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "python_api": "product_helpers:fake_mc_missing",
    }
    data["comparison"] = {
        "enabled": True,
        "reference_solver": "multi_term",
        "candidate_solvers": ["monte_carlo"],
        "compare_eedf": True,
        "required": True,
    }
    cfg = load_config(write_config(tmp_path, data))
    with pytest.raises(ValueError, match="same-angular PN vs MC"):
        run(cfg, write=False)
    with pytest.raises(ValueError, match="same-angular PN vs MC"):
        run(cfg, write=True)

    data["comparison"]["required"] = False
    cfg = load_config(write_config(tmp_path, data))
    run(cfg, write=True)
    summary = pd.read_csv(tmp_path / "prod_summary.csv")
    mc = summary[summary["solver"] == "monte_carlo"].iloc[0]
    assert mc["meta_angular_model"] == "unknown"
    assert mc["meta_angular_moment_source"] == "external_adapter"
    comparison = pd.read_csv(tmp_path / "prod_comparison_summary.csv")
    [row] = comparison.to_dict("records")
    assert row["angular_model_status"] == "unknown"
    assert bool(row["same_angular_model"]) is False
    assert row["angular_model_warning"] == "angular_model_unknown"


def test_monte_carlo_sampler_fallback_is_visible_in_summary(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["physics"]["angular_scattering"] = {
        "model": "momentum_power",
        "higher_moment_closure": "power",
    }
    data["solvers"]["monte_carlo"] = {
        "python_api": "product_helpers:fake_mc_matching",
        "angular_scattering": "same_as_physics",
    }
    data["feature_policy"] = {
        "unsupported": "approximate",
        "degraded": "warn",
        "allow_unsupported_fallback": True,
    }
    cfg = load_config(write_config(tmp_path, data))
    run(cfg, write=True)

    summary = pd.read_csv(tmp_path / "prod_summary.csv")
    [row] = summary.to_dict("records")
    assert row["meta_effective_angular_scattering"].endswith(
        "external_metadata_validation_fallback"
    )


def test_pn_mc_comparison_records_angular_mismatch_warning(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["multi_term", "monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "python_api": "product_helpers:fake_mc_mismatch",
    }
    data["comparison"] = {
        "enabled": True,
        "reference_solver": "multi_term",
        "candidate_solvers": ["monte_carlo"],
        "compare_eedf": True,
        "required": False,
    }
    cfg = load_config(write_config(tmp_path, data))
    run(cfg, write=True)

    comparison = pd.read_csv(tmp_path / "prod_comparison_summary.csv")
    [row] = comparison.to_dict("records")
    assert row["angular_model_status"] == "mismatch"
    assert bool(row["same_angular_model"]) is False
    assert row["angular_model_warning"] == "angular_model_mismatch"
    assert row["angular_model"] == "isotropic->mismatch"
