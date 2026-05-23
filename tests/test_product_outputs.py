from __future__ import annotations

from pathlib import Path

import pandas as pd

from electron_swarm import load_config, run

from product_helpers import base_product_config, write_config


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
        "meta_ordinary_integral_xs_closure",
        "meta_electron_electron_treatment",
        "meta_magnetic_field_treatment",
        "meta_tail_probability",
        "meta_energy_grid_tail_status",
    ]:
        assert column in summary.columns
    assert not any(col.startswith("meta_capability_") for col in summary.columns)
    mt = summary[summary["solver"] == "multi_term"].iloc[0]
    assert mt["meta_angular_model"] == "isotropic"
    assert bool(mt["meta_exact_dcs_based"]) is False
    assert bool(mt["meta_ordinary_integral_xs_closure"]) is True
    assert not any(col.startswith("meta_multiterm_") for col in summary.columns)

    comparison = pd.read_csv(tmp_path / "prod_comparison_summary.csv")
    assert list(comparison["candidate_solver"]) == ["multi_term"]
    assert "mean_energy_eV_relative_difference" in comparison.columns
    assert "eedf_l1_error" in comparison.columns

    rates = pd.read_csv(tmp_path / "prod_rates.csv")
    assert "tail_fraction" in rates.columns

    plan = pd.read_csv(tmp_path / "prod_solver_plan.csv")
    assert "capability_angular_scattering" in plan.columns
    assert "effective_rf_field" not in plan.columns
