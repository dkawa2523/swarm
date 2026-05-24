from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from electron_swarm import load_config, run

from product_helpers import base_product_config, write_config, write_moment_table


def test_two_term_smoke_run(tmp_path: Path) -> None:
    cfg = load_config(write_config(tmp_path, base_product_config(tmp_path, ["two_term"])))
    result = run(cfg, write=False)
    [case] = result.cases
    assert case.solver == "two_term"
    assert case.schema_version == "2"
    assert case.mean_energy_eV > 0.0
    assert case.metadata["energy_grid_tail_status"] in {"ok", "warning", "insufficient"}
    assert case.metadata["tail_probability"] >= 0.0
    assert np.all(np.isfinite(case.eedf))


def test_multi_term_surrogate_smoke_and_minimum_metadata(tmp_path: Path) -> None:
    cfg = load_config(write_config(tmp_path, base_product_config(tmp_path, ["multi_term"])))
    result = run(cfg, write=False)
    [case] = result.cases
    assert case.solver == "multi_term"
    assert case.metadata["physics_level"] == "surrogate"
    assert case.metadata["angular_model"] == "isotropic"
    assert case.metadata["angular_moment_source"] == "isotropic_closure"
    assert case.metadata["solver_method"] == "pn_closure_surrogate"
    assert case.metadata["lmax"] == 3
    assert case.metadata["exact_dcs_based"] is False
    assert case.metadata["ordinary_integral_xs_closure"] is True
    assert case.metadata["direct_pn_operator"] is False
    assert case.metadata["energy_grid_tail_status"] in {"ok", "warning", "insufficient"}
    assert case.mean_energy_eV > 0.0


@pytest.mark.parametrize(
    ("model", "closure"),
    [
        ("momentum_power", "power"),
        ("maxent_p1", "maxent"),
    ],
)
def test_multi_term_uses_selected_angular_closure(
    tmp_path: Path, model: str, closure: str
) -> None:
    data = base_product_config(tmp_path, ["multi_term"])
    data["physics"]["angular_scattering"] = {
        "model": model,
        "higher_moment_closure": closure,
    }
    cfg = load_config(write_config(tmp_path, data))
    result = run(cfg, write=False)
    [case] = result.cases

    assert case.metadata["angular_model"] == model
    assert case.metadata["angular_moment_source"] == "ordinary_integral_xs_closure"
    assert case.metadata["lmax"] == 3
    assert case.metadata["direct_pn_operator"] is False
    assert case.metadata["exact_dcs_based"] is False
    assert case.metadata["ordinary_integral_xs_closure"] is True


def test_multi_term_pn_dcs_uses_moment_table_metadata(tmp_path: Path) -> None:
    table = write_moment_table(tmp_path)
    data = base_product_config(tmp_path, ["multi_term"])
    data["solvers"]["multi_term"]["method"] = "pn_dcs"
    data["physics"]["angular_scattering"] = {
        "model": "moment_table",
        "moment_table": {
            "path": table.as_posix(),
            "format": "csv",
            "provenance": "dcs_derived",
            "extrapolation": "error",
        },
    }
    cfg = load_config(write_config(tmp_path, data))
    result = run(cfg, write=False)
    [case] = result.cases

    assert case.solver == "multi_term"
    assert case.metadata["solver_method"] == "pn_dcs"
    assert case.metadata["physics_level"] == "table_moments"
    assert case.metadata["angular_model"] == "moment_table"
    assert case.metadata["angular_moment_source"] == "moment_table"
    assert case.metadata["exact_dcs_based"] is True
    assert case.metadata["ordinary_integral_xs_closure"] is False
    assert case.metadata["direct_pn_operator"] is False
    assert np.isfinite(case.mean_energy_eV)


def test_multi_term_pn_dcs_moment_table_affects_solution(tmp_path: Path) -> None:
    cases = []
    for name, m1 in [("forward", 0.9), ("backward", -0.9)]:
        table = write_moment_table(tmp_path, name=f"{name}.csv", m1=m1)
        data = base_product_config(tmp_path, ["multi_term"])
        data["solvers"]["multi_term"]["method"] = "pn_dcs"
        data["physics"]["angular_scattering"] = {
            "model": "moment_table",
            "moment_table": {
                "path": table.as_posix(),
                "format": "csv",
                "provenance": "dcs_derived",
                "extrapolation": "error",
            },
        }
        cfg = load_config(write_config(tmp_path, data, name=f"{name}.yaml"))
        cases.append(run(cfg, write=False).cases[0])

    forward, backward = cases
    assert forward.metadata["angular_moment_source"] == "moment_table"
    assert backward.metadata["angular_moment_source"] == "moment_table"
    assert not np.isclose(forward.mean_energy_eV, backward.mean_energy_eV)
    assert not np.isclose(forward.drift_velocity_m_s, backward.drift_velocity_m_s)


def test_electron_electron_postprocess_marks_transport_stale(tmp_path: Path) -> None:
    baseline_cfg = load_config(
        write_config(
            tmp_path,
            base_product_config(tmp_path, ["two_term"]),
            name="baseline.yaml",
        )
    )
    baseline = run(baseline_cfg, write=False).cases[0]

    data = base_product_config(tmp_path, ["two_term"])
    data["physics"]["electron_electron"] = {
        "enabled": True,
        "model": "relaxation_postprocess",
        "relaxation_fraction": 0.5,
        "conserve_mean_energy": True,
    }
    cfg = load_config(write_config(tmp_path, data, name="ee.yaml"))
    result = run(cfg, write=False)
    [case] = result.cases
    assert not np.allclose(case.eedf, baseline.eedf)
    assert case.metadata["electron_electron_treatment"] == "relaxation_postprocess"
    assert case.metadata["electron_electron_affects_eedf"] is True
    assert case.metadata["electron_electron_affects_rates"] is True
    assert case.metadata["electron_electron_affects_transport"] is False
    assert case.metadata["electron_electron_transport_stale"] is True
    assert case.drift_velocity_m_s == pytest.approx(baseline.drift_velocity_m_s)
    assert case.diffusion_L_m2_s == pytest.approx(baseline.diffusion_L_m2_s)


def test_two_term_fp_energy_runs_inside_solver_and_recomputes_transport(
    tmp_path: Path,
) -> None:
    baseline_cfg = load_config(
        write_config(
            tmp_path,
            base_product_config(tmp_path, ["two_term"]),
            name="baseline_fp_two.yaml",
        )
    )
    baseline = run(baseline_cfg, write=False).cases[0]

    data = base_product_config(tmp_path, ["two_term"])
    data["physics"]["electron_electron"] = {
        "enabled": True,
        "model": "fp_energy",
        "strength_model": "simple_relaxation",
        "relaxation_fraction": 0.5,
        "conserve_mean_energy": False,
        "fallback_temperature_eV": 0.8,
    }
    cfg = load_config(write_config(tmp_path, data, name="fp_two.yaml"))
    [case] = run(cfg, write=False).cases

    assert not np.allclose(case.eedf, baseline.eedf)
    assert case.metadata["electron_electron_treatment"] == "fp_energy"
    assert case.metadata["electron_electron_affects_eedf"] is True
    assert case.metadata["electron_electron_affects_rates"] is True
    assert case.metadata["electron_electron_affects_transport"] is True
    assert case.metadata["electron_electron_transport_stale"] is False
    assert case.metadata["electron_electron_operator_scope"] == "f0_energy_only"


def test_multi_term_fp_energy_marks_transport_stale(tmp_path: Path) -> None:
    baseline_cfg = load_config(
        write_config(
            tmp_path,
            base_product_config(tmp_path, ["multi_term"]),
            name="baseline_fp_multi.yaml",
        )
    )
    baseline = run(baseline_cfg, write=False).cases[0]

    data = base_product_config(tmp_path, ["multi_term"])
    data["physics"]["electron_electron"] = {
        "enabled": True,
        "model": "fp_energy",
        "strength_model": "simple_relaxation",
        "relaxation_fraction": 0.5,
        "conserve_mean_energy": False,
        "fallback_temperature_eV": 0.8,
    }
    cfg = load_config(write_config(tmp_path, data, name="fp_multi.yaml"))
    [case] = run(cfg, write=False).cases

    assert not np.allclose(case.eedf, baseline.eedf)
    assert case.metadata["electron_electron_treatment"] == "fp_energy"
    assert case.metadata["electron_electron_affects_eedf"] is True
    assert case.metadata["electron_electron_affects_rates"] is True
    assert case.metadata["electron_electron_affects_transport"] is False
    assert case.metadata["electron_electron_transport_stale"] is True
    assert case.drift_velocity_m_s == pytest.approx(baseline.drift_velocity_m_s)
    assert case.diffusion_L_m2_s == pytest.approx(baseline.diffusion_L_m2_s)


def test_monte_carlo_same_as_physics_validates_reported_angular_metadata(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "python_api": "product_helpers:fake_mc_matching",
        "angular_scattering": "same_as_physics",
    }
    cfg = load_config(write_config(tmp_path, data))
    result = run(cfg, write=False)
    [case] = result.cases

    assert case.solver == "monte_carlo"
    assert case.metadata["angular_model"] == "isotropic"
    assert case.metadata["angular_moment_source"] == "isotropic_closure"
    assert case.metadata["ordinary_integral_xs_closure"] is True


def test_monte_carlo_same_as_physics_rejects_missing_or_mismatched_metadata(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "python_api": "product_helpers:fake_mc_missing",
        "angular_scattering": "same_as_physics",
    }
    with pytest.raises(ValueError, match="requires angular metadata"):
        run(load_config(write_config(tmp_path, data)), write=False)

    data["solvers"]["monte_carlo"]["python_api"] = "product_helpers:fake_mc_mismatch"
    with pytest.raises(ValueError, match="angular metadata mismatch"):
        run(load_config(write_config(tmp_path, data)), write=False)


def test_monte_carlo_external_missing_metadata_records_unknown(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "python_api": "product_helpers:fake_mc_missing",
    }
    cfg = load_config(write_config(tmp_path, data))
    result = run(cfg, write=False)
    [case] = result.cases

    assert case.metadata["angular_model"] == "unknown"
    assert case.metadata["angular_moment_source"] == "external_adapter"
    assert case.metadata["ordinary_integral_xs_closure"] is False


@pytest.mark.mc
def test_internal_monte_carlo_magnetic_smoke_run(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "backend": "internal",
        "angular_scattering": "same_as_physics",
        "particles": 16,
        "max_collisions": 8,
        "seed": 17,
    }
    data["physics"]["field"]["magnetic_field"] = {
        "enabled": True,
        "B_T": 0.01,
        "angle_EB_deg": 90.0,
    }
    cfg = load_config(write_config(tmp_path, data))
    [case] = run(cfg, write=False).cases

    assert case.solver == "monte_carlo"
    assert np.isfinite(case.mean_energy_eV)
    assert np.all(np.isfinite(case.eedf))
    assert case.metadata["magnetic_field_treatment"] == "boris_lorentz_push"
    assert case.metadata["field_integrator"] == "boris"
    assert case.metadata["magnetic_field_B_T"] == pytest.approx(0.01)
    assert case.metadata["magnetic_field_Bx_T"] == pytest.approx(0.01)
    assert case.metadata["magnetic_field_Bz_T"] == pytest.approx(0.0, abs=1.0e-14)
    assert case.metadata["reaction_rates_source"] == "eedf_convolution"
    assert case.metadata["attachment_trajectory_treatment"] == "rate_convolution_only"


@pytest.mark.mc
def test_internal_monte_carlo_zero_b_matches_disabled_with_same_seed(
    tmp_path: Path,
) -> None:
    base = base_product_config(tmp_path, ["monte_carlo"])
    base["solvers"]["monte_carlo"] = {
        "backend": "internal",
        "angular_scattering": "same_as_physics",
        "particles": 16,
        "max_collisions": 8,
        "seed": 19,
    }
    disabled = run(
        load_config(write_config(tmp_path, base, "mc_no_b.yaml")), write=False
    ).cases[0]

    enabled = base_product_config(tmp_path, ["monte_carlo"])
    enabled["solvers"]["monte_carlo"] = dict(base["solvers"]["monte_carlo"])
    enabled["physics"]["field"]["magnetic_field"] = {
        "enabled": True,
        "B_T": 0.0,
        "angle_EB_deg": 90.0,
    }
    zero_b = run(
        load_config(write_config(tmp_path, enabled, "mc_zero_b.yaml")), write=False
    ).cases[0]

    assert zero_b.mean_energy_eV == pytest.approx(disabled.mean_energy_eV)
    assert zero_b.drift_velocity_m_s == pytest.approx(disabled.drift_velocity_m_s)
    assert np.allclose(zero_b.eedf, disabled.eedf)
    assert disabled.metadata["field_integrator"] == "none"
    assert zero_b.metadata["field_integrator"] == "boris"


def test_external_monte_carlo_magnetic_requires_reported_metadata(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "python_api": "product_helpers:fake_mc_missing",
    }
    data["physics"]["field"]["magnetic_field"] = {
        "enabled": True,
        "B_T": 0.01,
        "angle_EB_deg": 90.0,
    }
    with pytest.raises(ValueError, match="magnetic_field requires MC output metadata"):
        run(load_config(write_config(tmp_path, data)), write=False)

    data["solvers"]["monte_carlo"]["python_api"] = "product_helpers:fake_mc_magnetic_matching"
    [case] = run(load_config(write_config(tmp_path, data)), write=False).cases
    assert case.metadata["magnetic_field_treatment"] == "external_lorentz_push"

    data["solvers"]["monte_carlo"]["python_api"] = "product_helpers:fake_mc_magnetic_mismatch"
    with pytest.raises(ValueError, match="magnetic_field metadata mismatch"):
        run(load_config(write_config(tmp_path, data)), write=False)
