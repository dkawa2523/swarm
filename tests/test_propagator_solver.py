from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from electron_swarm import load_config, run
from electron_swarm.orchestration.plan import build_solve_plan
from electron_swarm.solvers.propagator.steady import PropagatorConvergenceError

from product_helpers import base_product_config, write_config


def _propagator_config(
    tmp_path: Path,
    *,
    solvers: list[str] | None = None,
    max_iterations: int = 2000,
) -> dict:
    data = base_product_config(tmp_path, solvers or ["propagator"])
    data["run"]["e_over_n_Td"] = [10.0]
    data["solvers"]["propagator"] = {
        "method": "stationary_response",
        "energy_cells": 64,
        "polar_cells": 8,
        "max_iterations": max_iterations,
        "convergence_tolerance": 1.0e-8,
        "max_memory_mb": 128,
    }
    return data


def test_propagator_schema_plan_execution_and_canonical_outputs(
    tmp_path: Path,
) -> None:
    data = _propagator_config(tmp_path)
    config = load_config(write_config(tmp_path, data))
    [plan] = build_solve_plan(config)
    assert plan.solver == "propagator"
    assert plan.runnable is True
    assert plan.treatment("angular_scattering") == (
        "same_as_physics:isotropic:cell_kernel"
    )

    result = run(config, write=True)
    [case] = result.cases
    assert case.metadata["angular_scattering_treatment"] == (
        "same_as_physics:isotropic:cell_kernel"
    )
    assert case.metadata["angular_scattering_fidelity"] == (
        "integral_xs_closure"
    )
    assert case.metadata["angular_scattering_assumption"] == "isotropic:zero"
    diagnostic = case.diagnostics["propagator"]
    assert diagnostic["converged"] is True
    assert diagnostic["schema"] == "swarm.propagator_diagnostics.v2"
    assert diagnostic["acceleration_scheme"] == (
        "stationary_collision_coupled_shell_response"
    )
    assert diagnostic["steady_algorithm"] == (
        "positive_perron_cone_certified_safeguarded_"
        "secant_brent_number_balance"
    )
    timing = diagnostic["timings_s"]
    assert timing["response_build_blas_threads"] == 1.0
    assert timing["growth_response_system_builds"] >= 1.0
    assert timing["growth_secant_steps"] >= 0.0
    assert diagnostic["perron_branch_selection"] == (
        "staged_largest_real_cone_certificate_with_multiritz_fallback"
    )
    assert diagnostic["reaction_rate_estimator"] == (
        "xs_knot_threshold_partitioned_cell_quadrature"
    )
    assert diagnostic["inelastic_energy_transfer"] == (
        "positive_weak_projection_with_incident_and_daughter_energy_identity"
    )
    assert diagnostic["elastic_frequency_model"] == (
        "separate_cell_integrated_total_and_momentum"
    )
    assert diagnostic["elastic_equilibrium_measure"] == (
        "cell_integrated_maxwell_energy_mass"
    )
    assert diagnostic["shell_coefficient_model"] == (
        "exact_inverse_speed_average_with_bounded_refinement"
    )
    assert diagnostic["energy_domain_strategy"] == (
        "one_shot_core_plus_stretched_tail_to_configured_ceiling"
    )
    assert diagnostic["energy_core_coordinate"] == "sinh_stretched_speed"
    assert diagnostic["energy_core_coordinate_strength"] == pytest.approx(2.5)
    assert diagnostic["operator_residual_L1"] <= 1.0e-8
    assert diagnostic["number_balance_residual"] <= 1.0e-12
    assert diagnostic["growth_number_balance_relative_difference"] <= 1.0e-10
    assert diagnostic["negative_population_mass"] <= 1.0e-14
    assert diagnostic["iterations"] <= 2000
    assert diagnostic["estimated_peak_memory_bytes"] < 128 * 1024 * 1024
    timings = diagnostic["timings_s"]
    assert timings["response_ritz_pairs"] in {1.0, 4.0}
    assert timings["response_cone_compatible_ritz_pairs"] == 1.0
    assert timings["response_perron_real_part_gap"] > 0.0
    assert timings["response_shell_rational_error_estimate"] <= 2.0e-5
    assert timings["response_shell_coefficient_error_estimate"] <= 5.0e-4
    assert timings["response_shell_maximum_coefficient_segments"] <= 128
    assert case.diffusion_L_m2_s is None
    assert case.diffusion_T_m2_s is None
    assert case.reduced_electron_energy_mobility_m2_V_s_m3 is None
    assert case.energy_angle_distribution is not None
    distribution = case.energy_angle_distribution
    assert float(
        np.sum(
            distribution.density_eV_inv
            * distribution.energy_widths_eV[:, None]
            * distribution.mu_widths[None, :]
        )
    ) == pytest.approx(1.0, abs=2.0e-13)
    assert case.metadata["velocity_space_representation"] == (
        "axisymmetric_energy_theta_cells"
    )
    assert case.metadata["transport_components"] == "drift_velocity;mobility"
    assert case.metadata["inelastic_radial_transfer"] == (
        "xs_knot_threshold_partitioned_positive_weak_projection"
    )
    assert case.metadata["elastic_recoil_model"] == (
        "momentum_transfer_driven_finite_temperature_"
        "first_mass_ratio_reversible_sg"
    )
    elastic_loss = diagnostic["elastic_energy_loss"]
    assert elastic_loss["gas_temperature_terms_included"] is True
    assert elastic_loss["neutral_thermal_motion_model"] == (
        "finite_temperature_fokker_planck"
    )

    summary = pd.read_csv(tmp_path / "prod_summary.csv")
    assert pd.isna(summary.loc[0, "diffusion_L_m2_s"])
    assert pd.isna(summary.loc[0, "diffusion_T_m2_s"])
    angle = pd.read_csv(tmp_path / "prod_energy_angle_distribution.csv")
    assert len(angle) == len(distribution.energy_eV) * len(distribution.mu)
    assert float(
        np.sum(
            angle["energy_angle_density_eV_inv"]
            * angle["energy_width_eV"]
            * angle["mu_width"]
        )
    ) == pytest.approx(1.0, rel=1.0e-9)


def test_propagator_maxent_p1_runs_with_explicit_total_and_momentum_xs(
    tmp_path: Path,
) -> None:
    data = _propagator_config(tmp_path)
    source_path = Path(data["cross_sections"]["files"][0]["path"])
    source = pd.read_csv(source_path)
    momentum = source[source["type"] == "momentum"].copy()
    total = momentum.copy()
    total["process"] = "explicit elastic total"
    total["type"] = "elastic_total"
    total["cross_section_m2"] *= 2.0
    path = tmp_path / "maxent_cross_sections.csv"
    pd.concat([total, source], ignore_index=True).to_csv(path, index=False)
    data["cross_sections"]["files"][0]["path"] = path.as_posix()
    data["physics"]["angular_scattering"] = {
        "model": "maxent_p1",
        "higher_moment_closure": "maxent",
    }

    config = load_config(write_config(tmp_path, data, "maxent.yaml"))
    [plan] = build_solve_plan(config)
    assert plan.treatment("angular_scattering") == (
        "same_as_physics:maxent_p1:cell_kernel"
    )
    [case] = run(config, write=False).cases
    assert case.metadata["elastic_collision_xs_role"] == "explicit_total"
    assert case.metadata["elastic_total_xs_process_ids"]
    assert case.metadata["elastic_momentum_xs_process_ids"]
    diagnostic = case.diagnostics["propagator"]
    assert diagnostic["angular_kernel_source"] == (
        "maxent_p1_reversible_cell_joint_exact_p1"
    )
    assert diagnostic["operator_residual_L1"] <= 1.0e-8


def test_propagator_isotropic_closure_is_not_implementation_degradation(
    tmp_path: Path,
) -> None:
    data = _propagator_config(tmp_path)
    data["physics"]["energy_grid_policy"]["adaptive"] = False
    data["feature_policy"]["degraded"] = "fail"

    [plan] = build_solve_plan(load_config(write_config(tmp_path, data)))

    assert plan.runnable
    assert not plan.degraded
    assert plan.fidelity("angular_scattering") == "integral_xs_closure"
    assert plan.assumption("angular_scattering") == "isotropic:zero"


def test_propagator_enabled_zero_magnetic_field_is_runnable(
    tmp_path: Path,
) -> None:
    data = _propagator_config(tmp_path)
    data["physics"]["field"]["magnetic_field"] = {
        "enabled": True,
        "B_T": 0.0,
        "angle_EB_deg": 90.0,
    }
    config = load_config(write_config(tmp_path, data, "zero_magnetic_field.yaml"))

    [plan] = build_solve_plan(config)
    assert plan.runnable
    assert plan.treatment("magnetic_field") == "none"

    [case] = run(config, write=False).cases
    assert case.solver == "propagator"
    assert case.metadata["magnetic_field_treatment"] == "none"


def test_propagator_accepts_number_conserving_perron_root_at_growth_bound(
    tmp_path: Path,
) -> None:
    data = _propagator_config(tmp_path)
    data["run"]["e_over_n_Td"] = [0.05]
    data["cross_sections"] = {
        "format": "csv",
        "high_energy_extrapolation": "error",
        "files": [
            {
                "path": (
                    Path(__file__).resolve().parents[1]
                    / "examples"
                    / "cross_sections"
                    / "argon_magboltz_11_17_gas2_elastic.csv"
                ).as_posix(),
                "species": "Ar",
                "format": "csv",
            }
        ],
    }
    data["physics"]["angular_scattering"] = {
        "model": "maxent_p1",
        "higher_moment_closure": "maxent",
    }

    [case] = run(
        load_config(write_config(tmp_path, data, "conservative_endpoint.yaml")),
        write=False,
    ).cases
    diagnostic = case.diagnostics["propagator"]
    assert diagnostic["converged"] is True
    assert diagnostic["stop_reason"] == "converged"
    assert diagnostic["operator_residual_L1"] <= 1.0e-8
    assert diagnostic["number_balance_residual"] <= 1.0e-12
    assert diagnostic["tail_probability"] <= 1.0e-9
    assert case.metadata["elastic_collision_xs_role"] == "explicit_total"


def test_propagator_mixture_is_density_scaled_and_order_independent(
    tmp_path: Path,
) -> None:
    def mixture_case(
        name: str,
        *,
        pressure_Pa: float,
        reverse_components: bool = False,
    ):
        data = _propagator_config(tmp_path)
        data["run"]["e_over_n_Td"] = [0.1]
        components = [
            {"species": "Ar", "fraction": 0.9, "mass_amu": 39.948},
            {"species": "O2", "fraction": 0.1, "mass_amu": 31.998},
        ]
        data["conditions"] = {
            "gas_temperature_K": 300.0,
            "pressure_Pa": pressure_Pa,
            "gas_mixture": (
                list(reversed(components))
                if reverse_components
                else components
            ),
        }
        data["cross_sections"] = {
            "format": "csv",
            "high_energy_extrapolation": "hold",
            "files": [
                {
                    "path": (
                        Path(__file__).resolve().parents[1]
                        / "examples"
                        / "cross_sections"
                        / "ar_o2_minimal.csv"
                    ).as_posix(),
                    "format": "csv",
                }
            ],
        }
        data["physics"]["ionization"] = {
            "energy_sharing": "primary_secondary",
            "secondary_electron_energy_eV": 0.5,
        }
        data["physics"]["energy_grid_policy"][
            "threshold_refinement"
        ] = True
        return run(
            load_config(write_config(tmp_path, data, name)),
            write=False,
        ).cases[0]

    baseline = mixture_case("mixture_base.yaml", pressure_Pa=13.3)
    dense = mixture_case("mixture_dense.yaml", pressure_Pa=133.0)
    reordered = mixture_case(
        "mixture_reordered.yaml",
        pressure_Pa=13.3,
        reverse_components=True,
    )
    base_growth = baseline.diagnostics["propagator"]["growth_frequency_s_inv"]
    dense_growth = dense.diagnostics["propagator"]["growth_frequency_s_inv"]

    assert dense.energy_eV == pytest.approx(baseline.energy_eV, abs=0.0)
    assert dense.eedf == pytest.approx(baseline.eedf, rel=2.0e-8, abs=2.0e-12)
    assert dense.mean_energy_eV == pytest.approx(
        baseline.mean_energy_eV,
        rel=2.0e-8,
    )
    assert dense.drift_velocity_m_s == pytest.approx(
        baseline.drift_velocity_m_s,
        rel=2.0e-8,
    )
    assert dense_growth == pytest.approx(10.0 * base_growth, rel=2.0e-8)
    assert reordered.eedf == pytest.approx(baseline.eedf, abs=2.0e-13)
    assert reordered.mean_energy_eV == pytest.approx(
        baseline.mean_energy_eV,
        abs=2.0e-13,
    )
    assert reordered.drift_velocity_m_s == pytest.approx(
        baseline.drift_velocity_m_s,
        abs=2.0e-10,
    )


def test_obsolete_propagator_method_name_is_rejected(tmp_path: Path) -> None:
    data = _propagator_config(tmp_path)
    data["solvers"]["propagator"]["method"] = "relaxation_accelerated"
    with pytest.raises(ValueError, match="solvers.propagator.method"):
        load_config(write_config(tmp_path, data, "obsolete_method.yaml"))


def test_comparison_records_unavailable_propagator_diffusion(
    tmp_path: Path,
) -> None:
    data = _propagator_config(tmp_path, solvers=["two_term", "propagator"])
    data["comparison"] = {
        "enabled": True,
        "reference_solver": "two_term",
        "candidate_solvers": ["propagator"],
        "compare_eedf": True,
        "required": True,
    }
    result = run(load_config(write_config(tmp_path, data)), write=True)
    [row] = result.metadata["comparison_summary_rows"]
    assert row["diffusion_L_m2_s_status"] == "candidate_unavailable"
    assert row["diffusion_L_relative_difference"] is None
    assert row["mean_energy_eV_status"] == "available"
    assert np.isfinite(float(row["eedf_l1_error"]))


@pytest.mark.parametrize("feature", ["magnetic_field", "electron_electron"])
def test_propagator_unsupported_physics_obeys_feature_policy(
    tmp_path: Path,
    feature: str,
) -> None:
    data = _propagator_config(tmp_path)
    if feature == "magnetic_field":
        data["physics"]["field"]["magnetic_field"] = {
            "enabled": True,
            "B_T": 0.01,
            "angle_EB_deg": 90.0,
        }
    else:
        data["physics"]["electron_electron"] = {
            "enabled": True,
            "model": "fp_energy",
        }
    with pytest.raises(ValueError, match=feature.replace("_", " ")):
        build_solve_plan(load_config(write_config(tmp_path, data, "fail.yaml")))

    data["feature_policy"]["unsupported"] = "skip_solver"
    [item] = build_solve_plan(
        load_config(write_config(tmp_path, data, "skip.yaml"))
    )
    assert item.skipped is True
    assert item.runnable is False
    assert item.treatment(feature) in {"unsupported", "skipped"}


def test_propagator_hard_iteration_and_memory_limits_fail_closed(
    tmp_path: Path,
) -> None:
    data = _propagator_config(tmp_path, max_iterations=10)
    with pytest.raises(PropagatorConvergenceError) as error:
        run(load_config(write_config(tmp_path, data, "iteration.yaml")), write=False)
    assert error.value.diagnostics.iterations == 10
    assert error.value.diagnostics.stop_reason == "maximum_response_applications"

    memory = _propagator_config(tmp_path)
    memory["solvers"]["propagator"].update(
        {
            "energy_cells": 4096,
            "polar_cells": 360,
            "max_memory_mb": 128,
        }
    )
    with pytest.raises(MemoryError, match="preflight estimate"):
        run(load_config(write_config(tmp_path, memory, "memory.yaml")), write=False)


@pytest.mark.parametrize("field", ["energy_cells", "polar_cells", "max_iterations", "max_memory_mb"])
def test_propagator_integer_controls_reject_nonintegers(
    tmp_path: Path,
    field: str,
) -> None:
    data = _propagator_config(tmp_path)
    data["solvers"]["propagator"][field] = 64.5
    with pytest.raises(ValueError, match=f"solvers.propagator.{field} must be an integer"):
        load_config(write_config(tmp_path, data))
