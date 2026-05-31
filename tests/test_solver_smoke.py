from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from electron_swarm import load_config, run
from electron_swarm.core.cross_sections import ProcessType, load_cross_sections
from electron_swarm.diagnostics.eedf_compare import compare_eedf_cases
from electron_swarm.solvers.kinetic import (
    assemble_native_operator_blocks,
    cell_edges_from_centers,
    compute_rates_from_eedf,
    eepf_from_eedf,
    make_two_term_energy_grid,
    mean_energy_from_eedf,
    negative_mass_fraction,
    normalize_eedf,
)
from electron_swarm.solvers.multi_term.direct import build_higher_l_collision_damping
from electron_swarm.solvers.internal_monte_carlo import (
    _EnergyAudit,
    _mc_bin_relative_standard_error,
    _mc_effective_bin_counts,
    _mc_energy_edges,
    _mc_tail_uncertainty_metadata,
    _post_reaction_energy,
    _post_reaction_outcome,
    _validate_trial_collision_frequency,
)
from electron_swarm.solvers.two_term import TwoTermSolver

from product_helpers import ROOT, base_product_config, write_config, write_moment_table


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


def test_multi_term_default_direct_smoke_and_minimum_metadata(tmp_path: Path) -> None:
    cfg = load_config(write_config(tmp_path, base_product_config(tmp_path, ["multi_term"])))
    result = run(cfg, write=False)
    [case] = result.cases
    assert case.solver == "multi_term"
    assert case.metadata["physics_level"] == "angular_closure"
    assert case.metadata["angular_model"] == "isotropic"
    assert case.metadata["angular_moment_source"] == "isotropic_closure"
    assert case.metadata["solver_method"] == "pn_closure_direct"
    assert case.metadata["lmax"] == 3
    assert case.metadata["exact_dcs_based"] is False
    assert case.metadata["ordinary_integral_xs_closure"] is True
    assert case.metadata["direct_pn_operator"] is True
    assert case.metadata["energy_grid_tail_status"] in {"ok", "warning", "insufficient"}
    assert case.mean_energy_eV > 0.0


def test_multi_term_pn_closure_direct_lmax1_regresses_to_two_term(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["two_term", "multi_term"])
    data["solvers"]["multi_term"]["method"] = "pn_closure_direct"
    data["solvers"]["multi_term"]["lmax"] = 1
    cfg = load_config(write_config(tmp_path, data, name="direct_regression.yaml"))
    two_term, direct = run(cfg, write=False).cases

    assert direct.metadata["solver_method"] == "pn_closure_direct"
    assert direct.metadata["physics_level"] == "angular_closure"
    assert direct.metadata["direct_pn_operator"] is True
    assert direct.metadata["exact_dcs_based"] is False
    assert direct.metadata["ordinary_integral_xs_closure"] is True
    assert direct.metadata["lmax1_regression_target"] == "two_term"
    assert direct.metadata["negative_mass_fraction"] < 1.0e-8
    assert "direct_pn_iterations" not in direct.metadata
    assert "direct_pn_growth_frequency_s-1" not in direct.metadata
    assert np.all(np.isfinite(direct.eedf))
    widths = np.diff(
        np.r_[
            max(0.0, direct.energy_eV[0] - 0.5 * (direct.energy_eV[1] - direct.energy_eV[0])),
            0.5 * (direct.energy_eV[:-1] + direct.energy_eV[1:]),
            direct.energy_eV[-1] + 0.5 * (direct.energy_eV[-1] - direct.energy_eV[-2]),
        ]
    )
    assert float(np.sum(direct.eedf * widths)) == pytest.approx(1.0, abs=1.0e-12)
    assert float(np.sum(direct.energy_eV * direct.eedf * widths)) == pytest.approx(
        direct.mean_energy_eV,
        rel=1.0e-12,
    )

    comparison = compare_eedf_cases(two_term, direct)
    metrics = comparison.metrics
    assert metrics["eedf_relative_l1"] < 0.01
    assert metrics["mean_energy_relative_difference"] < 0.005
    assert metrics["drift_velocity_relative_difference"] < 0.01
    assert metrics["major_rate_relative_difference"] < 0.02
    assert metrics["normalization_error_candidate"] < 1.0e-8
    assert comparison.failures == []


def test_direct_pn_uses_shared_kinetic_block_not_two_term_solver(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["multi_term"])
    data["solvers"]["multi_term"]["method"] = "pn_closure_direct"
    data["solvers"]["multi_term"]["lmax"] = 1
    cfg = load_config(write_config(tmp_path, data, name="shared_kinetic.yaml"))
    cross_sections = load_cross_sections(cfg.cross_sections, cfg.conditions)
    grid = make_two_term_energy_grid(
        cfg.internal.two_term,
        cross_sections=cross_sections,
    )
    shared = assemble_native_operator_blocks(
        cfg,
        cross_sections,
        cfg.run.e_over_n_Td[0],
        grid,
    )
    via_solver = TwoTermSolver(cfg, cross_sections).assemble_native_operator_block(
        cfg.run.e_over_n_Td[0],
        grid.energy_eV,
        grid.edges_eV,
        grid.widths_eV,
    )

    assert np.array_equal(shared.energy_eV, via_solver.energy_eV)
    assert np.array_equal(shared.widths_eV, via_solver.widths_eV)
    assert np.allclose(
        shared.collisions.nu_m,
        via_solver.collisions.nu_m,
        rtol=0.0,
        atol=0.0,
    )
    assert shared.collisions.sigma_total_like.shape == shared.energy_eV.shape
    assert shared.collisions.inelastic_loss_frequency_s_inv.shape == shared.energy_eV.shape
    assert np.all(shared.collisions.sigma_total_like > 0.0)
    assert np.all(shared.collisions.inelastic_loss_frequency_s_inv >= 0.0)
    assert (shared.matrix - via_solver.matrix).nnz == 0

    direct_source = (
        ROOT / "electron_swarm" / "solvers" / "multi_term" / "direct.py"
    ).read_text(encoding="utf-8")
    assert "TwoTermSolver" not in direct_source


def test_direct_pn_higher_l_damping_gate_uses_shared_kinetic_data(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["multi_term"])
    data["solvers"]["multi_term"]["method"] = "pn_closure_direct"
    data["solvers"]["multi_term"]["lmax"] = 2
    cfg = load_config(write_config(tmp_path, data, name="higher_l_damping.yaml"))
    cross_sections = load_cross_sections(cfg.cross_sections, cfg.conditions)
    grid = make_two_term_energy_grid(cfg.internal.two_term, cross_sections=cross_sections)
    shared = assemble_native_operator_blocks(
        cfg,
        cross_sections,
        cfg.run.e_over_n_Td[0],
        grid,
    )
    moments = np.zeros((3, len(shared.energy_eV)))
    moments[0] = 1.0
    damping = build_higher_l_collision_damping(shared, moments)

    assert damping.shape == moments.shape
    assert np.allclose(damping[1], shared.collisions.nu_m)
    assert np.all(np.isfinite(damping[2]))
    assert np.all(damping[2] >= 0.0)

    shared.collisions.sigma_total_like = np.zeros_like(shared.collisions.sigma_total_like)
    with pytest.raises(ValueError, match="sigma_total_like must be positive"):
        build_higher_l_collision_damping(shared, moments)


def test_shared_kinetic_eedf_helpers_and_rate_convolution(tmp_path: Path) -> None:
    energy = np.array([0.5, 1.5, 2.5])
    widths = np.array([1.0, 1.0, 1.0])
    eedf = normalize_eedf(np.array([1.0, 2.0, 1.0]), widths)

    assert float(np.sum(eedf * widths)) == pytest.approx(1.0)
    assert mean_energy_from_eedf(energy, widths, eedf) == pytest.approx(1.5)
    assert np.all(np.isfinite(eepf_from_eedf(energy, eedf)))
    assert negative_mass_fraction(np.array([-0.1, 1.1, 0.0]), widths) > 0.0

    cfg = load_config(write_config(tmp_path, base_product_config(tmp_path, ["two_term"])))
    case = run(cfg, write=False).cases[0]
    case_widths = cell_edges_from_centers(case.energy_eV)[1]
    cross_sections = load_cross_sections(cfg.cross_sections, cfg.conditions)
    rates = compute_rates_from_eedf(
        cfg,
        cross_sections,
        case.energy_eV,
        case_widths,
        case.eedf,
        case_id=case.case_id,
        e_over_n_Td=case.e_over_n_Td,
        solver_name=case.solver,
    )

    expected = {
        (rate.species, rate.process, rate.process_type): rate.mixture_weighted_rate_m3_s
        for rate in case.rates
    }
    actual = {
        (rate.species, rate.process, rate.process_type): rate.mixture_weighted_rate_m3_s
        for rate in rates.rates
    }
    assert actual.keys() == expected.keys()
    for key, value in expected.items():
        assert actual[key] == pytest.approx(value)


def test_internal_monte_carlo_ionization_energy_sharing_modes(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "backend": "internal",
        "angular_scattering": "same_as_physics",
    }

    data["physics"]["ionization"] = {"energy_sharing": "equal"}
    cfg = load_config(write_config(tmp_path, data, "mc_equal.yaml"))
    rng = np.random.default_rng(1)
    assert _post_reaction_energy(
        cfg, ProcessType.IONIZATION, 10.0, 30.0, rng
    ) == pytest.approx(10.0)
    outcome = _post_reaction_outcome(
        cfg, ProcessType.IONIZATION, 10.0, 30.0, np.random.default_rng(1)
    )
    assert outcome.tracked_energy_eV == pytest.approx(10.0)
    assert outcome.ionization_threshold_loss_eV == pytest.approx(10.0)
    assert outcome.ionization_untracked_secondary_energy_eV == pytest.approx(10.0)
    assert (
        outcome.tracked_energy_eV
        + outcome.ionization_threshold_loss_eV
        + outcome.ionization_untracked_secondary_energy_eV
    ) == pytest.approx(30.0)

    data["physics"]["ionization"] = {
        "energy_sharing": "primary_secondary",
        "secondary_electron_energy_eV": 2.0,
    }
    cfg = load_config(write_config(tmp_path, data, "mc_primary_secondary.yaml"))
    rng = np.random.default_rng(3)
    samples = {
        _post_reaction_energy(cfg, ProcessType.IONIZATION, 10.0, 30.0, rng)
        for _ in range(40)
    }
    assert samples == {2.0, 18.0}
    outcome = _post_reaction_outcome(
        cfg, ProcessType.IONIZATION, 10.0, 30.0, np.random.default_rng(3)
    )
    assert outcome.tracked_energy_eV in {2.0, 18.0}
    assert outcome.ionization_untracked_secondary_energy_eV in {2.0, 18.0}
    assert outcome.tracked_energy_eV != outcome.ionization_untracked_secondary_energy_eV
    assert (
        outcome.tracked_energy_eV
        + outcome.ionization_threshold_loss_eV
        + outcome.ionization_untracked_secondary_energy_eV
    ) == pytest.approx(30.0)

    data["physics"]["ionization"] = {"energy_sharing": "loss_only"}
    cfg = load_config(write_config(tmp_path, data, "mc_loss_only.yaml"))
    rng = np.random.default_rng(1)
    assert _post_reaction_energy(
        cfg, ProcessType.IONIZATION, 10.0, 30.0, rng
    ) == pytest.approx(20.0)
    outcome = _post_reaction_outcome(
        cfg, ProcessType.IONIZATION, 10.0, 30.0, np.random.default_rng(1)
    )
    assert outcome.tracked_energy_eV == pytest.approx(20.0)
    assert outcome.ionization_untracked_secondary_energy_eV == pytest.approx(0.0)


def test_internal_monte_carlo_energy_audit_and_bin_uncertainty(
    tmp_path: Path,
) -> None:
    audit = _EnergyAudit()
    audit.tracked_particle_initial_energy_eV = 10.0
    audit.record_field_push(10.0, 16.0)
    audit.record_elastic_collision(16.0, 15.5)
    audit.record_reaction(
        _post_reaction_outcome(
            load_config(
                write_config(
                    tmp_path,
                    base_product_config(tmp_path, ["monte_carlo"]),
                    "tmp_audit.yaml",
                )
            ),
            ProcessType.IONIZATION,
            5.0,
            15.5,
            np.random.default_rng(0),
        )
    )
    audit.tracked_particle_final_energy_eV = 5.25
    metadata = audit.as_metadata()
    assert metadata["mc_energy_balance_status"] == "ok"
    assert metadata["mc_tracked_energy_balance_residual_fraction"] < 1.0e-12
    assert metadata["mc_physical_branching_gap_eV"] == pytest.approx(5.25)

    counts = np.array([0, 4, 25])
    rel = _mc_bin_relative_standard_error(counts)
    assert np.isnan(rel[0])
    assert rel[1] == pytest.approx(0.5)
    assert rel[2] == pytest.approx(0.2)
    tail = _mc_tail_uncertainty_metadata(
        np.array([1.0, 20.0, 40.0]),
        counts,
        15.0,
        min_count=20,
    )
    assert tail["mc_tail_uncertainty_status"] == "insufficient"
    assert tail["mc_min_tail_bin_count"] == 4
    assert tail["mc_tail_effective_sample_count_min"] == pytest.approx(4.0)
    assert tail["mc_tail_weak_probability_fraction"] > 0.05
    assert tail["mc_max_resolved_energy_eV"] == pytest.approx(40.0)

    weighted = _mc_effective_bin_counts(
        np.array([0.0, 2.0, 4.0]),
        np.array([0.0, 2.0, 4.0]),
    )
    assert weighted[0] == pytest.approx(0.0)
    assert weighted[1] == pytest.approx(2.0)
    assert weighted[2] == pytest.approx(4.0)


def test_internal_monte_carlo_null_collision_majorant_and_bins() -> None:
    _validate_trial_collision_frequency(9.0, 10.0)
    with pytest.raises(RuntimeError, match="null-collision majorant"):
        _validate_trial_collision_frequency(10.1, 10.0)

    edges = _mc_energy_edges(250.0)
    assert edges[0] == pytest.approx(0.0)
    assert edges[-1] == pytest.approx(250.0)
    assert np.all(np.diff(edges) > 0.0)


@pytest.mark.parametrize("lmax", [2, 4])
def test_multi_term_pn_closure_direct_lmax_gt_one_smoke(
    tmp_path: Path,
    lmax: int,
) -> None:
    data = base_product_config(tmp_path, ["multi_term"])
    data["solvers"]["multi_term"]["method"] = "pn_closure_direct"
    data["solvers"]["multi_term"]["lmax"] = lmax
    cfg = load_config(write_config(tmp_path, data, name=f"direct_lmax{lmax}.yaml"))
    [case] = run(cfg, write=False).cases

    assert case.metadata["solver_method"] == "pn_closure_direct"
    assert case.metadata["direct_pn_operator"] is True
    assert case.metadata["lmax"] == lmax
    assert case.metadata["higher_l_collision_model"] == "angular_closure_damping"
    assert case.metadata["higher_l_inelastic_model"] == "sink_only"
    assert case.metadata["lmax_convergence_status"] in {"ok", "warning"}
    assert case.metadata["negative_mass_fraction"] < 1.0e-7
    assert np.all(np.isfinite(case.eedf))
    assert case.mean_energy_eV > 0.0


def test_multi_term_pn_closure_direct_lmax_gt_one_fails_for_b_field(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["multi_term"])
    data["solvers"]["multi_term"]["method"] = "pn_closure_direct"
    data["solvers"]["multi_term"]["lmax"] = 2
    data["physics"]["field"]["magnetic_field"] = {
        "enabled": True,
        "B_T": 0.01,
        "angle_EB_deg": 90.0,
    }
    cfg = load_config(write_config(tmp_path, data, name="direct_lmax2_b.yaml"))
    with pytest.raises(ValueError, match="magnetic field"):
        run(cfg, write=False)


def test_multi_term_pn_closure_direct_lmax_gt_one_fails_for_moment_table(
    tmp_path: Path,
) -> None:
    table = write_moment_table(tmp_path)
    data = base_product_config(tmp_path, ["multi_term"])
    data["solvers"]["multi_term"]["method"] = "pn_closure_direct"
    data["solvers"]["multi_term"]["lmax"] = 2
    data["physics"]["angular_scattering"] = {
        "model": "moment_table",
        "moment_table": {
            "path": table.as_posix(),
            "format": "normalized_legendre_moments",
            "provenance": "model_derived",
            "extrapolation": "error",
        },
    }
    cfg = load_config(write_config(tmp_path, data, name="direct_lmax2_table.yaml"))
    with pytest.raises(NotImplementedError, match="does not use moment_table"):
        run(cfg, write=False)


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
    assert case.metadata["direct_pn_operator"] is True
    assert case.metadata["exact_dcs_based"] is False
    assert case.metadata["ordinary_integral_xs_closure"] is True


def test_multi_term_pn_dcs_with_moment_table_runs(
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
            "provenance": "dcs_derived",
            "extrapolation": "error",
        },
    }
    cfg = load_config(write_config(tmp_path, data))
    [case] = run(cfg, write=False).cases

    assert case.metadata["solver_method"] == "pn_dcs"
    assert case.metadata["physics_level"] == "dcs_moment_based"
    assert case.metadata["angular_model"] == "moment_table"
    assert case.metadata["angular_moment_source"] == "moment_table"
    assert case.metadata["moment_table_provenance"] == "dcs_derived"
    assert case.metadata["exact_dcs_based"] is True
    assert case.metadata["ordinary_integral_xs_closure"] is False
    assert case.metadata["direct_pn_operator"] is True
    assert case.metadata["higher_l_collision_model"] == "moment_table_damping"
    assert np.all(np.isfinite(case.eedf))
    assert case.mean_energy_eV > 0.0


def test_multi_term_pn_dcs_lmax1_sanity(tmp_path: Path) -> None:
    table = write_moment_table(tmp_path, m1=0.0)
    data = base_product_config(tmp_path, ["multi_term"])
    data["solvers"]["multi_term"]["method"] = "pn_dcs"
    data["solvers"]["multi_term"]["lmax"] = 1
    data["physics"]["angular_scattering"] = {
        "model": "moment_table",
        "moment_table": {
            "path": table.as_posix(),
            "format": "normalized_legendre_moments",
            "provenance": "model_derived",
            "extrapolation": "error",
        },
    }
    [case] = run(load_config(write_config(tmp_path, data)), write=False).cases

    assert case.metadata["solver_method"] == "pn_dcs"
    assert case.metadata["physics_level"] == "table_moment_based"
    assert case.metadata["moment_table_provenance"] == "model_derived"
    assert case.metadata["exact_dcs_based"] is False
    assert case.metadata["ordinary_integral_xs_closure"] is False
    assert case.metadata["direct_pn_operator"] is True
    assert np.all(np.isfinite(case.eedf))
    assert case.metadata["negative_mass_fraction"] < 1.0e-8


def test_multi_term_pn_dcs_moment_table_affects_solution(tmp_path: Path) -> None:
    cases = []
    for name, m1 in [("isotropic_like", 0.0), ("forward", 0.8)]:
        table = write_moment_table(tmp_path, name=f"{name}.csv", m1=m1)
        data = base_product_config(tmp_path, ["multi_term"])
        data["solvers"]["multi_term"]["method"] = "pn_dcs"
        data["solvers"]["multi_term"]["lmax"] = 3
        data["physics"]["angular_scattering"] = {
            "model": "moment_table",
            "moment_table": {
                "path": table.as_posix(),
                "format": "normalized_legendre_moments",
                "provenance": "model_derived",
                "extrapolation": "error",
            },
        }
        cfg = load_config(write_config(tmp_path, data, name=f"{name}.yaml"))
        cases.append(run(cfg, write=False).cases[0])

    first, second = cases
    assert first.metadata["angular_moment_source"] == "moment_table"
    assert second.metadata["angular_moment_source"] == "moment_table"
    assert not np.isclose(first.mean_energy_eV, second.mean_energy_eV)
    assert not np.isclose(first.drift_velocity_m_s, second.drift_velocity_m_s)


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
        "warmup_collisions": 3,
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
    assert case.metadata["ionization_source_treatment"] == "equal"
    assert case.metadata["secondary_electron_tracking"] is False
    assert case.metadata["mc_population_model"] == "fixed_particle_single_daughter"
    assert case.metadata["mc_warmup_collisions"] == 3
    assert case.metadata["mc_production_collisions"] == 8
    assert case.metadata["ionization_branching_model"] == "single_daughter_sampling"
    assert case.metadata["inelastic_angular_model"] == "isotropic_reset"
    assert case.metadata["eedf_estimator"] == "time_sampled_null_clock"
    assert case.metadata["mc_samples"] > 0
    assert case.metadata["mc_energy_balance_status"] == "ok"
    assert case.metadata["mc_tail_uncertainty_status"] in {"ok", "insufficient"}
    assert case.metadata["mc_tail_comparison_status"] in {
        "ok",
        "weak_tail_statistics",
        "energy_balance_warning",
    }
    assert case.eedf_counts is not None
    assert int(np.sum(case.eedf_counts)) == case.metadata["mc_samples"]


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
