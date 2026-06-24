from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from electron_swarm import load_config, run
from electron_swarm.core.cross_sections import ProcessType, load_cross_sections
from electron_swarm.core.solver_configs import build_internal_solver_configs
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
    _ParticleEnsemble,
    _flux_transport_estimates,
)
from electron_swarm.solvers._internal_mc.collisions import (
    _ionization_daughters,
    _moment_cross_sections,
    _post_reaction_energy,
    _post_reaction_outcome,
    _project_processes,
    _validate_maxent_p1_cross_sections,
    _validate_trial_collision_frequency,
)
from electron_swarm.solvers._internal_mc.histogram import (
    _build_eedf_histogram,
    _mc_bin_relative_standard_error,
    _mc_effective_bin_counts,
    _mc_energy_edges,
    _mc_tail_uncertainty_metadata,
)
from electron_swarm.physics.angular_scattering import build_angular_model
from electron_swarm.solvers.two_term import TwoTermSolver

from product_helpers import ROOT, base_product_config, write_config, write_moment_table


def _write_total_momentum_xs(
    tmp_path: Path,
    *,
    include_total: bool = True,
    include_momentum: bool = True,
) -> Path:
    rows = [
        "species,process,type,threshold_eV,mass_amu,energy_eV,cross_section_m2",
    ]
    for energy in (0.0, 10.0, 100.0, 1000.0):
        if include_total:
            rows.append(f"Ar,total_elastic,total_elastic,0,39.948,{energy},4e-20")
        if include_momentum:
            rows.append(f"Ar,momentum,momentum,0,39.948,{energy},1e-20")
    path = tmp_path / "total_momentum_xs.csv"
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def test_two_term_smoke_run(tmp_path: Path) -> None:
    cfg = load_config(write_config(tmp_path, base_product_config(tmp_path, ["two_term"])))
    result = run(cfg, write=False)
    [case] = result.cases
    assert case.solver == "two_term"
    assert case.schema_version == "2"
    assert case.mean_energy_eV > 0.0
    tail = case.diagnostics["tail_metrics"]
    assert tail["energy_grid_tail_status"] in {"ok", "warning", "insufficient"}
    assert tail["tail_probability"] >= 0.0
    assert np.all(np.isfinite(case.eedf))


def test_multi_term_default_direct_smoke_and_minimum_metadata(tmp_path: Path) -> None:
    cfg = load_config(write_config(tmp_path, base_product_config(tmp_path, ["multi_term"])))
    result = run(cfg, write=False)
    [case] = result.cases
    assert case.solver == "multi_term"
    assert case.metadata["angular_model"] == "isotropic"
    assert case.metadata["angular_moment_source"] == "isotropic_closure"
    assert case.metadata["solver_method"] == "pn_closure_direct"
    assert case.metadata["lmax"] == 3
    assert case.metadata["exact_dcs_based"] is False
    assert case.metadata["ordinary_integral_xs_closure"] is True
    assert case.metadata["direct_pn_operator"] is True
    tail = case.diagnostics["tail_metrics"]
    assert tail["energy_grid_tail_status"] in {"ok", "warning", "insufficient"}
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
    assert direct.metadata["direct_pn_operator"] is True
    assert direct.metadata["exact_dcs_based"] is False
    assert direct.metadata["ordinary_integral_xs_closure"] is True
    direct_diag = direct.diagnostics["multi_term"]
    assert direct_diag["lmax1_regression_target"] == "two_term"
    assert direct_diag["negative_mass_fraction"] < 1.0e-8
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
    internal = build_internal_solver_configs(cfg.solvers, cfg.physics)
    grid = make_two_term_energy_grid(
        internal.two_term,
        cross_sections=cross_sections,
    )
    shared = assemble_native_operator_blocks(
        cfg,
        cross_sections,
        cfg.run.e_over_n_Td[0],
        grid,
        internal.two_term,
    )
    via_solver = TwoTermSolver(
        cfg, cross_sections, internal.two_term
    ).assemble_native_operator_block(
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
    internal = build_internal_solver_configs(cfg.solvers, cfg.physics)
    grid = make_two_term_energy_grid(internal.two_term, cross_sections=cross_sections)
    shared = assemble_native_operator_blocks(
        cfg,
        cross_sections,
        cfg.run.e_over_n_Td[0],
        grid,
        internal.two_term,
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
    data["solvers"]["monte_carlo"] = {}

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


def test_internal_monte_carlo_ionization_daughters_and_resampling(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {}

    data["physics"]["ionization"] = {"energy_sharing": "equal"}
    cfg = load_config(write_config(tmp_path, data, "equal_daughters.yaml"))
    daughters = _ionization_daughters(cfg, 10.0, 30.0)
    assert daughters.threshold_loss_eV == pytest.approx(10.0)
    assert daughters.energies_eV == pytest.approx((10.0, 10.0))

    data["physics"]["ionization"] = {
        "energy_sharing": "primary_secondary",
        "secondary_electron_energy_eV": 2.0,
    }
    cfg = load_config(write_config(tmp_path, data, "primary_secondary_daughters.yaml"))
    daughters = _ionization_daughters(cfg, 10.0, 30.0)
    assert daughters.energies_eV == pytest.approx((18.0, 2.0))
    assert sum(daughters.energies_eV) + daughters.threshold_loss_eV == pytest.approx(
        30.0
    )

    data["physics"]["ionization"] = {"energy_sharing": "loss_only"}
    cfg = load_config(write_config(tmp_path, data, "loss_only_daughters.yaml"))
    daughters = _ionization_daughters(cfg, 10.0, 30.0)
    assert daughters.energies_eV == pytest.approx((20.0,))

    ensemble = _ParticleEnsemble(
        positions=np.zeros((4, 3), dtype=float),
        velocities=np.ones((4, 3), dtype=float),
        times=np.arange(4, dtype=float),
        weights=np.array([1.0, 2.0, 3.0, 4.0], dtype=float),
    )
    assert ensemble.systematic_resample(3, np.random.default_rng(0)) is True
    assert len(ensemble) == 3
    assert float(np.sum(ensemble.weights)) == pytest.approx(10.0)
    assert np.allclose(ensemble.weights, 10.0 / 3.0)


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
    assert 80.0 in edges
    assert 90.0 in edges
    assert not np.any((edges > 80.0) & (edges < 90.0))

    energy, widths, eedf, counts, effective_counts = _build_eedf_histogram(
        np.array([0.0, 1.0, 3.0]),
        np.array([2, 1]),
        np.array([2.0, 2.0]),
        np.array([2.0, 4.0]),
    )
    assert energy == pytest.approx([0.5, 2.0])
    assert widths == pytest.approx([1.0, 2.0])
    assert float(np.sum(eedf * widths)) == pytest.approx(1.0)
    assert counts == pytest.approx([2, 1])
    assert effective_counts == pytest.approx([2.0, 1.0])


def test_internal_mc_maxent_p1_uses_total_and_momentum_xs(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["cross_sections"]["files"] = [
        {
            "path": _write_total_momentum_xs(tmp_path).as_posix(),
            "species": "Ar",
            "format": "csv",
        }
    ]
    data["physics"]["angular_scattering"] = {
        "model": "maxent_p1",
        "higher_moment_closure": "maxent",
    }
    cfg = load_config(write_config(tmp_path, data, "maxent_p1.yaml"))
    cross_sections = load_cross_sections(cfg.cross_sections, cfg.conditions)
    projected = _project_processes(cfg, cross_sections)
    _validate_maxent_p1_cross_sections(projected)
    total, momentum = _moment_cross_sections(projected, "Ar", 5.0)
    assert total == pytest.approx(4.0e-20)
    assert momentum == pytest.approx(1.0e-20)

    rng = np.random.default_rng(12)
    angular = build_angular_model(cfg)
    samples = np.array(
        [
            angular.sample_mu(
                5.0,
                rng,
                sigma_total=total,
                sigma_momentum=momentum,
            )
            for _ in range(20000)
        ],
        dtype=float,
    )
    assert float(np.mean(samples)) == pytest.approx(0.75, abs=0.025)


def test_internal_mc_maxent_p1_requires_total_and_momentum_xs(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "particles": 4,
        "max_collisions": 2,
    }
    data["physics"]["angular_scattering"] = {
        "model": "maxent_p1",
        "higher_moment_closure": "maxent",
    }
    data["cross_sections"]["files"] = [
        {
            "path": _write_total_momentum_xs(
                tmp_path,
                include_total=False,
                include_momentum=True,
            ).as_posix(),
            "species": "Ar",
            "format": "csv",
        }
    ]
    with pytest.raises(NotImplementedError, match="requires both"):
        run(load_config(write_config(tmp_path, data, "maxent_missing_total.yaml")), write=False)


def test_internal_mc_flux_transport_uses_total_time_estimator() -> None:
    ensemble = _ParticleEnsemble(
        positions=np.array([[0.0, 0.0, -2.0], [1.0, 2.0, -8.0]], dtype=float),
        velocities=np.zeros((2, 3), dtype=float),
        times=np.array([1.0, 3.0], dtype=float),
        weights=np.array([1.0, 2.0], dtype=float),
    )

    drift, mobility, diffusion_l, diffusion_t, mean_time = (
        _flux_transport_estimates(ensemble, 10.0)
    )
    total_time = 1.0 * 1.0 + 2.0 * 3.0
    expected_drift = -((1.0 * -2.0) + (2.0 * -8.0)) / total_time
    expected_dz = np.array(
        [
            2.0 - expected_drift * 1.0,
            8.0 - expected_drift * 3.0,
        ]
    )
    expected_diff_l = float(
        np.sum(np.array([1.0, 2.0]) * expected_dz**2) / (2.0 * total_time)
    )
    expected_diff_t = float((1.0 * 0.0 + 2.0 * (1.0**2 + 2.0**2)) / (4.0 * total_time))
    assert drift == pytest.approx(expected_drift)
    assert mobility == pytest.approx(expected_drift / 10.0)
    assert diffusion_l == pytest.approx(expected_diff_l)
    assert diffusion_t == pytest.approx(expected_diff_t)
    assert mean_time == pytest.approx(total_time / 3.0)


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
    assert case.metadata["transport_definition"] == "f0_gradient_reconstruction"
    assert case.diagnostics["multi_term"]["negative_mass_fraction"] < 1.0e-7
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
    assert case.metadata["angular_model"] == "moment_table"
    assert case.metadata["angular_moment_source"] == "moment_table"
    assert case.metadata["moment_table_provenance"] == "dcs_derived"
    assert case.metadata["exact_dcs_based"] is True
    assert case.metadata["ordinary_integral_xs_closure"] is False
    assert case.metadata["direct_pn_operator"] is True
    assert "higher_l_collision_model" not in case.metadata
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
    assert case.metadata["moment_table_provenance"] == "model_derived"
    assert case.metadata["exact_dcs_based"] is False
    assert case.metadata["ordinary_integral_xs_closure"] is False
    assert case.metadata["direct_pn_operator"] is True
    assert np.all(np.isfinite(case.eedf))
    assert case.diagnostics["multi_term"]["negative_mass_fraction"] < 1.0e-8


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
    assert case.metadata["electron_electron_transport_stale"] is True
    assert case.drift_velocity_m_s == pytest.approx(baseline.drift_velocity_m_s)
    assert case.diffusion_L_m2_s == pytest.approx(baseline.diffusion_L_m2_s)


def test_two_term_fp_energy_postprocess_marks_transport_stale(
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
    assert case.metadata["electron_electron_transport_stale"] is True
    assert case.drift_velocity_m_s == pytest.approx(baseline.drift_velocity_m_s)
    assert case.diffusion_L_m2_s == pytest.approx(baseline.diffusion_L_m2_s)


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
    assert case.metadata["electron_electron_transport_stale"] is True
    assert case.drift_velocity_m_s == pytest.approx(baseline.drift_velocity_m_s)
    assert case.diffusion_L_m2_s == pytest.approx(baseline.diffusion_L_m2_s)


def test_monte_carlo_same_as_physics_validates_reported_angular_metadata(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "particles": 8,
        "max_collisions": 4,
        "seed": 11,
    }
    cfg = load_config(write_config(tmp_path, data))
    result = run(cfg, write=False)
    [case] = result.cases

    assert case.solver == "monte_carlo"
    assert case.metadata["angular_model"] == "isotropic"
    assert case.metadata["angular_moment_source"] == "isotropic_closure"
    assert case.metadata["ordinary_integral_xs_closure"] is True


@pytest.mark.mc
def test_internal_monte_carlo_magnetic_smoke_run(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {
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
    assert case.metadata["ionization_source_treatment"] == "equal"
    assert "mc_run" not in case.metadata
    assert "mc_energy_audit" not in case.metadata
    assert "mc_tail_audit" not in case.metadata
    assert "_dev_internal_monte_carlo_audit" not in case.metadata
    assert "internal_monte_carlo_audit" not in case.diagnostics
    assert case.eedf_counts is not None
    assert int(np.sum(case.eedf_counts)) > 0
    assert case.energy_widths_eV is not None
    assert float(np.sum(case.eedf * case.energy_widths_eV)) == pytest.approx(1.0)


@pytest.mark.mc
def test_internal_monte_carlo_weighted_branching_metadata(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["run"]["e_over_n_Td"] = [400.0]
    data["cross_sections"]["high_energy_extrapolation"] = "hold"
    data["solvers"]["monte_carlo"] = {
        "population_model": "weighted_branching",
        "particles": 12,
        "warmup_collisions": 2,
        "max_collisions": 12,
        "seed": 31,
    }
    cfg = load_config(write_config(tmp_path, data))
    [case] = run(cfg, write=False, collect_diagnostics=True).cases

    assert "_dev_internal_monte_carlo_audit" not in case.metadata
    dev_audit = case.diagnostics["internal_monte_carlo_audit"]
    run_meta = dev_audit["mc_run"]
    energy_meta = dev_audit["mc_energy_audit"]
    assert run_meta["mc_population_model"] == "weighted_branching"
    assert run_meta["secondary_electron_tracking"] is True
    assert (
        run_meta["swarm_population_treatment"]
        == "weighted_branching_resampled"
    )
    assert (
        run_meta["nonconservative_growth_treatment"]
        == "explicit_weighted_branching"
    )
    assert (
        run_meta["ionization_branching_model"]
        == "two_daughter_weighted_resampling"
    )
    assert (
        case.metadata["transport_definition"]
        == "mc_flux_particle_tracking_weighted_growth_population"
    )
    assert run_meta["transport_has_bulk"] is False
    for key in [
        "mc_secondary_electron_count",
        "mc_branching_resample_count",
        "mc_population_total_weight_final",
        "mc_population_log_growth_estimate_s_inv",
        "mc_population_weight_cv",
    ]:
        assert np.isfinite(run_meta[key])
    assert np.isfinite(energy_meta["mc_population_resampling_energy_adjustment_eV"])
    assert run_meta["mc_secondary_electron_count"] >= 0
    assert run_meta["mc_branching_resample_count"] >= 0
    assert run_meta["mc_population_total_weight_final"] > 0.0
    assert energy_meta["mc_physical_branching_gap_eV"] == pytest.approx(0.0)
    assert case.energy_widths_eV is not None
    assert float(np.sum(case.eedf * case.energy_widths_eV)) == pytest.approx(1.0)


@pytest.mark.mc
def test_internal_monte_carlo_zero_b_matches_disabled_with_same_seed(
    tmp_path: Path,
) -> None:
    base = base_product_config(tmp_path, ["monte_carlo"])
    base["solvers"]["monte_carlo"] = {
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
    assert disabled.metadata["magnetic_field_treatment"] == "none"
    assert zero_b.metadata["magnetic_field_treatment"] == "boris_lorentz_push"


def test_external_monte_carlo_product_fields_are_rejected(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "python_api": "product_helpers:fake_mc_missing",
    }
    with pytest.raises(ValueError, match="internal product backend only"):
        load_config(write_config(tmp_path, data))
