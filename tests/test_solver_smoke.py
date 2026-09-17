from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import electron_swarm.solvers.two_term.transport as kinetic_module
from electron_swarm import load_config, run
from electron_swarm.core.numerics import eedf_from_f0, eepf_from_eedf, f0_from_eedf
from electron_swarm.core.cross_sections import (
    CrossSectionProcess,
    ProcessType,
    load_cross_sections,
    prepare_active_mixture_inputs,
)
from electron_swarm.core.results import RateResult
from electron_swarm.core.result_metadata import PRODUCT_CASE_METADATA_KEYS
from electron_swarm.core.solver_configs import (
    ConvergenceConfig,
    build_internal_solver_configs,
)
from electron_swarm.solvers.boltzmann_common.grid import cell_edges_from_centers
from electron_swarm.solvers.boltzmann_common.observables import (
    compute_rates_from_eedf,
    mean_energy_from_eedf,
    negative_mass_fraction,
    normalize_eedf,
    weighted_integral,
)
from electron_swarm.solvers.boltzmann_common.operators import (
    assemble_energy_flux_operator,
)
from electron_swarm.solvers.multi_term.case import build_multiterm_case
from electron_swarm.solvers.multi_term.operator import assemble_pn_operator
from electron_swarm.solvers.multi_term.steady import solve_stationary_pn
from electron_swarm.solvers.monte_carlo.population import _ParticleEnsemble
from electron_swarm.solvers.monte_carlo.collisions import (
    _ionization_daughters,
    _moment_cross_sections,
    _post_reaction_outcome,
    _project_processes,
)
from electron_swarm.solvers.monte_carlo.histogram import (
    _build_eedf_histogram,
    _mc_energy_edges,
)
from electron_swarm.solvers.monte_carlo.evidence import (
    DIRECT_MC_TRANSPORT_LAG_PLANES,
    MC_EEDF_ESTIMATOR_SCHEMA_VERSION,
    MC_MAGNETIC_EEDF_ESTIMATOR_SCHEMA_VERSION,
    WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
    WEIGHTED_MC_TRANSPORT_MIN_EFFECTIVE_LINEAGE_COUNT,
    WEIGHTED_MC_TRANSPORT_MIN_EFFECTIVE_LINEAGE_FRACTION,
)
from electron_swarm.solvers.monte_carlo.direct_transport import (
    SynchronizedFluxTransportObserver,
    direct_flux_transport_snapshot,
    direct_transport_observation_times,
)
from electron_swarm.solvers.monte_carlo.weighted_transport import (
    SynchronizedWeightedGrowthFluxObserver,
    weighted_growth_flux_transport_snapshot,
    weighted_growth_transport_lag_plan,
)
from electron_swarm.physics.angular_scattering import build_angular_model
from electron_swarm.solvers.two_term.steady import assemble_native_operator_block

from product_helpers import ROOT, base_product_config, write_config, write_moment_table


_WEIGHTED_GROWTH_LAG_PLAN = weighted_growth_transport_lag_plan()


def _assert_canonical_transport(case) -> None:
    transport = case.transport
    density = transport.gas_number_density_m3
    assert density > 0.0
    assert case.mobility_m2_V_s == pytest.approx(
        case.reduced_mobility_m2_V_s_m3 / density
    )
    assert case.diffusion_L_m2_s == pytest.approx(
        case.reduced_diffusion_L_m2_s_m3 / density
    )
    assert case.diffusion_T_m2_s == pytest.approx(
        case.reduced_diffusion_T_m2_s_m3 / density
    )
    if case.solver == "two_term":
        assert case.reduced_electron_energy_mobility_m2_V_s_m3 is not None
        assert case.reduced_electron_energy_diffusion_m2_s_m3 is not None
        assert case.electron_energy_mobility_m2_V_s == pytest.approx(
            case.reduced_electron_energy_mobility_m2_V_s_m3 / density
        )
        assert case.electron_energy_diffusion_m2_s == pytest.approx(
            case.reduced_electron_energy_diffusion_m2_s_m3 / density
        )
    elif case.solver == "monte_carlo":
        assert transport.electron_energy_mobility_m2_V_s is not None
        assert transport.electron_energy_diffusion_m2_s is None
        assert transport.electron_energy_diffusion_L_m2_s is not None
        assert transport.electron_energy_diffusion_T_m2_s is not None


def _assert_eedf_normalized(case) -> None:
    widths = case.energy_widths_eV
    if widths is None:
        widths = cell_edges_from_centers(case.energy_eV)[1]
    assert float(np.sum(case.eedf * widths)) == pytest.approx(1.0)


def _relative_difference(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-300)


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
    cfg = load_config(
        write_config(tmp_path, base_product_config(tmp_path, ["two_term"]))
    )
    result = run(cfg, write=False)
    [case] = result.cases
    assert case.solver == "two_term"
    assert case.schema_version == "2"
    assert case.mean_energy_eV > 0.0
    tail = case.diagnostics["tail_metrics"]
    assert tail["energy_grid_tail_status"] in {"ok", "warning", "insufficient"}
    assert tail["tail_probability"] >= 0.0
    assert np.all(np.isfinite(case.eedf))
    _assert_eedf_normalized(case)
    _assert_canonical_transport(case)


def test_two_term_low_field_accepts_only_componentwise_roundoff_floor(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["two_term"])
    data["run"]["e_over_n_Td"] = [10.0]
    cfg = load_config(write_config(tmp_path, data, "two_term_10td.yaml"))

    [case] = run(cfg, write=False).cases
    diagnostic = case.diagnostics["two_term"]

    assert diagnostic["converged"] is True
    assert diagnostic["residual_roundoff_limited"] is True
    assert diagnostic["residual_L1"] > diagnostic["residual_requested_tolerance"]
    assert diagnostic["residual_L1"] <= diagnostic["residual_tolerance"]
    assert diagnostic["residual_backward_error"] <= 2.0e-15


def test_two_term_elastic_loss_is_same_discrete_operator_moment(
    tmp_path: Path,
) -> None:
    cfg = load_config(
        write_config(tmp_path, base_product_config(tmp_path, ["two_term"]))
    )
    [case] = run(cfg, write=False).cases
    cross_sections = prepare_active_mixture_inputs(
        load_cross_sections(cfg.cross_sections, cfg.conditions),
        cfg.conditions,
    )
    internal = build_internal_solver_configs(cfg.solvers, cfg.physics)
    edges, widths = cell_edges_from_centers(case.energy_eV)
    block = assemble_native_operator_block(
        cfg,
        cross_sections,
        internal.two_term,
        case.e_over_n_Td,
        case.energy_eV,
        edges,
        widths,
    )
    elastic_operator = assemble_energy_flux_operator(
        case.energy_eV,
        widths,
        0.0,
        block.collisions,
    )
    expected = (
        -weighted_integral(
            case.energy_eV * (elastic_operator @ case.eedf),
            widths,
        )
        / block.gas_number_density_m3
    )

    artifact = case.diagnostics["two_term"]["elastic_energy_loss"]
    assert artifact["rate_coefficient_eV_m3_s"] == pytest.approx(expected, rel=1.0e-13)
    assert artifact["source_eedf"] == "same_solved_eedf"
    assert artifact["gas_temperature_terms_included"] is True
    assert artifact["neutral_thermal_motion_model"] == (
        "finite_temperature_fokker_planck"
    )


def test_multi_term_default_direct_smoke_and_minimum_metadata(tmp_path: Path) -> None:
    cfg = load_config(
        write_config(tmp_path, base_product_config(tmp_path, ["multi_term"]))
    )
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
    _assert_eedf_normalized(case)
    _assert_canonical_transport(case)


def test_multi_term_pn_closure_direct_lmax1_is_physical_pn_sanity(
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
    assert direct_diag["negative_mass_fraction"] < 1.0e-8
    assert direct_diag["pn_converged"] is True
    assert direct_diag["pn_full_relative_residual"] < 1.0e-6
    assert direct_diag["field_coupling"] == ("bidirectional_adjacent_legendre_moments")
    assert direct_diag["drift_observable"] == "F1_velocity_moment"
    assert set(direct.metadata) <= PRODUCT_CASE_METADATA_KEYS
    assert np.all(np.isfinite(direct.eedf))
    widths = direct.energy_widths_eV
    assert widths is not None
    assert float(np.sum(direct.eedf * widths)) == pytest.approx(1.0, abs=1.0e-12)
    assert float(np.sum(direct.energy_eV * direct.eedf * widths)) == pytest.approx(
        direct.mean_energy_eV,
        rel=1.0e-12,
    )

    assert len(two_term.energy_eV) != len(direct.energy_eV)
    assert _relative_difference(two_term.mean_energy_eV, direct.mean_energy_eV) < 0.10
    assert (
        _relative_difference(
            two_term.drift_velocity_m_s,
            direct.drift_velocity_m_s,
        )
        < 0.10
    )
    two_rates = {rate.process: rate.rate_coefficient_m3_s for rate in two_term.rates}
    direct_rates = {rate.process: rate.rate_coefficient_m3_s for rate in direct.rates}
    common_rates = set(two_rates) & set(direct_rates)
    assert common_rates
    assert (
        max(
            _relative_difference(two_rates[process], direct_rates[process])
            for process in common_rates
        )
        < 0.15
    )


def test_direct_pn_uses_configured_grid_and_bidirectional_blocks(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["multi_term"])
    data["solvers"]["multi_term"]["method"] = "pn_closure_direct"
    data["solvers"]["multi_term"]["lmax"] = 2
    cfg = load_config(write_config(tmp_path, data, name="shared_kinetic.yaml"))
    cross_sections = load_cross_sections(cfg.cross_sections, cfg.conditions)
    internal = build_internal_solver_configs(cfg.solvers, cfg.physics)
    case = build_multiterm_case(
        cfg,
        cross_sections,
        internal.multi_term,
        cfg.run.e_over_n_Td[0],
    )
    operator = assemble_pn_operator(case)

    assert np.array_equal(operator.energy_eV, case.grid.centers_eV)
    assert np.array_equal(operator.widths_eV, case.grid.widths_eV)
    s0, s1, s2 = operator.moment_slices[:3]
    assert operator.matrix[s0, s1].nnz > 0
    assert operator.matrix[s1, s0].nnz > 0
    assert operator.matrix[s1, s2].nnz > 0
    assert operator.matrix[s2, s1].nnz > 0
    field_block = operator.matrix[s0, s1].toarray()
    conservation_defect = operator.widths_eV @ field_block
    column_scale = np.max(
        np.sum(np.abs(operator.widths_eV[:, None] * field_block), axis=0)
    )
    assert np.max(np.abs(conservation_defect)) / column_scale < 1.0e-14

    operator_source = (
        ROOT / "electron_swarm" / "solvers" / "multi_term" / "operator.py"
    ).read_text(encoding="utf-8")
    assert "TwoTermSolver" not in operator_source
    assert "assemble_native_operator_blocks" not in operator_source


def test_direct_pn_angular_damping_is_finite_on_solved_grid(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["multi_term"])
    data["solvers"]["multi_term"]["method"] = "pn_closure_direct"
    data["solvers"]["multi_term"]["lmax"] = 2
    cfg = load_config(write_config(tmp_path, data, name="higher_l_damping.yaml"))
    cross_sections = load_cross_sections(cfg.cross_sections, cfg.conditions)
    internal = build_internal_solver_configs(cfg.solvers, cfg.physics)
    case = build_multiterm_case(
        cfg,
        cross_sections,
        internal.multi_term,
        cfg.run.e_over_n_Td[0],
    )
    operator = assemble_pn_operator(case)

    assert operator.angular_damping_s_inv.shape == operator.angular_moments.shape
    assert np.all(np.isfinite(operator.angular_damping_s_inv[1:]))
    assert np.all(operator.angular_damping_s_inv[1:] > 0.0)


def test_shared_kinetic_eedf_helpers_and_rate_convolution(tmp_path: Path) -> None:
    energy = np.array([0.5, 1.5, 2.5])
    widths = np.array([1.0, 1.0, 1.0])
    eedf = normalize_eedf(np.array([1.0, 2.0, 1.0]), widths)

    assert float(np.sum(eedf * widths)) == pytest.approx(1.0)
    assert mean_energy_from_eedf(energy, widths, eedf) == pytest.approx(1.5)
    assert np.all(np.isfinite(eepf_from_eedf(energy, eedf)))
    f0 = f0_from_eedf(energy, eedf)
    assert float(np.sum(np.sqrt(energy) * f0 * widths)) == pytest.approx(1.0)
    assert np.allclose(eedf_from_f0(energy, f0), eedf)
    assert negative_mass_fraction(np.array([-0.1, 1.1, 0.0]), widths) > 0.0

    cfg = load_config(
        write_config(tmp_path, base_product_config(tmp_path, ["two_term"]))
    )
    case = run(cfg, write=False).cases[0]
    case_widths = cell_edges_from_centers(case.energy_eV)[1]
    cross_sections = load_cross_sections(cfg.cross_sections, cfg.conditions)
    internal = build_internal_solver_configs(cfg.solvers, cfg.physics)
    block = assemble_native_operator_block(
        cfg,
        cross_sections,
        internal.two_term,
        case.e_over_n_Td,
        case.energy_eV,
        *cell_edges_from_centers(case.energy_eV),
    )
    expected_transport = kinetic_module.transport_from_eedf(
        cfg,
        block,
        case.eedf,
        case.e_over_n_Td,
        internal.two_term,
        temporal_growth_frequency_s_inv=case.diagnostics["two_term"][
            "growth_frequency_s-1"
        ],
    )
    assert case.reduced_mobility_m2_V_s_m3 == pytest.approx(
        expected_transport.reduced_mobility_m2_V_s_m3
    )
    assert case.reduced_diffusion_L_m2_s_m3 == pytest.approx(
        expected_transport.reduced_diffusion_L_m2_s_m3
    )
    assert case.reduced_electron_energy_mobility_m2_V_s_m3 == pytest.approx(
        expected_transport.reduced_electron_energy_mobility_m2_V_s_m3
    )
    assert case.reduced_electron_energy_diffusion_m2_s_m3 == pytest.approx(
        expected_transport.reduced_electron_energy_diffusion_m2_s_m3
    )
    assert all(
        np.isfinite(v)
        for v in (
            expected_transport.reduced_mobility_m2_V_s_m3,
            expected_transport.reduced_diffusion_L_m2_s_m3,
            expected_transport.reduced_electron_energy_mobility_m2_V_s_m3,
            expected_transport.reduced_electron_energy_diffusion_m2_s_m3,
        )
    )
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


def test_two_term_energy_transport_density_scaling(tmp_path: Path) -> None:
    cases = []
    for index, density in enumerate((1.0e21, 2.0e21)):
        data = base_product_config(tmp_path, ["two_term"])
        data["conditions"].pop("pressure_Pa")
        data["conditions"]["gas_number_density_m3"] = density
        cfg = load_config(write_config(tmp_path, data, name=f"density_{index}.yaml"))
        cases.append(run(cfg, write=False).cases[0])

    low, high = cases
    for attr in (
        "reduced_mobility_m2_V_s_m3",
        "reduced_diffusion_L_m2_s_m3",
        "reduced_electron_energy_mobility_m2_V_s_m3",
        "reduced_electron_energy_diffusion_m2_s_m3",
    ):
        assert getattr(high, attr) == pytest.approx(getattr(low, attr), rel=1.0e-6)
    assert high.mobility_m2_V_s == pytest.approx(low.mobility_m2_V_s / 2.0, rel=1.0e-6)
    assert high.diffusion_L_m2_s == pytest.approx(
        low.diffusion_L_m2_s / 2.0, rel=1.0e-6
    )
    assert high.electron_energy_mobility_m2_V_s == pytest.approx(
        low.electron_energy_mobility_m2_V_s / 2.0, rel=1.0e-6
    )
    assert high.electron_energy_diffusion_m2_s == pytest.approx(
        low.electron_energy_diffusion_m2_s / 2.0, rel=1.0e-6
    )


def test_two_term_transport_refuses_nonpositive_kinetic_moment(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = load_config(
        write_config(tmp_path, base_product_config(tmp_path, ["two_term"]))
    )
    case = run(cfg, write=False).cases[0]
    cross_sections = load_cross_sections(cfg.cross_sections, cfg.conditions)
    internal = build_internal_solver_configs(cfg.solvers, cfg.physics)
    block = assemble_native_operator_block(
        cfg,
        cross_sections,
        internal.two_term,
        case.e_over_n_Td,
        case.energy_eV,
        *cell_edges_from_centers(case.energy_eV),
    )

    monkeypatch.setattr(
        kinetic_module,
        "f0_reduced_transport_from_eedf",
        lambda *_args, **_kwargs: (-1.0, 2.0, 3.0, 4.0),
    )

    with pytest.raises(
        FloatingPointError,
        match="refusing to replace kinetic output with a Drude fallback",
    ):
        kinetic_module.transport_from_eedf(
            cfg,
            block,
            case.eedf,
            case.e_over_n_Td,
            internal.two_term,
        )


def test_rate_result_energy_loss_signs_and_derived_values() -> None:
    base = {
        "solver": "two_term",
        "case_id": "case",
        "e_over_n_Td": 10.0,
        "species": "Ar",
        "rate_coefficient_m3_s": 2.0,
        "target_species_fraction": 0.25,
        "gas_number_density_m3": 3.0,
    }
    expected_losses = {
        "excitation": 4.0,
        "ionization": 4.0,
        "attachment": 0.0,
        "superelastic": -4.0,
    }
    for process_type, loss in expected_losses.items():
        rate = RateResult(
            **base,
            process=process_type,
            process_type=process_type,
            threshold_eV=4.0,
        )
        assert rate.energy_loss_eV == pytest.approx(loss)
        assert rate.mixture_weighted_rate_m3_s == pytest.approx(0.5)
        assert rate.frequency_s_inv == pytest.approx(1.5)
        assert rate.energy_loss_rate_coefficient_eV_m3_s == pytest.approx(2.0 * loss)
        assert rate.power_loss_eV_s == pytest.approx(1.5 * loss)


def test_internal_monte_carlo_ionization_energy_sharing_modes(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {}

    data["physics"]["ionization"] = {"energy_sharing": "equal"}
    cfg = load_config(write_config(tmp_path, data, "mc_equal.yaml"))
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
        _post_reaction_outcome(
            cfg, ProcessType.IONIZATION, 10.0, 30.0, rng
        ).tracked_energy_eV
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


def test_internal_monte_carlo_energy_bins() -> None:
    edges = _mc_energy_edges(250.0, thresholds_eV=(11.55, 15.76))
    assert edges[0] == pytest.approx(0.0)
    assert edges[-1] == pytest.approx(250.0)
    assert np.all(np.diff(edges) > 0.0)
    assert 11.55 in edges
    assert 15.76 in edges
    low_energy_edges = edges[edges <= 20.0]
    threshold_edges = np.isin(low_energy_edges, [11.55, 15.76])
    regular_edges = low_energy_edges[~threshold_edges]
    # The core is uniform in speed (sqrt energy), giving progressively finer
    # energy cells toward the origin while retaining a bounded width at 20 eV.
    assert np.diff(np.sqrt(regular_edges)) == pytest.approx(
        np.diff(np.sqrt(regular_edges))[0]
    )
    assert float(np.max(np.diff(low_energy_edges))) <= 0.025
    regular_tail_edges = edges[~np.isin(edges, [11.55, 15.76])]
    tail_widths = np.diff(regular_tail_edges[regular_tail_edges >= 20.0])
    assert len(tail_widths) > 10
    growth_widths = tail_widths[:-1]
    assert np.all(np.diff(growth_widths) >= -1.0e-12)
    assert np.max(growth_widths[1:] / growth_widths[:-1]) <= 1.04 + 1.0e-12

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


def test_internal_monte_carlo_energy_bins_coalesce_one_ulp_thresholds() -> None:
    threshold = 4.5
    near_regular_edge = np.nextafter(threshold, np.inf)
    edges = _mc_energy_edges(
        20.0,
        thresholds_eV=(threshold, near_regular_edge),
    )

    assert np.count_nonzero(np.isclose(edges, threshold, rtol=0.0, atol=1.0e-14)) == 1
    assert threshold in edges
    widths = np.diff(edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    reconstructed_left = centers - 0.5 * widths
    reconstructed_right = centers + 0.5 * widths
    assert np.allclose(reconstructed_right[:-1], reconstructed_left[1:], atol=1.0e-14)


def test_internal_monte_carlo_energy_bins_retain_rate_kernel_landmarks() -> None:
    process = CrossSectionProcess(
        species="X",
        process="narrow excitation",
        process_type=ProcessType.EXCITATION,
        threshold_eV=0.18,
        energy_eV=np.asarray([0.0, 0.18, 0.19, 0.20, 0.22, 0.40]),
        cross_section_m2=np.asarray([0.0, 0.0, 1.0e-20, 2.0e-20, 0.0, 0.0]),
    )

    edges = _mc_energy_edges(
        20.0,
        thresholds_eV=(0.18,),
        processes=(process,),
    )

    assert {0.18, 0.19, 0.20, 0.22} <= set(edges)
    local = edges[(edges >= 0.15) & (edges <= 0.25)]
    assert float(np.max(np.diff(local))) < 0.003


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
    cross_sections = prepare_active_mixture_inputs(
        load_cross_sections(cfg.cross_sections, cfg.conditions),
        cfg.conditions,
    )
    projected = _project_processes(cfg, cross_sections)
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
    with pytest.raises(ValueError, match="requires both"):
        run(
            load_config(write_config(tmp_path, data, "maxent_missing_total.yaml")),
            write=False,
        )


def test_direct_mc_block_helfand_recovers_anisotropic_diffusion() -> None:
    age = 2.0
    diffusion_l = 5.0
    diffusion_t = 3.0
    drift = 4.0
    transverse_scale = np.sqrt(6.0 * diffusion_t * age)
    longitudinal_scale = np.sqrt(6.0 * diffusion_l * age)
    fluctuations = np.array(
        [
            [transverse_scale, 0.0, 0.0],
            [-transverse_scale, 0.0, 0.0],
            [0.0, transverse_scale, 0.0],
            [0.0, -transverse_scale, 0.0],
            [0.0, 0.0, longitudinal_scale],
            [0.0, 0.0, -longitudinal_scale],
        ],
    )
    positions = fluctuations + np.array([0.0, 0.0, -drift * age])
    velocities = 0.5 * fluctuations / age + np.array([0.0, 0.0, -drift])
    estimate = direct_flux_transport_snapshot(
        positions_m=positions,
        velocities_m_s=velocities,
        ages_s=np.full(6, age),
        weights=np.ones(6),
        electric_field_V_m=np.array([0.0, 0.0, 2.0]),
    )

    assert estimate.drift_velocity_m_s == pytest.approx(drift)
    assert estimate.mobility_m2_V_s == pytest.approx(drift / 2.0)
    assert estimate.diffusion_L_m2_s == pytest.approx(diffusion_l)
    assert estimate.diffusion_T_m2_s == pytest.approx(diffusion_t)


def _deterministic_transport_observer(production_steps: int):
    times = direct_transport_observation_times(
        production_trial_steps=production_steps,
        trial_frequency_s_inv=1.0,
    )
    observer = SynchronizedFluxTransportObserver(
        observation_times_s=times,
        particles=2,
        electric_field_V_m=np.array([0.0, 0.0, 1.0]),
    )
    previous = 0.0
    for time_s in times:
        dt_s = float(time_s - previous)
        for particle, terminal_noise in enumerate((-1000.0, 1000.0)):
            observer.record_residence(
                displacement_m=np.array([0.0, 0.0, -2.0 * dt_s]),
                sample_energy_eV=3.0,
                dt_s=dt_s,
                weight=1.0,
            )
            assert observer.record_due(
                particle_index=particle,
                time_s=float(time_s),
                position_m=np.array([float(particle), 0.0, -2.0 * time_s]),
                velocity_m_s=np.array([terminal_noise, 0.0, -2.0]),
                weight=1.0,
            )
        previous = float(time_s)
    return observer.finalize()


def test_direct_mc_fixed_lag_does_not_grow_with_production() -> None:
    short, short_diagnostics = _deterministic_transport_observer(2560)
    long, long_diagnostics = _deterministic_transport_observer(5120)

    short_scan = short_diagnostics["lag_scan"]
    long_scan = long_diagnostics["lag_scan"]
    assert short_scan["lag_planes"] == list(DIRECT_MC_TRANSPORT_LAG_PLANES)
    assert long_scan["lag_s"] == pytest.approx(short_scan["lag_s"])
    assert long_scan["completed_origins"] == [
        value + 128 for value in short_scan["completed_origins"]
    ]
    assert short.mobility_m2_V_s == pytest.approx(2.0)
    assert long.mobility_m2_V_s == pytest.approx(2.0)
    assert short.energy_mobility_m2_V_s == pytest.approx(2.0)
    assert long.energy_mobility_m2_V_s == pytest.approx(2.0)
    assert short_diagnostics["production_estimates"]["drift_velocity_m_s"] == (
        pytest.approx(2.0)
    )
    assert (
        short_diagnostics["transport_sampling"]["production_length_controls_lag"]
        is False
    )
    assert short_diagnostics["estimator_schema_version"] == "direct_mc_transport.v4"
    assert short_diagnostics["component_estimators"] == {
        "drift_mobility": "production_trajectory_displacement_over_residence_time",
        "particle_diffusion": "block_helfand_mean_square_displacement",
        "energy_mobility": "production_residence_energy_current",
        "energy_diffusion": "restricted_density_packet_energy_flux_fixed_lag",
    }


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
    assert case.metadata["transport_definition"] == (
        "pn_f1_flux_drift_f0_gradient_diffusion"
    )
    assert case.diagnostics["multi_term"]["negative_mass_fraction"] < 1.0e-7
    assert case.diagnostics["multi_term"]["pn_converged"] is True
    assert case.diagnostics["multi_term"]["pn_full_relative_residual"] < 1.0e-6
    assert np.all(np.isfinite(case.eedf))
    assert case.mean_energy_eV > 0.0


def test_multi_term_lmax_changes_the_coupled_solution(tmp_path: Path) -> None:
    cases = []
    for lmax in (1, 4):
        data = base_product_config(tmp_path, ["multi_term"])
        data["solvers"]["multi_term"]["lmax"] = lmax
        cfg = load_config(write_config(tmp_path, data, name=f"lmax_effect_{lmax}.yaml"))
        cases.append(run(cfg, write=False).cases[0])

    low, high = cases
    assert not np.isclose(low.mean_energy_eV, high.mean_energy_eV, rtol=1.0e-4)
    assert not np.isclose(
        low.drift_velocity_m_s,
        high.drift_velocity_m_s,
        rtol=1.0e-4,
    )


def test_multi_term_stationary_solver_fails_closed_on_iteration_limit(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["multi_term"])
    cfg = load_config(write_config(tmp_path, data, name="pn_fail_closed.yaml"))
    cross_sections = load_cross_sections(cfg.cross_sections, cfg.conditions)
    internal = build_internal_solver_configs(cfg.solvers, cfg.physics)
    case = build_multiterm_case(
        cfg,
        cross_sections,
        internal.multi_term,
        cfg.run.e_over_n_Td[0],
    )
    operator = assemble_pn_operator(case)

    with pytest.raises(RuntimeError, match="did not converge"):
        solve_stationary_pn(operator, ConvergenceConfig(max_iterations=1))


def test_multi_term_full_pn_residual_covers_every_moment(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["multi_term"])
    data["solvers"]["multi_term"]["lmax"] = 4
    cfg = load_config(write_config(tmp_path, data, name="pn_residual.yaml"))
    cross_sections = load_cross_sections(cfg.cross_sections, cfg.conditions)
    internal = build_internal_solver_configs(cfg.solvers, cfg.physics)
    case = build_multiterm_case(
        cfg,
        cross_sections,
        internal.multi_term,
        cfg.run.e_over_n_Td[0],
    )
    operator = assemble_pn_operator(case)
    state = solve_stationary_pn(operator, internal.multi_term.convergence)
    flat = np.concatenate(state.raw_coefficients)
    applied = operator.matrix @ flat
    residual = applied - state.diagnostics.growth_frequency_s_inv * flat
    weights = np.concatenate(operator.moment_weights_eV)
    operator_norm = float(
        np.max(np.asarray(abs(operator.matrix).T @ weights).reshape(-1) / weights)
    )
    state_norm = float(np.sum(np.abs(flat) * weights))
    scale = (operator_norm + abs(state.diagnostics.growth_frequency_s_inv)) * state_norm

    block_residuals = [
        float(np.sum(np.abs(residual[block]) * weights) / scale)
        for block, weights in zip(
            operator.moment_slices,
            operator.moment_weights_eV,
            strict=True,
        )
    ]
    assert len(block_residuals) == 5
    assert max(block_residuals) < 1.0e-6
    assert sum(block_residuals) == pytest.approx(
        state.diagnostics.full_pn_relative_residual,
        rel=1.0e-12,
    )


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
    with pytest.raises(ValueError, match="use method 'pn_dcs' for moment_table"):
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
    data["cross_sections"]["files"] = [
        {
            "path": _write_total_momentum_xs(tmp_path).as_posix(),
            "species": "Ar",
            "format": "csv",
        }
    ]
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
    assert set(case.metadata) <= PRODUCT_CASE_METADATA_KEYS
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
    _assert_eedf_normalized(case)
    _assert_canonical_transport(case)


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
    provenance = case.diagnostics["internal_monte_carlo_transport"][
        "mc_run_provenance"
    ]
    assert provenance["eedf_estimator_schema_version"] == (
        MC_MAGNETIC_EEDF_ESTIMATOR_SCHEMA_VERSION
    )
    assert len(str(provenance["solver_source_sha256"])) == 64
    int(str(provenance["solver_source_sha256"]), 16)
    assert case.metadata["ionization_source_treatment"] == "equal"
    assert set(case.metadata) <= PRODUCT_CASE_METADATA_KEYS
    assert "internal_monte_carlo_audit" not in case.diagnostics
    assert case.eedf_counts is not None
    assert int(np.sum(case.eedf_counts)) > 0
    assert case.energy_widths_eV is not None
    assert float(np.sum(case.eedf * case.energy_widths_eV)) == pytest.approx(1.0)


@pytest.mark.mc
def test_internal_monte_carlo_weighted_branching_metadata(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    histogram_thresholds: list[tuple[float, ...]] = []
    histogram_process_counts: list[int] = []

    def capture_histogram_edges(max_energy_eV, thresholds_eV=None, processes=None):
        values = tuple(thresholds_eV or ())
        process_values = tuple(processes or ())
        histogram_thresholds.append(values)
        histogram_process_counts.append(len(process_values))
        return _mc_energy_edges(
            max_energy_eV,
            thresholds_eV=values,
            processes=process_values,
        )

    monkeypatch.setattr(
        "electron_swarm.solvers.monte_carlo.case_state._mc_energy_edges",
        capture_histogram_edges,
    )
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["run"]["e_over_n_Td"] = [400.0]
    data["cross_sections"]["high_energy_extrapolation"] = "hold"
    data["solvers"]["monte_carlo"] = {
        "population_model": "weighted_branching",
        "particles": 12,
        "warmup_collisions": 2,
        "max_collisions": 12,
        "seed": 31,
        "transport_correlation_lag_barriers": 8,
    }
    cfg = load_config(write_config(tmp_path, data))
    [case] = run(cfg, write=False, collect_diagnostics=True).cases

    assert len(histogram_thresholds) == 1
    assert set(histogram_thresholds[0]) == {0.0, 11.55, 15.76}
    assert histogram_process_counts == [3]
    assert set(case.metadata) <= PRODUCT_CASE_METADATA_KEYS
    dev_audit = case.diagnostics["internal_monte_carlo_audit"]
    run_meta = dev_audit["mc_run"]
    energy_meta = dev_audit["mc_energy_audit"]
    assert run_meta["mc_population_model"] == "weighted_branching"
    assert run_meta["secondary_electron_tracking"] is True
    assert run_meta["swarm_population_treatment"] == "weighted_branching_resampled"
    assert run_meta["nonconservative_growth_treatment"] == "explicit_weighted_branching"
    assert run_meta["ionization_branching_model"] == "two_daughter_weighted_resampling"
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
    transport_meta = case.diagnostics["internal_monte_carlo_transport"]
    assert transport_meta["mc_run_provenance"][
        "eedf_estimator_schema_version"
    ] == MC_EEDF_ESTIMATOR_SCHEMA_VERSION
    assert transport_meta["estimator_schema_version"] == (
        WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
    )
    assert transport_meta["transport_definition"] == "flux"
    assert transport_meta["bulk_transport_identified"] is False
    assert transport_meta["cross_gradient_response_identified"] is False
    assert transport_meta["block_lag_sampling"]["physical_time_barriers"] is True
    assert (
        transport_meta["block_lag_sampling"]["resampling_at_common_time_only"] is True
    )
    sampling = transport_meta["block_lag_sampling"]
    assert sampling["configured_correlation_lag_barriers"] == 8
    assert sampling["transport_block_barriers"] == 16
    assert sampling["lag_barriers"] == [2, 4, 8, 16]
    assert transport_meta["lag_scan"]["production_lag_barriers"] == 8
    assert transport_meta["lag_scan"]["hard_convergence_pair_barriers"] == [4, 8]
    assert transport_meta["lag_scan"]["supplemental_pair_barriers"] == [8, 16]
    assert transport_meta["lineage_sampling"]["qualification_lag_barriers"] == [4, 8]
    growth = transport_meta["population_growth_consistency"]
    assert growth["event_count_matches_secondary"] is True
    assert growth["event_exposure_consistent_across_processes"] is True
    assert transport_meta["energy_transport_status"] == "direct_mc_unqualified"
    assert case.transport.electron_energy_mobility_m2_V_s is not None
    assert case.transport.electron_energy_diffusion_L_m2_s is not None
    assert case.transport.electron_energy_diffusion_T_m2_s is not None


@pytest.mark.mc
def test_internal_monte_carlo_weighted_branching_rejects_loss_only(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "population_model": "weighted_branching",
        "particles": 8,
        "max_collisions": 2,
        "seed": 7,
    }
    data["physics"]["ionization"] = {"energy_sharing": "loss_only"}
    cfg = load_config(write_config(tmp_path, data, "weighted_loss_only.yaml"))

    with pytest.raises(NotImplementedError, match="energy_sharing=loss_only"):
        run(cfg, write=False)


@pytest.mark.mc
@pytest.mark.parametrize(
    "population_model",
    ["fixed_particle_single_daughter", "weighted_branching"],
)
def test_internal_monte_carlo_rejects_unimplemented_attachment_dynamics(
    tmp_path: Path,
    population_model: str,
) -> None:
    rows = [
        "species,process,type,threshold_eV,mass_amu,energy_eV,cross_section_m2",
    ]
    for energy in (0.0, 10.0, 100.0, 1000.0):
        rows.append(f"Ar,elastic,elastic,0,39.948,{energy},4e-20")
        rows.append(f"Ar,attachment,attachment,0,39.948,{energy},1e-22")
    xs_path = tmp_path / "attachment_xs.csv"
    xs_path.write_text("\n".join(rows) + "\n", encoding="utf-8")

    data = base_product_config(tmp_path, ["monte_carlo"])
    data["cross_sections"]["files"] = [
        {"path": xs_path.as_posix(), "species": "Ar", "format": "csv"}
    ]
    data["solvers"]["monte_carlo"] = {
        "population_model": population_model,
        "particles": 8,
        "max_collisions": 2,
        "seed": 7,
    }
    cfg = load_config(write_config(tmp_path, data, "attachment.yaml"))

    with pytest.raises(NotImplementedError, match="attachment particle loss"):
        run(cfg, write=False)


def test_weighted_growth_flux_snapshot_uses_position_velocity_covariance() -> None:
    positions = np.array(
        [
            [-1.0, 0.0, -1.0],
            [1.0, 0.0, 1.0],
            [0.0, -1.0, -1.0],
            [0.0, 1.0, 1.0],
        ]
    )
    velocities = positions + np.array([0.0, 0.0, -10.0])
    estimate = weighted_growth_flux_transport_snapshot(
        positions_m=positions,
        velocities_m_s=velocities,
        weights=np.ones(4),
        electric_field_V_m=np.array([0.0, 0.0, 2.0]),
    )

    assert estimate.drift_velocity_m_s == pytest.approx(10.0)
    assert estimate.mobility_m2_V_s == pytest.approx(5.0)
    assert estimate.diffusion_L_m2_s == pytest.approx(1.0)
    assert estimate.diffusion_T_m2_s == pytest.approx(0.5)
    assert np.isfinite(estimate.energy_mobility_m2_V_s)
    assert np.isfinite(estimate.energy_diffusion_L_m2_s)
    assert np.isfinite(estimate.energy_diffusion_T_m2_s)


def test_weighted_growth_v6_sampling_contract() -> None:
    assert WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION == (
        "direct_mc_weighted_growth_configured_lag.v6"
    )
    assert _WEIGHTED_GROWTH_LAG_PLAN.block_barriers == 128
    assert _WEIGHTED_GROWTH_LAG_PLAN.lag_barriers == (16, 32, 64, 128)
    assert _WEIGHTED_GROWTH_LAG_PLAN.production_lag_barriers == 64
    assert _WEIGHTED_GROWTH_LAG_PLAN.lineage_qualification_lag_barriers == (
        32,
        64,
    )
    assert WEIGHTED_MC_TRANSPORT_MIN_EFFECTIVE_LINEAGE_COUNT == 128.0
    assert WEIGHTED_MC_TRANSPORT_MIN_EFFECTIVE_LINEAGE_FRACTION == 0.05

    configured = weighted_growth_transport_lag_plan(256)
    assert configured.lag_barriers == (64, 128, 256, 512)
    assert configured.block_barriers == 512
    assert configured.production_lag_barriers == 256
    assert configured.hard_convergence_pair_barriers == (128, 256)
    assert configured.lineage_qualification_lag_barriers == (128, 256)
    assert configured.supplemental_pair_barriers == (256, 512)

    for invalid in (3, 12, True, 64.0):
        with pytest.raises(ValueError):
            weighted_growth_transport_lag_plan(invalid)


def test_weighted_growth_stationarity_uses_production_residence_ratios() -> None:
    observer = SynchronizedWeightedGrowthFluxObserver(
        particles=4,
        electric_field_V_m=np.array([0.0, 0.0, 2.0]),
        barrier_cadence_s=0.5,
    )
    positions = np.array(
        [
            [-1.0, 0.0, -1.0],
            [1.0, 0.0, 1.0],
            [0.0, -1.0, -1.0],
            [0.0, 1.0, 1.0],
        ]
    )
    for block in range(128):
        if block < 32:
            mobility, mean_energy, residence_time = 1.0, 2.0, 1.0
        elif block < 64:
            mobility, mean_energy, residence_time = 2.0, 4.0, 3.0
        else:
            mobility, mean_energy, residence_time = 3.0, 6.0, 1.0
        observer.record_residence(
            displacement_m=np.array([0.0, 0.0, -2.0 * mobility * residence_time]),
            sample_energy_eV=mean_energy,
            dt_s=residence_time,
            weight=1.0,
        )
        for lag in _WEIGHTED_GROWTH_LAG_PLAN.lag_barriers:
            lag_positions = (
                positions * 10.0
                if lag == _WEIGHTED_GROWTH_LAG_PLAN.block_barriers
                else positions
            )
            observer.record_plane(
                time_s=float(block * _WEIGHTED_GROWTH_LAG_PLAN.block_barriers + lag),
                positions_m=lag_positions,
                velocities_m_s=(lag_positions + np.array([0.0, 0.0, -10.0])),
                weights=np.ones(4),
                lineages=np.arange(4),
                lag_barriers=lag,
            )

    estimate, diagnostics = observer.finalize()

    assert estimate.mobility_m2_V_s == pytest.approx(13.0 / 6.0)
    assert estimate.energy_mobility_m2_V_s == pytest.approx(31.0 / 13.0)
    assert estimate.mean_energy_eV == pytest.approx(13.0 / 3.0)
    early = diagnostics["time_stationarity_window_estimate_early"]
    late = diagnostics["time_stationarity_window_estimate_late"]
    assert early["mobility_m2_V_s"] == pytest.approx(1.75)
    assert late["mobility_m2_V_s"] == pytest.approx(3.0)
    assert early["energy_mobility_m2_V_s"] == pytest.approx(13.0 / 7.0)
    assert late["energy_mobility_m2_V_s"] == pytest.approx(3.0)
    assert early["mean_energy_eV"] == pytest.approx(3.5)
    assert late["mean_energy_eV"] == pytest.approx(6.0)
    lag_scan = diagnostics["lag_scan"]
    production_index = lag_scan["lag_barriers"].index(
        _WEIGHTED_GROWTH_LAG_PLAN.production_lag_barriers
    )
    diffusion_fields = (
        "diffusion_L_m2_s",
        "diffusion_T_m2_s",
        "energy_diffusion_L_m2_s",
        "energy_diffusion_T_m2_s",
    )
    for field in diffusion_fields:
        production_value = lag_scan["estimates"][field][production_index]
        assert getattr(estimate, field) == pytest.approx(production_value)
        assert early[field] == pytest.approx(production_value)
        assert late[field] == pytest.approx(production_value)
    assert estimate.diffusion_L_m2_s != pytest.approx(
        lag_scan["estimates"]["diffusion_L_m2_s"][-1]
    )
    assert diagnostics["stationarity_sampling"] == {
        "window_definition": ("first_half_vs_second_half_complete_transport_blocks"),
        "state_moment_aggregation": ("sum_raw_weighted_residence_moments_before_ratio"),
        "state_estimator": "same_residence_ratio_formulas_as_production",
        "diffusion_estimator": (
            "mean_fixed_production_lag_snapshot_covariance_per_window"
        ),
        "diffusion_lag_barriers": (_WEIGHTED_GROWTH_LAG_PLAN.production_lag_barriers),
        "early_block_count": 64,
        "late_block_count": 64,
        "odd_center_block_excluded": False,
    }


def test_weighted_growth_observer_fails_closed_on_lineage_collapse() -> None:
    observer = SynchronizedWeightedGrowthFluxObserver(
        particles=4,
        electric_field_V_m=np.array([0.0, 0.0, 2.0]),
        barrier_cadence_s=0.5,
    )
    positions = np.array(
        [
            [-1.0, 0.0, -1.0],
            [1.0, 0.0, 1.0],
            [0.0, -1.0, -1.0],
            [0.0, 1.0, 1.0],
        ]
    )
    velocities = positions + np.array([0.0, 0.0, -10.0])
    observer.record_residence(
        displacement_m=np.array([0.0, 0.0, -1.0]),
        sample_energy_eV=2.0,
        dt_s=1.0,
        weight=1.0,
    )
    for block in range(128):
        if block > 0:
            observer.record_residence(
                displacement_m=np.array([0.0, 0.0, -1.0]),
                sample_energy_eV=2.0,
                dt_s=1.0,
                weight=1.0,
            )
        for lag in _WEIGHTED_GROWTH_LAG_PLAN.lag_barriers:
            observer.record_plane(
                time_s=float(block * _WEIGHTED_GROWTH_LAG_PLAN.block_barriers + lag),
                positions_m=positions,
                velocities_m_s=velocities,
                weights=np.ones(4),
                lineages=np.zeros(4, dtype=int),
                lag_barriers=lag,
            )
    _, diagnostics = observer.finalize()

    assert diagnostics["lineage_qualified"] is False
    lineage = diagnostics["lineage_sampling"]
    production_index = lineage["lag_barriers"].index(
        _WEIGHTED_GROWTH_LAG_PLAN.production_lag_barriers
    )
    assert lineage["minimum_effective_lineage_fraction"][
        production_index
    ] == pytest.approx(0.25)
    assert lineage["qualification_lag_barriers"] == [32, 64]
    assert lineage["required_effective_lineage_count"] == 128.0
    assert lineage["required_effective_lineage_fraction"] == 0.05
    assert diagnostics["energy_transport_status"] == "direct_mc_unqualified"


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
    assert zero_b.metadata["magnetic_field_treatment"] == "none"

    disabled_nonzero = base_product_config(tmp_path, ["monte_carlo"])
    disabled_nonzero["solvers"]["monte_carlo"] = dict(base["solvers"]["monte_carlo"])
    disabled_nonzero["physics"]["field"]["magnetic_field"] = {
        "enabled": False,
        "B_T": 0.2,
        "angle_EB_deg": 90.0,
    }
    disabled_b = run(
        load_config(
            write_config(tmp_path, disabled_nonzero, "mc_disabled_nonzero_b.yaml")
        ),
        write=False,
    ).cases[0]

    assert disabled_b.mean_energy_eV == disabled.mean_energy_eV
    assert disabled_b.drift_velocity_m_s == disabled.drift_velocity_m_s
    np.testing.assert_array_equal(disabled_b.eedf, disabled.eedf)
    assert disabled_b.metadata["magnetic_field_treatment"] == "none"
    provenance = disabled_b.diagnostics["internal_monte_carlo_transport"][
        "mc_run_provenance"
    ]
    assert provenance["magnetic_B_T"] == 0.0


def test_external_monte_carlo_product_fields_are_rejected(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "python_api": "product_helpers:fake_mc_missing",
    }
    with pytest.raises(ValueError, match="Unsupported solvers.monte_carlo fields"):
        load_config(write_config(tmp_path, data))
