from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

from electron_swarm import load_config, run
from electron_swarm.core.cross_sections import (
    CrossSectionProcess,
    CrossSectionSet,
    ProcessType,
    load_cross_sections,
)
from electron_swarm.core.transport import FluxTransport, BulkTransport, TransportSet
from electron_swarm.solvers.boltzmann_two_term import BoltzmannTwoTermSolver
from electron_swarm.solvers.multiterm_boltzmann import MultiTermBoltzmannSolver
from electron_swarm.solvers.multiterm_boltzmann import (
    EnergyGrid,
    _compute_rate_set,
    _elastic_power_loss_eV_s,
    _project_collision_data,
    build_multiterm_case,
)
from electron_swarm.solvers.multiterm_boltzmann.operator import (
    assemble_lmax1_native_operator_block,
    compare_lmax1_transport_to_reference,
    solve_lmax1_two_term_reference,
)

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "configs" / "unified" / "multiterm_boltzmann.yaml"
OPERATOR_LMAX1_CONFIG_PATH = (
    ROOT / "configs" / "unified" / "multiterm_operator_lmax1.yaml"
)
OPERATOR_CONFIG_PATH = (
    ROOT / "configs" / "unified" / "multiterm_operator_experimental.yaml"
)
ALL_CONFIG_PATH = ROOT / "configs" / "unified" / "all_template.yaml"
EXPORTER_ROOT = ROOT / "swarm_comsol_exporter"
if str(EXPORTER_ROOT) not in sys.path:
    sys.path.insert(0, str(EXPORTER_ROOT))
from swarm_comsol_export.hook import maybe_export_comsol  # noqa: E402
from swarm_comsol_export.discovery import discover_inputs  # noqa: E402


def _fast_multiterm_config(tmp_path: Path):
    cfg = load_config(CONFIG_PATH)
    cfg.output.directory = tmp_path
    cfg.output.write_plots = False
    cfg.run.e_over_n_Td = [50.0]
    cfg.multiterm_boltzmann.method = "moment_closure"
    cfg.multiterm_boltzmann.hydrodynamic = False
    cfg.multiterm_boltzmann.energy_grid.n = 80
    cfg.multiterm_boltzmann.energy_grid.max_eV = 60.0
    return cfg


def test_transport_requires_explicit_bulk():
    flux = FluxTransport(10.0, 2.0, 4.0)
    tset = TransportSet.from_flux_only(flux)
    with pytest.raises(RuntimeError):
        tset.require_bulk()
    bulk = BulkTransport(12.0, 2.4, 5.0)
    both = TransportSet.from_flux_bulk(flux, bulk, 7.0, 2.0)
    assert both.source.gradient_velocity_m_s == pytest.approx(2.0)
    assert both.source.curvature_diffusion_longitudinal_m2_s == pytest.approx(1.0)


def test_energy_grid_shapes_normalize():
    for grid in (
        EnergyGrid.linear(0.0, 10.0, 32),
        EnergyGrid.log(1.0e-3, 10.0, 32),
        EnergyGrid.log_linear(1.0e-3, 20.0, 48),
    ):
        assert np.all(np.diff(grid.edges_eV) > 0.0)
        pdf = grid.normalize_energy_pdf(np.exp(-grid.centers_eV))
        assert np.sum(pdf * grid.widths_eV) == pytest.approx(1.0)


def test_cross_section_conversion_from_existing_loader(tmp_path: Path):
    cfg = _fast_multiterm_config(tmp_path)
    cross_sections = load_cross_sections(cfg.cross_sections, cfg.conditions)
    case = build_multiterm_case(cfg, cross_sections, 50.0)
    assert case.grid.n_cells == 80
    assert case.cross_sections.processes
    assert case.gas_number_density_m3 > 0.0
    assert case.electric_field_V_m > 0.0


def test_superelastic_power_loss_is_energy_gain(tmp_path: Path):
    cfg = _fast_multiterm_config(tmp_path)
    proc = CrossSectionProcess(
        species="Ar",
        process="deexcitation",
        process_type=ProcessType.SUPERELASTIC,
        threshold_eV=11.5,
        energy_eV=np.array([0.0, 20.0]),
        cross_section_m2=np.array([1.0e-20, 1.0e-20]),
    )
    case = build_multiterm_case(cfg, CrossSectionSet([proc]), 50.0)
    F = case.grid.normalize_energy_pdf(np.exp(-case.grid.centers_eV))
    rates = _compute_rate_set(case, _project_collision_data(case), F, "super")
    assert rates.rates[0].power_loss_eV_s < 0.0


def test_elastic_power_loss_uses_cell_energy(tmp_path: Path):
    cfg = _fast_multiterm_config(tmp_path)
    proc = CrossSectionProcess(
        species="Ar",
        process="momentum",
        process_type=ProcessType.MOMENTUM,
        threshold_eV=None,
        mass_amu=39.948,
        energy_eV=np.array([0.0, 20.0]),
        cross_section_m2=np.array([1.0e-20, 3.0e-20]),
    )
    case = build_multiterm_case(cfg, CrossSectionSet([proc]), 50.0)
    collisions = _project_collision_data(case)
    F = case.grid.normalize_energy_pdf(1.0 + case.grid.centers_eV)
    thermal_mean = 0.5
    actual = _elastic_power_loss_eV_s(case, collisions, F, thermal_mean)
    mean_e = float(np.sum(case.grid.centers_eV * F * case.grid.widths_eV))
    old_style = float(
        np.sum(
            collisions.elastic_energy_frequency_s_inv
            * (mean_e - thermal_mean)
            * F
            * case.grid.widths_eV
        )
    )
    assert actual != pytest.approx(old_style)


def test_effective_cross_section_suppresses_elastic_momentum_double_count(tmp_path: Path):
    cfg = _fast_multiterm_config(tmp_path)
    energy = np.array([0.0, 10.0])
    effective = CrossSectionProcess(
        species="Ar",
        process="effective",
        process_type=ProcessType.EFFECTIVE,
        threshold_eV=None,
        energy_eV=energy,
        cross_section_m2=np.array([2.0e-20, 2.0e-20]),
    )
    elastic = CrossSectionProcess(
        species="Ar",
        process="elastic",
        process_type=ProcessType.ELASTIC,
        threshold_eV=None,
        energy_eV=energy,
        cross_section_m2=np.array([9.0e-20, 9.0e-20]),
    )
    case = build_multiterm_case(cfg, CrossSectionSet([effective, elastic]), 50.0)
    collisions = _project_collision_data(case)
    assert collisions.normalization.effective_species == ("Ar",)
    assert collisions.normalization.suppressed_momentum_processes == (
        "Ar:elastic:elastic",
    )
    np.testing.assert_allclose(
        collisions.momentum_frequency_s_inv,
        collisions.processes[0].frequency_s_inv,
    )


def test_lmax1_native_operator_block_is_public_and_sparse(tmp_path: Path):
    cfg = _fast_multiterm_config(tmp_path)
    cfg.boltzmann_two_term.energy_grid.n = 64
    cfg.boltzmann_two_term.energy_grid.max_eV = 35.0
    cfg.boltzmann_two_term.adaptive_grid.enabled = False
    cross_sections = load_cross_sections(cfg.cross_sections, cfg.conditions)

    block = assemble_lmax1_native_operator_block(cfg, cross_sections, 50.0)

    assert block.discretization == "finite_volume_scharfetter_gummel"
    assert block.matrix.shape == (64, 64)
    assert block.matrix.nnz > 64
    assert block.energy_eV.shape == block.widths_eV.shape
    assert block.edges_eV.shape == (65,)
    assert block.collisions.nu_m.shape == block.energy_eV.shape
    assert block.electric_field_V_m > 0.0
    assert block.gas_number_density_m3 > 0.0


def test_two_term_native_distribution_is_reusable_for_operator(tmp_path: Path):
    cfg = _fast_multiterm_config(tmp_path)
    cfg.boltzmann_two_term.energy_grid.n = 80
    cfg.boltzmann_two_term.energy_grid.max_eV = 40.0
    cfg.boltzmann_two_term.adaptive_grid.enabled = False
    cross_sections = load_cross_sections(cfg.cross_sections, cfg.conditions)
    solver = BoltzmannTwoTermSolver(cfg, cross_sections)
    block = solver.solve_native_distribution(50.0)
    case = solver.solve_native_reference_case(50.0, "ref")
    assert block.energy_eV.shape == block.eedf_eV_inv.shape
    assert block.widths_eV.shape == block.energy_eV.shape
    assert np.sum(block.eedf_eV_inv * block.widths_eV) == pytest.approx(1.0)
    assert block.metadata["discretization"] == "finite_volume_scharfetter_gummel"
    assert case.solver == "boltzmann_two_term"
    assert case.mean_energy_eV == pytest.approx(
        float(np.sum(block.energy_eV * block.eedf_eV_inv * block.widths_eV))
    )


def test_lmax1_two_term_reference_comparison_harness(tmp_path: Path):
    cfg = _fast_multiterm_config(tmp_path)
    cfg.multiterm_boltzmann.lmax = 1
    cfg.multiterm_boltzmann.method = "operator"
    cfg.boltzmann_two_term.energy_grid.n = 80
    cfg.boltzmann_two_term.energy_grid.max_eV = 40.0
    cfg.boltzmann_two_term.adaptive_grid.enabled = False
    cross_sections = load_cross_sections(cfg.cross_sections, cfg.conditions)
    candidate = MultiTermBoltzmannSolver(cfg, cross_sections).solve_case(
        50.0, "candidate"
    )
    assert candidate.metadata["multiterm_method_used"] == "operator_lmax1_two_term"
    assert candidate.metadata["backend"] == "native_operator_lmax1_two_term"
    assert candidate.metadata["physical_validity"] == "operator_lmax1_two_term_reference"
    assert candidate.transport is not None
    assert candidate.transport.bulk is None
    reference = solve_lmax1_two_term_reference(
        cfg, cross_sections, 50.0, "reference"
    )
    report = compare_lmax1_transport_to_reference(
        candidate,
        reference,
    )
    assert report.candidate_solver == "multiterm_boltzmann"
    assert report.reference_solver == "boltzmann_two_term"
    assert {item.metric for item in report.comparisons} == {
        "mean_energy_eV",
        "drift_velocity_m_s",
        "diffusion_L_m2_s",
        "net_ionization_frequency_s",
    }
    assert report.passed
    assert all(np.isfinite(item.relative_difference) for item in report.comparisons)
    metadata = report.as_metadata()
    assert "lmax1_reference_passed" in metadata


def test_operator_lmax1_config_runs_no_write(tmp_path: Path):
    cfg = load_config(OPERATOR_LMAX1_CONFIG_PATH)
    cfg.output.directory = tmp_path
    cfg.output.write_plots = False
    cfg.run.e_over_n_Td = [50.0]
    cfg.boltzmann_two_term.energy_grid.n = 80
    cfg.boltzmann_two_term.energy_grid.max_eV = 40.0
    cfg.boltzmann_two_term.adaptive_grid.enabled = False
    result = run(cfg, write=False)
    assert len(result.cases) == 1
    case = result.cases[0]
    assert case.solver == "multiterm_boltzmann"
    assert case.metadata["multiterm_method_used"] == "operator_lmax1_two_term"
    assert np.isfinite(case.mean_energy_eV)


def test_operator_config_runs_and_writes(tmp_path: Path):
    cfg = load_config(OPERATOR_CONFIG_PATH)
    cfg.output.directory = tmp_path
    cfg.output.write_plots = False
    result = run(cfg, write=True)
    assert len(result.cases) == 1
    case = result.cases[0]
    assert case.metadata["multiterm_method_used"] == "operator_reference_anchored_lmax_gt1"
    assert case.metadata["transport_definition"] == "flux"
    assert case.metadata["operator_coefficient_order"] == "legendre_major_energy_minor"
    assert (tmp_path / "argon_operator_experimental_summary.csv").exists()
    assert (tmp_path / "summary_multiterm.csv").exists()
    assert (tmp_path / "eedf_table_multiterm.csv").exists()


def test_operator_lmax_gt1_requires_experimental_opt_in(tmp_path: Path):
    cfg = _fast_multiterm_config(tmp_path)
    cfg.multiterm_boltzmann.method = "operator"
    cfg.multiterm_boltzmann.lmax = 2
    with pytest.raises(RuntimeError, match="allow_experimental_operator"):
        run(cfg, write=False)


def test_operator_runs_flux_only(tmp_path: Path):
    cfg = _fast_multiterm_config(tmp_path)
    cfg.multiterm_boltzmann.method = "operator"
    cfg.multiterm_boltzmann.allow_experimental_operator = True
    cfg.multiterm_boltzmann.lmax = 2
    cfg.multiterm_boltzmann.energy_grid.n = 80
    cfg.multiterm_boltzmann.energy_grid.max_eV = 60.0
    result = run(cfg, write=False)
    case = result.cases[0]
    assert case.metadata["multiterm_method_used"] == "operator_reference_anchored_lmax_gt1"
    assert (
        case.metadata["physical_validity"]
        == "reference_anchored_lmax_gt1_integral_cross_section_closure"
    )
    assert case.metadata["operator_validation_required"] is True
    assert case.metadata["operator_lmax"] == 2
    assert case.transport is not None
    assert case.transport.bulk is None
    assert np.isfinite(case.mean_energy_eV)
    assert np.isfinite(case.drift_velocity_m_s)
    assert np.isfinite(case.diffusion_L_m2_s)
    assert case.metadata["normalization_integral"] == pytest.approx(1.0)
    assert np.isfinite(case.metadata["operator_tail_rate_fraction"])
    assert np.isfinite(case.metadata["operator_highest_l_relative_l1"])
    assert isinstance(case.metadata["operator_lmax_convergence_ok"], bool)
    assert (
        case.metadata["operator_ionization_source_model"]
        == "two_term_native_equal_sharing"
    )
    assert (
        case.metadata["operator_l_gt_0_inelastic_model"]
        == "sink_only_isotropic_l0_source"
    )
    assert (
        case.metadata["operator_l_gt_0_elastic_model"]
        == "reference_anchored_integral_momentum_closure"
    )
    assert (
        case.metadata["operator_reference_gate_status"]
        == "anchored_to_native_two_term"
    )
    assert case.metadata["operator_reference_relerr_max"] == pytest.approx(0.0)
    assert case.metadata["operator_transport_reused_from"] == "native_two_term_reference"


def test_operator_synthetic_attachment_and_superelastic_are_finite(tmp_path: Path):
    cfg = _fast_multiterm_config(tmp_path)
    cfg.multiterm_boltzmann.method = "operator"
    cfg.multiterm_boltzmann.allow_experimental_operator = True
    cfg.multiterm_boltzmann.lmax = 2
    cfg.multiterm_boltzmann.energy_grid.n = 48
    cfg.multiterm_boltzmann.energy_grid.max_eV = 40.0
    energy = np.array([0.0, 3.0, 10.0, 40.0])
    momentum = CrossSectionProcess(
        species="Ar",
        process="momentum",
        process_type=ProcessType.MOMENTUM,
        threshold_eV=None,
        mass_amu=39.948,
        energy_eV=energy,
        cross_section_m2=np.array([1.0e-20, 1.0e-20, 1.2e-20, 1.2e-20]),
    )
    attachment = CrossSectionProcess(
        species="Ar",
        process="attachment",
        process_type=ProcessType.ATTACHMENT,
        threshold_eV=0.0,
        energy_eV=energy,
        cross_section_m2=np.array([0.0, 4.0e-22, 4.0e-22, 1.0e-22]),
    )
    excitation = CrossSectionProcess(
        species="Ar",
        process="excitation",
        process_type=ProcessType.EXCITATION,
        threshold_eV=1.5,
        energy_eV=energy,
        cross_section_m2=np.array([0.0, 1.0e-22, 1.0e-22, 1.0e-22]),
    )
    superelastic = CrossSectionProcess(
        species="Ar",
        process="deexcitation",
        process_type=ProcessType.SUPERELASTIC,
        threshold_eV=1.5,
        energy_eV=energy,
        cross_section_m2=np.array([5.0e-26, 5.0e-26, 5.0e-26, 5.0e-26]),
    )
    cases = [
        MultiTermBoltzmannSolver(
            cfg, CrossSectionSet([momentum, attachment])
        ).solve_case(40.0, "synthetic_attachment"),
        MultiTermBoltzmannSolver(
            cfg, CrossSectionSet([momentum, excitation, superelastic])
        ).solve_case(40.0, "synthetic_superelastic"),
    ]

    for case in cases:
        assert np.isfinite(case.mean_energy_eV)
        assert np.isfinite(case.drift_velocity_m_s)
        assert case.metadata["normalization_integral"] == pytest.approx(1.0)
        assert case.metadata["operator_l_gt_0_inelastic_model"] == (
            "sink_only_isotropic_l0_source"
        )
    assert np.isfinite(
        cases[0].metadata["convolution_attachment_rate_coefficient_m3_s"]
    )
    assert (
        cases[1].metadata["operator_ionization_source_model"]
        == "two_term_native_equal_sharing"
    )


def test_operator_hydrodynamic_request_is_marked_not_computed(tmp_path: Path):
    cfg = _fast_multiterm_config(tmp_path)
    cfg.multiterm_boltzmann.method = "operator"
    cfg.multiterm_boltzmann.allow_experimental_operator = True
    cfg.multiterm_boltzmann.hydrodynamic = True
    cfg.multiterm_boltzmann.lmax = 2
    cfg.multiterm_boltzmann.energy_grid.n = 80
    cfg.multiterm_boltzmann.energy_grid.max_eV = 60.0
    result = run(cfg, write=False)
    case = result.cases[0]
    assert case.metadata["multiterm_method_used"] == "operator_reference_anchored_lmax_gt1"
    assert case.metadata["transport_definition"] == "flux"
    assert case.metadata["physical_validity"] == (
        "reference_anchored_lmax_gt1_integral_cross_section_closure"
    )
    assert case.metadata["operator_validation_required"] is True
    assert case.metadata["operator_hydrodynamic_computed"] is False
    assert "operator_hydrodynamic_requested_but_not_computed" in case.metadata["warnings"]
    assert case.transport is not None
    assert case.transport.bulk is None
    assert np.isfinite(case.drift_velocity_m_s)
    assert np.isfinite(case.diffusion_L_m2_s)


def test_operator_hydrodynamic_rejects_nonphysical_field_scaling(tmp_path: Path):
    cfg = _fast_multiterm_config(tmp_path)
    cfg.multiterm_boltzmann.method = "operator"
    cfg.multiterm_boltzmann.allow_experimental_operator = True
    cfg.multiterm_boltzmann.hydrodynamic = True
    cfg.multiterm_boltzmann.lmax = 2
    cfg.multiterm_boltzmann.field_coupling_scale = 0.5
    result = run(cfg, write=False)
    case = result.cases[0]
    assert "operator_nonphysical_field_scaling" in case.metadata["warnings"]


def test_hybrid_method_fails_fast(tmp_path: Path):
    cfg = _fast_multiterm_config(tmp_path)
    cfg.multiterm_boltzmann.method = "hybrid"
    with pytest.raises(NotImplementedError):
        run(cfg, write=False)


def test_power_balance_failure_is_not_silent(tmp_path: Path):
    cfg = _fast_multiterm_config(tmp_path)
    with pytest.raises(RuntimeError):
        MultiTermBoltzmannSolver(cfg, CrossSectionSet([])).solve_case(50.0, "empty")


def test_multiterm_runs_and_writes_outputs(tmp_path: Path):
    cfg = _fast_multiterm_config(tmp_path)
    result = run(cfg, write=True)
    assert len(result.cases) == 1
    case = result.cases[0]
    assert case.solver == "multiterm_boltzmann"
    assert case.schema_version == "1.2"
    assert np.isfinite(case.mean_energy_eV) and case.mean_energy_eV > 0.0
    assert np.isfinite(case.drift_velocity_m_s)
    assert case.transport is not None
    assert case.transport.bulk is None
    assert case.metadata["multiterm_method_used"] == "moment_closure"
    assert case.metadata["transport_definition"] == "flux"
    assert "not_multiterm_operator" in case.metadata["physical_validity"]
    assert case.metadata["cross_section_high_energy_extrapolation"] == "zero"
    cross_sections = load_cross_sections(cfg.cross_sections, cfg.conditions)
    mt_case = build_multiterm_case(cfg, cross_sections, 50.0)
    assert np.sum(case.eedf * mt_case.grid.widths_eV) == pytest.approx(1.0)

    summary = pd.read_csv(tmp_path / "argon_multiterm_summary.csv")
    assert "multiterm_boltzmann" in set(summary["solver"])
    assert "meta_transport_definition" in summary.columns
    assert "meta_physical_validity" in summary.columns
    assert "meta_estimated_bulk_drift_velocity_m_s" in summary.columns
    assert (tmp_path / "summary_multiterm.csv").exists()
    assert (tmp_path / "eedf_table_multiterm.csv").exists()
    assert (tmp_path / "energy_table_multiterm.csv").exists()
    assert (tmp_path / "summary.csv").exists()


def test_all_mode_dispatches_two_term_and_multiterm(tmp_path: Path):
    cfg = load_config(ALL_CONFIG_PATH)
    cfg.output.directory = tmp_path
    cfg.output.write_plots = False
    cfg.run.e_over_n_Td = [50.0]
    cfg.monte_carlo.enabled = False
    cfg.boltzmann_two_term.energy_grid.n = 80
    cfg.boltzmann_two_term.energy_grid.max_eV = 40.0
    cfg.boltzmann_two_term.adaptive_grid.enabled = False
    cfg.multiterm_boltzmann.method = "operator"
    cfg.multiterm_boltzmann.allow_experimental_operator = True
    cfg.multiterm_boltzmann.lmax = 2
    cfg.multiterm_boltzmann.energy_grid.n = 60
    cfg.multiterm_boltzmann.energy_grid.max_eV = 50.0
    result = run(cfg, write=True)
    assert result.cases[0].transport is not None
    assert result.cases[0].transport.bulk is None
    assert result.cases[1].transport is not None
    assert result.cases[1].transport.bulk is None
    assert [case.solver for case in result.cases] == [
        "boltzmann_two_term",
        "multiterm_boltzmann",
    ]
    summary = pd.read_csv(tmp_path / "argon_all_summary.csv")
    assert {"boltzmann_two_term", "multiterm_boltzmann"} <= set(summary["solver"])
    mt = summary[summary["solver"] == "multiterm_boltzmann"].iloc[0]
    assert mt["meta_multiterm_method_used"] == "operator_reference_anchored_lmax_gt1"


def test_comsol_export_compat_multiterm(tmp_path: Path):
    cfg = _fast_multiterm_config(tmp_path)
    run(cfg, write=True)
    out = maybe_export_comsol(
        {
            "comsol_export": {
                "enabled": True,
                "input": {
                    "transport_csv": "summary.csv",
                    "rates_csv": "summary.csv",
                    "eedf": {
                        "mode": "stacked_csv",
                        "stacked_csv": "eedf_table.csv",
                    },
                },
            }
        },
        run_dir=tmp_path,
    )
    assert out is not None
    assert (out / "comsol_transport.csv").exists()
    assert (out / "comsol_rates.csv").exists()
    assert (out / "comsol_eedf.csv").exists()


def test_comsol_discovery_requires_primary_solver_for_ambiguous_fallback(tmp_path: Path):
    (tmp_path / "summary_boltzmann.csv").write_text("solver\n", encoding="utf-8")
    (tmp_path / "summary_multiterm.csv").write_text("solver\n", encoding="utf-8")
    (tmp_path / "eedf_table_multiterm.csv").write_text("energy\n", encoding="utf-8")
    with pytest.raises(ValueError, match="primary_solver"):
        discover_inputs(tmp_path, {})
    found = discover_inputs(tmp_path, {"primary_solver": "multiterm_boltzmann"})
    assert found.transport_csv == tmp_path / "summary_multiterm.csv"
