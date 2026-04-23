from pathlib import Path

import numpy as np

from electron_swarm import load_config, run
from electron_swarm.core.cross_sections import ProcessType, load_cross_sections, mixture_fraction
from electron_swarm.solvers.boltzmann_two_term import (
    BoltzmannTwoTermSolver,
    _electron_speed,
    _gas_number_density,
)

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "configs" / "unified" / "boltzmann_only.yaml"


def _fast_boltzmann_config(tmp_path: Path):
    cfg = load_config(CONFIG_PATH)
    cfg.output.directory = tmp_path
    cfg.output.write_plots = False
    cfg.run.e_over_n_Td = [50.0]
    cfg.boltzmann_two_term.backend = "native_bolsig"
    cfg.boltzmann_two_term.energy_grid.n = 120
    cfg.boltzmann_two_term.energy_grid.max_eV = 50.0
    cfg.boltzmann_two_term.adaptive_grid.enabled = False
    cfg.boltzmann_two_term.convergence.max_iterations = 80
    return cfg


def test_native_bolsig_runs_and_normalizes(tmp_path: Path):
    cfg = _fast_boltzmann_config(tmp_path)
    result = run(cfg, write=True)
    assert len(result.cases) == 1
    case = result.cases[0]
    assert case.metadata["backend"] == "native_bolsig"
    assert case.metadata["converged"] is True
    assert case.metadata["net_ionization_frequency_model"] == "convolution_effective_rate"
    assert case.mean_energy_eV > 0.0
    assert case.drift_velocity_m_s > 0.0
    assert case.reduced_mobility_m2_V_s_m3 > 0.0
    assert case.reduced_diffusion_L_m2_s_m3 > 0.0
    assert np.all(np.isfinite(case.eedf))
    assert np.trapezoid(case.eedf, case.energy_eV) > 0.8
    assert (tmp_path / "argon_summary.csv").exists()
    assert (tmp_path / "argon_eedf.csv").exists()
    assert (tmp_path / "argon_rates.csv").exists()
    assert (tmp_path / "summary_boltzmann.csv").exists()
    assert (tmp_path / "summary.csv").exists()


def test_e_over_n_scan_produces_monotonic_case_ids(tmp_path: Path):
    cfg = _fast_boltzmann_config(tmp_path)
    cfg.run.e_over_n_Td = [20.0, 100.0]
    cfg.boltzmann_two_term.energy_grid.n = 100
    result = run(cfg, write=False)
    assert [c.e_over_n_Td for c in result.cases] == [20.0, 100.0]
    assert all(c.solver == "boltzmann_two_term" for c in result.cases)
    assert result.cases[1].mean_energy_eV > result.cases[0].mean_energy_eV


def test_internal_backend_aliases_native_bolsig(tmp_path: Path):
    cfg = _fast_boltzmann_config(tmp_path)
    cfg.boltzmann_two_term.backend = "internal"
    cfg.boltzmann_two_term.energy_grid.n = 80
    result = run(cfg, write=False)
    assert result.cases[0].metadata["backend"] == "native_bolsig"


def test_effective_transport_frequency_includes_inelastic_cross_sections(tmp_path: Path):
    cfg = _fast_boltzmann_config(tmp_path)
    cross_sections = load_cross_sections(cfg.cross_sections, cfg.conditions)
    solver = BoltzmannTwoTermSolver(cfg, cross_sections)
    energy = np.array([15.0, 15.5])
    speed = _electron_speed(energy)
    coll = solver._effective_collision_data(energy, _gas_number_density(cfg))

    expected_sigma = np.zeros_like(energy)
    for proc in cross_sections.processes:
        frac = mixture_fraction(cfg.conditions, proc.species)
        if frac <= 0.0:
            continue
        if proc.process_type not in {
            ProcessType.MOMENTUM,
            ProcessType.ELASTIC,
            ProcessType.EXCITATION,
            ProcessType.IONIZATION,
            ProcessType.ATTACHMENT,
            ProcessType.SUPERELASTIC,
        }:
            continue
        expected_sigma += frac * np.maximum(
            proc.sigma(energy), cfg.boltzmann_two_term.min_momentum_cross_section_m2
        )

    assert np.all(expected_sigma > 0.0)
    assert np.allclose(coll.sigma_m, expected_sigma)
    assert np.allclose(coll.nu_m_over_N, expected_sigma * speed)
