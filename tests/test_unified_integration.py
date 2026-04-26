from pathlib import Path
import sys
import warnings

import numpy as np
import pandas as pd
import pytest

from electron_swarm import load_config, run
from electron_swarm.core.cross_sections import (
    ProcessType,
    load_cross_sections,
)
from electron_swarm.solvers.monte_carlo_adapter import MonteCarloAdapter
from EEDF_optimze import load_config as load_optimize_config
from EEDF_optimze import process_table
from swarm_mc.config import Config as McConfig
from swarm_mc.cross_section import InterpolatedCrossSectionSet
from swarm_mc.monte_carlo import MonteCarlo
from swarm_mc.unified_adapter import _build_mc_config

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "configs" / "unified" / "both_template.yaml"
EXPORTER_ROOT = ROOT / "swarm_comsol_exporter"
if str(EXPORTER_ROOT) not in sys.path:
    sys.path.insert(0, str(EXPORTER_ROOT))
from swarm_comsol_export.hook import maybe_export_comsol  # noqa: E402


def _fast_both_config(tmp_path: Path):
    cfg = load_config(CONFIG_PATH)
    cfg.output.directory = tmp_path
    cfg.output.write_plots = False
    cfg.run.e_over_n_Td = [50.0]
    cfg.boltzmann_two_term.energy_grid.n = 100
    cfg.boltzmann_two_term.energy_grid.max_eV = 40.0
    cfg.boltzmann_two_term.adaptive_grid.enabled = False
    cfg.boltzmann_two_term.convergence.max_iterations = 80
    base = cfg.monte_carlo.passthrough["base_config"]
    base["initial_state"]["num_e_initial"] = 1500
    base["simulation_settings"]["num_energy_bins"] = 128
    base["simulation_settings"]["num_e_max"] = 80000
    base["simulation_settings"]["seed"] = 1
    base["end_conditions"]["num_col_max"] = 12000
    base["end_conditions"]["w_tol"] = 0.15
    base["end_conditions"]["DN_tol"] = 0.15
    return cfg


@pytest.fixture(scope="module")
def integrated_outputs(tmp_path_factory):
    out_dir = tmp_path_factory.mktemp("unified_both")
    cfg = _fast_both_config(out_dir)
    result = run(cfg, write=True)
    return cfg, result, out_dir


def test_cross_section_txt_loader():
    cfg = load_config(ROOT / "configs" / "unified" / "boltzmann_only.yaml")
    loaded = load_cross_sections(cfg.cross_sections, cfg.conditions)
    assert loaded.processes
    assert cfg.cross_sections.high_energy_extrapolation == "zero"
    assert any(
        process.process_type
        in {ProcessType.MOMENTUM, ProcessType.ELASTIC, ProcessType.EFFECTIVE}
        for process in loaded.processes
    )
    assert any(process.process_type == ProcessType.IONIZATION for process in loaded.processes)
    proc = loaded.processes[0]
    assert proc.sigma(np.array([proc.energy_eV[-1] + 1.0]))[0] == pytest.approx(0.0)


def test_cross_section_high_energy_extrapolation_error_policy():
    cfg = load_config(ROOT / "configs" / "unified" / "boltzmann_only.yaml")
    cfg.cross_sections.high_energy_extrapolation = "error"
    loaded = load_cross_sections(cfg.cross_sections, cfg.conditions)
    proc = loaded.processes[0]
    with pytest.raises(ValueError, match="above the tabulated range"):
        proc.sigma(np.array([proc.energy_eV[-1] + 1.0]))


def test_lxcat_mass_ratio_rounding_does_not_warn():
    cases = [
        (ROOT / "cross_sections" / "Ar_Biagi.txt", "Ar"),
        (ROOT / "cross_sections" / "N2_Biagi.txt", "N2"),
        (ROOT / "cross_sections" / "Cl.txt", "Cl2"),
    ]

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        for path, species in cases:
            InterpolatedCrossSectionSet(50.0, str(path), species)

    unexpected = [
        warning
        for warning in caught
        if "Incorrect mass ratio." in str(warning.message)
    ]
    assert not unexpected


def test_mc_adapter_python_api(integrated_outputs):
    _, result, _ = integrated_outputs
    mc_cases = [case for case in result.cases if case.solver == "monte_carlo"]
    assert len(mc_cases) == 1
    case = mc_cases[0]
    assert np.isfinite(case.mean_energy_eV) and case.mean_energy_eV > 0.0
    assert np.isfinite(case.drift_velocity_m_s) and case.drift_velocity_m_s > 0.0
    assert np.isfinite(case.reduced_mobility_m2_V_s_m3)
    assert np.isfinite(case.diffusion_L_m2_s)
    assert np.isfinite(case.effective_townsend_m2)
    assert case.metadata["net_ionization_frequency_model"] == "convolution_effective_rate"
    assert np.isclose(
        case.net_ionization_frequency_s,
        case.metadata["convolution_net_ionization_frequency_s-1"],
        rtol=1e-12,
        atol=0.0,
    )


def test_command_adapter_infers_reduced_transport_from_density():
    cfg = load_config(CONFIG_PATH)
    cross_sections = load_cross_sections(cfg.cross_sections, cfg.conditions)
    adapter = MonteCarloAdapter(cfg, cross_sections)
    gas_density = cfg.conditions.pressure_Pa / (
        1.380649e-23 * cfg.conditions.gas_temperature_K
    )
    summary = pd.DataFrame(
        [
            {
                "run_label": "synthetic_0000",
                "E/N (Td)": 100.0,
                "mean energy (eV)": 2.0,
                "flux drift velocity (m.s-1)": 6.0,
                "mobility_m2_V_s": 2.0,
                "diffusion_L_m2_s": 3.0,
                "diffusion_T_m2_s": 4.0,
                "effective ionization rate coeff. (convolution) (m3.s-1)": 5.0,
                "meta_gas_number_density_m-3": gas_density,
            }
        ]
    )
    [case] = adapter._from_dataframes(summary, None)
    assert np.isclose(case.reduced_mobility_m2_V_s_m3, 2.0 * gas_density)
    assert np.isclose(case.reduced_diffusion_L_m2_s_m3, 3.0 * gas_density)
    assert np.isclose(case.reduced_diffusion_T_m2_s_m3, 4.0 * gas_density)
    assert np.isclose(case.net_ionization_frequency_s, 5.0 * gas_density)
    assert np.isclose(case.effective_townsend_m2, 5.0 / 6.0)


def test_monte_carlo_preserves_velocity_lut_selection():
    cfg = load_config(CONFIG_PATH)
    cross_sections = load_cross_sections(cfg.cross_sections, cfg.conditions)
    mc_cfg = _build_mc_config(
        cfg,
        cross_sections,
        case_id="lut_0000",
        e_over_n_Td=100.0,
        base_config=cfg.monte_carlo.passthrough["base_config"],
    )
    mc_cfg["simulation_settings"]["use_velocity_lut"] = True
    mc_cfg["simulation_settings"]["use_jit"] = False
    monte_carlo = MonteCarlo(McConfig(mc_cfg))
    assert monte_carlo._velocity_from_energy.__func__ is MonteCarlo._lut_velocity


def test_both_mode_dual_outputs(integrated_outputs):
    _, result, out_dir = integrated_outputs
    assert {case.solver for case in result.cases} == {
        "monte_carlo",
        "boltzmann_two_term",
    }

    summary = pd.read_csv(out_dir / "argon_summary.csv")
    assert {"monte_carlo", "boltzmann_two_term"} <= set(summary["solver"])

    summary_alias = pd.read_csv(out_dir / "summary.csv")
    assert summary_alias["solver"].eq("monte_carlo").all()
    assert (out_dir / "summary_mc.csv").exists()
    assert (out_dir / "summary_boltzmann.csv").exists()
    assert (out_dir / "eedf_table.csv").exists()
    assert (out_dir / "eedf_table_mc.csv").exists()
    assert (out_dir / "eedf_table_boltzmann.csv").exists()
    assert (out_dir / "energy_table.csv").exists()


def test_comsol_export_compat(integrated_outputs):
    _, _, out_dir = integrated_outputs
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
        run_dir=out_dir,
    )
    assert out is not None
    assert (out / "comsol_transport.csv").exists()
    assert (out / "comsol_rates.csv").exists()
    assert (out / "comsol_eedf.csv").exists()


def test_eedf_optimize_compat(integrated_outputs, tmp_path: Path):
    _, _, out_dir = integrated_outputs
    cfg = load_optimize_config(ROOT / "EEDF_optimize.yaml")
    results = process_table(
        out_dir / "eedf_table.csv",
        save_root=tmp_path,
        config=cfg,
        solver_filter="monte_carlo",
    )
    assert results
    fit_dir = tmp_path / out_dir.name
    assert (fit_dir / "fit_results.csv").exists()
