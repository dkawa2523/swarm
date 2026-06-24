from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from electron_swarm import load_config, run
from electron_swarm.core.cross_sections import (
    CrossSectionProcess,
    CrossSectionSet,
    ProcessType,
)
from electron_swarm.core.results import RateResult, SwarmCaseResult
from electron_swarm.diagnostics.tail import (
    attach_tail_metrics,
    rate_tail_fraction,
    resolve_tail_threshold_eV,
    tail_probability,
)

from product_helpers import base_product_config, write_config


def _constant_process(process_type: ProcessType = ProcessType.IONIZATION) -> CrossSectionProcess:
    return CrossSectionProcess(
        species="Ar",
        process=process_type.value,
        process_type=process_type,
        threshold_eV=2.0 if process_type != ProcessType.MOMENTUM else None,
        energy_eV=np.array([0.0, 1.0, 2.0, 3.0]),
        cross_section_m2=np.ones(4),
    )


def test_tail_probability_calculation() -> None:
    energy = np.array([0.0, 1.0, 2.0, 3.0])
    eedf = np.array([0.1, 0.2, 0.3, 0.4])
    widths = np.ones_like(energy)
    assert tail_probability(energy, eedf, 2.0, widths) == pytest.approx(0.7)


def test_rate_tail_fraction_calculation() -> None:
    energy = np.array([0.0, 1.0, 2.0, 3.0])
    eedf = np.ones_like(energy)
    process = _constant_process()
    fraction = rate_tail_fraction(energy, eedf, process, 2.0, np.ones_like(energy))
    expected = (np.sqrt(2.0) + np.sqrt(3.0)) / (
        np.sqrt(1.0) + np.sqrt(2.0) + np.sqrt(3.0)
    )
    assert fraction == pytest.approx(expected)


def test_tail_threshold_fallback_without_reaction_threshold(tmp_path: Path) -> None:
    data = base_product_config(tmp_path)
    data["physics"]["energy_grid_policy"]["tail_threshold_eV"] = 12.0
    cfg = load_config(write_config(tmp_path, data))
    xs = CrossSectionSet([_constant_process(ProcessType.MOMENTUM)])
    assert resolve_tail_threshold_eV(cfg, xs) == pytest.approx(12.0)

    data["physics"]["energy_grid_policy"]["tail_threshold_eV"] = None
    cfg = load_config(write_config(tmp_path, data, name="default_tail.yaml"))
    assert resolve_tail_threshold_eV(cfg, xs) == pytest.approx(20.0)


def test_tail_status_insufficient_for_cutoff_rate_contribution(tmp_path: Path) -> None:
    data = base_product_config(tmp_path)
    data["physics"]["energy_grid_policy"]["tail_threshold_eV"] = 2.0
    data["physics"]["energy_grid_policy"]["tail_rate_warning_fraction"] = 0.05
    cfg = load_config(write_config(tmp_path, data))
    process = CrossSectionProcess(
        species="Ar",
        process="ionization",
        process_type=ProcessType.IONIZATION,
        threshold_eV=2.0,
        energy_eV=np.array([1.0, 2.0, 9.0, 10.0]),
        cross_section_m2=np.array([0.0, 0.0, 1.0, 1.0]),
    )
    case = SwarmCaseResult(
        solver="two_term",
        case_id="tail",
        e_over_n_Td=100.0,
        mean_energy_eV=1.0,
        drift_velocity_m_s=1.0,
        mobility_m2_V_s=1.0,
        reduced_mobility_m2_V_s_m3=1.0,
        diffusion_L_m2_s=1.0,
        diffusion_T_m2_s=1.0,
        reduced_diffusion_L_m2_s_m3=1.0,
        reduced_diffusion_T_m2_s_m3=1.0,
        net_ionization_frequency_s=1.0,
        effective_townsend_m2=1.0,
        energy_eV=np.array([1.0, 2.0, 9.0, 10.0]),
        eedf=np.array([0.0, 0.0, 0.0, 1.0]),
        eepf=np.array([0.0, 0.0, 0.0, 1.0]),
        rates=[
            RateResult(
                solver="two_term",
                case_id="tail",
                e_over_n_Td=100.0,
                species="Ar",
                process="ionization",
                process_type="ionization",
                threshold_eV=2.0,
                rate_coefficient_m3_s=1.0,
                mixture_weighted_rate_m3_s=1.0,
            )
        ],
    )
    attach_tail_metrics(case, cfg, CrossSectionSet([process]))

    tail = case.diagnostics["tail_metrics"]
    assert tail["energy_grid_tail_status"] == "insufficient"
    assert tail["high_energy_cutoff_rate_fraction"] > 0.05
    assert case.rates[0].tail_fraction == pytest.approx(1.0)


def test_tail_metrics_disabled_keeps_output_shape(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["two_term"])
    data["physics"]["energy_grid_policy"]["tail_metrics"] = False
    cfg = load_config(write_config(tmp_path, data))
    result = run(cfg, write=True)
    [case] = result.cases
    assert case.diagnostics["tail_metrics"]["tail_metrics_enabled"] is False
    assert all(rate.tail_fraction is None for rate in case.rates)
    summary = pd.read_csv(tmp_path / "prod_summary.csv")
    rates = pd.read_csv(tmp_path / "prod_rates.csv")
    assert "meta_tail_refinement_treatment" in summary.columns
    assert "meta_tail_rate_fraction_max" not in summary.columns
    assert "meta_tail_probability" not in summary.columns
    assert summary["meta_tail_refinement_treatment"].iloc[0] == "approximate"
    assert "tail_fraction" in rates.columns
