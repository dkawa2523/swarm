from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from electron_swarm import load_config, run
from electron_swarm.core.solver_configs import build_internal_solver_configs
from electron_swarm.solvers.multi_term.case import build_multiterm_case
from electron_swarm.core.cross_sections import load_cross_sections, prepare_active_mixture_inputs

from product_helpers import base_product_config, write_config


@pytest.mark.parametrize("adaptive", [False, True])
def test_energy_ceiling_bounds_both_solver_cells_and_comparison_outputs(
    tmp_path: Path, adaptive: bool,
) -> None:
    data = base_product_config(tmp_path, ["two_term", "multi_term"])
    data["solvers"]["multi_term"]["lmax"] = 1
    data["physics"]["energy_grid_policy"].update(max_eV_limit=30.0, adaptive=adaptive)
    data["comparison"] = {
        "enabled": True, "required": True,
        "reference_solver": "two_term", "candidate_solvers": ["multi_term"],
    }
    result = run(load_config(write_config(tmp_path, data)), write=True)
    for case in result.cases:
        assert case.energy_eV[-1] <= 30.0
        assert np.sum(case.energy_widths_eV) <= 30.0 + 1.0e-12
        assert np.sum(case.eedf * case.energy_widths_eV) == pytest.approx(1.0)
    comparison = pd.read_csv(tmp_path / "prod_comparison_summary.csv")
    assert np.isfinite(comparison.loc[0, "eedf_l1_error"])
    pn = result.by_solver()["multi_term"][0]
    assert pn.diagnostics["multi_term"]["energy_grid_limit_eV"] == 30.0
    assert pn.diagnostics["multi_term"]["grid_refinement_cycles"] == 1
    assert pn.diagnostics["multi_term"]["tail_target_met"] == (
        pn.diagnostics["multi_term"]["tail_fraction"]
        <= data["physics"]["energy_grid_policy"]["tail_probability_target"]
    )


def test_multiterm_adaptive_request_re_solves_instead_of_only_labelling_it(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["multi_term"])
    data["solvers"]["multi_term"]["lmax"] = 1
    data["physics"]["energy_grid_policy"].update(
        max_eV_limit=350.0, adaptive=False,
    )
    config = load_config(write_config(tmp_path, data))
    fixed = run(config, write=False).cases[0]
    measured_tail = fixed.diagnostics["multi_term"]["tail_fraction"]
    assert measured_tail > 0.0
    config.physics.energy_grid_policy.adaptive = True
    config.physics.energy_grid_policy.tail_probability_target = measured_tail * 0.5
    refined = run(config, write=False).cases[0]
    detail = refined.diagnostics["multi_term"]
    assert detail["grid_refinement_cycles"] > 1
    assert refined.energy_eV[-1] > fixed.energy_eV[-1]
    assert np.sum(refined.energy_widths_eV) <= 350.0 + 1.0e-12
    assert detail["tail_target_met"] == (
        detail["tail_fraction"] <= config.physics.energy_grid_policy.tail_probability_target
    )


def test_multiterm_case_rejects_limit_below_its_minimum_energy(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["multi_term"])
    data["physics"]["energy_grid_policy"]["max_eV_limit"] = 1.0e-6
    config = load_config(write_config(tmp_path, data))
    inputs = prepare_active_mixture_inputs(
        load_cross_sections(config.cross_sections, config.conditions), config.conditions,
    )
    with pytest.raises(ValueError, match="max_eV must be larger than min_eV"):
        build_multiterm_case(
            config, inputs,
            build_internal_solver_configs(config.solvers, config.physics).multi_term,
            50.0,
        )
