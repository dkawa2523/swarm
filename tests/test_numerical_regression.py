from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from electron_swarm import load_config, run

from product_helpers import base_product_config, write_config


BASELINE = Path(__file__).parent / "regression_baselines" / "argon_two_term_multi_term.json"


def _assert_relative(actual: float, expected: float, tolerance: float) -> None:
    assert math.isfinite(actual)
    floor = max(abs(expected) * tolerance, 1.0e-40)
    assert abs(actual - expected) <= floor


@pytest.mark.regression
def test_argon_two_term_multi_term_scalar_regression(tmp_path: Path) -> None:
    baseline = json.loads(BASELINE.read_text(encoding="utf-8"))
    data = base_product_config(tmp_path, ["two_term", "multi_term"])
    data["comparison"] = {"enabled": False}
    cfg = load_config(write_config(tmp_path, data, "regression.yaml"))
    result = run(cfg, write=True)

    by_solver = {case.solver: case for case in result.cases}
    assert set(by_solver) == {"two_term", "multi_term"}
    for solver, expected_values in baseline["solvers"].items():
        case = by_solver[solver]
        assert case.mean_energy_eV > 0.0
        assert np.all(np.isfinite(case.eedf))
        assert np.isclose(np.trapezoid(case.eedf, case.energy_eV), 1.0, rtol=2.0e-2)
        transport_tol = baseline["tolerances"]["transport_relative"]
        rate_tol = baseline["tolerances"]["rate_relative"]
        for metric in ["mean_energy_eV", "drift_velocity_m_s", "mobility_m2_V_s"]:
            _assert_relative(float(getattr(case, metric)), expected_values[metric], transport_tol)
        for metric in ["net_ionization_frequency_s", "effective_townsend_m2"]:
            _assert_relative(float(getattr(case, metric)), expected_values[metric], rate_tol)

    summary = pd.read_csv(tmp_path / "prod_summary.csv")
    for column in [
        "solver",
        "solver_method",
        "mean_energy_eV",
        "drift_velocity_m_s",
        "mobility_m2_V_s",
        "meta_transport_definition",
        "meta_tail_refinement_treatment",
    ]:
        assert column in summary.columns
