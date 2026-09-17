from __future__ import annotations

from typing import Any, cast

import numpy as np
import pytest

from electron_swarm import load_config
from electron_swarm.core.cross_sections import load_cross_sections
from electron_swarm.core.solver_configs import build_internal_solver_configs
from electron_swarm.solvers.two_term import TwoTermSolver
from electron_swarm.solvers.two_term.models import (
    NativeDistributionResult,
    NativeSolveDiagnostics,
    TimePeriodicDistributionResult,
)

from product_helpers import base_product_config, write_config


def _solver(tmp_path, data: dict[str, Any]) -> TwoTermSolver:
    config = load_config(write_config(tmp_path, data))
    cross_sections = load_cross_sections(config.cross_sections, config.conditions)
    internal = build_internal_solver_configs(config.solvers, config.physics)
    return TwoTermSolver(config, cross_sections, internal.two_term)


def test_stationary_two_term_refuses_last_unconverged_iterate(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    solver = _solver(tmp_path, base_product_config(tmp_path, ["two_term"]))
    energy = np.array([0.5, 1.5])
    result = NativeDistributionResult(
        energy_eV=energy,
        edges_eV=np.array([0.0, 1.0, 2.0]),
        widths_eV=np.ones(2),
        eedf_eV_inv=np.array([0.5, 0.5]),
        diagnostics=NativeSolveDiagnostics(
            converged=False,
            iterations=1,
            residual=1.0,
            residual_requested_tolerance=1.0e-8,
            residual_effective_tolerance=1.0e-8,
            residual_roundoff_bound=0.0,
            residual_backward_error=1.0,
            residual_roundoff_limited=False,
            growth_frequency_s=0.0,
            tail_probability=0.0,
            edge_to_peak=0.0,
            grid_max_eV=2.0,
            regrid_cycles=1,
        ),
        metadata={},
        operator_block=cast(Any, None),
    )
    monkeypatch.setattr(solver, "solve_native_distribution", lambda _field: result)

    with pytest.raises(RuntimeError, match="stationary distribution did not converge"):
        solver.solve_case(10.0, "case")


def test_periodic_two_term_refuses_unconverged_cycle(
    tmp_path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    data = base_product_config(tmp_path, ["two_term"])
    data["physics"]["field"] = {
        "type": "time_dependent",
        "magnetic_field": {
            "enabled": False,
            "B_T": 0.0,
            "angle_EB_deg": 0.0,
        },
        "time_dependent": {
            "waveform": "sinusoidal",
            "frequency_Hz": 13.56e6,
            "amplitude_definition": "rms",
            "momentum_response": "instantaneous",
            "phase_steps": 8,
            "max_periods": 4,
            "periodic_tolerance": 0.5,
        },
    }
    solver = _solver(tmp_path, data)
    result = TimePeriodicDistributionResult(
        energy_eV=np.array([0.5, 1.5]),
        widths_eV=np.ones(2),
        cycle_averaged_eedf_eV_inv=np.array([0.5, 0.5]),
        phase=cast(Any, None),
        diagnostics={"converged": False},
        operator_block=cast(Any, None),
    )
    monkeypatch.setattr(
        solver,
        "solve_time_periodic_distribution",
        lambda _field, *, case_id: result,
    )

    with pytest.raises(RuntimeError, match="time-periodic distribution did not converge"):
        solver.solve_case(10.0, "case")
