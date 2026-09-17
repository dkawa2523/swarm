"""Execute one internal Monte Carlo field case with its derived RNG stream."""

from __future__ import annotations

from copy import deepcopy

import numpy as np

import electron_swarm.solvers.monte_carlo.case_phases as _phases
import electron_swarm.solvers.monte_carlo.case_result as _result
import electron_swarm.solvers.monte_carlo.case_state as _state
from electron_swarm.solvers.monte_carlo.field_parity import (
    apply_field_parity_response,
)
import electron_swarm.solvers.monte_carlo.setup as _setup
from electron_swarm.core.results import SwarmCaseResult


def _run_field_leg(
    setup: _setup._MonteCarloRunSetup,
    rng: np.random.Generator,
    *,
    e_over_n_Td: float,
    case_id: str,
    case_seed: int,
    field_polarity: int,
    run_tail: bool,
) -> SwarmCaseResult:
    state = _state._initialize_monte_carlo_case(
        setup,
        rng,
        e_over_n_Td=e_over_n_Td,
        case_id=case_id,
        case_seed=case_seed,
        field_polarity=field_polarity,
    )
    _phases._run_fixed_particle_phase(setup, state, rng)
    if setup.branching_active:
        if run_tail:
            _phases._run_weighted_phase(setup, state, rng)
        else:
            _phases._run_weighted_transport_phase(setup, state, rng)
    return _result._build_monte_carlo_case_result(setup, state)


def _run_monte_carlo_case(
    setup: _setup._MonteCarloRunSetup,
    rng: np.random.Generator,
    *,
    e_over_n_Td: float,
    case_id: str,
    case_seed: int,
) -> SwarmCaseResult:
    paired = setup.cfg.transport_estimator == "paired_field_parity"
    # Cloning does not consume the primary stream.  Consequently EEDF, rates,
    # and the primary case stream remain bitwise the same as single-field
    # execution for the same positive-field budget.
    mirror_rng = deepcopy(rng) if paired else None
    primary = _run_field_leg(
        setup,
        rng,
        e_over_n_Td=e_over_n_Td,
        case_id=case_id,
        case_seed=case_seed,
        field_polarity=1,
        run_tail=True,
    )
    if mirror_rng is None:
        return primary
    mirror = _run_field_leg(
        setup,
        mirror_rng,
        e_over_n_Td=e_over_n_Td,
        case_id=f"{case_id}__negative_field_leg",
        case_seed=case_seed,
        field_polarity=-1,
        run_tail=False,
    )
    return apply_field_parity_response(primary, mirror)
