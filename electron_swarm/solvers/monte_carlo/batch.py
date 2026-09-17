"""Order-independent batch adapter for internal Monte Carlo cases."""

from __future__ import annotations

import numpy as np

import electron_swarm.solvers.monte_carlo.case as _case
import electron_swarm.solvers.monte_carlo.setup as _setup
from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.cross_sections import ActiveMixtureInputs
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.core.solver_configs import MonteCarloAdapterConfig
from electron_swarm.solvers.monte_carlo.seeding import derive_case_seed


def run_internal_monte_carlo(
    config: SwarmConfig,
    cross_sections: ActiveMixtureInputs,
    solver_config: MonteCarloAdapterConfig,
    *,
    collect_audit: bool = False,
) -> list[SwarmCaseResult]:
    if solver_config.seed is None:
        raise ValueError(
            "solvers.monte_carlo.seed is required for a direct Monte Carlo run; "
            "workflow campaigns supply each replica seed from mc.base_seed"
        )
    if isinstance(solver_config.seed, bool) or not isinstance(solver_config.seed, int):
        raise ValueError("solvers.monte_carlo.seed must be a nonnegative integer")
    if solver_config.seed < 0:
        raise ValueError("solvers.monte_carlo.seed must be a nonnegative integer")
    setup = _setup._prepare_monte_carlo_run(
        config,
        cross_sections,
        solver_config,
        collect_audit=collect_audit,
    )
    results: list[SwarmCaseResult] = []
    for index, e_over_n_Td in enumerate(config.run.e_over_n_Td):
        case_seed = derive_case_seed(
            base_seed=solver_config.seed,
            e_over_n_Td=float(e_over_n_Td),
        )
        rng = np.random.default_rng(case_seed)
        results.append(
            _case._run_monte_carlo_case(
                setup,
                rng,
                e_over_n_Td=e_over_n_Td,
                case_id=f"{config.run.case_prefix}_{index:04d}",
                case_seed=case_seed,
            )
        )
    return results
