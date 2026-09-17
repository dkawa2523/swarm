"""Public adapter for the internal Monte Carlo solver."""

from __future__ import annotations

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.cross_sections import ActiveMixtureInputs
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.core.solver_configs import MonteCarloAdapterConfig
from electron_swarm.solvers.base import SwarmSolver
import electron_swarm.solvers.monte_carlo.batch as _batch


class MonteCarloSolver(SwarmSolver):
    name = "monte_carlo"

    def __init__(
        self,
        config: SwarmConfig,
        cross_sections: ActiveMixtureInputs,
        solver_config: MonteCarloAdapterConfig,
    ) -> None:
        super().__init__(config, cross_sections)
        self.solver_config = solver_config

    def solve_all(self) -> list[SwarmCaseResult]:
        return _batch.run_internal_monte_carlo(
            self.config,
            self.cross_sections,
            self.solver_config,
            collect_audit=self.solver_config.collect_audit,
        )
