"""Product Monte Carlo solver entry point.

The product schema exposes Monte Carlo as the limited internal particle backend.
External MC/BOLSIG/MCIG execution and ingest live in benchmark tools.
"""

from __future__ import annotations

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.core.solver_configs import MonteCarloAdapterConfig

from .base import SwarmSolver
from .internal_monte_carlo import run_internal_monte_carlo


class MonteCarloAdapter(SwarmSolver):
    name = "monte_carlo"

    def __init__(
        self,
        config: SwarmConfig,
        cross_sections: CrossSectionSet,
        solver_config: MonteCarloAdapterConfig,
    ) -> None:
        super().__init__(config, cross_sections)
        self.solver_config = solver_config

    def solve_all(self) -> list[SwarmCaseResult]:
        return run_internal_monte_carlo(
            self.config,
            self.cross_sections,
            self.solver_config,
            collect_audit=self.solver_config.collect_audit,
        )

    def solve_case(self, e_over_n_Td: float, case_id: str) -> SwarmCaseResult:
        raise NotImplementedError("monte_carlo executes through solve_all()")
