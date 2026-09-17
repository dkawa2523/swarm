"""Common execution contracts for product solvers."""

from __future__ import annotations

from abc import ABC, abstractmethod

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.cross_sections import (
    ActiveMixtureInputs,
    prepare_active_mixture_inputs,
)
from electron_swarm.core.results import SwarmCaseResult


class SwarmSolver(ABC):
    """A solver that executes the configured E/N batch.

    Batch execution is the product boundary used by orchestration.  Solvers
    whose cases are independent can inherit :class:`IndependentCaseSolver`
    instead of reimplementing the small batch loop.
    """

    name: str

    def __init__(
        self,
        config: SwarmConfig,
        cross_sections: ActiveMixtureInputs,
    ) -> None:
        self.config = config
        self.cross_sections: ActiveMixtureInputs = prepare_active_mixture_inputs(
            cross_sections,
            config.conditions,
        )

    @abstractmethod
    def solve_all(self) -> list[SwarmCaseResult]:
        raise NotImplementedError


class IndependentCaseSolver(SwarmSolver):
    """Adapter for solvers that can evaluate each configured case alone."""

    @abstractmethod
    def solve_case(self, e_over_n_Td: float, case_id: str) -> SwarmCaseResult:
        raise NotImplementedError

    def solve_all(self) -> list[SwarmCaseResult]:
        return [
            self.solve_case(e_over_n_Td, f"{self.config.run.case_prefix}_{i:04d}")
            for i, e_over_n_Td in enumerate(self.config.run.e_over_n_Td)
        ]
