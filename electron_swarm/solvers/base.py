"""Common solver protocol."""

from __future__ import annotations

from abc import ABC, abstractmethod

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.core.results import SwarmCaseResult


class SwarmSolver(ABC):
    name: str

    def __init__(self, config: SwarmConfig, cross_sections: CrossSectionSet) -> None:
        self.config = config
        self.cross_sections = cross_sections

    @abstractmethod
    def solve_case(self, e_over_n_Td: float, case_id: str) -> SwarmCaseResult:
        raise NotImplementedError

    def solve_all(self) -> list[SwarmCaseResult]:
        return [
            self.solve_case(e_over_n_Td, f"{self.config.run.case_prefix}_{i:04d}")
            for i, e_over_n_Td in enumerate(self.config.run.e_over_n_Td)
        ]
