"""Thin product entry point for the axisymmetric multi-term solver."""

from __future__ import annotations

from time import perf_counter

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.cross_sections import ActiveMixtureInputs
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.core.solver_configs import MultiTermInternalConfig
from electron_swarm.solvers.base import IndependentCaseSolver

from .refinement import solve_refined_case
from .result import to_case_result


MULTITERM_SOLVER_NAME = "multi_term"


class MultiTermSolver(IndependentCaseSolver):
    """Solve a homogeneous DC axisymmetric PN system on an energy grid."""

    name = MULTITERM_SOLVER_NAME

    def __init__(
        self,
        config: SwarmConfig,
        cross_sections: ActiveMixtureInputs,
        solver_config: MultiTermInternalConfig,
    ) -> None:
        super().__init__(config, cross_sections)
        self.solver_config = solver_config

    def solve_case(self, e_over_n_Td: float, case_id: str) -> SwarmCaseResult:
        started = perf_counter()
        if self.solver_config.method not in {"pn_closure_direct", "pn_dcs"}:
            raise NotImplementedError(
                f"Unsupported multi_term method {self.solver_config.method!r}"
            )
        case, solution = solve_refined_case(
            self.config,
            self.cross_sections,
            self.solver_config,
            e_over_n_Td,
            case_id,
            self.name,
        )
        return to_case_result(
            case,
            solution,
            case_id,
            perf_counter() - started,
            self.name,
        )


__all__ = ["MultiTermSolver"]
