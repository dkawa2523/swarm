"""Execute runnable solver-plan rows and attach product metadata."""

from __future__ import annotations

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.core.solver_configs import InternalSolverConfigs
from electron_swarm.core.solver_registry import build_solver
from electron_swarm.orchestration.metadata import attach_product_metadata
from electron_swarm.orchestration.plan import SolverPlanItem


def execute_solve_plan(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    plan: list[SolverPlanItem],
    internal: InternalSolverConfigs,
) -> list[SwarmCaseResult]:
    cases: list[SwarmCaseResult] = []
    for item in plan:
        if not item.runnable:
            continue
        solver = build_solver(config, cross_sections, internal, item.solver)
        solver_cases = solver.solve_all()
        for case in solver_cases:
            cases.append(attach_product_metadata(case, config=config, item=item))
    return cases
