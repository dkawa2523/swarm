"""Execute runnable solver-plan rows and attach product metadata."""

from __future__ import annotations

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.cross_sections import ActiveMixtureInputs
from electron_swarm.core.results import SwarmCaseResult, validate_case_result
from electron_swarm.core.solver_configs import InternalSolverConfigs
from electron_swarm.core.solver_registry import build_solver
from electron_swarm.orchestration.metadata import attach_product_metadata
from electron_swarm.orchestration.plan import SolverPlanItem


def execute_solve_plan(
    config: SwarmConfig,
    cross_sections: ActiveMixtureInputs,
    plan: list[SolverPlanItem],
    internal: InternalSolverConfigs,
) -> list[SwarmCaseResult]:
    cases: list[SwarmCaseResult] = []
    for item in plan:
        if not item.runnable:
            continue
        solver = build_solver(config, cross_sections, internal, item.solver)
        solver_cases = solver.solve_all()
        expected = [
            (f"{config.run.case_prefix}_{index:04d}", float(e_over_n_Td))
            for index, e_over_n_Td in enumerate(config.run.e_over_n_Td)
        ]
        if len(solver_cases) != len(expected):
            raise ValueError(
                f"{item.solver} returned {len(solver_cases)} cases; "
                f"expected {len(expected)}"
            )
        for case, (expected_id, expected_e_over_n) in zip(
            solver_cases,
            expected,
            strict=True,
        ):
            if case.solver != item.solver:
                raise ValueError(
                    f"{item.solver} returned a result labelled {case.solver!r}"
                )
            if case.case_id != expected_id or not abs(
                float(case.e_over_n_Td) - expected_e_over_n
            ) <= 1.0e-12 * max(1.0, abs(expected_e_over_n)):
                raise ValueError(
                    f"{item.solver} returned an unexpected case identity"
                )
            validate_case_result(case)
            cases.append(attach_product_metadata(case, config=config, item=item))
    return cases
