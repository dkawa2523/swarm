"""Bounded energy-domain refinement for a complete stationary PN case."""

from __future__ import annotations

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.core.solver_configs import (
    AdaptiveGridConfig,
    MultiTermInternalConfig,
)

from .case import build_multiterm_case
from .models import MultiTermCase, MultiTermSolution
from .observables import build_solution
from .operator import assemble_pn_operator
from .steady import solve_stationary_pn


def solve_refined_case(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    solver_config: MultiTermInternalConfig,
    e_over_n_Td: float,
    case_id: str,
    solver_name: str,
) -> tuple[MultiTermCase, MultiTermSolution]:
    """Re-solve all PN moments when the requested energy-tail target requires it.

    Refinement is limited by both the product energy ceiling and a finite
    numerical cycle budget. A bound-limited answer retains its measured tail
    and target status; it never claims that the requested target was attained.
    """

    policy = config.physics.energy_grid_policy
    controls = AdaptiveGridConfig()
    maximum = min(float(solver_config.energy_grid.max_eV), policy.max_eV_limit)
    cycle_budget = controls.max_cycles if policy.adaptive else 1
    for cycle in range(cycle_budget):
        case = build_multiterm_case(
            config,
            cross_sections,
            solver_config,
            e_over_n_Td,
            max_eV_override=maximum,
        )
        operator = assemble_pn_operator(case)
        state = solve_stationary_pn(operator, solver_config.convergence)
        solution = build_solution(case, operator, state, case_id, solver_name)
        target_met = bool(
            solution.diagnostics.eedf_tail_fraction <= policy.tail_probability_target
        )
        solution.metadata.update(
            {
                "energy_grid_limit_eV": float(policy.max_eV_limit),
                "grid_refinement_cycles": cycle + 1,
                "tail_probability_target": float(policy.tail_probability_target),
                "tail_target_met": target_met,
            }
        )
        if not policy.adaptive or target_met:
            break
        next_maximum = min(
            policy.max_eV_limit,
            max(
                maximum * 1.25,
                controls.mean_energy_multiplier * solution.mean_energy_eV,
            ),
        )
        if next_maximum <= maximum * (1.0 + 1.0e-12):
            break
        maximum = next_maximum
    return case, solution
