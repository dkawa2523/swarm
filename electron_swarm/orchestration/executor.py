"""Execute runnable solver-plan rows and attach product metadata."""

from __future__ import annotations

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.orchestration.plan import SolverPlanItem
from electron_swarm.physics.angular_scattering import build_angular_model
from electron_swarm.solvers.boltzmann_two_term import BoltzmannTwoTermSolver
from electron_swarm.solvers.monte_carlo_adapter import MonteCarloAdapter
from electron_swarm.solvers.multiterm_boltzmann import MultiTermBoltzmannSolver


SOLVER_FACTORIES = {
    "two_term": BoltzmannTwoTermSolver,
    "multi_term": MultiTermBoltzmannSolver,
    "monte_carlo": MonteCarloAdapter,
}


def _solver_method(config: SwarmConfig, solver: str) -> str:
    if solver == "two_term":
        return str(config.solvers.two_term.backend)
    if solver == "multi_term":
        return str(config.solvers.multi_term.method)
    if solver == "monte_carlo":
        if config.solvers.monte_carlo.python_api:
            return "python_api"
        if config.solvers.monte_carlo.command:
            return "command"
        return "adapter_unconfigured"
    return ""


def _attach_product_metadata(
    case: SwarmCaseResult,
    *,
    config: SwarmConfig,
    item: SolverPlanItem,
) -> SwarmCaseResult:
    case.solver = item.solver
    case.schema_version = "2"
    for rate in case.rates:
        rate.solver = item.solver

    angular_model = build_angular_model(config)
    angular_metadata = angular_model.metadata()
    physics_level = "surrogate" if item.solver == "multi_term" else "approximate"
    metadata = {
        "schema_version": config.schema_version,
        "solver_id": item.solver,
        "solver_method": _solver_method(config, item.solver),
        "solver_label": item.solver,
        "physics_level": physics_level,
        **angular_metadata,
        "exact_dcs_based": False,
        "ordinary_integral_xs_closure": "",
        "direct_pn_operator": False,
        "electron_electron_treatment": item.effective_physics.get(
            "electron_electron", "none"
        ),
        "electron_electron_affects_eedf": False,
        "electron_electron_affects_rates": False,
        "electron_electron_affects_transport": False,
        "electron_electron_transport_stale": False,
        "tail_refinement_treatment": item.effective_physics.get(
            "tail_refinement", "none"
        ),
        "feature_degraded": bool(item.degraded),
        "feature_warnings": "; ".join(item.warnings),
        "solve_plan_runnable": bool(item.runnable),
        "solve_plan_degraded": bool(item.degraded),
        **{f"effective_{key}": value for key, value in item.effective_physics.items()},
    }
    if item.solver == "multi_term":
        metadata.update(
            {
                "lmax": config.solvers.multi_term.lmax,
                "direct_pn_operator": False,
                "angular_moment_source": angular_metadata["angular_moment_source"],
                "exact_dcs_based": False,
                "ordinary_integral_xs_closure": True,
                "physics_level": "surrogate",
            }
        )
    magnetic = config.physics.field.magnetic_field
    metadata.update(
        {
            "magnetic_field_requested": bool(magnetic.enabled),
            "magnetic_field_B_T": float(magnetic.B_T),
            "magnetic_field_angle_EB_deg": float(magnetic.angle_EB_deg),
        }
    )
    metadata["magnetic_field_treatment"] = item.effective_physics.get(
        "magnetic_field", "none"
    )
    case.metadata.update(metadata)
    return case


def execute_solve_plan(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    plan: list[SolverPlanItem],
) -> list[SwarmCaseResult]:
    cases: list[SwarmCaseResult] = []
    for item in plan:
        if not item.runnable:
            continue
        factory = SOLVER_FACTORIES[item.solver]
        solver_cases = factory(config, cross_sections).solve_all()
        for case in solver_cases:
            cases.append(_attach_product_metadata(case, config=config, item=item))
    return cases
