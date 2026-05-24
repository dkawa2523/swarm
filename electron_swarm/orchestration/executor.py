"""Execute runnable solver-plan rows and attach product metadata."""

from __future__ import annotations

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.orchestration.plan import SolverPlanItem
from electron_swarm.physics.angular_scattering import (
    ANGULAR_METADATA_KEYS,
    expected_angular_metadata,
)
from electron_swarm.solvers.two_term import TwoTermSolver
from electron_swarm.solvers.monte_carlo_adapter import MonteCarloAdapter
from electron_swarm.solvers.multi_term import MultiTermSolver


SOLVER_FACTORIES = {
    "two_term": TwoTermSolver,
    "multi_term": MultiTermSolver,
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


def _transport_definition(case: SwarmCaseResult, solver: str) -> str:
    if solver == "monte_carlo":
        if case.metadata.get("monte_carlo_backend") == "internal":
            return "mc_particle_tracking"
        return "external_adapter"
    if case.transport is not None and case.transport.metadata is not None:
        return str(case.transport.metadata.coefficient_definition)
    return "flux"


def _attach_product_metadata(
    case: SwarmCaseResult,
    *,
    config: SwarmConfig,
    item: SolverPlanItem,
) -> SwarmCaseResult:
    solver_ee_metadata = {
        key: value
        for key, value in case.metadata.items()
        if key.startswith("electron_electron_")
    }
    solver_magnetic_metadata = {
        key: value
        for key, value in case.metadata.items()
        if key.startswith("magnetic_field_") or key == "field_integrator"
    }
    case.solver = item.solver
    case.schema_version = "2"
    for rate in case.rates:
        rate.solver = item.solver

    expected_angular = expected_angular_metadata(config)
    if item.solver == "monte_carlo":
        angular_metadata = {
            key: case.metadata.get(key)
            for key in ANGULAR_METADATA_KEYS
            if key in case.metadata
        }
        if config.solvers.monte_carlo.angular_scattering == "same_as_physics":
            angular_metadata = expected_angular
    else:
        angular_metadata = expected_angular
    physics_level = "surrogate" if item.solver == "multi_term" else "approximate"
    metadata = {
        "schema_version": config.schema_version,
        "solver_id": item.solver,
        "solver_method": _solver_method(config, item.solver),
        "solver_label": item.solver,
        "physics_level": physics_level,
        **angular_metadata,
        "direct_pn_operator": False,
        "electron_electron_treatment": item.effective_physics.get(
            "electron_electron", "none"
        ),
        "electron_electron_affects_eedf": False,
        "electron_electron_affects_rates": False,
        "electron_electron_affects_transport": False,
        "electron_electron_transport_stale": False,
        "ionization_source_model": config.physics.ionization.energy_sharing,
        "ionization_source_treatment": item.effective_physics.get(
            "ionization_source", "none"
        ),
        "ionization_secondary_electron_energy_eV": (
            config.physics.ionization.secondary_electron_energy_eV
        ),
        "transport_definition": _transport_definition(case, item.solver),
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
        method = config.solvers.multi_term.method
        is_pn_dcs = method == "pn_dcs"
        table = config.physics.angular_scattering.moment_table
        exact_dcs_based = bool(
            is_pn_dcs and table is not None and table.provenance == "dcs_derived"
        )
        metadata.update(
            {
                "lmax": config.solvers.multi_term.lmax,
                "direct_pn_operator": False,
                "angular_moment_source": str(angular_metadata["angular_moment_source"]),
                "exact_dcs_based": exact_dcs_based,
                "ordinary_integral_xs_closure": not is_pn_dcs,
                "physics_level": "table_moments" if is_pn_dcs else "surrogate",
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
    case.metadata.update(solver_ee_metadata)
    case.metadata.update(solver_magnetic_metadata)
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
