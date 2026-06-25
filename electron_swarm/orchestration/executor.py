"""Execute runnable solver-plan rows and attach product metadata."""

from __future__ import annotations

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.core.result_metadata import (
    TRANSPORT_FLUX,
    TRANSPORT_MC_FIXED_POPULATION,
)
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.core.solver_configs import InternalSolverConfigs
from electron_swarm.core.solver_registry import solver_method
from electron_swarm.orchestration.plan import SolverPlanItem
from electron_swarm.physics.angular_scattering import expected_angular_metadata
from electron_swarm.solvers.base import SwarmSolver
from electron_swarm.solvers.internal_monte_carlo import MonteCarloSolver
from electron_swarm.solvers.two_term import TwoTermSolver
from electron_swarm.solvers.multi_term import MultiTermSolver


def _build_solver(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    internal: InternalSolverConfigs,
    solver: str,
) -> SwarmSolver:
    match solver:
        case "two_term":
            return TwoTermSolver(config, cross_sections, internal.two_term)
        case "multi_term":
            return MultiTermSolver(
                config,
                cross_sections,
                internal.multi_term,
                internal.two_term,
            )
        case "monte_carlo":
            return MonteCarloSolver(config, cross_sections, internal.monte_carlo)
        case _:
            raise ValueError(f"Unsupported solver id: {solver}")


def _transport_definition(
    case: SwarmCaseResult,
    config: SwarmConfig,
    solver: str,
) -> str:
    if solver == "monte_carlo":
        return str(
            case.metadata.get(
                "transport_definition",
                TRANSPORT_MC_FIXED_POPULATION,
            )
        )
    if case.transport is not None:
        return case.transport.definition
    return TRANSPORT_FLUX


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

    expected_angular = expected_angular_metadata(config)
    angular_metadata = expected_angular
    metadata = {
        "schema_version": config.schema_version,
        "solver_method": solver_method(config, item.solver),
        **angular_metadata,
        "direct_pn_operator": False,
        "electron_electron_treatment": item.treatment("electron_electron"),
        "electron_electron_transport_stale": False,
        "ionization_source_model": config.physics.ionization.energy_sharing,
        "ionization_source_treatment": item.treatment("ionization_source"),
        "ionization_secondary_electron_energy_eV": (
            config.physics.ionization.secondary_electron_energy_eV
        ),
        "transport_definition": _transport_definition(case, config, item.solver),
        "tail_refinement_treatment": item.treatment("tail_refinement"),
    }
    if item.solver == "multi_term":
        method = config.solvers.multi_term.method
        is_pn_dcs = method == "pn_dcs"
        is_pn_direct = method == "pn_closure_direct"
        table = config.physics.angular_scattering.moment_table
        exact_dcs_based = bool(
            is_pn_dcs and table is not None and table.provenance == "dcs_derived"
        )
        metadata.update(
            {
                "lmax": config.solvers.multi_term.lmax,
                "direct_pn_operator": bool(
                    (is_pn_direct or is_pn_dcs)
                    and case.metadata.get("direct_pn_operator") is True
                ),
                "angular_moment_source": str(angular_metadata["angular_moment_source"]),
                "moment_table_provenance": (
                    str(table.provenance) if is_pn_dcs and table is not None else ""
                ),
                "exact_dcs_based": exact_dcs_based,
                "ordinary_integral_xs_closure": not is_pn_dcs,
            }
        )
    metadata["magnetic_field_treatment"] = str(
        case.metadata.get(
            "magnetic_field_treatment",
            item.treatment("magnetic_field"),
        )
    )
    case.metadata.update(metadata)
    return case


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
        solver = _build_solver(config, cross_sections, internal, item.solver)
        solver_cases = solver.solve_all()
        for case in solver_cases:
            cases.append(_attach_product_metadata(case, config=config, item=item))
    return cases
