"""Product-visible metadata assembly for executed solver-plan cases."""

from __future__ import annotations

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.result_metadata import (
    TRANSPORT_FLUX,
    TRANSPORT_F0_GRADIENT_RECONSTRUCTION,
    TRANSPORT_MC_FIXED_POPULATION,
)
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.core.solver_registry import solver_method
from electron_swarm.orchestration.plan import SolverPlanItem
from electron_swarm.physics.angular_scattering import expected_angular_metadata


def _transport_definition(
    case: SwarmCaseResult,
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


def _multi_term_metadata(
    case: SwarmCaseResult,
    config: SwarmConfig,
    angular_metadata: dict[str, object],
) -> dict[str, object]:
    method = config.solvers.multi_term.method
    is_pn_dcs = method == "pn_dcs"
    is_pn_direct = method == "pn_closure_direct"
    table = config.physics.angular_scattering.moment_table
    exact_dcs_based = bool(
        is_pn_dcs and table is not None and table.provenance == "dcs_derived"
    )
    return {
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
        "transport_definition": TRANSPORT_F0_GRADIENT_RECONSTRUCTION,
    }


def attach_product_metadata(
    case: SwarmCaseResult,
    *,
    config: SwarmConfig,
    item: SolverPlanItem,
) -> SwarmCaseResult:
    """Attach the compact product metadata contract to one solver result."""

    case.solver = item.solver
    case.schema_version = "2"
    for rate in case.rates:
        rate.solver = item.solver

    angular_metadata = expected_angular_metadata(config)
    metadata: dict[str, object] = {
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
        "transport_definition": _transport_definition(case, item.solver),
        "tail_refinement_treatment": item.treatment("tail_refinement"),
        "magnetic_field_treatment": str(
            case.metadata.get(
                "magnetic_field_treatment",
                item.treatment("magnetic_field"),
            )
        ),
    }
    if item.solver == "multi_term":
        metadata.update(_multi_term_metadata(case, config, angular_metadata))

    case.metadata.update(metadata)
    return case
