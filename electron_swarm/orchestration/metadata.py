"""Product-visible metadata assembly for executed solver-plan cases."""

from __future__ import annotations

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.core.solver_registry import solver_method
from electron_swarm.orchestration.plan import SolverPlanItem
from electron_swarm.orchestration.solver_product_contracts import (
    solver_product_contract,
)
from electron_swarm.physics.angular_scattering import expected_angular_metadata


def attach_product_metadata(
    case: SwarmCaseResult,
    *,
    config: SwarmConfig,
    item: SolverPlanItem,
) -> SwarmCaseResult:
    """Attach the compact product metadata contract to one solver result."""

    angular_metadata = expected_angular_metadata(config)
    contract_metadata = solver_product_contract(item.solver).product_metadata(
        case,
        config=config,
        item=item,
        selected_angular_metadata=angular_metadata,
    )

    case.solver = item.solver
    case.schema_version = "2"
    for rate in case.rates:
        rate.solver = item.solver

    metadata: dict[str, object] = {
        "schema_version": config.schema_version,
        "solver_method": solver_method(config, item.solver),
        **angular_metadata,
        "angular_scattering_treatment": item.treatment("angular_scattering"),
        "angular_scattering_fidelity": item.fidelity("angular_scattering"),
        "angular_scattering_assumption": item.assumption("angular_scattering"),
        "direct_pn_operator": False,
        "electron_electron_treatment": item.treatment("electron_electron"),
        "electron_electron_transport_stale": False,
        "ionization_source_model": config.physics.ionization.energy_sharing,
        "ionization_source_treatment": item.treatment("ionization_source"),
        "ionization_secondary_electron_energy_eV": (
            config.physics.ionization.secondary_electron_energy_eV
        ),
        "transport_definition": contract_metadata["transport_definition"],
        "tail_refinement_treatment": item.treatment("tail_refinement"),
        "magnetic_field_treatment": item.treatment("magnetic_field"),
        "rf_field_treatment": item.treatment("rf_field"),
        "rf_frequency_Hz": (
            config.physics.field.time_dependent.frequency_Hz
            if config.physics.field.type == "time_dependent"
            else None
        ),
        "rf_amplitude_definition": (
            config.physics.field.time_dependent.amplitude_definition
            if config.physics.field.type == "time_dependent"
            else "none"
        ),
    }
    metadata.update(contract_metadata)

    case.metadata.update(metadata)
    return case
