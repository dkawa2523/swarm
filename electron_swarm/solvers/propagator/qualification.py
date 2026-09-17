"""P1 scope and input qualification."""

from __future__ import annotations

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.cross_sections import (
    CrossSectionSet,
    ProcessType,
    mixture_fraction,
)
from electron_swarm.core.scattering import resolve_particle_scattering


def qualify_p1_inputs(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
) -> dict[str, object]:
    if config.physics.field.type != "dc":
        raise NotImplementedError("propagator P1 supports homogeneous DC fields only")
    magnetic = config.physics.field.magnetic_field
    if magnetic.enabled and magnetic.B_T != 0.0:
        raise NotImplementedError("propagator P1 requires B=0")
    if config.physics.electron_electron.enabled:
        raise NotImplementedError("propagator P1 does not implement e-e collisions")
    if config.physics.finite_k.enabled:
        raise NotImplementedError("propagator P1 does not implement finite-k swarms")
    angular = config.physics.angular_scattering.model
    if angular not in {"isotropic", "maxent_p1"}:
        raise NotImplementedError(
            f"propagator P1 does not implement angular model {angular!r}"
        )
    unsupported = sorted(
        {
            process.process_type.value
            for process in cross_sections.processes
            if mixture_fraction(config.conditions, process.species) > 0.0
            if process.process_type == ProcessType.SUPERELASTIC
        }
    )
    if unsupported:
        raise NotImplementedError(
            "propagator P1 does not implement active process types "
            f"{unsupported}"
        )
    channels = resolve_particle_scattering(
        cross_sections,
        angular_model=angular,
        active_species=(
            component.species
            for component in config.conditions.gas_mixture
            if component.fraction > 0.0
        ),
    )
    return {
        "field": "homogeneous_dc",
        "magnetic_field": "zero",
        "angular_model": angular,
        "scattering_species": [channel.species for channel in channels],
        "elastic_collision_xs_role": sorted(
            {channel.collision_total_model for channel in channels}
        ),
        "inelastic_angular_model": "isotropic_integral_xs_closure",
    }
