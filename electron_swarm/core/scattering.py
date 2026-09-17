"""Consumer-specific resolution of integral scattering cross sections."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal

from electron_swarm.core.cross_sections import (
    CrossSectionProcess,
    CrossSectionSet,
    ScatteringRole,
)


AngularCollisionModel = Literal["isotropic", "maxent_p1"]


@dataclass(frozen=True, slots=True)
class SpeciesScatteringChannels:
    """Resolved scattering inputs for one active gas species."""

    species: str
    collision_total: tuple[CrossSectionProcess, ...]
    momentum_transfer: tuple[CrossSectionProcess, ...]
    effective_momentum: tuple[CrossSectionProcess, ...]
    collision_total_model: str


def _role_processes(
    cross_sections: CrossSectionSet,
    species: str,
    role: ScatteringRole,
) -> tuple[CrossSectionProcess, ...]:
    return tuple(
        process
        for process in cross_sections.by_species(species)
        if process.scattering_role == role
    )


def resolve_species_scattering(
    cross_sections: CrossSectionSet,
    species: str,
    *,
    angular_model: AngularCollisionModel,
) -> SpeciesScatteringChannels:
    """Resolve collision and moment integrals without double counting.

    Particle methods consume a physical total collision cross section.
    Isotropic scattering has the exact identity sigma_total=sigma_m and may
    therefore use a momentum-transfer input as the total. The maxent-P1
    closure needs two independent integrals and never consumes an effective
    momentum cross section as a collision event.
    """

    totals = _role_processes(
        cross_sections,
        species,
        ScatteringRole.ELASTIC_TOTAL,
    )
    momentum = _role_processes(
        cross_sections,
        species,
        ScatteringRole.ELASTIC_MOMENTUM_TRANSFER,
    )
    effective = _role_processes(
        cross_sections,
        species,
        ScatteringRole.EFFECTIVE_MOMENTUM_TRANSFER,
    )
    if angular_model == "maxent_p1":
        missing: list[str] = []
        if not totals:
            missing.append("elastic_total")
        if not momentum:
            missing.append("elastic_momentum_transfer")
        if missing:
            raise NotImplementedError(
                f"{species}: maxent_p1 requires both explicit elastic_total "
                "and elastic_momentum_transfer cross sections; missing "
                + ", ".join(missing)
            )
        total_model = "explicit_total"
        collision_total = totals
    elif angular_model == "isotropic":
        if totals:
            collision_total = totals
            total_model = "explicit_total"
        elif momentum:
            collision_total = momentum
            total_model = "momentum_as_total_isotropic"
        elif effective:
            raise NotImplementedError(
                f"{species}: effective momentum cross sections do not define "
                "a particle collision frequency"
            )
        else:
            raise ValueError(f"{species}: no elastic scattering cross section")
    else:
        raise NotImplementedError(
            f"particle scattering model {angular_model!r} is not implemented"
        )
    return SpeciesScatteringChannels(
        species=species,
        collision_total=collision_total,
        momentum_transfer=momentum,
        effective_momentum=effective,
        collision_total_model=total_model,
    )


def resolve_particle_scattering(
    cross_sections: CrossSectionSet,
    *,
    angular_model: AngularCollisionModel,
    active_species: Iterable[str] | None = None,
) -> tuple[SpeciesScatteringChannels, ...]:
    """Resolve all species that contain at least one scattering integral."""

    roles = {
        ScatteringRole.ELASTIC_TOTAL,
        ScatteringRole.ELASTIC_MOMENTUM_TRANSFER,
        ScatteringRole.EFFECTIVE_MOMENTUM_TRANSFER,
    }
    species = {
        process.species
        for process in cross_sections.processes
        if process.scattering_role in roles
    }
    if active_species is not None:
        requested = {str(item) for item in active_species}
        missing = sorted(requested - species)
        if missing:
            raise ValueError(
                "Active gas species lack an elastic scattering cross section: "
                + ", ".join(missing)
            )
        species &= requested
    ordered_species = sorted(species)
    if not ordered_species:
        raise ValueError("At least one scattering cross section is required")
    return tuple(
        resolve_species_scattering(
            cross_sections,
            item,
            angular_model=angular_model,
        )
        for item in ordered_species
    )


def transport_momentum_processes(
    cross_sections: CrossSectionSet,
    species: str,
) -> tuple[tuple[CrossSectionProcess, ...], bool, str]:
    """Return the elastic contribution to a Boltzmann momentum frequency.

    The boolean indicates that an effective input already includes inelastic
    momentum loss and therefore suppresses separate inelastic additions.
    """

    effective = _role_processes(
        cross_sections,
        species,
        ScatteringRole.EFFECTIVE_MOMENTUM_TRANSFER,
    )
    if effective:
        return effective, True, "effective_momentum_includes_inelastic"
    momentum = _role_processes(
        cross_sections,
        species,
        ScatteringRole.ELASTIC_MOMENTUM_TRANSFER,
    )
    if momentum:
        return momentum, False, "explicit_elastic_momentum_transfer"
    totals = _role_processes(
        cross_sections,
        species,
        ScatteringRole.ELASTIC_TOTAL,
    )
    if totals:
        return totals, False, "elastic_total_isotropic_momentum_closure"
    return (), False, "missing"
