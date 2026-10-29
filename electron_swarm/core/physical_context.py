"""Resolved physical problem identity, independent of numerical solver choice."""

from dataclasses import asdict, dataclass
from hashlib import sha256
from typing import TYPE_CHECKING

from .gas import gas_number_density
from .identity import content_digest

if TYPE_CHECKING:
    from .config import SwarmConfig

PHYSICAL_CONTEXT_SCHEMA = "swarm_physical_context.v2"


@dataclass(frozen=True, slots=True)
class AngularPhysics:
    model: str
    higher_moment_closure: str
    moment_table_sha256: str | None
    moment_table_provenance: str | None
    moment_table_extrapolation: str | None


@dataclass(frozen=True, slots=True)
class ElectronElectronPhysics:
    enabled: bool
    model: str
    strength_model: str | None
    relaxation_fraction: float | None
    conserve_mean_energy: bool | None
    target_temperature_eV: float | None


def angular_physics(config: "SwarmConfig") -> AngularPhysics:
    angular = config.physics.angular_scattering
    table = angular.moment_table
    return AngularPhysics(
        model=angular.model,
        higher_moment_closure=angular.higher_moment_closure,
        moment_table_sha256=sha256(table.path.read_bytes()).hexdigest()
        if table
        else None,
        moment_table_provenance=table.provenance if table else None,
        moment_table_extrapolation=table.extrapolation if table else None,
    )


def electron_electron_physics(config: "SwarmConfig") -> ElectronElectronPhysics:
    ee = config.physics.electron_electron
    return ElectronElectronPhysics(
        enabled=ee.enabled,
        model=ee.model if ee.enabled else "none",
        strength_model=ee.strength_model if ee.enabled else None,
        relaxation_fraction=ee.relaxation_fraction if ee.enabled else None,
        conserve_mean_energy=ee.conserve_mean_energy if ee.enabled else None,
        target_temperature_eV=ee.fallback_temperature_eV
        if ee.enabled and not ee.conserve_mean_energy
        else None,
    )


def _field_context(config: "SwarmConfig") -> dict:
    field = asdict(config.physics.field)
    magnetic = field["magnetic_field"]
    if not magnetic["enabled"] or magnetic["B_T"] == 0:
        field["magnetic_field"] = {"enabled": False, "B_T": 0.0, "angle_EB_deg": 0.0}
    if field["type"] == "dc":
        field.pop("time_dependent")
    else:
        for name in ("phase_steps", "max_periods", "periodic_tolerance"):
            field["time_dependent"].pop(name)
    return field


def physical_features(config: "SwarmConfig") -> dict:
    """One normalized physics view shared by direct and workflow results."""
    ionization = asdict(config.physics.ionization)
    if ionization["energy_sharing"] != "primary_secondary":
        ionization["secondary_electron_energy_eV"] = 0.0
    finite_k = config.physics.finite_k
    return {
        "field": _field_context(config),
        "angular_scattering": asdict(angular_physics(config)),
        "electron_electron": asdict(electron_electron_physics(config)),
        "ionization": ionization,
        "finite_k": {
            "enabled": finite_k.enabled,
            "k_m_inv": finite_k.k_m_inv if finite_k.enabled else None,
        },
    }


def workflow_physical_context(config: "SwarmConfig") -> dict:
    return {
        "schema": PHYSICAL_CONTEXT_SCHEMA,
        "gas_temperature_K": config.conditions.gas_temperature_K,
        "gas_number_density_m3": gas_number_density(config),
        **physical_features(config),
        "high_energy_extrapolation": config.cross_sections.high_energy_extrapolation,
    }


@dataclass(frozen=True, slots=True)
class PhysicalContext:
    gas_temperature_K: float
    gas_number_density_m3: float
    species: tuple[tuple[str, float, float], ...]
    collision_data_sha256: str
    cross_section_extrapolation: str
    field_type: str
    e_over_n_Td: float
    field_frequency_Hz: float | None
    magnetic_field_T: float
    angle_EB_deg: float
    growth_model: str
    ionization_energy_sharing: str
    secondary_electron_energy_eV: float
    angular_scattering: AngularPhysics
    electron_electron: ElectronElectronPhysics


def collision_data_identity(inputs) -> str:
    processes = sorted(
        inputs.processes, key=lambda p: (p.species, p.process, p.process_type.value)
    )
    return content_digest(
        [
            (
                p.species,
                p.process,
                p.process_type.value,
                p.threshold_eV,
                p.mass_amu,
                p.scattering_role.value if p.scattering_role is not None else None,
                p.energy_eV,
                p.cross_section_m2,
            )
            for p in processes
        ]
    )


def physical_context(
    config: "SwarmConfig", collision_identity: str, e_over_n_Td: float
) -> PhysicalContext:
    physics = config.physics
    field = physics.field
    magnetic = field.magnetic_field
    ionization = physics.ionization
    return PhysicalContext(
        gas_temperature_K=config.conditions.gas_temperature_K,
        gas_number_density_m3=gas_number_density(config),
        species=tuple(
            sorted(
                (gas.species, gas.fraction, gas.mass_amu)
                for gas in config.conditions.gas_mixture
                if gas.fraction > 0
            )
        ),
        collision_data_sha256=collision_identity,
        cross_section_extrapolation=config.cross_sections.high_energy_extrapolation,
        field_type=field.type,
        e_over_n_Td=float(e_over_n_Td),
        field_frequency_Hz=field.time_dependent.frequency_Hz
        if field.type == "time_dependent"
        else None,
        magnetic_field_T=magnetic.B_T if magnetic.enabled else 0.0,
        angle_EB_deg=magnetic.angle_EB_deg
        if magnetic.enabled and magnetic.B_T
        else 0.0,
        growth_model="temporal_growth"
        if ionization.population_model == "weighted_branching"
        else "fixed_population",
        ionization_energy_sharing=ionization.energy_sharing,
        secondary_electron_energy_eV=ionization.secondary_electron_energy_eV
        if ionization.energy_sharing == "primary_secondary"
        else 0.0,
        angular_scattering=angular_physics(config),
        electron_electron=electron_electron_physics(config),
    )
