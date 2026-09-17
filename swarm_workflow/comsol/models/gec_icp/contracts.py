"""Typed mapping and offline MPH contracts for the GEC ICP model."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


class GecIcpContractError(RuntimeError):
    """Raised when a GEC ICP mapping or local MPH violates its contract."""


@dataclass(frozen=True, slots=True)
class GecIcpModelSpec:
    input_mph: Path
    output_mph: Path
    component: str
    geometry: str
    mesh: str
    geometry_dimension: int
    axisymmetric: bool
    plasma_physics: str
    plasma_feature: str
    magnetic_physics: str
    coil_feature: str
    plasma_conductivity_coupling: str
    electron_heat_source_coupling: str
    study: str
    study_feature: str
    solution: str
    dataset: str


@dataclass(frozen=True, slots=True)
class GecIcpParameterTags:
    power: str
    gas_temperature: str
    pressure: str


@dataclass(frozen=True, slots=True)
class GecIcpSpeciesTags:
    ground: str
    excited: str
    ion: str


@dataclass(frozen=True, slots=True)
class GecIcpReactionTags:
    elastic: str
    excitation: str
    superelastic: str
    ionization: str
    stepwise_ionization: str


@dataclass(frozen=True, slots=True)
class GecIcpMapping:
    path: Path
    root: Path
    model: GecIcpModelSpec
    parameters: GecIcpParameterTags
    species: GecIcpSpeciesTags
    reactions: GecIcpReactionTags


@dataclass(frozen=True, slots=True)
class MphComponentContract:
    tag: str
    axisymmetric: bool


@dataclass(frozen=True, slots=True)
class MphGeometryContract:
    tag: str
    dimension: int


@dataclass(frozen=True, slots=True)
class MphPhysicsContract:
    tag: str
    operation: str
    geometry: str | None


@dataclass(frozen=True, slots=True)
class MphPhysicsFeatureContract:
    physics: str
    tag: str
    operation: str
    enabled: bool


@dataclass(frozen=True, slots=True)
class MphCouplingContract:
    tag: str
    operation: str
    component: str | None
    enabled: bool


@dataclass(frozen=True, slots=True)
class MphStudyFeatureContract:
    study: str
    tag: str
    operation: str
    enabled: bool


@dataclass(frozen=True, slots=True)
class MphSolutionContract:
    tag: str
    study: str | None


@dataclass(frozen=True, slots=True)
class MphDatasetContract:
    tag: str
    operation: str
    solution: str | None


@dataclass(frozen=True, slots=True)
class MphParameterContract:
    name: str
    expression: str


@dataclass(frozen=True, slots=True)
class MphSpeciesContract:
    physics: str
    tag: str
    species_type: str
    enabled: bool


@dataclass(frozen=True, slots=True)
class MphReactionContract:
    physics: str
    tag: str
    formula: str
    collision_type: str
    energy_loss_eV: float
    specification: str
    energy_data_eV: tuple[float, ...]
    cross_section_data_m2: tuple[float, ...]
    enabled: bool


@dataclass(frozen=True, slots=True)
class GecIcpMphContract:
    path: Path
    comsol_version: str | None
    model_title: str | None
    components: tuple[MphComponentContract, ...]
    geometries: tuple[MphGeometryContract, ...]
    meshes: tuple[str, ...]
    physics: tuple[MphPhysicsContract, ...]
    physics_features: tuple[MphPhysicsFeatureContract, ...]
    couplings: tuple[MphCouplingContract, ...]
    studies: tuple[str, ...]
    study_features: tuple[MphStudyFeatureContract, ...]
    solutions: tuple[MphSolutionContract, ...]
    datasets: tuple[MphDatasetContract, ...]
    parameters: tuple[MphParameterContract, ...]
    species: tuple[MphSpeciesContract, ...]
    reactions: tuple[MphReactionContract, ...]
