"""Shared, explicit gas conditions used by every numerical method."""

from dataclasses import dataclass, field

from .constants import BOLTZMANN_J_K


@dataclass(slots=True)
class GasComponent:
    species: str
    fraction: float
    mass_amu: float


@dataclass(slots=True)
class ConditionsConfig:
    gas_temperature_K: float = 300.0
    pressure_Pa: float | None = None
    gas_number_density_m3: float | None = None
    gas_mixture: list[GasComponent] = field(default_factory=list)


def gas_number_density(config: object) -> float:
    conditions = config.conditions
    if conditions.gas_number_density_m3 is not None:
        return float(conditions.gas_number_density_m3)
    if conditions.pressure_Pa is None:
        raise ValueError("gas pressure or number density is required")
    return float(
        conditions.pressure_Pa / (BOLTZMANN_J_K * conditions.gas_temperature_K)
    )
