"""Internal data models for multi-term Boltzmann backends."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.core.results import RateResult
from electron_swarm.core.transport import BulkTransport, TransportSet

from .diagnostics import SolverDiagnostics
from .grid import EnergyGrid


@dataclass(frozen=True, slots=True)
class MultiTermCase:
    config: SwarmConfig
    cross_sections: CrossSectionSet
    grid: EnergyGrid
    e_over_n_Td: float
    gas_number_density_m3: float
    electric_field_V_m: float


@dataclass(frozen=True, slots=True)
class RateSet:
    rates: tuple[RateResult, ...]
    ionization_frequency_s_inv: float
    attachment_frequency_s_inv: float

    @property
    def effective_growth_frequency_s_inv(self) -> float:
        return self.ionization_frequency_s_inv - self.attachment_frequency_s_inv


@dataclass(frozen=True, slots=True)
class MultiTermSolution:
    energy_eV: np.ndarray
    widths_eV: np.ndarray
    coefficients: np.ndarray
    eedf_eV_inv: np.ndarray
    rates: RateSet
    transport: TransportSet
    estimated_bulk: BulkTransport | None
    diagnostics: SolverDiagnostics
    method_used: str
    metadata: dict[str, object] = field(default_factory=dict)
    angular_moments: np.ndarray | None = None

    @property
    def mean_energy_eV(self) -> float:
        return float(np.sum(self.energy_eV * self.eedf_eV_inv * self.widths_eV))
