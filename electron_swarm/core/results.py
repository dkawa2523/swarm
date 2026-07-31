"""Shared result models for Monte Carlo and Boltzmann solvers."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from electron_swarm.core.numerics import eepf_from_eedf
from electron_swarm.core.transport import ElectronTransport


@dataclass(slots=True)
class RateResult:
    """Rate coefficient for one collision/reaction process."""

    solver: str
    case_id: str
    e_over_n_Td: float
    species: str
    process: str
    process_type: str
    threshold_eV: float | None
    rate_coefficient_m3_s: float
    target_species_fraction: float = 1.0
    energy_loss_eV: float | None = None
    gas_number_density_m3: float | None = None
    tail_fraction: float | None = None

    def __post_init__(self) -> None:
        if self.energy_loss_eV is None:
            process_type = self.process_type.lower()
            threshold = float(self.threshold_eV or 0.0)
            if process_type in {"excitation", "ionization"}:
                self.energy_loss_eV = threshold
            elif process_type == "superelastic":
                self.energy_loss_eV = -abs(threshold)
            else:
                self.energy_loss_eV = 0.0

    @property
    def mixture_weighted_rate_m3_s(self) -> float:
        return self.rate_coefficient_m3_s * self.target_species_fraction

    @property
    def energy_loss_rate_coefficient_eV_m3_s(self) -> float:
        return self.rate_coefficient_m3_s * float(self.energy_loss_eV or 0.0)

    @property
    def frequency_s_inv(self) -> float | None:
        if self.gas_number_density_m3 is None:
            return None
        return self.gas_number_density_m3 * self.mixture_weighted_rate_m3_s

    @property
    def power_loss_eV_s(self) -> float | None:
        if self.gas_number_density_m3 is None:
            return None
        return (
            self.gas_number_density_m3
            * self.target_species_fraction
            * self.energy_loss_rate_coefficient_eV_m3_s
        )


@dataclass(slots=True)
class SwarmCaseResult:
    """A single solver result for one E/N case."""

    solver: str
    case_id: str
    e_over_n_Td: float
    mean_energy_eV: float
    net_ionization_frequency_s: float
    effective_townsend_m2: float
    transport: ElectronTransport
    energy_eV: np.ndarray = field(repr=False)
    eedf: np.ndarray = field(repr=False)  # normalized energy distribution, 1/eV
    energy_widths_eV: np.ndarray | None = field(default=None, repr=False)
    eedf_counts: np.ndarray | None = field(default=None, repr=False)
    eedf_effective_counts: np.ndarray | None = field(default=None, repr=False)
    rates: list[RateResult] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)
    diagnostics: dict[str, Any] = field(default_factory=dict)
    schema_version: str = "2"

    @property
    def drift_velocity_m_s(self) -> float:
        return self.transport.drift_velocity_m_s

    @property
    def mobility_m2_V_s(self) -> float:
        return self.transport.mobility_m2_V_s

    @property
    def reduced_mobility_m2_V_s_m3(self) -> float:
        return self.transport.reduced_mobility_m2_V_s_m3

    @property
    def diffusion_L_m2_s(self) -> float:
        return self.transport.diffusion_L_m2_s

    @property
    def diffusion_T_m2_s(self) -> float:
        return self.transport.diffusion_T_m2_s

    @property
    def reduced_diffusion_L_m2_s_m3(self) -> float:
        return self.transport.reduced_diffusion_L_m2_s_m3

    @property
    def reduced_diffusion_T_m2_s_m3(self) -> float:
        return self.transport.reduced_diffusion_T_m2_s_m3

    @property
    def reduced_electron_energy_mobility_m2_V_s_m3(self) -> float | None:
        return self.transport.reduced_electron_energy_mobility_m2_V_s_m3

    @property
    def reduced_electron_energy_diffusion_m2_s_m3(self) -> float | None:
        return self.transport.reduced_electron_energy_diffusion_m2_s_m3

    @property
    def electron_energy_mobility_m2_V_s(self) -> float | None:
        return self.transport.electron_energy_mobility_m2_V_s

    @property
    def electron_energy_diffusion_m2_s(self) -> float | None:
        return self.transport.electron_energy_diffusion_m2_s

    @property
    def eepf(self) -> np.ndarray:
        return eepf_from_eedf(self.energy_eV, self.eedf)


@dataclass(slots=True)
class SwarmRunResult:
    """All results from a YAML run."""

    cases: list[SwarmCaseResult]
    metadata: dict[str, Any] = field(default_factory=dict)

    def by_solver(self) -> dict[str, list[SwarmCaseResult]]:
        grouped: dict[str, list[SwarmCaseResult]] = {}
        for case in self.cases:
            grouped.setdefault(case.solver, []).append(case)
        return grouped
