"""Shared result models for Monte Carlo and Boltzmann solvers."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from electron_swarm.core.transport import TransportSet


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
    mixture_weighted_rate_m3_s: float
    frequency_s_inv: float | None = None
    power_loss_eV_s: float | None = None
    tail_fraction: float | None = None


@dataclass(slots=True)
class SwarmCaseResult:
    """A single solver result for one E/N case."""

    solver: str
    case_id: str
    e_over_n_Td: float
    mean_energy_eV: float
    drift_velocity_m_s: float
    mobility_m2_V_s: float
    reduced_mobility_m2_V_s_m3: float
    diffusion_L_m2_s: float
    diffusion_T_m2_s: float
    reduced_diffusion_L_m2_s_m3: float
    reduced_diffusion_T_m2_s_m3: float
    net_ionization_frequency_s: float
    effective_townsend_m2: float
    energy_eV: np.ndarray = field(repr=False)
    eedf: np.ndarray = field(repr=False)  # normalized energy distribution, 1/eV
    eepf: np.ndarray = field(repr=False)  # EEPF-like eedf/sqrt(eV), eV^-3/2
    energy_widths_eV: np.ndarray | None = field(default=None, repr=False)
    eedf_counts: np.ndarray | None = field(default=None, repr=False)
    eedf_effective_counts: np.ndarray | None = field(default=None, repr=False)
    rates: list[RateResult] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)
    diagnostics: dict[str, Any] = field(default_factory=dict)
    transport: TransportSet | None = None
    schema_version: str = "2"


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
