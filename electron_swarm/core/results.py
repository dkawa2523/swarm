"""Shared result models for Monte Carlo and Boltzmann solvers."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from electron_swarm.core.transport import TransportSet

SUMMARY_METADATA_KEYS = (
    "physics_level",
    "solver_method",
    "angular_model",
    "angular_moment_source",
    "exact_dcs_based",
    "ordinary_integral_xs_closure",
    "direct_pn_operator",
    "lmax",
    "pn_residual",
    "negative_mass_fraction",
    "lmax1_regression_target",
    "effective_angular_scattering",
    "ionization_source_model",
    "ionization_source_treatment",
    "electron_electron_treatment",
    "electron_electron_transport_stale",
    "transport_definition",
    "magnetic_field_treatment",
    "magnetic_field_B_T",
    "magnetic_field_angle_EB_deg",
    "field_integrator",
    "tail_probability",
    "tail_rate_fraction_max",
    "dominant_tail_process",
    "energy_grid_tail_status",
)


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
    rates: list[RateResult] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)
    transport: TransportSet | None = None
    schema_version: str = "1.2"

    def summary_dict(self) -> dict[str, Any]:
        return {
            "solver": self.solver,
            "schema_version": self.schema_version,
            "case_id": self.case_id,
            "E_over_N_Td": self.e_over_n_Td,
            "mean_energy_eV": self.mean_energy_eV,
            "drift_velocity_m_s": self.drift_velocity_m_s,
            "mobility_m2_V_s": self.mobility_m2_V_s,
            "reduced_mobility_m2_V_s_m3": self.reduced_mobility_m2_V_s_m3,
            "diffusion_L_m2_s": self.diffusion_L_m2_s,
            "diffusion_T_m2_s": self.diffusion_T_m2_s,
            "reduced_diffusion_L_m2_s_m3": self.reduced_diffusion_L_m2_s_m3,
            "reduced_diffusion_T_m2_s_m3": self.reduced_diffusion_T_m2_s_m3,
            "net_ionization_frequency_s": self.net_ionization_frequency_s,
            "effective_townsend_m2": self.effective_townsend_m2,
            **{
                f"meta_{key}": self.metadata.get(key, "")
                for key in SUMMARY_METADATA_KEYS
            },
        }


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
