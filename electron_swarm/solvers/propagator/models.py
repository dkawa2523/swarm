"""Typed numerical models for the DC energy-angle propagator."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import sparse

from electron_swarm.core.cross_sections import CrossSectionProcess


@dataclass(frozen=True, slots=True)
class PropagatorGrid:
    energy_edges_eV: np.ndarray
    energy_centers_eV: np.ndarray
    energy_widths_eV: np.ndarray
    speed_edges_m_s: np.ndarray
    speed_centers_m_s: np.ndarray
    theta_edges_rad: np.ndarray
    theta_centers_rad: np.ndarray
    mu_centers: np.ndarray
    mu_widths: np.ndarray
    solid_angles_sr: np.ndarray
    cell_volumes_v3: np.ndarray

    @property
    def energy_cells(self) -> int:
        return int(self.energy_centers_eV.size)

    @property
    def polar_cells(self) -> int:
        return int(self.theta_centers_rad.size)

    @property
    def cells(self) -> int:
        return self.energy_cells * self.polar_cells


@dataclass(frozen=True, slots=True)
class ElasticTransfer:
    species: str
    collision_frequency_s_inv: np.ndarray
    momentum_frequency_s_inv: np.ndarray
    mean_cosine: np.ndarray
    angular_kernel: np.ndarray | None
    isotropic: bool
    target_mass_amu: float
    xs_role: str
    collision_processes: tuple[CrossSectionProcess, ...]
    momentum_processes: tuple[CrossSectionProcess, ...]
    density_scale_m3: float


@dataclass(frozen=True, slots=True)
class InelasticTransfer:
    species: str
    process: str
    process_type: str
    frequency_s_inv: np.ndarray
    incident_energy_rate_eV_s_inv: np.ndarray
    daughter_energy_rate_eV_s_inv: np.ndarray
    gain_matrix_s_inv: sparse.csr_matrix
    daughter_count: int
    number_change_per_event: int
    threshold_eV: float


@dataclass(frozen=True, slots=True)
class CollisionOperatorData:
    elastic: tuple[ElasticTransfer, ...]
    inelastic: tuple[InelasticTransfer, ...]
    elastic_energy_generator_s_inv: sparse.csr_matrix
    elastic_energy_gain_s_inv: sparse.csr_matrix
    elastic_energy_outflow_s_inv: np.ndarray
    inelastic_outflow_s_inv: np.ndarray
    nonlocal_outflow_s_inv: np.ndarray
    elastic_A_eV_s: np.ndarray
    elastic_D_eV2_s: np.ndarray
    elastic_equilibrium_mass: np.ndarray
    isotropic_weights: np.ndarray
    memory_bytes: int
    elastic_xs_roles: tuple[str, ...]
    elastic_total_process_ids: tuple[str, ...]
    elastic_momentum_process_ids: tuple[str, ...]


@dataclass(slots=True)
class PropagatorConvergenceDiagnostics:
    converged: bool = False
    iterations: int = 0
    shape_change_L1: float = float("inf")
    growth_relative_change: float = float("inf")
    mean_energy_relative_change: float = float("inf")
    drift_relative_change: float = float("inf")
    operator_residual_L1: float = float("inf")
    normalization_error: float = float("inf")
    negative_population_mass: float = float("inf")
    tail_probability: float = float("inf")
    outer_acceleration_flux_fraction: float = float("inf")
    number_balance_residual: float = float("inf")
    stop_reason: str = "not_started"
    timings_s: dict[str, float] = field(default_factory=dict)

    def as_dict(self) -> dict[str, object]:
        return {
            "converged": self.converged,
            "iterations": self.iterations,
            "shape_change_L1": self.shape_change_L1,
            "growth_relative_change": self.growth_relative_change,
            "mean_energy_relative_change": self.mean_energy_relative_change,
            "drift_relative_change": self.drift_relative_change,
            "operator_residual_L1": self.operator_residual_L1,
            "normalization_error": self.normalization_error,
            "negative_population_mass": self.negative_population_mass,
            "tail_probability": self.tail_probability,
            "outer_acceleration_flux_fraction": (
                self.outer_acceleration_flux_fraction
            ),
            "number_balance_residual": self.number_balance_residual,
            "stop_reason": self.stop_reason,
            "timings_s": dict(self.timings_s),
        }


@dataclass(frozen=True, slots=True)
class PropagatorSteadySolution:
    population: np.ndarray
    growth_frequency_s_inv: float
    diagnostics: PropagatorConvergenceDiagnostics
