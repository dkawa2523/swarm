"""Shared schemas and typed contracts for workflow table materialization."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

FORMAT_VERSION = 1
TD_TO_V_M2 = 1.0e-21
TWO_TERM_SOLVER = "two_term"
PROPAGATOR_SOLVER = "propagator"
MC_SOLVER = "monte_carlo"
SOURCE_CHOICES = (TWO_TERM_SOLVER, PROPAGATOR_SOLVER, MC_SOLVER)
MC_QUALIFICATION_FULL_TRANSPORT = "full_transport"
MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA = "function_eedf_restricted_lmea"
MC_QUALIFICATION_PROFILES = (
    MC_QUALIFICATION_FULL_TRANSPORT,
    MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA,
)
RATE_EVIDENCE_TABLE = "rate_evidence.csv"
ELASTIC_ENERGY_LOSS_TABLE = "elastic_energy_loss_vs_mean_energy.csv"
COLLISION_RATE_KERNEL_TABLE = "collision_rate_kernels.csv"
COLLISION_RATE_KERNEL_COLUMNS = (
    "species",
    "process",
    "process_type",
    "threshold_eV",
    "electron_energy_eV",
    "cross_section_m2",
    "high_energy_extrapolation",
)
TWO_TERM_TEMPORAL_GROWTH_TRANSPORT_COLUMNS = (
    "gas_number_density_m3",
    "temporal_growth_frequency_s_inv",
    "reduced_temporal_growth_frequency_m3_s",
)
TWO_TERM_MOMENTUM_CROSS_SECTION_FLOOR_M2 = 1.0e-24
ZERO_EVENT_CONFIDENCE = 0.95

CASE_COLUMNS = (
    "E_over_N_Td",
    "E_over_N_V_m2",
    "mean_energy_eV",
    "drift_velocity_m_s",
    "reduced_mobility_m2_V_s_m3",
    "reduced_diffusion_L_m2_s_m3",
    "reduced_diffusion_T_m2_s_m3",
    "reduced_electron_energy_mobility_m2_V_s_m3",
    "reduced_electron_energy_diffusion_m2_s_m3",
    "reduced_electron_energy_diffusion_L_m2_s_m3",
    "reduced_electron_energy_diffusion_T_m2_s_m3",
)
RATE_COLUMNS = (
    "E_over_N_Td",
    "E_over_N_V_m2",
    "species",
    "process",
    "process_type",
    "threshold_eV",
    "target_species_fraction",
    "rate_coefficient_m3_s",
    "mixture_weighted_rate_m3_s",
    "reduced_townsend_m2",
    "mixture_weighted_reduced_townsend_m2",
)
ENERGY_LOSS_COLUMNS = (
    "E_over_N_Td",
    "E_over_N_V_m2",
    "species",
    "process",
    "process_type",
    "threshold_eV",
    "target_species_fraction",
    "energy_loss_eV",
    "energy_loss_rate_coefficient_eV_m3_s",
)
EEDF_COLUMNS = (
    "electron_energy_eV",
    "energy_width_eV",
    "E_over_N_Td",
    "E_over_N_V_m2",
    "mean_energy_eV",
    "eedf",
    "pooled_effective_sample_count",
)


RATE_EVIDENCE_COLUMNS = (
    "E_over_N_Td",
    "E_over_N_V_m2",
    "mean_energy_eV",
    "species",
    "process",
    "process_type",
    "threshold_eV",
    "target_species_fraction",
    "rate_coefficient_mean_m3_s",
    "rate_coefficient_standard_error_m3_s",
    "rate_coefficient_relative_standard_error",
    "rate_coefficient_ci95_low_m3_s",
    "rate_coefficient_ci95_high_m3_s",
    "ci95_critical_value",
    "estimate_status",
    "valid_replicates",
    "uncertainty_available",
    "pooled_event_count",
    "pooled_target_exposure_s_m3",
    "pooled_all_zero_upper_95_m3_s",
    "pooled_zero_event_status",
)
ELASTIC_ENERGY_LOSS_COLUMNS = (
    "mean_energy_eV",
    "E_over_N_Td",
    "E_over_N_V_m2",
    "elastic_energy_loss_rate_coefficient_eV_m3_s",
    "elastic_energy_loss_standard_error_eV_m3_s",
    "elastic_energy_loss_relative_standard_error",
    "elastic_energy_loss_ci95_low_eV_m3_s",
    "elastic_energy_loss_ci95_high_eV_m3_s",
    "ci95_critical_value",
    "estimate_status",
    "valid_replicates",
    "uncertainty_available",
)

UNITS = {
    "E_over_N_Td": "Td",
    "E_over_N_V_m2": "V m^2",
    "mean_energy_eV": "eV",
    "electron_energy_eV": "eV",
    "energy_width_eV": "eV",
    "pooled_effective_sample_count": "1",
    "drift_velocity_m_s": "m/s",
    "reduced_mobility_m2_V_s_m3": "1/(V m s)",
    "reduced_diffusion_L_m2_s_m3": "1/(m s)",
    "reduced_diffusion_T_m2_s_m3": "1/(m s)",
    "reduced_electron_energy_mobility_m2_V_s_m3": "1/(V m s)",
    "reduced_electron_energy_diffusion_m2_s_m3": "1/(m s)",
    "reduced_electron_energy_diffusion_L_m2_s_m3": "1/(m s)",
    "reduced_electron_energy_diffusion_T_m2_s_m3": "1/(m s)",
    "gas_number_density_m3": "1/m^3",
    "temporal_growth_frequency_s_inv": "1/s",
    "reduced_temporal_growth_frequency_m3_s": "m^3/s",
    "effective_townsend_m2": "m^2",
    "rate_coefficient_m3_s": "m^3/s",
    "mixture_weighted_rate_m3_s": "m^3/s",
    "reduced_townsend_m2": "m^2",
    "mixture_weighted_reduced_townsend_m2": "m^2",
    "threshold_eV": "eV",
    "energy_loss_eV": "eV",
    "energy_loss_rate_coefficient_eV_m3_s": "eV m^3/s",
    "elastic_energy_loss_rate_coefficient_eV_m3_s": "eV m^3/s",
    "elastic_energy_loss_standard_error_eV_m3_s": "eV m^3/s",
    "elastic_energy_loss_relative_standard_error": "1",
    "elastic_energy_loss_ci95_low_eV_m3_s": "eV m^3/s",
    "elastic_energy_loss_ci95_high_eV_m3_s": "eV m^3/s",
    "eedf": "eV^-1",
    "cross_section_m2": "m^2",
    "rate_coefficient_mean_m3_s": "m^3/s",
    "rate_coefficient_standard_error_m3_s": "m^3/s",
    "rate_coefficient_ci95_low_m3_s": "m^3/s",
    "rate_coefficient_ci95_high_m3_s": "m^3/s",
    "pooled_target_exposure_s_m3": "s/m^3",
    "pooled_all_zero_upper_95_m3_s": "m^3/s",
    "solver_population_growth_pooled_exposure_s": "s",
    "solver_population_growth_pooled_frequency_ci95_low_s_inv": "1/s",
    "solver_population_growth_pooled_frequency_ci95_high_s_inv": "1/s",
    "required_rate_min_process_peak_fraction": "1",
}


class TableBuildError(RuntimeError):
    """Raised when workflow tables cannot be built safely."""


class MonteCarloQualificationError(TableBuildError):
    """MC statistics do not yet provide an exportable interpolation support."""


@dataclass(frozen=True, slots=True)
class MonteCarloEvidence:
    quality: list[dict[str, Any]]
    estimator_schema: str | None
    aggregate_eligible: set[float]
    full_transport_eligible: set[float]
    active_failure_reasons: dict[float, list[str]]
    rate_failures: dict[float, list[str]]
    transport_eligible: set[float]
    eedf_rate_eligible: set[float]


@dataclass(frozen=True, slots=True)
class TableBuildSummary:
    database_path: Path
    output_directory: Path
    source: str
    mixtures: int
    unqualified_mixtures: int = 0


@dataclass(slots=True)
class TableDataset:
    source: str
    mixture_id: int
    mixture: list[dict[str, str | float]]
    cases: list[dict[str, Any]]
    rates: list[dict[str, Any]]
    eedf: list[dict[str, Any]]
    rate_evidence: list[dict[str, Any]] | None
    elastic_energy_loss: list[dict[str, Any]] | None
    elastic_energy_loss_metadata: dict[str, Any] | None
    quality: list[dict[str, Any]]
    source_policy: dict[str, Any]
