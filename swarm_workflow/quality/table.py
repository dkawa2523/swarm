"""Source-specific schemas for the canonical quality table.

Every solver exports the same qualification outcome fields.  Numerical
evidence is intentionally solver-specific: deterministic convergence evidence
does not belong in a Monte Carlo table, and Monte Carlo ensemble evidence does
not belong in deterministic solver tables.
"""

from __future__ import annotations

from dataclasses import dataclass


COMMON_QUALITY_COLUMNS = (
    "E_over_N_Td",
    "E_over_N_V_m2",
    "passed",
    "aggregate_quality_passed",
    "failure_reasons_json",
    "aggregate_failure_reasons_json",
    "eedf_normalization_error",
    "valid_replicates",
    "uncertainty_available",
    "required_rate_min_process_peak_fraction",
    "quality_policy_reevaluated",
    "solver_diagnostics_available",
    "solver_diagnostics_passed",
    "quality_source",
)

DETERMINISTIC_QUALITY_EVIDENCE_COLUMNS = (
    "solver_converged",
    "solver_iterations_max",
    "solver_residual_L1_max",
    "solver_residual_tolerance",
    "solver_tail_probability_max",
    "solver_tail_probability_target",
    "solver_edge_to_peak_max",
    "solver_edge_to_peak_target",
    "solver_grid_max_eV_max",
    "solver_grid_max_limit_eV",
    "solver_grid_limit_hit",
)

MONTE_CARLO_QUALITY_EVIDENCE_COLUMNS = (
    "failure_axes_json",
    "qualification_profile",
    "active_closure_quality_passed",
    "active_closure_failure_reasons_json",
    "active_closure_failure_axes_json",
    "mobility_rse",
    "diffusion_L_rse",
    "diffusion_T_rse",
    "energy_mobility_rse",
    "energy_diffusion_L_rse",
    "energy_diffusion_T_rse",
    "max_major_rate_rse",
    "solver_transport_qualified",
    "solver_transport_replicates",
    "solver_origin_stationarity_limiting_field",
    "solver_origin_stationarity_limiting_relative_ci95_bound",
    "solver_mean_energy_stationarity_relative_ci95_bound",
    "solver_mobility_stationarity_absolute_log_drift",
    "solver_lag_convergence_max_relative_ci95_bound",
    "solver_lag_supplemental_replicates",
    "solver_lag_supplemental_max_relative_ci95_bound",
    "solver_transport_mean_energy_max_relative_ci95_bound",
    "solver_population_growth_gate_mode",
    "solver_population_growth_max_relative_ci95_bound",
    "solver_population_growth_poisson_interval_max_ratio",
    "solver_population_growth_sparse_max_metric",
    "solver_population_growth_pooled_event_count",
    "solver_population_growth_pooled_exposure_s",
    "solver_population_growth_pooled_frequency_ci95_low_s_inv",
    "solver_population_growth_pooled_frequency_ci95_high_s_inv",
    "solver_population_growth_pooled_direct_within_ci95",
    "solver_population_growth_pooled_population_within_ci95",
    "solver_origin_stationarity_limiting_relative_tolerance",
    "solver_lag_convergence_relative_tolerance",
    "solver_transport_mean_energy_relative_tolerance",
    "solver_population_growth_relative_tolerance",
)


@dataclass(frozen=True, slots=True)
class QualityTableSchema:
    """Column contract for one canonical solver source."""

    source: str
    evidence_kind: str
    evidence_columns: tuple[str, ...]

    @property
    def columns(self) -> tuple[str, ...]:
        return (*COMMON_QUALITY_COLUMNS, *self.evidence_columns)


_QUALITY_SCHEMAS = {
    "two_term": QualityTableSchema(
        source="two_term",
        evidence_kind="deterministic_convergence",
        evidence_columns=DETERMINISTIC_QUALITY_EVIDENCE_COLUMNS,
    ),
    "propagator": QualityTableSchema(
        source="propagator",
        evidence_kind="deterministic_convergence",
        evidence_columns=DETERMINISTIC_QUALITY_EVIDENCE_COLUMNS,
    ),
    "monte_carlo": QualityTableSchema(
        source="monte_carlo",
        evidence_kind="independent_replica_statistics",
        evidence_columns=MONTE_CARLO_QUALITY_EVIDENCE_COLUMNS,
    ),
}


def quality_table_schema(source: str) -> QualityTableSchema:
    """Return the canonical quality-table schema for ``source``."""

    try:
        return _QUALITY_SCHEMAS[source]
    except KeyError as exc:
        raise ValueError(f"unsupported quality-table source: {source}") from exc
