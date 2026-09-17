from __future__ import annotations

import csv
from hashlib import sha256
import json
import math
from pathlib import Path
import sqlite3

import pytest

from electron_swarm.solvers.monte_carlo.evidence import (
    DIRECT_MC_TRANSPORT_CADENCE_TRIAL_PERIODS,
    DIRECT_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
    DIRECT_MC_TRANSPORT_LAG_PLANES,
    DIRECT_MC_TRANSPORT_OBSERVATION_PLANES,
    WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
    WEIGHTED_MC_TRANSPORT_MIN_EFFECTIVE_LINEAGE_COUNT,
    WEIGHTED_MC_TRANSPORT_MIN_EFFECTIVE_LINEAGE_FRACTION,
    WEIGHTED_MC_TRANSPORT_MIN_BLOCKS,
)
from electron_swarm.solvers.monte_carlo.weighted_transport import (
    weighted_growth_transport_lag_plan,
)
from swarm_workflow.campaign.aggregate import (
    aggregate_workflow_results,
)
from swarm_workflow.quality.table import (
    COMMON_QUALITY_COLUMNS,
    quality_table_schema,
)
from swarm_workflow.quality.monte_carlo.direct_transport import (
    summarize_mc_transport_diagnostics,
)
from swarm_workflow.quality.monte_carlo.contracts import (
    MC_TRANSPORT_LAG_FIELDS,
    WEIGHTED_MC_TRANSPORT_DIFFUSION_IDENTIFICATION_RELATIVE_TOLERANCE,
    mc_function_eedf_restricted_lmea_failure_reasons,
)
from swarm_workflow.quality.monte_carlo.weighted_transport import (
    summarize_weighted_mc_transport_diagnostics,
)
from swarm_workflow.quality.policy import QualityThresholds, quality_thresholds_json
from swarm_workflow.tables.contracts import (
    ELASTIC_ENERGY_LOSS_TABLE,
    MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA,
    MonteCarloQualificationError,
    TableBuildError,
    TWO_TERM_TEMPORAL_GROWTH_TRANSPORT_COLUMNS,
)
from swarm_workflow.tables.energy_loss import _load_elastic_energy_loss
from swarm_workflow.tables.monte_carlo import (
    _mc_eedf_rate_consistency_failures,
    _validate_mc_censored_rate_relevance,
)
from swarm_workflow.tables.repository import _validate_mc_sampling_plan_against_cases
from swarm_workflow.tables import build_tables
from swarm_workflow.campaign.store import (
    MC_SAMPLING_PLAN_METADATA_KEY,
    WorkflowSchemaError,
    WorkflowStore,
    canonical_mc_sampling_plan_json,
)
from product_helpers import (
    insert_workflow_case,
    insert_workflow_eedf,
    insert_workflow_rate,
)


_WEIGHTED_GROWTH_LAG_PLAN = weighted_growth_transport_lag_plan()


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fp:
        return list(csv.DictReader(fp))


def _mc_rate_diagnostics(
    *,
    event_count: int,
    residence_time_s: float,
) -> dict[str, object]:
    return {
        "internal_monte_carlo_reaction_rates": {
            "estimator": "trajectory_time_average_sigma_v",
            "rates": [
                {
                    "species": "Ar",
                    "process": "ionization",
                    "process_type": "ionization",
                    "threshold_eV": 15.0,
                    "target_species_fraction": 1.0,
                    "event_sampling_enabled": True,
                    "event_count": event_count,
                    "event_observation_residence_time_s": residence_time_s,
                    "zero_event_confidence": 0.95,
                }
            ],
        }
    }


def _mc_transport_diagnostic(
    *,
    seed: int = 1,
    stationarity_change: float = 0.01,
    lag_change: float = 0.02,
) -> dict[str, object]:
    production = {
        "drift_velocity_m_s": 1.0,
        "mobility_m2_V_s": 1.0,
        "energy_mobility_m2_V_s": 4.0,
        "mean_energy_eV": 7.0,
    }
    diffusion = {
        "diffusion_L_m2_s": 2.0,
        "diffusion_T_m2_s": 3.0,
        "energy_diffusion_L_m2_s": 5.0,
        "energy_diffusion_T_m2_s": 6.0,
    }
    lag_estimates = {
        name: [
            value * (1.0 - 3.0 * lag_change),
            value * (1.0 - 2.0 * lag_change),
            value * (1.0 - lag_change),
            value,
        ]
        for name, value in diffusion.items()
    }
    cadence_s = 1.0e-9
    return {
        "estimator": "direct_mc_fixed_lag_block_helfand_energy_flux",
        "estimator_schema_version": DIRECT_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
        "energy_transport_status": ("direct_mc_candidate_for_ensemble_qualification"),
        "finite": True,
        "strictly_positive": True,
        "multiple_time_origins": True,
        "complete_observation": True,
        "origin_stationarity_window_mean_early": {
            "mean_energy_eV": production["mean_energy_eV"]
            * (1.0 - stationarity_change),
        },
        "origin_stationarity_window_mean_late": {
            "mean_energy_eV": production["mean_energy_eV"],
        },
        "component_estimators": {
            "drift_mobility": (
                "production_trajectory_displacement_over_residence_time"
            ),
            "particle_diffusion": "block_helfand_mean_square_displacement",
            "energy_mobility": "production_residence_energy_current",
            "energy_diffusion": ("restricted_density_packet_energy_flux_fixed_lag"),
        },
        "production_estimates": production,
        "lag_scan": {
            "lag_planes": list(DIRECT_MC_TRANSPORT_LAG_PLANES),
            "lag_s": [value * cadence_s for value in DIRECT_MC_TRANSPORT_LAG_PLANES],
            "estimates": lag_estimates,
            "completed_origins": [
                DIRECT_MC_TRANSPORT_OBSERVATION_PLANES - value
                for value in DIRECT_MC_TRANSPORT_LAG_PLANES
            ],
            "selected_lag_index": len(DIRECT_MC_TRANSPORT_LAG_PLANES) - 1,
        },
        "transport_sampling": {
            "lag_grid_source": ("fixed_physical_time_independent_of_production_length"),
            "origin_cadence_s": cadence_s,
            "origin_cadence_trial_periods": (DIRECT_MC_TRANSPORT_CADENCE_TRIAL_PERIODS),
            "rolling_origins": True,
            "production_length_controls_lag": False,
            "warmup_accumulators_reset": True,
            "observation_planes_requested": DIRECT_MC_TRANSPORT_OBSERVATION_PLANES,
            "observation_planes_complete": DIRECT_MC_TRANSPORT_OBSERVATION_PLANES,
            "minimum_observation_planes": DIRECT_MC_TRANSPORT_OBSERVATION_PLANES,
        },
        "mc_run_provenance": {
            "seed": seed,
            "particles": 128,
            "warmup_collisions": 100,
            "production_collisions": 2000,
            "tail_max_collisions": 2000,
            "tail_rate_rse_trigger": 0.25,
            "tail_collisions_executed": 0,
            "population_model": "fixed_particle_single_daughter",
            "transport_estimator": "single_field",
            "gas_number_density_m3": 1.0e22,
            "electric_field_V_m": 1.0e3,
            "trial_collision_frequency_s_inv": 1.0e9,
            "collision_clock": "global",
            "magnetic_enabled": False,
            "magnetic_B_T": 0.0,
            "magnetic_angle_EB_deg": 0.0,
            "angular_scattering_model": "isotropic",
            "ionization_source_model": "equal",
            "high_energy_extrapolation": "zero",
            "max_energy_limit_eV": 1000.0,
        },
        "_case_mean_energy_eV": production["mean_energy_eV"],
        "_case_reported_transport": {
            **production,
            **diffusion,
        },
    }


def _weighted_mc_transport_diagnostic(
    seed: int,
    *,
    correlation_lag: int = 64,
) -> dict[str, object]:
    lag_barriers = [
        correlation_lag // 4,
        correlation_lag // 2,
        correlation_lag,
        2 * correlation_lag,
    ]
    transport_block_barriers = 2 * correlation_lag
    production = {
        "drift_velocity_m_s": 1.0,
        "mobility_m2_V_s": 2.0,
        "diffusion_L_m2_s": 3.0,
        "diffusion_T_m2_s": 4.0,
        "energy_mobility_m2_V_s": 5.0,
        "energy_diffusion_L_m2_s": 6.0,
        "energy_diffusion_T_m2_s": 7.0,
        "mean_energy_eV": 8.0,
    }
    stationarity = {
        key: value for key, value in production.items() if key != "drift_velocity_m_s"
    }
    return {
        "estimator": "direct_mc_weighted_growth_block_lag_flux",
        "estimator_schema_version": (WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION),
        "energy_transport_status": ("direct_mc_candidate_for_ensemble_qualification"),
        "finite": True,
        "strictly_positive": True,
        "multiple_time_blocks": True,
        "complete_observation": True,
        "lineage_qualified": True,
        "time_stationarity_window_estimate_early": dict(stationarity),
        "time_stationarity_window_estimate_late": dict(stationarity),
        "component_estimators": {
            "drift_mobility": "normalized_weight_residence_velocity_moment",
            "particle_diffusion": (
                "synchronized_block_lag_position_velocity_flux_covariance"
            ),
            "energy_mobility": "normalized_weight_residence_energy_current",
            "energy_diffusion": (
                "synchronized_block_lag_restricted_density_packet_"
                "energy_current_covariance"
            ),
        },
        "production_estimates": production,
        "block_lag_sampling": {
            "physical_time_barriers": True,
            "resampling_at_common_time_only": True,
            "barrier_cadence_s": 1.0e-9,
            "configured_correlation_lag_barriers": correlation_lag,
            "transport_block_barriers": transport_block_barriers,
            "lag_barriers": lag_barriers,
            "complete_blocks": WEIGHTED_MC_TRANSPORT_MIN_BLOCKS,
            "minimum_blocks": WEIGHTED_MC_TRANSPORT_MIN_BLOCKS,
            "positions_reset_each_block": True,
            "warmup_accumulators_reset": True,
        },
        "stationarity_sampling": {
            "window_definition": (
                "first_half_vs_second_half_complete_transport_blocks"
            ),
            "state_moment_aggregation": (
                "sum_raw_weighted_residence_moments_before_ratio"
            ),
            "state_estimator": "same_residence_ratio_formulas_as_production",
            "diffusion_estimator": (
                "mean_fixed_production_lag_snapshot_covariance_per_window"
            ),
            "diffusion_lag_barriers": (correlation_lag),
            "early_block_count": WEIGHTED_MC_TRANSPORT_MIN_BLOCKS // 2,
            "late_block_count": WEIGHTED_MC_TRANSPORT_MIN_BLOCKS // 2,
            "odd_center_block_excluded": False,
        },
        "lag_scan": {
            "lag_barriers": lag_barriers,
            "lag_s": [lag * 1.0e-9 for lag in lag_barriers],
            "estimates": {
                key: [value] * len(lag_barriers)
                for key, value in production.items()
                if "diffusion" in key
            },
            "completed_blocks": [WEIGHTED_MC_TRANSPORT_MIN_BLOCKS] * len(lag_barriers),
            "production_lag_barriers": correlation_lag,
            "hard_convergence_pair_barriers": [
                correlation_lag // 2,
                correlation_lag,
            ],
            "supplemental_pair_barriers": [
                correlation_lag,
                transport_block_barriers,
            ],
        },
        "lineage_sampling": {
            "lineage_horizon": "each_sampled_lag",
            "lag_barriers": lag_barriers,
            "minimum_effective_lineage_count": [128.0] * len(lag_barriers),
            "minimum_effective_lineage_fraction": [0.25] * len(lag_barriers),
            "production_lag_barriers": correlation_lag,
            "qualification_lag_barriers": [
                correlation_lag // 2,
                correlation_lag,
            ],
            "required_effective_lineage_count": (
                WEIGHTED_MC_TRANSPORT_MIN_EFFECTIVE_LINEAGE_COUNT
            ),
            "required_effective_lineage_fraction": (
                WEIGHTED_MC_TRANSPORT_MIN_EFFECTIVE_LINEAGE_FRACTION
            ),
        },
        "transport_definition": "flux",
        "bulk_transport_identified": False,
        "cross_gradient_response_identified": False,
        "population_growth_consistency": {
            "model": "explicit_ionization_branching_no_attachment",
            "production_elapsed_time_s": 1.0e-6,
            "cumulative_log_growth": 10.0,
            "population_log_growth_frequency_s_inv": 1.0e7,
            "direct_ionization_frequency_s_inv": 1.0e7,
            "event_evidence_model": "raw_macro_poisson",
            "event_count_ionization_frequency_s_inv": 1.0e7,
            "ionization_event_count": 1000,
            "weighted_ionization_event_count": 1000.0,
            "secondary_electron_events": 1000,
            "event_count_matches_secondary": True,
            "event_observation_residence_time_s": 1.0e-4,
            "event_exposure_consistent_across_processes": True,
            "zero_event_frequency_upper_95_s_inv": None,
        },
        "mc_run_provenance": {
            "seed": seed,
            "particles": 512,
            "warmup_collisions": 100,
            "production_collisions": (
                transport_block_barriers * WEIGHTED_MC_TRANSPORT_MIN_BLOCKS
            ),
            "tail_max_collisions": (
                transport_block_barriers * WEIGHTED_MC_TRANSPORT_MIN_BLOCKS
            ),
            "tail_rate_rse_trigger": 0.25,
            "tail_collisions_executed": 0,
            "population_model": "weighted_branching",
            "transport_estimator": "single_field",
            "gas_number_density_m3": 1.0e22,
            "electric_field_V_m": 1.0e3,
            "trial_collision_frequency_s_inv": 1.0e9,
            "collision_clock": "global",
            "magnetic_enabled": False,
            "magnetic_B_T": 0.0,
            "magnetic_angle_EB_deg": 0.0,
            "angular_scattering_model": "isotropic",
            "ionization_source_model": "equal",
            "high_energy_extrapolation": "zero",
            "max_energy_limit_eV": 1000.0,
        },
        "_case_mean_energy_eV": production["mean_energy_eV"],
        "_case_reported_transport": {
            key: value for key, value in production.items() if key != "mean_energy_eV"
        },
    }


def test_weighted_mc_transport_qualification_accepts_synchronized_flux() -> None:
    summary = summarize_weighted_mc_transport_diagnostics(
        [_weighted_mc_transport_diagnostic(seed) for seed in (11, 12, 13)]
    )

    assert summary["solver_transport_qualified"] == 1
    assert summary["solver_failure_reasons"] == []
    assert summary["solver_population_growth_gate_mode"] == ("positive_paired_log")
    assert summary["solver_lag_convergence_max_relative_ci95_bound"] == 0.0
    assert summary["solver_lag_convergence_relative_tolerance"] == (
        WEIGHTED_MC_TRANSPORT_DIFFUSION_IDENTIFICATION_RELATIVE_TOLERANCE
    )
    assert summary["solver_lag_supplemental_replicates"] == 3
    assert summary["solver_lag_supplemental_max_relative_ci95_bound"] == 0.0
    assert summary["solver_population_growth_max_relative_ci95_bound"] == 0.0
    assert summary["solver_population_growth_poisson_interval_max_ratio"] < 1.0
    assert summary["solver_population_growth_sparse_max_metric"] > 0.1
    assert summary["solver_population_growth_pooled_event_count"] == 3000
    assert summary["solver_population_growth_pooled_exposure_s"] == pytest.approx(
        3.0e-4
    )
    assert (
        summary["solver_population_growth_pooled_frequency_ci95_low_s_inv"]
        < 1.0e7
        < summary["solver_population_growth_pooled_frequency_ci95_high_s_inv"]
    )
    assert summary["solver_population_growth_pooled_direct_within_ci95"] == 1
    assert summary["solver_population_growth_pooled_population_within_ci95"] == 1
    assert summary["solver_population_growth_relative_tolerance"] == 0.1
    assert summary["solver_origin_stationarity_limiting_field"] == ("mobility_m2_V_s")


def test_restricted_lmea_ignores_only_inactive_mc_transport_failures() -> None:
    row = {
        "failure_reasons_json": json.dumps(
            [
                "reduced_diffusion_L_rse_unavailable_or_exceeds_threshold",
                ("reduced_electron_energy_diffusion_L_m2_s_m3_not_finite_nonnegative"),
                "mc_weighted_growth_lag_not_converged:diffusion_L_m2_s",
                (
                    "mc_weighted_growth_stationarity_not_converged:"
                    "energy_mobility_m2_V_s"
                ),
            ]
        ),
        "valid_replicates": 4,
        "solver_transport_replicates": 4,
        "solver_diagnostics_available": 1,
        "mobility_rse": 0.02,
        "eedf_normalization_error": 1.0e-12,
        "max_major_rate_rse": 0.03,
        "solver_mean_energy_stationarity_relative_ci95_bound": 0.02,
        "solver_mobility_stationarity_absolute_log_drift": 0.02,
        "solver_transport_mean_energy_max_relative_ci95_bound": 0.02,
        "solver_transport_mean_energy_relative_tolerance": 0.1,
        "solver_population_growth_relative_tolerance": 0.1,
        "solver_population_growth_gate_mode": "positive_paired_log",
        "solver_population_growth_max_relative_ci95_bound": 0.02,
        "quality_source": "monte_carlo_aggregate_quality",
    }
    thresholds = QualityThresholds(
        mobility_rse=0.08,
        diffusion_rse=0.12,
        major_rate_rse=0.12,
    )

    assert MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA == (
        "function_eedf_restricted_lmea"
    )
    assert (
        mc_function_eedf_restricted_lmea_failure_reasons(
            row,
            thresholds=thresholds,
            estimator_schema=WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
        )
        == []
    )

    row["failure_reasons_json"] = json.dumps(
        ["mc_weighted_growth_stationarity_not_converged:mobility_m2_V_s"]
    )
    assert (
        mc_function_eedf_restricted_lmea_failure_reasons(
            row,
            thresholds=thresholds,
            estimator_schema=WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
        )
        == []
    )

    row["solver_mobility_stationarity_absolute_log_drift"] = 0.11
    assert mc_function_eedf_restricted_lmea_failure_reasons(
        row,
        thresholds=thresholds,
        estimator_schema=WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
    ) == ["solver_mobility_stationarity_absolute_log_drift"]

    row["solver_mobility_stationarity_absolute_log_drift"] = 0.02
    row["solver_mean_energy_stationarity_relative_ci95_bound"] = 0.11
    assert mc_function_eedf_restricted_lmea_failure_reasons(
        row,
        thresholds=thresholds,
        estimator_schema=WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
    ) == ["solver_mean_energy_stationarity_relative_ci95_bound"]

    row["solver_mean_energy_stationarity_relative_ci95_bound"] = 0.02
    row["mobility_rse"] = 0.04
    assert mc_function_eedf_restricted_lmea_failure_reasons(
        row,
        thresholds=thresholds,
        estimator_schema=WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
    ) == ["mobility_relative_ci95_half_width"]

    row["mobility_rse"] = 0.02
    row["solver_mean_energy_stationarity_relative_ci95_bound"] = None
    assert mc_function_eedf_restricted_lmea_failure_reasons(
        row,
        thresholds=thresholds,
        estimator_schema=WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
    ) == ["solver_mean_energy_stationarity_relative_ci95_bound"]


def test_weighted_mc_transport_accepts_configured_lag_256() -> None:
    diagnostics = [
        _weighted_mc_transport_diagnostic(seed, correlation_lag=256)
        for seed in (11, 12, 13)
    ]

    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 1
    assert diagnostics[0]["block_lag_sampling"]["lag_barriers"] == [
        64,
        128,
        256,
        512,
    ]


def test_weighted_mc_transport_rejects_partial_production_block() -> None:
    diagnostics = [_weighted_mc_transport_diagnostic(seed) for seed in (11, 12, 13)]
    for item in diagnostics:
        item["mc_run_provenance"]["production_collisions"] += 1

    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 0
    assert (
        "mc_weighted_growth_ensemble_evidence_missing"
        in summary["solver_failure_reasons"]
    )


def test_weighted_mc_transport_rejects_reported_plan_relation_mismatch() -> None:
    diagnostics = [
        _weighted_mc_transport_diagnostic(seed, correlation_lag=256)
        for seed in (11, 12, 13)
    ]
    diagnostics[0]["lag_scan"]["hard_convergence_pair_barriers"] = [32, 256]

    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 0
    assert (
        "mc_weighted_growth_ensemble_evidence_missing"
        in summary["solver_failure_reasons"]
    )


def test_weighted_mc_sparse_growth_qualifies_transport_only() -> None:
    diagnostics = [_weighted_mc_transport_diagnostic(seed) for seed in (11, 12, 13)]
    for item in diagnostics:
        growth = item["population_growth_consistency"]
        growth["cumulative_log_growth"] = 0.0
        growth["population_log_growth_frequency_s_inv"] = 0.0
        growth["direct_ionization_frequency_s_inv"] = 1000.0
        growth["event_count_ionization_frequency_s_inv"] = 0.0
        growth["ionization_event_count"] = 0
        growth["secondary_electron_events"] = 0
        growth["zero_event_frequency_upper_95_s_inv"] = (
            -math.log(0.05) / growth["event_observation_residence_time_s"]
        )

    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 1
    assert summary["solver_population_growth_gate_mode"] == (
        "sparse_negligible_transport"
    )
    assert summary["solver_population_growth_max_relative_ci95_bound"] is None
    assert summary["solver_population_growth_poisson_interval_max_ratio"] < 1.0
    assert summary["solver_population_growth_sparse_max_metric"] < 0.1
    assert summary["solver_population_growth_pooled_event_count"] == 0
    assert summary[
        "solver_population_growth_pooled_frequency_ci95_high_s_inv"
    ] == pytest.approx(-math.log(0.05) / 3.0e-4)

    for item in diagnostics:
        item["population_growth_consistency"]["direct_ionization_frequency_s_inv"] = (
            1.0e6
        )
    rejected = summarize_weighted_mc_transport_diagnostics(diagnostics)
    assert rejected["solver_transport_qualified"] == 0
    assert (
        "mc_weighted_growth_sparse_poisson_inconsistent"
        in rejected["solver_failure_reasons"]
    )
    assert summary["solver_origin_stationarity_limiting_relative_tolerance"] == 0.1


def test_weighted_ensemble_growth_uses_weighted_events_without_poisson() -> None:
    diagnostics = [_weighted_mc_transport_diagnostic(seed) for seed in (11, 12, 13)]
    for item in diagnostics:
        growth = item["population_growth_consistency"]
        exposure = float(growth["event_observation_residence_time_s"])
        growth["event_evidence_model"] = "weighted_ensemble_correlated_events"
        growth["weighted_ionization_event_count"] = 1000.0
        growth["event_count_ionization_frequency_s_inv"] = 1000.0 / exposure
        growth["zero_event_frequency_upper_95_s_inv"] = None

    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 1
    assert summary["solver_population_growth_gate_mode"] == (
        "weighted_ensemble_positive_paired_log"
    )
    assert summary["solver_population_growth_poisson_interval_max_ratio"] is None
    assert summary["solver_population_growth_max_relative_ci95_bound"] == 0.0


def test_weighted_mc_positive_growth_uses_paired_log_not_poisson_gate() -> None:
    diagnostics = [_weighted_mc_transport_diagnostic(seed) for seed in (11, 12, 13)]
    for item in diagnostics:
        growth = item["population_growth_consistency"]
        growth["cumulative_log_growth"] = 10.8
        growth["population_log_growth_frequency_s_inv"] = 1.08e7
        growth["direct_ionization_frequency_s_inv"] = 1.08e7

    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 1
    assert summary["solver_population_growth_gate_mode"] == "positive_paired_log"
    assert summary["solver_population_growth_max_relative_ci95_bound"] == pytest.approx(
        abs(math.log(1.0e7 / 1.08e7))
    )
    assert summary["solver_population_growth_pooled_direct_within_ci95"] == 0
    assert summary["solver_population_growth_pooled_population_within_ci95"] == 0


def test_weighted_mc_sparse_mixed_counts_use_pooled_interval() -> None:
    diagnostics = [_weighted_mc_transport_diagnostic(seed) for seed in (11, 12, 13)]
    for index, item in enumerate(diagnostics):
        growth = item["population_growth_consistency"]
        if index == 0:
            growth["cumulative_log_growth"] = 0.0
            growth["population_log_growth_frequency_s_inv"] = 0.0
            growth["direct_ionization_frequency_s_inv"] = 5000.0
            growth["event_count_ionization_frequency_s_inv"] = 0.0
            growth["ionization_event_count"] = 0
            growth["secondary_electron_events"] = 0
            growth["zero_event_frequency_upper_95_s_inv"] = (
                -math.log(0.05) / growth["event_observation_residence_time_s"]
            )
        else:
            growth["cumulative_log_growth"] = 0.01
            growth["population_log_growth_frequency_s_inv"] = 1.0e4
            growth["direct_ionization_frequency_s_inv"] = 7500.0
            growth["event_count_ionization_frequency_s_inv"] = 1.0e4
            growth["ionization_event_count"] = 1
            growth["secondary_electron_events"] = 1
            growth["zero_event_frequency_upper_95_s_inv"] = None

    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 1
    assert summary["solver_population_growth_gate_mode"] == (
        "sparse_negligible_transport"
    )
    assert summary["solver_population_growth_pooled_event_count"] == 2
    assert summary["solver_population_growth_pooled_direct_within_ci95"] == 1
    assert summary["solver_population_growth_pooled_population_within_ci95"] == 1
    assert summary["solver_population_growth_max_relative_ci95_bound"] is None


def test_weighted_mc_sparse_positive_counts_do_not_gate_paired_log() -> None:
    diagnostics = [_weighted_mc_transport_diagnostic(seed) for seed in (11, 12, 13)]
    for item in diagnostics:
        growth = item["population_growth_consistency"]
        growth["cumulative_log_growth"] = 0.02
        growth["population_log_growth_frequency_s_inv"] = 2.0e4
        growth["direct_ionization_frequency_s_inv"] = 2.0e4
        growth["event_count_ionization_frequency_s_inv"] = 1.0e4
        growth["ionization_event_count"] = 1
        growth["secondary_electron_events"] = 1
        growth["zero_event_frequency_upper_95_s_inv"] = None

    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 1
    assert summary["solver_population_growth_gate_mode"] == (
        "sparse_negligible_transport"
    )
    assert summary["solver_population_growth_sparse_max_metric"] < 0.1
    assert summary["solver_population_growth_max_relative_ci95_bound"] == pytest.approx(
        math.log(2.0)
    )


def test_weighted_mc_zero_growth_separates_branching_negligibility() -> None:
    diagnostics = [_weighted_mc_transport_diagnostic(seed) for seed in (11, 12, 13)]
    for item in diagnostics:
        growth = item["population_growth_consistency"]
        growth["cumulative_log_growth"] = 0.0
        growth["population_log_growth_frequency_s_inv"] = 0.0
        growth["direct_ionization_frequency_s_inv"] = 1000.0
        growth["event_count_ionization_frequency_s_inv"] = 0.0
        growth["ionization_event_count"] = 0
        growth["secondary_electron_events"] = 0
        growth["zero_event_frequency_upper_95_s_inv"] = (
            -math.log(0.05) / growth["event_observation_residence_time_s"]
        )
        sampling = item["block_lag_sampling"]
        sampling["barrier_cadence_s"] = 1.0e-5
        item["lag_scan"]["lag_s"] = [
            lag * 1.0e-5 for lag in _WEIGHTED_GROWTH_LAG_PLAN.lag_barriers
        ]
        item["mc_run_provenance"]["trial_collision_frequency_s_inv"] = 1.0e5

    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 0
    assert summary["solver_population_growth_poisson_interval_max_ratio"] < 1.0
    assert summary["solver_population_growth_sparse_max_metric"] > 0.1
    assert (
        "mc_weighted_growth_branching_not_negligible"
        in summary["solver_failure_reasons"]
    )


def test_weighted_mc_stationarity_reports_limiting_field_contract() -> None:
    diagnostics = [_weighted_mc_transport_diagnostic(seed) for seed in (11, 12, 13)]
    for item in diagnostics:
        early = item["time_stationarity_window_estimate_early"]
        late = item["time_stationarity_window_estimate_late"]
        late["diffusion_L_m2_s"] = early["diffusion_L_m2_s"] * math.exp(0.115)

    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 1
    assert summary["solver_origin_stationarity_limiting_field"] == ("diffusion_L_m2_s")
    assert summary[
        "solver_origin_stationarity_limiting_relative_ci95_bound"
    ] == pytest.approx(0.115)
    assert summary["solver_origin_stationarity_limiting_relative_tolerance"] == 0.25


def test_weighted_mc_stationarity_selects_worst_normalized_margin() -> None:
    diagnostics = [_weighted_mc_transport_diagnostic(seed) for seed in (11, 12, 13)]
    for item in diagnostics:
        early = item["time_stationarity_window_estimate_early"]
        late = item["time_stationarity_window_estimate_late"]
        late["diffusion_L_m2_s"] = early["diffusion_L_m2_s"] * math.exp(0.2)
        late["mobility_m2_V_s"] = early["mobility_m2_V_s"] * math.exp(0.11)

    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 0
    assert summary["solver_origin_stationarity_limiting_field"] == ("mobility_m2_V_s")
    assert summary[
        "solver_origin_stationarity_limiting_relative_ci95_bound"
    ] == pytest.approx(0.11)
    assert summary["solver_origin_stationarity_limiting_relative_tolerance"] == 0.1


def test_weighted_mc_transport_qualification_rejects_lag64_lineage_collapse() -> None:
    diagnostics = [_weighted_mc_transport_diagnostic(seed) for seed in (11, 12, 13)]
    diagnostics[0]["lineage_qualified"] = False
    lineage = diagnostics[0]["lineage_sampling"]
    assert isinstance(lineage, dict)
    production_index = list(_WEIGHTED_GROWTH_LAG_PLAN.lag_barriers).index(
        _WEIGHTED_GROWTH_LAG_PLAN.production_lag_barriers
    )
    lineage["minimum_effective_lineage_count"][production_index] = 127.0
    lineage["minimum_effective_lineage_fraction"][production_index] = 127.0 / 512.0
    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 0
    assert "mc_weighted_growth_lineage_degenerate" in summary["solver_failure_reasons"]


def test_weighted_mc_transport_qualification_rejects_lag32_lineage_collapse() -> None:
    diagnostics = [_weighted_mc_transport_diagnostic(seed) for seed in (11, 12, 13)]
    diagnostics[0]["lineage_qualified"] = False
    lineage = diagnostics[0]["lineage_sampling"]
    assert isinstance(lineage, dict)
    lag32_index = list(_WEIGHTED_GROWTH_LAG_PLAN.lag_barriers).index(32)
    lineage["minimum_effective_lineage_count"][lag32_index] = 127.0
    lineage["minimum_effective_lineage_fraction"][lag32_index] = 127.0 / 512.0

    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 0
    assert "mc_weighted_growth_lineage_degenerate" in summary["solver_failure_reasons"]


def test_weighted_mc_transport_qualification_rejects_fraction_degeneracy() -> None:
    diagnostics = [_weighted_mc_transport_diagnostic(seed) for seed in (11, 12, 13)]
    diagnostics[0]["lineage_qualified"] = False
    diagnostics[0]["mc_run_provenance"]["particles"] = 4096
    lineage = diagnostics[0]["lineage_sampling"]
    assert isinstance(lineage, dict)
    lineage["minimum_effective_lineage_fraction"] = [
        count / 4096.0 for count in lineage["minimum_effective_lineage_count"]
    ]

    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 0
    assert "mc_weighted_growth_lineage_degenerate" in summary["solver_failure_reasons"]


def test_weighted_mc_transport_qualification_scales_with_absolute_lineage_ess() -> None:
    diagnostics = [_weighted_mc_transport_diagnostic(seed) for seed in (11, 12, 13)]
    for item in diagnostics:
        item["mc_run_provenance"]["particles"] = 2048
        lineage = item["lineage_sampling"]
        assert isinstance(lineage, dict)
        lineage["minimum_effective_lineage_count"] = [160.0] * len(
            _WEIGHTED_GROWTH_LAG_PLAN.lag_barriers
        )
        lineage["minimum_effective_lineage_fraction"] = [160.0 / 2048.0] * len(
            _WEIGHTED_GROWTH_LAG_PLAN.lag_barriers
        )

    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 1
    assert summary["solver_failure_reasons"] == []


def test_weighted_mc_transport_qualification_rejects_lag_drift() -> None:
    diagnostics = [_weighted_mc_transport_diagnostic(seed) for seed in (11, 12, 13)]
    for item in diagnostics:
        item["lag_scan"]["estimates"]["diffusion_L_m2_s"][1] = 1.0

    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 0
    assert (
        "mc_weighted_growth_lag_not_converged:diffusion_L_m2_s"
        in summary["solver_failure_reasons"]
    )


def test_weighted_mc_lag128_drift_gates_when_lineage_is_qualified() -> None:
    diagnostics = [_weighted_mc_transport_diagnostic(seed) for seed in (11, 12, 13)]
    for item in diagnostics:
        item["lag_scan"]["estimates"]["diffusion_L_m2_s"][3] *= 4.0

    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 0
    assert summary["solver_lag_supplemental_replicates"] == 3
    assert summary["solver_lag_supplemental_max_relative_ci95_bound"] == pytest.approx(
        math.log(4.0)
    )
    assert (
        "mc_weighted_growth_lag_not_converged:diffusion_L_m2_s"
        in summary["solver_failure_reasons"]
    )


def test_weighted_mc_supplemental_128_requires_qualified_lineage() -> None:
    diagnostics = [_weighted_mc_transport_diagnostic(seed) for seed in (11, 12, 13)]
    for item in diagnostics:
        item["lag_scan"]["estimates"]["diffusion_L_m2_s"][3] *= 4.0
        lineage = item["lineage_sampling"]
        lineage["minimum_effective_lineage_count"][3] = 127.0
        lineage["minimum_effective_lineage_fraction"][3] = 127.0 / 512.0

    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 1
    assert summary["solver_lag_supplemental_replicates"] == 0
    assert summary["solver_lag_supplemental_max_relative_ci95_bound"] is None


@pytest.mark.parametrize("invalid_value", [0.0, -1.0, float("nan")])
def test_weighted_mc_invalid_128_moment_rejects_when_lineage_is_qualified(
    invalid_value: float,
) -> None:
    diagnostics = [_weighted_mc_transport_diagnostic(seed) for seed in (11, 12, 13)]
    diagnostics[0]["lag_scan"]["estimates"]["diffusion_L_m2_s"][3] = invalid_value

    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 0
    assert summary["solver_lag_convergence_max_relative_ci95_bound"] == 0.0
    assert summary["solver_lag_supplemental_replicates"] == 2
    assert summary["solver_lag_supplemental_max_relative_ci95_bound"] is None
    assert (
        "mc_weighted_growth_supplemental_moment_evidence_missing"
        in summary["solver_failure_reasons"]
    )


@pytest.mark.parametrize("field", MC_TRANSPORT_LAG_FIELDS)
def test_weighted_mc_all_production_diffusion_fields_are_lag64(
    field: str,
) -> None:
    diagnostics = [_weighted_mc_transport_diagnostic(seed) for seed in (11, 12, 13)]
    diagnostics[0]["production_estimates"][field] *= 1.01
    diagnostics[0]["_case_reported_transport"][field] *= 1.01

    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 0
    assert (
        "mc_weighted_growth_ensemble_evidence_missing"
        in summary["solver_failure_reasons"]
    )


def test_weighted_mc_transport_qualification_rejects_growth_rate_mismatch() -> None:
    diagnostics = [_weighted_mc_transport_diagnostic(seed) for seed in (11, 12, 13)]
    for item in diagnostics:
        growth = item["population_growth_consistency"]
        growth["cumulative_log_growth"] = 2.0
        growth["population_log_growth_frequency_s_inv"] = 2.0e6

    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 0
    assert (
        "mc_weighted_growth_population_rate_inconsistent"
        in summary["solver_failure_reasons"]
    )


def test_mc_transport_quality_uses_paired_ensemble_evidence() -> None:
    summary = summarize_mc_transport_diagnostics(
        [_mc_transport_diagnostic(seed=index) for index in range(4)]
    )

    assert summary["solver_diagnostics_passed"] == 1
    assert summary["solver_converged"] is None
    assert summary["solver_transport_qualified"] == 1
    assert summary["solver_transport_replicates"] == 4
    assert summary[
        "solver_origin_stationarity_limiting_relative_ci95_bound"
    ] == pytest.approx(-math.log(1.0 - 0.01))
    assert summary["solver_origin_stationarity_limiting_field"] == ("mean_energy_eV")
    assert summary["solver_lag_convergence_max_relative_ci95_bound"] == pytest.approx(
        -math.log(1.0 - 0.02)
    )


def test_mc_transport_quality_rejects_unconverged_ensemble_lag() -> None:
    summary = summarize_mc_transport_diagnostics(
        [_mc_transport_diagnostic(seed=index, lag_change=0.3) for index in range(4)]
    )

    assert summary["solver_diagnostics_passed"] == 0
    assert any(
        reason.startswith("mc_transport_lag_plateau_not_converged:")
        for reason in summary["solver_failure_reasons"]
    )


def test_mc_transport_origin_stationarity_tracks_state_not_diffusion_noise() -> None:
    diagnostics = [_mc_transport_diagnostic(seed=index) for index in range(4)]
    for item in diagnostics:
        item["lag_scan"]["estimates"]["energy_diffusion_L_m2_s"][0] *= 0.25

    summary = summarize_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 1
    assert summary[
        "solver_origin_stationarity_limiting_relative_ci95_bound"
    ] == pytest.approx(-math.log(1.0 - 0.01))
    assert summary["solver_origin_stationarity_limiting_field"] == ("mean_energy_eV")
    assert summary["solver_origin_stationarity_limiting_relative_tolerance"] == 0.1
    assert summary["solver_lag_convergence_relative_tolerance"] == 0.25


def test_mc_transport_quality_fails_closed_without_lag_scan() -> None:
    legacy = {
        "energy_transport_status": "direct_mc_qualified",
        "finite": True,
        "strictly_positive": True,
        "multiple_time_origins": True,
    }
    summary = summarize_mc_transport_diagnostics([legacy] * 4)

    assert summary["solver_diagnostics_passed"] == 0
    assert "mc_transport_ensemble_evidence_missing" in summary["solver_failure_reasons"]


def test_mc_transport_quality_rejects_pre_v4_estimator_metadata() -> None:
    diagnostics = [_mc_transport_diagnostic(seed=index) for index in range(4)]
    diagnostics[-1]["estimator_schema_version"] = "direct_mc_transport.v3"

    summary = summarize_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 0
    assert (
        "mc_transport_estimator_provenance_mismatch"
        in summary["solver_failure_reasons"]
    )


def test_mc_transport_quality_rejects_report_not_backed_by_production() -> None:
    diagnostics = [_mc_transport_diagnostic(seed=index) for index in range(4)]
    diagnostics[-1]["_case_reported_transport"]["mobility_m2_V_s"] *= 1.01

    summary = summarize_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 0
    assert (
        "mc_transport_report_not_production_estimate"
        in summary["solver_failure_reasons"]
    )


def test_mc_transport_quality_rejects_component_estimator_substitution() -> None:
    diagnostics = [_mc_transport_diagnostic(seed=index) for index in range(4)]
    diagnostics[-1]["component_estimators"]["particle_diffusion"] = "einstein_relation"

    summary = summarize_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 0
    assert (
        "mc_transport_component_estimator_contract_mismatch"
        in summary["solver_failure_reasons"]
    )


@pytest.mark.parametrize(
    ("mutation", "reason"),
    [
        (
            lambda item: item["lag_scan"].__setitem__(
                "lag_s", [8.0e-9, 16.0e-9, 32.0e-9, 65.0e-9]
            ),
            "mc_transport_positive_ensemble_evidence_required",
        ),
        (
            lambda item: item["transport_sampling"].__setitem__(
                "observation_planes_complete", 127
            ),
            "mc_transport_positive_ensemble_evidence_required",
        ),
        (
            lambda item: item["origin_stationarity_window_mean_early"].__setitem__(
                "mean_energy_eV", -1.0
            ),
            "mc_transport_positive_ensemble_evidence_required",
        ),
        (
            lambda item: item["mc_run_provenance"].__setitem__("particles", 129),
            "mc_transport_replica_settings_mismatch",
        ),
        (
            lambda item: item["mc_run_provenance"].__setitem__(
                "tail_rate_rse_trigger", 0.1
            ),
            "mc_transport_replica_settings_mismatch",
        ),
    ],
)
def test_mc_transport_quality_fails_closed_on_mixed_or_invalid_evidence(
    mutation: object,
    reason: str,
) -> None:
    diagnostics = [_mc_transport_diagnostic(seed=index) for index in range(4)]
    mutation(diagnostics[-1])

    summary = summarize_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 0
    assert reason in summary["solver_failure_reasons"]


def test_mc_transport_run_provenance_allows_additive_execution_metadata() -> None:
    diagnostics = [_mc_transport_diagnostic(seed=index) for index in range(4)]
    diagnostics[-1]["mc_run_provenance"]["numeric_kernel_used"] = "numba"

    summary = summarize_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 1


def test_weighted_mc_tail_activation_may_differ_between_replicas() -> None:
    diagnostics = [_weighted_mc_transport_diagnostic(seed) for seed in (11, 12, 13)]
    provenance = diagnostics[-1]["mc_run_provenance"]
    provenance["tail_collisions_executed"] = provenance["tail_max_collisions"]

    summary = summarize_weighted_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 1


def test_mc_eedf_rate_consistency_uses_same_trajectory_evidence(
    tmp_path: Path,
) -> None:
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        for replicate in range(3):
            insert_workflow_case(
                store.connection,
                e_over_n=100.0,
                replicate=replicate,
            )
            insert_workflow_rate(
                store.connection,
                e_over_n=100.0,
                replicate=replicate,
                value=1.0e-15,
            )
        row = store.connection.execute(
            """
            SELECT diagnostics_json FROM cases
            WHERE solver = 'monte_carlo' AND e_over_n_Td = 100.0
              AND replicate = 0
            """
        ).fetchone()
        diagnostics = json.loads(str(row[0]))
        diagnostics["internal_monte_carlo_reaction_rates"]["rates"][0][
            "histogram_convolution_rate_coefficient_m3_s"
        ] = 0.5e-15
        store.connection.execute(
            """
            UPDATE cases SET diagnostics_json = ?
            WHERE solver = 'monte_carlo' AND e_over_n_Td = 100.0
              AND replicate = 0
            """,
            (json.dumps(diagnostics),),
        )
        store.connection.row_factory = sqlite3.Row
        failures = _mc_eedf_rate_consistency_failures(
            store.connection,
            0,
            allowed_e={100.0},
        )
        diagnostics["internal_monte_carlo_reaction_rates"]["rates"][0][
            "rate_coefficient_m3_s"
        ] = 0.0
        diagnostics["internal_monte_carlo_reaction_rates"]["rates"][0][
            "histogram_convolution_rate_coefficient_m3_s"
        ] = 0.0
        store.connection.execute(
            """
            UPDATE cases SET diagnostics_json = ?
            WHERE solver = 'monte_carlo' AND e_over_n_Td = 100.0
              AND replicate = 0
            """,
            (json.dumps(diagnostics),),
        )
        unresolved_failures = _mc_eedf_rate_consistency_failures(
            store.connection,
            0,
            allowed_e={100.0},
        )

    assert failures == {100.0: ("mc_eedf_rate_consistency_failed",)}
    assert unresolved_failures == {}


def test_mc_transport_quality_rejects_duplicate_seed() -> None:
    diagnostics = [_mc_transport_diagnostic(seed=index) for index in range(4)]
    diagnostics[-1]["mc_run_provenance"]["seed"] = 0

    summary = summarize_mc_transport_diagnostics(diagnostics)

    assert summary["solver_transport_qualified"] == 0
    assert (
        "mc_transport_replica_seeds_missing_or_not_distinct"
        in summary["solver_failure_reasons"]
    )


def _probability_masses_for_mean(
    mean_energy_eV: float,
    *,
    bin_count: int,
) -> tuple[float, ...]:
    centers = [0.5 + index for index in range(bin_count)]
    if not centers[0] <= mean_energy_eV <= centers[-1]:
        raise ValueError("mean energy is outside the test EEDF grid")
    masses = [0.0] * bin_count
    if mean_energy_eV == centers[-1]:
        masses[-1] = 1.0
        return tuple(masses)
    lower = max(
        index for index, center in enumerate(centers) if center <= mean_energy_eV
    )
    upper = min(lower + 1, bin_count - 1)
    if lower == upper:
        masses[lower] = 1.0
    else:
        upper_fraction = (mean_energy_eV - centers[lower]) / (
            centers[upper] - centers[lower]
        )
        masses[lower] = 1.0 - upper_fraction
        masses[upper] = upper_fraction
    return tuple(masses)


def _probability_masses_with_tail(
    mean_energy_eV: float,
    *,
    tail_mass: float = 0.05,
) -> tuple[float, ...]:
    tail_energy = 2.5
    bulk_mean = (mean_energy_eV - tail_mass * tail_energy) / (1.0 - tail_mass)
    bulk = _probability_masses_for_mean(bulk_mean, bin_count=2)
    return (
        (1.0 - tail_mass) * bulk[0],
        (1.0 - tail_mass) * bulk[1],
        tail_mass,
    )


def _insert_case_set(
    store: WorkflowStore,
    *,
    solver: str,
    mixture_id: int = 0,
    e_values: tuple[float, ...] = (10.0, 20.0),
    replicate: int = 0,
    mean_values: tuple[float, ...] = (1.0, 2.0),
    mobility_values: tuple[float, ...] = (10.0, 11.0),
) -> None:
    for index, e_over_n in enumerate(e_values):
        insert_workflow_case(
            store.connection,
            mixture_id=mixture_id,
            solver=solver,
            e_over_n=e_over_n,
            replicate=replicate,
            mean_energy=mean_values[index],
            drift=5.0,
            reduced_mobility=mobility_values[index],
            reduced_diffusion_l=20.0 + index,
            reduced_diffusion_t=30.0 + index,
            effective_townsend=1.0e-16,
        )
        insert_workflow_rate(
            store.connection,
            mixture_id=mixture_id,
            solver=solver,
            e_over_n=e_over_n,
            replicate=replicate,
            value=1.0e-15 * (index + 1),
            target_fraction=0.5,
            process="ionization",
            process_type="ionization",
            threshold_eV=15.0,
        )
        insert_workflow_eedf(
            store.connection,
            mixture_id=mixture_id,
            solver=solver,
            e_over_n=e_over_n,
            replicate=replicate,
        )


def _two_term_elastic_loss_diagnostic(value: float) -> dict[str, object]:
    return {
        "two_term": {
            "elastic_energy_loss": {
                "schema": "swarm.elastic_energy_loss.v1",
                "status": "available",
                "symbol": "K_epsilon_el",
                "rate_coefficient_eV_m3_s": value,
                "estimator": ("same_discrete_elastic_collision_operator_energy_moment"),
                "operator": (
                    "native_finite_volume_scharfetter_gummel_elastic_A_D_zero_field"
                ),
                "operator_isolation": (
                    "zero_field_reassembly_from_same_elastic_A_D_coefficients"
                ),
                "source_eedf": "same_solved_eedf",
                "gas_temperature_K": 300.0,
                "neutral_thermal_motion_model": ("finite_temperature_fokker_planck"),
                "gas_temperature_terms_included": True,
                "sign_convention": "positive_is_net_electron_energy_loss",
                "uncertainty_status": "deterministic_kinetic_solver",
            }
        }
    }


def _two_term_temporal_growth_diagnostic(
    growth_frequency_s_inv: float,
    *,
    applied: bool = True,
) -> dict[str, object]:
    return {
        "two_term": {
            "converged": True,
            "iterations": 2,
            "residual_L1": 1.0e-8,
            "residual_tolerance": 1.0e-6,
            "tail_probability": 1.0e-12,
            "tail_probability_target": 1.0e-9,
            "edge_to_peak": 1.0e-12,
            "edge_to_peak_target": 1.0e-10,
            "grid_max_eV": 100.0,
            "grid_max_limit_eV": 1000.0,
            "growth_frequency_s-1": growth_frequency_s_inv,
            "temporal_growth_momentum_correction": {
                "applied": applied,
                "model": ("nu_m_eff(epsilon)=nu_m(epsilon)+growth_frequency"),
                "growth_frequency_s-1": growth_frequency_s_inv,
                "frequency_source": "converged_temporal_growth_eigenvalue",
            },
        }
    }


def test_two_term_temporal_growth_transport_provenance_is_exported(
    tmp_path: Path,
) -> None:
    database = tmp_path / "workflow.sqlite"
    with WorkflowStore(database) as store:
        _insert_case_set(store, solver="two_term")
        for field, growth in ((10.0, 2.0e5), (20.0, 4.0e5)):
            store.connection.execute(
                """
                UPDATE cases SET diagnostics_json = ?
                WHERE solver = 'two_term' AND e_over_n_Td = ?
                """,
                (
                    json.dumps(_two_term_temporal_growth_diagnostic(growth)),
                    field,
                ),
            )
        store.connection.commit()

    build_tables(database, tmp_path / "tables", source="two_term")
    mixture = tmp_path / "tables" / "mixture_0000"
    rows = _read_csv(mixture / "transport_vs_mean_energy.csv")
    assert (
        tuple(
            name
            for name in TWO_TERM_TEMPORAL_GROWTH_TRANSPORT_COLUMNS
            if name in rows[0]
        )
        == TWO_TERM_TEMPORAL_GROWTH_TRANSPORT_COLUMNS
    )
    assert [
        float(row["temporal_growth_frequency_s_inv"]) for row in rows
    ] == pytest.approx([2.0e5, 4.0e5])
    assert [
        float(row["reduced_temporal_growth_frequency_m3_s"]) for row in rows
    ] == pytest.approx([2.0e-17, 4.0e-17])
    manifest = json.loads((mixture / "manifest.json").read_text())
    contract = manifest["source_policy"]["two_term_transport_kernel"]
    assert contract["schema"] == "two_term_temporal_growth_transport.v1"
    assert contract["correction_applied_every_anchor"] is True
    assert contract["momentum_cross_section_floor_policy"] == (
        "sum_raw_process_cross_sections_then_single_floor"
    )


def test_two_term_temporal_growth_transport_provenance_is_fail_closed(
    tmp_path: Path,
) -> None:
    database = tmp_path / "workflow.sqlite"
    with WorkflowStore(database) as store:
        _insert_case_set(store, solver="two_term")
        for field, applied in ((10.0, True), (20.0, False)):
            store.connection.execute(
                """
                UPDATE cases SET diagnostics_json = ?
                WHERE solver = 'two_term' AND e_over_n_Td = ?
                """,
                (
                    json.dumps(
                        _two_term_temporal_growth_diagnostic(
                            2.0e5,
                            applied=applied,
                        )
                    ),
                    field,
                ),
            )
        store.connection.commit()

    with pytest.raises(TableBuildError, match="provenance is partial"):
        build_tables(database, tmp_path / "tables", source="two_term")


def _mc_elastic_loss_diagnostic(
    value: float,
    *,
    seed: int = 1,
) -> dict[str, object]:
    return {
        "internal_monte_carlo_transport": {
            "mc_run_provenance": {"seed": seed, "gas_temperature_K": 300.0}
        },
        "internal_monte_carlo_reaction_rates": {
            "rates": [
                {
                    "species": "Ar",
                    "process": "elastic",
                    "process_type": "elastic",
                    "target_species_fraction": 0.5,
                    "event_sampling_enabled": True,
                    "energy_loss": {
                        "status": "direct_event_estimator",
                        "estimator": (
                            "trajectory_event_energy_change_per_target_density_"
                            "residence_time"
                        ),
                        "model": ("maxwellian_target_exact_binary_collision_isotropic"),
                        "rate_coefficient_eV_m3_s": value,
                        "neutral_thermal_motion_model": (
                            "maxwellian_relative_speed_exact_binary_collision"
                        ),
                        "gas_temperature_terms_included": True,
                    },
                }
            ]
        },
    }


def test_two_term_elastic_energy_loss_table_has_operator_contract(
    tmp_path: Path,
) -> None:
    database = tmp_path / "workflow.sqlite"
    with WorkflowStore(database) as store:
        _insert_case_set(store, solver="two_term")
        for e_over_n, value in ((10.0, 1.0e-18), (20.0, 2.0e-18)):
            store.connection.execute(
                """
                UPDATE cases SET diagnostics_json = ?
                WHERE solver = 'two_term' AND e_over_n_Td = ?
                """,
                (json.dumps(_two_term_elastic_loss_diagnostic(value)), e_over_n),
            )
        store.connection.commit()

    build_tables(database, tmp_path / "tables", source="two_term")

    mixture = tmp_path / "tables" / "mixture_0000"
    rows = _read_csv(mixture / ELASTIC_ENERGY_LOSS_TABLE)
    assert [float(row["mean_energy_eV"]) for row in rows] == [1.0, 2.0]
    assert [
        float(row["elastic_energy_loss_rate_coefficient_eV_m3_s"]) for row in rows
    ] == pytest.approx([1.0e-18, 2.0e-18])
    manifest = json.loads((mixture / "manifest.json").read_text())
    table = manifest["tables"][ELASTIC_ENERGY_LOSS_TABLE]
    assert table["artifact_role"] == ("canonical_comsol_elastic_energy_loss_input")
    assert table["physics_contract"]["gas_temperature_terms_included"] is True


def test_mc_elastic_energy_loss_aggregates_direct_replicas_with_uncertainty(
    tmp_path: Path,
) -> None:
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        for replicate, value in enumerate((1.0e-18, 2.0e-18, 3.0e-18, 4.0e-18)):
            insert_workflow_case(
                store.connection,
                e_over_n=10.0,
                replicate=replicate,
                mean_energy=2.0,
                diagnostics=_mc_elastic_loss_diagnostic(value, seed=100 + replicate),
            )
        store.connection.commit()
        store.connection.row_factory = sqlite3.Row
        rows, contract = _load_elastic_energy_loss(
            store.connection,
            "monte_carlo",
            0,
            [{"E_over_N_Td": 10.0, "mean_energy_eV": 2.0}],
            allowed_e={10.0},
        )

    assert rows is not None and contract is not None
    [row] = rows
    assert row["elastic_energy_loss_rate_coefficient_eV_m3_s"] == pytest.approx(
        1.25e-18
    )
    assert row["elastic_energy_loss_standard_error_eV_m3_s"] is not None
    assert row["uncertainty_available"] == 1
    assert contract["neutral_thermal_motion_model"] == (
        "maxwellian_relative_speed_exact_binary_collision"
    )
    assert contract["gas_temperature_terms_included"] is True
    assert contract["gas_temperature_K"] == pytest.approx(300.0)


def test_elastic_energy_loss_fails_closed_on_partial_replica_artifact(
    tmp_path: Path,
) -> None:
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        insert_workflow_case(
            store.connection,
            solver="two_term",
            e_over_n=10.0,
            replicate=0,
            diagnostics=_two_term_elastic_loss_diagnostic(1.0e-18),
        )
        insert_workflow_case(
            store.connection,
            solver="two_term",
            e_over_n=10.0,
            replicate=1,
            diagnostics={},
        )
        store.connection.commit()
        store.connection.row_factory = sqlite3.Row
        with pytest.raises(TableBuildError, match="only part"):
            _load_elastic_energy_loss(
                store.connection,
                "two_term",
                0,
                [{"E_over_N_Td": 10.0, "mean_energy_eV": 2.0}],
                allowed_e={10.0},
            )


def test_build_tables_two_term_writes_two_mixture_directories(tmp_path: Path) -> None:
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        _insert_case_set(store, solver="two_term", mixture_id=0)
        _insert_case_set(store, solver="two_term", mixture_id=1)
        store.connection.commit()

    summary = build_tables(
        tmp_path / "workflow.sqlite",
        tmp_path / "tables",
        source="two_term",
    )

    assert summary.mixtures == 2
    root_manifest = json.loads((tmp_path / "tables" / "manifest.json").read_text())
    assert [item["mixture_id"] for item in root_manifest["mixtures"]] == [0, 1]
    for mixture_id in [0, 1]:
        mixture_dir = tmp_path / "tables" / f"mixture_{mixture_id:04d}"
        rows = _read_csv(mixture_dir / "transport_vs_en.csv")
        assert len(rows) == 2
        assert float(rows[0]["E_over_N_V_m2"]) == pytest.approx(
            float(rows[0]["E_over_N_Td"]) * 1.0e-21
        )
        for table_name in ["transport_vs_en.csv", "rates_vs_en.csv"]:
            for row in _read_csv(mixture_dir / table_name):
                for value in row.values():
                    if value:
                        number = float(value) if _looks_numeric(value) else None
                        if number is not None:
                            assert math.isfinite(number)
        manifest = json.loads((mixture_dir / "manifest.json").read_text())
        assert manifest["table_argument"]["secondary"] == "mean_energy_eV"
        assert (
            manifest["source_policy"]["transport_definition"] == "fixed_particle_flux"
        )


def test_rebuild_removes_stale_mean_energy_tables_when_nonmonotonic(
    tmp_path: Path,
) -> None:
    database = tmp_path / "workflow.sqlite"
    output = tmp_path / "tables"
    with WorkflowStore(database) as store:
        _insert_case_set(store, solver="two_term")
        store.connection.commit()

    build_tables(database, output, source="two_term")
    mixture_dir = output / "mixture_0000"
    assert (mixture_dir / "transport_vs_mean_energy.csv").exists()
    assert (mixture_dir / "rates_vs_mean_energy.csv").exists()
    (mixture_dir / "unlisted_stale.csv").write_text("stale\n", encoding="utf-8")

    with sqlite3.connect(database) as connection:
        connection.execute(
            "UPDATE cases SET mean_energy_eV = 1.0 WHERE e_over_n_Td = 20.0"
        )
        connection.commit()
        aggregate_workflow_results(connection)

    build_tables(database, output, source="two_term")

    assert not (mixture_dir / "transport_vs_mean_energy.csv").exists()
    assert not (mixture_dir / "rates_vs_mean_energy.csv").exists()
    assert not (mixture_dir / "unlisted_stale.csv").exists()
    manifest = json.loads((mixture_dir / "manifest.json").read_text())
    assert manifest["table_argument"]["secondary"] is None
    assert manifest["monotonicity"]["mean_energy_strictly_monotonic"] is False


def test_build_tables_refreshes_partial_aggregate_after_new_case_commit(
    tmp_path: Path,
) -> None:
    database = tmp_path / "workflow.sqlite"
    output = tmp_path / "tables"
    with WorkflowStore(database) as store:
        _insert_case_set(
            store,
            solver="two_term",
            e_values=(10.0,),
            mean_values=(1.0,),
            mobility_values=(10.0,),
        )
        store.connection.commit()
    build_tables(database, output, source="two_term")

    with WorkflowStore(database) as store:
        _insert_case_set(
            store,
            solver="two_term",
            e_values=(20.0,),
            mean_values=(2.0,),
            mobility_values=(11.0,),
        )
        store.connection.commit()

    build_tables(database, output, source="two_term")

    rows = _read_csv(output / "mixture_0000" / "transport_vs_en.csv")
    assert [float(row["E_over_N_Td"]) for row in rows] == [10.0, 20.0]


def test_build_tables_rejects_aggregate_rows_with_different_quality_policy(
    tmp_path: Path,
) -> None:
    database = tmp_path / "workflow.sqlite"
    with WorkflowStore(database) as store:
        _insert_case_set(store, solver="two_term")
        store.connection.commit()
    build_tables(database, tmp_path / "tables", source="two_term")

    with sqlite3.connect(database) as connection:
        connection.execute(
            "UPDATE aggregate_quality SET thresholds_json = ?",
            (quality_thresholds_json(QualityThresholds(mobility_rse=0.1)),),
        )
        connection.commit()

    with pytest.raises(
        WorkflowSchemaError,
        match="aggregate_quality source/evaluation policy provenance",
    ):
        build_tables(database, tmp_path / "tables", source="two_term")


def test_two_term_quality_includes_and_enforces_solver_diagnostics(
    tmp_path: Path,
) -> None:
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        _insert_case_set(
            store,
            solver="two_term",
            e_values=(1.0,),
            mean_values=(1.0,),
            mobility_values=(10.0,),
        )
        diagnostics = {
            "two_term": {
                "converged": False,
                "iterations": 600,
                "residual_L1": 2.0e-7,
                "residual_tolerance": 1.0e-7,
                "tail_probability": 1.0e-20,
                "tail_probability_target": 1.0e-9,
                "edge_to_peak": 1.0e-30,
                "edge_to_peak_target": 1.0e-10,
                "grid_max_eV": 100.0,
                "grid_max_limit_eV": 20000.0,
            }
        }
        store.connection.execute(
            """
            UPDATE cases
            SET diagnostics_json = ?
            WHERE solver = 'two_term' AND e_over_n_Td = 1.0
            """,
            (json.dumps(diagnostics),),
        )
        store.connection.commit()

    build_tables(
        tmp_path / "workflow.sqlite",
        tmp_path / "tables",
        source="two_term",
    )
    [quality] = _read_csv(tmp_path / "tables" / "mixture_0000" / "quality.csv")
    assert quality["passed"] == "0"
    assert quality["solver_diagnostics_available"] == "1"
    assert quality["solver_converged"] == "0"
    assert quality["solver_diagnostics_passed"] == "0"
    assert "mobility_rse" not in quality
    assert "qualification_profile" not in quality
    assert "solver_transport_qualified" not in quality
    assert "solver_not_converged" in json.loads(quality["failure_reasons_json"])
    manifest = json.loads(
        (tmp_path / "tables" / "mixture_0000" / "manifest.json").read_text()
    )
    assert manifest["quality_summary"]["passed"] is False
    quality_entry = manifest["tables"]["quality.csv"]
    deterministic_schema = quality_table_schema("two_term")
    assert quality_entry["columns"] == list(deterministic_schema.columns)
    assert quality_entry["common_columns"] == list(COMMON_QUALITY_COLUMNS)
    assert quality_entry["solver_evidence"] == {
        "kind": "deterministic_convergence",
        "columns": list(deterministic_schema.evidence_columns),
    }


@pytest.mark.parametrize("removed_source", ["monte_carlo_smoothed", "hybrid"])
def test_build_tables_rejects_removed_noncanonical_sources(
    tmp_path: Path,
    removed_source: str,
) -> None:
    with pytest.raises(TableBuildError, match="unsupported table source"):
        build_tables(
            tmp_path / "workflow.sqlite",
            tmp_path / "tables",
            source=removed_source,
        )


def test_build_tables_nonmonotonic_mean_energy_skips_mean_energy_tables(
    tmp_path: Path,
) -> None:
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        _insert_case_set(
            store,
            solver="two_term",
            e_values=(10.0, 20.0, 30.0),
            mean_values=(1.0, 3.0, 2.0),
            mobility_values=(10.0, 11.0, 12.0),
        )
        store.connection.commit()

    build_tables(tmp_path / "workflow.sqlite", tmp_path / "tables", source="two_term")

    mixture_dir = tmp_path / "tables" / "mixture_0000"
    manifest = json.loads((mixture_dir / "manifest.json").read_text())
    assert manifest["table_argument"]["secondary"] is None
    assert manifest["monotonicity"]["reason"] == "mean_energy_not_strictly_monotonic"
    assert not (mixture_dir / "transport_vs_mean_energy.csv").exists()


def test_build_tables_monte_carlo_requires_four_qualified_anchors(
    tmp_path: Path,
) -> None:
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        _insert_case_set(store, solver="monte_carlo", e_values=(10.0,))
        store.connection.commit()

    with pytest.raises(TableBuildError, match="at least four aggregate-"):
        build_tables(
            tmp_path / "workflow.sqlite",
            tmp_path / "tables",
            source="monte_carlo",
        )


def test_unqualified_mc_evidence_survives_without_interpolation_tables(
    tmp_path: Path,
) -> None:
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        _insert_case_set(store, solver="monte_carlo", e_values=(10.0,))
        store.connection.commit()
    build_tables(
        tmp_path / "workflow.sqlite",
        tmp_path / "tables",
        source="monte_carlo",
        allow_unqualified_mc=True,
    )
    directory = tmp_path / "tables/mixture_0000"
    manifest = json.loads((directory / "manifest.json").read_text())
    database_digest = sha256((tmp_path / "workflow.sqlite").read_bytes()).hexdigest()
    assert manifest["hashes"]["source_database_sha256"] == database_digest
    root_manifest = json.loads((tmp_path / "tables/manifest.json").read_text())
    assert root_manifest["hashes"]["source_database_sha256"] == database_digest
    assert manifest["status"] == "unqualified"
    assert manifest["mc_qualification"]["all_planned_anchors"] == 1
    assert manifest["mc_qualification"]["coefficient_tables_available"] is False
    assert (directory / "mc_qualification.csv").is_file()
    assert not (directory / "transport_vs_mean_energy.csv").exists()
    from swarm_workflow.comsol.input import ComsolExportError, export_comsol_bundle

    with pytest.raises(ComsolExportError, match="no exportable coefficient"):
        export_comsol_bundle(directory, tmp_path / "bundle")


def test_mc_table_rejects_legacy_v2_database_without_sampling_plan(
    tmp_path: Path,
) -> None:
    database = tmp_path / "workflow.sqlite"
    with WorkflowStore(database) as store:
        insert_workflow_case(
            store.connection,
            solver="monte_carlo",
            e_over_n=10.0,
        )
        store.connection.execute(
            "DELETE FROM metadata WHERE key = ?",
            (MC_SAMPLING_PLAN_METADATA_KEY,),
        )
        store.connection.commit()

    with pytest.raises(WorkflowSchemaError, match="missing mc_sampling_plan_json"):
        build_tables(database, tmp_path / "tables", source="monte_carlo")


def test_mc_table_rejects_sampling_plan_that_disagrees_with_run_provenance(
    tmp_path: Path,
) -> None:
    database = tmp_path / "workflow.sqlite"
    with WorkflowStore(database) as store:
        insert_workflow_case(
            store.connection,
            solver="monte_carlo",
            e_over_n=10.0,
        )
        raw = store.connection.execute(
            "SELECT value FROM metadata WHERE key = ?",
            (MC_SAMPLING_PLAN_METADATA_KEY,),
        ).fetchone()[0]
        plan = json.loads(str(raw))
        plan[0]["transport_estimator"] = "paired_field_parity"
        store.connection.execute(
            "UPDATE metadata SET value = ? WHERE key = ?",
            (
                canonical_mc_sampling_plan_json(plan),
                MC_SAMPLING_PLAN_METADATA_KEY,
            ),
        )
        store.connection.commit()

    with pytest.raises(WorkflowSchemaError, match="controls do not match"):
        build_tables(database, tmp_path / "tables", source="monte_carlo")


def test_v6_sampling_plan_validates_anchor_specific_reported_lags(
    tmp_path: Path,
) -> None:
    database = tmp_path / "workflow.sqlite"
    with WorkflowStore(database) as store:
        for e_over_n, lag in ((10.0, 64), (20.0, 256)):
            insert_workflow_case(
                store.connection,
                solver="monte_carlo",
                e_over_n=e_over_n,
                diagnostics={
                    "internal_monte_carlo_transport": (
                        _weighted_mc_transport_diagnostic(
                            int(e_over_n),
                            correlation_lag=lag,
                        )
                    )
                },
            )
        store.connection.execute(
            "UPDATE metadata SET value = ? WHERE key = ?",
            (
                WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
                "mc_transport_estimator_schema_version",
            ),
        )
        store.connection.commit()
        metadata = store.metadata()
        store.connection.row_factory = sqlite3.Row

        sampling = _validate_mc_sampling_plan_against_cases(
            store.connection,
            metadata,
        )

    assert [
        entry["transport_correlation_lag_barriers"] for entry in sampling["entries"]
    ] == [64, 256]


def test_v6_sampling_plan_rejects_case_reported_lag_mismatch(
    tmp_path: Path,
) -> None:
    database = tmp_path / "workflow.sqlite"
    with WorkflowStore(database) as store:
        insert_workflow_case(
            store.connection,
            solver="monte_carlo",
            e_over_n=10.0,
            diagnostics={
                "internal_monte_carlo_transport": (
                    _weighted_mc_transport_diagnostic(10, correlation_lag=256)
                )
            },
        )
        raw = store.connection.execute(
            "SELECT value FROM metadata WHERE key = ?",
            (MC_SAMPLING_PLAN_METADATA_KEY,),
        ).fetchone()[0]
        plan = json.loads(str(raw))
        plan[0]["transport_correlation_lag_barriers"] = 64
        store.connection.execute(
            "UPDATE metadata SET value = ? WHERE key = ?",
            (
                canonical_mc_sampling_plan_json(plan),
                MC_SAMPLING_PLAN_METADATA_KEY,
            ),
        )
        store.connection.execute(
            "UPDATE metadata SET value = ? WHERE key = ?",
            (
                WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
                "mc_transport_estimator_schema_version",
            ),
        )
        store.connection.commit()
        store.connection.row_factory = sqlite3.Row

        with pytest.raises(WorkflowSchemaError, match="controls do not match"):
            _validate_mc_sampling_plan_against_cases(
                store.connection,
                store.metadata(),
            )


def test_mc_table_fails_closed_when_raw_aggregate_passes_but_diagnostic_missing(
    tmp_path: Path,
) -> None:
    database = tmp_path / "workflow.sqlite"
    with WorkflowStore(database) as store:
        for e_over_n in (10.0, 20.0, 30.0, 40.0):
            for replicate in range(3):
                insert_workflow_case(
                    store.connection,
                    solver="monte_carlo",
                    e_over_n=e_over_n,
                    replicate=replicate,
                    mean_energy=1.0 + e_over_n / 100.0,
                )
                insert_workflow_rate(
                    store.connection,
                    solver="monte_carlo",
                    e_over_n=e_over_n,
                    replicate=replicate,
                )
                insert_workflow_eedf(
                    store.connection,
                    solver="monte_carlo",
                    e_over_n=e_over_n,
                    replicate=replicate,
                )
        aggregate_workflow_results(store.connection)
        assert (
            int(
                store.connection.execute(
                    """
            SELECT passed FROM aggregate_quality
            WHERE solver = 'monte_carlo' AND e_over_n_Td = 10.0
            """
                ).fetchone()[0]
            )
            == 1
        )
        diagnostic_row = store.connection.execute(
            """
            SELECT diagnostics_json FROM cases
            WHERE solver = 'monte_carlo' AND e_over_n_Td = 10.0
              AND replicate = 0
            """
        ).fetchone()
        diagnostics = json.loads(str(diagnostic_row[0]))
        del diagnostics["internal_monte_carlo_transport"][
            "origin_stationarity_window_mean_early"
        ]
        store.connection.execute(
            """
            UPDATE cases SET diagnostics_json = ?
            WHERE solver = 'monte_carlo' AND e_over_n_Td = 10.0
              AND replicate = 0
            """,
            (json.dumps(diagnostics),),
        )
        store.connection.commit()

    with pytest.raises(TableBuildError, match="at least four aggregate-"):
        build_tables(database, tmp_path / "tables", source="monte_carlo")


def test_build_tables_monte_carlo_uses_only_qualified_mc_evidence(
    tmp_path: Path,
) -> None:
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        for index, e_over_n in enumerate((75.0, 100.0, 150.0, 200.0, 300.0)):
            decoy_two_term_mean = 0.8 + 0.15 * index
            insert_workflow_case(
                store.connection,
                solver="two_term",
                e_over_n=e_over_n,
                mean_energy=decoy_two_term_mean,
                reduced_mobility=10.0 + index,
            )
            insert_workflow_rate(
                store.connection,
                solver="two_term",
                e_over_n=e_over_n,
                value=(index + 1) * 0.8e-16,
            )
            insert_workflow_eedf(
                store.connection,
                solver="two_term",
                e_over_n=e_over_n,
                probability_masses=_probability_masses_with_tail(decoy_two_term_mean),
            )
            target_mean = 0.9 + 0.1 * index
            for replicate, scale in enumerate((0.99, 1.0, 1.01)):
                mc_probability_masses = (
                    _probability_masses_with_tail(
                        target_mean,
                        tail_mass=(0.0, 0.04, 0.08)[replicate],
                    )
                    if index == 1
                    else _probability_masses_for_mean(
                        target_mean,
                        bin_count=3,
                    )
                )
                insert_workflow_case(
                    store.connection,
                    solver="monte_carlo",
                    e_over_n=e_over_n,
                    replicate=replicate,
                    mean_energy=target_mean,
                    reduced_mobility=(12.0 - index) * scale,
                    reduced_diffusion_l=(20.0 + index) * scale,
                    reduced_diffusion_t=(30.0 + index) * scale,
                    reduced_energy_mobility=(40.0 + index) * scale,
                    reduced_energy_diffusion=None,
                    reduced_energy_diffusion_l=(50.0 + index) * scale,
                    reduced_energy_diffusion_t=(60.0 + index) * scale,
                )
                insert_workflow_rate(
                    store.connection,
                    solver="monte_carlo",
                    e_over_n=e_over_n,
                    replicate=replicate,
                    value=(index + 1) * 1.0e-16 * scale,
                )
                insert_workflow_eedf(
                    store.connection,
                    solver="monte_carlo",
                    e_over_n=e_over_n,
                    replicate=replicate,
                    probability_masses=mc_probability_masses,
                )
        store.connection.commit()

    build_tables(
        tmp_path / "workflow.sqlite",
        tmp_path / "tables",
        source="monte_carlo",
    )

    mixture_dir = tmp_path / "tables" / "mixture_0000"
    transport = _read_csv(mixture_dir / "transport_vs_mean_energy.csv")
    rates = _read_csv(mixture_dir / "rates_vs_mean_energy.csv")
    assert all(float(row["reduced_mobility_m2_V_s_m3"]) > 0.0 for row in transport)
    direct_energy_fields = (
        "reduced_electron_energy_mobility_m2_V_s_m3",
        "reduced_electron_energy_diffusion_L_m2_s_m3",
        "reduced_electron_energy_diffusion_T_m2_s_m3",
    )
    assert all(
        float(row[field]) > 0.0 for row in transport for field in direct_energy_fields
    )
    assert all(
        row["reduced_electron_energy_diffusion_m2_s_m3"] == "" for row in transport
    )
    assert all(float(row["rate_coefficient_m3_s"]) > 0.0 for row in rates)
    assert all(
        right > left
        for left, right in zip(
            [float(row["mean_energy_eV"]) for row in transport],
            [float(row["mean_energy_eV"]) for row in transport][1:],
        )
    )
    manifest = json.loads((mixture_dir / "manifest.json").read_text())
    policy = manifest["source_policy"]
    assert manifest["source"] == "monte_carlo"
    sampling = manifest["mc_sampling_plan"]
    assert sampling["status"] == "validated_against_mc_cases"
    assert sampling["mixtures"] == 1
    assert len(sampling["entries"]) == 5
    sampling_json = manifest["hashes"][MC_SAMPLING_PLAN_METADATA_KEY]
    assert json.loads(sampling_json) == sampling["entries"]
    assert (
        manifest["hashes"]["mc_sampling_plan_sha256"]
        == sha256(sampling_json.encode("utf-8")).hexdigest()
    )
    root_manifest = json.loads((tmp_path / "tables" / "manifest.json").read_text())
    assert root_manifest["mc_sampling_plan"] == sampling
    assert policy["source"] == "monte_carlo"
    assert policy["postprocess"] == "none"
    assert policy["transport_eligible_anchor_points"] == 5
    assert policy["eedf_rate_eligible_anchor_points"] == 5
    assert policy["smoothing"] == {
        "transport_and_rates": "none",
        "eedf": (
            "shape_preserving_local_reconstruction_only; no_external_distribution_prior"
        ),
    }
    assert policy["transport_scope"] == {
        "standard_energy_density_gradient_closure": (
            "restricted_direct_density_packet_energy_flux_correlation"
        ),
        "full_cross_gradient_response_identified": False,
    }
    assert policy["reaction_rate_policy"] == (
        "direct_mc_trajectory_rate_coefficients_as_active_input;_"
        "raw_rate_evidence_as_qualification"
    )
    assert manifest["tables"]["eedf.csv"]["artifact_role"] == (
        "raw_swarm_eedf_audit_and_canonicalization_source"
    )
    quality_rows = _read_csv(mixture_dir / "quality.csv")
    assert quality_rows
    assert "solver_converged" not in quality_rows[0]
    assert "solver_residual_L1_max" not in quality_rows[0]
    assert "qualification_profile" in quality_rows[0]
    assert "solver_transport_qualified" in quality_rows[0]
    mc_schema = quality_table_schema("monte_carlo")
    quality_entry = manifest["tables"]["quality.csv"]
    assert quality_entry["columns"] == list(mc_schema.columns)
    assert quality_entry["solver_evidence"] == {
        "kind": "independent_replica_statistics",
        "columns": list(mc_schema.evidence_columns),
    }


def test_build_tables_monte_carlo_rejects_censored_tail_without_substitution(
    tmp_path: Path,
) -> None:
    e_values = (20.0, 35.0, 50.0, 75.0, 100.0)
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        for index, e_over_n in enumerate(e_values):
            mean_energy = 2.0 + index
            decoy_two_term_mean = 1.5 + 1.25 * index
            insert_workflow_case(
                store.connection,
                solver="two_term",
                e_over_n=e_over_n,
                mean_energy=decoy_two_term_mean,
                reduced_mobility=10.0 + index,
            )
            insert_workflow_rate(
                store.connection,
                solver="two_term",
                e_over_n=e_over_n,
                value=1.0e-30 * 100.0**index,
                process="ionization",
                process_type="ionization",
                threshold_eV=15.0,
            )
            insert_workflow_eedf(
                store.connection,
                solver="two_term",
                e_over_n=e_over_n,
                probability_masses=_probability_masses_for_mean(
                    decoy_two_term_mean,
                    bin_count=8,
                ),
            )
            for replicate, scale in enumerate((0.99, 1.0, 1.01)):
                insert_workflow_case(
                    store.connection,
                    solver="monte_carlo",
                    e_over_n=e_over_n,
                    replicate=replicate,
                    mean_energy=mean_energy * scale,
                    reduced_mobility=(12.0 + index) * scale,
                    reduced_diffusion_l=(20.0 + index) * scale,
                    reduced_diffusion_t=(30.0 + index) * scale,
                    diagnostics=_mc_rate_diagnostics(
                        event_count=0 if index < 2 else 1,
                        residence_time_s=(replicate + 1) * 1.0e-6,
                    ),
                )
                insert_workflow_rate(
                    store.connection,
                    solver="monte_carlo",
                    e_over_n=e_over_n,
                    replicate=replicate,
                    value=(0.0 if index < 2 else index * 1.0e-18 * scale),
                    process="ionization",
                    process_type="ionization",
                    threshold_eV=15.0,
                )
                insert_workflow_eedf(
                    store.connection,
                    solver="monte_carlo",
                    e_over_n=e_over_n,
                    replicate=replicate,
                    probability_masses=_probability_masses_for_mean(
                        mean_energy * scale,
                        bin_count=8,
                    ),
                )
        store.connection.commit()

    with pytest.raises(TableBuildError, match="instead of substituting"):
        build_tables(
            tmp_path / "workflow.sqlite",
            tmp_path / "tables",
            source="monte_carlo",
        )


def test_build_tables_monte_carlo_preserves_evidence_when_qualification_fails(
    tmp_path: Path,
) -> None:
    e_values = (20.0, 35.0, 50.0, 75.0, 100.0)
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        for index, e_over_n in enumerate(e_values):
            for replicate, scale in enumerate((0.99, 1.0, 1.01)):
                insert_workflow_case(
                    store.connection,
                    solver="monte_carlo",
                    e_over_n=e_over_n,
                    replicate=replicate,
                    mean_energy=(2.0 + index) * scale,
                    reduced_mobility=(12.0 + index) * scale,
                    reduced_diffusion_l=(20.0 + index) * scale,
                    reduced_diffusion_t=(30.0 + index) * scale,
                    diagnostics=_mc_rate_diagnostics(
                        event_count=0 if index < 2 else 1,
                        residence_time_s=(replicate + 1) * 1.0e-6,
                    ),
                )
                insert_workflow_rate(
                    store.connection,
                    solver="monte_carlo",
                    e_over_n=e_over_n,
                    replicate=replicate,
                    value=(0.0 if index < 2 else index * 1.0e-18 * scale),
                    process="ionization",
                    process_type="ionization",
                    threshold_eV=15.0,
                )
                insert_workflow_eedf(
                    store.connection,
                    solver="monte_carlo",
                    e_over_n=e_over_n,
                    replicate=replicate,
                    probability_masses=_probability_masses_for_mean(
                        (2.0 + index) * scale,
                        bin_count=8,
                    ),
                )
        store.connection.commit()

    summary = build_tables(
        tmp_path / "workflow.sqlite",
        tmp_path / "tables",
        source="monte_carlo",
        allow_unqualified_mc=True,
    )

    mixture_dir = summary.output_directory / "mixture_0000"
    qualification_path = mixture_dir / "mc_qualification.csv"
    assert qualification_path.is_file()
    assert not (mixture_dir / "rate_coefficients.csv").exists()
    qualification_rows = list(
        csv.DictReader(qualification_path.open(encoding="utf-8"))
    )
    unresolved_tail_rows = qualification_rows[:2]
    assert all(
        "rare_tail" in json.loads(row["active_closure_failure_axes_json"])
        for row in unresolved_tail_rows
    )
    manifest = json.loads((mixture_dir / "manifest.json").read_text())
    assert manifest["mc_qualification"]["coefficient_tables_available"] is False
    assert manifest["mc_qualification"]["table_build_failure"]


def test_censored_rate_relevance_is_a_qualification_failure() -> None:
    process = {
        "species": "Ar",
        "process": "ionization",
        "process_type": "ionization",
        "threshold_eV": 15.0,
    }
    rates = [
        {**process, "E_over_N_Td": 1.0, "rate_coefficient_m3_s": 0.0},
        {**process, "E_over_N_Td": 100.0, "rate_coefficient_m3_s": 1.0e-16},
    ]
    upper = 2.0e-20
    exposure = -math.log(1.0 - 0.95) / upper
    evidence = [
        {
            **process,
            "E_over_N_Td": 1.0,
            "rate_coefficient_mean_m3_s": 0.0,
            "estimate_status": "censored_all_zero",
            "pooled_event_count": 0.0,
            "pooled_target_exposure_s_m3": exposure,
            "pooled_all_zero_upper_95_m3_s": upper,
            "pooled_zero_event_status": "all_zero_upper_95",
            "uncertainty_available": 0,
        },
        {
            **process,
            "E_over_N_Td": 100.0,
            "rate_coefficient_mean_m3_s": 1.0e-16,
            "estimate_status": "observed",
        },
    ]

    with pytest.raises(
        MonteCarloQualificationError,
        match="censored-rate pooled 95% upper bound",
    ):
        _validate_mc_censored_rate_relevance(
            rates,
            evidence,
            minimum_fraction=1.0e-4,
        )


def _insert_qualified_mc_with_decoy_two_term(
    store: WorkflowStore,
    *,
    mc_means: tuple[float, ...],
    decoy_two_term_means: tuple[float, ...],
) -> None:
    e_values = tuple(20.0 + 20.0 * index for index in range(len(mc_means)))
    for index, (e_over_n, mc_mean, decoy_two_term_mean) in enumerate(
        zip(e_values, mc_means, decoy_two_term_means, strict=True)
    ):
        insert_workflow_case(
            store.connection,
            solver="two_term",
            e_over_n=e_over_n,
            mean_energy=decoy_two_term_mean,
            reduced_mobility=10.0 + index,
        )
        insert_workflow_rate(
            store.connection,
            solver="two_term",
            e_over_n=e_over_n,
            value=(index + 1) * 1.0e-16,
            process="ionization",
            process_type="ionization",
            threshold_eV=15.0,
        )
        insert_workflow_eedf(
            store.connection,
            solver="two_term",
            e_over_n=e_over_n,
            probability_masses=_probability_masses_for_mean(
                decoy_two_term_mean,
                bin_count=8,
            ),
        )
        for replicate, scale in enumerate((0.99, 1.0, 1.01)):
            replicate_mean = mc_mean * scale
            insert_workflow_case(
                store.connection,
                solver="monte_carlo",
                e_over_n=e_over_n,
                replicate=replicate,
                mean_energy=replicate_mean,
                reduced_mobility=(12.0 + index) * scale,
                reduced_diffusion_l=(20.0 + index) * scale,
                reduced_diffusion_t=(30.0 + index) * scale,
            )
            insert_workflow_rate(
                store.connection,
                solver="monte_carlo",
                e_over_n=e_over_n,
                replicate=replicate,
                value=(index + 1) * 1.0e-16 * scale,
                process="ionization",
                process_type="ionization",
                threshold_eV=15.0,
            )
            insert_workflow_eedf(
                store.connection,
                solver="monte_carlo",
                e_over_n=e_over_n,
                replicate=replicate,
                probability_masses=_probability_masses_for_mean(
                    replicate_mean,
                    bin_count=8,
                ),
            )


def test_qualified_monte_carlo_rejects_internal_failed_transport_anchor(
    tmp_path: Path,
) -> None:
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        _insert_qualified_mc_with_decoy_two_term(
            store,
            mc_means=(1.0, 2.0, 3.0, 4.0, 5.0),
            decoy_two_term_means=(0.8, 1.5, 2.5, 4.0, 6.0),
        )
        store.connection.execute(
            """
            UPDATE cases
            SET reduced_electron_energy_diffusion_L_m2_s_m3 = 0.0
            WHERE solver = 'monte_carlo' AND e_over_n_Td = 60.0
            """
        )
        store.connection.commit()

    with pytest.raises(
        TableBuildError,
        match="crosses failed internal planned E/N anchors: 60",
    ):
        build_tables(
            tmp_path / "workflow.sqlite",
            tmp_path / "tables",
            source="monte_carlo",
        )


def test_mc_transport_lag_failure_retains_independently_qualified_eedf_and_rates(
    tmp_path: Path,
) -> None:
    database = tmp_path / "workflow.sqlite"
    with WorkflowStore(database) as store:
        _insert_qualified_mc_with_decoy_two_term(
            store,
            mc_means=(1.0, 2.0, 3.0, 4.0, 5.0),
            decoy_two_term_means=(0.8, 1.5, 2.5, 4.0, 6.0),
        )
        rows = store.connection.execute(
            """
            SELECT replicate, diagnostics_json FROM cases
            WHERE solver = 'monte_carlo' AND e_over_n_Td = 20.0
            """
        ).fetchall()
        for row in rows:
            replicate, diagnostics_json = row
            diagnostics = json.loads(str(diagnostics_json))
            estimates = diagnostics["internal_monte_carlo_transport"]["lag_scan"][
                "estimates"
            ]["energy_diffusion_L_m2_s"]
            estimates[-2] = float(estimates[-1]) * 0.5
            store.connection.execute(
                """
                UPDATE cases SET diagnostics_json = ?
                WHERE solver = 'monte_carlo' AND e_over_n_Td = 20.0
                  AND replicate = ?
                """,
                (json.dumps(diagnostics), int(replicate)),
            )
        store.connection.commit()

    build_tables(database, tmp_path / "tables", source="monte_carlo")

    mixture_dir = tmp_path / "tables" / "mixture_0000"
    transport_e = {
        float(row["E_over_N_Td"])
        for row in _read_csv(mixture_dir / "transport_vs_mean_energy.csv")
    }
    eedf_e = {float(row["E_over_N_Td"]) for row in _read_csv(mixture_dir / "eedf.csv")}
    rate_e = {
        float(row["E_over_N_Td"])
        for row in _read_csv(mixture_dir / "rates_vs_mean_energy.csv")
    }
    assert transport_e == {40.0, 60.0, 80.0, 100.0}
    assert eedf_e == {20.0, 40.0, 60.0, 80.0, 100.0}
    assert rate_e == eedf_e

    policy = json.loads((mixture_dir / "manifest.json").read_text())["source_policy"]
    assert policy["transport_eligible_anchor_points"] == 4
    assert policy["eedf_rate_eligible_anchor_points"] == 5
    assert policy["eedf_rate_excluded_E_over_N_Td"] == {}
    assert policy["transport_excluded_E_over_N_Td"]["20"] == [
        "mc_transport_lag_plateau_not_converged:energy_diffusion_L_m2_s"
    ]


def test_qualified_monte_carlo_rejects_nonmonotonic_mean_energy(
    tmp_path: Path,
) -> None:
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        _insert_qualified_mc_with_decoy_two_term(
            store,
            mc_means=(1.0, 4.0, 1.5, 5.0, 6.0),
            decoy_two_term_means=(0.8, 2.0, 3.5, 5.5, 7.0),
        )
        store.connection.commit()

    with pytest.raises(TableBuildError, match="strictly increasing"):
        build_tables(
            tmp_path / "workflow.sqlite",
            tmp_path / "tables",
            source="monte_carlo",
        )


def test_pure_monte_carlo_ignores_decoy_two_term_coverage(
    tmp_path: Path,
) -> None:
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        _insert_qualified_mc_with_decoy_two_term(
            store,
            mc_means=(1.0, 2.0, 3.0, 4.0, 6.0),
            decoy_two_term_means=(0.8, 1.5, 2.5, 3.5, 4.5),
        )
        store.connection.commit()

    build_tables(
        tmp_path / "workflow.sqlite",
        tmp_path / "tables",
        source="monte_carlo",
    )
    manifest = json.loads(
        (tmp_path / "tables" / "mixture_0000" / "manifest.json").read_text()
    )
    assert manifest["source_policy"]["quality_policy"]["failed_points"] == (
        "excluded_per_component_only_outside_contiguous_transport_support_not_repaired"
    )


def _looks_numeric(value: str) -> bool:
    try:
        float(value)
    except ValueError:
        return False
    return True
