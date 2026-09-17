"""Quality evaluation for weighted-growth Monte Carlo transport."""

from __future__ import annotations

from dataclasses import dataclass, field
import json
import math
from typing import Any

import electron_swarm.solvers.monte_carlo.evidence as _mc_evidence

from .contracts import (
    MC_POPULATION_GROWTH_RELATIVE_TOLERANCE,
    MC_TRANSPORT_DIFFUSION_LAG_RELATIVE_TOLERANCE,
    MC_TRANSPORT_LAG_FIELDS,
    MC_TRANSPORT_MIN_ENSEMBLE_REPLICATES,
    MC_TRANSPORT_STATE_RELATIVE_TOLERANCE,
    WEIGHTED_MC_TRANSPORT_COMPONENT_ESTIMATORS,
    WEIGHTED_MC_TRANSPORT_DIFFUSION_IDENTIFICATION_RELATIVE_TOLERANCE,
    WEIGHTED_MC_TRANSPORT_STATIONARITY_FIELDS,
    _normalized_mc_run_settings,
    _paired_log_ratio_absolute_mean,
    _poisson_frequency_interval_95,
    _worst_stationarity_evidence,
    paired_log_ratio_ci95_bound,
)


_PRODUCTION_FIELDS = (
    "drift_velocity_m_s",
    "mobility_m2_V_s",
    "diffusion_L_m2_s",
    "diffusion_T_m2_s",
    "energy_mobility_m2_V_s",
    "energy_diffusion_L_m2_s",
    "energy_diffusion_T_m2_s",
    "mean_energy_eV",
)
_SAMPLING_FIELDS = {
    "physical_time_barriers",
    "resampling_at_common_time_only",
    "barrier_cadence_s",
    "configured_correlation_lag_barriers",
    "transport_block_barriers",
    "lag_barriers",
    "complete_blocks",
    "minimum_blocks",
    "positions_reset_each_block",
    "warmup_accumulators_reset",
}
_STATIONARITY_SAMPLING_FIELDS = {
    "window_definition",
    "state_moment_aggregation",
    "state_estimator",
    "diffusion_estimator",
    "diffusion_lag_barriers",
    "early_block_count",
    "late_block_count",
    "odd_center_block_excluded",
}
_LINEAGE_FIELDS = {
    "lineage_horizon",
    "lag_barriers",
    "minimum_effective_lineage_count",
    "minimum_effective_lineage_fraction",
    "production_lag_barriers",
    "qualification_lag_barriers",
    "required_effective_lineage_count",
    "required_effective_lineage_fraction",
}
_LAG_SCAN_FIELDS = {
    "lag_barriers",
    "lag_s",
    "estimates",
    "completed_blocks",
    "production_lag_barriers",
    "hard_convergence_pair_barriers",
    "supplemental_pair_barriers",
}
_GROWTH_FIELDS = {
    "model",
    "production_elapsed_time_s",
    "cumulative_log_growth",
    "population_log_growth_frequency_s_inv",
    "direct_ionization_frequency_s_inv",
    "event_evidence_model",
    "event_count_ionization_frequency_s_inv",
    "ionization_event_count",
    "weighted_ionization_event_count",
    "secondary_electron_events",
    "event_count_matches_secondary",
    "event_observation_residence_time_s",
    "event_exposure_consistent_across_processes",
    "zero_event_frequency_upper_95_s_inv",
}


@dataclass(slots=True)
class _WeightedTransportEvidence:
    """Normalized ensemble evidence shared by weighted-MC quality stages."""

    reasons: list[str] = field(default_factory=list)
    early_values: dict[str, list[float]] = field(
        default_factory=lambda: {
            name: [] for name in WEIGHTED_MC_TRANSPORT_STATIONARITY_FIELDS
        }
    )
    late_values: dict[str, list[float]] = field(
        default_factory=lambda: {
            name: [] for name in WEIGHTED_MC_TRANSPORT_STATIONARITY_FIELDS
        }
    )
    lag_reference: dict[str, list[float]] = field(
        default_factory=lambda: {name: [] for name in MC_TRANSPORT_LAG_FIELDS}
    )
    lag_candidate: dict[str, list[float]] = field(
        default_factory=lambda: {name: [] for name in MC_TRANSPORT_LAG_FIELDS}
    )
    supplemental_lag_reference: dict[str, list[float]] = field(
        default_factory=lambda: {name: [] for name in MC_TRANSPORT_LAG_FIELDS}
    )
    supplemental_lag_candidate: dict[str, list[float]] = field(
        default_factory=lambda: {name: [] for name in MC_TRANSPORT_LAG_FIELDS}
    )
    supplemental_lineage_replicates: int = 0
    growth_direct: list[float] = field(default_factory=list)
    growth_population: list[float] = field(default_factory=list)
    growth_event: list[float] = field(default_factory=list)
    growth_event_models: list[str] = field(default_factory=list)
    growth_event_counts: list[int] = field(default_factory=list)
    growth_weighted_event_counts: list[float] = field(default_factory=list)
    growth_event_exposures_s: list[float] = field(default_factory=list)
    growth_lag_horizons_s: list[float] = field(default_factory=list)
    case_mean_energy: list[float] = field(default_factory=list)
    transport_mean_energy: list[float] = field(default_factory=list)
    seeds: list[int] = field(default_factory=list)
    common_settings: str | None = None
    common_sampling: str | None = None
    complete: bool = True

    def reject(self, reason: str) -> None:
        if reason not in self.reasons:
            self.reasons.append(reason)

    def incomplete(self, reason: str) -> None:
        self.reject(reason)
        self.complete = False


@dataclass(frozen=True, slots=True)
class _WeightedReplicaContext:
    particles: int
    production_lag: int
    lag_barriers: list[int]
    block_barriers: int
    lineage_qualification_lags: list[int]
    cadence_s: float


@dataclass(frozen=True, slots=True)
class _WeightedMomentEvidence:
    production_values: dict[str, float]
    lag_seconds: list[float]
    production_index: int


@dataclass(slots=True)
class _WeightedPrecisionEvidence:
    stationarity_bounds: dict[str, float] = field(default_factory=dict)
    lag_bounds: dict[str, float] = field(default_factory=dict)
    supplemental_lag_bounds: dict[str, float] = field(default_factory=dict)
    supplemental_replicates: int = 0
    mean_energy_bound: float | None = None
    stationarity_worst_field: str | None = None
    stationarity_max: float | None = None
    stationarity_tolerance: float | None = None


@dataclass(slots=True)
class _WeightedGrowthEvidence:
    bounds: dict[str, float] = field(default_factory=dict)
    mode: str | None = None
    poisson_ratio: float | None = None
    sparse_metric: float | None = None
    pooled_event_count: int | None = None
    pooled_exposure_s: float | None = None
    pooled_frequency_low: float | None = None
    pooled_frequency_high: float | None = None
    pooled_direct_within: int | None = None
    pooled_population_within: int | None = None


def _normalize_weighted_replica_input(
    item: dict[str, Any], evidence: _WeightedTransportEvidence
) -> _WeightedReplicaContext:
    if not isinstance(item, dict) or item.get("_diagnostic_invalid"):
        raise ValueError("invalid diagnostic")
    if (
        item.get("estimator") != "direct_mc_weighted_growth_block_lag_flux"
        or item.get("estimator_schema_version")
        != _mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
    ):
        evidence.incomplete("mc_transport_estimator_provenance_mismatch")
    if (
        item.get("energy_transport_status")
        != "direct_mc_candidate_for_ensemble_qualification"
        or item.get("finite") is not True
        or item.get("strictly_positive") is not True
        or item.get("multiple_time_blocks") is not True
        or item.get("complete_observation") is not True
        or item.get("lineage_qualified") is not True
        or item.get("transport_definition") != "flux"
        or item.get("bulk_transport_identified") is not False
        or item.get("cross_gradient_response_identified") is not False
    ):
        evidence.incomplete("mc_weighted_growth_transport_basic_unqualified")
    if item.get("component_estimators") != (WEIGHTED_MC_TRANSPORT_COMPONENT_ESTIMATORS):
        evidence.incomplete("mc_transport_component_estimator_contract_mismatch")

    provenance = item.get("mc_run_provenance")
    if (
        not isinstance(provenance, dict)
        or provenance.get("population_model") != "weighted_branching"
    ):
        raise ValueError("invalid weighted-growth provenance")
    settings = _normalized_mc_run_settings(provenance)
    seed = provenance.get("seed")
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError("weighted-growth seed is invalid")
    if any(
        isinstance(provenance.get(name), bool)
        or not isinstance(provenance.get(name), int)
        for name in (
            "particles",
            "warmup_collisions",
            "production_collisions",
            "tail_max_collisions",
            "tail_collisions_executed",
        )
    ):
        raise ValueError("weighted-growth integer settings are invalid")
    if any(
        not isinstance(provenance.get(name), str) or not provenance.get(name)
        for name in (
            "population_model",
            "collision_clock",
            "angular_scattering_model",
            "ionization_source_model",
            "high_energy_extrapolation",
        )
    ) or not isinstance(provenance.get("magnetic_enabled"), bool):
        raise ValueError("weighted-growth model settings are invalid")
    numeric_names = (
        "gas_number_density_m3",
        "electric_field_V_m",
        "trial_collision_frequency_s_inv",
        "tail_rate_rse_trigger",
        "magnetic_B_T",
        "magnetic_angle_EB_deg",
        "max_energy_limit_eV",
    )
    numeric = [float(provenance[name]) for name in numeric_names]
    if (
        any(not math.isfinite(value) for value in numeric)
        or float(provenance["gas_number_density_m3"]) <= 0.0
        or float(provenance["electric_field_V_m"]) <= 0.0
        or float(provenance["trial_collision_frequency_s_inv"]) <= 0.0
        or int(provenance["particles"]) < 2
        or int(provenance["warmup_collisions"]) < 0
        or int(provenance["production_collisions"]) <= 0
        or int(provenance["tail_max_collisions"]) < 0
        or int(provenance["tail_collisions_executed"])
        not in {0, int(provenance["tail_max_collisions"])}
        or not (0.0 < float(provenance["tail_rate_rse_trigger"]) <= 1.0)
        or provenance.get("magnetic_enabled") is not False
    ):
        raise ValueError("weighted-growth run settings are invalid")
    evidence.seeds.append(seed)
    if evidence.common_settings is None:
        evidence.common_settings = settings
    elif settings != evidence.common_settings:
        evidence.incomplete("mc_transport_replica_settings_mismatch")

    sampling = item.get("block_lag_sampling")
    if not isinstance(sampling, dict) or set(sampling) != _SAMPLING_FIELDS:
        raise ValueError("weighted-growth sampling evidence is missing")
    configured_lag = sampling["configured_correlation_lag_barriers"]
    if isinstance(configured_lag, bool) or not isinstance(configured_lag, int):
        raise ValueError("weighted-growth configured lag is invalid")
    production_lag = int(configured_lag)
    if production_lag < 4 or production_lag & (production_lag - 1):
        raise ValueError("weighted-growth configured lag is invalid")
    lag_barriers = [
        production_lag // 4,
        production_lag // 2,
        production_lag,
        2 * production_lag,
    ]
    block_barriers = 2 * production_lag
    lineage_qualification_lags = [production_lag // 2, production_lag]
    cadence_s = float(sampling["barrier_cadence_s"])
    expected_cadence = 1.0 / float(provenance["trial_collision_frequency_s_inv"])
    if (
        sampling["physical_time_barriers"] is not True
        or sampling["resampling_at_common_time_only"] is not True
        or sampling["positions_reset_each_block"] is not True
        or sampling["warmup_accumulators_reset"] is not True
        or not math.isclose(
            cadence_s,
            expected_cadence,
            rel_tol=1.0e-12,
            abs_tol=0.0,
        )
        or sampling["transport_block_barriers"] != block_barriers
        or int(provenance["warmup_collisions"]) <= 0
        or int(provenance["production_collisions"]) % block_barriers != 0
        or sampling["lag_barriers"] != lag_barriers
        or sampling["minimum_blocks"] != _mc_evidence.WEIGHTED_MC_TRANSPORT_MIN_BLOCKS
        or int(sampling["complete_blocks"])
        < _mc_evidence.WEIGHTED_MC_TRANSPORT_MIN_BLOCKS
    ):
        raise ValueError("weighted-growth synchronization is invalid")
    sampling_contract = json.dumps(
        sampling,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    if evidence.common_sampling is None:
        evidence.common_sampling = sampling_contract
    elif sampling_contract != evidence.common_sampling:
        evidence.incomplete("mc_transport_replica_sampling_mismatch")

    stationarity_sampling = item.get("stationarity_sampling")
    if (
        not isinstance(stationarity_sampling, dict)
        or set(stationarity_sampling) != _STATIONARITY_SAMPLING_FIELDS
        or stationarity_sampling["window_definition"]
        != "first_half_vs_second_half_complete_transport_blocks"
        or stationarity_sampling["state_moment_aggregation"]
        != "sum_raw_weighted_residence_moments_before_ratio"
        or stationarity_sampling["state_estimator"]
        != "same_residence_ratio_formulas_as_production"
        or stationarity_sampling["diffusion_estimator"]
        != "mean_fixed_production_lag_snapshot_covariance_per_window"
        or int(stationarity_sampling["diffusion_lag_barriers"]) != production_lag
        or int(stationarity_sampling["early_block_count"])
        != int(sampling["complete_blocks"]) // 2
        or int(stationarity_sampling["late_block_count"])
        != int(sampling["complete_blocks"]) // 2
        or not isinstance(stationarity_sampling["odd_center_block_excluded"], bool)
        or stationarity_sampling["odd_center_block_excluded"]
        != bool(
            int(sampling["complete_blocks"]) > 1
            and int(sampling["complete_blocks"]) % 2
        )
    ):
        raise ValueError("weighted-growth stationarity estimator is invalid")
    return _WeightedReplicaContext(
        particles=int(provenance["particles"]),
        production_lag=production_lag,
        lag_barriers=lag_barriers,
        block_barriers=block_barriers,
        lineage_qualification_lags=lineage_qualification_lags,
        cadence_s=cadence_s,
    )


def _collect_weighted_lineage(
    item: dict[str, Any],
    replica: _WeightedReplicaContext,
    evidence: _WeightedTransportEvidence,
) -> bool:
    lineage = item.get("lineage_sampling")
    if not isinstance(lineage, dict) or set(lineage) != _LINEAGE_FIELDS:
        raise ValueError("weighted-growth lineage evidence is missing")
    lineage_lags = [int(value) for value in lineage["lag_barriers"]]
    lineage_counts = [
        float(value) for value in lineage["minimum_effective_lineage_count"]
    ]
    lineage_fractions = [
        float(value) for value in lineage["minimum_effective_lineage_fraction"]
    ]
    qualification_lags = [int(value) for value in lineage["qualification_lag_barriers"]]
    required_count = float(lineage["required_effective_lineage_count"])
    required_fraction = float(lineage["required_effective_lineage_fraction"])
    if (
        lineage["lineage_horizon"] != "each_sampled_lag"
        or lineage_lags != replica.lag_barriers
        or len(lineage_counts) != len(lineage_lags)
        or len(lineage_fractions) != len(lineage_lags)
        or int(lineage["production_lag_barriers"]) != replica.production_lag
        or qualification_lags != replica.lineage_qualification_lags
        or not math.isclose(
            required_count,
            _mc_evidence.WEIGHTED_MC_TRANSPORT_MIN_EFFECTIVE_LINEAGE_COUNT,
            rel_tol=0.0,
            abs_tol=0.0,
        )
        or not math.isclose(
            required_fraction,
            _mc_evidence.WEIGHTED_MC_TRANSPORT_MIN_EFFECTIVE_LINEAGE_FRACTION,
            rel_tol=0.0,
            abs_tol=0.0,
        )
        or any(not math.isfinite(count) or count <= 0.0 for count in lineage_counts)
        or any(
            not math.isfinite(fraction)
            or fraction <= 0.0
            or not math.isclose(
                fraction,
                count / replica.particles,
                rel_tol=1.0e-12,
                abs_tol=0.0,
            )
            for count, fraction in zip(
                lineage_counts,
                lineage_fractions,
                strict=True,
            )
        )
    ):
        raise ValueError("weighted-growth lineage evidence is invalid")
    lineage_fraction_by_lag = dict(zip(lineage_lags, lineage_fractions, strict=True))
    lineage_count_by_lag = dict(zip(lineage_lags, lineage_counts, strict=True))
    hard_lineage_ok = all(
        lineage_count_by_lag[lag] >= required_count
        and lineage_fraction_by_lag[lag] >= required_fraction
        for lag in qualification_lags
    )
    if not hard_lineage_ok:
        evidence.incomplete("mc_weighted_growth_lineage_degenerate")
    supplemental_lineage_ok = bool(
        lineage_count_by_lag[replica.block_barriers] >= required_count
        and lineage_fraction_by_lag[replica.block_barriers] >= required_fraction
    )
    evidence.supplemental_lineage_replicates += int(supplemental_lineage_ok)
    return supplemental_lineage_ok


def _collect_weighted_moments(
    item: dict[str, Any],
    replica: _WeightedReplicaContext,
    evidence: _WeightedTransportEvidence,
    *,
    supplemental_lineage_ok: bool,
) -> _WeightedMomentEvidence:
    production = item.get("production_estimates")
    early = item.get("time_stationarity_window_estimate_early")
    late = item.get("time_stationarity_window_estimate_late")
    if (
        not isinstance(production, dict)
        or set(production) != set(_PRODUCTION_FIELDS)
        or not isinstance(early, dict)
        or set(early) != set(WEIGHTED_MC_TRANSPORT_STATIONARITY_FIELDS)
        or not isinstance(late, dict)
        or set(late) != set(WEIGHTED_MC_TRANSPORT_STATIONARITY_FIELDS)
    ):
        raise ValueError("weighted-growth moment evidence is missing")
    production_values = {name: float(production[name]) for name in _PRODUCTION_FIELDS}
    if any(
        not math.isfinite(value) or value <= 0.0 for value in production_values.values()
    ):
        raise ValueError("weighted-growth production moment is invalid")
    for name in WEIGHTED_MC_TRANSPORT_STATIONARITY_FIELDS:
        early_value = float(early[name])
        late_value = float(late[name])
        if (
            not math.isfinite(early_value)
            or not math.isfinite(late_value)
            or early_value <= 0.0
            or late_value <= 0.0
        ):
            raise ValueError("weighted-growth stationarity value is invalid")
        evidence.early_values[name].append(early_value)
        evidence.late_values[name].append(late_value)

    lag_scan = item.get("lag_scan")
    if not isinstance(lag_scan, dict) or set(lag_scan) != _LAG_SCAN_FIELDS:
        raise ValueError("weighted-growth lag evidence is missing")
    lag_barriers = [int(value) for value in lag_scan["lag_barriers"]]
    lag_seconds = [float(value) for value in lag_scan["lag_s"]]
    completed_blocks = [int(value) for value in lag_scan["completed_blocks"]]
    production_lag = int(lag_scan["production_lag_barriers"])
    hard_pair = [int(value) for value in lag_scan["hard_convergence_pair_barriers"]]
    supplemental_pair = [int(value) for value in lag_scan["supplemental_pair_barriers"]]
    lag_estimates = lag_scan["estimates"]
    if (
        lag_barriers != replica.lag_barriers
        or len(lag_seconds) != len(lag_barriers)
        or len(completed_blocks) != len(lag_barriers)
        or production_lag != replica.production_lag
        or hard_pair != [replica.production_lag // 2, production_lag]
        or supplemental_pair != [production_lag, replica.block_barriers]
        or not isinstance(lag_estimates, dict)
        or set(lag_estimates) != set(MC_TRANSPORT_LAG_FIELDS)
        or any(
            value < _mc_evidence.WEIGHTED_MC_TRANSPORT_MIN_BLOCKS
            for value in completed_blocks
        )
        or any(
            not math.isclose(
                value,
                lag * replica.cadence_s,
                rel_tol=1.0e-12,
                abs_tol=0.0,
            )
            for value, lag in zip(lag_seconds, lag_barriers, strict=True)
        )
    ):
        raise ValueError("weighted-growth lag schedule is invalid")
    hard_reference_index = lag_barriers.index(hard_pair[0])
    production_index = lag_barriers.index(production_lag)
    supplemental_index = lag_barriers.index(supplemental_pair[1])
    parsed_lag_estimates: dict[str, list[float]] = {}
    supplemental_moments_valid = True
    for name in MC_TRANSPORT_LAG_FIELDS:
        values = [float(value) for value in lag_estimates[name]]
        if len(values) != len(lag_barriers) or any(
            not math.isfinite(values[index]) or values[index] <= 0.0
            for index in range(production_index + 1)
        ):
            raise ValueError("weighted-growth lag moment is invalid")
        if not math.isclose(
            production_values[name],
            values[production_index],
            rel_tol=1.0e-12,
            abs_tol=0.0,
        ):
            raise ValueError("weighted-growth production lag is inconsistent")
        evidence.lag_reference[name].append(values[hard_reference_index])
        evidence.lag_candidate[name].append(values[production_index])
        parsed_lag_estimates[name] = values
        supplemental_moments_valid = bool(
            supplemental_moments_valid
            and math.isfinite(values[supplemental_index])
            and values[supplemental_index] > 0.0
        )
    if supplemental_lineage_ok and supplemental_moments_valid:
        for name in MC_TRANSPORT_LAG_FIELDS:
            values = parsed_lag_estimates[name]
            evidence.supplemental_lag_reference[name].append(values[production_index])
            evidence.supplemental_lag_candidate[name].append(values[supplemental_index])
    return _WeightedMomentEvidence(
        production_values=production_values,
        lag_seconds=lag_seconds,
        production_index=production_index,
    )


def _collect_weighted_growth(
    item: dict[str, Any],
    moments: _WeightedMomentEvidence,
    evidence: _WeightedTransportEvidence,
) -> None:
    growth = item.get("population_growth_consistency")
    if not isinstance(growth, dict) or set(growth) != _GROWTH_FIELDS:
        raise ValueError("weighted-growth population evidence is missing")
    if (
        growth["model"] != "explicit_ionization_branching_no_attachment"
        or growth["event_count_matches_secondary"] is not True
        or growth["event_exposure_consistent_across_processes"] is not True
    ):
        raise ValueError("weighted-growth population model is inconsistent")
    elapsed = float(growth["production_elapsed_time_s"])
    log_growth = float(growth["cumulative_log_growth"])
    population_frequency = float(growth["population_log_growth_frequency_s_inv"])
    direct_frequency = float(growth["direct_ionization_frequency_s_inv"])
    event_evidence_model = str(growth["event_evidence_model"])
    event_frequency = float(growth["event_count_ionization_frequency_s_inv"])
    event_exposure = float(growth["event_observation_residence_time_s"])
    ionization_events = growth["ionization_event_count"]
    weighted_ionization_events = float(growth["weighted_ionization_event_count"])
    secondary_events = growth["secondary_electron_events"]
    if (
        any(
            not math.isfinite(value) or value < 0.0
            for value in (
                log_growth,
                population_frequency,
                direct_frequency,
                event_frequency,
                event_exposure,
                weighted_ionization_events,
            )
        )
        or not math.isfinite(elapsed)
        or elapsed <= 0.0
        or isinstance(ionization_events, bool)
        or isinstance(secondary_events, bool)
        or int(ionization_events) < 0
        or int(secondary_events) < 0
        or int(ionization_events) != int(secondary_events)
        or not math.isclose(
            population_frequency,
            log_growth / elapsed,
            rel_tol=1.0e-12,
            abs_tol=0.0,
        )
    ):
        raise ValueError("weighted-growth population evidence is invalid")
    if event_exposure > 0.0:
        expected_event_frequency = (
            weighted_ionization_events / event_exposure
            if event_evidence_model == "weighted_ensemble_correlated_events"
            else int(ionization_events) / event_exposure
        )
        if not math.isclose(
            event_frequency,
            expected_event_frequency,
            rel_tol=1.0e-12,
            abs_tol=0.0,
        ):
            raise ValueError("weighted-growth event frequency is invalid")
    elif int(ionization_events) != 0 or event_frequency != 0.0:
        raise ValueError("weighted-growth event exposure is invalid")
    zero_upper = growth["zero_event_frequency_upper_95_s_inv"]
    if event_evidence_model == "weighted_ensemble_correlated_events":
        if zero_upper is not None:
            raise ValueError("weighted-ensemble events cannot provide a Poisson bound")
    elif event_evidence_model != "raw_macro_poisson":
        raise ValueError("weighted-growth event evidence model is invalid")
    elif int(ionization_events) == 0 and event_exposure > 0.0:
        if zero_upper is None or not math.isclose(
            float(zero_upper),
            -math.log(0.05) / event_exposure,
            rel_tol=1.0e-12,
            abs_tol=0.0,
        ):
            raise ValueError("weighted-growth zero-event bound is invalid")
    elif int(ionization_events) > 0 and zero_upper is not None:
        raise ValueError("weighted-growth zero-event bound is invalid")
    elif (
        int(ionization_events) == 0
        and event_exposure == 0.0
        and (zero_upper is None or float(zero_upper) != 0.0)
    ):
        raise ValueError("weighted-growth zero-event bound is invalid")
    count = int(ionization_events)
    if event_evidence_model == "raw_macro_poisson":
        _, poisson_upper = _poisson_frequency_interval_95(count, event_exposure)
        if (
            count == 0
            and zero_upper is not None
            and not math.isclose(
                poisson_upper,
                float(zero_upper),
                rel_tol=1.0e-12,
                abs_tol=0.0,
            )
        ):
            raise ValueError("weighted-growth Poisson evidence is inconsistent")
    evidence.growth_direct.append(direct_frequency)
    evidence.growth_population.append(population_frequency)
    evidence.growth_event.append(event_frequency)
    evidence.growth_event_models.append(event_evidence_model)
    evidence.growth_event_counts.append(count)
    evidence.growth_weighted_event_counts.append(weighted_ionization_events)
    evidence.growth_event_exposures_s.append(event_exposure)
    evidence.growth_lag_horizons_s.append(moments.lag_seconds[moments.production_index])


def _collect_weighted_report(
    item: dict[str, Any],
    moments: _WeightedMomentEvidence,
    evidence: _WeightedTransportEvidence,
) -> None:
    reported = item.get("_case_reported_transport")
    if not isinstance(reported, dict):
        raise ValueError("weighted-growth case transport is missing")
    for name in _PRODUCTION_FIELDS[:-1]:
        if not math.isclose(
            float(reported[name]),
            moments.production_values[name],
            rel_tol=1.0e-12,
            abs_tol=0.0,
        ):
            evidence.incomplete("mc_transport_report_not_production_estimate")
    evidence.case_mean_energy.append(float(item["_case_mean_energy_eV"]))
    evidence.transport_mean_energy.append(moments.production_values["mean_energy_eV"])


def _collect_weighted_evidence(
    diagnostics: list[dict[str, Any]],
) -> _WeightedTransportEvidence:
    evidence = _WeightedTransportEvidence()
    if len(diagnostics) < MC_TRANSPORT_MIN_ENSEMBLE_REPLICATES:
        evidence.reject("mc_transport_ensemble_replicates_insufficient")
    for item in diagnostics:
        try:
            replica = _normalize_weighted_replica_input(item, evidence)
            supplemental_lineage_ok = _collect_weighted_lineage(item, replica, evidence)
            moments = _collect_weighted_moments(
                item,
                replica,
                evidence,
                supplemental_lineage_ok=supplemental_lineage_ok,
            )
            _collect_weighted_growth(item, moments, evidence)
            _collect_weighted_report(item, moments, evidence)
        except (KeyError, TypeError, ValueError, OverflowError):
            evidence.incomplete("mc_weighted_growth_ensemble_evidence_missing")
    if len(evidence.seeds) != len(diagnostics) or len(set(evidence.seeds)) != len(
        diagnostics
    ):
        evidence.incomplete("mc_transport_replica_seeds_missing_or_not_distinct")
    return evidence


def _evaluate_weighted_coefficient_precision(
    evidence: _WeightedTransportEvidence,
) -> _WeightedPrecisionEvidence:
    precision = _WeightedPrecisionEvidence()
    if evidence.complete:
        for name in WEIGHTED_MC_TRANSPORT_STATIONARITY_FIELDS:
            bound = paired_log_ratio_ci95_bound(
                evidence.early_values[name], evidence.late_values[name]
            )
            if bound is None:
                evidence.complete = False
                break
            precision.stationarity_bounds[name] = bound
    if evidence.complete:
        for name in MC_TRANSPORT_LAG_FIELDS:
            bound = paired_log_ratio_ci95_bound(
                evidence.lag_reference[name], evidence.lag_candidate[name]
            )
            if bound is None:
                evidence.complete = False
                break
            precision.lag_bounds[name] = bound
    precision.supplemental_replicates = len(
        evidence.supplemental_lag_reference[MC_TRANSPORT_LAG_FIELDS[0]]
    )
    if (
        evidence.supplemental_lineage_replicates >= MC_TRANSPORT_MIN_ENSEMBLE_REPLICATES
        and precision.supplemental_replicates < MC_TRANSPORT_MIN_ENSEMBLE_REPLICATES
    ):
        evidence.reject("mc_weighted_growth_supplemental_moment_evidence_missing")
    if precision.supplemental_replicates >= MC_TRANSPORT_MIN_ENSEMBLE_REPLICATES:
        for name in MC_TRANSPORT_LAG_FIELDS:
            bound = paired_log_ratio_ci95_bound(
                evidence.supplemental_lag_reference[name],
                evidence.supplemental_lag_candidate[name],
            )
            if bound is not None:
                precision.supplemental_lag_bounds[name] = bound
    return precision


def _growth_arrays_complete(
    evidence: _WeightedTransportEvidence, replicate_count: int
) -> bool:
    return not (
        len(evidence.growth_direct) != replicate_count
        or len(evidence.growth_population) != replicate_count
        or len(evidence.growth_event) != replicate_count
        or len(evidence.growth_event_models) != replicate_count
        or len(evidence.growth_event_counts) != replicate_count
        or len(evidence.growth_weighted_event_counts) != replicate_count
        or len(evidence.growth_event_exposures_s) != replicate_count
        or len(evidence.growth_lag_horizons_s) != replicate_count
    )


def _paired_growth_bounds(
    evidence: _WeightedTransportEvidence, *, weighted_events: bool
) -> dict[str, float]:
    paired_growth_defined = all(
        min(direct, population, event) > 0.0
        for direct, population, event in zip(
            evidence.growth_direct,
            evidence.growth_population,
            evidence.growth_event,
            strict=True,
        )
    )
    if not paired_growth_defined:
        return {}
    population_bound = paired_log_ratio_ci95_bound(
        evidence.growth_direct,
        evidence.growth_population,
    )
    event_bound = paired_log_ratio_ci95_bound(
        evidence.growth_direct,
        evidence.growth_event,
    )
    if population_bound is None or event_bound is None:
        return {}
    return {
        "population_vs_direct": population_bound,
        (
            "weighted_event_vs_direct" if weighted_events else "event_vs_direct"
        ): event_bound,
    }


def _evaluate_raw_poisson_growth(
    evidence: _WeightedTransportEvidence,
    growth: _WeightedGrowthEvidence,
) -> None:
    assert growth.pooled_event_count is not None
    assert growth.pooled_exposure_s is not None
    growth.pooled_frequency_low, growth.pooled_frequency_high = (
        _poisson_frequency_interval_95(
            growth.pooled_event_count,
            growth.pooled_exposure_s,
        )
    )
    pooled_direct = float(
        sum(
            value * exposure
            for value, exposure in zip(
                evidence.growth_direct,
                evidence.growth_event_exposures_s,
                strict=True,
            )
        )
        / growth.pooled_exposure_s
    )
    pooled_population = float(
        sum(
            value * exposure
            for value, exposure in zip(
                evidence.growth_population,
                evidence.growth_event_exposures_s,
                strict=True,
            )
        )
        / growth.pooled_exposure_s
    )
    growth.pooled_direct_within = int(
        growth.pooled_frequency_low <= pooled_direct <= growth.pooled_frequency_high
    )
    growth.pooled_population_within = int(
        growth.pooled_frequency_low <= pooled_population <= growth.pooled_frequency_high
    )

    def interval_ratio(value: float) -> float:
        upper_ratio = value / growth.pooled_frequency_high
        lower_ratio = (
            growth.pooled_frequency_low / value
            if value > 0.0
            else (0.0 if growth.pooled_frequency_low == 0.0 else math.inf)
        )
        return max(upper_ratio, lower_ratio)

    growth.poisson_ratio = max(
        interval_ratio(pooled_direct),
        interval_ratio(pooled_population),
    )
    growth.sparse_metric = -math.expm1(
        -growth.pooled_frequency_high * evidence.growth_lag_horizons_s[0]
    )
    growth.bounds = _paired_growth_bounds(evidence, weighted_events=False)
    if growth.sparse_metric <= MC_POPULATION_GROWTH_RELATIVE_TOLERANCE:
        growth.mode = "sparse_negligible_transport"
        if not (growth.pooled_direct_within and growth.pooled_population_within):
            evidence.reject("mc_weighted_growth_sparse_poisson_inconsistent")
    elif all(count > 0 for count in evidence.growth_event_counts) and growth.bounds:
        growth.mode = "positive_paired_log"
        if max(growth.bounds.values()) > MC_POPULATION_GROWTH_RELATIVE_TOLERANCE:
            evidence.reject("mc_weighted_growth_population_rate_inconsistent")
    else:
        evidence.reject("mc_weighted_growth_branching_not_negligible")
        evidence.reject("mc_weighted_growth_population_rate_unresolved")


def _evaluate_correlated_weighted_growth(
    evidence: _WeightedTransportEvidence,
    growth: _WeightedGrowthEvidence,
) -> None:
    growth.sparse_metric = max(
        -math.expm1(-max(direct, population) * lag)
        for direct, population, lag in zip(
            evidence.growth_direct,
            evidence.growth_population,
            evidence.growth_lag_horizons_s,
            strict=True,
        )
    )
    growth.bounds = _paired_growth_bounds(evidence, weighted_events=True)
    if growth.sparse_metric <= MC_POPULATION_GROWTH_RELATIVE_TOLERANCE:
        growth.mode = "weighted_ensemble_sparse_negligible_transport"
    elif growth.bounds:
        growth.mode = "weighted_ensemble_positive_paired_log"
        if max(growth.bounds.values()) > MC_POPULATION_GROWTH_RELATIVE_TOLERANCE:
            evidence.reject("mc_weighted_growth_population_rate_inconsistent")
    else:
        evidence.reject("mc_weighted_growth_branching_not_negligible")
        evidence.reject("mc_weighted_growth_population_rate_unresolved")


def _evaluate_weighted_growth_precision(
    evidence: _WeightedTransportEvidence, replicate_count: int
) -> _WeightedGrowthEvidence:
    growth = _WeightedGrowthEvidence()
    if evidence.complete and not _growth_arrays_complete(evidence, replicate_count):
        evidence.complete = False
    if not evidence.complete:
        return growth
    growth.pooled_event_count = sum(evidence.growth_event_counts)
    growth.pooled_exposure_s = float(sum(evidence.growth_event_exposures_s))
    if growth.pooled_exposure_s <= 0.0:
        evidence.complete = False
    elif all(model == "raw_macro_poisson" for model in evidence.growth_event_models):
        _evaluate_raw_poisson_growth(evidence, growth)
    elif all(
        model == "weighted_ensemble_correlated_events"
        for model in evidence.growth_event_models
    ):
        _evaluate_correlated_weighted_growth(evidence, growth)
    else:
        evidence.reject("mc_weighted_growth_event_evidence_models_mixed")
        evidence.complete = False
    return growth


def _finish_weighted_precision(
    evidence: _WeightedTransportEvidence,
    precision: _WeightedPrecisionEvidence,
) -> None:
    if evidence.complete:
        precision.mean_energy_bound = paired_log_ratio_ci95_bound(
            evidence.transport_mean_energy,
            evidence.case_mean_energy,
        )
        if precision.mean_energy_bound is None:
            evidence.complete = False
    if not evidence.complete:
        evidence.reject("mc_weighted_growth_ensemble_evidence_missing")
    (
        precision.stationarity_worst_field,
        precision.stationarity_max,
        precision.stationarity_tolerance,
    ) = _worst_stationarity_evidence(precision.stationarity_bounds)


def _apply_weighted_acceptance(
    evidence: _WeightedTransportEvidence,
    precision: _WeightedPrecisionEvidence,
) -> None:
    for name, bound in precision.stationarity_bounds.items():
        limit = (
            MC_TRANSPORT_DIFFUSION_LAG_RELATIVE_TOLERANCE
            if "diffusion" in name
            else MC_TRANSPORT_STATE_RELATIVE_TOLERANCE
        )
        if bound > limit:
            evidence.reject(f"mc_weighted_growth_stationarity_not_converged:{name}")
    for name, bound in precision.lag_bounds.items():
        if bound > WEIGHTED_MC_TRANSPORT_DIFFUSION_IDENTIFICATION_RELATIVE_TOLERANCE:
            evidence.reject(f"mc_weighted_growth_lag_not_converged:{name}")
    for name, bound in precision.supplemental_lag_bounds.items():
        if bound > WEIGHTED_MC_TRANSPORT_DIFFUSION_IDENTIFICATION_RELATIVE_TOLERANCE:
            evidence.reject(f"mc_weighted_growth_lag_not_converged:{name}")
    if (
        precision.mean_energy_bound is not None
        and precision.mean_energy_bound > MC_TRANSPORT_STATE_RELATIVE_TOLERANCE
    ):
        evidence.reject("mc_transport_case_mean_energy_inconsistent")


def _weighted_summary(
    diagnostics: list[dict[str, Any]],
    evidence: _WeightedTransportEvidence,
    precision: _WeightedPrecisionEvidence,
    growth: _WeightedGrowthEvidence,
) -> dict[str, Any]:
    qualified = not evidence.reasons
    return {
        "solver_diagnostics_available": 1,
        "solver_converged": None,
        "solver_transport_qualified": int(qualified),
        "solver_iterations_max": None,
        "solver_residual_L1_max": None,
        "solver_residual_tolerance": None,
        "solver_tail_probability_max": None,
        "solver_tail_probability_target": None,
        "solver_edge_to_peak_max": None,
        "solver_edge_to_peak_target": None,
        "solver_grid_max_eV_max": None,
        "solver_grid_max_limit_eV": None,
        "solver_grid_limit_hit": None,
        "solver_diagnostics_passed": int(qualified),
        "solver_transport_replicates": len(diagnostics),
        "solver_origin_stationarity_limiting_field": (
            precision.stationarity_worst_field
        ),
        "solver_origin_stationarity_limiting_relative_ci95_bound": (
            precision.stationarity_max
        ),
        "solver_mean_energy_stationarity_relative_ci95_bound": (
            precision.stationarity_bounds.get("mean_energy_eV")
        ),
        "solver_mobility_stationarity_absolute_log_drift": (
            _paired_log_ratio_absolute_mean(
                evidence.early_values["mobility_m2_V_s"],
                evidence.late_values["mobility_m2_V_s"],
            )
        ),
        "solver_lag_convergence_max_relative_ci95_bound": max(
            precision.lag_bounds.values(), default=None
        ),
        "solver_lag_supplemental_replicates": precision.supplemental_replicates,
        "solver_lag_supplemental_max_relative_ci95_bound": max(
            precision.supplemental_lag_bounds.values(), default=None
        ),
        "solver_transport_mean_energy_max_relative_ci95_bound": (
            precision.mean_energy_bound
        ),
        "solver_population_growth_gate_mode": growth.mode,
        "solver_population_growth_max_relative_ci95_bound": (
            max(growth.bounds.values()) if growth.bounds else None
        ),
        "solver_population_growth_poisson_interval_max_ratio": (growth.poisson_ratio),
        "solver_population_growth_sparse_max_metric": growth.sparse_metric,
        "solver_population_growth_pooled_event_count": growth.pooled_event_count,
        "solver_population_growth_pooled_exposure_s": growth.pooled_exposure_s,
        "solver_population_growth_pooled_frequency_ci95_low_s_inv": (
            growth.pooled_frequency_low
        ),
        "solver_population_growth_pooled_frequency_ci95_high_s_inv": (
            growth.pooled_frequency_high
        ),
        "solver_population_growth_pooled_direct_within_ci95": (
            growth.pooled_direct_within
        ),
        "solver_population_growth_pooled_population_within_ci95": (
            growth.pooled_population_within
        ),
        "solver_origin_stationarity_limiting_relative_tolerance": (
            precision.stationarity_tolerance
        ),
        "solver_lag_convergence_relative_tolerance": (
            WEIGHTED_MC_TRANSPORT_DIFFUSION_IDENTIFICATION_RELATIVE_TOLERANCE
        ),
        "solver_transport_mean_energy_relative_tolerance": (
            MC_TRANSPORT_STATE_RELATIVE_TOLERANCE
        ),
        "solver_population_growth_relative_tolerance": (
            MC_POPULATION_GROWTH_RELATIVE_TOLERANCE
        ),
        "solver_failure_reasons": evidence.reasons,
    }


def summarize_weighted_mc_transport_diagnostics(
    diagnostics: list[dict[str, Any]],
) -> dict[str, Any]:
    """Qualify block-lag weighted-growth flux moments fail closed."""

    evidence = _collect_weighted_evidence(diagnostics)
    precision = _evaluate_weighted_coefficient_precision(evidence)
    growth = _evaluate_weighted_growth_precision(evidence, len(diagnostics))
    _finish_weighted_precision(evidence, precision)
    _apply_weighted_acceptance(evidence, precision)
    return _weighted_summary(diagnostics, evidence, precision, growth)
