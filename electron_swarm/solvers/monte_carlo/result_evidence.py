"""Reaction-rate and population evidence assembly for internal Monte Carlo."""

from __future__ import annotations

import numpy as np

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.cross_sections import ProcessType
from electron_swarm.core.result_metadata import (
    TRANSPORT_MC_FIXED_POPULATION,
    TRANSPORT_MC_WEIGHTED_GROWTH,
)
from electron_swarm.core.results import RateResult
from electron_swarm.physics.kinetics import cell_integrated_rate_coefficient
from electron_swarm.solvers.monte_carlo.collisions import (
    _ProjectedProcess,
    _energy_loss_eV,
)
from electron_swarm.solvers.monte_carlo.reaction_rates import ReactionRateEstimate


def _weighted_branching_active(config: SwarmConfig, population_model: str) -> bool:
    return (
        population_model == "weighted_branching"
        and config.physics.ionization.energy_sharing != "loss_only"
    )


def _ionization_branching_model(config: SwarmConfig, population_model: str) -> str:
    if config.physics.ionization.energy_sharing == "loss_only":
        return "none"
    if population_model == "weighted_branching":
        return "two_daughter_weighted_resampling"
    return "single_daughter_sampling"


def _population_metadata(
    config: SwarmConfig, population_model: str
) -> dict[str, str | bool]:
    if _weighted_branching_active(config, population_model):
        return {
            "swarm_population_treatment": "weighted_branching_resampled",
            "nonconservative_growth_treatment": "explicit_weighted_branching",
            "secondary_electron_tracking": True,
            "ionization_branching_model": _ionization_branching_model(
                config, population_model
            ),
            "transport_definition": TRANSPORT_MC_WEIGHTED_GROWTH,
            "transport_has_bulk": False,
        }
    return {
        "swarm_population_treatment": "fixed_population_no_growth",
        "nonconservative_growth_treatment": "not_tracked",
        "secondary_electron_tracking": False,
        "ionization_branching_model": _ionization_branching_model(
            config, population_model
        ),
        "transport_definition": TRANSPORT_MC_FIXED_POPULATION,
        "transport_has_bulk": False,
    }


def _rate_results(
    config: SwarmConfig,
    projected: list[_ProjectedProcess],
    estimates: list[ReactionRateEstimate],
    energy: np.ndarray,
    eedf: np.ndarray,
    widths: np.ndarray,
    case_id: str,
    e_over_n_Td: float,
    gas_number_density_m3: float,
    event_evidence_estimates: list[ReactionRateEstimate] | None = None,
) -> tuple[list[RateResult], float, list[dict[str, object]]]:
    if len(projected) != len(estimates):
        raise ValueError("trajectory reaction-rate estimates do not match processes")
    density = float(gas_number_density_m3)
    rates: list[RateResult] = []
    reaction_metadata: list[dict[str, object]] = []
    ion_rate = 0.0
    if event_evidence_estimates is not None and len(event_evidence_estimates) != len(
        estimates
    ):
        raise ValueError("reaction event evidence does not match rate estimates")
    event_keys = (
        "event_count",
        "weighted_event_count",
        "event_sampling_enabled",
        "event_observation_status",
        "event_count_estimator",
        "event_count_rate_coefficient_m3_s",
        "zero_event_rate_upper_95_m3_s",
        "zero_event_confidence",
        "zero_event_upper_bound_method",
        "event_observation_residence_time_s",
    )
    for process_index, (item, estimate) in enumerate(
        zip(projected, estimates, strict=True)
    ):
        proc = item.process
        frac = float(item.fraction)
        k = float(estimate.rate_coefficient_m3_s)
        histogram_rate = cell_integrated_rate_coefficient(
            energy,
            widths,
            eedf,
            proc,
        )
        kmix = frac * k
        if proc.process_type == ProcessType.IONIZATION:
            ion_rate += kmix
        rates.append(
            RateResult(
                solver="monte_carlo",
                case_id=case_id,
                e_over_n_Td=e_over_n_Td,
                species=proc.species,
                process=proc.process,
                process_type=proc.process_type.value,
                threshold_eV=proc.threshold_eV,
                rate_coefficient_m3_s=k,
                target_species_fraction=frac,
                energy_loss_eV=_energy_loss_eV(proc.process_type, proc.threshold_eV),
                gas_number_density_m3=density,
            )
        )
        estimate_metadata = estimate.metadata()
        if event_evidence_estimates is not None:
            evidence_metadata = event_evidence_estimates[process_index].metadata()
            estimate_metadata.update(
                {key: evidence_metadata[key] for key in event_keys}
            )
            estimate_metadata["event_evidence_source"] = (
                "independent_transport_ensemble"
            )
        else:
            estimate_metadata["event_evidence_source"] = "rate_ensemble"
        reaction_metadata.append(
            {
                "species": proc.species,
                "process": proc.process,
                "process_type": proc.process_type.value,
                "threshold_eV": proc.threshold_eV,
                "target_species_fraction": frac,
                **estimate_metadata,
                "histogram_convolution_rate_coefficient_m3_s": histogram_rate,
                "direct_minus_histogram_rate_coefficient_m3_s": (k - histogram_rate),
            }
        )
    net_freq = density * ion_rate
    return rates, net_freq, reaction_metadata


def _weighted_growth_rate_consistency(
    *,
    projected: list[_ProjectedProcess],
    estimates: list[ReactionRateEstimate],
    event_estimates: list[ReactionRateEstimate] | None = None,
    gas_number_density_m3: float,
    elapsed_time_s: float,
    cumulative_log_growth: float,
    secondary_electron_events: int,
) -> dict[str, object]:
    """Cross-check branching growth against two same-trajectory rate estimates."""

    if event_estimates is not None and len(event_estimates) != len(estimates):
        raise ValueError("weighted-growth event evidence does not match rates")
    ionization = [
        (item, estimate)
        for item, estimate in zip(projected, estimates, strict=True)
        if item.process.process_type == ProcessType.IONIZATION
    ]
    direct_frequency = float(gas_number_density_m3) * sum(
        float(item.fraction) * float(estimate.rate_coefficient_m3_s)
        for item, estimate in ionization
    )
    evidence_source = estimates if event_estimates is None else event_estimates
    event_ionization = [
        (item, estimate)
        for item, estimate in zip(projected, evidence_source, strict=True)
        if item.process.process_type == ProcessType.IONIZATION
    ]
    event_count = sum(int(estimate.event_count) for _, estimate in event_ionization)
    weighted_event_count = sum(
        float(estimate.weighted_event_count) for _, estimate in event_ionization
    )
    weighted_ensemble_events = any(
        estimate.event_observation_status
        == "weighted_ensemble_correlated_events_audit_only"
        for _, estimate in event_ionization
    )
    exposures = [
        float(
            estimate.weighted_residence_time_s
            if weighted_ensemble_events
            else estimate.event_observation_residence_time_s
        )
        for _, estimate in event_ionization
    ]
    event_exposure = exposures[0] if exposures else 0.0
    exposure_consistent = bool(
        not exposures
        or all(
            np.isclose(value, event_exposure, rtol=1.0e-12, atol=0.0)
            for value in exposures
        )
    )
    if exposures and (not exposure_consistent or event_exposure <= 0.0):
        event_frequency: float | None = None
        zero_event_upper: float | None = None
    elif exposures:
        event_frequency = float(
            (weighted_event_count if weighted_ensemble_events else event_count)
            / event_exposure
        )
        zero_event_upper = (
            None
            if weighted_ensemble_events
            else (float(-np.log(0.05) / event_exposure) if event_count == 0 else None)
        )
    else:
        event_frequency = 0.0
        zero_event_upper = 0.0
    elapsed = float(elapsed_time_s)
    population_frequency = (
        float(cumulative_log_growth) / elapsed if elapsed > 0.0 else None
    )
    return {
        "model": "explicit_ionization_branching_no_attachment",
        "production_elapsed_time_s": elapsed,
        "cumulative_log_growth": float(cumulative_log_growth),
        "population_log_growth_frequency_s_inv": population_frequency,
        "direct_ionization_frequency_s_inv": direct_frequency,
        "event_evidence_model": (
            "weighted_ensemble_correlated_events"
            if weighted_ensemble_events
            else "raw_macro_poisson"
        ),
        "event_count_ionization_frequency_s_inv": event_frequency,
        "ionization_event_count": int(event_count),
        "weighted_ionization_event_count": float(weighted_event_count),
        "secondary_electron_events": int(secondary_electron_events),
        "event_count_matches_secondary": bool(
            event_count == int(secondary_electron_events)
        ),
        "event_observation_residence_time_s": float(event_exposure),
        "event_exposure_consistent_across_processes": exposure_consistent,
        "zero_event_frequency_upper_95_s_inv": zero_event_upper,
    }
