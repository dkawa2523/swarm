"""Observable finalization and result assembly for one Monte Carlo case."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any

import numpy as np

import electron_swarm.solvers.monte_carlo.evidence as _mc_evidence
import electron_swarm.solvers.monte_carlo.result_evidence as _result_evidence
from electron_swarm.core.results import RateResult, SwarmCaseResult
from electron_swarm.core.transport import ElectronTransport
from electron_swarm.physics.angular_scattering import expected_angular_metadata
from electron_swarm.solvers.monte_carlo.case_state import _MonteCarloCaseState
from electron_swarm.solvers.monte_carlo.compiled_kernel import (
    COMPILED_MC_KERNEL_SCHEMA_VERSION,
)
from electron_swarm.solvers.monte_carlo.histogram import (
    _build_eedf_histogram,
    _mc_tail_comparison_status,
    _mc_tail_uncertainty_metadata_from_effective_counts,
)
from electron_swarm.solvers.monte_carlo.reaction_rates import ReactionRateEstimate
from electron_swarm.solvers.monte_carlo.setup import _MonteCarloRunSetup
from electron_swarm.solvers.monte_carlo.weighted_transport import (
    SynchronizedWeightedGrowthFluxObserver,
)


logger = logging.getLogger("electron_swarm.solvers.monte_carlo.solver")


def _tail_refinement_treatment(
    setup: _MonteCarloRunSetup,
    state: _MonteCarloCaseState,
) -> str:
    if setup.tail_collisions <= 0:
        return "disabled"
    if state.tail_sampling_used:
        return "executed"
    return "configured_not_triggered"


def _eedf_estimator_schema(setup: _MonteCarloRunSetup) -> str:
    if setup.magnetic_enabled:
        return _mc_evidence.MC_MAGNETIC_EEDF_ESTIMATOR_SCHEMA_VERSION
    return _mc_evidence.MC_EEDF_ESTIMATOR_SCHEMA_VERSION


@dataclass(frozen=True, slots=True)
class _MonteCarloCaseStatistics:
    energy: np.ndarray
    widths: np.ndarray
    eedf: np.ndarray
    counts: np.ndarray
    effective_counts: np.ndarray
    nonzero_bins: np.ndarray
    max_nonzero_energy: float
    mean_energy: float
    mean_time: float
    drift_velocity: float
    mobility: float
    diffusion_l: float | None
    diffusion_t: float | None
    energy_mobility: float | None
    energy_diffusion_l: float | None
    energy_diffusion_t: float | None
    transport_diagnostics: dict[str, Any]
    reaction_rate_estimates: list[ReactionRateEstimate]


@dataclass(frozen=True, slots=True)
class _MonteCarloCaseEvidence:
    rates: list[RateResult]
    net_frequency: float
    effective_townsend: float
    population_metadata: dict[str, Any]
    mc_run_details: dict[str, Any]
    diagnostics: dict[str, Any]


def _finalize_case_statistics(
    setup: _MonteCarloRunSetup,
    state: _MonteCarloCaseState,
) -> _MonteCarloCaseStatistics:
    state.audit.tracked_particle_final_energy_eV = (
        state.ensemble.total_weighted_energy_eV()
    )
    energy, widths, eedf, counts, effective_counts = _build_eedf_histogram(
        state.edges,
        state.counts,
        state.weighted_hist,
        state.weighted_square_hist,
    )
    nonzero_bins = counts > 0
    max_nonzero_energy = (
        float(np.max(energy[nonzero_bins])) if np.any(nonzero_bins) else 0.0
    )
    mean_energy = float(np.sum(energy * eedf * widths))
    mean_time = float(
        np.sum(state.ensemble.weights * state.ensemble.times)
        / max(float(np.sum(state.ensemble.weights)), 1.0e-300)
    )
    energy_mobility: float | None = None
    energy_diffusion_l: float | None = None
    energy_diffusion_t: float | None = None
    if setup.branching_active:
        assert isinstance(
            state.transport_observer,
            SynchronizedWeightedGrowthFluxObserver,
        )
    direct_estimate, transport_diagnostics = state.transport_observer.finalize()
    drift_velocity = direct_estimate.drift_velocity_m_s
    mobility = direct_estimate.mobility_m2_V_s
    diffusion_l = direct_estimate.diffusion_L_m2_s
    diffusion_t = direct_estimate.diffusion_T_m2_s
    if setup.magnetic_enabled:
        transport_diagnostics["energy_transport_status"] = (
            "unsupported_magnetic_tensor_not_identified"
        )
    else:
        energy_mobility = direct_estimate.energy_mobility_m2_V_s
        energy_diffusion_l = direct_estimate.energy_diffusion_L_m2_s
        energy_diffusion_t = direct_estimate.energy_diffusion_T_m2_s

    # Persist qualification-critical settings independently of optional audit.
    tail_treatment = _tail_refinement_treatment(setup, state)
    transport_diagnostics["mc_run_provenance"] = {
        "seed": setup.cfg.seed,
        "case_seed": int(state.case_seed),
        "case_seed_derivation": _mc_evidence.MC_CASE_SEED_DERIVATION,
        "solver_source_sha256": setup.solver_source_sha256,
        "particles": int(setup.particles),
        "warmup_collisions": int(setup.warmup_collisions),
        "production_collisions": int(setup.collisions),
        "tail_max_collisions": int(setup.tail_collisions),
        "tail_rate_rse_trigger": float(setup.cfg.tail_rate_rse_trigger),
        "tail_collisions_executed": (
            int(setup.tail_collisions) if state.tail_sampling_used else 0
        ),
        "tail_refinement_treatment": tail_treatment,
        "population_model": setup.population_model,
        "gas_number_density_m3": float(setup.density),
        "electric_field_V_m": float(state.electric_field_scalar_V_m),
        "trial_collision_frequency_s_inv": float(setup.trial_frequency),
        "collision_clock": "global",
        "numeric_kernel_requested": setup.cfg.numeric_kernel,
        "numeric_kernel_used": (
            "numba" if setup.compiled_plan is not None else "python"
        ),
        "numeric_kernel_schema_version": (
            COMPILED_MC_KERNEL_SCHEMA_VERSION
            if setup.compiled_plan is not None
            else None
        ),
        "numeric_kernel_fallback_reason": setup.compiled_fallback_reason,
        "transport_estimator": setup.cfg.transport_estimator,
        "eedf_estimator_schema_version": (
            _eedf_estimator_schema(setup)
        ),
        "transport_field_legs": (
            2 if setup.cfg.transport_estimator == "paired_field_parity" else 1
        ),
        "transport_budget_interpretation": "configured_budget_per_field_leg",
        "synchronization_frequency_s_inv": float(setup.barrier_frequency),
        "barrier_cadence_s": float(setup.barrier_dt),
        "magnetic_enabled": setup.magnetic_enabled,
        "magnetic_B_T": setup.effective_B_T,
        "magnetic_angle_EB_deg": float(setup.magnetic.angle_EB_deg),
        "angular_scattering_model": setup.angular_name,
        "elastic_collision_clock_model": (
            "marked_species_maxwellian_relative_speed_global"
        ),
        "elastic_collision_kinematics": "exact_center_of_mass_binary_collision",
        "neutral_thermal_motion_model": (
            "maxwellian_relative_speed_exact_binary_collision"
        ),
        "gas_temperature_K": float(setup.config.conditions.gas_temperature_K),
        "inelastic_target_motion_model": "stationary_target_lab_energy",
        "ionization_source_model": setup.config.physics.ionization.energy_sharing,
        "high_energy_extrapolation": (
            setup.config.cross_sections.high_energy_extrapolation
        ),
        "max_energy_limit_eV": float(setup.max_energy_limit),
    }
    reaction_rate_estimates = state.reaction_rate_accumulator.estimates()
    if setup.branching_active:
        transport_diagnostics["population_growth_consistency"] = (
            _result_evidence._weighted_growth_rate_consistency(
                projected=setup.projected,
                estimates=reaction_rate_estimates,
                event_estimates=state.event_evidence_estimates,
                gas_number_density_m3=setup.density,
                elapsed_time_s=(
                    state.production_mean_time
                    if state.production_mean_time is not None
                    else mean_time
                ),
                cumulative_log_growth=state.cumulative_log_growth,
                secondary_electron_events=state.production_secondary_events,
            )
        )
    return _MonteCarloCaseStatistics(
        energy=energy,
        widths=widths,
        eedf=eedf,
        counts=counts,
        effective_counts=effective_counts,
        nonzero_bins=nonzero_bins,
        max_nonzero_energy=max_nonzero_energy,
        mean_energy=mean_energy,
        mean_time=mean_time,
        drift_velocity=drift_velocity,
        mobility=mobility,
        diffusion_l=diffusion_l,
        diffusion_t=diffusion_t,
        energy_mobility=energy_mobility,
        energy_diffusion_l=energy_diffusion_l,
        energy_diffusion_t=energy_diffusion_t,
        transport_diagnostics=transport_diagnostics,
        reaction_rate_estimates=reaction_rate_estimates,
    )


def _assemble_case_evidence(
    setup: _MonteCarloRunSetup,
    state: _MonteCarloCaseState,
    statistics: _MonteCarloCaseStatistics,
) -> _MonteCarloCaseEvidence:
    rates, net_frequency, reaction_rate_metadata = _result_evidence._rate_results(
        setup.config,
        setup.projected,
        statistics.reaction_rate_estimates,
        statistics.energy,
        statistics.eedf,
        statistics.widths,
        state.case_id,
        state.e_over_n_Td,
        gas_number_density_m3=setup.density,
        event_evidence_estimates=state.event_evidence_estimates,
    )
    effective_townsend = net_frequency / max(
        abs(statistics.drift_velocity) * setup.density,
        1.0e-300,
    )
    population_metadata = _result_evidence._population_metadata(
        setup.config,
        setup.population_model,
    )
    residence_integration = (
        "endpoint_average"
        if setup.magnetic_enabled
        else "exact_bin_crossings_and_piecewise_linear_rate_integrals"
    )
    mc_run_details = {
        "adapter": "internal_monte_carlo",
        "monte_carlo_backend": "internal",
        "mc_seed": setup.cfg.seed,
        "mc_case_seed": int(state.case_seed),
        "mc_case_seed_derivation": _mc_evidence.MC_CASE_SEED_DERIVATION,
        "mc_solver_source_sha256": setup.solver_source_sha256,
        "mc_particles": int(setup.particles),
        "mc_population_model": setup.population_model,
        "mc_collision_clock": "global",
        "mc_tail_max_collisions": int(setup.tail_collisions),
        "mc_tail_rate_rse_trigger": float(setup.cfg.tail_rate_rse_trigger),
        "mc_tail_collisions_executed": (
            int(setup.tail_collisions) if state.tail_sampling_used else 0
        ),
        "tail_refinement_treatment": _tail_refinement_treatment(setup, state),
        "mc_numeric_kernel_requested": setup.cfg.numeric_kernel,
        "mc_numeric_kernel_used": (
            "numba" if setup.compiled_plan is not None else "python"
        ),
        "mc_numeric_kernel_schema_version": (
            COMPILED_MC_KERNEL_SCHEMA_VERSION
            if setup.compiled_plan is not None
            else None
        ),
        "mc_numeric_kernel_fallback_reason": setup.compiled_fallback_reason,
        "mc_transport_estimator": setup.cfg.transport_estimator,
        "mc_transport_field_legs": (
            2 if setup.cfg.transport_estimator == "paired_field_parity" else 1
        ),
        "mc_transport_budget_interpretation": "configured_budget_per_field_leg",
        "mc_synchronization_frequency_s_inv": float(setup.barrier_frequency),
        "magnetic_field_treatment": (
            "boris_lorentz_push" if setup.magnetic_enabled else "none"
        ),
        "magnetic_field_B_T": setup.effective_B_T,
        "magnetic_field_angle_EB_deg": float(setup.magnetic.angle_EB_deg),
        "gas_number_density_m-3": float(setup.density),
        "electric_field_V_m": float(state.electric_field_scalar_V_m),
        "reaction_rates_source": "trajectory_time_average_sigma_v",
        "elastic_collision_clock_model": (
            "marked_species_maxwellian_relative_speed_global"
        ),
        "elastic_collision_kinematics": "exact_center_of_mass_binary_collision",
        "neutral_thermal_motion_model": (
            "maxwellian_relative_speed_exact_binary_collision"
        ),
        "gas_temperature_K": float(setup.config.conditions.gas_temperature_K),
        "inelastic_target_motion_model": "stationary_target_lab_energy",
        "ionization_source_model": setup.config.physics.ionization.energy_sharing,
        "ionization_source_treatment": (setup.config.physics.ionization.energy_sharing),
        "eedf_estimator": (
            f"pooled_main_and_tail_time_residence_{residence_integration}"
            if state.tail_sampling_used
            else f"time_residence_{residence_integration}"
        ),
        "eedf_estimator_schema_version": (
            _eedf_estimator_schema(setup)
        ),
        "tail_sampling_model": (
            "reaction_kernel_weighted_ensemble"
            if state.tail_sampling_used
            else (
                "ordinary_trajectory_sampling_resolved"
                if setup.tail_strata_edges is not None
                else "ordinary_trajectory_sampling"
            )
        ),
        "tail_estimator_schema_version": (
            _mc_evidence.WEIGHTED_MC_TAIL_ESTIMATOR_SCHEMA_VERSION
            if setup.tail_strata_edges is not None
            else None
        ),
        "tail_sampling_strata_edges_eV": (
            [float(value) for value in setup.tail_strata_edges]
            if setup.tail_strata_edges is not None
            else []
        ),
        **population_metadata,
    }
    diagnostics: dict[str, Any] = {
        "internal_monte_carlo_transport": statistics.transport_diagnostics,
        "internal_monte_carlo_reaction_rates": {
            "estimator": "trajectory_time_average_sigma_v",
            "sampling_model": mc_run_details["tail_sampling_model"],
            "sampling_strata_edges_eV": mc_run_details["tail_sampling_strata_edges_eV"],
            "stratum_weight_conservation": (
                "exact_per_occupied_stratum"
                if state.tail_sampling_used
                else "not_applicable"
            ),
            "raw_event_count_role": (
                "independent_transport_ensemble_poisson_evidence"
                if state.tail_sampling_used
                else "independent_poisson_evidence"
            ),
            "energy_loss_estimator": (
                "process_specific_residence_or_exact_thermal_event"
            ),
            "rates": reaction_rate_metadata,
        },
    }
    statistics.transport_diagnostics["tail_sampling"] = {
        "state": mc_run_details["tail_refinement_treatment"],
        "estimator_schema_version": mc_run_details["tail_estimator_schema_version"],
        "model": mc_run_details["tail_sampling_model"],
        "strata_edges_eV": mc_run_details["tail_sampling_strata_edges_eV"],
        "resampling_cadence": (
            "independent_rate_ensemble_fixed_barriers"
            if state.tail_sampling_used
            else "none"
        ),
        "resampling_interval_barriers": setup.tail_resample_interval,
        "configured_max_collisions": int(setup.tail_collisions),
        "executed_collisions": (
            int(setup.tail_collisions) if state.tail_sampling_used else 0
        ),
        "activation_rate_rse_trigger": float(setup.cfg.tail_rate_rse_trigger),
        "transport_ensemble_resampled_for_tail": False,
        "reported_eedf_includes_main_production": True,
        "reported_rates_include_main_production": True,
        "stratum_weight_conservation": (
            "exact_per_occupied_stratum"
            if state.tail_sampling_used
            else "not_applicable"
        ),
    }
    return _MonteCarloCaseEvidence(
        rates=rates,
        net_frequency=net_frequency,
        effective_townsend=effective_townsend,
        population_metadata=population_metadata,
        mc_run_details=mc_run_details,
        diagnostics=diagnostics,
    )


def _attach_audit_diagnostics_and_log(
    setup: _MonteCarloRunSetup,
    state: _MonteCarloCaseState,
    statistics: _MonteCarloCaseStatistics,
    evidence: _MonteCarloCaseEvidence,
) -> None:
    if setup.collect_audit:
        state.run_audit.set_population_state(
            state.ensemble.weights,
            statistics.mean_time,
            state.audit_cumulative_log_growth,
        )
        audit_metadata = state.audit.as_metadata()
        run_audit_metadata = state.run_audit.as_metadata()
        tail_uncertainty = _mc_tail_uncertainty_metadata_from_effective_counts(
            statistics.energy,
            statistics.effective_counts,
            setup.tail_threshold,
            bin_probability=statistics.eedf * statistics.widths,
        )
        tail_comparison_status = _mc_tail_comparison_status(
            energy_balance_status=str(audit_metadata["mc_energy_balance_status"]),
            tail_uncertainty_status=str(tail_uncertainty["mc_tail_uncertainty_status"]),
        )
        evidence.diagnostics["internal_monte_carlo_audit"] = {
            "mc_run": {
                **run_audit_metadata,
                **evidence.mc_run_details,
                "mc_samples": int(np.sum(statistics.counts)),
                "mc_nonzero_bins": int(np.count_nonzero(statistics.nonzero_bins)),
                "mc_max_nonzero_energy_eV": statistics.max_nonzero_energy,
                "mc_energy_bin_max_eV": float(state.edges[-1]),
            },
            "mc_energy_audit": audit_metadata,
            "mc_tail_audit": {
                **tail_uncertainty,
                "mc_tail_comparison_status": tail_comparison_status,
            },
        }
        direct_event_loss = float(
            sum(
                estimate.energy_loss.event_weighted_energy_loss_eV
                for estimate in statistics.reaction_rate_estimates
            )
        )
        event_complete = all(
            estimate.event_count == estimate.energy_loss.event_loss_sample_count
            for estimate in statistics.reaction_rate_estimates
            if estimate.event_sampling_enabled
        )
        audit_physical_loss = float(
            audit_metadata["mc_elastic_energy_loss_eV"]
            + audit_metadata["mc_inelastic_energy_loss_eV"]
            + audit_metadata["mc_ionization_threshold_loss_eV"]
        )
        loss_difference = direct_event_loss - audit_physical_loss
        loss_scale = max(
            abs(direct_event_loss),
            abs(audit_metadata["mc_elastic_energy_loss_eV"])
            + abs(audit_metadata["mc_inelastic_energy_loss_eV"])
            + abs(audit_metadata["mc_ionization_threshold_loss_eV"]),
            1.0e-300,
        )
        matched = event_complete and abs(loss_difference) <= 1.0e-12 * loss_scale
        evidence.diagnostics["internal_monte_carlo_reaction_rates"][
            "energy_audit_consistency"
        ] = {
            "status": "matched" if matched else "mismatch",
            "direct_event_physical_energy_loss_eV": direct_event_loss,
            "mc_energy_audit_physical_collision_loss_eV": audit_physical_loss,
            "difference_eV": float(loss_difference),
        }
        logger.info(
            "internal MC case %s: base_seed=%s case_seed=%d "
            "particles=%d samples=%d "
            "mean_energy_eV=%.6g max_sampled_energy_eV=%.6g "
            "xs_above_fraction=%.3g",
            state.case_id,
            setup.cfg.seed,
            state.case_seed,
            setup.particles,
            int(run_audit_metadata["mc_histogram_samples"] or 0),
            statistics.mean_energy,
            float(run_audit_metadata["mc_max_sampled_energy_eV"] or 0.0),
            float(run_audit_metadata["mc_energy_samples_above_xs_max_fraction"] or 0.0),
        )
    else:
        logger.info(
            "internal MC case %s: base_seed=%s case_seed=%d particles=%d "
            "samples=%d mean_energy_eV=%.6g",
            state.case_id,
            setup.cfg.seed,
            state.case_seed,
            setup.particles,
            int(np.sum(statistics.counts)),
            statistics.mean_energy,
        )


def _assemble_swarm_case_result(
    setup: _MonteCarloRunSetup,
    state: _MonteCarloCaseState,
    statistics: _MonteCarloCaseStatistics,
    evidence: _MonteCarloCaseEvidence,
) -> SwarmCaseResult:
    metadata = {
        "magnetic_field_treatment": (
            "boris_lorentz_push" if setup.magnetic_enabled else "none"
        ),
        "magnetic_field_B_T": setup.effective_B_T,
        "magnetic_field_angle_EB_deg": float(setup.magnetic.angle_EB_deg),
        "ionization_source_model": setup.config.physics.ionization.energy_sharing,
        "ionization_source_treatment": setup.config.physics.ionization.energy_sharing,
        "tail_refinement_treatment": _tail_refinement_treatment(setup, state),
        "monte_carlo_base_seed": int(setup.cfg.seed),
        "monte_carlo_case_seed": int(state.case_seed),
        "transport_definition": evidence.population_metadata["transport_definition"],
        **expected_angular_metadata(setup.config),
    }
    transport = ElectronTransport.from_actual(
        definition=str(evidence.population_metadata["transport_definition"]),
        gas_number_density_m3=setup.density,
        drift_velocity_m_s=statistics.drift_velocity,
        mobility_m2_V_s=statistics.mobility,
        diffusion_L_m2_s=statistics.diffusion_l,
        diffusion_T_m2_s=statistics.diffusion_t,
        electron_energy_mobility_m2_V_s=statistics.energy_mobility,
        electron_energy_diffusion_L_m2_s=statistics.energy_diffusion_l,
        electron_energy_diffusion_T_m2_s=statistics.energy_diffusion_t,
    )
    return SwarmCaseResult(
        solver="monte_carlo",
        case_id=state.case_id,
        e_over_n_Td=state.e_over_n_Td,
        mean_energy_eV=statistics.mean_energy,
        net_ionization_frequency_s=evidence.net_frequency,
        effective_townsend_m2=evidence.effective_townsend,
        transport=transport,
        energy_eV=statistics.energy,
        eedf=statistics.eedf,
        energy_widths_eV=statistics.widths,
        eedf_counts=statistics.counts.astype(int),
        eedf_effective_counts=statistics.effective_counts,
        rates=evidence.rates,
        metadata=metadata,
        diagnostics=evidence.diagnostics,
    )


def _build_monte_carlo_case_result(
    setup: _MonteCarloRunSetup,
    state: _MonteCarloCaseState,
) -> SwarmCaseResult:
    statistics = _finalize_case_statistics(setup, state)
    evidence = _assemble_case_evidence(setup, state, statistics)
    _attach_audit_diagnostics_and_log(setup, state, statistics, evidence)
    return _assemble_swarm_case_result(setup, state, statistics, evidence)
