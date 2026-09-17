"""Weighted-branching phase runtime for the internal Monte Carlo solver."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.constants import AMU_KG, E_CHARGE_C, ELECTRON_MASS_KG
from electron_swarm.core.cross_sections import ProcessType, gas_mass_amu
from electron_swarm.core.solver_configs import MonteCarloAdapterConfig
import electron_swarm.solvers.monte_carlo.audits as _audits
import electron_swarm.solvers.monte_carlo.population as _population
from electron_swarm.solvers.monte_carlo.collisions import (
    SCATTERING_TYPES,
    TRAJECTORY_REACTION_TYPES,
    PreparedCollisionSampler,
    _PostReactionOutcome,
    _ProjectedProcess,
    _ionization_daughters,
    _is_scattering_event_process,
    _moment_cross_sections,
    _post_reaction_outcome,
    _thermal_elastic_velocity_after_collision,
)
from electron_swarm.solvers.monte_carlo.compiled_kernel import (
    IONIZATION_EQUAL,
    IONIZATION_PRIMARY_SECONDARY,
    KERNEL_CROSS_SECTION_RANGE_EXCEEDED,
    KERNEL_ENERGY_LIMIT_EXCEEDED,
    KERNEL_ERROR_EXTRAPOLATION,
    KERNEL_FLIGHT_INTEGRATION_FAILED,
    KERNEL_INVALID_RESIDENCE_SAMPLE,
    KERNEL_MAJORANT_EXCEEDED,
    KERNEL_OK,
    PROCESS_EXCITATION,
    PROCESS_IONIZATION,
    PROCESS_NONE,
    PROCESS_SCATTERING,
    PROCESS_SUPERELASTIC,
    advance_weighted_barrier,
)
from electron_swarm.solvers.monte_carlo.cross_section_table import NUMBA_AVAILABLE
from electron_swarm.solvers.monte_carlo.kinematics import (
    energy_from_velocity_eV,
    random_direction,
    speed_from_energy_m_s,
)
from electron_swarm.solvers.monte_carlo.orbit import _advance_to_trial_event
from electron_swarm.solvers.monte_carlo.reaction_rates import (
    TrajectoryReactionRateAccumulator,
)
from electron_swarm.solvers.monte_carlo.weighted_transport import (
    SynchronizedWeightedGrowthFluxObserver,
    WeightedGrowthTransportLagPlan,
)


@dataclass(frozen=True, slots=True)
class _CompiledKernelPlan:
    process_kind: np.ndarray
    process_threshold_eV: np.ndarray
    process_target_mass_fraction: np.ndarray
    ionization_model: int


def _compiled_kernel_plan(
    *,
    config: SwarmConfig,
    solver_config: MonteCarloAdapterConfig,
    projected: list[_ProjectedProcess],
    angular_name: str,
    branching_active: bool,
    magnetic_enabled: bool,
    collect_audit: bool,
) -> tuple[_CompiledKernelPlan | None, str | None]:
    """Resolve the compiled kernel without changing requested physics."""

    if solver_config.numeric_kernel == "python":
        return None, "python_kernel_requested"
    unsupported: list[str] = []
    if not NUMBA_AVAILABLE:
        unsupported.append("numba_unavailable")
    if not branching_active:
        unsupported.append("population_model_not_weighted_branching")
    if magnetic_enabled:
        unsupported.append("magnetic_orbit")
    if angular_name != "isotropic":
        unsupported.append("angular_model_not_isotropic")
    if config.physics.ionization.energy_sharing not in {
        "equal",
        "primary_secondary",
    }:
        unsupported.append("ionization_energy_sharing")
    # The detailed audit records every individual numerical event.  Running
    # its reference path is intentional and is reported in provenance.
    if collect_audit:
        unsupported.append("detailed_event_audit")
    if unsupported:
        reason = ",".join(unsupported)
        if solver_config.numeric_kernel == "numba" and not collect_audit:
            raise NotImplementedError(
                "monte_carlo numeric_kernel=numba does not support this execution: "
                f"{reason}"
            )
        return None, reason

    kinds: list[int] = []
    thresholds: list[float] = []
    target_mass_fractions: list[float] = []
    for item in projected:
        process = item.process
        if _is_scattering_event_process(item, projected, angular_name):
            kind = PROCESS_SCATTERING
        elif process.process_type == ProcessType.EXCITATION:
            kind = PROCESS_EXCITATION
        elif process.process_type == ProcessType.IONIZATION:
            kind = PROCESS_IONIZATION
        elif process.process_type == ProcessType.SUPERELASTIC:
            kind = PROCESS_SUPERELASTIC
        else:
            kind = PROCESS_NONE
        kinds.append(kind)
        thresholds.append(float(process.threshold_eV or 0.0))
        if kind == PROCESS_SCATTERING:
            mass_amu = process.mass_amu or gas_mass_amu(
                config.conditions,
                process.species,
            )
            target_mass = float(mass_amu) * AMU_KG
            target_mass_fractions.append(target_mass / (ELECTRON_MASS_KG + target_mass))
        else:
            target_mass_fractions.append(0.0)
    ionization_model = (
        IONIZATION_EQUAL
        if config.physics.ionization.energy_sharing == "equal"
        else IONIZATION_PRIMARY_SECONDARY
    )
    return (
        _CompiledKernelPlan(
            process_kind=np.asarray(kinds, dtype=np.int64),
            process_threshold_eV=np.asarray(thresholds, dtype=float),
            process_target_mass_fraction=np.asarray(
                target_mass_fractions,
                dtype=float,
            ),
            ionization_model=ionization_model,
        ),
        None,
    )


def _sample_weighted_branching_collision(
    *,
    config: SwarmConfig,
    ensemble: _population._ParticleEnsemble,
    particle_index: int,
    energy_eV: float,
    projected: list[_ProjectedProcess],
    collision_sampler: PreparedCollisionSampler,
    angular_name: str,
    angular_model: object,
    rng: np.random.Generator,
    audit: _audits._EnergyAudit | _audits._NullEnergyAudit,
    run_audit: _audits._MonteCarloRunAudit | _audits._NullMonteCarloRunAudit,
    reaction_rate_accumulator: TrajectoryReactionRateAccumulator,
    reaction_rate_block_index: int | None,
    sampling_enabled: bool,
) -> bool:
    """Apply one null/physical event in the synchronized branching kernel."""

    particle = int(particle_index)
    weight = float(ensemble.weights[particle])
    trial = collision_sampler.sample_trial(
        ensemble.velocities[particle],
        rng,
    )
    audit.record_trial_ratio(trial.acceptance_probability)
    if trial.process_index is None:
        if sampling_enabled:
            run_audit.record_null_collision()
        return False

    choice = int(trial.process_index)
    proc = projected[choice].process
    if sampling_enabled:
        run_audit.record_collision(proc.process_type)
    event_energy_loss_eV: float | None = None
    if proc.process_type in SCATTERING_TYPES:
        if trial.target_velocity_m_s is None or trial.relative_energy_eV is None:
            raise RuntimeError("thermal elastic collision lacks target state")
        relative_energy_eV = float(trial.relative_energy_eV)
        if angular_name == "maxent_p1":
            total, momentum = _moment_cross_sections(
                projected,
                proc.species,
                relative_energy_eV,
            )
            mu = float(
                angular_model.sample_mu(
                    relative_energy_eV,
                    rng,
                    sigma_total=total,
                    sigma_momentum=momentum,
                )
            )
        else:
            mu = float(angular_model.sample_mu(relative_energy_eV, rng))
        post_velocity = _thermal_elastic_velocity_after_collision(
            config,
            proc,
            ensemble.velocities[particle],
            trial.target_velocity_m_s,
            mu,
            rng,
        )
        post_energy = float(energy_from_velocity_eV(post_velocity))
        if post_energy > collision_sampler.max_energy_eV * (1.0 + 1.0e-12):
            raise RuntimeError(
                "internal_monte_carlo thermal collision exceeded "
                "physics.energy_grid_policy.max_eV_limit"
            )
        event_energy_loss_eV = float(energy_eV - post_energy)
        audit.record_elastic_collision(energy_eV, post_energy, weight)
        ensemble.velocities[particle] = post_velocity
    elif proc.process_type in TRAJECTORY_REACTION_TYPES:
        if proc.process_type == ProcessType.IONIZATION:
            daughters = _ionization_daughters(
                config,
                proc.threshold_eV,
                energy_eV,
            )
            outcome = _PostReactionOutcome(
                tracked_energy_eV=float(sum(daughters.energies_eV)),
                ionization_threshold_loss_eV=daughters.threshold_loss_eV,
            )
            event_energy_loss_eV = float(daughters.threshold_loss_eV)
            audit.record_reaction(outcome, weight)
            energies = daughters.energies_eV
            ensemble.velocities[particle] = random_direction(rng) * speed_from_energy_m_s(
                max(energies[0], 1.0e-4)
            )
            if len(energies) > 1:
                ensemble.append_particle(
                    position=ensemble.positions[particle].copy(),
                    velocity=(
                        random_direction(rng)
                        * speed_from_energy_m_s(max(energies[1], 1.0e-4))
                    ),
                    time_s=float(ensemble.times[particle]),
                    weight=weight,
                    lineage=int(ensemble.lineages[particle]),
                )
                if sampling_enabled:
                    run_audit.record_secondary_electron()
        else:
            outcome = _post_reaction_outcome(
                config,
                proc.process_type,
                proc.threshold_eV,
                energy_eV,
                rng,
            )
            event_energy_loss_eV = float(
                outcome.inelastic_energy_loss_eV + outcome.ionization_threshold_loss_eV
            )
            audit.record_reaction(outcome, weight)
            ensemble.velocities[particle] = random_direction(rng) * speed_from_energy_m_s(
                outcome.tracked_energy_eV
            )

    if sampling_enabled:
        if event_energy_loss_eV is None or reaction_rate_block_index is None:
            raise RuntimeError(
                "sampled weighted-growth collision lacks physical rate evidence"
            )
        reaction_rate_accumulator.record_event(
            choice,
            weight,
            energy_loss_eV=event_energy_loss_eV,
            block_index=reaction_rate_block_index,
        )
    return bool(
        proc.process_type == ProcessType.IONIZATION
        and config.physics.ionization.energy_sharing != "loss_only"
    )


def _advance_compiled_weighted_barrier(
    *,
    config: SwarmConfig,
    ensemble: _population._ParticleEnsemble,
    barrier_time_s: float,
    trial_frequency: float,
    density: float,
    electric_acceleration_m_s2: np.ndarray,
    collision_sampler: PreparedCollisionSampler,
    compiled_plan: _CompiledKernelPlan,
    rng: np.random.Generator,
    edges: np.ndarray,
    counts: np.ndarray,
    weighted_hist: np.ndarray,
    weighted_square_hist: np.ndarray,
    reaction_rate_accumulator: TrajectoryReactionRateAccumulator,
    reaction_rate_block_index: int,
    transport_observer: SynchronizedWeightedGrowthFluxObserver | None,
    sampling_enabled: bool,
    zero_high_energy_policy: bool,
    max_cross_section_energy_eV: float,
    max_energy_limit_eV: float,
) -> int:
    """Run one dense barrier and commit its sufficient statistics."""

    collision_table = collision_sampler.table
    rate_table = reaction_rate_accumulator._prepared_cross_sections
    storage = ensemble.compiled_storage()
    result = advance_weighted_barrier(
        *storage,
        float(barrier_time_s),
        float(trial_frequency),
        float(density),
        electric_acceleration_m_s2,
        collision_table.grid_eV,
        collision_table.values_m2,
        collision_table.slopes_m2_eV,
        collision_table.process_min_eV,
        collision_table.process_max_eV,
        collision_table.right_values_m2,
        collision_table.error_extrapolation,
        collision_sampler.process_collision_group,
        float(collision_sampler.cold_bound_m3_s),
        collision_sampler.thermal_group_bounds_m3_s,
        collision_sampler.thermal_group_component_std_m_s,
        collision_sampler.thermal_group_mean_speed_m_s,
        collision_sampler.thermal_group_proposal_mode,
        compiled_plan.process_kind,
        compiled_plan.process_threshold_eV,
        compiled_plan.process_target_mass_fraction,
        compiled_plan.ionization_model,
        float(config.physics.ionization.secondary_electron_energy_eV),
        bool(sampling_enabled),
        bool(sampling_enabled and transport_observer is not None),
        int(reaction_rate_block_index),
        edges,
        counts,
        weighted_hist,
        weighted_square_hist,
        reaction_rate_accumulator._flight_integration_breakpoints_eV,
        rate_table.grid_eV,
        rate_table.values_m2,
        rate_table.slopes_m2_eV,
        rate_table.process_min_eV,
        rate_table.process_max_eV,
        rate_table.right_values_m2,
        rate_table.error_extrapolation,
        reaction_rate_accumulator._constant_energy_loss_eV,
        reaction_rate_accumulator._recoil_mass_ratio,
        reaction_rate_accumulator._weighted_time,
        reaction_rate_accumulator._integrals,
        reaction_rate_accumulator._energy_loss_integrals,
        reaction_rate_accumulator._event_counts,
        reaction_rate_accumulator._weighted_event_counts,
        reaction_rate_accumulator._event_loss_sample_counts,
        reaction_rate_accumulator._weighted_event_energy_loss_eV,
        reaction_rate_accumulator._weighted_event_counts_by_block,
        reaction_rate_accumulator._weighted_event_energy_loss_eV_by_block,
        bool(zero_high_energy_policy),
        float(max_cross_section_energy_eV),
        float(max_energy_limit_eV),
        rng,
    )
    (
        positions,
        velocities,
        times,
        weights,
        lineages,
        size,
        event_observation_time_s,
        residence_totals,
        secondary_events,
        error_code,
        error_value,
        error_process_index,
    ) = result
    ensemble.adopt_compiled_storage(
        positions,
        velocities,
        times,
        weights,
        lineages,
        size,
    )
    if sampling_enabled:
        reaction_rate_accumulator.add_compiled_event_observation_time(
            event_observation_time_s
        )
        if transport_observer is not None and residence_totals[0] > 0.0:
            transport_observer.record_aggregate_residence(
                weighted_time_s=float(residence_totals[0]),
                weighted_displacement_m=residence_totals[1:4],
                weighted_energy_time_eV_s=float(residence_totals[4]),
                weighted_energy_displacement_eV_m=residence_totals[5:8],
            )
    if error_code == KERNEL_OK:
        return int(secondary_events)
    if error_code == KERNEL_MAJORANT_EXCEEDED:
        raise RuntimeError(
            "internal_monte_carlo null-collision majorant was exceeded in the "
            f"compiled kernel; collision_frequency_s_inv={error_value:.17g}, "
            f"candidate_majorant_s_inv={trial_frequency:.17g}"
        )
    if error_code == KERNEL_CROSS_SECTION_RANGE_EXCEEDED:
        raise RuntimeError(
            "internal_monte_carlo particle exceeded the cross-section energy "
            "range while high_energy_extrapolation=zero; use hold extrapolation "
            f"or extend the cross-section table; energy_eV={error_value:.17g}"
        )
    if error_code == KERNEL_ENERGY_LIMIT_EXCEEDED:
        raise RuntimeError(
            "internal_monte_carlo particle exceeded "
            "physics.energy_grid_policy.max_eV_limit; "
            f"energy_eV={error_value:.17g}"
        )
    if error_code == KERNEL_ERROR_EXTRAPOLATION:
        index = int(error_process_index)
        process = reaction_rate_accumulator.processes[index]
        raise ValueError(
            "Cross-section interpolation requested above the tabulated range "
            f"for {process.species}:{process.process}; energy_eV={error_value:.17g}"
        )
    if error_code == KERNEL_INVALID_RESIDENCE_SAMPLE:
        raise ValueError(
            "compiled weighted-growth residence sample is invalid; "
            f"energy_eV={error_value:.17g}"
        )
    if error_code == KERNEL_FLIGHT_INTEGRATION_FAILED:
        raise RuntimeError(
            "MC B=0 flight reaction-rate integration did not converge; "
            f"maximum_energy_eV={error_value:.17g}"
        )
    raise RuntimeError(f"unknown compiled Monte Carlo kernel status {error_code}")


def _run_weighted_branching_phase(
    *,
    config: SwarmConfig,
    gas_number_density_m3: float,
    ensemble: _population._ParticleEnsemble,
    barriers: int,
    particles: int,
    trial_frequency: float,
    barrier_dt_s: float,
    electric_field_V_m: np.ndarray,
    magnetic_field_T: np.ndarray,
    magnetic_enabled: bool,
    magnetic_B_T: float,
    projected: list[_ProjectedProcess],
    collision_sampler: PreparedCollisionSampler,
    compiled_plan: _CompiledKernelPlan | None,
    angular_name: str,
    angular_model: object,
    rng: np.random.Generator,
    audit: _audits._EnergyAudit | _audits._NullEnergyAudit,
    run_audit: _audits._MonteCarloRunAudit | _audits._NullMonteCarloRunAudit,
    edges: np.ndarray,
    counts: np.ndarray,
    weighted_hist: np.ndarray,
    weighted_square_hist: np.ndarray,
    reaction_rate_accumulator: TrajectoryReactionRateAccumulator,
    transport_observer: SynchronizedWeightedGrowthFluxObserver | None,
    transport_lag_plan: WeightedGrowthTransportLagPlan,
    sampling_enabled: bool,
    zero_high_energy_policy: bool,
    max_cross_section_energy_eV: float,
    max_energy_limit_eV: float,
    tail_strata_edges_eV: np.ndarray | None,
    tail_reaction_importance: np.ndarray | None = None,
    tail_resample_interval_barriers: int | None = None,
) -> tuple[float, int]:
    """Advance a weighted ensemble between synchronized time barriers."""

    barrier_count = int(barriers)
    if barrier_count <= 0:
        return 0.0, 0
    barrier_dt = float(barrier_dt_s)
    if not np.isfinite(barrier_dt) or barrier_dt <= 0.0:
        raise ValueError("weighted-growth barrier cadence must be positive")
    cumulative_log_growth = 0.0
    secondary_events = 0
    density = float(gas_number_density_m3)
    electric_acceleration = (
        -E_CHARGE_C / ELECTRON_MASS_KG * np.asarray(electric_field_V_m, dtype=float)
    )
    for barrier in range(barrier_count):
        barrier_time = (barrier + 1) * barrier_dt
        block_index = min(
            int(
                (barrier + 0.5) * reaction_rate_accumulator.block_count / barrier_count
            ),
            reaction_rate_accumulator.block_count - 1,
        )
        if compiled_plan is not None:
            secondary_events += _advance_compiled_weighted_barrier(
                config=config,
                ensemble=ensemble,
                barrier_time_s=barrier_time,
                trial_frequency=trial_frequency,
                density=density,
                electric_acceleration_m_s2=electric_acceleration,
                collision_sampler=collision_sampler,
                compiled_plan=compiled_plan,
                rng=rng,
                edges=edges,
                counts=counts,
                weighted_hist=weighted_hist,
                weighted_square_hist=weighted_square_hist,
                reaction_rate_accumulator=reaction_rate_accumulator,
                reaction_rate_block_index=block_index,
                transport_observer=(transport_observer if sampling_enabled else None),
                sampling_enabled=sampling_enabled,
                zero_high_energy_policy=zero_high_energy_policy,
                max_cross_section_energy_eV=max_cross_section_energy_eV,
                max_energy_limit_eV=max_energy_limit_eV,
            )
        else:
            particle = 0
            while particle < len(ensemble):
                tolerance = max(1.0e-12 * barrier_time, 1.0e-300)
                while float(ensemble.times[particle]) < barrier_time - tolerance:
                    remaining = barrier_time - float(ensemble.times[particle])
                    sampled_dt = float(rng.exponential(1.0 / trial_frequency))
                    collision_due = sampled_dt < remaining
                    flight_dt = sampled_dt if collision_due else remaining
                    energy_eV = _advance_to_trial_event(
                        ensemble=ensemble,
                        particle_index=particle,
                        trial_dt_s=flight_dt,
                        electric_field_V_m=electric_field_V_m,
                        magnetic_field_T=magnetic_field_T,
                        magnetic_enabled=magnetic_enabled,
                        magnetic_B_T=magnetic_B_T,
                        audit=audit,
                        run_audit=run_audit,
                        edges=edges,
                        counts=counts,
                        weighted_hist=weighted_hist,
                        weighted_square_hist=weighted_square_hist,
                        reaction_rate_accumulator=reaction_rate_accumulator,
                        reaction_rate_block_index=(
                            block_index if sampling_enabled else None
                        ),
                        transport_observer=(
                            transport_observer if sampling_enabled else None
                        ),
                        sampling_enabled=sampling_enabled,
                        zero_high_energy_policy=zero_high_energy_policy,
                        max_cross_section_energy_eV=max_cross_section_energy_eV,
                        max_energy_limit_eV=max_energy_limit_eV,
                    )
                    if collision_due:
                        secondary_created = _sample_weighted_branching_collision(
                            config=config,
                            ensemble=ensemble,
                            particle_index=particle,
                            energy_eV=energy_eV,
                            projected=projected,
                            collision_sampler=collision_sampler,
                            angular_name=angular_name,
                            angular_model=angular_model,
                            rng=rng,
                            audit=audit,
                            run_audit=run_audit,
                            reaction_rate_accumulator=reaction_rate_accumulator,
                            reaction_rate_block_index=(
                                block_index if sampling_enabled else None
                            ),
                            sampling_enabled=sampling_enabled,
                        )
                        if secondary_created:
                            secondary_events += 1
                ensemble.times[particle] = barrier_time
                particle += 1

        if compiled_plan is None:
            if not np.allclose(
                ensemble.times,
                barrier_time,
                rtol=1.0e-12,
                atol=max(1.0e-15 * barrier_dt, 1.0e-300),
            ):
                raise RuntimeError(
                    "weighted-growth ensemble missed a common-time barrier"
                )
        lag_barriers = (barrier + 1) % transport_lag_plan.block_barriers
        if lag_barriers == 0:
            lag_barriers = transport_lag_plan.block_barriers
        observation_due = bool(
            sampling_enabled
            and transport_observer is not None
            and (
                lag_barriers in transport_lag_plan.lag_barriers
                or barrier + 1 == barrier_count
            )
        )
        if observation_due:
            # At the maximum lag this also closes the raw residence-moment
            # block, before population resampling changes the weights.
            transport_observer.record_plane(
                time_s=barrier_time,
                positions_m=ensemble.positions,
                velocities_m_s=ensemble.velocities,
                weights=ensemble.weights,
                lineages=ensemble.lineages,
                lag_barriers=lag_barriers,
            )

        total_weight = float(np.sum(ensemble.weights))
        if not np.isfinite(total_weight) or total_weight <= 0.0:
            raise RuntimeError("weighted-growth population weight is invalid")
        cumulative_log_growth += float(
            np.log(total_weight / max(float(particles), 1.0e-300))
        )
        population_adjustment_due = bool(
            len(ensemble) > particles or total_weight != float(particles)
        )
        if population_adjustment_due:
            before_resample_energy = (
                ensemble.total_weighted_energy_eV() if audit.enabled else 0.0
            )
            resampled = ensemble.systematic_resample(particles, rng)
            ensemble.normalize_total_weight(float(particles))
            if audit.enabled and (
                resampled or not np.isclose(total_weight, float(particles))
            ):
                audit.record_population_resampling(
                    before_resample_energy,
                    ensemble.total_weighted_energy_eV(),
                )
            if sampling_enabled and (
                resampled or not np.isclose(total_weight, float(particles))
            ):
                run_audit.record_branching_resample()
        tail_resample_due = bool(
            tail_strata_edges_eV is not None
            and tail_resample_interval_barriers is not None
            and tail_resample_interval_barriers > 0
            and (barrier + 1) % tail_resample_interval_barriers == 0
        )
        if tail_resample_due:
            before_tail_energy = (
                ensemble.total_weighted_energy_eV() if audit.enabled else 0.0
            )
            tail_plan = ensemble.threshold_weighted_resample(
                particles,
                tail_strata_edges_eV,
                rng,
                reaction_importance=tail_reaction_importance,
            )
            if tail_plan is not None and audit.enabled:
                audit.record_population_resampling(
                    before_tail_energy,
                    ensemble.total_weighted_energy_eV(),
                )
                if sampling_enabled:
                    run_audit.record_tail_weighted_resample(tail_plan)
        if lag_barriers == transport_lag_plan.block_barriers:
            # The gas and field are homogeneous, so a translation of every
            # trajectory origin leaves the transport process unchanged.  A
            # fresh block origin makes the sampled displacement horizon equal
            # to the lineage horizon used by the qualification gate.
            ensemble.positions[:] = 0.0
            ensemble.lineages = np.arange(len(ensemble), dtype=np.int64)
    return cumulative_log_growth, secondary_events
