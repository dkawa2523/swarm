"""State-transition phases for one internal Monte Carlo case."""

from __future__ import annotations

import numpy as np

import electron_swarm.solvers.monte_carlo.population as _population
import electron_swarm.solvers.monte_carlo.weighted_runtime as _weighted_runtime
from electron_swarm.core.cross_sections import ProcessType
from electron_swarm.solvers.monte_carlo.case_state import (
    _MonteCarloCaseState,
    _new_energy_audit,
    _new_run_audit,
)
from electron_swarm.solvers.monte_carlo.collisions import (
    SCATTERING_TYPES,
    TRAJECTORY_REACTION_TYPES,
    _PostReactionOutcome,
    _ionization_daughters,
    _moment_cross_sections,
    _post_reaction_outcome,
    _thermal_elastic_velocity_after_collision,
)
from electron_swarm.solvers.monte_carlo.kinematics import (
    energy_from_velocity_eV,
    random_direction,
    speed_from_energy_m_s,
)
from electron_swarm.solvers.monte_carlo.orbit import _advance_to_trial_event
from electron_swarm.solvers.monte_carlo.reaction_rates import (
    TrajectoryReactionRateAccumulator,
)
from electron_swarm.solvers.monte_carlo.setup import _MonteCarloRunSetup
from electron_swarm.solvers.monte_carlo.weighted_transport import (
    SynchronizedWeightedGrowthFluxObserver,
)


def _run_fixed_particle_phase(
    setup: _MonteCarloRunSetup,
    state: _MonteCarloCaseState,
    rng: np.random.Generator,
) -> None:
    """Advance the direct fixed-population warmup and production window."""

    total_collision_steps = (
        0 if setup.branching_active else setup.warmup_collisions + setup.collisions
    )
    for step in range(total_collision_steps):
        sampling_enabled = step >= setup.warmup_collisions
        if step == setup.warmup_collisions and setup.warmup_collisions:
            state.ensemble.positions[:] = 0.0
            state.ensemble.times[:] = 0.0
            if setup.branching_active:
                state.ensemble.normalize_total_weight(float(setup.particles))
            # Keep the detailed energy audit on the same production
            # window as the direct residence and event estimators.
            state.audit = _new_energy_audit(setup, state.ensemble)
        particles_this_round = len(state.ensemble)
        for particle_index in range(particles_this_round):
            if particle_index >= len(state.ensemble):
                break
            trial_dt = float(rng.exponential(1.0 / setup.trial_frequency))
            weight = float(state.ensemble.weights[particle_index])
            reaction_rate_block_index = (
                state.reaction_rate_accumulator.block_index_for_trial(
                    production_step=step - setup.warmup_collisions,
                    production_steps=setup.collisions,
                    particle_index=particle_index,
                    particles_this_round=particles_this_round,
                )
                if sampling_enabled
                else None
            )
            energy_eV = _advance_to_trial_event(
                ensemble=state.ensemble,
                particle_index=particle_index,
                trial_dt_s=trial_dt,
                electric_field_V_m=state.electric_field_V_m,
                magnetic_field_T=setup.B,
                magnetic_enabled=setup.magnetic_enabled,
                magnetic_B_T=setup.effective_B_T,
                audit=state.audit,
                run_audit=state.run_audit,
                edges=state.edges,
                counts=state.counts,
                weighted_hist=state.weighted_hist,
                weighted_square_hist=state.weighted_square_hist,
                reaction_rate_accumulator=state.reaction_rate_accumulator,
                reaction_rate_block_index=reaction_rate_block_index,
                transport_observer=state.transport_observer,
                sampling_enabled=sampling_enabled,
                zero_high_energy_policy=setup.zero_high_energy_policy,
                max_cross_section_energy_eV=setup.max_cross_section_energy,
                max_energy_limit_eV=setup.max_energy_limit,
            )
            trial = setup.collision_sampler.sample_trial(
                state.ensemble.velocities[particle_index],
                rng,
            )
            state.audit.record_trial_ratio(trial.acceptance_probability)
            if trial.process_index is None:
                if sampling_enabled:
                    state.run_audit.record_null_collision()
                continue
            choice = int(trial.process_index)
            proc = setup.projected[choice].process
            if sampling_enabled:
                state.run_audit.record_collision(proc.process_type)
            event_energy_loss_eV: float | None = None
            if proc.process_type in SCATTERING_TYPES:
                if (
                    trial.target_velocity_m_s is None
                    or trial.relative_energy_eV is None
                ):
                    raise RuntimeError("thermal elastic collision lacks target state")
                relative_energy_eV = float(trial.relative_energy_eV)
                if setup.angular_name == "maxent_p1":
                    total, momentum = _moment_cross_sections(
                        setup.projected,
                        proc.species,
                        relative_energy_eV,
                    )
                    mu = float(
                        setup.angular_model.sample_mu(
                            relative_energy_eV,
                            rng,
                            sigma_total=total,
                            sigma_momentum=momentum,
                        )
                    )
                else:
                    mu = float(setup.angular_model.sample_mu(relative_energy_eV, rng))
                post_velocity = _thermal_elastic_velocity_after_collision(
                    setup.config,
                    proc,
                    state.ensemble.velocities[particle_index],
                    trial.target_velocity_m_s,
                    mu,
                    rng,
                )
                post_energy = float(energy_from_velocity_eV(post_velocity))
                if post_energy > setup.max_energy_limit * (1.0 + 1.0e-12):
                    raise RuntimeError(
                        "internal_monte_carlo thermal collision exceeded "
                        "physics.energy_grid_policy.max_eV_limit"
                    )
                event_energy_loss_eV = float(energy_eV - post_energy)
                state.audit.record_elastic_collision(energy_eV, post_energy, weight)
                state.ensemble.velocities[particle_index] = post_velocity
            elif proc.process_type in TRAJECTORY_REACTION_TYPES:
                if (
                    proc.process_type == ProcessType.IONIZATION
                    and setup.branching_active
                ):
                    daughters = _ionization_daughters(
                        setup.config,
                        proc.threshold_eV,
                        energy_eV,
                    )
                    outcome = _PostReactionOutcome(
                        tracked_energy_eV=float(sum(daughters.energies_eV)),
                        ionization_threshold_loss_eV=daughters.threshold_loss_eV,
                    )
                    event_energy_loss_eV = float(daughters.threshold_loss_eV)
                    state.audit.record_reaction(outcome, weight)
                    energies = daughters.energies_eV
                    state.ensemble.velocities[particle_index] = random_direction(
                        rng
                    ) * speed_from_energy_m_s(max(energies[0], 1.0e-4))
                    if len(energies) > 1:
                        state.ensemble.append_particle(
                            position=state.ensemble.positions[particle_index].copy(),
                            velocity=random_direction(rng)
                            * speed_from_energy_m_s(max(energies[1], 1.0e-4)),
                            time_s=float(state.ensemble.times[particle_index]),
                            weight=weight,
                            lineage=int(state.ensemble.lineages[particle_index]),
                        )
                        if sampling_enabled:
                            state.run_audit.record_secondary_electron()
                else:
                    outcome = _post_reaction_outcome(
                        setup.config,
                        proc.process_type,
                        proc.threshold_eV,
                        energy_eV,
                        rng,
                    )
                    event_energy_loss_eV = float(
                        outcome.inelastic_energy_loss_eV
                        + outcome.ionization_threshold_loss_eV
                    )
                    state.audit.record_reaction(outcome, weight)
                    state.ensemble.velocities[particle_index] = random_direction(
                        rng
                    ) * speed_from_energy_m_s(outcome.tracked_energy_eV)
            if sampling_enabled:
                if event_energy_loss_eV is None:
                    raise RuntimeError(
                        "sampled MC collision has no physical energy-loss audit"
                    )
                state.reaction_rate_accumulator.record_event(
                    choice,
                    weight,
                    energy_loss_eV=event_energy_loss_eV,
                    block_index=reaction_rate_block_index,
                )
        if setup.branching_active and len(state.ensemble) > 2 * setup.particles:
            before_resample_energy = state.ensemble.total_weighted_energy_eV()
            if state.ensemble.systematic_resample(setup.particles, rng):
                after_resample_energy = state.ensemble.total_weighted_energy_eV()
                state.audit.record_population_resampling(
                    before_resample_energy,
                    after_resample_energy,
                )
                if sampling_enabled:
                    state.run_audit.record_branching_resample()


def _advance_weighted_phase(
    setup: _MonteCarloRunSetup,
    state: _MonteCarloCaseState,
    rng: np.random.Generator,
    *,
    ensemble: _population._ParticleEnsemble,
    barriers: int,
    reaction_rate_accumulator: TrajectoryReactionRateAccumulator,
    transport_observer: SynchronizedWeightedGrowthFluxObserver | None,
    sampling_enabled: bool,
    tail_strata_edges_eV: np.ndarray | None,
    tail_reaction_importance: np.ndarray | None = None,
    tail_resample_interval_barriers: int | None = None,
) -> tuple[float, int]:
    return _weighted_runtime._run_weighted_branching_phase(
        config=setup.config,
        gas_number_density_m3=setup.density,
        ensemble=ensemble,
        barriers=barriers,
        particles=setup.particles,
        trial_frequency=setup.trial_frequency,
        barrier_dt_s=setup.barrier_dt,
        electric_field_V_m=state.electric_field_V_m,
        magnetic_field_T=setup.B,
        magnetic_enabled=setup.magnetic_enabled,
        magnetic_B_T=setup.effective_B_T,
        projected=setup.projected,
        collision_sampler=setup.collision_sampler,
        compiled_plan=setup.compiled_plan,
        angular_name=setup.angular_name,
        angular_model=setup.angular_model,
        rng=rng,
        audit=state.audit,
        run_audit=state.run_audit,
        edges=state.edges,
        counts=state.counts,
        weighted_hist=state.weighted_hist,
        weighted_square_hist=state.weighted_square_hist,
        reaction_rate_accumulator=reaction_rate_accumulator,
        transport_observer=transport_observer,
        transport_lag_plan=setup.transport_lag_plan,
        sampling_enabled=sampling_enabled,
        zero_high_energy_policy=setup.zero_high_energy_policy,
        max_cross_section_energy_eV=setup.max_cross_section_energy,
        max_energy_limit_eV=setup.max_energy_limit,
        tail_strata_edges_eV=tail_strata_edges_eV,
        tail_reaction_importance=tail_reaction_importance,
        tail_resample_interval_barriers=tail_resample_interval_barriers,
    )


def _run_tail_phase(
    setup: _MonteCarloRunSetup,
    state: _MonteCarloCaseState,
    rng: np.random.Generator,
    *,
    tail_ensemble: _population._ParticleEnsemble,
) -> None:
    """Run the independent rare-tail ensemble and merge its estimators."""

    main_reaction_rate_accumulator = state.reaction_rate_accumulator
    main_counts = state.counts.copy()
    main_weighted_hist = state.weighted_hist.copy()
    main_weighted_square_hist = state.weighted_square_hist.copy()
    state.ensemble = tail_ensemble
    state.counts = np.zeros_like(state.energy, dtype=int)
    state.weighted_hist = np.zeros_like(state.energy, dtype=float)
    state.weighted_square_hist = np.zeros_like(state.energy, dtype=float)
    state.reaction_rate_accumulator = TrajectoryReactionRateAccumulator(
        processes=tuple(item.process for item in setup.projected),
        fractions=tuple(item.fraction for item in setup.projected),
        event_sampled=setup.event_sampled,
        gas_number_density_m3=setup.density,
        angular_scattering_model=setup.angular_name,
        poisson_event_evidence=False,
    )
    state.audit = _new_energy_audit(setup, state.ensemble)
    state.run_audit = _new_run_audit(
        setup,
        production_collisions=setup.tail_collisions,
    )
    state.audit_cumulative_log_growth, _ = _advance_weighted_phase(
        setup,
        state,
        rng,
        ensemble=state.ensemble,
        barriers=setup.tail_collisions,
        reaction_rate_accumulator=state.reaction_rate_accumulator,
        transport_observer=None,
        sampling_enabled=True,
        tail_strata_edges_eV=setup.tail_strata_edges,
        tail_reaction_importance=setup.tail_reaction_importance,
        tail_resample_interval_barriers=setup.tail_resample_interval,
    )
    state.counts += main_counts
    state.weighted_hist += main_weighted_hist
    state.weighted_square_hist += main_weighted_square_hist
    main_reaction_rate_accumulator.merge_residence_statistics(
        state.reaction_rate_accumulator
    )
    state.reaction_rate_accumulator = main_reaction_rate_accumulator


def _run_weighted_transport_phase(
    setup: _MonteCarloRunSetup,
    state: _MonteCarloCaseState,
    rng: np.random.Generator,
) -> _population._ParticleEnsemble:
    """Run weighted warmup and production without a secondary tail phase."""

    assert isinstance(
        state.transport_observer,
        SynchronizedWeightedGrowthFluxObserver,
    )
    if setup.warmup_collisions:
        _advance_weighted_phase(
            setup,
            state,
            rng,
            ensemble=state.ensemble,
            barriers=setup.warmup_collisions,
            reaction_rate_accumulator=state.reaction_rate_accumulator,
            transport_observer=None,
            sampling_enabled=False,
            tail_strata_edges_eV=None,
        )
        state.ensemble.positions[:] = 0.0
        state.ensemble.times[:] = 0.0
        state.ensemble.normalize_total_weight(float(setup.particles))
        state.ensemble.lineages = np.arange(len(state.ensemble), dtype=np.int64)
        state.audit = _new_energy_audit(setup, state.ensemble)
    tail_ensemble = state.ensemble.copy()
    (
        state.cumulative_log_growth,
        state.production_secondary_events,
    ) = _advance_weighted_phase(
        setup,
        state,
        rng,
        ensemble=state.ensemble,
        barriers=setup.collisions,
        reaction_rate_accumulator=state.reaction_rate_accumulator,
        transport_observer=state.transport_observer,
        sampling_enabled=True,
        tail_strata_edges_eV=None,
    )
    state.production_mean_time = float(
        np.sum(state.ensemble.weights * state.ensemble.times)
        / max(float(np.sum(state.ensemble.weights)), 1.0e-300)
    )
    state.audit_cumulative_log_growth = state.cumulative_log_growth
    return tail_ensemble


def _run_weighted_phase(
    setup: _MonteCarloRunSetup,
    state: _MonteCarloCaseState,
    rng: np.random.Generator,
) -> None:
    """Run weighted transport and activate independent rare-tail sampling."""

    tail_ensemble = _run_weighted_transport_phase(setup, state, rng)
    main_reaction_rate_estimates = state.reaction_rate_accumulator.estimates()
    state.tail_sampling_used = bool(
        setup.tail_strata_edges is not None
        and any(
            estimate.rate_coefficient_m3_s <= 0.0
            or estimate.relative_standard_error is None
            or estimate.relative_standard_error > setup.cfg.tail_rate_rse_trigger
            for item, estimate in zip(
                setup.projected,
                main_reaction_rate_estimates,
                strict=True,
            )
            if item.process.process_type in TRAJECTORY_REACTION_TYPES
        )
    )
    if state.tail_sampling_used:
        # Transport and rare-tail sampling have different statistical jobs.
        # The copied ensemble owns tail refinement; the production observer is
        # never resampled, and histogram/rate residence evidence is pooled.
        state.event_evidence_estimates = main_reaction_rate_estimates
        _run_tail_phase(
            setup,
            state,
            rng,
            tail_ensemble=tail_ensemble,
        )
    else:
        state.audit_cumulative_log_growth = state.cumulative_log_growth
