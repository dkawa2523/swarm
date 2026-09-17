"""Mutable state and initialization for one internal Monte Carlo case."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

import electron_swarm.solvers.monte_carlo.audits as _audits
import electron_swarm.solvers.monte_carlo.population as _population
from electron_swarm.core.constants import TOWNSEND
from electron_swarm.solvers.monte_carlo.histogram import _mc_energy_edges
from electron_swarm.solvers.monte_carlo.kinematics import speed_from_energy_m_s
from electron_swarm.solvers.monte_carlo.reaction_rates import (
    ReactionRateEstimate,
    TrajectoryReactionRateAccumulator,
)
from electron_swarm.solvers.monte_carlo.setup import _MonteCarloRunSetup
from electron_swarm.solvers.monte_carlo.direct_transport import (
    SynchronizedFluxTransportObserver,
    direct_transport_observation_times,
)
from electron_swarm.solvers.monte_carlo.weighted_transport import (
    SynchronizedWeightedGrowthFluxObserver,
)


_EnergyAudit = _audits._EnergyAudit | _audits._NullEnergyAudit
_RunAudit = _audits._MonteCarloRunAudit | _audits._NullMonteCarloRunAudit
_TransportObserver = (
    SynchronizedFluxTransportObserver | SynchronizedWeightedGrowthFluxObserver
)


@dataclass(slots=True)
class _MonteCarloCaseState:
    """State mutated in simulation order by the phases of one field case."""

    e_over_n_Td: float
    case_id: str
    case_seed: int
    electric_field_scalar_V_m: float
    electric_field_V_m: np.ndarray
    audit: _EnergyAudit
    run_audit: _RunAudit
    reaction_rate_accumulator: TrajectoryReactionRateAccumulator
    transport_observer: _TransportObserver
    ensemble: _population._ParticleEnsemble
    edges: np.ndarray
    energy: np.ndarray
    widths: np.ndarray
    counts: np.ndarray
    weighted_hist: np.ndarray
    weighted_square_hist: np.ndarray
    cumulative_log_growth: float = 0.0
    audit_cumulative_log_growth: float = 0.0
    production_secondary_events: int = 0
    event_evidence_estimates: list[ReactionRateEstimate] | None = None
    tail_sampling_used: bool = False
    production_mean_time: float | None = None


def _new_energy_audit(
    setup: _MonteCarloRunSetup,
    ensemble: _population._ParticleEnsemble,
) -> _EnergyAudit:
    audit: _EnergyAudit = (
        _audits._EnergyAudit() if setup.collect_audit else _audits._NullEnergyAudit()
    )
    audit.tracked_particle_initial_energy_eV = ensemble.total_weighted_energy_eV()
    return audit


def _new_run_audit(
    setup: _MonteCarloRunSetup,
    *,
    production_collisions: int,
) -> _RunAudit:
    if not setup.collect_audit:
        return _audits._NullMonteCarloRunAudit()
    return _audits._MonteCarloRunAudit(
        seed=setup.cfg.seed,
        particles=setup.particles,
        population_model=setup.population_model,
        warmup_collisions=setup.warmup_collisions,
        production_collisions=production_collisions,
        trial_collision_frequency_s_inv=setup.trial_frequency,
        max_cross_section_energy_eV=setup.max_cross_section_energy,
        tail_threshold_eV=setup.tail_threshold,
    )


def _initialize_monte_carlo_case(
    setup: _MonteCarloRunSetup,
    rng: np.random.Generator,
    *,
    e_over_n_Td: float,
    case_id: str,
    case_seed: int,
    field_polarity: int = 1,
) -> _MonteCarloCaseState:
    if field_polarity not in {-1, 1}:
        raise ValueError("Monte Carlo field polarity must be +1 or -1")
    electric_field_scalar = (
        float(field_polarity) * float(e_over_n_Td) * TOWNSEND * setup.density
    )
    electric_field = np.array([0.0, 0.0, electric_field_scalar], dtype=float)
    audit: _EnergyAudit = (
        _audits._EnergyAudit() if setup.collect_audit else _audits._NullEnergyAudit()
    )
    run_audit = _new_run_audit(
        setup,
        production_collisions=setup.collisions,
    )
    reaction_rate_accumulator = TrajectoryReactionRateAccumulator(
        processes=tuple(item.process for item in setup.projected),
        fractions=tuple(item.fraction for item in setup.projected),
        event_sampled=setup.event_sampled,
        gas_number_density_m3=setup.density,
        angular_scattering_model=setup.angular_name,
    )
    if setup.branching_active:
        transport_observer: _TransportObserver = SynchronizedWeightedGrowthFluxObserver(
            particles=setup.particles,
            electric_field_V_m=electric_field,
            barrier_cadence_s=setup.barrier_dt,
            lag_plan=setup.transport_lag_plan,
        )
    else:
        transport_observer = SynchronizedFluxTransportObserver(
            observation_times_s=direct_transport_observation_times(
                production_trial_steps=setup.collisions,
                trial_frequency_s_inv=setup.trial_frequency,
            ),
            particles=setup.particles,
            electric_field_V_m=electric_field,
        )
    initial_speed = speed_from_energy_m_s(1.0)
    ensemble = _population._ParticleEnsemble.initialize(
        setup.particles,
        initial_speed,
        rng,
    )
    audit.tracked_particle_initial_energy_eV = ensemble.total_weighted_energy_eV()
    edges = _mc_energy_edges(
        setup.max_energy_limit,
        thresholds_eV=setup.reaction_thresholds,
        processes=(item.process for item in setup.projected),
    )
    widths = np.diff(edges)
    energy = 0.5 * (edges[:-1] + edges[1:])
    return _MonteCarloCaseState(
        e_over_n_Td=float(e_over_n_Td),
        case_id=case_id,
        case_seed=int(case_seed),
        electric_field_scalar_V_m=electric_field_scalar,
        electric_field_V_m=electric_field,
        audit=audit,
        run_audit=run_audit,
        reaction_rate_accumulator=reaction_rate_accumulator,
        transport_observer=transport_observer,
        ensemble=ensemble,
        edges=edges,
        energy=energy,
        widths=widths,
        counts=np.zeros_like(energy, dtype=int),
        weighted_hist=np.zeros_like(energy, dtype=float),
        weighted_square_hist=np.zeros_like(energy, dtype=float),
    )
