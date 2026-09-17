"""Orbit integration helpers for the internal Monte Carlo backend."""

from __future__ import annotations

import numpy as np

from electron_swarm.core.constants import E_CHARGE_C, ELECTRON_MASS_KG
from electron_swarm.solvers.monte_carlo.histogram import _record_histogram_flight
from electron_swarm.solvers.monte_carlo.flight_integration import (
    accumulate_dc_flight_histogram,
    integrate_dc_flight_rates,
)
from electron_swarm.solvers.monte_carlo.reaction_rates import (
    TrajectoryReactionRateAccumulator,
)
from electron_swarm.solvers.monte_carlo.kinematics import (
    boris_push,
    energy_from_velocity_eV,
)
from electron_swarm.solvers.monte_carlo.direct_transport import (
    SynchronizedFluxTransportObserver,
)
from electron_swarm.solvers.monte_carlo.weighted_transport import (
    SynchronizedWeightedGrowthFluxObserver,
)


def _magnetic_flight_energy_eV(
    before_velocity_m_s: np.ndarray,
    after_velocity_m_s: np.ndarray,
) -> float:
    """Preserve magnetic no-work behavior in the fallback orbit quadrature."""

    before = np.asarray(before_velocity_m_s, dtype=float)
    after = np.asarray(after_velocity_m_s, dtype=float)
    return 0.5 * (
        energy_from_velocity_eV(before) + energy_from_velocity_eV(after)
    )


def _dc_flight_quadrature(
    before_velocity_m_s: np.ndarray,
    after_velocity_m_s: np.ndarray,
    dt_s: float,
) -> tuple[tuple[float, float, np.ndarray], ...]:
    """Three-point Gauss integration samples for a constant-acceleration flight."""

    before = np.asarray(before_velocity_m_s, dtype=float)
    after = np.asarray(after_velocity_m_s, dtype=float)
    duration = float(dt_s)
    root = np.sqrt(3.0 / 5.0)
    fractions = (0.5 * (1.0 - root), 0.5, 0.5 * (1.0 + root))
    weights = (5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0)
    delta = after - before
    return tuple(
        (
            energy_from_velocity_eV(before + fraction * delta),
            duration * quadrature_weight,
            (before + fraction * delta) * duration * quadrature_weight,
        )
        for fraction, quadrature_weight in zip(fractions, weights, strict=True)
    )


def _record_dc_flight_residence(
    *,
    before_velocity_m_s: np.ndarray,
    after_velocity_m_s: np.ndarray,
    dt_s: float,
    weight: float,
    edges: np.ndarray,
    counts: np.ndarray,
    weighted_hist: np.ndarray,
    weighted_square_hist: np.ndarray,
    reaction_rate_accumulator: TrajectoryReactionRateAccumulator,
    reaction_rate_block_index: int,
    transport_observer: (
        SynchronizedFluxTransportObserver
        | SynchronizedWeightedGrowthFluxObserver
        | None
    ),
    run_audit,
) -> None:
    """Record exact-bin residence and resolved rate integrals for one B=0 flight."""

    if transport_observer is not None:
        # E is quadratic and E*v is cubic under constant acceleration, so this
        # established three-point rule is exact for the transport moments.
        for sample_energy, sample_dt, sample_displacement in _dc_flight_quadrature(
            before_velocity_m_s,
            after_velocity_m_s,
            dt_s,
        ):
            transport_observer.record_residence(
                displacement_m=sample_displacement,
                sample_energy_eV=sample_energy,
                dt_s=sample_dt,
                weight=weight,
            )

    table = reaction_rate_accumulator._prepared_cross_sections
    maximum_energy = max(
        energy_from_velocity_eV(before_velocity_m_s),
        energy_from_velocity_eV(after_velocity_m_s),
    )
    invalid = table.error_extrapolation & (maximum_energy > table.process_max_eV)
    if np.any(invalid):
        process = reaction_rate_accumulator.processes[int(np.flatnonzero(invalid)[0])]
        raise ValueError(
            "Cross-section interpolation requested above the "
            f"tabulated range for {process.species}:{process.process}; "
            "increase the energy grid or set "
            "cross_sections.high_energy_extrapolation."
        )
    (
        rate_integrals,
        loss_integrals,
        quadrature_workspace,
        interval_workspace,
        depth_workspace,
    ) = reaction_rate_accumulator._flight_rate_workspace
    converged = integrate_dc_flight_rates(
        np.asarray(before_velocity_m_s, dtype=float),
        np.asarray(after_velocity_m_s, dtype=float),
        float(dt_s),
        reaction_rate_accumulator._flight_integration_breakpoints_eV,
        table.grid_eV,
        table.values_m2,
        table.slopes_m2_eV,
        table.process_min_eV,
        table.process_max_eV,
        table.right_values_m2,
        reaction_rate_accumulator._constant_energy_loss_eV,
        reaction_rate_accumulator._recoil_mass_ratio,
        rate_integrals,
        loss_integrals,
        quadrature_workspace,
        interval_workspace,
        depth_workspace,
    )
    if not converged:
        raise RuntimeError("MC B=0 flight reaction-rate integration did not converge")
    reaction_rate_accumulator.record_integrated_residence(
        block_index=reaction_rate_block_index,
        weight=weight,
        dt_s=dt_s,
        rate_integrals_m3=rate_integrals,
        energy_loss_integrals_eV_m3=loss_integrals,
    )

    scratch = getattr(reaction_rate_accumulator, "_flight_histogram_scratch", None)
    if scratch is None or scratch.shape != counts.shape:
        scratch = np.zeros_like(weighted_hist, dtype=float)
        reaction_rate_accumulator._flight_histogram_scratch = scratch
        reaction_rate_accumulator._flight_histogram_touched = np.empty(
            counts.size, dtype=np.int64
        )
    touched_bins = reaction_rate_accumulator._flight_histogram_touched
    touched_count = accumulate_dc_flight_histogram(
        np.asarray(before_velocity_m_s, dtype=float),
        np.asarray(after_velocity_m_s, dtype=float),
        float(dt_s),
        float(weight),
        edges,
        counts,
        weighted_hist,
        weighted_square_hist,
        scratch,
        touched_bins,
    )
    # The audit records exact extrema and one representative per occupied bin;
    # it does not participate in the estimator.
    run_audit.record_energy_sample(energy_from_velocity_eV(before_velocity_m_s))
    run_audit.record_energy_sample(energy_from_velocity_eV(after_velocity_m_s))
    for touched_index in range(touched_count):
        bin_index = int(touched_bins[touched_index])
        run_audit.record_histogram_sample(
            0.5 * (float(edges[bin_index]) + float(edges[bin_index + 1]))
        )


def _advance_to_trial_event(
    *,
    ensemble,
    particle_index: int,
    trial_dt_s: float,
    electric_field_V_m: np.ndarray,
    magnetic_field_T: np.ndarray,
    magnetic_enabled: bool,
    magnetic_B_T: float,
    audit,
    run_audit,
    edges: np.ndarray,
    counts: np.ndarray,
    weighted_hist: np.ndarray,
    weighted_square_hist: np.ndarray,
    reaction_rate_accumulator: TrajectoryReactionRateAccumulator,
    reaction_rate_block_index: int | None,
    transport_observer: (
        SynchronizedFluxTransportObserver
        | SynchronizedWeightedGrowthFluxObserver
        | None
    ),
    sampling_enabled: bool,
    zero_high_energy_policy: bool,
    max_cross_section_energy_eV: float,
    max_energy_limit_eV: float,
) -> float:
    remaining = float(trial_dt_s)
    if remaining <= 0.0 or not np.isfinite(remaining):
        raise ValueError(
            "internal_monte_carlo trial event time must be finite and positive"
        )
    max_step = remaining
    if magnetic_enabled:
        gyro = E_CHARGE_C * float(magnetic_B_T) / ELECTRON_MASS_KG
        max_step = min(max_step, 0.2 / max(gyro, 1.0))
    weight = float(ensemble.weights[particle_index])
    final_energy = max(
        energy_from_velocity_eV(ensemble.velocities[particle_index]), 1.0e-6
    )
    while remaining > 0.0:
        dt = min(remaining, max_step)
        if sampling_enabled and transport_observer is not None:
            dt = transport_observer.cap_step(
                particle_index,
                float(ensemble.times[particle_index]),
                dt,
            )
            if dt <= 0.0:
                if not transport_observer.record_due(
                    particle_index=particle_index,
                    time_s=float(ensemble.times[particle_index]),
                    position_m=ensemble.positions[particle_index],
                    velocity_m_s=ensemble.velocities[particle_index],
                    weight=weight,
                ):
                    raise RuntimeError("MC transport observation time did not advance")
                continue
        before_push_velocity = ensemble.velocities[particle_index].copy()
        before_push_energy = energy_from_velocity_eV(before_push_velocity)
        ensemble.velocities[particle_index] = boris_push(
            ensemble.velocities[particle_index],
            electric_field_V_m,
            magnetic_field_T,
            dt,
        )
        after_push_energy = energy_from_velocity_eV(
            ensemble.velocities[particle_index]
        )
        audit.record_field_push(before_push_energy, after_push_energy, weight)
        if magnetic_enabled:
            # Keep the established Boris-orbit position rule for magnetic
            # trajectories; full Lorentz-orbit accuracy is a separate scope.
            position_velocity = ensemble.velocities[particle_index]
        else:
            # With B=0, acceleration is constant over a free flight.  The
            # trapezoidal velocity is therefore the exact displacement and
            # avoids a systematic field-direction bias in diffusion moments.
            position_velocity = 0.5 * (
                before_push_velocity + ensemble.velocities[particle_index]
            )
        ensemble.positions[particle_index] += position_velocity * dt
        ensemble.times[particle_index] += dt
        if sampling_enabled and transport_observer is not None:
            transport_observer.record_due(
                particle_index=particle_index,
                time_s=float(ensemble.times[particle_index]),
                position_m=ensemble.positions[particle_index],
                velocity_m_s=ensemble.velocities[particle_index],
                weight=weight,
            )
        final_energy = max(after_push_energy, 1.0e-6)
        if zero_high_energy_policy and final_energy > max_cross_section_energy_eV:
            raise RuntimeError(
                "internal_monte_carlo particle exceeded the cross-section "
                "energy range while high_energy_extrapolation=zero; use "
                "hold extrapolation or extend the cross-section table"
            )
        if final_energy > max_energy_limit_eV:
            raise RuntimeError(
                "internal_monte_carlo particle exceeded "
                "physics.energy_grid_policy.max_eV_limit"
            )
        run_audit.record_orbit_substep()
        if sampling_enabled:
            if reaction_rate_block_index is None:
                raise RuntimeError("production reaction-rate sample has no batch index")
            if magnetic_enabled:
                midpoint_energy = _magnetic_flight_energy_eV(
                    before_push_velocity,
                    ensemble.velocities[particle_index],
                )
                quadrature = ((midpoint_energy, dt, position_velocity * dt),)
                histogram_energy: list[float] = []
                histogram_contribution: list[float] = []
                for sample_energy, sample_dt, sample_displacement in quadrature:
                    if transport_observer is not None:
                        transport_observer.record_residence(
                            displacement_m=sample_displacement,
                            sample_energy_eV=sample_energy,
                            dt_s=sample_dt,
                            weight=weight,
                        )
                    run_audit.record_energy_sample(sample_energy)
                    reaction_rate_accumulator.record_residence(
                        block_index=reaction_rate_block_index,
                        energy_eV=sample_energy,
                        weight=weight,
                        dt_s=sample_dt,
                    )
                    histogram_energy.append(sample_energy)
                    histogram_contribution.append(weight * sample_dt)
                _record_histogram_flight(
                    edges=edges,
                    counts=counts,
                    weighted_hist=weighted_hist,
                    weighted_square_hist=weighted_square_hist,
                    run_audit=run_audit,
                    energy_eV=np.asarray(histogram_energy, dtype=float),
                    contributions=np.asarray(histogram_contribution, dtype=float),
                )
            else:
                _record_dc_flight_residence(
                    before_velocity_m_s=before_push_velocity,
                    after_velocity_m_s=ensemble.velocities[particle_index],
                    dt_s=dt,
                    weight=weight,
                    edges=edges,
                    counts=counts,
                    weighted_hist=weighted_hist,
                    weighted_square_hist=weighted_square_hist,
                    reaction_rate_accumulator=reaction_rate_accumulator,
                    reaction_rate_block_index=reaction_rate_block_index,
                    transport_observer=transport_observer,
                    run_audit=run_audit,
                )
        remaining -= dt
        if remaining <= max(1.0e-15 * trial_dt_s, 1.0e-300):
            break
    return final_energy
