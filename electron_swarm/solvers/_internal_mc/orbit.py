"""Orbit integration helpers for the internal Monte Carlo backend."""

from __future__ import annotations

import numpy as np

from electron_swarm.core.constants import E_CHARGE_C, ELECTRON_MASS_KG, EV_TO_J
from electron_swarm.solvers._internal_mc.histogram import _record_histogram_residence


def magnetic_field_vector(B_T: float, angle_EB_deg: float) -> np.ndarray:
    if B_T < 0.0 or not np.isfinite(B_T) or not np.isfinite(angle_EB_deg):
        raise ValueError("magnetic field magnitude and angle must be finite")
    if angle_EB_deg < 0.0 or angle_EB_deg > 180.0:
        raise ValueError("magnetic field angle_EB_deg must be in [0, 180]")
    angle = np.deg2rad(float(angle_EB_deg))
    return float(B_T) * np.array([np.sin(angle), 0.0, np.cos(angle)], dtype=float)

def boris_push(
    velocity_m_s: np.ndarray,
    electric_field_V_m: np.ndarray,
    magnetic_field_T: np.ndarray,
    dt_s: float,
) -> np.ndarray:
    """Advance electron velocity with the Boris algorithm."""

    v = np.asarray(velocity_m_s, dtype=float)
    E = np.asarray(electric_field_V_m, dtype=float)
    B = np.asarray(magnetic_field_T, dtype=float)
    qmdt2 = -E_CHARGE_C / ELECTRON_MASS_KG * float(dt_s) * 0.5
    v_minus = v + qmdt2 * E
    t = qmdt2 * B
    s = 2.0 * t / (1.0 + float(np.dot(t, t)))
    v_prime = v_minus + np.cross(v_minus, t)
    v_plus = v_minus + np.cross(v_prime, s)
    return v_plus + qmdt2 * E

def _energy_from_velocity(velocity: np.ndarray) -> float:
    return float(0.5 * ELECTRON_MASS_KG * np.dot(velocity, velocity) / EV_TO_J)

def _speed_from_energy(energy_eV: float) -> float:
    return float(np.sqrt(max(2.0 * energy_eV * EV_TO_J / ELECTRON_MASS_KG, 0.0)))

def _random_direction(rng: np.random.Generator) -> np.ndarray:
    mu = rng.uniform(-1.0, 1.0)
    phi = rng.uniform(0.0, 2.0 * np.pi)
    sint = np.sqrt(max(1.0 - mu * mu, 0.0))
    return np.array([sint * np.cos(phi), sint * np.sin(phi), mu], dtype=float)

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
    sampling_enabled: bool,
    zero_high_energy_policy: bool,
    max_cross_section_energy_eV: float,
    max_energy_limit_eV: float,
) -> float:
    remaining = float(trial_dt_s)
    if remaining <= 0.0 or not np.isfinite(remaining):
        raise ValueError("internal_monte_carlo trial event time must be finite and positive")
    max_step = remaining
    if magnetic_enabled:
        gyro = E_CHARGE_C * float(magnetic_B_T) / ELECTRON_MASS_KG
        max_step = min(max_step, 0.2 / max(gyro, 1.0))
    weight = float(ensemble.weights[particle_index])
    final_energy = max(_energy_from_velocity(ensemble.velocities[particle_index]), 1.0e-6)
    while remaining > 0.0:
        dt = min(remaining, max_step)
        before_push_energy = _energy_from_velocity(ensemble.velocities[particle_index])
        ensemble.velocities[particle_index] = boris_push(
            ensemble.velocities[particle_index],
            electric_field_V_m,
            magnetic_field_T,
            dt,
        )
        after_push_energy = _energy_from_velocity(ensemble.velocities[particle_index])
        audit.record_field_push(before_push_energy, after_push_energy, weight)
        ensemble.positions[particle_index] += ensemble.velocities[particle_index] * dt
        ensemble.times[particle_index] += dt
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
            midpoint_energy = 0.5 * (before_push_energy + after_push_energy)
            run_audit.record_energy_sample(midpoint_energy)
            _record_histogram_residence(
                edges=edges,
                counts=counts,
                weighted_hist=weighted_hist,
                weighted_square_hist=weighted_square_hist,
                run_audit=run_audit,
                energy_eV=midpoint_energy,
                contribution=weight * dt,
            )
        remaining -= dt
        if remaining <= max(1.0e-15 * trial_dt_s, 1.0e-300):
            break
    return final_energy
