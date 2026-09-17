"""Low-level particle kinematics shared across Monte Carlo phases."""

from __future__ import annotations

import numpy as np

from electron_swarm.core.constants import E_CHARGE_C, ELECTRON_MASS_KG, EV_TO_J


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
    if not np.any(B):
        return v_minus + qmdt2 * E
    t = qmdt2 * B
    s = 2.0 * t / (1.0 + float(np.dot(t, t)))
    v_prime = v_minus + np.cross(v_minus, t)
    v_plus = v_minus + np.cross(v_prime, s)
    return v_plus + qmdt2 * E


def energy_from_velocity_eV(velocity_m_s: np.ndarray) -> float:
    return float(
        0.5
        * ELECTRON_MASS_KG
        * np.dot(velocity_m_s, velocity_m_s)
        / EV_TO_J
    )


def speed_from_energy_m_s(energy_eV: float) -> float:
    return float(
        np.sqrt(max(2.0 * energy_eV * EV_TO_J / ELECTRON_MASS_KG, 0.0))
    )


def random_direction(rng: np.random.Generator) -> np.ndarray:
    mu = rng.uniform(-1.0, 1.0)
    phi = rng.uniform(0.0, 2.0 * np.pi)
    sint = np.sqrt(max(1.0 - mu * mu, 0.0))
    return np.array([sint * np.cos(phi), sint * np.sin(phi), mu], dtype=float)
