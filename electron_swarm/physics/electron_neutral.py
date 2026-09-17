"""Shared electron-neutral kinematics used by particle-like solvers."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from electron_swarm.core.constants import (
    AMU_KG,
    BOLTZMANN_J_K,
    ELECTRON_MASS_KG,
    EV_TO_J,
)


def speed_from_energy_m_s(energy_eV: float | np.ndarray) -> float | np.ndarray:
    energy = np.asarray(energy_eV, dtype=float)
    if np.any(~np.isfinite(energy)) or np.any(energy < 0.0):
        raise ValueError("electron energy must be finite and nonnegative")
    speed = np.sqrt(2.0 * EV_TO_J * energy / ELECTRON_MASS_KG)
    return float(speed) if energy.ndim == 0 else speed


def energy_from_speed_eV(speed_m_s: float | np.ndarray) -> float | np.ndarray:
    speed = np.asarray(speed_m_s, dtype=float)
    if np.any(~np.isfinite(speed)) or np.any(speed < 0.0):
        raise ValueError("electron speed must be finite and nonnegative")
    energy = 0.5 * ELECTRON_MASS_KG * speed * speed / EV_TO_J
    return float(energy) if speed.ndim == 0 else energy


def first_order_elastic_recoil_energy_eV(
    incident_energy_eV: float | np.ndarray,
    scattering_cosine: float | np.ndarray,
    target_mass_amu: float,
) -> float | np.ndarray:
    """Stationary-target first-order mass-ratio elastic kinematics."""

    energy = np.asarray(incident_energy_eV, dtype=float)
    cosine = np.asarray(scattering_cosine, dtype=float)
    mass = float(target_mass_amu)
    if (
        np.any(~np.isfinite(energy))
        or np.any(energy < 0.0)
        or np.any(~np.isfinite(cosine))
        or np.any(cosine < -1.0)
        or np.any(cosine > 1.0)
        or not math.isfinite(mass)
        or mass <= 0.0
    ):
        raise ValueError("invalid elastic recoil inputs")
    ratio = ELECTRON_MASS_KG / (mass * AMU_KG)
    result = energy * (1.0 - 2.0 * ratio * (1.0 - cosine))
    result = np.maximum(result, 0.0)
    return float(result) if result.ndim == 0 else result


def maxwell_velocity_component_std_m_s(
    gas_temperature_K: float,
    target_mass_amu: float,
) -> float:
    """Return the one-component thermal-velocity standard deviation."""

    temperature = float(gas_temperature_K)
    mass_amu = float(target_mass_amu)
    if (
        not math.isfinite(temperature)
        or temperature <= 0.0
        or not math.isfinite(mass_amu)
        or mass_amu <= 0.0
    ):
        raise ValueError("thermal target temperature and mass must be positive")
    return math.sqrt(BOLTZMANN_J_K * temperature / (mass_amu * AMU_KG))


def elastic_binary_collision_velocity_m_s(
    electron_velocity_m_s: np.ndarray,
    target_velocity_m_s: np.ndarray,
    scattering_cosine: float,
    azimuth_rad: float,
    target_mass_amu: float,
) -> np.ndarray:
    """Return the exact post-collision electron velocity for an elastic pair.

    The scattering angle rotates the relative velocity in the center-of-mass
    frame.  This conserves the pair momentum and kinetic energy before the
    untracked neutral is discarded.
    """

    electron = np.asarray(electron_velocity_m_s, dtype=float)
    target = np.asarray(target_velocity_m_s, dtype=float)
    cosine = float(scattering_cosine)
    azimuth = float(azimuth_rad)
    mass_amu = float(target_mass_amu)
    if (
        electron.shape != (3,)
        or target.shape != (3,)
        or np.any(~np.isfinite(electron))
        or np.any(~np.isfinite(target))
        or not math.isfinite(cosine)
        or cosine < -1.0
        or cosine > 1.0
        or not math.isfinite(azimuth)
        or not math.isfinite(mass_amu)
        or mass_amu <= 0.0
    ):
        raise ValueError("invalid elastic binary-collision inputs")

    target_mass = mass_amu * AMU_KG
    relative = electron - target
    relative_speed = float(np.linalg.norm(relative))
    center_of_mass = (
        ELECTRON_MASS_KG * electron + target_mass * target
    ) / (ELECTRON_MASS_KG + target_mass)
    if relative_speed <= 0.0:
        return center_of_mass.copy()

    axis = relative / relative_speed
    reference = (
        np.array([1.0, 0.0, 0.0])
        if abs(float(axis[0])) < 0.9
        else np.array([0.0, 1.0, 0.0])
    )
    transverse_1 = np.cross(axis, reference)
    transverse_1 /= max(float(np.linalg.norm(transverse_1)), 1.0e-300)
    transverse_2 = np.cross(axis, transverse_1)
    sine = math.sqrt(max(1.0 - cosine * cosine, 0.0))
    rotated_relative = relative_speed * (
        cosine * axis
        + sine
        * (
            math.cos(azimuth) * transverse_1
            + math.sin(azimuth) * transverse_2
        )
    )
    return center_of_mass + (
        target_mass / (ELECTRON_MASS_KG + target_mass)
    ) * rotated_relative


@dataclass(frozen=True, slots=True)
class IonizationDaughters:
    threshold_loss_eV: float
    primary_energy_eV: float
    secondary_energy_eV: float | None

    @property
    def energies_eV(self) -> tuple[float, ...]:
        if self.secondary_energy_eV is None:
            return (self.primary_energy_eV,)
        return (self.primary_energy_eV, self.secondary_energy_eV)


def ionization_daughter_energy_arrays(
    incident_energy_eV: np.ndarray,
    threshold_eV: float,
    *,
    model: str,
    secondary_electron_energy_eV: float = 0.0,
) -> tuple[np.ndarray, ...]:
    """Vectorized daughter energies for deterministic transfer operators."""

    incident = np.asarray(incident_energy_eV, dtype=float)
    threshold = float(threshold_eV)
    secondary = float(secondary_electron_energy_eV)
    if (
        np.any(~np.isfinite(incident))
        or np.any(incident < 0.0)
        or not math.isfinite(threshold)
        or threshold < 0.0
        or not math.isfinite(secondary)
        or secondary < 0.0
    ):
        raise ValueError("invalid ionization kinematics inputs")
    available = np.maximum(incident - threshold, 0.0)
    if model == "equal":
        half = 0.5 * available
        return half, half.copy()
    if model == "primary_secondary":
        secondary_energy = np.minimum(secondary, available)
        return available - secondary_energy, secondary_energy
    if model == "loss_only":
        return (available,)
    raise ValueError(f"unknown ionization energy-sharing model {model!r}")


def ionization_daughters(
    incident_energy_eV: float,
    threshold_eV: float,
    *,
    model: str,
    secondary_electron_energy_eV: float = 0.0,
) -> IonizationDaughters:
    """Return deterministic daughter energies for the configured source model."""

    threshold = float(threshold_eV)
    energies = ionization_daughter_energy_arrays(
        np.asarray([incident_energy_eV], dtype=float),
        threshold,
        model=model,
        secondary_electron_energy_eV=secondary_electron_energy_eV,
    )
    primary = float(energies[0][0])
    secondary = float(energies[1][0]) if len(energies) == 2 else None
    return IonizationDaughters(threshold, primary, secondary)
