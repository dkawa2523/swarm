"""Electron-electron relaxation target helpers for product schema v2."""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq


_EPS = 1.0e-300


def prepare_relaxation_input(
    energy_eV: np.ndarray,
    widths_eV: np.ndarray,
    eedf_eV_inv: np.ndarray,
) -> np.ndarray:
    """Remove bounded PN roundoff while preserving probability and energy.

    The native PN solution remains untouched. Only a relaxation input with
    negative probability and energy fractions below 1e-8 may be projected;
    its total probability change must remain below 1e-7.
    """
    energy = np.asarray(energy_eV, dtype=float)
    widths = np.asarray(widths_eV, dtype=float)
    values = np.asarray(eedf_eV_inv, dtype=float)
    total = float(np.dot(values, widths))
    if not np.isfinite(total) or total <= 0.0:
        raise ValueError(
            "electron_electron relaxation requires nonzero finite probability"
        )
    values = values / total
    if np.all(values >= 0.0):
        return values
    mass = values * widths
    negative = np.minimum(mass, 0.0)
    target = float(np.dot(energy, mass))
    if (
        target <= 0.0
        or -np.sum(negative) > 1e-8
        or -np.dot(energy, negative) > 1e-8 * target
    ):
        raise ValueError(
            "electron_electron negative input exceeds numerical roundoff tolerance"
        )
    projected = np.maximum(mass, 0.0)
    projected /= np.sum(projected)
    mean = float(np.dot(energy, projected))
    difference = target - mean
    if difference != 0.0:
        support = np.flatnonzero(projected > 0.0)
        index = int(support[-1] if difference > 0.0 else support[0])
        denominator = energy[index] - mean
        if denominator == 0.0 or not 0.0 <= difference / denominator <= 1.0:
            raise ValueError(
                "electron_electron cannot preserve input energy under positivity"
            )
        fraction = difference / denominator
        projected *= 1.0 - fraction
        projected[index] += fraction
    if (
        np.sum(np.abs(projected - mass)) > 1e-7
        or abs(np.dot(energy, projected) - target) > 2e-12 * target
    ):
        raise ValueError(
            "electron_electron input positivity correction exceeds roundoff tolerance"
        )
    return projected / widths


def _normalize(values: np.ndarray, widths: np.ndarray) -> np.ndarray:
    total = float(np.sum(values * widths))
    if total <= 0.0 or not np.isfinite(total):
        raise ValueError("Cannot normalize non-positive EEDF target")
    return values / total


def maxwellian_at_temperature(
    energy_eV, widths_eV, temperature_eV: float
) -> np.ndarray:
    """Finite-domain Maxwell shape at the specified kinetic temperature."""
    if not np.isfinite(temperature_eV) or temperature_eV <= 0:
        raise ValueError("electron_electron temperature must be finite and positive")
    energy = np.asarray(energy_eV, dtype=float)
    shape = np.sqrt(energy) * np.exp(-(energy - energy.min()) / temperature_eV)
    return _normalize(shape, np.asarray(widths_eV, dtype=float))


def mean_energy_eV(
    energy_eV: np.ndarray, widths_eV: np.ndarray, eedf_eV_inv: np.ndarray
) -> float:
    energy = np.asarray(energy_eV, dtype=float)
    widths = np.asarray(widths_eV, dtype=float)
    eedf = np.asarray(eedf_eV_inv, dtype=float)
    norm = max(float(np.sum(eedf * widths)), _EPS)
    return float(np.sum(energy * eedf * widths) / norm)


def maxwellian_eedf_like(
    energy_eV: np.ndarray,
    widths_eV: np.ndarray,
    mean_energy_target_eV: float,
) -> np.ndarray:
    """Return a normalized Maxwellian-like energy PDF matching mean energy."""

    energy = np.asarray(energy_eV, dtype=float)
    widths = np.asarray(widths_eV, dtype=float)
    target_mean = float(mean_energy_target_eV)
    limiting = _normalize(np.sqrt(energy), widths)
    maximum_mean = float(np.dot(energy * widths, limiting))
    if not np.isfinite(target_mean) or not energy.min() <= target_mean <= maximum_mean:
        raise ValueError(
            "mean energy cannot be conserved by a positive-temperature Maxwell target on this energy domain"
        )
    tolerance = 2e-12 * max(target_mean, float(energy.min()))
    if abs(maximum_mean - target_mean) <= tolerance:
        return limiting
    if abs(target_mean - energy.min()) <= tolerance:
        result = np.zeros_like(energy)
        result[np.argmin(energy)] = 1.0 / widths[np.argmin(energy)]
        return result

    def distribution(beta: float) -> np.ndarray:
        return _normalize(
            np.sqrt(energy) * np.exp(-beta * (energy - energy.min())), widths
        )

    def balance(beta: float) -> float:
        return float(np.dot(energy * widths, distribution(beta))) - target_mean

    high = 1.0 / target_mean
    for _ in range(60):
        if balance(high) <= 0:
            break
        high *= 2.0
    else:
        raise ValueError("cannot bracket the energy-conserving Maxwell target")
    beta = brentq(balance, 0.0, high, xtol=1e-14, rtol=1e-13)
    result = distribution(beta)
    if abs(np.dot(energy * widths, result) - target_mean) > tolerance:
        raise FloatingPointError(
            "Maxwell target did not conserve the specified mean energy"
        )
    return result


def relaxation_target(
    energy_eV: np.ndarray,
    widths_eV: np.ndarray,
    eedf_eV_inv: np.ndarray,
    *,
    conserve_mean_energy: bool = True,
    fallback_temperature_eV: float = 2.0,
) -> np.ndarray:
    """Return the EEDF target for the relaxation_postprocess e-e treatment."""

    if conserve_mean_energy:
        target_mean = mean_energy_eV(energy_eV, widths_eV, eedf_eV_inv)
        return maxwellian_eedf_like(energy_eV, widths_eV, target_mean)
    return maxwellian_at_temperature(energy_eV, widths_eV, fallback_temperature_eV)
