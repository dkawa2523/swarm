"""Bounded, number- and energy-conserving isotropic FP relaxation.

The positive SG generator relaxes toward a Maxwell shape. When energy is to
be conserved, its inverse temperature is determined by the implicit step's
energy balance, rather than repairing the distribution after the step.
This is a relaxation model, not a density-dependent Landau collision model.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import solve_banded
from scipy.optimize import brentq

from .electron_electron import prepare_relaxation_input


def _bernoulli(values: np.ndarray) -> np.ndarray:
    result = np.empty_like(values)
    small = np.abs(values) < 1e-6
    positive = values > 50.0
    negative = values < -50.0
    regular = ~(small | positive | negative)
    result[small] = 1.0 - values[small] / 2.0 + values[small] ** 2 / 12.0
    result[positive] = values[positive] * np.exp(-values[positive])
    result[negative] = -values[negative]
    result[regular] = values[regular] / np.expm1(values[regular])
    return result


def apply_fp_energy_operator(
    energy_eV: np.ndarray,
    widths_eV: np.ndarray,
    eedf_eV_inv: np.ndarray,
    *,
    relaxation_fraction: float,
    conserve_mean_energy: bool = True,
    fallback_temperature_eV: float = 2.0,
) -> np.ndarray:
    energy = np.asarray(energy_eV, dtype=float)
    widths = np.asarray(widths_eV, dtype=float)
    old = np.asarray(eedf_eV_inv, dtype=float)
    if energy.ndim != 1 or widths.shape != energy.shape or old.shape != energy.shape:
        raise ValueError("electron_electron fp_energy arrays must share one shape")
    if (
        len(energy) < 2
        or np.any(np.diff(energy) <= 0.0)
        or np.any(energy <= 0.0)
        or np.any(widths <= 0.0)
        or not np.all(np.isfinite(energy + widths + old))
    ):
        raise ValueError(
            "electron_electron fp_energy requires positive energy cells and a finite nonnegative EEDF"
        )
    alpha = float(relaxation_fraction)
    if not 0.0 <= alpha <= 1.0:
        raise ValueError("electron_electron relaxation_fraction must be in [0, 1]")
    old = prepare_relaxation_input(energy, widths, old)
    if alpha == 0.0:
        return old
    target_mean = float(np.dot(energy * widths, old))
    step = min(alpha / max(1.0 - alpha, 1e-6), 1e6) * float(np.median(widths) ** 2)
    separation = np.diff(energy)
    log_speed_ratio = 0.5 * np.log(energy[1:] / energy[:-1])

    def advance(beta: float) -> np.ndarray:
        delta = beta * separation - log_speed_ratio
        forward = step * _bernoulli(delta) / separation
        backward = step * _bernoulli(-delta) / separation
        diagonal = np.ones(len(energy))
        diagonal[:-1] += forward / widths[:-1]
        diagonal[1:] += backward / widths[1:]
        bands = np.zeros((3, len(energy)))
        bands[0, 1:] = -backward / widths[:-1]
        bands[1] = diagonal
        bands[2, :-1] = -forward / widths[1:]
        result = solve_banded((1, 1), bands, old)
        if np.any(result < 0.0) or not np.all(np.isfinite(result)):
            raise FloatingPointError("electron_electron FP step lost positivity")
        return result / np.sum(result * widths)

    if not conserve_mean_energy:
        temperature = float(fallback_temperature_eV)
        if not np.isfinite(temperature) or temperature <= 0.0:
            raise ValueError("electron_electron fallback temperature must be positive")
        return advance(1.0 / temperature)

    energy_tolerance = 2e-12 * max(target_mean, float(energy[0]))

    def balance(beta: float) -> float:
        return float(np.dot(energy * widths, advance(beta))) - target_mean

    low = 0.0
    if balance(low) < -energy_tolerance:
        raise ValueError("energy-conserving FP relaxation needs a larger energy domain")
    high = 1.5 / max(target_mean, 1e-12)
    for _ in range(60):
        if balance(high) <= 0.0:
            break
        high *= 2.0
    else:
        raise ValueError(
            "energy-conserving FP relaxation could not bracket a positive temperature"
        )
    beta = brentq(balance, low, high, xtol=1e-13, rtol=1e-13)
    result = advance(beta)
    if abs(float(np.dot(energy * widths, result)) - target_mean) > energy_tolerance:
        raise FloatingPointError("electron_electron FP energy balance did not converge")
    return result
