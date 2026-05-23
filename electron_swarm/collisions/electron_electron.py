"""Electron-electron relaxation target helpers for product schema v2."""

from __future__ import annotations

import numpy as np


_EPS = 1.0e-300


def _normalize(values: np.ndarray, widths: np.ndarray) -> np.ndarray:
    total = float(np.sum(values * widths))
    if total <= 0.0 or not np.isfinite(total):
        raise ValueError("Cannot normalize non-positive EEDF target")
    return values / total


def _maxwellian_shape(energy_eV: np.ndarray, kT_eV: float) -> np.ndarray:
    energy = np.asarray(energy_eV, dtype=float)
    kT = max(float(kT_eV), 1.0e-12)
    return np.sqrt(np.maximum(energy, 0.0)) * np.exp(-np.maximum(energy, 0.0) / kT)


def mean_energy_eV(energy_eV: np.ndarray, widths_eV: np.ndarray, eedf_eV_inv: np.ndarray) -> float:
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
    target_mean = max(float(mean_energy_target_eV), 0.0)

    def distribution(kT_eV: float) -> np.ndarray:
        return _normalize(np.clip(_maxwellian_shape(energy, kT_eV), 0.0, None), widths)

    lo = 1.0e-8
    hi = max(2.0 * target_mean / 3.0, lo)
    for _ in range(80):
        if mean_energy_eV(energy, widths, distribution(hi)) >= target_mean:
            break
        hi *= 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if mean_energy_eV(energy, widths, distribution(mid)) < target_mean:
            lo = mid
        else:
            hi = mid
    return distribution(hi)


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
    else:
        target_mean = 1.5 * float(fallback_temperature_eV)
    return maxwellian_eedf_like(energy_eV, widths_eV, target_mean)
