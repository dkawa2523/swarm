"""Small numerical helpers shared by runtime, diagnostics, and tools."""

from __future__ import annotations

import numpy as np


def widths_from_centers(energy_eV: np.ndarray) -> np.ndarray:
    """Return finite-volume cell widths from increasing cell centers."""

    energy = np.asarray(energy_eV, dtype=float)
    if energy.ndim != 1 or len(energy) == 0:
        return np.asarray([], dtype=float)
    if len(energy) == 1:
        return np.ones_like(energy)
    edges = np.empty(len(energy) + 1, dtype=float)
    edges[1:-1] = 0.5 * (energy[:-1] + energy[1:])
    edges[0] = max(0.0, energy[0] - 0.5 * (energy[1] - energy[0]))
    edges[-1] = energy[-1] + 0.5 * (energy[-1] - energy[-2])
    return np.diff(edges)


def eepf_from_eedf(energy_eV: np.ndarray, eedf: np.ndarray) -> np.ndarray:
    return np.asarray(eedf, dtype=float) / np.sqrt(
        np.maximum(np.asarray(energy_eV, dtype=float), 1.0e-30)
    )


def f0_from_eedf(energy_eV: np.ndarray, eedf: np.ndarray) -> np.ndarray:
    return eepf_from_eedf(energy_eV, eedf)


def eedf_from_f0(energy_eV: np.ndarray, f0: np.ndarray) -> np.ndarray:
    return np.asarray(f0, dtype=float) * np.sqrt(
        np.maximum(np.asarray(energy_eV, dtype=float), 1.0e-30)
    )
