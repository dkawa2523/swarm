"""Solver-neutral energy-grid geometry used by Boltzmann solvers."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from electron_swarm.core.constants import ELECTRON_MASS_KG, EV_TO_J


@dataclass(frozen=True, slots=True)
class KineticGrid:
    energy_eV: np.ndarray
    edges_eV: np.ndarray
    widths_eV: np.ndarray
    speed_m_s: np.ndarray


def electron_speed_m_s(energy_eV: np.ndarray) -> np.ndarray:
    return np.sqrt(np.maximum(2.0 * EV_TO_J * energy_eV / ELECTRON_MASS_KG, 0.0))


def cell_edges_from_centers(centers: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    edges = np.empty(len(centers) + 1)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0] = max(0.0, centers[0] - 0.5 * (centers[1] - centers[0]))
    edges[-1] = centers[-1] + 0.5 * (centers[-1] - centers[-2])
    widths = np.diff(edges)
    if np.any(widths <= 0.0):
        raise ValueError("Energy grid must be strictly increasing")
    return edges, widths


__all__ = ["KineticGrid", "cell_edges_from_centers", "electron_speed_m_s"]
