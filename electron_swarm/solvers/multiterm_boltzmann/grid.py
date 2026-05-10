"""Energy-grid utilities for multi-term Boltzmann backends."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.constants import ELECTRON_MASS_KG, EV_TO_J
from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.grids.energy import build_energy_grid as build_shared_energy_grid


def electron_speed_m_s(energy_eV: np.ndarray | float) -> np.ndarray:
    energy = np.asarray(energy_eV, dtype=float)
    return np.sqrt(np.maximum(2.0 * EV_TO_J * energy / ELECTRON_MASS_KG, 0.0))


@dataclass(frozen=True, slots=True)
class EnergyGrid:
    centers_eV: np.ndarray
    edges_eV: np.ndarray

    def __post_init__(self) -> None:
        centers = np.asarray(self.centers_eV, dtype=float)
        edges = np.asarray(self.edges_eV, dtype=float)
        if centers.ndim != 1 or edges.ndim != 1:
            raise ValueError("Energy grid arrays must be one-dimensional")
        if len(edges) != len(centers) + 1:
            raise ValueError("Energy grid edges length must be centers length + 1")
        if np.any(np.diff(edges) <= 0.0):
            raise ValueError("Energy grid edges must be strictly increasing")
        if np.any(centers <= edges[:-1]) or np.any(centers >= edges[1:]):
            raise ValueError("Energy grid centers must lie inside their cells")
        object.__setattr__(self, "centers_eV", centers)
        object.__setattr__(self, "edges_eV", edges)

    @property
    def widths_eV(self) -> np.ndarray:
        return np.diff(self.edges_eV)

    @property
    def n_cells(self) -> int:
        return len(self.centers_eV)

    @property
    def speeds_m_s(self) -> np.ndarray:
        return electron_speed_m_s(self.centers_eV)

    @classmethod
    def linear(cls, emin_eV: float, emax_eV: float, n_cells: int) -> "EnergyGrid":
        edges = np.linspace(float(emin_eV), float(emax_eV), int(n_cells) + 1)
        centers = 0.5 * (edges[:-1] + edges[1:])
        return cls(centers, edges)

    @classmethod
    def log(cls, emin_eV: float, emax_eV: float, n_cells: int) -> "EnergyGrid":
        if emin_eV <= 0.0:
            raise ValueError("log grid requires positive min_eV")
        edges = np.geomspace(float(emin_eV), float(emax_eV), int(n_cells) + 1)
        centers = np.sqrt(edges[:-1] * edges[1:])
        return cls(centers, edges)

    @classmethod
    def log_linear(
        cls,
        emin_eV: float,
        emax_eV: float,
        n_cells: int,
        linear_until_eV: float = 2.0,
    ) -> "EnergyGrid":
        if not (emin_eV < linear_until_eV < emax_eV):
            return cls.log(max(emin_eV, 1.0e-8), emax_eV, n_cells)
        n_cells = int(n_cells)
        n_lin = max(
            8,
            min(
                n_cells // 3,
                int(n_cells * linear_until_eV / max(emax_eV, linear_until_eV)),
            ),
        )
        n_log = n_cells - n_lin
        e1 = np.linspace(float(emin_eV), float(linear_until_eV), n_lin + 1)
        e2 = np.geomspace(float(linear_until_eV), float(emax_eV), n_log + 1)[1:]
        edges = np.concatenate([e1, e2])
        centers = 0.5 * (edges[:-1] + edges[1:])
        return cls(centers, edges)

    def integrate_energy_pdf(self, values: np.ndarray) -> float:
        return float(np.sum(np.asarray(values, dtype=float) * self.widths_eV))

    def normalize_energy_pdf(self, values: np.ndarray) -> np.ndarray:
        values = np.asarray(values, dtype=float).copy()
        if not np.all(np.isfinite(values)):
            raise ValueError("cannot normalize non-finite EEDF")
        if np.any(values < -1.0e-300):
            raise ValueError("cannot normalize EEDF with negative values")
        values = np.clip(values, 0.0, None)
        norm = self.integrate_energy_pdf(values)
        if not np.isfinite(norm) or norm <= 0.0:
            raise ValueError("cannot normalize non-positive EEDF")
        return values / norm


def make_energy_grid(
    config: SwarmConfig,
    cross_sections: CrossSectionSet | None = None,
) -> EnergyGrid:
    grid = config.multiterm_boltzmann.energy_grid
    if grid.refine.enabled:
        shared = build_shared_energy_grid(
            min_eV=grid.min_eV,
            max_eV=grid.max_eV,
            n=grid.n,
            spacing=grid.spacing,
            linear_until_eV=grid.linear_until_eV,
            cross_sections=cross_sections,
            refine=True,
            threshold_padding_eV=grid.refine.threshold_padding_eV,
            points_per_threshold=grid.refine.points_per_threshold,
            max_extra_points=grid.refine.max_extra_points,
        )
        return EnergyGrid(shared.centers_eV, shared.edges_eV)
    if grid.spacing == "linear":
        return EnergyGrid.linear(grid.min_eV, grid.max_eV, grid.n)
    if grid.spacing == "log":
        return EnergyGrid.log(grid.min_eV, grid.max_eV, grid.n)
    return EnergyGrid.log_linear(
        grid.min_eV,
        grid.max_eV,
        grid.n,
        linear_until_eV=grid.linear_until_eV,
    )
