"""Energy-grid construction and initialization for the two-term solver."""

from __future__ import annotations

import numpy as np

from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.core.solver_configs import TwoTermInternalConfig
from electron_swarm.grids.energy import build_energy_grid
from electron_swarm.solvers.boltzmann_common.grid import cell_edges_from_centers


def make_energy_grid(
    cfg: TwoTermInternalConfig,
    *,
    max_eV_override: float | None = None,
    n_override: int | None = None,
    cross_sections: CrossSectionSet | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    grid = cfg.energy_grid
    n = int(n_override or grid.n)
    emin = max(float(grid.min_eV), 0.0)
    emax = float(max_eV_override if max_eV_override is not None else grid.max_eV)
    if n < 8:
        raise ValueError("Boltzmann energy grid requires at least 8 cells")
    if emax <= emin:
        raise ValueError("Boltzmann energy grid max_eV must be larger than min_eV")

    if grid.refine.enabled:
        refined = build_energy_grid(
            min_eV=emin,
            max_eV=emax,
            n=n,
            spacing=grid.spacing,
            cross_sections=cross_sections,
            refine=True,
            threshold_padding_eV=grid.refine.threshold_padding_eV,
            points_per_threshold=grid.refine.points_per_threshold,
            max_extra_points=grid.refine.max_extra_points,
        )
        energy = refined.centers_eV
        edges = refined.edges_eV
        widths = refined.widths_eV
    elif grid.spacing == "linear":
        energy = np.linspace(emin, emax, n)
        edges, widths = cell_edges_from_centers(energy)
    elif grid.spacing == "quadratic":
        x = np.linspace(0.0, 1.0, n)
        energy = emin + (emax - emin) * x * x
        edges, widths = cell_edges_from_centers(energy)
    elif grid.spacing == "log":
        energy = np.geomspace(max(emin, 1.0e-8), emax, n)
        edges, widths = cell_edges_from_centers(energy)
    else:
        raise ValueError(f"Unsupported energy spacing: {grid.spacing}")
    return energy, edges, widths


def grid_metadata(
    cfg: TwoTermInternalConfig, energy: np.ndarray
) -> dict[str, bool | float | int | str]:
    grid = cfg.energy_grid
    return {
        "grid_n_cells": int(len(energy)),
        "grid_min_eV": float(np.min(energy)) if len(energy) else float("nan"),
        "grid_max_eV": float(np.max(energy)) if len(energy) else float("nan"),
        "grid_spacing": grid.spacing,
        "threshold_refined": bool(grid.refine.enabled and len(energy) > grid.n),
    }


def maxwell_eedf(energy: np.ndarray, kT_eV: float) -> np.ndarray:
    """Return a Maxwell-Boltzmann EEDF seed before grid normalization."""

    kT = max(float(kT_eV), 1.0e-8)
    f = np.sqrt(np.maximum(energy, 0.0)) * np.exp(-np.maximum(energy, 0.0) / kT)
    return np.clip(f, 0.0, None)
