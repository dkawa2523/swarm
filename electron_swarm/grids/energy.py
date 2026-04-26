"""Shared energy-grid construction utilities.

The existing solvers keep their validated grid builders.  This module provides a
small common helper for new code paths and benchmark tools, including optional
threshold-aware refinement without changing solver contracts.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

import numpy as np

from electron_swarm.core.cross_sections import CrossSectionSet, ProcessType


@dataclass(frozen=True, slots=True)
class EnergyGrid:
    centers_eV: np.ndarray
    edges_eV: np.ndarray
    widths_eV: np.ndarray
    metadata: dict[str, float | int | bool | str] = field(default_factory=dict)


def cell_edges_from_centers(centers_eV: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    centers = np.asarray(centers_eV, dtype=float)
    if centers.ndim != 1 or len(centers) < 2:
        raise ValueError("Energy grid requires at least two centers")
    if np.any(np.diff(centers) <= 0.0):
        raise ValueError("Energy grid centers must be strictly increasing")
    edges = np.empty(len(centers) + 1, dtype=float)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0] = max(0.0, centers[0] - 0.5 * (centers[1] - centers[0]))
    edges[-1] = centers[-1] + 0.5 * (centers[-1] - centers[-2])
    widths = np.diff(edges)
    if np.any(widths <= 0.0):
        raise ValueError("Energy grid widths must be positive")
    return edges, widths


def _base_grid(
    min_eV: float,
    max_eV: float,
    n: int,
    spacing: str,
    *,
    linear_until_eV: float | None = None,
) -> np.ndarray:
    if n < 8:
        raise ValueError("Energy grid requires at least 8 cells")
    emin = float(min_eV)
    emax = float(max_eV)
    if emax <= emin:
        raise ValueError("max_eV must be larger than min_eV")
    spacing = spacing.lower()
    if spacing == "linear":
        return np.linspace(emin, emax, n)
    if spacing == "quadratic":
        x = np.linspace(0.0, 1.0, n)
        return emin + (emax - emin) * x * x
    if spacing == "log":
        return np.geomspace(max(emin, 1.0e-12), emax, n)
    if spacing == "log_linear":
        split = min(max(linear_until_eV or 2.0, emin), emax)
        n_linear = max(3, min(n - 3, int(0.35 * n)))
        n_log = n - n_linear + 1
        low = np.linspace(emin, split, n_linear)
        high = np.geomspace(max(split, 1.0e-12), emax, n_log)
        return np.unique(np.concatenate([low, high[1:]]))
    raise ValueError(f"Unsupported energy grid spacing: {spacing}")


def process_thresholds(
    cross_sections: CrossSectionSet | None,
    include_types: Iterable[ProcessType] | None = None,
) -> tuple[float, ...]:
    if cross_sections is None:
        return ()
    active = set(
        include_types
        or (
            ProcessType.EXCITATION,
            ProcessType.IONIZATION,
            ProcessType.ATTACHMENT,
            ProcessType.SUPERELASTIC,
        )
    )
    values = []
    for process in cross_sections.processes:
        if process.process_type not in active:
            continue
        if process.threshold_eV is not None and np.isfinite(process.threshold_eV):
            values.append(float(process.threshold_eV))
    return tuple(sorted({value for value in values if value >= 0.0}))


def refine_with_thresholds(
    centers_eV: np.ndarray,
    thresholds_eV: Iterable[float],
    *,
    padding_eV: float = 0.15,
    points_per_threshold: int = 8,
    max_extra_points: int = 120,
) -> np.ndarray:
    centers = np.asarray(centers_eV, dtype=float)
    additions = []
    for threshold in thresholds_eV:
        th = float(threshold)
        lo = max(0.0, th - padding_eV)
        hi = th + padding_eV
        if hi < centers[0] or lo > centers[-1]:
            continue
        additions.append(np.linspace(lo, hi, max(2, points_per_threshold)))
    if not additions:
        return centers
    extra = np.unique(np.concatenate(additions))
    if len(extra) > max_extra_points:
        idx = np.linspace(0, len(extra) - 1, max_extra_points).astype(int)
        extra = extra[idx]
    merged = np.unique(np.concatenate([centers, extra]))
    return merged[(merged >= centers[0]) & (merged <= centers[-1])]


def build_energy_grid(
    *,
    min_eV: float,
    max_eV: float,
    n: int,
    spacing: str = "linear",
    linear_until_eV: float | None = None,
    cross_sections: CrossSectionSet | None = None,
    refine: bool = False,
    threshold_padding_eV: float = 0.15,
    points_per_threshold: int = 8,
    max_extra_points: int = 120,
) -> EnergyGrid:
    centers = _base_grid(
        min_eV,
        max_eV,
        n,
        spacing,
        linear_until_eV=linear_until_eV,
    )
    thresholds = process_thresholds(cross_sections)
    if refine and thresholds:
        centers = refine_with_thresholds(
            centers,
            thresholds,
            padding_eV=threshold_padding_eV,
            points_per_threshold=points_per_threshold,
            max_extra_points=max_extra_points,
        )
    edges, widths = cell_edges_from_centers(centers)
    return EnergyGrid(
        centers_eV=centers,
        edges_eV=edges,
        widths_eV=widths,
        metadata={
            "spacing": spacing,
            "n_cells": int(len(centers)),
            "min_eV": float(centers[0]),
            "max_eV": float(centers[-1]),
            "threshold_refined": bool(refine and thresholds),
            "n_thresholds": int(len(thresholds)),
        },
    )
