"""Validate and normalize finite-volume EEDF source slices."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from electron_swarm.grids.energy import reported_cell_edges

from .contracts import FunctionEedfError


@dataclass(frozen=True, slots=True)
class _SourceSlice:
    mean_energy_eV: float
    left_edges_eV: np.ndarray
    right_edges_eV: np.ndarray
    bin_f0_eV_m32: np.ndarray
    source_mean_relative_error: float


def _prepare_source_slice(
    mean_energy_eV: float,
    values: list[tuple[float, float, float]],
) -> _SourceSlice:
    centers = np.asarray([item[0] for item in values], dtype=float)
    widths = np.asarray([item[1] for item in values], dtype=float)
    density = np.asarray([item[2] for item in values], dtype=float)
    if len(centers) < 2 or np.any(np.diff(centers) <= 0.0):
        raise FunctionEedfError("source energy grid must be strictly increasing")
    try:
        edges = reported_cell_edges(centers, widths)
    except ValueError as exc:
        raise FunctionEedfError(f"source EEDF has an invalid cell grid: {exc}") from exc
    if np.any(~np.isfinite(density)) or np.any(density < 0.0):
        raise FunctionEedfError("source EEDF must be finite and nonnegative")
    left = edges[:-1]
    right = edges[1:]
    masses = density * widths
    source_norm = float(np.sum(masses))
    if not math.isfinite(source_norm) or source_norm <= 0.0:
        raise FunctionEedfError("source EEDF has zero mass")
    source_mean = float(np.sum(centers * masses) / source_norm)
    masses /= source_norm
    measure = (2.0 / 3.0) * (np.power(right, 1.5) - np.power(left, 1.5))
    return _SourceSlice(
        mean_energy_eV=float(mean_energy_eV),
        left_edges_eV=left,
        right_edges_eV=right,
        bin_f0_eV_m32=masses / measure,
        source_mean_relative_error=abs(source_mean - mean_energy_eV)
        / max(abs(mean_energy_eV), 1.0e-12),
    )
