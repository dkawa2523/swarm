"""Validate and normalize finite-volume EEDF source slices."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

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
    if np.any(~np.isfinite(widths)) or np.any(widths <= 0.0):
        raise FunctionEedfError("source energy widths must be positive")
    if np.any(~np.isfinite(density)) or np.any(density < 0.0):
        raise FunctionEedfError("source EEDF must be finite and nonnegative")
    left = np.maximum(0.0, centers - 0.5 * widths)
    right = centers + 0.5 * widths
    overlap = right[:-1] - left[1:]
    scale = np.maximum.reduce(
        [np.ones_like(overlap), np.abs(right[:-1]), np.abs(left[1:])]
    )
    if np.any(overlap > 1.0e-10 * scale):
        raise FunctionEedfError("source EEDF bins overlap")
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
