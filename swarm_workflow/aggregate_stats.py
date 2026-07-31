"""Statistics and conservative EEDF-grid operations for workflow aggregation."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import numpy as np


@dataclass(frozen=True, slots=True)
class ScalarStats:
    mean: float | None
    sample_stddev: float | None
    standard_error: float | None
    relative_standard_error: float | None
    ci95_low: float | None
    ci95_high: float | None
    valid_replicates: int
    uncertainty_available: bool


def summarize_replicates(values: Iterable[object]) -> ScalarStats:
    finite = np.asarray(
        [
            float(value)
            for value in values
            if value is not None and math.isfinite(float(value))
        ],
        dtype=float,
    )
    count = int(len(finite))
    if count == 0:
        return ScalarStats(None, None, None, None, None, None, 0, False)
    mean = float(np.mean(finite))
    if count == 1:
        return ScalarStats(mean, None, None, None, None, None, 1, False)
    sample_stddev = float(np.std(finite, ddof=1))
    standard_error = float(sample_stddev / math.sqrt(count))
    relative_standard_error = (
        0.0 if mean == 0.0 and standard_error == 0.0 else None
    )
    if mean != 0.0:
        relative_standard_error = float(abs(standard_error / mean))
    half_width = 1.96 * standard_error
    return ScalarStats(
        mean=mean,
        sample_stddev=sample_stddev,
        standard_error=standard_error,
        relative_standard_error=relative_standard_error,
        ci95_low=float(mean - half_width),
        ci95_high=float(mean + half_width),
        valid_replicates=count,
        uncertainty_available=True,
    )


def conservative_rebin_probability_mass(
    source_edges: Iterable[float],
    source_probability_mass: Iterable[float],
    target_edges: Iterable[float],
) -> np.ndarray:
    source_edges_arr = np.asarray(list(source_edges), dtype=float)
    target_edges_arr = np.asarray(list(target_edges), dtype=float)
    source_mass = np.asarray(list(source_probability_mass), dtype=float)
    if len(source_edges_arr) != len(source_mass) + 1:
        raise ValueError("source_edges must have one more entry than source mass")
    if len(target_edges_arr) < 2:
        raise ValueError("target_edges must contain at least two edges")
    if np.any(np.diff(source_edges_arr) <= 0.0):
        raise ValueError("source_edges must be strictly increasing")
    if np.any(np.diff(target_edges_arr) <= 0.0):
        raise ValueError("target_edges must be strictly increasing")
    if np.any(~np.isfinite(source_mass)):
        raise ValueError("source probability mass must be finite")

    rebinned = np.zeros(len(target_edges_arr) - 1, dtype=float)
    source_index = 0
    target_index = 0
    while source_index < len(source_mass) and target_index < len(rebinned):
        source_left = source_edges_arr[source_index]
        source_right = source_edges_arr[source_index + 1]
        target_left = target_edges_arr[target_index]
        target_right = target_edges_arr[target_index + 1]
        overlap = min(source_right, target_right) - max(source_left, target_left)
        if overlap > 0.0:
            rebinned[target_index] += (
                source_mass[source_index] * overlap / (source_right - source_left)
            )
        if source_right <= target_right:
            source_index += 1
        else:
            target_index += 1
    return rebinned


def energy_edges_from_centers(centers: Iterable[float]) -> np.ndarray:
    values = np.asarray(list(centers), dtype=float)
    if len(values) == 0:
        return np.asarray([], dtype=float)
    if len(values) == 1:
        return np.asarray([max(0.0, values[0] - 0.5), values[0] + 0.5], dtype=float)
    edges = np.empty(len(values) + 1, dtype=float)
    edges[1:-1] = 0.5 * (values[:-1] + values[1:])
    edges[0] = max(0.0, values[0] - 0.5 * (values[1] - values[0]))
    edges[-1] = values[-1] + 0.5 * (values[-1] - values[-2])
    return edges
