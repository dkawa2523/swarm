"""Statistics and conservative EEDF-grid operations for workflow aggregation."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import numpy as np
from scipy.stats import t as student_t


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
    estimate_status: str
    ci95_critical_value: float | None


def summarize_replicates(
    values: Iterable[object],
    *,
    zero_is_censored: bool = False,
) -> ScalarStats:
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
        return ScalarStats(
            None,
            None,
            None,
            None,
            None,
            None,
            0,
            False,
            "unavailable",
            None,
        )
    mean = float(np.mean(finite))
    if count == 1:
        return ScalarStats(
            mean,
            None,
            None,
            None,
            None,
            None,
            1,
            False,
            "censored_all_zero" if zero_is_censored and mean == 0.0 else "estimate",
            None,
        )
    sample_stddev = float(np.std(finite, ddof=1))
    if zero_is_censored and bool(np.all(finite == 0.0)):
        # A zero event count does not establish a zero physical rate.  The
        # workflow database has no collision exposure from which a numerical
        # Poisson upper bound can be reconstructed, so retain the observed
        # mean but make uncertainty/CI unavailable and force the quality gate
        # to treat the rate as censored.
        return ScalarStats(
            mean=mean,
            sample_stddev=sample_stddev,
            standard_error=None,
            relative_standard_error=None,
            ci95_low=None,
            ci95_high=None,
            valid_replicates=count,
            uncertainty_available=False,
            estimate_status="censored_all_zero",
            ci95_critical_value=None,
        )
    standard_error = float(sample_stddev / math.sqrt(count))
    relative_standard_error = (
        0.0 if mean == 0.0 and standard_error == 0.0 else None
    )
    if mean != 0.0:
        relative_standard_error = float(abs(standard_error / mean))
    critical_value = float(student_t.ppf(0.975, df=count - 1))
    half_width = critical_value * standard_error
    return ScalarStats(
        mean=mean,
        sample_stddev=sample_stddev,
        standard_error=standard_error,
        relative_standard_error=relative_standard_error,
        ci95_low=float(mean - half_width),
        ci95_high=float(mean + half_width),
        valid_replicates=count,
        uncertainty_available=True,
        estimate_status="estimate",
        ci95_critical_value=critical_value,
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


def energy_edges_from_cells(
    centers: Iterable[float],
    widths: Iterable[float],
) -> np.ndarray:
    """Return exact finite-volume edges encoded by cell centers and widths.

    Histogram densities are cell averages, not point samples.  Reconstructing
    their edges from distances between neighboring centers changes the cell
    volume whenever the grid spacing changes.  The stored width is therefore
    part of the physical definition of every EEDF bin.
    """

    center_values = np.asarray(list(centers), dtype=float)
    width_values = np.asarray(list(widths), dtype=float)
    if center_values.ndim != 1 or width_values.ndim != 1:
        raise ValueError("EEDF cell centers and widths must be one-dimensional")
    if len(center_values) == 0 or len(center_values) != len(width_values):
        raise ValueError(
            "EEDF cell centers and widths must have equal nonzero length"
        )
    if np.any(~np.isfinite(center_values)) or np.any(~np.isfinite(width_values)):
        raise ValueError("EEDF cell centers and widths must be finite")
    if np.any(width_values <= 0.0):
        raise ValueError("EEDF cell widths must be positive")
    if len(center_values) > 1 and np.any(np.diff(center_values) <= 0.0):
        raise ValueError("EEDF cell centers must be strictly increasing")

    left = center_values - 0.5 * width_values
    right = center_values + 0.5 * width_values
    scale = np.maximum.reduce(
        (
            np.ones_like(center_values),
            np.abs(center_values),
            width_values,
        )
    )
    tolerance = 64.0 * np.finfo(float).eps * scale
    if left[0] < -tolerance[0]:
        raise ValueError("EEDF cells must not extend below zero energy")
    if left[0] < 0.0:
        left[0] = 0.0

    if len(center_values) > 1:
        gap = left[1:] - right[:-1]
        boundary_tolerance = np.maximum(tolerance[1:], tolerance[:-1])
        if np.any(np.abs(gap) > boundary_tolerance):
            raise ValueError("EEDF cells must be contiguous and nonoverlapping")
        shared = 0.5 * (right[:-1] + left[1:])
    else:
        shared = np.asarray([], dtype=float)

    edges = np.concatenate(([left[0]], shared, [right[-1]]))
    if np.any(np.diff(edges) <= 0.0):
        raise ValueError("EEDF cell edges must be strictly increasing")
    return edges


def energy_edges_from_nodes(nodes: Iterable[float]) -> np.ndarray:
    """Return control-volume edges for a nodal EEDF representation."""

    values = np.asarray(list(nodes), dtype=float)
    if values.ndim != 1 or len(values) == 0:
        raise ValueError("EEDF nodes must be a nonempty one-dimensional array")
    if np.any(~np.isfinite(values)):
        raise ValueError("EEDF nodes must be finite")
    if len(values) > 1 and np.any(np.diff(values) <= 0.0):
        raise ValueError("EEDF nodes must be strictly increasing")
    if len(values) == 1:
        return np.asarray([max(0.0, values[0] - 0.5), values[0] + 0.5])
    edges = np.empty(len(values) + 1, dtype=float)
    edges[1:-1] = 0.5 * (values[:-1] + values[1:])
    edges[0] = max(0.0, values[0] - 0.5 * (values[1] - values[0]))
    edges[-1] = values[-1] + 0.5 * (values[-1] - values[-2])
    return edges
