"""Construct and project dimensionless C1 Function-EEDF families."""

from __future__ import annotations

import math

import numpy as np

from .adaptive_grid import (
    adapt_energy_axis,
    audit_energy_axis,
    closure_energy_support_max,
    initial_energy_axis,
    mean_axis_refinement,
    project_rows,
    rate_error_scales,
    table_projection_error,
)
from .contracts import (
    SOURCE_MEAN_RELATIVE_ERROR_LIMIT,
    C1FunctionEedf,
    CollisionRateKernel,
    ComsolFunctionEedfGrid,
    FunctionEedfError,
)
from .moments import (
    _pchip_minimum,
    _tilt_to_unit_mean,
    pchip_weighted_moments,
    piecewise_linear_weighted_moments,
)
from .source import _SourceSlice, _prepare_source_slice


SHAPE_TOTAL_VARIATION_TOLERANCE = 5.0e-3
RATE_SCALED_ERROR_TOLERANCE = 5.0e-3
RATE_IMPORTANCE_FRACTION = 1.0e-12
# Rates below this fraction of a kernel's family-wide peak are controlled in
# absolute, rather than relative, terms while refining the mean-energy axis.
# A relative norm in the vanishing-rate limit needlessly recreates a dense
# uniform table without improving the plasma source terms.
MEAN_AXIS_RATE_IMPORTANCE_FRACTION = 1.0e-3


def build_c1_function_eedf(
    grouped_bins: dict[float, list[tuple[float, float, float]]],
) -> C1FunctionEedf:
    """Build nonnegative C1 shapes from ``(center, width, EEDF)`` bins.

    For ``x=E/m`` and ``h=m**(3/2) f0``, each shape is constrained to
    ``integral(sqrt(x) h dx)=integral(x**(3/2) h dx)=1``.  A convex C1
    blend of these shapes therefore remains normalized and has mean ``m``.
    """

    slices = [
        _prepare_source_slice(mean, sorted(values))
        for mean, values in sorted(grouped_bins.items())
    ]
    if len(slices) < 2:
        raise FunctionEedfError("Function EEDF requires at least two mean energies")
    means = np.asarray([slice_.mean_energy_eV for slice_ in slices], dtype=float)
    if np.any(~np.isfinite(means)) or np.any(means <= 0.0):
        raise FunctionEedfError("mean-energy anchors must be positive and finite")
    if np.any(np.diff(means) <= 0.0):
        raise FunctionEedfError("mean-energy anchors must be strictly increasing")
    source_mean_error = max(slice_.source_mean_relative_error for slice_ in slices)
    if source_mean_error > SOURCE_MEAN_RELATIVE_ERROR_LIMIT:
        raise FunctionEedfError(
            "source EEDF mean energy differs from its case mean: "
            f"relative error {source_mean_error:.3e} exceeds "
            f"{SOURCE_MEAN_RELATIVE_ERROR_LIMIT:.3e}"
        )

    x_grid = _dimensionless_grid(slices)
    shapes: list[np.ndarray] = []
    norm_error = 0.0
    mean_error = 0.0
    minimum = math.inf
    for slice_ in slices:
        mean = slice_.mean_energy_eV
        left = slice_.left_edges_eV / mean
        right = slice_.right_edges_eV / mean
        bin_shape = (mean**1.5) * slice_.bin_f0_eV_m32
        sampled = _sample_bins(left, right, bin_shape, x_grid)
        matched = _tilt_to_unit_mean(x_grid, sampled)
        normalization, first_moment = pchip_weighted_moments(x_grid, matched)
        norm_error = max(norm_error, abs(normalization - 1.0))
        mean_error = max(
            mean_error,
            abs(first_moment / max(normalization, 1.0e-300) - 1.0),
        )
        minimum = min(minimum, _pchip_minimum(x_grid, matched))
        shapes.append(matched)

    return C1FunctionEedf(
        mean_energies_eV=means,
        dimensionless_energy=x_grid,
        shapes=tuple(shapes),
        source_energy_max_eV=max(float(slice_.right_edges_eV[-1]) for slice_ in slices),
        source_mean_relative_error_max=source_mean_error,
        normalization_error_max=norm_error,
        unit_mean_error_max=mean_error,
        minimum_value=minimum,
    )


def project_c1_function_eedf_to_comsol_grid(
    representation: C1FunctionEedf,
    *,
    rate_kernels: tuple[CollisionRateKernel, ...] = (),
    shape_total_variation_tolerance: float = SHAPE_TOTAL_VARIATION_TOLERANCE,
    rate_scaled_error_tolerance: float = RATE_SCALED_ERROR_TOLERANCE,
) -> ComsolFunctionEedfGrid:
    """Adapt the shared energy axis until physical moments and rates are preserved.

    The continuous C1 family is the projection reference.  Both tensor axes
    are refined by probability-distribution total variation and, when
    supplied, every physical ``sigma*v`` kernel.  The mean-energy axis retains
    every solver anchor and adds only the C1 rows needed by COMSOL's linear
    interpolator.  Each materialized row has exact zeroth and first moments,
    so interpolation between rows preserves normalization and the requested
    mean.  This replaces the old fixed ``2049 x 257`` sampling and applies
    identically to every solver source.
    """

    if (
        not math.isfinite(shape_total_variation_tolerance)
        or shape_total_variation_tolerance <= 0.0
        or not math.isfinite(rate_scaled_error_tolerance)
        or rate_scaled_error_tolerance <= 0.0
    ):
        raise FunctionEedfError("adaptive Function-EEDF tolerances must be positive")
    x_values = representation.dimensionless_energy
    positive_x = x_values[x_values > 0.0]
    anchors = representation.mean_energies_eV
    if not len(positive_x) or not len(anchors):
        raise FunctionEedfError("Function-EEDF support is empty")
    minimum_energy = float(np.min(anchors[:, None] * positive_x[0]))
    # The common dimensionless grid contains the largest E/m value from every
    # source slice. Multiplying that global value by the largest mean energy
    # can therefore create a large, entirely zero physical tail. Bound the
    # active COMSOL table by the source grid while retaining any larger
    # nonzero support introduced by interpolation between adjacent shapes.
    maximum_energy = max(
        representation.source_energy_max_eV,
        closure_energy_support_max(representation),
    )
    if (
        not math.isfinite(minimum_energy)
        or not math.isfinite(maximum_energy)
        or minimum_energy <= 0.0
        or maximum_energy <= minimum_energy
    ):
        raise FunctionEedfError("Function-EEDF physical energy range is invalid")
    audit_energy = audit_energy_axis(
        minimum_energy,
        maximum_energy,
        rate_kernels,
    )
    rate_scales = rate_error_scales(
        representation,
        audit_energy,
        anchors,
        rate_kernels,
        importance_fraction=RATE_IMPORTANCE_FRACTION,
    )
    mean_rate_scales = rate_error_scales(
        representation,
        audit_energy,
        anchors,
        rate_kernels,
        importance_fraction=MEAN_AXIS_RATE_IMPORTANCE_FRACTION,
    )
    means = anchors.copy()
    energies = initial_energy_axis(
        minimum_energy,
        maximum_energy,
        rate_kernels,
    )
    reference_cache: dict[float, np.ndarray] = {}
    mean_shape_error = math.inf
    mean_rate_error: float | None = math.inf if rate_kernels else None
    for _iteration in range(32):
        energies = adapt_energy_axis(
            representation,
            energies,
            audit_energy,
            means,
            rate_kernels,
            rate_scales,
            shape_tolerance=shape_total_variation_tolerance,
            rate_tolerance=rate_scaled_error_tolerance,
        )
        values = project_rows(representation, energies, means)
        additions, mean_shape_error, mean_rate_error = mean_axis_refinement(
            representation,
            energies,
            audit_energy,
            means,
            values,
            rate_kernels,
            mean_rate_scales,
            reference_cache,
            shape_tolerance=shape_total_variation_tolerance,
            rate_tolerance=rate_scaled_error_tolerance,
        )
        if not additions:
            break
        means = np.unique(np.concatenate((means, np.asarray(additions))))
        if len(means) > 513:
            raise FunctionEedfError(
                "adaptive Function-EEDF mean-energy grid exceeded 513 points"
            )
    else:
        raise FunctionEedfError(
            "adaptive Function-EEDF mean-energy refinement did not converge "
            f"({len(means)} points, TV={mean_shape_error:.3e}, "
            f"rate={mean_rate_error if mean_rate_error is not None else 'not_requested'})"
        )
    shape_error, rate_error = table_projection_error(
        representation,
        energies,
        means,
        values,
        audit_energy,
        rate_kernels,
        rate_scales,
    )
    if shape_error > shape_total_variation_tolerance or (
        rate_error is not None and rate_error > rate_scaled_error_tolerance
    ):
        raise FunctionEedfError(
            "adaptive Function-EEDF projection did not meet its shape/rate tolerance"
        )

    normalization_error = 0.0
    mean_error = 0.0
    minimum = math.inf
    for mean, projected in zip(means, values, strict=True):
        normalization, first_moment = piecewise_linear_weighted_moments(
            energies,
            projected,
        )
        normalization_error = max(normalization_error, abs(normalization - 1.0))
        mean_error = max(
            mean_error,
            abs(first_moment / normalization - mean) / mean,
        )
        minimum = min(minimum, float(np.min(projected)))
    if (
        np.any(~np.isfinite(values))
        or minimum < 0.0
        or normalization_error > 1.0e-8
        or mean_error > 1.0e-8
    ):
        raise FunctionEedfError(
            "projected COMSOL Function EEDF failed its moment contract"
        )
    return ComsolFunctionEedfGrid(
        electron_energies_eV=energies,
        mean_energies_eV=means,
        values_eV_m32=values,
        normalization_error_max=normalization_error,
        mean_energy_relative_error_max=mean_error,
        minimum_value=minimum,
        shape_total_variation_error_max=shape_error,
        rate_scaled_error_max=rate_error,
        rate_kernel_count=len(rate_kernels),
        mean_axis_shape_total_variation_error_max=mean_shape_error,
        mean_axis_rate_scaled_error_max=mean_rate_error,
    )


def _dimensionless_grid(slices: list[_SourceSlice]) -> np.ndarray:
    left_edges = [slice_.left_edges_eV / slice_.mean_energy_eV for slice_ in slices]
    right_edges = [slice_.right_edges_eV / slice_.mean_energy_eV for slice_ in slices]
    max_x = max(float(right[-1]) for right in right_edges)
    positive = np.concatenate(
        [
            np.concatenate((left[left > 0.0], right[right > 0.0]))
            for left, right in zip(left_edges, right_edges, strict=True)
        ]
    )
    if max_x <= 1.0 or not len(positive):
        raise FunctionEedfError("source support must extend beyond its mean")
    point_count = min(2400, max(600, 2 * max(len(left) for left in left_edges)))
    first = max(float(np.min(positive)) * 1.0e-3, max_x * 1.0e-14)
    grid = set(np.geomspace(first, max_x, point_count - 1).tolist())
    grid.add(0.0)
    terminal_widths: list[float] = []
    for left, right in zip(left_edges, right_edges, strict=True):
        edge = float(right[-1])
        width = float(right[-1] - left[-1])
        terminal_widths.append(width)
        grid.add(edge)
    # Source support ends at the final finite-volume face.  A zero guard after
    # that face makes zero extrapolation explicit without deleting any source
    # bin based on solver type or sample count.
    grid.add(max_x + max(max(terminal_widths), max_x * 1.0e-3))
    result = np.asarray(sorted(grid), dtype=float)
    if np.any(np.diff(result) <= 0.0):
        raise FunctionEedfError("dimensionless EEDF grid is not increasing")
    return result


def _sample_bins(
    left: np.ndarray,
    right: np.ndarray,
    values: np.ndarray,
    points: np.ndarray,
) -> np.ndarray:
    indices = np.searchsorted(right, points, side="right")
    result = np.zeros_like(points)
    valid = indices < len(values)
    valid_indices = indices[valid]
    valid[valid] &= points[valid] >= left[valid_indices]
    result[valid] = values[indices[valid]]
    return result
