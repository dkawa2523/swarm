"""Adaptive tensor-grid construction for a continuous Function-EEDF family."""

from __future__ import annotations

import math

import numpy as np

from .contracts import C1FunctionEedf, CollisionRateKernel, FunctionEedfError
from .kernels import collision_rate_coefficient
from .moments import (
    _tilt_linear_row_to_mean,
    evaluate_c1_function_eedf,
)


_MEAN_AXIS_AUDIT_FRACTIONS = (0.25, 0.5, 0.75)


def mean_axis_refinement(
    representation: C1FunctionEedf,
    energies: np.ndarray,
    audit_energy: np.ndarray,
    means: np.ndarray,
    rows: np.ndarray,
    kernels: tuple[CollisionRateKernel, ...],
    rate_scales: tuple[float, ...],
    reference_cache: dict[float, np.ndarray],
    *,
    shape_tolerance: float,
    rate_tolerance: float,
) -> tuple[list[float], float, float | None]:
    """Find C1 mean rows needed to bound the actual tensor-table error."""

    def reference(mean: float) -> np.ndarray:
        if mean not in reference_cache:
            reference_cache[mean] = np.maximum(
                evaluate_c1_function_eedf(
                    representation,
                    audit_energy,
                    mean,
                ),
                0.0,
            )
        return reference_cache[mean]

    additions: list[float] = []
    worst_shape = 0.0
    worst_rate: float | None = 0.0 if kernels else None
    for index, (left, right) in enumerate(zip(means[:-1], means[1:], strict=True)):
        left_row = rows[index]
        right_row = rows[index + 1]
        interval_score = 0.0
        interval_point: float | None = None
        for fraction in _MEAN_AXIS_AUDIT_FRACTIONS:
            mean = float(left + fraction * (right - left))
            table_row = (1.0 - fraction) * left_row + fraction * right_row
            candidate = np.interp(audit_energy, energies, table_row)
            shape_error, rate_error = row_projection_error(
                audit_energy,
                reference(mean),
                candidate,
                kernels,
                rate_scales,
            )
            worst_shape = max(worst_shape, shape_error)
            if worst_rate is not None and rate_error is not None:
                worst_rate = max(worst_rate, rate_error)
            score = max(
                shape_error / shape_tolerance,
                rate_error / rate_tolerance if rate_error is not None else 0.0,
            )
            if score > interval_score:
                interval_score = score
                interval_point = mean
        if interval_score > 1.0 and interval_point is not None:
            additions.append(interval_point)
    return additions, worst_shape, worst_rate


def initial_energy_axis(
    minimum_energy: float,
    maximum_energy: float,
    kernels: tuple[CollisionRateKernel, ...],
) -> np.ndarray:
    points = {0.0, maximum_energy}
    points.update(np.geomspace(minimum_energy, maximum_energy, 65).tolist())
    for kernel in kernels:
        points.add(float(np.clip(kernel.threshold_eV, 0.0, maximum_energy)))
        points.update(
            float(value)
            for value in kernel.electron_energies_eV
            if 0.0 < value < maximum_energy
        )
    return np.asarray(sorted(points), dtype=float)


def audit_energy_axis(
    minimum_energy: float,
    maximum_energy: float,
    kernels: tuple[CollisionRateKernel, ...],
) -> np.ndarray:
    points = {0.0, maximum_energy}
    points.update(np.geomspace(minimum_energy, maximum_energy, 2049).tolist())
    for kernel in kernels:
        points.add(float(np.clip(kernel.threshold_eV, 0.0, maximum_energy)))
        points.update(
            float(value)
            for value in kernel.electron_energies_eV
            if 0.0 < value < maximum_energy
        )
    return np.asarray(sorted(points), dtype=float)


def project_rows(
    representation: C1FunctionEedf,
    energies: np.ndarray,
    means: np.ndarray,
) -> np.ndarray:
    rows = []
    for mean in means:
        sampled = np.maximum(
            evaluate_c1_function_eedf(representation, energies, float(mean)),
            0.0,
        )
        rows.append(_tilt_linear_row_to_mean(energies, sampled, float(mean)))
    return np.asarray(rows, dtype=float)


def adapt_energy_axis(
    representation: C1FunctionEedf,
    initial: np.ndarray,
    audit_energy: np.ndarray,
    means: np.ndarray,
    kernels: tuple[CollisionRateKernel, ...],
    rate_scales: tuple[float, ...],
    *,
    shape_tolerance: float,
    rate_tolerance: float,
) -> np.ndarray:
    energies = initial
    references = [
        np.maximum(
            evaluate_c1_function_eedf(representation, audit_energy, float(mean)),
            0.0,
        )
        for mean in means
    ]
    worst_shape = math.inf
    worst_rate: float | None = math.inf if kernels else None
    for _iteration in range(32):
        rows = project_rows(representation, energies, means)
        global_score = np.zeros_like(audit_energy)
        needs_refinement = False
        worst_shape = 0.0
        worst_rate = 0.0 if kernels else None
        for _mean, row, reference in zip(means, rows, references, strict=True):
            candidate = np.interp(audit_energy, energies, row)
            shape_error, rate_error = row_projection_error(
                audit_energy,
                reference,
                candidate,
                kernels,
                rate_scales,
            )
            worst_shape = max(worst_shape, shape_error)
            if worst_rate is not None and rate_error is not None:
                worst_rate = max(worst_rate, rate_error)
            if shape_error <= shape_tolerance and (
                rate_error is None or rate_error <= rate_tolerance
            ):
                continue
            needs_refinement = True
            score = np.sqrt(audit_energy) * np.abs(candidate - reference)
            for kernel, rate_scale in zip(kernels, rate_scales, strict=True):
                reference_rate = collision_rate_coefficient(
                    audit_energy, reference, kernel
                )
                if reference_rate <= 0.0 and rate_scale <= 0.0:
                    continue
                right = (
                    0.0
                    if kernel.high_energy_extrapolation == "zero"
                    else float(kernel.cross_sections_m2[-1])
                )
                sigma = np.interp(
                    audit_energy,
                    kernel.electron_energies_eV,
                    kernel.cross_sections_m2,
                    left=0.0,
                    right=right,
                )
                score = np.maximum(
                    score,
                    audit_energy
                    * sigma
                    * np.abs(candidate - reference)
                    / max(reference_rate, rate_scale),
                )
            global_score = np.maximum(global_score, score)
        if not needs_refinement:
            return energies
        additions = interval_refinement_points(
            energies,
            audit_energy,
            global_score,
            maximum_points=64,
        )
        if not additions:
            break
        energies = np.unique(np.concatenate((energies, np.asarray(sorted(additions)))))
        if len(energies) > 4097:
            raise FunctionEedfError(
                "adaptive Function-EEDF energy grid exceeded 4097 points"
            )
    raise FunctionEedfError(
        "adaptive Function-EEDF energy refinement did not converge "
        f"({len(energies)} points, TV={worst_shape:.3e}, "
        f"rate={worst_rate if worst_rate is not None else 'not_requested'})"
    )


def interval_refinement_points(
    energies: np.ndarray,
    audit_energy: np.ndarray,
    score: np.ndarray,
    *,
    maximum_points: int,
) -> set[float]:
    """Choose at most one highest-error audit point from each table interval."""

    interval = np.searchsorted(energies, audit_energy, side="right") - 1
    available = (
        (interval >= 0)
        & (interval < len(energies) - 1)
        & ~np.isin(audit_energy, energies)
        & (score > 0.0)
    )
    if not np.any(available):
        return set()
    interval_score = np.zeros(len(energies) - 1, dtype=float)
    np.maximum.at(interval_score, interval[available], score[available])
    selected_intervals = np.argsort(interval_score)[-maximum_points:]
    additions: set[float] = set()
    for selected in selected_intervals:
        if interval_score[selected] <= 0.0:
            continue
        local = available & (interval == selected)
        index = int(np.argmax(np.where(local, score, -1.0)))
        additions.add(float(audit_energy[index]))
    return additions


def row_projection_error(
    energies: np.ndarray,
    reference: np.ndarray,
    candidate: np.ndarray,
    kernels: tuple[CollisionRateKernel, ...],
    rate_scales: tuple[float, ...],
) -> tuple[float, float | None]:
    shape_error = 0.5 * float(
        np.trapezoid(np.sqrt(energies) * np.abs(candidate - reference), energies)
    )
    if not kernels:
        return shape_error, None
    rate_error = 0.0
    for kernel, rate_scale in zip(kernels, rate_scales, strict=True):
        expected = collision_rate_coefficient(energies, reference, kernel)
        actual = collision_rate_coefficient(energies, candidate, kernel)
        scale = max(expected, rate_scale)
        if scale > 0.0:
            rate_error = max(rate_error, abs(actual - expected) / scale)
        elif actual != 0.0:
            rate_error = math.inf
    return shape_error, rate_error


def table_projection_error(
    representation: C1FunctionEedf,
    energies: np.ndarray,
    means: np.ndarray,
    rows: np.ndarray,
    audit_energy: np.ndarray,
    kernels: tuple[CollisionRateKernel, ...],
    rate_scales: tuple[float, ...],
) -> tuple[float, float | None]:
    shape_error = 0.0
    rate_error = 0.0 if kernels else None
    for mean, row in zip(means, rows, strict=True):
        reference = np.maximum(
            evaluate_c1_function_eedf(representation, audit_energy, float(mean)),
            0.0,
        )
        candidate = np.interp(audit_energy, energies, row)
        local_shape, local_rate = row_projection_error(
            audit_energy, reference, candidate, kernels, rate_scales
        )
        shape_error = max(shape_error, local_shape)
        if rate_error is not None and local_rate is not None:
            rate_error = max(rate_error, local_rate)
    return shape_error, rate_error


def rate_error_scales(
    representation: C1FunctionEedf,
    energies: np.ndarray,
    means: np.ndarray,
    kernels: tuple[CollisionRateKernel, ...],
    *,
    importance_fraction: float,
) -> tuple[float, ...]:
    if not kernels:
        return ()
    peaks = [0.0] * len(kernels)
    for mean in means:
        reference = np.maximum(
            evaluate_c1_function_eedf(representation, energies, float(mean)),
            0.0,
        )
        for index, kernel in enumerate(kernels):
            peaks[index] = max(
                peaks[index],
                collision_rate_coefficient(energies, reference, kernel),
            )
    return tuple(peak * importance_fraction for peak in peaks)


def closure_energy_support_max(representation: C1FunctionEedf) -> float:
    """Return the physical support of every local C1 mean-energy blend."""

    support_x: list[float] = []
    x_values = representation.dimensionless_energy
    for shape in representation.shapes:
        positive = np.flatnonzero(shape > 0.0)
        if not len(positive):
            raise FunctionEedfError("Function-EEDF shape has empty support")
        zero_index = min(int(positive[-1]) + 1, len(x_values) - 1)
        support_x.append(float(x_values[zero_index]))
    supports = np.asarray(support_x, dtype=float)
    anchors = representation.mean_energies_eV
    candidates = (anchors * supports).tolist()
    candidates.extend(
        float(right_mean * max(left_support, right_support))
        for right_mean, left_support, right_support in zip(
            anchors[1:],
            supports[:-1],
            supports[1:],
            strict=True,
        )
    )
    maximum = max(candidates)
    if not math.isfinite(maximum) or maximum <= 0.0:
        raise FunctionEedfError("Function-EEDF support is invalid")
    return maximum
