"""Moment integration, interpolation, and conservative EEDF tilting."""

from __future__ import annotations

import math

import numpy as np
from scipy.interpolate import PchipInterpolator

from .contracts import C1FunctionEedf, ComsolFunctionEedfGrid, FunctionEedfError


_GAUSS_NODES, _GAUSS_WEIGHTS = np.polynomial.legendre.leggauss(8)


def evaluate_comsol_function_eedf_grid(
    grid: ComsolFunctionEedfGrid,
    electron_energy_eV: np.ndarray | float,
    mean_energy_eV: float,
) -> np.ndarray:
    """Evaluate COMSOL's linear 2D table with constant extrapolation."""

    mean = float(mean_energy_eV)
    energy = np.asarray(electron_energy_eV, dtype=float)
    if not math.isfinite(mean) or mean <= 0.0 or not np.all(np.isfinite(energy)):
        raise ValueError(
            "mean energy must be positive and all Function-EEDF arguments finite"
        )
    anchors = grid.mean_energies_eV
    if mean <= anchors[0]:
        row = grid.values_eV_m32[0]
    elif mean >= anchors[-1]:
        row = grid.values_eV_m32[-1]
    else:
        left = int(np.searchsorted(anchors, mean, side="right") - 1)
        fraction = (mean - anchors[left]) / (anchors[left + 1] - anchors[left])
        row = (1.0 - fraction) * grid.values_eV_m32[
            left
        ] + fraction * grid.values_eV_m32[left + 1]
    flat = np.interp(
        energy.reshape(-1),
        grid.electron_energies_eV,
        row,
        left=float(row[0]),
        right=float(row[-1]),
    )
    return flat.reshape(energy.shape)


def evaluate_c1_function_eedf(
    representation: C1FunctionEedf,
    electron_energy_eV: np.ndarray | float,
    mean_energy_eV: float,
) -> np.ndarray:
    """Evaluate ``f0(E,m)`` with local smoothstep blending."""

    mean = float(mean_energy_eV)
    if not math.isfinite(mean) or mean <= 0.0:
        raise ValueError("mean_energy_eV must be positive and finite")
    energy = np.asarray(electron_energy_eV, dtype=float)
    x_value = energy / mean
    anchors = representation.mean_energies_eV
    if mean <= anchors[0]:
        shape = _evaluate_shape(representation, 0, x_value)
    elif mean >= anchors[-1]:
        shape = _evaluate_shape(representation, len(anchors) - 1, x_value)
    else:
        left_index = int(np.searchsorted(anchors, mean, side="right") - 1)
        right_index = left_index + 1
        fraction = (mean - anchors[left_index]) / (
            anchors[right_index] - anchors[left_index]
        )
        weight = fraction * fraction * (3.0 - 2.0 * fraction)
        shape = (1.0 - weight) * _evaluate_shape(
            representation, left_index, x_value
        ) + weight * _evaluate_shape(representation, right_index, x_value)
    return np.asarray(shape, dtype=float) / (mean**1.5)


def scaled_shape_moments(
    representation: C1FunctionEedf,
    mean_energy_eV: float,
) -> tuple[float, float]:
    """Integrate a blended shape exactly as the closure evaluates it."""

    x_values = representation.dimensionless_energy
    left = x_values[:-1, None]
    right = x_values[1:, None]
    half_width = 0.5 * (right - left)
    points = 0.5 * (right + left) + half_width * _GAUSS_NODES
    physical_f0 = evaluate_c1_function_eedf(
        representation, float(mean_energy_eV) * points, mean_energy_eV
    )
    shape = (float(mean_energy_eV) ** 1.5) * physical_f0
    weighted = half_width * _GAUSS_WEIGHTS
    return (
        float(np.sum(weighted * np.sqrt(points) * shape)),
        float(np.sum(weighted * np.power(points, 1.5) * shape)),
    )


def piecewise_linear_weighted_moments(
    electron_energies_eV: np.ndarray,
    values_eV_m32: np.ndarray,
) -> tuple[float, float]:
    """Integrate the two EEDF moments exactly for piecewise-linear ``f0``."""

    energies = np.asarray(electron_energies_eV, dtype=float)
    values = np.asarray(values_eV_m32, dtype=float)
    if (
        energies.ndim != 1
        or values.shape != energies.shape
        or len(energies) < 2
        or energies[0] != 0.0
        or np.any(np.diff(energies) <= 0.0)
        or np.any(~np.isfinite(values))
    ):
        raise FunctionEedfError("piecewise-linear Function-EEDF row is invalid")
    return (
        _piecewise_linear_weighted_integral(energies, values, 0.5),
        _piecewise_linear_weighted_integral(energies, values, 1.5),
    )


def pchip_weighted_moments(
    x_values: np.ndarray,
    y_values: np.ndarray,
) -> tuple[float, float]:
    interpolator = _pchip(x_values, y_values)
    left = x_values[:-1, None]
    right = x_values[1:, None]
    half_width = 0.5 * (right - left)
    points = 0.5 * (right + left) + half_width * _GAUSS_NODES
    values = interpolator(points)
    if np.any(~np.isfinite(values)):
        raise FunctionEedfError("cubic Function EEDF is nonfinite")
    weighted = half_width * _GAUSS_WEIGHTS
    return (
        float(np.sum(weighted * np.sqrt(points) * values)),
        float(np.sum(weighted * np.power(points, 1.5) * values)),
    )


def _tilt_to_unit_mean(x_values: np.ndarray, values: np.ndarray) -> np.ndarray:
    positive = values > 0.0
    if not np.any(positive) or x_values[positive][-1] < 1.0:
        raise FunctionEedfError("cannot mean-match the Function EEDF")

    def tilted(lam: float) -> tuple[np.ndarray, float]:
        exponent = lam * (x_values - 1.0)
        exponent -= float(np.max(exponent[positive]))
        candidate = values * np.exp(np.clip(exponent, -745.0, 0.0))
        normalization, first_moment = pchip_weighted_moments(x_values, candidate)
        if not math.isfinite(normalization) or normalization <= 0.0:
            raise FunctionEedfError("mean matching lost normalization")
        candidate /= normalization
        return candidate, first_moment / normalization

    result, mean = tilted(0.0)
    if abs(mean - 1.0) <= 1.0e-12:
        return result
    low, high = -1.0, 1.0
    _unused, low_mean = tilted(low)
    _unused, high_mean = tilted(high)
    for _index in range(80):
        if low_mean <= 1.0 <= high_mean:
            break
        if low_mean > 1.0:
            low *= 2.0
            _unused, low_mean = tilted(low)
        if high_mean < 1.0:
            high *= 2.0
            _unused, high_mean = tilted(high)
    else:
        raise FunctionEedfError("could not bracket EEDF mean matching")
    for _index in range(80):
        middle = 0.5 * (low + high)
        result, mean = tilted(middle)
        if abs(mean - 1.0) <= 1.0e-12:
            break
        if mean < 1.0:
            low = middle
        else:
            high = middle
    return result


def _tilt_linear_row_to_mean(
    energies: np.ndarray,
    values: np.ndarray,
    requested_mean: float,
) -> np.ndarray:
    positive = values > 0.0
    if not np.any(positive):
        raise FunctionEedfError("cannot project an empty Function-EEDF row")

    def tilted(lam: float) -> tuple[np.ndarray, float]:
        exponent = lam * (energies / requested_mean - 1.0)
        exponent -= float(np.max(exponent[positive]))
        candidate = values * np.exp(np.clip(exponent, -745.0, 0.0))
        normalization, first_moment = piecewise_linear_weighted_moments(
            energies,
            candidate,
        )
        if not math.isfinite(normalization) or normalization <= 0.0:
            raise FunctionEedfError("COMSOL grid projection lost normalization")
        candidate /= normalization
        return candidate, first_moment / normalization

    result, mean = tilted(0.0)
    if abs(mean - requested_mean) <= 1.0e-13 * requested_mean:
        return result
    low, high = -1.0, 1.0
    _unused, low_mean = tilted(low)
    _unused, high_mean = tilted(high)
    for _index in range(80):
        if low_mean <= requested_mean <= high_mean:
            break
        if low_mean > requested_mean:
            low *= 2.0
            _unused, low_mean = tilted(low)
        if high_mean < requested_mean:
            high *= 2.0
            _unused, high_mean = tilted(high)
    else:
        raise FunctionEedfError("could not bracket COMSOL EEDF grid projection")
    for _index in range(100):
        middle = 0.5 * (low + high)
        result, mean = tilted(middle)
        if abs(mean - requested_mean) <= 1.0e-13 * requested_mean:
            break
        if mean < requested_mean:
            low = middle
        else:
            high = middle
    normalization, _unused = piecewise_linear_weighted_moments(energies, result)
    return result / normalization


def _piecewise_linear_weighted_integral(
    energies: np.ndarray,
    values: np.ndarray,
    power: float,
) -> float:
    left = energies[:-1]
    right = energies[1:]
    width = right - left
    left_weight = np.empty_like(width)
    right_weight = np.empty_like(width)

    at_origin = left == 0.0
    if np.any(at_origin):
        scale = np.power(width[at_origin], power + 1.0)
        left_weight[at_origin] = scale / ((power + 1.0) * (power + 2.0))
        right_weight[at_origin] = scale / (power + 2.0)

    narrow = (~at_origin) & (width / np.maximum(left, 1.0e-300) < 0.25)
    if np.any(narrow):
        # Vectorized binomial expansion avoids cancellation for dense adaptive
        # knots while keeping the moment projection cheap inside root solves.
        z = width[narrow] / left[narrow]
        coefficient = 1.0
        z_power = np.ones_like(z)
        a0 = np.zeros_like(z)
        a1 = np.zeros_like(z)
        for order in range(40):
            term = coefficient * z_power
            a0 += term / (order + 1.0)
            a1 += term / (order + 2.0)
            coefficient *= (power - order) / (order + 1.0)
            z_power *= z
        scale = width[narrow] * np.power(left[narrow], power)
        left_weight[narrow] = scale * (a0 - a1)
        right_weight[narrow] = scale * a1

    regular = ~(at_origin | narrow)
    if np.any(regular):
        a = left[regular]
        b = right[regular]
        w = width[regular]
        base = (np.power(b, power + 1.0) - np.power(a, power + 1.0)) / (power + 1.0)
        raised = (np.power(b, power + 2.0) - np.power(a, power + 2.0)) / (power + 2.0)
        right_weight[regular] = (raised - a * base) / w
        left_weight[regular] = base - right_weight[regular]
    return float(np.dot(left_weight, values[:-1]) + np.dot(right_weight, values[1:]))


def _evaluate_shape(
    representation: C1FunctionEedf,
    index: int,
    x_value: np.ndarray,
) -> np.ndarray:
    interpolator = _pchip(
        representation.dimensionless_energy,
        representation.shapes[index],
    )
    maximum = representation.dimensionless_energy[-1]
    result = np.asarray(interpolator(np.clip(x_value, 0.0, maximum)), dtype=float)
    return np.where(x_value > maximum, 0.0, result)


def _pchip_minimum(x_values: np.ndarray, values: np.ndarray) -> float:
    interpolator = _pchip(x_values, values)
    left = x_values[:-1, None]
    right = x_values[1:, None]
    points = 0.5 * (right + left) + 0.5 * (right - left) * _GAUSS_NODES
    return min(float(np.min(values)), float(np.min(interpolator(points))))


def _pchip(x_values: np.ndarray, values: np.ndarray) -> PchipInterpolator:
    # Compact EEDF tails legitimately span hundreds of orders of magnitude.
    # SciPy's harmonic-slope construction can overflow in an intermediate
    # reciprocal even when its monotone result is finite. The positivity and
    # finiteness audits below remain authoritative.
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        return PchipInterpolator(x_values, values, extrapolate=False)
