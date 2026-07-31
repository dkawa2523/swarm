"""Numerical helpers for workflow table construction."""

from __future__ import annotations

import math
from typing import Any, Iterable

import numpy as np


def strictly_monotonic(values: Iterable[float]) -> bool:
    items = [float(value) for value in values]
    if len(items) < 2:
        return True
    diffs = [right - left for left, right in zip(items, items[1:])]
    return all(diff > 0.0 for diff in diffs) or all(diff < 0.0 for diff in diffs)


def finite_range(values: Iterable[object]) -> list[float] | None:
    finite = [
        number
        for number in (_float_or_none(value) for value in values)
        if number is not None
    ]
    if not finite:
        return None
    return [min(finite), max(finite)]


def positive_log_interpolate(
    rows: list[dict[str, Any]],
    scalar: str,
    x: float,
) -> float | None:
    points = [
        (float(row["E_over_N_Td"]), _finite_positive(row.get(scalar)))
        for row in rows
    ]
    points = [(point_x, value) for point_x, value in points if value is not None]
    if not points:
        return None
    points.sort()
    xs = [point[0] for point in points]
    logs = [math.log(point[1]) for point in points]
    if x < xs[0] or x > xs[-1]:
        return None
    if len(xs) == 1:
        return math.exp(logs[0]) if math.isclose(x, xs[0]) else None
    return math.exp(float(np.interp(x, xs, logs)))


def interpolate_no_extrapolate(
    x: float,
    known_x: list[float],
    known_delta: list[float],
) -> float:
    if not known_x:
        return 0.0
    paired = sorted(zip(known_x, known_delta))
    xs = [item[0] for item in paired]
    deltas = [item[1] for item in paired]
    if len(xs) == 1:
        return deltas[0] if math.isclose(x, xs[0]) else 0.0
    if x < xs[0] or x > xs[-1]:
        return 0.0
    return float(np.interp(x, xs, deltas))


def _finite_positive(value: object) -> float | None:
    number = _float_or_none(value)
    if number is None or number <= 0.0:
        return None
    return number


def _float_or_none(value: object) -> float | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None
