"""Numerical helpers for workflow table construction."""

from __future__ import annotations

import math
from typing import Iterable


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


def _float_or_none(value: object) -> float | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None
