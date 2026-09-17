"""Prepared piecewise-linear cross-section tables for Monte Carlo kernels.

The product cross sections remain the source of truth.  This module only
projects their existing interpolation and endpoint policies onto one union
grid so every process shares a single interval lookup at run time.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np

from electron_swarm.core.cross_sections import CrossSectionProcess

try:  # The scalar fallback keeps the reference backend usable without Numba.
    from numba import njit

    NUMBA_AVAILABLE = True
except ImportError:  # pragma: no cover - exercised only in minimal installs
    NUMBA_AVAILABLE = False

    def njit(*_args, **_kwargs):  # type: ignore[no-untyped-def]
        def decorate(function):  # type: ignore[no-untyped-def]
            return function

        return decorate


@njit(cache=True, fastmath=False)
def evaluate_prepared_cross_sections(
    energy_eV: float,
    scale: float,
    grid_eV: np.ndarray,
    values_m2: np.ndarray,
    slopes_m2_eV: np.ndarray,
    process_min_eV: np.ndarray,
    process_max_eV: np.ndarray,
    right_values_m2: np.ndarray,
    out: np.ndarray,
) -> None:
    """Fill ``out`` with the prepared table value for one energy.

    Error extrapolation is checked by :class:`PreparedCrossSectionTable`
    before entering this numerical primitive.  Keeping that diagnostic in
    Python preserves the process name in the exception while the hot path is
    compiled.
    """

    energy = float(energy_eV)
    multiplier = float(scale)
    process_count = values_m2.shape[0]
    if not math.isfinite(energy) or not math.isfinite(multiplier):
        for process_index in range(process_count):
            out[process_index] = 0.0
        return

    grid_count = grid_eV.size
    if energy < grid_eV[0]:
        for process_index in range(process_count):
            out[process_index] = 0.0
        return

    if energy > grid_eV[grid_count - 1]:
        for process_index in range(process_count):
            out[process_index] = max(right_values_m2[process_index], 0.0) * multiplier
        return

    if energy == grid_eV[grid_count - 1]:
        interval = grid_count - 1
        for process_index in range(process_count):
            if energy < process_min_eV[process_index]:
                value = 0.0
            elif energy > process_max_eV[process_index]:
                value = right_values_m2[process_index]
            else:
                value = values_m2[process_index, interval]
            out[process_index] = max(value, 0.0) * multiplier
        return

    # Manual upper-bound search avoids allocating a temporary array and works
    # in both the compiled and reference implementations.
    lower = 0
    upper = grid_count
    while lower < upper:
        middle = (lower + upper) // 2
        if energy < grid_eV[middle]:
            upper = middle
        else:
            lower = middle + 1
    interval = max(0, min(lower - 1, grid_count - 2))
    delta = energy - grid_eV[interval]
    for process_index in range(process_count):
        if energy < process_min_eV[process_index]:
            value = 0.0
        elif energy > process_max_eV[process_index]:
            value = right_values_m2[process_index]
        else:
            value = (
                values_m2[process_index, interval]
                + slopes_m2_eV[process_index, interval] * delta
            )
        out[process_index] = max(value, 0.0) * multiplier


@dataclass(slots=True)
class PreparedCrossSectionTable:
    """Immutable interpolation data plus one caller-owned result buffer."""

    processes: tuple[CrossSectionProcess, ...]
    grid_eV: np.ndarray
    values_m2: np.ndarray
    slopes_m2_eV: np.ndarray
    process_min_eV: np.ndarray
    process_max_eV: np.ndarray
    right_values_m2: np.ndarray
    error_extrapolation: np.ndarray
    buffer: np.ndarray

    @classmethod
    def build(
        cls,
        processes: Sequence[CrossSectionProcess],
        *,
        multipliers: Sequence[float] | None = None,
    ) -> "PreparedCrossSectionTable":
        items = tuple(processes)
        if not items:
            raise ValueError(
                "prepared cross-section table requires at least one process"
            )
        factors = np.ones(len(items), dtype=float)
        if multipliers is not None:
            factors = np.asarray(tuple(multipliers), dtype=float)
            if factors.shape != (len(items),):
                raise ValueError("cross-section table multipliers must match processes")
            if np.any(~np.isfinite(factors)) or np.any(factors < 0.0):
                raise ValueError(
                    "cross-section table multipliers must be finite and nonnegative"
                )

        grid_parts = [item.energy_eV for item in items]
        grid_parts.extend(
            np.asarray([threshold], dtype=float)
            for item in items
            if (threshold := item.incident_threshold_eV) is not None
        )
        grid = np.unique(np.concatenate(grid_parts)).astype(float, copy=False)
        if grid.size < 2 or np.any(~np.isfinite(grid)) or np.any(np.diff(grid) <= 0.0):
            raise ValueError(
                "prepared cross-section energy grid must strictly increase"
            )
        policies = tuple(
            str(item.metadata.get("high_energy_extrapolation", "hold")).lower()
            for item in items
        )
        unsupported = sorted(set(policies) - {"zero", "hold", "error"})
        if unsupported:
            raise ValueError(
                "Unsupported high-energy extrapolation policy for prepared "
                f"cross sections: {unsupported}"
            )
        right = (
            np.asarray(
                [
                    0.0 if policy == "zero" else float(item.cross_section_m2[-1])
                    for item, policy in zip(items, policies, strict=True)
                ],
                dtype=float,
            )
            * factors
        )
        values = np.stack(
            [
                item.sigma(
                    grid,
                    left=0.0,
                    right=(
                        0.0 if policy == "zero" else float(item.cross_section_m2[-1])
                    ),
                )
                * factor
                for item, policy, factor in zip(
                    items,
                    policies,
                    factors,
                    strict=True,
                )
            ]
        )
        slopes = np.diff(values, axis=1) / np.diff(grid)[None, :]
        return cls(
            processes=items,
            grid_eV=np.ascontiguousarray(grid),
            values_m2=np.ascontiguousarray(values),
            slopes_m2_eV=np.ascontiguousarray(slopes),
            process_min_eV=np.asarray(
                [
                    max(
                        float(item.energy_eV[0]),
                        float(item.incident_threshold_eV or 0.0),
                    )
                    for item in items
                ],
                dtype=float,
            ),
            process_max_eV=np.asarray(
                [item.energy_eV[-1] for item in items], dtype=float
            ),
            right_values_m2=np.ascontiguousarray(right),
            error_extrapolation=np.asarray(
                [policy == "error" for policy in policies], dtype=bool
            ),
            buffer=np.zeros(len(items), dtype=float),
        )

    def evaluate(self, energy_eV: float, *, scale: float = 1.0) -> np.ndarray:
        energy = float(energy_eV)
        invalid = self.error_extrapolation & (energy > self.process_max_eV)
        if np.any(invalid):
            index = int(np.flatnonzero(invalid)[0])
            process = self.processes[index]
            raise ValueError(
                "Cross-section interpolation requested above the "
                f"tabulated range for {process.species}:{process.process}; "
                "increase the energy grid or set "
                "cross_sections.high_energy_extrapolation."
            )
        evaluate_prepared_cross_sections(
            energy,
            float(scale),
            self.grid_eV,
            self.values_m2,
            self.slopes_m2_eV,
            self.process_min_eV,
            self.process_max_eV,
            self.right_values_m2,
            self.buffer,
        )
        return self.buffer
