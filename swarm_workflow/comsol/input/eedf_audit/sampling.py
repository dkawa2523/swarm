"""Deterministic native-grid sampling plans for Function-EEDF audits."""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from ..function_eedf import ComsolFunctionEedfGrid
from .contracts import (
    _PROJECTION_SIGNIFICANCE_FRACTION,
    ComsolEedfAuditError,
)


def _audit_queries(
    grid: ComsolFunctionEedfGrid, *, moment_mean_count: int
) -> list[dict[str, Any]]:
    if moment_mean_count < 3:
        raise ComsolEedfAuditError("moment_mean_count must be at least three")
    energies = grid.electron_energies_eV
    means = grid.mean_energies_eV
    rows: list[dict[str, Any]] = []

    def append(
        kind: str,
        group: str,
        position: str,
        energy: float,
        mean: float,
        expected_a: float,
        expected_b: float | None = None,
    ) -> None:
        rows.append(
            {
                "point_id": f"q{len(rows):07d}",
                "kind": kind,
                "group": group,
                "position": position,
                "electron_energy_eV": energy,
                "mean_energy_eV": mean,
                "expected_a": expected_a,
                "expected_b": expected_a if expected_b is None else expected_b,
            }
        )

    anchor_count = (moment_mean_count + 1) // 2
    midpoint_count = moment_mean_count - anchor_count
    moment_specs: list[tuple[str, float, np.ndarray]] = []
    for index in _representative_indices(len(means), anchor_count):
        moment_specs.append(
            (f"moment_anchor_{index:04d}", means[index], grid.values_eV_m32[index])
        )
    for index in _representative_indices(len(means) - 1, midpoint_count):
        mean = 0.5 * (means[index] + means[index + 1])
        values = 0.5 * (grid.values_eV_m32[index] + grid.values_eV_m32[index + 1])
        moment_specs.append((f"moment_mid_{index:04d}", mean, values))
    for group, mean, values in moment_specs:
        for energy_index, (energy, value) in enumerate(
            zip(energies, values, strict=True)
        ):
            append("moment", group, str(energy_index), energy, mean, value)

    for mean_index, mean in enumerate(means):
        for energy_index in _representative_significant_indices(
            grid.values_eV_m32[mean_index], 13
        ):
            append(
                "anchor",
                f"anchor_{mean_index:04d}",
                str(energy_index),
                energies[energy_index],
                mean,
                grid.values_eV_m32[mean_index, energy_index],
            )

    mean_indices = _representative_indices(len(means), 7)
    for mean_index in mean_indices:
        midpoint_values = 0.5 * (
            grid.values_eV_m32[mean_index, :-1] + grid.values_eV_m32[mean_index, 1:]
        )
        for energy_index in _representative_significant_indices(midpoint_values, 13):
            append(
                "energy_midpoint",
                f"energy_mid_{mean_index:04d}",
                str(energy_index),
                0.5 * (energies[energy_index] + energies[energy_index + 1]),
                means[mean_index],
                0.5
                * (
                    grid.values_eV_m32[mean_index, energy_index]
                    + grid.values_eV_m32[mean_index, energy_index + 1]
                ),
            )

    for mean_index in _representative_indices(len(means) - 1, 7):
        cell_values = 0.25 * (
            grid.values_eV_m32[mean_index, :-1]
            + grid.values_eV_m32[mean_index, 1:]
            + grid.values_eV_m32[mean_index + 1, :-1]
            + grid.values_eV_m32[mean_index + 1, 1:]
        )
        for energy_index in _representative_significant_indices(cell_values, 13):
            append(
                "cell_center",
                f"cell_{mean_index:04d}_{energy_index:04d}",
                "center",
                0.5 * (energies[energy_index] + energies[energy_index + 1]),
                0.5 * (means[mean_index] + means[mean_index + 1]),
                float(cell_values[energy_index]),
            )

    for mean_index in _representative_indices(len(means) - 2, 7) + 1:
        left_gap = means[mean_index] - means[mean_index - 1]
        right_gap = means[mean_index + 1] - means[mean_index]
        delta = 1.0e-4 * min(left_gap, right_gap)
        left_slopes = (
            grid.values_eV_m32[mean_index] - grid.values_eV_m32[mean_index - 1]
        ) / left_gap
        right_slopes = (
            grid.values_eV_m32[mean_index + 1] - grid.values_eV_m32[mean_index]
        ) / right_gap
        slope_activity = np.maximum(np.abs(left_slopes), np.abs(right_slopes))
        for energy_index in _representative_significant_indices(slope_activity, 13):
            group = f"derivative_{mean_index:04d}_{energy_index:04d}"
            center_value = grid.values_eV_m32[mean_index, energy_index]
            left_slope = left_slopes[energy_index]
            right_slope = right_slopes[energy_index]
            append(
                "derivative",
                group,
                "left",
                energies[energy_index],
                means[mean_index] - delta,
                center_value - delta * left_slope,
            )
            append(
                "derivative",
                group,
                "center",
                energies[energy_index],
                means[mean_index],
                center_value,
            )
            append(
                "derivative",
                group,
                "right",
                energies[energy_index],
                means[mean_index] + delta,
                center_value + delta * right_slope,
            )

    mid_mean_index = len(means) // 2
    active_mid_indices = _representative_significant_indices(
        grid.values_eV_m32[mid_mean_index], 3
    )
    mid_energy_index = int(active_mid_indices[len(active_mid_indices) // 2])
    outside_energies = (-0.05 * max(energies[-1], 1.0), 1.05 * energies[-1])
    outside_means = (0.5 * means[0], 1.05 * means[-1])
    for side, energy, source_index in (
        ("energy_low", outside_energies[0], 0),
        ("energy_high", outside_energies[1], len(energies) - 1),
    ):
        append(
            "outside",
            side,
            "constant",
            energy,
            means[mid_mean_index],
            grid.values_eV_m32[mid_mean_index, source_index],
        )
    for side, mean, source_index in (
        ("mean_low", outside_means[0], 0),
        ("mean_high", outside_means[1], len(means) - 1),
    ):
        append(
            "outside",
            side,
            "constant",
            energies[mid_energy_index],
            mean,
            grid.values_eV_m32[source_index, mid_energy_index],
        )
    return rows


def _representative_indices(size: int, count: int) -> np.ndarray:
    if size <= 0 or count <= 0:
        return np.asarray([], dtype=int)
    return np.unique(np.rint(np.linspace(0, size - 1, min(size, count))).astype(int))


def _representative_significant_indices(
    values: np.ndarray,
    count: int,
) -> np.ndarray:
    """Sample the active support without repeatedly landing on exact zeros."""

    magnitude = np.abs(np.asarray(values, dtype=float))
    if magnitude.ndim != 1 or magnitude.size == 0 or count <= 0:
        return np.asarray([], dtype=int)
    peak = float(np.max(magnitude))
    if not math.isfinite(peak) or peak <= 0.0:
        return _representative_indices(magnitude.size, count)
    candidates = np.flatnonzero(magnitude >= peak * _PROJECTION_SIGNIFICANCE_FRACTION)
    positions = _representative_indices(len(candidates), count)
    return candidates[positions]
