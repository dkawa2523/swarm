"""Diagnostics for multi-term solver quality and validity."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, slots=True)
class AngularConvergence:
    current_lmax: int
    previous_lmax: int | None
    relative_changes: dict[str, float]
    converged: bool


@dataclass(frozen=True, slots=True)
class SolverDiagnostics:
    warnings: tuple[str, ...]
    mean_energy_eV: float
    eedf_tail_fraction: float
    power_balance_relative_residual: float | None = None
    angular_convergence: AngularConvergence | None = None
    fit_residual: float | None = None
    n_samples: int | None = None


def compute_lmax_convergence(
    current_lmax: int,
    current: dict[str, float],
    previous_lmax: int | None,
    previous: dict[str, float] | None,
    tolerance: float = 0.02,
) -> AngularConvergence:
    if previous_lmax is None or previous is None:
        return AngularConvergence(current_lmax, None, {}, False)
    rel: dict[str, float] = {}
    for key, val in current.items():
        if key in previous:
            denom = max(abs(float(val)), 1.0e-300)
            rel[key] = abs(float(val) - float(previous[key])) / denom
    converged = bool(rel) and all(v <= tolerance for v in rel.values() if np.isfinite(v))
    return AngularConvergence(current_lmax, previous_lmax, rel, converged)
