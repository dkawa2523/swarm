"""Numerical evidence produced by the stationary PN solve."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class SolverDiagnostics:
    iterations: int
    growth_frequency_s_inv: float
    full_pn_relative_residual: float
    shape_change: float
    eigenvalue_change: float
    mean_energy_eV: float
    eedf_tail_fraction: float
