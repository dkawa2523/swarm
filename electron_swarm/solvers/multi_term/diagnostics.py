"""Compact diagnostics carried by multi-term solver results."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class SolverDiagnostics:
    warnings: tuple[str, ...]
    mean_energy_eV: float
    eedf_tail_fraction: float
    power_balance_relative_residual: float | None = None
