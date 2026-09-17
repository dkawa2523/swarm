"""Typed numerical results produced by the native two-term solver."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from electron_swarm.core.results import RfPhaseResult
from electron_swarm.solvers.boltzmann_common.operators import KineticOperatorBlock


@dataclass(slots=True)
class NativeSolveDiagnostics:
    converged: bool
    iterations: int
    residual: float
    residual_requested_tolerance: float
    residual_effective_tolerance: float
    residual_roundoff_bound: float
    residual_backward_error: float
    residual_roundoff_limited: bool
    growth_frequency_s: float
    tail_probability: float
    edge_to_peak: float
    grid_max_eV: float
    regrid_cycles: int


@dataclass(slots=True)
class NativeDistributionResult:
    """Solved stationary distribution and its reusable case assembly."""

    energy_eV: np.ndarray
    edges_eV: np.ndarray
    widths_eV: np.ndarray
    eedf_eV_inv: np.ndarray
    diagnostics: NativeSolveDiagnostics
    metadata: dict[str, object]
    operator_block: KineticOperatorBlock


@dataclass(slots=True)
class TimePeriodicDistributionResult:
    """Cycle-averaged and phase-resolved homogeneous EEDF solution."""

    energy_eV: np.ndarray
    widths_eV: np.ndarray
    cycle_averaged_eedf_eV_inv: np.ndarray
    phase: RfPhaseResult
    diagnostics: dict[str, object]
    operator_block: KineticOperatorBlock
