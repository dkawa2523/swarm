"""Solver capability matrix for product schema v2."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class SupportLevel(str, Enum):
    EXACT = "exact"
    APPROXIMATE = "approximate"
    UNSUPPORTED = "unsupported"


@dataclass(frozen=True, slots=True)
class SolverCapabilities:
    solver: str
    angular_scattering: SupportLevel
    ionization_source: SupportLevel
    electron_electron: SupportLevel
    magnetic_field: SupportLevel
    tail_refinement: SupportLevel


def get_solver_capabilities(solver: str) -> SolverCapabilities:
    from electron_swarm.core.solver_registry import solver_descriptor

    return solver_descriptor(solver).capabilities
