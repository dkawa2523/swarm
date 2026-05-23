"""Solver capability matrix for product schema v2."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class SupportLevel(str, Enum):
    EXACT = "exact"
    APPROXIMATE = "approximate"
    POSTPROCESS = "postprocess"
    SURROGATE = "surrogate"
    DIAGNOSTIC = "diagnostic"
    UNSUPPORTED = "unsupported"


@dataclass(frozen=True, slots=True)
class SolverCapabilities:
    solver: str
    electron_neutral: SupportLevel
    angular_scattering: SupportLevel
    electron_electron: SupportLevel
    magnetic_field: SupportLevel
    tail_refinement: SupportLevel
    bulk_transport: SupportLevel


CAPABILITIES_BY_SOLVER: dict[str, SolverCapabilities] = {
    "two_term": SolverCapabilities(
        solver="two_term",
        electron_neutral=SupportLevel.APPROXIMATE,
        angular_scattering=SupportLevel.APPROXIMATE,
        electron_electron=SupportLevel.POSTPROCESS,
        magnetic_field=SupportLevel.UNSUPPORTED,
        tail_refinement=SupportLevel.APPROXIMATE,
        bulk_transport=SupportLevel.APPROXIMATE,
    ),
    "multi_term": SolverCapabilities(
        solver="multi_term",
        electron_neutral=SupportLevel.SURROGATE,
        angular_scattering=SupportLevel.SURROGATE,
        electron_electron=SupportLevel.POSTPROCESS,
        magnetic_field=SupportLevel.UNSUPPORTED,
        tail_refinement=SupportLevel.APPROXIMATE,
        bulk_transport=SupportLevel.UNSUPPORTED,
    ),
    "monte_carlo": SolverCapabilities(
        solver="monte_carlo",
        electron_neutral=SupportLevel.APPROXIMATE,
        angular_scattering=SupportLevel.APPROXIMATE,
        electron_electron=SupportLevel.UNSUPPORTED,
        magnetic_field=SupportLevel.UNSUPPORTED,
        tail_refinement=SupportLevel.APPROXIMATE,
        bulk_transport=SupportLevel.APPROXIMATE,
    ),
}


def get_solver_capabilities(solver: str) -> SolverCapabilities:
    try:
        return CAPABILITIES_BY_SOLVER[solver]
    except KeyError as exc:
        raise ValueError(f"Unknown canonical solver id: {solver!r}") from exc
