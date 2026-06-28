"""Canonical solver descriptors for product solver modes."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Callable

from electron_swarm.core.capabilities import SolverCapabilities, SupportLevel

if TYPE_CHECKING:  # pragma: no cover
    from electron_swarm.core.config import SwarmConfig
    from electron_swarm.core.cross_sections import CrossSectionSet
    from electron_swarm.core.solver_configs import InternalSolverConfigs
    from electron_swarm.solvers.base import SwarmSolver

SolverFactory = Callable[
    ["SwarmConfig", "CrossSectionSet", "InternalSolverConfigs"],
    "SwarmSolver",
]
MethodLabel = Callable[["SwarmConfig"], str]


@dataclass(frozen=True, slots=True)
class SolverDescriptor:
    id: str
    capabilities: SolverCapabilities
    method_label: MethodLabel
    factory: SolverFactory


def _two_term_factory(
    config: "SwarmConfig",
    cross_sections: "CrossSectionSet",
    internal: "InternalSolverConfigs",
) -> "SwarmSolver":
    from electron_swarm.solvers.two_term import TwoTermSolver

    return TwoTermSolver(config, cross_sections, internal.two_term)


def _multi_term_factory(
    config: "SwarmConfig",
    cross_sections: "CrossSectionSet",
    internal: "InternalSolverConfigs",
) -> "SwarmSolver":
    from electron_swarm.solvers.multi_term import MultiTermSolver

    return MultiTermSolver(
        config,
        cross_sections,
        internal.multi_term,
        internal.two_term,
    )


def _monte_carlo_factory(
    config: "SwarmConfig",
    cross_sections: "CrossSectionSet",
    internal: "InternalSolverConfigs",
) -> "SwarmSolver":
    from electron_swarm.solvers.internal_monte_carlo import MonteCarloSolver

    return MonteCarloSolver(config, cross_sections, internal.monte_carlo)


SOLVER_DESCRIPTORS: dict[str, SolverDescriptor] = {
    "two_term": SolverDescriptor(
        id="two_term",
        capabilities=SolverCapabilities(
            solver="two_term",
            angular_scattering=SupportLevel.APPROXIMATE,
            ionization_source=SupportLevel.EXACT,
            electron_electron=SupportLevel.APPROXIMATE,
            magnetic_field=SupportLevel.UNSUPPORTED,
            tail_refinement=SupportLevel.APPROXIMATE,
        ),
        method_label=lambda config: str(config.solvers.two_term.backend),
        factory=_two_term_factory,
    ),
    "multi_term": SolverDescriptor(
        id="multi_term",
        capabilities=SolverCapabilities(
            solver="multi_term",
            angular_scattering=SupportLevel.APPROXIMATE,
            ionization_source=SupportLevel.UNSUPPORTED,
            electron_electron=SupportLevel.APPROXIMATE,
            magnetic_field=SupportLevel.UNSUPPORTED,
            tail_refinement=SupportLevel.APPROXIMATE,
        ),
        method_label=lambda config: str(config.solvers.multi_term.method),
        factory=_multi_term_factory,
    ),
    "monte_carlo": SolverDescriptor(
        id="monte_carlo",
        capabilities=SolverCapabilities(
            solver="monte_carlo",
            angular_scattering=SupportLevel.APPROXIMATE,
            ionization_source=SupportLevel.UNSUPPORTED,
            electron_electron=SupportLevel.UNSUPPORTED,
            magnetic_field=SupportLevel.APPROXIMATE,
            tail_refinement=SupportLevel.APPROXIMATE,
        ),
        method_label=lambda config: "internal",
        factory=_monte_carlo_factory,
    ),
}

CANONICAL_SOLVER_IDS = tuple(SOLVER_DESCRIPTORS)


def solver_descriptor(solver: str) -> SolverDescriptor:
    try:
        return SOLVER_DESCRIPTORS[solver]
    except KeyError as exc:
        raise ValueError(f"Unknown canonical solver id: {solver!r}") from exc


def solver_method(config: "SwarmConfig", solver: str) -> str:
    """Return the user-facing method label for a canonical solver id."""

    return solver_descriptor(solver).method_label(config)


def build_solver(
    config: "SwarmConfig",
    cross_sections: "CrossSectionSet",
    internal: "InternalSolverConfigs",
    solver: str,
) -> "SwarmSolver":
    return solver_descriptor(solver).factory(config, cross_sections, internal)
