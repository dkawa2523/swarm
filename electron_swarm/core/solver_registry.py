"""Canonical solver descriptors for product solver modes."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Callable

from electron_swarm.core.capabilities import SolverCapabilities, SupportLevel
from electron_swarm.core.feature_resolution import (
    FeatureTreatment,
    monte_carlo_feature_resolver,
    multi_term_feature_resolver,
    propagator_feature_resolver,
    two_term_feature_resolver,
)
from electron_swarm.core.solver_ids import CANONICAL_SOLVER_IDS

if TYPE_CHECKING:  # pragma: no cover
    from electron_swarm.core.config import SwarmConfig
    from electron_swarm.core.cross_sections import (
        ActiveMixtureInputs,
        CrossSectionSet,
    )
    from electron_swarm.core.solver_configs import InternalSolverConfigs
    from electron_swarm.solvers.base import SwarmSolver

SolverFactory = Callable[
    ["SwarmConfig", "ActiveMixtureInputs", "InternalSolverConfigs"],
    "SwarmSolver",
]
MethodLabel = Callable[["SwarmConfig"], str]
FeatureResolver = Callable[
    ["SwarmConfig", "CrossSectionSet | None", SolverCapabilities],
    dict[str, FeatureTreatment],
]


@dataclass(frozen=True, slots=True)
class SolverDescriptor:
    id: str
    capabilities: SolverCapabilities
    method_label: MethodLabel
    factory: SolverFactory
    feature_resolver: FeatureResolver


def _two_term_factory(
    config: "SwarmConfig",
    cross_sections: "ActiveMixtureInputs",
    internal: "InternalSolverConfigs",
) -> "SwarmSolver":
    from electron_swarm.solvers.two_term import TwoTermSolver

    return TwoTermSolver(config, cross_sections, internal.two_term)


def _multi_term_factory(
    config: "SwarmConfig",
    cross_sections: "ActiveMixtureInputs",
    internal: "InternalSolverConfigs",
) -> "SwarmSolver":
    from electron_swarm.solvers.multi_term import MultiTermSolver

    return MultiTermSolver(
        config,
        cross_sections,
        internal.multi_term,
    )


def _monte_carlo_factory(
    config: "SwarmConfig",
    cross_sections: "ActiveMixtureInputs",
    internal: "InternalSolverConfigs",
) -> "SwarmSolver":
    from electron_swarm.solvers.monte_carlo import MonteCarloSolver

    return MonteCarloSolver(config, cross_sections, internal.monte_carlo)


def _propagator_factory(
    config: "SwarmConfig",
    cross_sections: "ActiveMixtureInputs",
    internal: "InternalSolverConfigs",
) -> "SwarmSolver":
    from electron_swarm.solvers.propagator import PropagatorSolver

    return PropagatorSolver(config, cross_sections, internal.propagator)


SOLVER_DESCRIPTORS: dict[str, SolverDescriptor] = {
    "two_term": SolverDescriptor(
        id="two_term",
        capabilities=SolverCapabilities(
            solver="two_term",
            angular_scattering=SupportLevel.EXACT,
            ionization_source=SupportLevel.EXACT,
            electron_electron=SupportLevel.APPROXIMATE,
            magnetic_field=SupportLevel.UNSUPPORTED,
            rf_field=SupportLevel.APPROXIMATE,
            tail_refinement=SupportLevel.APPROXIMATE,
        ),
        method_label=lambda config: str(config.solvers.two_term.backend),
        factory=_two_term_factory,
        feature_resolver=two_term_feature_resolver,
    ),
    "multi_term": SolverDescriptor(
        id="multi_term",
        capabilities=SolverCapabilities(
            solver="multi_term",
            angular_scattering=SupportLevel.EXACT,
            ionization_source=SupportLevel.UNSUPPORTED,
            electron_electron=SupportLevel.APPROXIMATE,
            magnetic_field=SupportLevel.UNSUPPORTED,
            rf_field=SupportLevel.UNSUPPORTED,
            tail_refinement=SupportLevel.APPROXIMATE,
        ),
        method_label=lambda config: str(config.solvers.multi_term.method),
        factory=_multi_term_factory,
        feature_resolver=multi_term_feature_resolver,
    ),
    "monte_carlo": SolverDescriptor(
        id="monte_carlo",
        capabilities=SolverCapabilities(
            solver="monte_carlo",
            angular_scattering=SupportLevel.EXACT,
            ionization_source=SupportLevel.UNSUPPORTED,
            electron_electron=SupportLevel.UNSUPPORTED,
            magnetic_field=SupportLevel.APPROXIMATE,
            rf_field=SupportLevel.UNSUPPORTED,
            tail_refinement=SupportLevel.APPROXIMATE,
        ),
        method_label=lambda config: "internal",
        factory=_monte_carlo_factory,
        feature_resolver=monte_carlo_feature_resolver,
    ),
    "propagator": SolverDescriptor(
        id="propagator",
        capabilities=SolverCapabilities(
            solver="propagator",
            angular_scattering=SupportLevel.EXACT,
            ionization_source=SupportLevel.EXACT,
            electron_electron=SupportLevel.UNSUPPORTED,
            magnetic_field=SupportLevel.UNSUPPORTED,
            rf_field=SupportLevel.UNSUPPORTED,
            tail_refinement=SupportLevel.APPROXIMATE,
        ),
        method_label=lambda config: str(config.solvers.propagator.method),
        factory=_propagator_factory,
        feature_resolver=propagator_feature_resolver,
    ),
}

if tuple(SOLVER_DESCRIPTORS) != CANONICAL_SOLVER_IDS:
    raise RuntimeError("solver descriptor registry does not match canonical solver ids")


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
    cross_sections: "ActiveMixtureInputs",
    internal: "InternalSolverConfigs",
    solver: str,
) -> "SwarmSolver":
    return solver_descriptor(solver).factory(config, cross_sections, internal)
