"""Shared solver identifiers and presentation metadata."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class SolverInfo:
    name: str
    legacy_tag: str
    plot_label: str
    plot_color: str


SOLVER_REGISTRY: dict[str, SolverInfo] = {
    "monte_carlo": SolverInfo(
        name="monte_carlo",
        legacy_tag="mc",
        plot_label="Monte Carlo",
        plot_color="#c62828",
    ),
    "boltzmann_two_term": SolverInfo(
        name="boltzmann_two_term",
        legacy_tag="boltzmann",
        plot_label="Boltzmann",
        plot_color="#1565c0",
    ),
    "multiterm_boltzmann": SolverInfo(
        name="multiterm_boltzmann",
        legacy_tag="multiterm",
        plot_label="Multi-term closure",
        plot_color="#2e7d32",
    ),
}

RUN_MODES = (
    "monte_carlo",
    "boltzmann_two_term",
    "multiterm_boltzmann",
    "both",
    "all",
)
PRIMARY_SOLVERS = tuple(SOLVER_REGISTRY)


def solver_legacy_tag(solver: str) -> str:
    info = SOLVER_REGISTRY.get(solver)
    return info.legacy_tag if info is not None else solver


def solver_plot_label(solver: str) -> str:
    info = SOLVER_REGISTRY.get(solver)
    return info.plot_label if info is not None else solver


def solver_plot_color(solver: str) -> str | None:
    info = SOLVER_REGISTRY.get(solver)
    return info.plot_color if info is not None else None
