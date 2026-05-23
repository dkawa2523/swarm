"""Canonical solver identifiers and presentation metadata."""

from __future__ import annotations

from dataclasses import dataclass

CANONICAL_SOLVER_IDS = ("two_term", "multi_term", "monte_carlo")


@dataclass(frozen=True, slots=True)
class SolverInfo:
    name: str
    plot_label: str
    plot_color: str


SOLVER_REGISTRY: dict[str, SolverInfo] = {
    "two_term": SolverInfo(
        name="two_term",
        plot_label="Two-term",
        plot_color="#1565c0",
    ),
    "multi_term": SolverInfo(
        name="multi_term",
        plot_label="Multi-term closure",
        plot_color="#2e7d32",
    ),
    "monte_carlo": SolverInfo(
        name="monte_carlo",
        plot_label="Monte Carlo",
        plot_color="#c62828",
    ),
}


def solver_plot_label(solver: str) -> str:
    info = SOLVER_REGISTRY.get(solver)
    return info.plot_label if info is not None else solver


def solver_plot_color(solver: str) -> str | None:
    info = SOLVER_REGISTRY.get(solver)
    return info.plot_color if info is not None else None
