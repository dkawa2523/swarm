"""Canonical solver identifiers and compact product metadata."""

from __future__ import annotations

from typing import TYPE_CHECKING

CANONICAL_SOLVER_IDS = ("two_term", "multi_term", "monte_carlo")

if TYPE_CHECKING:  # pragma: no cover
    from electron_swarm.core.config import SwarmConfig


def solver_method(config: "SwarmConfig", solver: str) -> str:
    """Return the user-facing method label for a canonical solver id."""

    if solver == "two_term":
        return str(config.solvers.two_term.backend)
    if solver == "multi_term":
        return str(config.solvers.multi_term.method)
    if solver == "monte_carlo":
        return "internal"
    raise ValueError(f"Unknown canonical solver id: {solver!r}")
