"""Canonical product solver identifiers shared by config and execution."""

from typing import Final, Literal


SolverId = Literal["two_term", "multi_term", "monte_carlo", "propagator"]

CANONICAL_SOLVER_IDS: Final[tuple[SolverId, ...]] = (
    "two_term",
    "multi_term",
    "monte_carlo",
    "propagator",
)
