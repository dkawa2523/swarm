from __future__ import annotations

from importlib import import_module
from pathlib import Path

import pytest

from electron_swarm.solvers.base import IndependentCaseSolver, SwarmSolver


ROOT = Path(__file__).resolve().parents[1]
SOLVERS = {
    "two_term": "TwoTermSolver",
    "multi_term": "MultiTermSolver",
    "monte_carlo": "MonteCarloSolver",
    "propagator": "PropagatorSolver",
}


@pytest.mark.parametrize(("solver_id", "class_name"), SOLVERS.items())
def test_canonical_solver_has_one_package_entry_point(
    solver_id: str,
    class_name: str,
) -> None:
    package = import_module(f"electron_swarm.solvers.{solver_id}")

    assert package.__all__ == [class_name]
    assert issubclass(getattr(package, class_name), SwarmSolver)


def test_solver_execution_contract_matches_each_solver_granularity() -> None:
    assert SwarmSolver.__abstractmethods__ == frozenset({"solve_all"})
    assert IndependentCaseSolver.__abstractmethods__ == frozenset({"solve_case"})

    for solver_id in ("two_term", "multi_term", "propagator"):
        package = import_module(f"electron_swarm.solvers.{solver_id}")
        assert issubclass(
            getattr(package, SOLVERS[solver_id]),
            IndependentCaseSolver,
        )

    monte_carlo = import_module("electron_swarm.solvers.monte_carlo")
    monte_carlo_type = getattr(monte_carlo, SOLVERS["monte_carlo"])
    assert issubclass(monte_carlo_type, SwarmSolver)
    assert not issubclass(monte_carlo_type, IndependentCaseSolver)
    assert "solve_case" not in monte_carlo_type.__dict__


def test_obsolete_solver_module_paths_are_removed() -> None:
    solver_root = ROOT / "electron_swarm" / "solvers"
    for relative in (
        "two_term.py",
        "internal_monte_carlo.py",
        "monte_carlo_evidence.py",
        "kinetic.py",
        "_internal_mc/__init__.py",
    ):
        assert not (solver_root / relative).exists()


def test_propagator_package_does_not_depend_on_other_solver_implementations() -> None:
    propagator_root = ROOT / "electron_swarm" / "solvers" / "propagator"
    forbidden = (
        "electron_swarm.solvers.two_term",
        "electron_swarm.solvers.multi_term",
        "electron_swarm.solvers.monte_carlo",
        "electron_swarm.solvers.boltzmann_common",
    )
    for path in propagator_root.glob("*.py"):
        source = path.read_text(encoding="utf-8")
        assert all(name not in source for name in forbidden)
