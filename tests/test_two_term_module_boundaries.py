from __future__ import annotations

import ast
from pathlib import Path
from types import ModuleType

import electron_swarm.solvers.two_term as package
import electron_swarm.solvers.boltzmann_common.collisions as common_collisions
import electron_swarm.solvers.boltzmann_common.grid as common_grid
import electron_swarm.solvers.boltzmann_common.observables as common_observables
import electron_swarm.solvers.boltzmann_common.operators as common_operators
import electron_swarm.solvers.two_term.grid as grid
import electron_swarm.solvers.two_term.models as models
import electron_swarm.solvers.two_term.observables as observables
import electron_swarm.solvers.two_term.solver as solver
import electron_swarm.solvers.two_term.steady as steady
import electron_swarm.solvers.two_term.time_periodic as time_periodic
import electron_swarm.solvers.two_term.transport as transport


def _tree(module: ModuleType) -> ast.Module:
    return ast.parse(Path(module.__file__).read_text(encoding="utf-8"))


def _defined_names(module: ModuleType) -> set[str]:
    return {
        node.name
        for node in _tree(module).body
        if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
    }


def _imported_modules(module: ModuleType) -> set[str]:
    imported: set[str] = set()
    for node in ast.walk(_tree(module)):
        if isinstance(node, ast.Import):
            imported.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            imported.add(node.module)
    return imported


def test_two_term_public_package_only_exports_solver() -> None:
    assert package.__all__ == ["TwoTermSolver"]
    assert package.TwoTermSolver is solver.TwoTermSolver


def test_two_term_implementation_has_one_responsibility_owner() -> None:
    assert {
        "NativeSolveDiagnostics",
        "NativeDistributionResult",
        "TimePeriodicDistributionResult",
    } <= _defined_names(models)
    assert {"make_energy_grid", "grid_metadata", "maxwell_eedf"} <= _defined_names(grid)
    assert {
        "solve_native_distribution",
        "assemble_native_operator_block",
        "solve_native_on_grid",
        "solve_normalized",
    } <= _defined_names(steady)
    assert {"solve_time_periodic_distribution"} <= _defined_names(time_periodic)
    assert {"build_case_result"} <= _defined_names(observables)
    assert "transport_from_eedf" not in _defined_names(observables)
    assert {
        "transport_from_eedf",
        "temporal_growth_effective_momentum_frequency",
        "discrete_elastic_energy_loss_rate_coefficient",
    } <= _defined_names(transport)
    assert _defined_names(solver) == {"TwoTermSolver"}


def test_boltzmann_common_has_explicit_solver_neutral_owners() -> None:
    assert {"KineticGrid", "cell_edges_from_centers", "electron_speed_m_s"} <= (
        _defined_names(common_grid)
    )
    assert {
        "EffectiveCollisionData",
        "build_effective_collision_data",
        "assemble_collision_operator",
    } <= _defined_names(common_collisions)
    assert {"KineticOperatorBlock", "assemble_energy_flux_operator"} <= (
        _defined_names(common_operators)
    )
    assert {
        "normalize_eedf",
        "mean_energy_from_eedf",
        "compute_rates_from_eedf",
    } <= _defined_names(common_observables)
    for module in (
        common_grid,
        common_collisions,
        common_operators,
        common_observables,
    ):
        imports = _imported_modules(module)
        assert all(
            not name.startswith("electron_swarm.solvers.two_term")
            and not name.startswith("electron_swarm.solvers.multi_term")
            for name in imports
        )


def test_two_term_internal_dependencies_do_not_point_back_to_adapter() -> None:
    adapter = "electron_swarm.solvers.two_term.solver"
    for module in (models, grid, steady, time_periodic, transport, observables):
        assert adapter not in _imported_modules(module)


def test_two_term_obsolete_private_forwarders_are_removed() -> None:
    defined = set().union(
        *(
            _defined_names(module)
            for module in (grid, steady, time_periodic, transport, observables, solver)
        )
    )
    assert "_gas_number_density" not in defined
    assert "_cell_edges_from_centers" not in defined


def test_two_term_adapter_remains_thin() -> None:
    lines = Path(solver.__file__).read_text(encoding="utf-8").splitlines()
    assert len(lines) <= 220
