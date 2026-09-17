from __future__ import annotations

import ast
from pathlib import Path
from types import ModuleType

import electron_swarm.solvers.monte_carlo.audits as audits
import electron_swarm.solvers.monte_carlo as package
import electron_swarm.solvers.monte_carlo.batch as batch
import electron_swarm.solvers.monte_carlo.case as case
import electron_swarm.solvers.monte_carlo.case_phases as case_phases
import electron_swarm.solvers.monte_carlo.case_result as case_result
import electron_swarm.solvers.monte_carlo.case_state as case_state
import electron_swarm.solvers.monte_carlo.direct_transport as direct_transport
import electron_swarm.solvers.monte_carlo.population as population
import electron_swarm.solvers.monte_carlo.result_evidence as result_evidence
import electron_swarm.solvers.monte_carlo.setup as setup
import electron_swarm.solvers.monte_carlo.transport_common as transport_common
import electron_swarm.solvers.monte_carlo.weighted_transport as weighted_transport
import electron_swarm.solvers.monte_carlo.weighted_runtime as weighted_runtime
import electron_swarm.solvers.monte_carlo.solver as solver


MOVED_NAMES = {
    "_EnergyAudit",
    "_NullEnergyAudit",
    "_MonteCarloRunAudit",
    "_NullMonteCarloRunAudit",
    "_ParticleEnsemble",
    "_weighted_branching_active",
    "_ionization_branching_model",
    "_population_metadata",
    "_rate_results",
    "_weighted_growth_rate_consistency",
    "_CompiledKernelPlan",
    "_compiled_kernel_plan",
    "_sample_weighted_branching_collision",
    "_advance_compiled_weighted_barrier",
    "_run_weighted_branching_phase",
    "_MonteCarloRunSetup",
    "_prepare_monte_carlo_run",
    "_MonteCarloCaseState",
    "_initialize_monte_carlo_case",
    "_new_energy_audit",
    "_new_run_audit",
    "_run_fixed_particle_phase",
    "_advance_weighted_phase",
    "_run_tail_phase",
    "_run_weighted_phase",
    "_MonteCarloCaseStatistics",
    "_MonteCarloCaseEvidence",
    "_finalize_case_statistics",
    "_assemble_case_evidence",
    "_attach_audit_diagnostics_and_log",
    "_assemble_swarm_case_result",
    "_build_monte_carlo_case_result",
    "_run_monte_carlo_case",
    "DirectFluxTransport",
    "_weighted_mean",
    "direct_transport_observation_times",
    "_validate_snapshot_arrays",
    "direct_flux_transport_snapshot",
    "DirectFluxTransportAccumulator",
    "_PendingPlane",
    "_CompletePlane",
    "SynchronizedFluxTransportObserver",
    "WeightedGrowthTransportLagPlan",
    "weighted_growth_transport_lag_plan",
    "weighted_growth_flux_transport_snapshot",
    "_effective_lineage_count",
    "SynchronizedWeightedGrowthFluxObserver",
    "run_internal_monte_carlo",
}
OWNER_NAMES = {
    audits: {
        "_EnergyAudit",
        "_NullEnergyAudit",
        "_MonteCarloRunAudit",
        "_NullMonteCarloRunAudit",
    },
    population: {"_ParticleEnsemble"},
    result_evidence: {
        "_weighted_branching_active",
        "_ionization_branching_model",
        "_population_metadata",
        "_rate_results",
        "_weighted_growth_rate_consistency",
    },
    weighted_runtime: {
        "_CompiledKernelPlan",
        "_compiled_kernel_plan",
        "_sample_weighted_branching_collision",
        "_advance_compiled_weighted_barrier",
        "_run_weighted_branching_phase",
    },
    setup: {
        "_MonteCarloRunSetup",
        "_prepare_monte_carlo_run",
    },
    case_state: {
        "_MonteCarloCaseState",
        "_initialize_monte_carlo_case",
        "_new_energy_audit",
        "_new_run_audit",
    },
    case_phases: {
        "_run_fixed_particle_phase",
        "_advance_weighted_phase",
        "_run_tail_phase",
        "_run_weighted_phase",
    },
    case_result: {
        "_MonteCarloCaseStatistics",
        "_MonteCarloCaseEvidence",
        "_finalize_case_statistics",
        "_assemble_case_evidence",
        "_attach_audit_diagnostics_and_log",
        "_assemble_swarm_case_result",
        "_build_monte_carlo_case_result",
    },
    transport_common: {
        "DirectFluxTransport",
        "_weighted_mean",
    },
    direct_transport: {
        "direct_transport_observation_times",
        "_validate_snapshot_arrays",
        "direct_flux_transport_snapshot",
        "DirectFluxTransportAccumulator",
        "_PendingPlane",
        "_CompletePlane",
        "SynchronizedFluxTransportObserver",
    },
    weighted_transport: {
        "WeightedGrowthTransportLagPlan",
        "weighted_growth_transport_lag_plan",
        "weighted_growth_flux_transport_snapshot",
        "_effective_lineage_count",
        "SynchronizedWeightedGrowthFluxObserver",
    },
    case: {"_run_monte_carlo_case"},
    batch: {"run_internal_monte_carlo"},
}


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


def _function_lengths(module: ModuleType) -> dict[str, int]:
    return {
        node.name: node.end_lineno - node.lineno + 1
        for node in _tree(module).body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.end_lineno is not None
    }


def test_internal_mc_helpers_have_one_implementation_owner() -> None:
    assert _defined_names(solver).isdisjoint(MOVED_NAMES)
    for module, names in OWNER_NAMES.items():
        assert names <= _defined_names(module)


def test_internal_mc_solver_does_not_reexport_moved_helpers() -> None:
    assert MOVED_NAMES.isdisjoint(vars(solver))
    assert (MOVED_NAMES - OWNER_NAMES[batch]).isdisjoint(vars(batch))
    assert (MOVED_NAMES - OWNER_NAMES[case]).isdisjoint(vars(case))


def test_transport_owners_do_not_reexport_each_others_private_surface() -> None:
    common_names = OWNER_NAMES[transport_common]
    direct_names = OWNER_NAMES[direct_transport]
    weighted_names = OWNER_NAMES[weighted_transport]
    assert common_names.isdisjoint(vars(direct_transport))
    assert common_names.isdisjoint(vars(weighted_transport))
    assert direct_names.isdisjoint(vars(weighted_transport))
    assert weighted_names.isdisjoint(vars(direct_transport))


def test_internal_mc_public_solver_api_remains_available() -> None:
    assert package.__all__ == ["MonteCarloSolver"]
    assert package.MonteCarloSolver is solver.MonteCarloSolver
    assert solver.MonteCarloSolver.name == "monte_carlo"
    assert callable(batch.run_internal_monte_carlo)
    assert "electron_swarm.physics.kinetics" in _imported_modules(setup)


def test_internal_mc_helper_dependencies_do_not_point_back_to_solver() -> None:
    for module in OWNER_NAMES:
        assert "electron_swarm.solvers.monte_carlo.solver" not in (
            _imported_modules(module)
        )
        assert "electron_swarm.solvers.monte_carlo" not in _imported_modules(module)


def test_internal_mc_orchestrator_and_runtime_stay_bounded() -> None:
    assert len(Path(solver.__file__).read_text(encoding="utf-8").splitlines()) <= 200
    assert len(Path(batch.__file__).read_text(encoding="utf-8").splitlines()) <= 100
    assert len(Path(setup.__file__).read_text(encoding="utf-8").splitlines()) <= 300
    assert len(Path(case.__file__).read_text(encoding="utf-8").splitlines()) <= 100
    assert len(Path(case_state.__file__).read_text(encoding="utf-8").splitlines()) <= 250
    assert len(Path(case_phases.__file__).read_text(encoding="utf-8").splitlines()) <= 500
    assert len(Path(case_result.__file__).read_text(encoding="utf-8").splitlines()) <= 600
    assert (
        len(Path(transport_common.__file__).read_text(encoding="utf-8").splitlines())
        <= 100
    )
    assert (
        len(Path(direct_transport.__file__).read_text(encoding="utf-8").splitlines())
        <= 650
    )
    assert (
        len(Path(weighted_transport.__file__).read_text(encoding="utf-8").splitlines())
        <= 750
    )
    assert (
        len(Path(weighted_runtime.__file__).read_text(encoding="utf-8").splitlines())
        <= 1_200
    )


def test_internal_mc_case_stage_functions_stay_cohesive() -> None:
    for module in (case, case_state, case_phases, case_result):
        assert max(_function_lengths(module).values()) <= 250


def test_legacy_transport_moments_owner_is_removed() -> None:
    package_dir = Path(package.__file__).resolve().parent
    assert not (package_dir / "transport_moments.py").exists()
