from __future__ import annotations

import ast
from pathlib import Path
from types import ModuleType

import swarm_workflow.quality.monte_carlo as mc_quality
import swarm_workflow.quality.monte_carlo.contracts as contracts
import swarm_workflow.quality.monte_carlo.direct_transport as direct_quality
import swarm_workflow.quality.monte_carlo.weighted_transport as weighted_quality
import swarm_workflow.tables.repository as repository


ROOT = Path(__file__).resolve().parents[1]
REMOVED_MODULE = "mc_transport_quality"
DIRECT_SUMMARY = "summarize_mc_transport_diagnostics"
WEIGHTED_SUMMARY = "summarize_weighted_mc_transport_diagnostics"
CONTRACT_NAMES = {
    "MC_TRANSPORT_STATE_RELATIVE_TOLERANCE",
    "MC_TRANSPORT_DIFFUSION_LAG_RELATIVE_TOLERANCE",
    "WEIGHTED_MC_TRANSPORT_DIFFUSION_IDENTIFICATION_RELATIVE_TOLERANCE",
    "MC_EEDF_RATE_ENSEMBLE_RELATIVE_TOLERANCE",
    "MC_POPULATION_GROWTH_RELATIVE_TOLERANCE",
    "MC_TRANSPORT_MIN_ENSEMBLE_REPLICATES",
    "MC_TRANSPORT_STATE_STATIONARITY_FIELDS",
    "MC_TRANSPORT_LAG_FIELDS",
    "MC_TRANSPORT_PRODUCTION_FIELDS",
    "MC_TRANSPORT_COMPONENT_ESTIMATORS",
    "WEIGHTED_MC_TRANSPORT_STATIONARITY_FIELDS",
    "WEIGHTED_MC_TRANSPORT_COMPONENT_ESTIMATORS",
    "MC_TRANSPORT_SAMPLING_FIELDS",
    "MC_TRANSPORT_RUN_PROVENANCE_FIELDS",
    "quality_failure_reasons",
    "mc_function_eedf_restricted_lmea_failure_reasons",
    "paired_log_ratio_ci95_bound",
}


def _tree(module: ModuleType) -> ast.Module:
    return ast.parse(Path(module.__file__).read_text(encoding="utf-8"))


def _imports(module: ModuleType) -> set[str]:
    names = {
        node.module
        for node in ast.walk(_tree(module))
        if isinstance(node, ast.ImportFrom) and node.module is not None
    }
    names.update(
        alias.name
        for node in ast.walk(_tree(module))
        if isinstance(node, ast.Import)
        for alias in node.names
    )
    return names


def test_mc_quality_symbols_have_one_owner() -> None:
    assert CONTRACT_NAMES <= vars(contracts).keys()
    assert direct_quality.summarize_mc_transport_diagnostics.__module__ == (
        direct_quality.__name__
    )
    assert weighted_quality.summarize_weighted_mc_transport_diagnostics.__module__ == (
        weighted_quality.__name__
    )

    assert DIRECT_SUMMARY not in vars(contracts)
    assert WEIGHTED_SUMMARY not in vars(contracts)
    assert WEIGHTED_SUMMARY not in vars(direct_quality)
    assert DIRECT_SUMMARY not in vars(weighted_quality)
    assert DIRECT_SUMMARY not in vars(mc_quality)
    assert WEIGHTED_SUMMARY not in vars(mc_quality)
    repository_imports = _imports(repository)
    assert "quality.monte_carlo.direct_transport" in repository_imports
    assert "quality.monte_carlo.weighted_transport" in repository_imports


def test_removed_mc_quality_module_has_no_importers() -> None:
    assert not (ROOT / "swarm_workflow" / f"{REMOVED_MODULE}.py").exists()
    for path in (ROOT / "swarm_workflow").rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        imported = {
            node.module
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom) and node.module is not None
        }
        imported.update(
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.Import)
            for alias in node.names
        )
        assert all(not name.endswith(REMOVED_MODULE) for name in imported), path


def test_mc_quality_dependencies_point_from_summaries_to_contracts() -> None:
    contract_imports = _imports(contracts)
    direct_imports = _imports(direct_quality)
    weighted_imports = _imports(weighted_quality)

    assert "contracts" in direct_imports
    assert "contracts" in weighted_imports
    assert "direct_transport" not in contract_imports
    assert "weighted_transport" not in contract_imports
    assert "weighted_transport" not in direct_imports
    assert "direct_transport" not in weighted_imports
    for imported in contract_imports | direct_imports | weighted_imports:
        assert (
            imported == "electron_swarm.solvers.monte_carlo.evidence"
            or not imported.startswith("electron_swarm.solvers.monte_carlo.")
        )


def test_mc_quality_modules_stay_below_architecture_size_limit() -> None:
    for module in (contracts, direct_quality, weighted_quality):
        physical_lines = len(
            Path(module.__file__).read_text(encoding="utf-8").splitlines()
        )
        assert physical_lines <= 1_200, module.__name__


def test_mc_quality_stage_functions_stay_bounded() -> None:
    summary_limits = {
        DIRECT_SUMMARY: 150,
        WEIGHTED_SUMMARY: 150,
    }
    for module in (direct_quality, weighted_quality):
        functions = {
            node.name: node.end_lineno - node.lineno + 1
            for node in ast.walk(_tree(module))
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        }
        for name, length in functions.items():
            assert length <= summary_limits.get(name, 250), (
                module.__name__,
                name,
                length,
            )
