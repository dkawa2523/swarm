from __future__ import annotations

import ast
from pathlib import Path

import swarm_workflow.quality as quality
import swarm_workflow.quality.monte_carlo as mc_quality
import swarm_workflow.tables as tables
import swarm_workflow.tables.builder as table_builder


ROOT = Path(__file__).resolve().parents[1]
OLD_FLAT_MODULES = {
    "swarm_workflow.mc_direct_transport_quality",
    "swarm_workflow.mc_policy",
    "swarm_workflow.mc_quality_contracts",
    "swarm_workflow.mc_weighted_transport_quality",
    "swarm_workflow.quality_policy",
    "swarm_workflow.quality_table",
    "swarm_workflow.solver_qualification",
    "swarm_workflow.table_contracts",
    "swarm_workflow.table_energy_loss",
    "swarm_workflow.table_math",
    "swarm_workflow.table_mc_evidence",
    "swarm_workflow.table_repository",
}


def _module_name(path: Path) -> str:
    relative = path.relative_to(ROOT).with_suffix("")
    parts = list(relative.parts)
    if parts[-1] == "__init__":
        parts.pop()
    return ".".join(parts)


def _resolved_imports(path: Path) -> set[str]:
    module = _module_name(path)
    package = module if path.name == "__init__.py" else module.rpartition(".")[0]
    imports: set[str] = set()
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if isinstance(node, ast.Import):
            imports.update(alias.name for alias in node.names)
            continue
        if not isinstance(node, ast.ImportFrom):
            continue
        if node.level == 0:
            if node.module is not None:
                imports.add(node.module)
            continue
        package_parts = package.split(".") if package else []
        keep = max(0, len(package_parts) - node.level + 1)
        base = ".".join(package_parts[:keep])
        if node.module is not None:
            imports.add(".".join(part for part in (base, node.module) if part))
        else:
            imports.update(
                ".".join(part for part in (base, alias.name) if part)
                for alias in node.names
            )
    return imports


def test_old_flat_quality_and_table_modules_are_absent() -> None:
    for module in OLD_FLAT_MODULES:
        assert not ROOT.joinpath(*module.split(".")).with_suffix(".py").exists()
    assert not (ROOT / "swarm_workflow" / "tables.py").exists()


def test_no_python_import_uses_an_old_flat_module() -> None:
    for source_root in ("electron_swarm", "swarm_workflow", "tests", "tools"):
        for path in (ROOT / source_root).rglob("*.py"):
            assert _resolved_imports(path).isdisjoint(OLD_FLAT_MODULES), path


def test_quality_does_not_depend_on_tables_or_comsol() -> None:
    for path in (ROOT / "swarm_workflow" / "quality").rglob("*.py"):
        for imported in _resolved_imports(path):
            assert imported != "swarm_workflow.tables"
            assert not imported.startswith("swarm_workflow.tables.")
            assert imported != "swarm_workflow.comsol"
            assert not imported.startswith("swarm_workflow.comsol.")


def test_tables_do_not_depend_on_comsol() -> None:
    for path in (ROOT / "swarm_workflow" / "tables").rglob("*.py"):
        for imported in _resolved_imports(path):
            assert imported != "swarm_workflow.comsol"
            assert not imported.startswith("swarm_workflow.comsol.")


def test_new_packages_expose_only_genuine_entrypoints() -> None:
    assert set(quality.__all__) == {"QualityThresholds", "RequiredRateRse"}
    assert set(mc_quality.__all__) == {
        "FailureAxis",
        "MonteCarloConvergencePolicy",
        "MonteCarloPolicyError",
        "PolicyDecisionSummary",
        "SamplingPlanEntry",
        "decide_monte_carlo_closure",
        "failure_axes_for_reasons",
    }
    assert tables.__all__ == [
        "CompositeTableSummary",
        "build_tables",
        "compose_mc_anchor_fallback_tables",
    ]
    assert tables.build_tables is table_builder.build_tables
