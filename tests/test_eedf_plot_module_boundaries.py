from __future__ import annotations

import ast
from pathlib import Path
from types import ModuleType

import swarm_workflow.plots.eedf_comparison as facade
import swarm_workflow.plots.eedf_contracts as contracts
import swarm_workflow.plots.eedf_data as data
import swarm_workflow.plots.eedf_metrics as metrics
import swarm_workflow.plots.eedf_provenance as provenance
import swarm_workflow.plots.eedf_render as render


def _imports(module: ModuleType) -> set[str]:
    tree = ast.parse(Path(module.__file__).read_text(encoding="utf-8"))
    names = {
        str(node.module)
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module is not None
    }
    names.update(
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    )
    return names


def test_eedf_plot_dependencies_point_toward_owners() -> None:
    prefix = "swarm_workflow.plots."
    assert not any(name.startswith(prefix) for name in _imports(contracts))
    assert f"{prefix}eedf_render" not in _imports(data)
    assert f"{prefix}eedf_comparison" not in _imports(data)
    assert f"{prefix}eedf_data" not in _imports(provenance)
    assert f"{prefix}eedf_render" not in _imports(provenance)
    assert f"{prefix}eedf_comparison" not in _imports(provenance)
    assert f"{prefix}eedf_render" not in _imports(metrics)
    assert f"{prefix}eedf_comparison" not in _imports(metrics)
    assert f"{prefix}eedf_comparison" not in _imports(render)


def test_eedf_plot_modules_remain_bounded_responsibility_units() -> None:
    for module in (contracts, provenance, data, metrics, render, facade):
        lines = Path(module.__file__).read_text(encoding="utf-8").splitlines()
        assert len(lines) <= 500, module.__name__
