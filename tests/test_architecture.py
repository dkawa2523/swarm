"""Protect dependency direction without freezing private names or file sizes."""

import ast
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _imports(path: Path) -> set[str]:
    package = list(path.relative_to(ROOT).with_suffix("").parts[:-1])
    result = set()
    for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
        if isinstance(node, ast.Import):
            result.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            prefix = package[:len(package) - node.level + 1] if node.level else []
            module = ".".join([*prefix, *([node.module] if node.module else [])])
            result.add(module)
            if node.module is None:
                result.update(f"{module}.{alias.name}" for alias in node.names)
    return result


@pytest.mark.parametrize("directory, forbidden", [
    ("electron_swarm", ("swarm_workflow", "sqlite3", "subprocess", "comsol", "mph")),
    ("electron_swarm/solvers", ("electron_swarm.orchestration", "electron_swarm.io")),
    ("swarm_workflow/campaign", ("swarm_workflow.comsol",)),
    ("swarm_workflow/selection", ("swarm_workflow.comsol",)),
    ("swarm_workflow/tables", ("swarm_workflow.comsol",)),
    ("swarm_workflow/quality", ("swarm_workflow.comsol", "swarm_workflow.tables")),
    ("swarm_workflow/comsol/input", ("swarm_workflow.comsol.models", "swarm_workflow.comsol.runtime")),
    ("swarm_workflow/comsol/runtime", ("swarm_workflow.comsol.models",)),
])
def test_dependency_direction(directory: str, forbidden: tuple[str, ...]) -> None:
    for path in (ROOT / directory).rglob("*.py"):
        violations = {name for name in _imports(path)
                      if any(name == prefix or name.startswith(prefix + ".") for prefix in forbidden)}
        assert not violations, (path.relative_to(ROOT), violations)


def test_numerical_solvers_do_not_import_each_other() -> None:
    modes = ("two_term", "multi_term", "monte_carlo", "propagator")
    for mode in modes:
        for path in (ROOT / "electron_swarm" / "solvers" / mode).rglob("*.py"):
            for name in _imports(path):
                assert not any(name.startswith(f"electron_swarm.solvers.{other}.")
                               or name == f"electron_swarm.solvers.{other}"
                               for other in modes if other != mode), (path, name)
