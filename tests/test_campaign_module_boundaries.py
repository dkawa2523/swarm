from __future__ import annotations

import ast
from pathlib import Path

import swarm_workflow
import swarm_workflow.campaign as campaign
import swarm_workflow.campaign.config as campaign_config
import swarm_workflow.campaign.sweep as campaign_sweep


ROOT = Path(__file__).resolve().parents[1]
OLD_CAMPAIGN_MODULES = {
    "swarm_workflow.aggregate",
    "swarm_workflow.aggregate_stats",
    "swarm_workflow.store",
    "swarm_workflow.sweep",
    "swarm_workflow.workflow_config",
    "swarm_workflow.workflow_repository",
}
PUBLIC_CAMPAIGN_ENTRYPOINTS = {
    "DeterministicExecutionConfig",
    "MeanEnergySupportConfig",
    "MixtureSpec",
    "SweepSummary",
    "WorkflowConfig",
    "load_workflow",
    "run_sweep",
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
    tree = ast.parse(path.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
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


def test_legacy_campaign_module_paths_are_removed() -> None:
    for module in OLD_CAMPAIGN_MODULES:
        path = ROOT.joinpath(*module.split(".")).with_suffix(".py")
        assert not path.exists(), path


def test_no_python_import_uses_a_legacy_campaign_module() -> None:
    for source_root in ("electron_swarm", "swarm_workflow", "tests", "tools"):
        for path in (ROOT / source_root).rglob("*.py"):
            assert _resolved_imports(path).isdisjoint(OLD_CAMPAIGN_MODULES), path


def test_campaign_package_does_not_depend_on_comsol() -> None:
    for path in (ROOT / "swarm_workflow" / "campaign").glob("*.py"):
        assert not any(
            imported == "swarm_workflow.comsol"
            or imported.startswith("swarm_workflow.comsol.")
            for imported in _resolved_imports(path)
        ), path


def test_campaign_package_exports_only_public_entrypoints() -> None:
    assert set(campaign.__all__) == PUBLIC_CAMPAIGN_ENTRYPOINTS
    assert campaign.DeterministicExecutionConfig is (
        campaign_config.DeterministicExecutionConfig
    )
    assert campaign.MeanEnergySupportConfig is campaign_config.MeanEnergySupportConfig
    assert campaign.MixtureSpec is campaign_config.MixtureSpec
    assert campaign.WorkflowConfig is campaign_config.WorkflowConfig
    assert campaign.load_workflow is campaign_config.load_workflow
    assert campaign.SweepSummary is campaign_sweep.SweepSummary
    assert campaign.run_sweep is campaign_sweep.run_sweep
    for name in PUBLIC_CAMPAIGN_ENTRYPOINTS:
        assert getattr(swarm_workflow, name) is getattr(campaign, name)
