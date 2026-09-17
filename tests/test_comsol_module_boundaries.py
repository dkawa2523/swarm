from __future__ import annotations

import ast
from pathlib import Path

import swarm_workflow.comsol.input as comsol_input
import swarm_workflow.comsol.input.eedf_audit as eedf_audit
import swarm_workflow.comsol.input.eedf_audit.analysis as eedf_audit_analysis
import swarm_workflow.comsol.input.eedf_audit.contracts as eedf_audit_contracts
import swarm_workflow.comsol.input.eedf_audit.io as eedf_audit_io
import swarm_workflow.comsol.input.eedf_audit.java as eedf_audit_java
import swarm_workflow.comsol.input.eedf_audit.planning as eedf_audit_planning
import swarm_workflow.comsol.input.export as input_export
import swarm_workflow.comsol.input.export.bundle as export_bundle
import swarm_workflow.comsol.input.function_eedf as function_eedf
import swarm_workflow.comsol.input.function_eedf.c1 as eedf_c1
import swarm_workflow.comsol.input.function_eedf.io as eedf_io
import swarm_workflow.comsol.input.function_eedf.kernels as eedf_kernels
import swarm_workflow.comsol.input.function_eedf.moments as eedf_moments
import swarm_workflow.comsol.runtime as runtime
import swarm_workflow.selection as selection


ROOT = Path(__file__).resolve().parents[1]
OLD_INPUT_MODULES = {
    "swarm_workflow.closure_selection",
    "swarm_workflow.comsol_eedf_audit",
    "swarm_workflow.comsol_export",
    "swarm_workflow.function_eedf",
    "swarm_workflow.mean_energy",
    "swarm_workflow.comsol.input.anchor_fallback",
    "swarm_workflow.comsol.input.selection",
}
PUBLIC_INPUT_ENTRYPOINTS = {
    "ComsolExportError",
    "ComsolExportSummary",
    "export_comsol_bundle",
}
PUBLIC_EXPORT_ENTRYPOINTS = {
    "ComsolExportError",
    "ComsolExportSummary",
    "export_comsol_bundle",
}
PUBLIC_FUNCTION_EEDF_ENTRYPOINTS = {
    "COMSOL_EEDF_COLUMNS",
    "MEAN_AXIS_RATE_IMPORTANCE_FRACTION",
    "SOURCE_MEAN_RELATIVE_ERROR_LIMIT",
    "RATE_IMPORTANCE_FRACTION",
    "RATE_SCALED_ERROR_TOLERANCE",
    "SHAPE_TOTAL_VARIATION_TOLERANCE",
    "C1FunctionEedf",
    "COLLISION_RATE_KERNEL_COLUMNS",
    "COLLISION_RATE_KERNEL_TABLE",
    "CollisionRateKernel",
    "ComsolEedfImportContract",
    "ComsolFunctionEedfGrid",
    "FunctionEedfError",
    "build_c1_function_eedf",
    "collision_rate_coefficient",
    "evaluate_c1_function_eedf",
    "evaluate_comsol_function_eedf_grid",
    "pchip_weighted_moments",
    "piecewise_linear_weighted_moments",
    "project_c1_function_eedf_to_comsol_grid",
    "read_c1_function_eedf",
    "read_collision_rate_kernels",
    "read_comsol_function_eedf_grid",
    "scaled_shape_moments",
}
PUBLIC_EEDF_AUDIT_ENTRYPOINTS = {
    "ComsolEedfAuditError",
    "ComsolEedfAuditPlan",
    "analyze_comsol_eedf_audit",
    "extract_comsol_eedf_audit_log",
    "prepare_comsol_eedf_audit",
    "read_native_eedf_grid",
    "render_comsol_eedf_audit_java",
}
REMOVED_INPUT_MONOLITHS = (
    ROOT / "swarm_workflow" / "comsol" / "input" / "eedf_audit.py",
    ROOT / "swarm_workflow" / "comsol" / "input" / "export.py",
    ROOT / "swarm_workflow" / "comsol" / "input" / "function_eedf.py",
)


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


def test_comsol_runtime_does_not_depend_on_model_packages() -> None:
    runtime_dir = Path(runtime.__file__).parent
    imported_modules: set[str] = set()
    for path in runtime_dir.glob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        imported_modules.update(
            node.module
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom) and node.module is not None
        )
        imported_modules.update(
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.Import)
            for alias in node.names
        )

    assert not any("comsol.models" in name for name in imported_modules)


def test_selection_package_is_comsol_independent() -> None:
    selection_dir = Path(selection.__file__).parent
    for path in selection_dir.glob("*.py"):
        assert not any(
            imported == "swarm_workflow.comsol"
            or imported.startswith("swarm_workflow.comsol.")
            for imported in _resolved_imports(path)
        ), path


def test_legacy_comsol_input_module_paths_are_removed() -> None:
    for module in OLD_INPUT_MODULES:
        path = ROOT.joinpath(*module.split(".")).with_suffix(".py")
        assert not path.exists(), path


def test_no_python_import_uses_a_legacy_comsol_input_module() -> None:
    for source_root in ("electron_swarm", "swarm_workflow", "tests", "tools"):
        for path in (ROOT / source_root).rglob("*.py"):
            assert _resolved_imports(path).isdisjoint(OLD_INPUT_MODULES), path


def test_comsol_input_does_not_depend_on_models_or_runtime() -> None:
    forbidden = (
        "swarm_workflow.comsol.models",
        "swarm_workflow.comsol.runtime",
    )
    for path in (ROOT / "swarm_workflow" / "comsol" / "input").rglob("*.py"):
        assert not any(
            imported == prefix or imported.startswith(f"{prefix}.")
            for prefix in forbidden
            for imported in _resolved_imports(path)
        ), path


def test_comsol_input_exports_only_public_entrypoints() -> None:
    assert set(comsol_input.__all__) == PUBLIC_INPUT_ENTRYPOINTS
    assert comsol_input.ComsolExportError is input_export.ComsolExportError
    assert comsol_input.ComsolExportSummary is input_export.ComsolExportSummary
    assert comsol_input.export_comsol_bundle is input_export.export_comsol_bundle


def test_input_monoliths_are_replaced_by_real_packages() -> None:
    assert all(not path.exists() for path in REMOVED_INPUT_MONOLITHS)
    for name in ("eedf_audit", "export", "function_eedf"):
        package = ROOT / "swarm_workflow" / "comsol" / "input" / name
        assert package.is_dir()
        assert (package / "__init__.py").is_file()


def test_nested_input_packages_export_only_public_entrypoints() -> None:
    assert set(eedf_audit.__all__) == PUBLIC_EEDF_AUDIT_ENTRYPOINTS
    assert set(input_export.__all__) == PUBLIC_EXPORT_ENTRYPOINTS
    assert set(function_eedf.__all__) == PUBLIC_FUNCTION_EEDF_ENTRYPOINTS
    assert eedf_audit.ComsolEedfAuditError is (
        eedf_audit_contracts.ComsolEedfAuditError
    )
    assert eedf_audit.ComsolEedfAuditPlan is eedf_audit_contracts.ComsolEedfAuditPlan
    assert eedf_audit.analyze_comsol_eedf_audit is (
        eedf_audit_analysis.analyze_comsol_eedf_audit
    )
    assert eedf_audit.extract_comsol_eedf_audit_log is (
        eedf_audit_java.extract_comsol_eedf_audit_log
    )
    assert eedf_audit.prepare_comsol_eedf_audit is (
        eedf_audit_planning.prepare_comsol_eedf_audit
    )
    assert eedf_audit.read_native_eedf_grid is eedf_audit_io.read_native_eedf_grid
    assert eedf_audit.render_comsol_eedf_audit_java is (
        eedf_audit_java.render_comsol_eedf_audit_java
    )
    assert input_export.export_comsol_bundle is export_bundle.export_comsol_bundle
    assert function_eedf.build_c1_function_eedf is eedf_c1.build_c1_function_eedf
    assert function_eedf.read_comsol_function_eedf_grid is (
        eedf_io.read_comsol_function_eedf_grid
    )
    assert function_eedf.piecewise_linear_weighted_moments is (
        eedf_moments.piecewise_linear_weighted_moments
    )
    assert function_eedf.read_collision_rate_kernels is (
        eedf_kernels.read_collision_rate_kernels
    )


def test_input_implementation_files_remain_bounded() -> None:
    input_root = ROOT / "swarm_workflow" / "comsol" / "input"
    for package_name in ("eedf_audit", "export", "function_eedf"):
        for path in (input_root / package_name).glob("*.py"):
            assert len(path.read_text(encoding="utf-8").splitlines()) <= 650, path
