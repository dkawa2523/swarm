from __future__ import annotations

import ast
from pathlib import Path
from types import ModuleType

import swarm_workflow.comsol.models.gec_ccp.plots.contracts as contracts
import swarm_workflow.comsol.models.gec_ccp.plots.data as data
import swarm_workflow.comsol.models.gec_ccp.plots.eedf as eedf
import swarm_workflow.comsol.models.gec_ccp.plots.spatial as spatial
import swarm_workflow.comsol.models.gec_ccp.plots.workflow as facade


DATA_HELPERS = (
    "_comparison_input_mph_sha256",
    "_comsol_column_names",
    "_comsol_domain_data",
    "_comsol_export_dimension",
    "_comsol_field_matrix",
    "_comsol_numeric_matrix",
    "_comsol_radial_midplane_data",
    "_comsol_sorted_phase_field",
    "_comsol_spatial_data",
    "_comsol_waveform_data",
    "_nearest",
    "_normalized_coordinate",
    "_numeric_rows",
    "_positive",
    "_rows",
    "_save",
    "_spatial_comparison_metrics",
    "_validated_plot_provenance",
    "_validated_result_artifacts",
    "_waveform_comparison_metrics",
)
SPATIAL_HELPERS = (
    "_draw_gec_boundaries",
    "_gec_domain_geometry",
    "_gec_domain_polygon",
    "_gec_domain_triangulation",
    "_plot_domain_2d_comparison",
    "_plot_domain_2d_difference",
    "_plot_gec_geometry_2d",
    "_plot_phase_axial_distributions",
    "_plot_radial_midplane_solver_comparison",
    "_plot_spatial_cut_overview",
    "_plot_waveform_comparison",
    "_style_comsol_surface_axis",
)
EEDF_HELPERS = (
    "_eedf_slice_metrics",
    "_eedf_tail_curve",
    "_function_eedf_slice",
    "_function_eedf_tail_curve",
    "_interpolate_tail_curve",
    "_nearest_eedf_slice",
    "_plot_operating_eedf",
    "_positive_log_interpolate",
    "_tail_probability",
    "_tail_survival",
    "_weighted_quantiles",
    "plot_closure_overview",
    "plot_operating_eedf_if_available",
)


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


def _function_line_count(function: object) -> int:
    module = __import__(function.__module__, fromlist=["*"])
    tree = ast.parse(Path(module.__file__).read_text(encoding="utf-8"))
    definition = next(
        node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == function.__name__
    )
    return definition.end_lineno - definition.lineno + 1


def test_plot_helpers_have_one_implementation_owner() -> None:
    modules = (facade, data, spatial, eedf)
    for owner, names in (
        (data, DATA_HELPERS),
        (spatial, SPATIAL_HELPERS),
        (eedf, EEDF_HELPERS),
    ):
        for name in names:
            assert name in owner.__dict__
            assert all(
                name not in module.__dict__
                for module in modules
                if module is not owner
            )


def test_plot_consumers_follow_their_dependency_direction() -> None:
    package = "swarm_workflow.comsol.models.gec_ccp.plots"

    def imports(module: ModuleType, *names: str) -> bool:
        expected = {f"{package}.{name}" for name in names}
        return bool(_imports(module) & expected)

    assert not imports(contracts, "data", "eedf", "spatial", "workflow")
    assert not imports(data, "eedf", "spatial", "workflow")
    assert not imports(spatial, "eedf", "workflow")
    assert not imports(eedf, "spatial", "workflow")


def test_model_layers_do_not_depend_on_plot_consumers() -> None:
    model_root = Path(facade.__file__).resolve().parents[1]
    for path in model_root.rglob("*.py"):
        if "plots" in path.parts:
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"))
        imports = {
            str(node.module)
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom) and node.module is not None
        }
        assert not any(
            name.startswith("swarm_workflow.comsol.models.gec_ccp.plots")
            for name in imports
        )


def test_plot_contracts_remain_available_from_public_facade() -> None:
    for name in contracts.__all__:
        assert getattr(facade, name) is getattr(contracts, name)


def test_plot_modules_remain_bounded_cohesive_units() -> None:
    for module in (contracts, data, spatial, eedf, facade):
        lines = Path(module.__file__).read_text(encoding="utf-8").splitlines()
        assert len(lines) <= 1_200, module.__name__
    assert _function_line_count(facade.plot_gec_ccp_results) <= 150
    assert _function_line_count(facade._plot_comsol_result_comparisons) <= 150
    assert _function_line_count(eedf.plot_closure_overview) <= 150
    assert _function_line_count(eedf.plot_operating_eedf_if_available) <= 150
