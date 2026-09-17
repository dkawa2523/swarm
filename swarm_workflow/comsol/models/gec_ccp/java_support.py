"""Shared expressions and interpolation primitives for GEC-CCP Java sources."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Literal

from swarm_workflow.comsol.java import java_path, java_string


def java_header(class_name: str, input_mph: Path, model_name: str) -> list[str]:
    return [
        "import com.comsol.model.*;",
        "import com.comsol.model.util.*;",
        "",
        f"public class {class_name} {{",
        "  public static void main(String[] args) throws Exception {",
        "    Model model = ModelUtil.loadCopy("
        + java_string(model_name)
        + ", "
        + java_string(java_path(input_mph))
        + ");",
    ]


def java_boolean(value: bool) -> str:
    return "true" if value else "false"


def inline_interpolation_lines(
    tag: str,
    x_values: list[float],
    y_values: list[float],
) -> list[str]:
    table = ", ".join(
        "new String[]{"
        + java_string(f"{x_value:.17e}")
        + ", "
        + java_string(f"{y_value:.17e}")
        + "}"
        for x_value, y_value in zip(x_values, y_values, strict=True)
    )
    qtag = java_string(tag)
    return [
        f"    replaceInterpolation(model, {qtag});",
        f'    model.func({qtag}).set("source", "table");',
        f'    model.func({qtag}).set("nargs", 1);',
        f'    model.func({qtag}).set("argunit", "1");',
        f'    model.func({qtag}).set("fununit", "1");',
        f'    model.func({qtag}).set("interp", "piecewisecubic");',
        f'    model.func({qtag}).set("extrap", "const");',
        f'    model.func({qtag}).set("table", new String[][]{{{table}}});',
    ]


def data_export_lines(
    tag: str,
    dataset: str,
    expressions: Iterable[str],
    units: Iterable[str],
    output_csv: Path,
) -> list[str]:
    expr = ", ".join(java_string(value) for value in expressions)
    unit = ", ".join(java_string(value) for value in units)
    qtag = java_string(tag)
    return [
        f"    try {{ model.result().export().remove({qtag}); }} catch (Exception ex) {{}}",
        f'    model.result().export().create({qtag}, "Data");',
        f'    model.result().export({qtag}).set("data", {java_string(dataset)});',
        f'    model.result().export({qtag}).set("expr", new String[]{{{expr}}});',
        f'    model.result().export({qtag}).set("unit", new String[]{{{unit}}});',
        f'    model.result().export({qtag}).set("filename", '
        + java_string(java_path(output_csv))
        + ");",
        f"    model.result().export({qtag}).run();",
    ]


def numerical_table_export_lines(
    tag: str,
    numerical_type: Literal["IntSurface", "IntLine", "EvalGlobal"],
    dataset: str,
    expressions: Iterable[str],
    units: Iterable[str],
    descriptions: Iterable[str],
    entities: Iterable[int] | None,
    axisymmetric_measure: Literal["intvolume", "intsurface"] | None,
    output_csv: Path,
) -> list[str]:
    """Generate one auditable COMSOL Derived Values table export."""

    expr = ", ".join(java_string(value) for value in expressions)
    unit = ", ".join(java_string(value) for value in units)
    descr = ", ".join(java_string(value) for value in descriptions)
    selection = (
        ", ".join(str(int(value)) for value in entities) if entities is not None else ""
    )
    numerical_tag = java_string(tag)
    table_tag = java_string(f"{tag}Table")
    export_tag = java_string(f"{tag}Export")
    selection_lines = (
        [
            f"    model.result().numerical({numerical_tag}).selection().set("
            f"new int[]{{{selection}}});"
        ]
        if entities is not None
        else []
    )
    measure_lines = (
        [
            f"    model.result().numerical({numerical_tag}).set("
            f"{java_string(axisymmetric_measure)}, true);"
        ]
        if axisymmetric_measure is not None
        else []
    )
    return [
        f"    try {{ model.result().numerical().remove({numerical_tag}); }} "
        "catch (Exception ex) {}",
        f"    try {{ model.result().table().remove({table_tag}); }} "
        "catch (Exception ex) {}",
        f"    try {{ model.result().export().remove({export_tag}); }} "
        "catch (Exception ex) {}",
        f'    model.result().table().create({table_tag}, "Table");',
        f"    model.result().numerical().create({numerical_tag}, "
        f"{java_string(numerical_type)});",
        f"    model.result().numerical({numerical_tag}).set("
        f'"data", {java_string(dataset)});',
        f"    model.result().numerical({numerical_tag}).set("
        f'"expr", new String[]{{{expr}}});',
        f"    model.result().numerical({numerical_tag}).set("
        f'"unit", new String[]{{{unit}}});',
        f"    model.result().numerical({numerical_tag}).set("
        f'"descr", new String[]{{{descr}}});',
        *selection_lines,
        *measure_lines,
        f'    model.result().numerical({numerical_tag}).set("table", {table_tag});',
        f"    model.result().numerical({numerical_tag}).setResult();",
        f'    model.result().export().create({export_tag}, "Table");',
        f'    model.result().export({export_tag}).set("table", {table_tag});',
        f'    model.result().export({export_tag}).set("filename", '
        + java_string(java_path(output_csv))
        + ");",
        f"    model.result().export({export_tag}).run();",
    ]
