"""Verify COMSOL interpolation functions against mapped Swarm CSV tables."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any, Iterable

from ._io import write_csv as _write_csv
from ._io import write_json as _write_json
from .comsol_adapter import ComsolExecutionSummary, execute_generated_comsol_java
from .comsol_java import java_path, java_string
from .comsol_mapping import (
    ComsolModelMapping,
    FunctionMapping,
    load_comsol_mapping,
)


EXPECTED_CSV = "function_values_expected.csv"
COMSOL_CSV = "function_values_comsol.csv"
SUMMARY_CSV = "function_verify_summary.csv"
SUMMARY_JSON = "function_verify_summary.json"
VERIFY_MANIFEST = "verify_manifest.json"
VERIFY_CLASS_NAME = "SwarmComsolVerify"


class ComsolVerifyError(RuntimeError):
    """Raised when COMSOL function verification cannot be planned or fails."""


@dataclass(frozen=True, slots=True)
class FunctionControlPoint:
    function_name: str
    function_tag: str
    file: str
    column: str
    argument_column: str
    argument: float
    expected_value: float
    argunit: str
    fununit: str
    row_index: int
    reason: str


@dataclass(frozen=True, slots=True)
class FeatureTableControlPoint:
    """One value read back from a COMSOL physics-feature lookup array."""

    quantity: str
    component: str
    physics: str
    feature: str
    x_property: str
    y_property: str
    file: str
    column: str
    argument_column: str
    argument: float
    expected_value: float
    argunit: str
    fununit: str
    row_index: int
    table_length: int
    reason: str = "feature_table_row"

    @property
    def function_name(self) -> str:
        return self.quantity

    @property
    def function_tag(self) -> str:
        return f"{self.feature}.{self.y_property}"


@dataclass(frozen=True, slots=True)
class VerifyComsolPlan:
    mapping: ComsolModelMapping
    points: tuple[FunctionControlPoint, ...]
    feature_points: tuple[FeatureTableControlPoint, ...]
    output_dir: Path
    expected_csv: Path
    comsol_csv: Path
    summary_csv: Path
    summary_json: Path
    manifest_path: Path
    java_path: Path
    relative_tolerance: float
    absolute_tolerance: float


@dataclass(frozen=True, slots=True)
class VerifyComsolSummary:
    plan: VerifyComsolPlan
    passed: int
    failed: int
    expected_csv: Path
    comsol_csv: Path
    summary_csv: Path
    summary_json: Path
    java_path: Path
    execution: ComsolExecutionSummary | Any | None


def plan_verify_comsol_functions(mapping_path: str | Path) -> VerifyComsolPlan:
    mapping = load_comsol_mapping(mapping_path)
    _validate_verify_sources(mapping)
    output_dir = mapping.verify.output_path
    return VerifyComsolPlan(
        mapping=mapping,
        points=tuple(_control_points(mapping)),
        feature_points=tuple(_feature_table_points(mapping)),
        output_dir=output_dir,
        expected_csv=output_dir / EXPECTED_CSV,
        comsol_csv=output_dir / COMSOL_CSV,
        summary_csv=output_dir / SUMMARY_CSV,
        summary_json=output_dir / SUMMARY_JSON,
        manifest_path=output_dir / VERIFY_MANIFEST,
        java_path=output_dir / f"{VERIFY_CLASS_NAME}.java",
        relative_tolerance=mapping.verify.relative_tolerance,
        absolute_tolerance=mapping.verify.absolute_tolerance,
    )


def execute_verify_comsol_functions(
    mapping_path: str | Path,
    *,
    comsol_executable: str | Path | None = None,
) -> VerifyComsolSummary:
    plan = plan_verify_comsol_functions(mapping_path)
    if not plan.mapping.model.output_mph.exists():
        write_verify_manifest(plan, status="failed", failure_reason="missing_output_mph")
        raise ComsolVerifyError(
            f"COMSOL output .mph does not exist: {plan.mapping.model.output_mph}"
        )
    plan.output_dir.mkdir(parents=True, exist_ok=True)
    _clear_previous_verify_outputs(plan)
    write_expected_values_csv(plan)
    write_verify_manifest(plan, status="running")
    source = generate_verify_java_source(plan)
    plan.java_path.write_text(source, encoding="utf-8")
    execution = execute_generated_comsol_java(
        plan.mapping,
        plan.java_path,
        operation="verify",
        comsol_executable=comsol_executable,
    )
    if not plan.comsol_csv.exists():
        _extract_comsol_csv_from_stdout(execution.stdout_paths[-1], plan.comsol_csv)
    if not plan.comsol_csv.exists():
        raise ComsolVerifyError(
            "COMSOL verification did not produce "
            f"{plan.comsol_csv}; see log_dir: {execution.log_dir}"
        )
    summary = compare_function_values(plan)
    write_verify_manifest(
        plan,
        status="passed" if summary.failed == 0 else "failed",
        passed=summary.passed,
        failed=summary.failed,
    )
    if summary.failed:
        raise ComsolVerifyError(
            "COMSOL function verification failed: "
            f"{summary.failed} failed point(s); summary: {summary.summary_csv}"
        )
    return VerifyComsolSummary(
        plan=plan,
        passed=summary.passed,
        failed=summary.failed,
        expected_csv=plan.expected_csv,
        comsol_csv=plan.comsol_csv,
        summary_csv=plan.summary_csv,
        summary_json=plan.summary_json,
        java_path=plan.java_path,
        execution=execution,
    )


def _clear_previous_verify_outputs(plan: VerifyComsolPlan) -> None:
    for path in (
        plan.comsol_csv,
        plan.summary_csv,
        plan.summary_json,
    ):
        if path.exists() and path.is_file():
            path.unlink()


def write_verify_manifest(
    plan: VerifyComsolPlan,
    *,
    status: str,
    passed: int | None = None,
    failed: int | None = None,
    failure_reason: str | None = None,
) -> None:
    functions = []
    for function in plan.mapping.functions:
        function_points = [
            point
            for point in plan.points
            if point.function_name == function.name and point.function_tag == function.tag
        ]
        functions.append(
            {
                "name": function.name,
                "tag": function.tag,
                "file": function.file,
                "column": (
                    function.column
                    if function.column is not None
                    else (function_points[0].column if function_points else None)
                ),
                "nargs": function.nargs,
                "argunit": function.argunit,
                "fununit": function.fununit,
                "interp": function.interp,
                "extrap": function.extrap,
                "control_points": len(function_points),
            }
        )
    closure_items = []
    for quantity in _ordered_unique(
        point.quantity for point in plan.feature_points
    ):
        quantity_points = [
            point for point in plan.feature_points if point.quantity == quantity
        ]
        first = quantity_points[0]
        closure_items.append(
            {
                "quantity": quantity,
                "verification_scope": "model_feature_table_all_rows",
                "component": first.component,
                "physics": first.physics,
                "feature": first.feature,
                "x_property": first.x_property,
                "y_property": first.y_property,
                "file": first.file,
                "column": first.column,
                "control_points": len(quantity_points),
            }
        )
    payload: dict[str, Any] = {
        "status": status,
        "mapping": str(plan.mapping.path),
        "bundle": str(plan.mapping.bundle.path),
        "output_mph": str(plan.mapping.model.output_mph),
        "expected_csv": str(plan.expected_csv),
        "comsol_csv": str(plan.comsol_csv),
        "summary_csv": str(plan.summary_csv),
        "summary_json": str(plan.summary_json),
        "java_path": str(plan.java_path),
        "relative_tolerance": plan.relative_tolerance,
        "absolute_tolerance": plan.absolute_tolerance,
        "control_points": len(plan.points) + len(plan.feature_points),
        "function_control_points": len(plan.points),
        "feature_table_control_points": len(plan.feature_points),
        "functions": functions,
        "closure_items": closure_items,
    }
    if passed is not None:
        payload["passed"] = passed
    if failed is not None:
        payload["failed"] = failed
    if failure_reason is not None:
        payload["failure_reason"] = failure_reason
    _write_json(plan.manifest_path, payload)


def write_expected_values_csv(plan: VerifyComsolPlan) -> None:
    plan.expected_csv.parent.mkdir(parents=True, exist_ok=True)
    _write_csv(
        plan.expected_csv,
        (
            "function_name",
            "function_tag",
            "verification_kind",
            "quantity",
            "file",
            "column",
            "argument_column",
            "argument",
            "expected_value",
            "expected_unit",
            "argunit",
            "row_index",
            "reason",
            "component",
            "physics",
            "feature",
            "x_property",
            "y_property",
            "relative_tolerance",
            "absolute_tolerance",
        ),
        _expected_rows(plan),
    )


def _expected_rows(plan: VerifyComsolPlan) -> Iterable[dict[str, Any]]:
    for point in plan.points:
        yield {
            "function_name": point.function_name,
            "function_tag": point.function_tag,
            "verification_kind": "interpolation_function",
            "quantity": point.function_name,
            "file": point.file,
            "column": point.column,
            "argument_column": point.argument_column,
            "argument": point.argument,
            "expected_value": point.expected_value,
            "expected_unit": point.fununit,
            "argunit": point.argunit,
            "row_index": point.row_index,
            "reason": point.reason,
            "component": "",
            "physics": "",
            "feature": "",
            "x_property": "",
            "y_property": "",
            "relative_tolerance": plan.relative_tolerance,
            "absolute_tolerance": plan.absolute_tolerance,
        }
    for point in plan.feature_points:
        yield {
            "function_name": point.function_name,
            "function_tag": point.function_tag,
            "verification_kind": "feature_table",
            "quantity": point.quantity,
            "file": point.file,
            "column": point.column,
            "argument_column": point.argument_column,
            "argument": point.argument,
            "expected_value": point.expected_value,
            "expected_unit": point.fununit,
            "argunit": point.argunit,
            "row_index": point.row_index,
            "reason": point.reason,
            "component": point.component,
            "physics": point.physics,
            "feature": point.feature,
            "x_property": point.x_property,
            "y_property": point.y_property,
            "relative_tolerance": plan.relative_tolerance,
            "absolute_tolerance": plan.absolute_tolerance,
        }


def compare_function_values(plan: VerifyComsolPlan) -> VerifyComsolSummary:
    expected = _read_csv_dicts(plan.expected_csv)
    comsol_rows = _read_csv_dicts(plan.comsol_csv)
    if len(expected) != len(comsol_rows):
        raise ComsolVerifyError(
            f"COMSOL verification row count mismatch: expected {len(expected)}, "
            f"got {len(comsol_rows)}"
        )
    quantity_scales: dict[tuple[str, str], float] = {}
    for row in expected:
        key = (row.get("verification_kind", ""), row.get("quantity", ""))
        value = abs(_required_float(row, "expected_value"))
        quantity_scales[key] = max(quantity_scales.get(key, 0.0), value)
    summary_rows: list[dict[str, Any]] = []
    passed = 0
    failed = 0
    for index, (expected_row, comsol_row) in enumerate(zip(expected, comsol_rows)):
        _assert_same_identity(index, expected_row, comsol_row)
        expected_value = _required_float(expected_row, "expected_value")
        comsol_value = _required_float(comsol_row, "comsol_value")
        negative_or_nan = expected_value < 0.0 or comsol_value < 0.0
        absolute_error = abs(comsol_value - expected_value)
        denominator = abs(expected_value)
        relative_error = (
            0.0 if denominator == 0.0 and absolute_error == 0.0 else absolute_error / denominator
        )
        scale_key = (
            expected_row.get("verification_kind", ""),
            expected_row.get("quantity", ""),
        )
        scale = quantity_scales.get(scale_key, denominator)
        effective_absolute_tolerance = min(
            plan.absolute_tolerance,
            max(scale * plan.relative_tolerance, 1.0e-300),
        )
        ok = (
            not negative_or_nan
            and (
                absolute_error <= effective_absolute_tolerance
                or relative_error <= plan.relative_tolerance
            )
        )
        if ok:
            passed += 1
        else:
            failed += 1
        summary_rows.append(
            {
                "passed": int(ok),
                "function_name": expected_row["function_name"],
                "function_tag": expected_row["function_tag"],
                "verification_kind": expected_row["verification_kind"],
                "quantity": expected_row["quantity"],
                "verification_scope": (
                    "model_feature_table_all_rows"
                    if expected_row["verification_kind"] == "feature_table"
                    else "imported_interpolation_function_control_points"
                ),
                "file": expected_row["file"],
                "column": expected_row["column"],
                "argument": expected_row["argument"],
                "expected_value": expected_value,
                "comsol_value": comsol_value,
                "absolute_error": absolute_error,
                "relative_error": relative_error,
                "effective_absolute_tolerance": effective_absolute_tolerance,
                "negative_or_nan": int(negative_or_nan),
                "expected_unit": expected_row["expected_unit"],
                "argunit": expected_row["argunit"],
                "reason": expected_row["reason"],
            }
        )
    _write_csv(
        plan.summary_csv,
        (
            "passed",
            "function_name",
            "function_tag",
            "verification_kind",
            "quantity",
            "verification_scope",
            "file",
            "column",
            "argument",
            "expected_value",
            "comsol_value",
            "absolute_error",
            "relative_error",
            "effective_absolute_tolerance",
            "negative_or_nan",
            "expected_unit",
            "argunit",
            "reason",
        ),
        summary_rows,
    )
    max_absolute_error = max(
        (float(row["absolute_error"]) for row in summary_rows),
        default=0.0,
    )
    max_relative_error = max(
        (float(row["relative_error"]) for row in summary_rows),
        default=0.0,
    )
    _write_json(
        plan.summary_json,
        {
            "status": "passed" if failed == 0 else "failed",
            "passed": passed,
            "failed": failed,
            "control_points": len(summary_rows),
            "relative_tolerance": plan.relative_tolerance,
            "absolute_tolerance": plan.absolute_tolerance,
            "max_relative_error": max_relative_error,
            "max_absolute_error": max_absolute_error,
            "functions": sorted(
                {
                    str(row["function_tag"])
                    for row in summary_rows
                }
            ),
            "closure_items": _closure_item_summary(summary_rows),
            "rows": summary_rows,
        },
    )
    return VerifyComsolSummary(
        plan=plan,
        passed=passed,
        failed=failed,
        expected_csv=plan.expected_csv,
        comsol_csv=plan.comsol_csv,
        summary_csv=plan.summary_csv,
        summary_json=plan.summary_json,
        java_path=plan.java_path,
        execution=None,
    )


def _closure_item_summary(
    rows: Iterable[dict[str, Any]],
) -> list[dict[str, Any]]:
    feature_rows = [
        row for row in rows if row.get("verification_kind") == "feature_table"
    ]
    summaries: list[dict[str, Any]] = []
    for quantity in _ordered_unique(str(row["quantity"]) for row in feature_rows):
        quantity_rows = [
            row for row in feature_rows if str(row["quantity"]) == quantity
        ]
        summaries.append(
            {
                "quantity": quantity,
                "verification_scope": "model_feature_table_all_rows",
                "status": (
                    "passed"
                    if all(bool(row["passed"]) for row in quantity_rows)
                    else "failed"
                ),
                "control_points": len(quantity_rows),
                "max_relative_error": max(
                    (float(row["relative_error"]) for row in quantity_rows),
                    default=0.0,
                ),
                "max_absolute_error": max(
                    (float(row["absolute_error"]) for row in quantity_rows),
                    default=0.0,
                ),
            }
        )
    return summaries


def _extract_comsol_csv_from_stdout(stdout_path: Path, output_csv: Path) -> None:
    if not stdout_path.exists():
        return
    lines = stdout_path.read_text(encoding="utf-8", errors="replace").splitlines()
    try:
        start = lines.index("SWARM_COMSOL_VERIFY_CSV_BEGIN") + 1
        end = lines.index("SWARM_COMSOL_VERIFY_CSV_END", start)
    except ValueError:
        return
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_csv.write_text("\n".join(lines[start:end]) + "\n", encoding="utf-8")


def generate_verify_java_source(plan: VerifyComsolPlan) -> str:
    lines = [
        "import com.comsol.model.*;",
        "import com.comsol.model.util.*;",
        "import java.io.*;",
        "",
        f"public class {VERIFY_CLASS_NAME} {{",
        "  public static void main(String[] args) throws Exception {",
        "    Model model = ModelUtil.loadCopy("
        + java_string("swarmVerifyModel")
        + ", "
        + java_string(java_path(plan.mapping.model.output_mph))
        + ");",
        '    System.out.println("SWARM_COMSOL_VERIFY_CSV_BEGIN");',
        '    System.out.println("function_name,function_tag,argument,comsol_value");',
    ]
    for index, point in enumerate(plan.points):
        lines.extend(
            [
                f"    model.func({java_string(point.function_tag)});",
                "    writeValue("
                + java_string(point.function_name)
                + ", "
                + java_string(point.function_tag)
                + ", "
                + _java_double(point.argument)
                + ", tableValue(model, "
                + java_string(point.function_tag)
                + ", "
                + _java_double(point.argument)
                + "));",
            ]
        )
    for point in plan.feature_points:
        lines.extend(
            [
                "    writeValue("
                + java_string(point.function_name)
                + ", "
                + java_string(point.function_tag)
                + ", "
                + _java_double(point.argument)
                + ", featureTableValue(model, "
                + java_string(point.component)
                + ", "
                + java_string(point.physics)
                + ", "
                + java_string(point.feature)
                + ", "
                + java_string(point.x_property)
                + ", "
                + java_string(point.y_property)
                + ", "
                + str(point.row_index)
                + ", "
                + str(point.table_length)
                + ", "
                + _java_double(point.argument)
                + "));",
            ]
        )
    lines.extend(
        [
            '    System.out.println("SWARM_COMSOL_VERIFY_CSV_END");',
            "  }",
            "",
            "  private static double tableValue(Model model, String tag, double argument) {",
            "    String table = model.func(tag).getString(\"table\");",
            "    String[] rows = table.split(\";\");",
            "    double[] xs = new double[rows.length];",
            "    double[] ys = new double[rows.length];",
            "    int n = 0;",
            "    for (String row : rows) {",
            "      String trimmed = row.trim();",
            "      if (trimmed.length() == 0) continue;",
            "      String[] parts = trimmed.split(\",\");",
            "      if (parts.length < 2) continue;",
            "      xs[n] = Double.parseDouble(parts[0].trim());",
            "      ys[n] = Double.parseDouble(parts[1].trim());",
            "      n++;",
            "    }",
            "    if (n == 0) throw new RuntimeException(\"empty function table: \" + tag);",
            "    if (argument <= xs[0]) return ys[0];",
            "    for (int i = 1; i < n; i++) {",
            "      if (argument <= xs[i]) {",
            "        double span = xs[i] - xs[i - 1];",
            "        if (span == 0.0) return ys[i];",
            "        double t = (argument - xs[i - 1]) / span;",
            "        return ys[i - 1] + t * (ys[i] - ys[i - 1]);",
            "      }",
            "    }",
            "    return ys[n - 1];",
            "  }",
            "",
            "  private static double featureTableValue(Model model, String component,",
            "      String physics, String feature, String xProperty, String yProperty,",
            "      int index, int expectedLength, double expectedArgument) {",
            "    double[] xs = model.component(component).physics(physics).feature(feature)",
            "        .getDoubleArray(xProperty);",
            "    double[] ys = model.component(component).physics(physics).feature(feature)",
            "        .getDoubleArray(yProperty);",
            "    if (xs.length != ys.length || xs.length != expectedLength) {",
            "      throw new RuntimeException(\"feature table length mismatch: \"",
            "          + feature + \".\" + xProperty + \"/\" + yProperty);",
            "    }",
            "    if (index < 0 || index >= xs.length) {",
            "      throw new RuntimeException(\"feature table row missing: \"",
            "          + feature + \".\" + yProperty + \" index=\" + index);",
            "    }",
            "    double tolerance = Math.max(1e-15, Math.abs(expectedArgument) * 1e-12);",
            "    if (Math.abs(xs[index] - expectedArgument) > tolerance) {",
            "      throw new RuntimeException(\"feature table argument mismatch: \"",
            "          + feature + \".\" + xProperty + \" index=\" + index);",
            "    }",
            "    return ys[index];",
            "  }",
            "",
            "  private static void writeValue(String name, "
            + "String tag, double argument, double value) {",
            "    System.out.println(csv(name) + \",\" + csv(tag) + \",\" + "
            + "Double.toString(argument) + \",\" + Double.toString(value));",
            "  }",
            "",
            "  private static String csv(String value) {",
            "    return \"\\\"\" + value.replace(\"\\\"\", \"\\\"\\\"\") + \"\\\"\";",
            "  }",
            "}",
            "",
        ]
    )
    return "\n".join(lines)


def format_verify_plan(plan: VerifyComsolPlan) -> str:
    lines = [
        "COMSOL function verify dry-run",
        f"mapping: {plan.mapping.path}",
        f"output_mph: {plan.mapping.model.output_mph}",
        f"bundle: {plan.mapping.bundle.path}",
        f"output_dir: {plan.output_dir}",
        f"expected_csv: {plan.expected_csv}",
        f"manifest: {plan.manifest_path}",
        f"relative_tolerance: {plan.relative_tolerance:g}",
        f"absolute_tolerance: {plan.absolute_tolerance:g}",
        "control_points:",
    ]
    for point in plan.points:
        lines.append(
            "  - "
            f"{point.function_tag}({point.argument:.17g} {point.argunit}) "
            f"expected={point.expected_value:.17g} {point.fununit} "
            f"file={point.file} column={point.column} reason={point.reason}"
        )
    lines.append("feature_table_points:")
    for point in plan.feature_points:
        lines.append(
            "  - "
            f"{point.quantity}: "
            f"{point.physics}.{point.feature}.{point.y_property}"
            f"[{point.row_index}] expected={point.expected_value:.17g} "
            f"{point.fununit} at {point.argument:.17g} {point.argunit}"
        )
    lines.append("COMSOL command: not executed")
    return "\n".join(lines)


def format_verify_summary(summary: VerifyComsolSummary) -> str:
    return "\n".join(
        [
            "COMSOL function verification completed",
            f"passed: {summary.passed}",
            f"failed: {summary.failed}",
            f"expected_csv: {summary.expected_csv}",
            f"comsol_csv: {summary.comsol_csv}",
            f"summary_csv: {summary.summary_csv}",
            f"summary_json: {summary.summary_json}",
            f"manifest: {summary.plan.manifest_path}",
            f"java_path: {summary.java_path}",
        ]
    )


def _validate_verify_sources(mapping: ComsolModelMapping) -> None:
    if not mapping.bundle.path.exists() or not mapping.bundle.path.is_dir():
        raise ComsolVerifyError(
            f"COMSOL table bundle does not exist: {mapping.bundle.path}. "
            + _missing_bundle_hint(mapping.bundle.path)
        )
    manifest = _read_manifest(mapping.bundle.path / "manifest.json")
    tables = manifest.get("tables", {})
    if not isinstance(tables, dict):
        tables = {}
    for function in mapping.functions:
        if function.nargs != 1:
            raise ComsolVerifyError(
                f"COMSOL function verification supports nargs=1 only: "
                f"{function.name} has nargs={function.nargs}"
            )
        if not function.path.exists():
            raise ComsolVerifyError(
                f"mapped function CSV does not exist for {function.name}: "
                f"{function.path}"
            )
        header = _csv_header(function.path)
        column = _resolve_value_column(function, header)
        argument_column = _resolve_argument_column(function, header)
        if column not in header:
            raise ComsolVerifyError(
                f"mapped function column {column!r} not found in {function.file}"
            )
        table_meta = tables.get(function.file)
        if isinstance(table_meta, dict):
            columns = table_meta.get("columns")
            if isinstance(columns, list) and column not in columns:
                raise ComsolVerifyError(
                    f"mapped function column {column!r} is not listed in "
                    f"manifest table {function.file}"
                )
        units = _unit_metadata(table_meta, manifest)
        if column not in units:
            raise ComsolVerifyError(
                f"mapped function column {column!r} has no unit metadata in "
                f"bundle manifest for {function.file}"
            )
        if argument_column not in units:
            raise ComsolVerifyError(
                f"mapped function argument column {argument_column!r} has no unit "
                f"metadata in bundle manifest for {function.file}"
            )
        _require_unit_compatible(
            function.argunit,
            units[argument_column],
            column=argument_column,
            function_name=function.name,
        )
        _require_unit_compatible(
            function.fununit,
            units[column],
            column=column,
            function_name=function.name,
        )


def _control_points(mapping: ComsolModelMapping) -> list[FunctionControlPoint]:
    points: list[FunctionControlPoint] = []
    for function in mapping.functions:
        rows = _read_csv_dicts(function.path)
        if not rows:
            raise ComsolVerifyError(f"mapped function CSV is empty: {function.path}")
        header = tuple(rows[0])
        argument_column = _resolve_argument_column(function, header)
        value_column = _resolve_value_column(function, header)
        numeric_rows: list[tuple[int, float, float]] = []
        for index, row in enumerate(rows):
            argument = _float_or_none(row.get(argument_column))
            value = _float_or_none(row.get(value_column))
            if argument is None or value is None:
                continue
            if argument < 0.0:
                raise ComsolVerifyError(
                    f"mapped function {function.name} has negative argument "
                    f"in {function.file} row {index}"
                )
            if value < 0.0:
                raise ComsolVerifyError(
                    f"mapped function {function.name} has negative expected value "
                    f"in {function.file} row {index}"
                )
            numeric_rows.append((index, argument, value))
        if not numeric_rows:
            raise ComsolVerifyError(
                f"mapped function {function.name} has no finite control values "
                f"in {function.file}"
            )
        numeric_rows.sort(key=lambda item: item[1])
        for row_index, reason in _selected_indices(numeric_rows):
            source_index, argument, value = numeric_rows[row_index]
            points.append(
                FunctionControlPoint(
                    function_name=function.name,
                    function_tag=function.tag,
                    file=function.file,
                    column=value_column,
                    argument_column=argument_column,
                    argument=argument,
                    expected_value=value,
                    argunit=function.argunit,
                    fununit=function.fununit,
                    row_index=source_index,
                    reason=reason,
                )
            )
    return points


def _feature_table_points(
    mapping: ComsolModelMapping,
) -> list[FeatureTableControlPoint]:
    """Build all-row checks for the seven values injected into COMSOL."""

    closure = mapping.closure
    points: list[FeatureTableControlPoint] = []
    points.extend(
        _points_from_table(
            quantity="mean_energy",
            component=mapping.model.component,
            physics=closure.physics,
            feature=closure.feature,
            x_property="enrgXdata",
            y_property="enrgYdata",
            path=closure.mean_energy.path,
            file=closure.mean_energy.table,
            argument_column=closure.mean_energy.x_column,
            value_column=closure.mean_energy.y_column,
            argunit="Td",
            fununit="eV",
        )
    )
    transport_path = mapping.bundle.path / "transport_vs_mean_energy.csv"
    for quantity, x_property, y_property, column, unit in (
        (
            "reduced_mobility",
            "muNXdata",
            "muNYdata",
            "reduced_mobility_m2_V_s_m3",
            "1/(V*m*s)",
        ),
        (
            "reduced_longitudinal_diffusion",
            "deNXdata",
            "deNYdata",
            "reduced_diffusion_L_m2_s_m3",
            "1/(m*s)",
        ),
        (
            "reduced_energy_mobility",
            "mueNXdata",
            "mueNYdata",
            "reduced_electron_energy_mobility_m2_V_s_m3",
            "1/(V*m*s)",
        ),
        (
            "reduced_energy_diffusion",
            "denNXdata",
            "denNYdata",
            "reduced_electron_energy_diffusion_m2_s_m3",
            "1/(m*s)",
        ),
    ):
        points.extend(
            _points_from_table(
                quantity=quantity,
                component=mapping.model.component,
                physics=closure.physics,
                feature=closure.feature,
                x_property=x_property,
                y_property=y_property,
                path=transport_path,
                file=transport_path.name,
                argument_column="mean_energy_eV",
                value_column=column,
                argunit="eV",
                fununit=unit,
            )
        )
    for lookup in mapping.reaction_lookups:
        points.extend(
            _points_from_table(
                quantity=f"{lookup.process_type}_townsend",
                component=mapping.model.component,
                physics=lookup.physics,
                feature=lookup.feature,
                x_property="xtownratedata",
                y_property="ytownratedata",
                path=lookup.path,
                file=lookup.table,
                argument_column=lookup.x_column,
                value_column=lookup.y_column,
                argunit="eV",
                fununit="m^2",
                process_type=lookup.process_type,
            )
        )
    return points


def _points_from_table(
    *,
    quantity: str,
    component: str,
    physics: str,
    feature: str,
    x_property: str,
    y_property: str,
    path: Path,
    file: str,
    argument_column: str,
    value_column: str,
    argunit: str,
    fununit: str,
    process_type: str | None = None,
) -> list[FeatureTableControlPoint]:
    if not path.exists():
        raise ComsolVerifyError(
            f"required closure verification table does not exist: {path}"
        )
    values: dict[float, float] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        for source_index, row in enumerate(csv.DictReader(stream)):
            if process_type is not None and row.get("process_type") != process_type:
                continue
            argument = _float_or_none(row.get(argument_column))
            value = _float_or_none(row.get(value_column))
            if argument is None or value is None:
                raise ComsolVerifyError(
                    f"required closure quantity {quantity} has a nonnumeric row "
                    f"in {file} row {source_index}"
                )
            if argument < 0.0 or value < 0.0:
                raise ComsolVerifyError(
                    f"required closure quantity {quantity} has a negative value "
                    f"in {file} row {source_index}"
                )
            values[argument] = value
    if len(values) < 2:
        raise ComsolVerifyError(
            f"required closure quantity {quantity} needs at least two rows in {file}"
        )
    ordered = sorted(values.items())
    table_length = len(ordered)
    return [
        FeatureTableControlPoint(
            quantity=quantity,
            component=component,
            physics=physics,
            feature=feature,
            x_property=x_property,
            y_property=y_property,
            file=file,
            column=value_column,
            argument_column=argument_column,
            argument=argument,
            expected_value=value,
            argunit=argunit,
            fununit=fununit,
            row_index=index,
            table_length=table_length,
        )
        for index, (argument, value) in enumerate(ordered)
    ]


def _ordered_unique(values: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(values))


def _selected_indices(rows: list[tuple[int, float, float]]) -> list[tuple[int, str]]:
    raw = [(0, "min"), (len(rows) // 2, "mid"), (len(rows) - 1, "max")]
    max_abs = max(abs(value) for _index, _argument, value in rows)
    if max_abs > 0.0:
        threshold = max_abs * 1.0e-2
        for index, (_source_index, _argument, value) in enumerate(rows):
            if abs(value) >= threshold:
                raw.append((index, "rise"))
                break
    seen: set[int] = set()
    selected: list[tuple[int, str]] = []
    for index, reason in raw:
        if index not in seen:
            seen.add(index)
            selected.append((index, reason))
    return selected


def _resolve_argument_column(function: FunctionMapping, header: Iterable[str]) -> str:
    columns = tuple(header)
    normalized = function.argunit.replace(" ", "")
    if normalized in {"V*m^2", "V*m2"} and "E_over_N_V_m2" in columns:
        return "E_over_N_V_m2"
    if normalized == "Td" and "E_over_N_Td" in columns:
        return "E_over_N_Td"
    if normalized == "eV" and "mean_energy_eV" in columns:
        return "mean_energy_eV"
    if "E_over_N_V_m2" in columns:
        return "E_over_N_V_m2"
    if not columns:
        raise ComsolVerifyError(f"mapped function CSV is empty: {function.file}")
    return columns[0]


def _resolve_value_column(function: FunctionMapping, header: Iterable[str]) -> str:
    columns = tuple(header)
    if function.column is not None:
        return function.column
    argument_column = _resolve_argument_column(function, columns)
    metadata_columns = {
        argument_column,
        "E_over_N_Td",
        "E_over_N_V_m2",
        "electron_energy_eV",
        "energy_width_eV",
    }
    if argument_column == "mean_energy_eV":
        metadata_columns.add("mean_energy_eV")
    candidates = [column for column in columns if column not in metadata_columns]
    if len(candidates) == 1:
        return candidates[0]
    raise ComsolVerifyError(
        f"functions.{function.name}.column is required because {function.file} "
        "has multiple value columns"
    )


def _assert_same_identity(
    index: int,
    expected_row: dict[str, str],
    comsol_row: dict[str, str],
) -> None:
    for column in ("function_name", "function_tag"):
        if expected_row.get(column) != comsol_row.get(column):
            raise ComsolVerifyError(
                f"COMSOL verification row {index} does not match expected "
                f"{column}: {expected_row.get(column)!r} != "
                f"{comsol_row.get(column)!r}"
            )
    expected_argument = _required_float(expected_row, "argument")
    comsol_argument = _required_float(comsol_row, "argument")
    if not math.isclose(expected_argument, comsol_argument, rel_tol=1.0e-12, abs_tol=0.0):
        raise ComsolVerifyError(
            f"COMSOL verification row {index} does not match expected argument: "
            f"{expected_argument!r} != {comsol_argument!r}"
        )


def _read_manifest(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise ComsolVerifyError(
            f"COMSOL table bundle manifest.json does not exist: {path}. "
            + _missing_bundle_hint(path.parent)
        )
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:  # pragma: no cover - exact JSON error is not useful here.
        raise ComsolVerifyError(f"invalid bundle manifest JSON: {path}") from exc
    if not isinstance(data, dict):
        raise ComsolVerifyError(f"bundle manifest must be a mapping: {path}")
    return data


def _missing_bundle_hint(path: Path) -> str:
    return (
        "Generate the Ar-only bundle first, for example: "
        "swarm-workflow sweep examples/workflow_argon_comsol.yaml; "
        "swarm-workflow aggregate outputs/argon_comsol/swarm.sqlite --output "
        "outputs/argon_comsol/aggregate; "
        "swarm-workflow build-tables outputs/argon_comsol/swarm.sqlite --output "
        "outputs/argon_comsol/tables --source two_term; "
        "swarm-workflow export-comsol outputs/argon_comsol/tables/mixture_0000 "
        "--output "
        f"{path.as_posix()}"
    )


def _csv_header(path: Path) -> tuple[str, ...]:
    with path.open("r", encoding="utf-8", newline="") as fp:
        reader = csv.reader(fp)
        try:
            return tuple(next(reader))
        except StopIteration as exc:
            raise ComsolVerifyError(f"mapped function CSV is empty: {path}") from exc


def _read_csv_dicts(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fp:
        return list(csv.DictReader(fp))


def _unit_metadata(table_meta: object, manifest: dict[str, Any]) -> dict[str, str]:
    units: dict[str, str] = {}
    table_units: object = {}
    if isinstance(table_meta, dict):
        table_units = table_meta.get("units", {})
    manifest_units = manifest.get("units", {})
    if isinstance(manifest_units, dict):
        units.update({str(key): str(value) for key, value in manifest_units.items()})
    if isinstance(table_units, dict):
        units.update({str(key): str(value) for key, value in table_units.items()})
    return units


def _require_unit_compatible(
    mapping_unit: str,
    manifest_unit: str,
    *,
    column: str,
    function_name: str,
) -> None:
    mapping_normalized = _normalize_unit(mapping_unit)
    manifest_normalized = _normalize_unit(manifest_unit)
    if mapping_normalized == manifest_normalized:
        return
    aliases = {_normalize_unit(unit) for unit in _UNIT_ALIASES.get(column, ())}
    if mapping_normalized in aliases and manifest_normalized in aliases:
        return
    raise ComsolVerifyError(
        f"unit mismatch for functions.{function_name} column {column}: "
        f"mapping has {mapping_unit!r}, bundle manifest has {manifest_unit!r}"
    )


def _normalize_unit(value: str) -> str:
    return (
        value.lower()
        .replace(" ", "")
        .replace("*", "")
        .replace("·", "")
        .replace("^", "")
        .replace("(", "")
        .replace(")", "")
    )


_UNIT_ALIASES = {
    "E_over_N_V_m2": {"V*m^2", "V m^2", "V*m2", "V m2"},
    "E_over_N_Td": {"Td"},
    "mean_energy_eV": {"eV"},
    "reduced_mobility_m2_V_s_m3": {
        "1/(V*m*s)",
        "m^2/(V s) m^3",
        "m^2/(V*s)*1/m^3",
    },
    "reduced_diffusion_L_m2_s_m3": {"1/(m*s)", "m^2/s m^3"},
    "reduced_diffusion_T_m2_s_m3": {"1/(m*s)", "m^2/s m^3"},
    "reduced_electron_energy_mobility_m2_V_s_m3": {
        "1/(V*m*s)",
        "m^2/(V s) m^3",
        "m^2/(V*s)*1/m^3",
    },
    "reduced_electron_energy_diffusion_m2_s_m3": {
        "1/(m*s)",
        "m^2/s m^3",
    },
    "rate_coefficient_m3_s": {"m^3/s", "m^3 s^-1"},
    "mixture_weighted_rate_m3_s": {"m^3/s", "m^3 s^-1"},
    "reduced_townsend_m2": {"m^2"},
    "mixture_weighted_reduced_townsend_m2": {"m^2"},
    "effective_townsend_m2": {"m^2"},
    "energy_loss_rate_coefficient_eV_m3_s": {"eV*m^3/s", "eV m^3/s"},
}


def _required_float(row: dict[str, str], key: str) -> float:
    value = _float_or_none(row.get(key))
    if value is None:
        raise ComsolVerifyError(f"missing finite {key} in COMSOL verification CSV")
    return value


def _float_or_none(value: object) -> float | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _java_double(value: float) -> str:
    return format(float(value), ".15g")
