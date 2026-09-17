"""COMSOL Java rendering and sentinel-log extraction for EEDF audits."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ...._io import write_csv
from .contracts import (
    _QUERY_COLUMNS,
    _VALUE_COLUMNS,
    ComsolEedfAuditError,
    ComsolEedfAuditPlan,
)
from .io import _finite_float, _read_rows, _validate_audit_inputs


def render_comsol_eedf_audit_java(
    *,
    class_name: str,
    model_path: Path,
    component: str,
    physics: str,
    function_tag: str,
    audit_rows: list[dict[str, Any]],
    table_sha256: str,
    query_sha256: str,
) -> str:
    """Render a COMSOL program that inspects and evaluates the saved Function."""

    values = {
        "model": _java_string(model_path.as_posix()),
        "component": _java_string(component),
        "physics": _java_string(physics),
        "function": _java_string(function_tag),
    }
    point_methods: list[str] = []
    point_calls: list[str] = []
    for batch_index, start in enumerate(range(0, len(audit_rows), 400)):
        points = ",\n".join(
            "      new double[] {"
            f"{float(row['electron_energy_eV']):.17g}, "
            f"{float(row['mean_energy_eV']):.17g}"
            "}"
            for row in audit_rows[start : start + 400]
        )
        point_methods.append(
            f"  private static double[][] points{batch_index}() {{\n"
            "    return new double[][] {\n"
            f"{points}\n"
            "    };\n"
            "  }"
        )
        point_calls.append(f"points{batch_index}()")
    methods = "\n\n".join(point_methods)
    calls = ", ".join(point_calls)
    return f"""import com.comsol.model.FunctionFeature;
import com.comsol.model.Model;
import com.comsol.model.util.ModelUtil;

public class {class_name} {{
  private static String joined(String[] values) {{
    return values == null ? "" : String.join("|", values);
  }}

  private static void property(String key, String value) {{
    System.out.println("SWARM_EEDF_CONTRACT\\t" + key + "\\t"
        + (value == null ? "" : value));
  }}

{methods}

  public static void main(String[] ignored) throws Exception {{
    String componentTag = {values["component"]};
    String physicsTag = {values["physics"]};
    String functionTag = {values["function"]};
    Model model = ModelUtil.loadCopy("swarmEedfAudit", {values["model"]});
    property("audit_table_sha256", {_java_string(table_sha256)});
    property("audit_queries_sha256", {_java_string(query_sha256)});
    FunctionFeature function = model.component(componentTag).func(functionTag);
    property("function_names", joined(function.functionNames()));
    property("source", function.getString("source"));
    property("struct", function.getString("struct"));
    property("interp", function.getString("interp"));
    property("extrap", function.getString("extrap"));
    property("fununit", joined(function.getStringArray("fununit")));
    property("argunit", joined(function.getStringArray("argunit")));
    property("imported_name", function.getString("importedname"));
    property("filename", function.getString("filename"));
    String[] problems = function.problem().tags();
    property("problem_count", Integer.toString(problems == null ? 0 : problems.length));
    property("problems", joined(problems));
    property(
        "eedf_selection",
        model.component(componentTag).physics(physicsTag)
            .prop("EEDFSettings").getString("eedf"));

    double[][][] batches = new double[][][] {{{calls}}};
    int index = 0;
    for (double[][] points : batches) {{
      for (double[] point : points) {{
        double energy = point[0];
        double mean = point[1];
        String expression = componentTag + "." + functionTag
            + "(" + Double.toString(energy) + "[eV],"
            + Double.toString(mean) + "[eV])";
        model.param().set("sw_eedf_audit_value", expression);
        double value = model.param().evaluate("sw_eedf_audit_value");
        if (!Double.isFinite(value)) throw new IllegalStateException("nonfinite EEDF value");
        System.out.println(
            "SWARM_EEDF_VALUE\\t" + Integer.toString(index) + "\\t"
                + Double.toString(value));
        index++;
      }}
    }}
    ModelUtil.remove("swarmEedfAudit");
  }}
}}
"""


def extract_comsol_eedf_audit_log(
    plan: ComsolEedfAuditPlan,
    batch_log_path: str | Path,
) -> None:
    """Extract sentinel records without granting Java filesystem write access."""

    _validate_audit_inputs(plan)
    contract: dict[str, str] = {}
    values: dict[int, float] = {}
    for line in (
        Path(batch_log_path)
        .read_text(encoding="utf-8-sig", errors="replace")
        .splitlines()
    ):
        if line.startswith("SWARM_EEDF_CONTRACT\t"):
            unused_marker, key, value = line.split("\t", 2)
            if key in contract:
                raise ComsolEedfAuditError(f"duplicate COMSOL contract key: {key}")
            contract[key] = value
        elif line.startswith("SWARM_EEDF_VALUE\t"):
            unused_marker, raw_index, raw_value = line.split("\t", 2)
            index = int(raw_index)
            if index in values:
                raise ComsolEedfAuditError(f"duplicate COMSOL audit index: {index}")
            values[index] = _finite_float(raw_value)
    queries = _read_rows(plan.query_path, _QUERY_COLUMNS)
    if set(values) != set(range(len(queries))) or not contract:
        raise ComsolEedfAuditError("incomplete COMSOL EEDF audit batch log")
    if (
        contract.get("audit_table_sha256") != plan.table_sha256
        or contract.get("audit_queries_sha256") != plan.query_sha256
    ):
        raise ComsolEedfAuditError(
            "COMSOL audit log belongs to a different input or query plan"
        )
    plan.contract_path.write_text(
        "".join(f"{key}\t{value}\n" for key, value in contract.items()),
        encoding="utf-8",
    )
    write_csv(
        plan.values_path,
        _VALUE_COLUMNS,
        (
            {**query, "comsol_value": values[index]}
            for index, query in enumerate(queries)
        ),
    )


def _java_string(value: str) -> str:
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'
