"""Generate Java source for applying mapped CSV tables to COMSOL models."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Iterable

from ...java import java_path, java_string
from .lookup import read_lookup
from .config import (
    ClosureMapping,
    ComsolModelMapping,
    FunctionMapping,
    ReactionLookupMapping,
)


DEFAULT_CLASS_NAME = "SwarmComsolApply"
VOLTAGE_LOG_PREFIX = "SWARM_COMSOL_VOLTAGE_V="
PRESSURE_LOG_PREFIX = "SWARM_COMSOL_PRESSURE_PA="
GAS_TEMPERATURE_LOG_PREFIX = "SWARM_COMSOL_GAS_TEMPERATURE_K="
MESH_ELEMENTS_LOG_PREFIX = "SWARM_COMSOL_MESH_ELEMENTS="
MEAN_ENERGY_MODEL_LOG_PREFIX = "SWARM_COMSOL_MEAN_ENERGY_MODEL="


def generate_apply_java_source(
    mapping: ComsolModelMapping,
    *,
    class_name: str = DEFAULT_CLASS_NAME,
    pressure_Pa: float | None = None,
    gas_temperature_K: float | None = None,
    mesh_elements: int | None = None,
) -> str:
    """Return standalone Java source that updates mapped interpolation functions."""

    _validate_java_class_name(class_name)
    lines: list[str] = [
        "import com.comsol.model.*;",
        "import com.comsol.model.util.*;",
        "",
        f"public class {class_name} {{",
        "  public static void main(String[] args) throws Exception {",
        "    Model model = ModelUtil.loadCopy("
        + java_string("swarmModel")
        + ", "
        + java_string(java_path(mapping.model.input_mph))
        + ");",
    ]
    if any(
        value is not None for value in (pressure_Pa, gas_temperature_K, mesh_elements)
    ):
        lines.append("    configureModelConditions(model);")
    lines.extend(
        [
            "    apply(model);",
            "    model.save(" + java_string(java_path(mapping.model.output_mph)) + ");",
            "  }",
            "",
            "  public static void apply(Model model) {",
        ]
    )
    for function in mapping.functions:
        lines.extend(_function_lines(function))
    lines.extend(_closure_lines(mapping, mapping.closure))
    for lookup in mapping.reaction_lookups:
        lines.extend(_reaction_lookup_lines(mapping, lookup))
    condition_method = _model_condition_method_lines(
        mapping,
        pressure_Pa=pressure_Pa,
        gas_temperature_K=gas_temperature_K,
        mesh_elements=mesh_elements,
    )
    lines.extend(
        [
            "  }",
            "",
            *condition_method,
            "  private static void ensureInterpolation(Model model, String tag) {",
            "    try {",
            "      model.func(tag);",
            "    } catch (Exception ex) {",
            '      model.func().create(tag, "Interpolation");',
            "    }",
            "  }",
            "}",
            "",
        ]
    )
    return "\n".join(lines)


def _model_condition_method_lines(
    mapping: ComsolModelMapping,
    *,
    pressure_Pa: float | None,
    gas_temperature_K: float | None,
    mesh_elements: int | None,
) -> list[str]:
    if all(value is None for value in (pressure_Pa, gas_temperature_K, mesh_elements)):
        return []
    component = java_string(mapping.model.component)
    physics = java_string(mapping.model.physics)
    plasma_feature = java_string(mapping.closure.feature)
    lines = [
        "  private static void configureModelConditions(Model model) {",
    ]
    if pressure_Pa is not None:
        pressure = _finite_positive_java_number(pressure_Pa, "pressure_Pa")
        lines.extend(
            [
                f'    model.param().set("p0", "{pressure}[Pa]");',
                (
                    f"    model.component({component}).physics({physics})"
                    f'.feature({plasma_feature}).set("pA_src", "userdef");'
                ),
                (
                    f"    model.component({component}).physics({physics})"
                    f'.feature({plasma_feature}).set("pA", "p0");'
                ),
                f'    System.out.println("{PRESSURE_LOG_PREFIX}{pressure}");',
            ]
        )
    if gas_temperature_K is not None:
        temperature = _finite_positive_java_number(
            gas_temperature_K,
            "gas_temperature_K",
        )
        lines.extend(
            [
                (
                    f"    model.component({component}).physics({physics})"
                    f'.feature({plasma_feature}).set("T_src", "userdef");'
                ),
                (
                    f"    model.component({component}).physics({physics})"
                    f'.feature({plasma_feature}).set("T", '
                    f'"{temperature}[K]");'
                ),
                (
                    "    String configuredGasTemperature = "
                    f"model.component({component}).physics({physics})"
                    f'.feature({plasma_feature}).getString("T");'
                ),
                (
                    f'    System.out.println("{GAS_TEMPERATURE_LOG_PREFIX}"'
                    " + configuredGasTemperature);"
                ),
            ]
        )
    if mesh_elements is not None:
        if mesh_elements not in (200, 400):
            raise ValueError("mesh_elements must be 200 or 400")
        lines.extend(
            [
                (
                    f'    model.component({component}).mesh("mesh1")'
                    '.feature("edg1").feature("dis1")'
                    f'.set("elemcount", {mesh_elements});'
                ),
                f'    model.component({component}).mesh("mesh1").run();',
                (
                    "    String configuredMeshElements = "
                    f'model.component({component}).mesh("mesh1")'
                    '.feature("edg1").feature("dis1")'
                    '.getString("elemcount");'
                ),
                (
                    f'    System.out.println("{MESH_ELEMENTS_LOG_PREFIX}"'
                    " + configuredMeshElements);"
                ),
            ]
        )
    lines.extend(["  }", ""])
    return lines


def _finite_positive_java_number(value: float, name: str) -> str:
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return format(number, ".17g")


def generate_continuation_java_source(
    mapping: ComsolModelMapping,
    *,
    class_name: str,
) -> str:
    """Run the target voltage directly, with continuation only as fallback."""

    _validate_java_class_name(class_name)
    component = java_string(mapping.model.component)
    physics = java_string(mapping.model.physics)
    feature = java_string(mapping.run.voltage_feature)
    voltage_property = java_string(mapping.run.voltage_property)
    voltage_parameter = java_string(mapping.run.voltage_parameter)
    voltage_expression = java_string(mapping.run.voltage_parameter)
    lines = [
        "import com.comsol.model.*;",
        "import com.comsol.model.util.*;",
        "",
        f"public class {class_name} {{",
        "  public static void main(String[] args) throws Exception {",
        "    Model model = ModelUtil.loadCopy("
        + java_string("swarmPositiveColumn")
        + ", "
        + java_string(java_path(mapping.model.output_mph))
        + ");",
        "    configureVoltage(model);",
        "    try {",
        "      runAtVoltage(model, "
        + java_string(f"{mapping.run.voltages_V[-1]:g}[V]")
        + ");",
        "    } catch (Exception directFailure) {",
        '      System.out.println("Direct target-voltage solve failed; using continuation.");',
        '      ModelUtil.remove("swarmPositiveColumn");',
        "      model = ModelUtil.loadCopy("
        + java_string("swarmPositiveColumn")
        + ", "
        + java_string(java_path(mapping.model.output_mph))
        + ");",
        "      configureVoltage(model);",
    ]
    for voltage in mapping.run.voltages_V:
        value = java_string(f"{voltage:g}[V]")
        lines.append(f"      runAtVoltage(model, {value});")
    lines.extend(
        [
            "    }",
            "  }",
            "",
            "  private static void configureVoltage(Model model) {",
            f"    model.component({component}).physics({physics}).feature({feature})"
            f".set({voltage_property}, {voltage_expression});",
            "  }",
            "",
            (
                "  private static void runAtVoltage(Model model, String voltage) "
                "throws Exception {"
            ),
            (
                "    System.out.println("
                + java_string(VOLTAGE_LOG_PREFIX)
                + " + voltage);"
            ),
            f"    model.param().set({voltage_parameter}, voltage);",
            f"    model.study({java_string(mapping.model.study)}).run();",
            "    model.save(" + java_string(java_path(mapping.model.output_mph)) + ");",
            "  }",
            "}",
            "",
        ]
    )
    return "\n".join(lines)


def _function_lines(function: FunctionMapping) -> list[str]:
    tag = java_string(function.tag)
    return [
        f"    ensureInterpolation(model, {tag});",
        f'    model.func({tag}).set("source", "file");',
        f'    model.func({tag}).set("filename", '
        + java_string(java_path(function.path))
        + ");",
        f'    model.func({tag}).set("nargs", {function.nargs});',
        f'    model.func({tag}).set("argunit", {java_string(function.argunit)});',
        f'    model.func({tag}).set("fununit", {java_string(function.fununit)});',
        f'    model.func({tag}).set("interp", {java_string(function.interp)});',
        f'    model.func({tag}).set("extrap", {java_string(function.extrap)});',
        f"    model.func({tag}).importData();",
        "",
    ]


def _closure_lines(
    mapping: ComsolModelMapping,
    closure: ClosureMapping,
) -> list[str]:
    component = java_string(mapping.model.component)
    physics = java_string(closure.physics)
    feature = java_string(closure.feature)
    data = read_lookup(
        closure.mean_energy.path,
        closure.mean_energy.x_column,
        closure.mean_energy.y_column,
    )
    x_values, y_values = data.x, data.y
    transport = mapping.bundle.path / "transport_vs_mean_energy.csv"
    transport_tables = (
        ("muNXdata", "muNYdata", "reduced_mobility_m2_V_s_m3"),
        ("deNXdata", "deNYdata", "reduced_diffusion_L_m2_s_m3"),
        (
            "mueNXdata",
            "mueNYdata",
            "reduced_electron_energy_mobility_m2_V_s_m3",
        ),
        (
            "denNXdata",
            "denNYdata",
            "reduced_electron_energy_diffusion_m2_s_m3",
        ),
    )
    target = f"model.component({component}).physics({physics}).feature({feature})"
    formulation = closure.mean_energy_formulation
    physics_target = f"model.component({component}).physics({physics})"
    formulation_target = (
        f"{physics_target}.prop({java_string(formulation.property_group)})"
    )
    expected_formulation = java_string(formulation.comsol_value)
    lines = [
        (
            "    // Explicit mean-energy formulation plus written Swarm table "
            f"audit ({formulation.mode})."
        ),
        (
            f"    {formulation_target}.set("
            f"{java_string(formulation.property)}, {expected_formulation});"
        ),
        (
            "    String configuredMeanEnergyModel = "
            f"{formulation_target}.getString({java_string(formulation.property)});"
        ),
        (f"    if (!{expected_formulation}.equals(configuredMeanEnergyModel)) {{"),
        (
            "      throw new IllegalStateException("
            '"COMSOL mean-energy formulation read-back mismatch: expected " + '
            f'{expected_formulation} + ", got " + configuredMeanEnergyModel);'
        ),
        "    }",
        (
            f'    System.out.println("{MEAN_ENERGY_MODEL_LOG_PREFIX}"'
            " + configuredMeanEnergyModel);"
        ),
        f'    {target}.set("SpecifyElectronDensityAndEnergy", "UseLookupTables");',
        f'    {target}.set("SpecifyMeanElectronEnergy", "MeanEnergyTable");',
        f'    {target}.set("enrgXdata", new double[]{{{_java_double_array(x_values)}}});',
        f'    {target}.set("enrgYdata", new double[]{{{_java_double_array(y_values)}}});',
    ]
    for x_property, y_property, column in transport_tables:
        data = read_lookup(transport, "mean_energy_eV", column)
        x_transport, y_transport = data.x, data.y
        lines.extend(
            [
                f"    {target}.set({java_string(x_property)}, "
                f"new double[]{{{_java_double_array(x_transport)}}});",
                f"    {target}.set({java_string(y_property)}, "
                f"new double[]{{{_java_double_array(y_transport)}}});",
            ]
        )
    lines.append("")
    return lines


def _reaction_lookup_lines(
    mapping: ComsolModelMapping,
    lookup: ReactionLookupMapping,
) -> list[str]:
    if lookup.form != "townsend":
        raise ValueError(f"unsupported reaction lookup form: {lookup.form}")
    data = read_lookup(
        lookup.path, lookup.x_column, lookup.y_column, reaction=lookup.identity
    )
    x_values, y_values = data.x, data.y
    component = java_string(mapping.model.component)
    physics = java_string(lookup.physics)
    feature = java_string(lookup.feature)
    return [
        f"    // Mapping-driven Townsend lookup for {lookup.name} "
        f"({lookup.source_model}).",
        f"    model.component({component}).physics({physics}).feature({feature})"
        '.set("SpecifyReactionUsing", "UseLookupTable");',
        f"    model.component({component}).physics({physics}).feature({feature})"
        '.set("RateConstantForm", "UseTownsend");',
        f"    model.component({component}).physics({physics}).feature({feature})"
        '.set("UseTownsendTwoTermsBoltzmann", false);',
        f"    model.component({component}).physics({physics}).feature({feature})"
        f'.set("xtownratedata", new double[]{{{_java_double_array(x_values)}}});',
        f"    model.component({component}).physics({physics}).feature({feature})"
        f'.set("ytownratedata", new double[]{{{_java_double_array(y_values)}}});',
        "",
    ]


def _java_double_array(values: Iterable[float]) -> str:
    return ", ".join(f"{float(value):.17e}" for value in values)


def generate_field_export_java_source(
    mapping: ComsolModelMapping,
    *,
    class_name: str,
    output_csv: Path,
) -> str:
    """Export selected saved solutions without running a COMSOL study."""
    _validate_java_class_name(class_name)
    results = mapping.results
    if results is None:
        raise ValueError("field export requires mapping results")
    expressions = ", ".join(java_string(field.expression) for field in results.fields)
    units = ", ".join(java_string(field.unit) for field in results.fields)
    header = (
        "dataset,solution_index,point_index,time_s,parameter_name,parameter_unit,parameter_value,"
        + ",".join(field.name for field in results.fields)
    )
    first = "0" if results.solutions == "all" else "solutionNumbers.length - 1"
    return f"""import com.comsol.model.*;
import com.comsol.model.util.*;
import java.nio.file.*;
import java.nio.charset.StandardCharsets;
import java.io.BufferedWriter;
import java.util.Arrays;

public class {class_name} {{
  private static String csv(String value) {{
    if (value == null) return "";
    return "\\\"" + value.replace("\\\"", "\\\"\\\"") + "\\\"";
  }}
  public static void main(String[] args) throws Exception {{
    Model model = ModelUtil.loadCopy("swarmPositiveColumnResults", {java_string(java_path(mapping.model.output_mph))});
    String dataset = {java_string(results.dataset)};
    String solution = model.result().dataset(dataset).getString("solution");
    SolutionInfo info = model.sol(solution).getSolutioninfo();
    int[] outer = info.getOuterSolnum();
    if (outer.length != 1 || outer[0] != 1) {{
      throw new IllegalStateException("positive-column export does not support outer parameter sweeps");
    }}
    int[] solutionNumbers = info.getSolnum(1, true);
    String[] parameterNames = model.sol(solution).getPNames();
    if (parameterNames.length > 1) {{
      throw new IllegalStateException("positive-column export supports at most one inner solution parameter");
    }}
    String parameterName = parameterNames.length == 0 ? "" : parameterNames[0];
    String parameterUnit = parameterName.isEmpty() ? "" : info.getUnit(parameterName);
    double[] parameterValues = parameterName.isEmpty() ? new double[0] : model.sol(solution).getPVals();
    if (solutionNumbers.length == 0 || (!parameterName.isEmpty() && parameterValues.length != solutionNumbers.length)) {{
      throw new IllegalStateException("saved solution parameters do not match solution numbers");
    }}
    boolean transientSolution = parameterName.equals("t");
    String[] expressions = new String[]{{{expressions}}};
    String[] units = new String[]{{{units}}};
    int fieldCount = expressions.length;
    if (transientSolution) {{
      expressions = Arrays.copyOf(expressions, fieldCount + 1);
      units = Arrays.copyOf(units, fieldCount + 1);
      expressions[fieldCount] = "t";
      units[fieldCount] = "s";
    }}
    String evalTag = "swpc_fields";
    try {{ model.result().numerical().remove(evalTag); }} catch (Exception ignored) {{}}
    model.result().numerical().create(evalTag, "Eval");
    model.result().numerical(evalTag).set("data", dataset);
    model.result().numerical(evalTag).set("matherr", "on");
    model.result().numerical(evalTag).set("expr", expressions);
    model.result().numerical(evalTag).set("unit", units);
    double[][][] values = model.result().numerical(evalTag).getData();
    if (model.result().numerical(evalTag).isComplex()) {{
      throw new IllegalStateException("positive-column fields must be real; refusing to discard imaginary values");
    }}
    String[] evaluatedUnits = model.result().numerical(evalTag).getStringArray("unit");
    if (evaluatedUnits.length != units.length) {{
      throw new IllegalStateException("evaluated units do not match the requested fields");
    }}
    for (int field = 0; field < units.length; field++) {{
      if (!units[field].replaceAll("\\\\s+", "").equals(evaluatedUnits[field].replaceAll("\\\\s+", ""))) {{
        throw new IllegalStateException("evaluated unit differs for " + expressions[field]
            + ": requested=" + units[field] + ", evaluated=" + evaluatedUnits[field]);
      }}
    }}
    int points = model.result().numerical(evalTag).getNData();
    if (values.length != expressions.length || points == 0) {{
      throw new IllegalStateException("empty or inconsistent field export");
    }}
    for (double[][] expression : values) {{
      if (expression.length != solutionNumbers.length) {{
        throw new IllegalStateException("field export omitted saved solutions");
      }}
      for (double[] data : expression) {{
        if (data.length != points) throw new IllegalStateException("inconsistent spatial field grid");
        for (double value : data) {{
          if (!Double.isFinite(value)) throw new IllegalStateException("nonfinite exported field");
        }}
      }}
    }}
    Path output = Paths.get({java_string(java_path(output_csv))});
    Files.createDirectories(output.toAbsolutePath().getParent());
    Path temporary = Files.createTempFile(output.toAbsolutePath().getParent(), ".fields-", ".csv");
    try {{
      try (BufferedWriter writer = Files.newBufferedWriter(temporary, StandardCharsets.UTF_8)) {{
        writer.write({java_string(header)});
        writer.newLine();
        for (int s = {first}; s < solutionNumbers.length; s++) {{
          for (int point = 0; point < points; point++) {{
            String time = transientSolution ? Double.toString(values[fieldCount][s][point]) : "";
            String parameter = parameterName.isEmpty() ? "" : Double.toString(parameterValues[s]);
            writer.write(csv(dataset) + "," + solutionNumbers[s] + "," + point + "," + time
                + "," + csv(parameterName) + "," + csv(parameterUnit) + "," + parameter);
            for (int field = 0; field < fieldCount; field++) {{
              writer.write("," + Double.toString(values[field][s][point]));
            }}
            writer.newLine();
          }}
        }}
      }}
      Files.move(temporary, output, StandardCopyOption.REPLACE_EXISTING, StandardCopyOption.ATOMIC_MOVE);
    }} finally {{
      Files.deleteIfExists(temporary);
    }}
  }}
}}
"""


def _validate_java_class_name(value: str) -> None:
    if not value or not value.replace("_", "a").isalnum() or value[0].isdigit():
        raise ValueError(f"invalid Java class name: {value}")
