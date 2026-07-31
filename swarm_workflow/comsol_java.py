"""Generate Java source for applying mapped CSV tables to COMSOL models."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable

from .comsol_mapping import (
    ClosureMapping,
    ComsolModelMapping,
    FunctionMapping,
    MeanEnergyLookup,
    ReactionLookupMapping,
)


DEFAULT_CLASS_NAME = "SwarmComsolApply"
ProbeSpec = tuple[str, str, str]
PROBE_CSV_BEGIN = "SWARM_COMSOL_PROBE_CSV_BEGIN"
PROBE_CSV_END = "SWARM_COMSOL_PROBE_CSV_END"
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
        value is not None
        for value in (pressure_Pa, gas_temperature_K, mesh_elements)
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
    if all(
        value is None
        for value in (pressure_Pa, gas_temperature_K, mesh_elements)
    ):
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
                    f".feature({plasma_feature}).set(\"pA_src\", \"userdef\");"
                ),
                (
                    f"    model.component({component}).physics({physics})"
                    f".feature({plasma_feature}).set(\"pA\", \"p0\");"
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
                    f".feature({plasma_feature}).set(\"T_src\", \"userdef\");"
                ),
                (
                    f"    model.component({component}).physics({physics})"
                    f".feature({plasma_feature}).set(\"T\", "
                    f"\"{temperature}[K]\");"
                ),
                (
                    "    String configuredGasTemperature = "
                    f"model.component({component}).physics({physics})"
                    f".feature({plasma_feature}).getString(\"T\");"
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
                    f"    model.component({component}).mesh(\"mesh1\")"
                    ".feature(\"edg1\").feature(\"dis1\")"
                    f'.set("elemcount", {mesh_elements});'
                ),
                f'    model.component({component}).mesh("mesh1").run();',
                (
                    "    String configuredMeshElements = "
                    f"model.component({component}).mesh(\"mesh1\")"
                    ".feature(\"edg1\").feature(\"dis1\")"
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
            "    model.save("
            + java_string(java_path(mapping.model.output_mph))
            + ");",
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
        f"    model.func({tag}).set(\"source\", \"file\");",
        f"    model.func({tag}).set(\"filename\", "
        + java_string(java_path(function.path))
        + ");",
        f"    model.func({tag}).set(\"nargs\", {function.nargs});",
        f"    model.func({tag}).set(\"argunit\", {java_string(function.argunit)});",
        f"    model.func({tag}).set(\"fununit\", {java_string(function.fununit)});",
        f"    model.func({tag}).set(\"interp\", {java_string(function.interp)});",
        f"    model.func({tag}).set(\"extrap\", {java_string(function.extrap)});",
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
    x_values, y_values = _csv_values(
        closure.mean_energy.path,
        closure.mean_energy.x_column,
        closure.mean_energy.y_column,
        "mean-energy lookup",
    )
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
        (
            f"    if (!{expected_formulation}.equals(configuredMeanEnergyModel)) {{"
        ),
        (
            "      throw new IllegalStateException("
            '"COMSOL mean-energy formulation read-back mismatch: expected " + '
            f"{expected_formulation} + \", got \" + configuredMeanEnergyModel);"
        ),
        "    }",
        (
            f'    System.out.println("{MEAN_ENERGY_MODEL_LOG_PREFIX}"'
            " + configuredMeanEnergyModel);"
        ),
        f"    {target}.set(\"SpecifyElectronDensityAndEnergy\", \"UseLookupTables\");",
        f"    {target}.set(\"SpecifyMeanElectronEnergy\", \"MeanEnergyTable\");",
        f"    {target}.set(\"enrgXdata\", new double[]{{{_java_double_array(x_values)}}});",
        f"    {target}.set(\"enrgYdata\", new double[]{{{_java_double_array(y_values)}}});",
    ]
    for x_property, y_property, column in transport_tables:
        x_transport, y_transport = _csv_values(
            transport,
            "mean_energy_eV",
            column,
            column,
        )
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


def _csv_lookup_values(
    table: MeanEnergyLookup,
) -> tuple[list[float], list[float]]:
    return _csv_values(
        table.path,
        table.x_column,
        table.y_column,
        "mean-energy lookup",
    )


def _csv_values(
    path: Path,
    x_column: str,
    y_column: str,
    name: str,
) -> tuple[list[float], list[float]]:
    rows: dict[float, float] = {}
    with path.open("r", encoding="utf-8", newline="") as fp:
        for row in csv.DictReader(fp):
            try:
                x_value = float(row[x_column])
                y_value = float(row[y_column])
            except (KeyError, ValueError) as exc:
                raise ValueError(f"invalid numeric lookup row in {path}") from exc
            if not math.isfinite(x_value) or not math.isfinite(y_value):
                raise ValueError(f"{name} has nonfinite values")
            rows[x_value] = y_value
    if len(rows) < 2:
        raise ValueError(f"{name} needs at least two points")
    x_values = sorted(rows)
    return x_values, [rows[value] for value in x_values]


def _reaction_lookup_lines(
    mapping: ComsolModelMapping,
    lookup: ReactionLookupMapping,
) -> list[str]:
    if lookup.form != "townsend":
        raise ValueError(f"unsupported reaction lookup form: {lookup.form}")
    x_values, y_values = _reaction_lookup_values(lookup)
    component = java_string(mapping.model.component)
    physics = java_string(lookup.physics)
    feature = java_string(lookup.feature)
    return [
        f"    // Mapping-driven Townsend lookup for {lookup.name} "
        f"({lookup.source_model}).",
        f"    model.component({component}).physics({physics}).feature({feature})"
        ".set(\"SpecifyReactionUsing\", \"UseLookupTable\");",
        f"    model.component({component}).physics({physics}).feature({feature})"
        ".set(\"RateConstantForm\", \"UseTownsend\");",
        f"    model.component({component}).physics({physics}).feature({feature})"
        ".set(\"UseTownsendTwoTermsBoltzmann\", false);",
        f"    model.component({component}).physics({physics}).feature({feature})"
        f".set(\"xtownratedata\", new double[]{{{_java_double_array(x_values)}}});",
        f"    model.component({component}).physics({physics}).feature({feature})"
        f".set(\"ytownratedata\", new double[]{{{_java_double_array(y_values)}}});",
        "",
    ]


def _reaction_lookup_values(
    lookup: ReactionLookupMapping,
) -> tuple[list[float], list[float]]:
    rows: dict[float, float] = {}
    with lookup.path.open("r", encoding="utf-8", newline="") as fp:
        for row in csv.DictReader(fp):
            if row.get("process_type") != lookup.process_type:
                continue
            try:
                x_value = float(row[lookup.x_column])
                y_value = float(row[lookup.y_column])
            except (KeyError, ValueError) as exc:
                raise ValueError(
                    f"invalid numeric reaction lookup row in {lookup.path}"
                ) from exc
            if not math.isfinite(x_value) or not math.isfinite(y_value):
                raise ValueError(
                    f"reaction lookup {lookup.name} has nonfinite values"
                )
            if y_value < 0.0:
                raise ValueError(
                    f"reaction lookup {lookup.name} has negative values"
                )
            rows[x_value] = y_value
    if not rows:
        raise ValueError(
            f"reaction lookup {lookup.name} has no rows for process_type "
            f"{lookup.process_type!r}"
        )
    x_values = sorted(rows)
    if len(x_values) < 2:
        raise ValueError(f"reaction lookup {lookup.name} needs at least two points")
    return x_values, [rows[value] for value in x_values]


def _java_double_array(values: Iterable[float]) -> str:
    return ", ".join(f"{float(value):.17e}" for value in values)


def generate_probe_table_export_java_source(
    *,
    class_name: str,
    model_name: str,
    input_mph: Path,
    output_csv: Path,
    eval_tag: str,
    probes: Iterable[ProbeSpec],
    pre_export_lines: Iterable[str] = (),
) -> str:
    """Return Java source that evaluates mapped result probes to a CSV file."""

    _validate_java_class_name(class_name)
    probe_list = list(probes)
    expressions = ", ".join(java_string(expression) for _, expression, _ in probe_list)
    units = ", ".join(java_string(unit) for _, _, unit in probe_list)
    names = ",".join(name for name, _, _ in probe_list)
    lines = [
        "import com.comsol.model.*;",
        "import com.comsol.model.util.*;",
        "",
        f"public class {class_name} {{",
        "  // output_csv: " + java_path(output_csv),
        "  public static void main(String[] args) throws Exception {",
        "    Model model = ModelUtil.loadCopy("
        + java_string(model_name)
        + ", "
        + java_string(java_path(input_mph))
        + ");",
        *pre_export_lines,
        "    String evalTag = " + java_string(eval_tag) + ";",
        "    try {",
        "      model.result().numerical().remove(evalTag);",
        "    } catch (Exception ex) {",
        "    }",
        '    model.result().numerical().create(evalTag, "Eval");',
        "    model.result().numerical(evalTag).set(\"expr\", new String[]{"
        + expressions
        + "});",
        "    model.result().numerical(evalTag).set(\"unit\", new String[]{"
        + units
        + "});",
        "    double[][][] values = model.result().numerical(evalTag).getData();",
        "    double[][] coords = model.result().numerical(evalTag).getCoordinates();",
        "    int n = model.result().numerical(evalTag).getNData();",
        "    System.out.println(" + java_string(PROBE_CSV_BEGIN) + ");",
        "    System.out.println(" + java_string(names) + ");",
        "    for (int i = 0; i < n; i++) {",
        "      StringBuilder row = new StringBuilder();",
        f"      for (int j = 0; j < {len(probe_list)}; j++) {{",
        "        if (j > 0) row.append(',');",
        "        double[][] exprValues = values[j];",
        "        double value = exprValues[exprValues.length - 1][i];",
        "        row.append(Double.toString(value));",
        "      }",
        "      System.out.println(row.toString());",
        "    }",
        "    System.out.println(" + java_string(PROBE_CSV_END) + ");",
        "  }",
        "}",
        "",
    ]
    return "\n".join(lines)


def extract_probe_csv_from_stdout(stdout_path: Path, output_csv: Path) -> None:
    """Write probe CSV rows emitted between stdout markers, if present."""

    if not stdout_path.exists():
        return
    lines = stdout_path.read_text(encoding="utf-8", errors="replace").splitlines()
    try:
        start = lines.index(PROBE_CSV_BEGIN) + 1
        end = lines.index(PROBE_CSV_END, start)
    except ValueError:
        return
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    output_csv.write_text("\n".join(lines[start:end]) + "\n", encoding="utf-8")


def java_path(value: object) -> str:
    """Return a path string that is stable inside generated Java source."""

    return str(value).replace("\\", "/")


def java_string(value: object) -> str:
    """Return ``value`` as an escaped Java string literal."""

    text = str(value)
    escaped = (
        text.replace("\\", "\\\\")
        .replace('"', '\\"')
        .replace("\n", "\\n")
        .replace("\r", "\\r")
    )
    return f'"{escaped}"'


def _validate_java_class_name(value: str) -> None:
    if not value or not value.replace("_", "a").isalnum() or value[0].isdigit():
        raise ValueError(f"invalid Java class name: {value}")
