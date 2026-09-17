"""Generate Java for the GEC-ICP Function-EEDF apply and cold solve."""

from __future__ import annotations

import csv
import math
from pathlib import Path

from swarm_workflow.comsol.input.mean_energy import MeanEnergyArgument
from swarm_workflow.comsol.java import java_path, java_string

from .run_contracts import GecIcpRunMapping, GecIcpWorkflowError


APPLY_CLASS = "SwarmGecIcpApply"
CLOSURE_CLASS = "SwarmGecIcpClosure"
SOLVE_CLASS = "SwarmGecIcpSolve"
NATIVE_EEDF_AUDIT_CLASS = "SwarmGecIcpNativeEedfAudit"
TRANSPORT_TABLE = "transport_vs_mean_energy.csv"


def read_mobility_table(path: str | Path) -> tuple[list[float], list[float]]:
    """Read the one active external transport coefficient."""

    table = Path(path)
    try:
        with table.open(encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream))
    except (OSError, UnicodeError, csv.Error) as exc:
        raise GecIcpWorkflowError(f"cannot read mobility table: {table}") from exc
    energies: list[float] = []
    mobilities: list[float] = []
    try:
        for row in rows:
            energies.append(float(row["mean_energy_eV"]))
            mobilities.append(float(row["reduced_mobility_m2_V_s_m3"]))
    except (KeyError, TypeError, ValueError) as exc:
        raise GecIcpWorkflowError(
            f"invalid mobility columns or values in {table}"
        ) from exc
    if len(energies) < 2:
        raise GecIcpWorkflowError("mobility table requires at least two rows")
    if any(not math.isfinite(value) or value <= 0.0 for value in energies + mobilities):
        raise GecIcpWorkflowError("mobility table values must be positive and finite")
    if any(right <= left for left, right in zip(energies, energies[1:])):
        raise GecIcpWorkflowError("mobility mean-energy axis must increase strictly")
    return energies, mobilities


def generate_apply_java(mapping: GecIcpRunMapping) -> str:
    """Bind source EEDF and mobility while retaining COMSOL chemistry."""

    model = mapping.model_mapping.model
    table_path = mapping.bundle.path / TRANSPORT_TABLE
    energies, mobilities = read_mobility_table(table_path)
    support = MeanEnergyArgument(energies[0], energies[-1])
    component = java_string(model.component)
    physics = java_string(model.plasma_physics)
    plasma_feature = java_string(model.plasma_feature)
    target = (
        f"model.component({component}).physics({physics}).feature({plasma_feature})"
    )
    physics_target = f"model.component({component}).physics({physics})"
    lines = _java_header(
        APPLY_CLASS,
        model.input_mph,
        "swarmGecIcpApply",
    )
    lines.extend(
        [
            "    // Restricted local-mean-energy closure shared by all",
            "    // Swarm layers: external mu*N and f0(E,<E>); COMSOL",
            "    // retains diffusion, energy transport, RF and chemistry.",
            f'    {physics_target}.prop("ElectronProperties")'
            '.set("ReducedProps", true);',
            f'    {physics_target}.prop("ElectronProperties")'
            '.set("TensorElectronProps", false);',
            f'    {physics_target}.prop("ElectronProperties")'
            '.set("IncludeThermalDiffusion", false);',
            f'    {physics_target}.prop("ElectronProperties")'
            '.set("MeanElectronEnergyModel", "LocalEnergyApproximationE");',
            f'    {target}.set("SpecifyElectronDensityAndEnergy", "SpecifyMueOnly");',
            f"    try {{ model.component({component})"
            '.variable("swIcpClosureVars"); } catch (Exception missing) {',
            f'      model.component({component}).variable().create("swIcpClosureVars");',
            "    }",
            f'    model.component({component}).variable("swIcpClosureVars")'
            '.set("sw_icp_logeps", ' + java_string(support.expression("En-Ne")) + ");",
        ]
    )
    lines.extend(_mobility_function_lines(energies, mobilities))
    mobility = "exp(sw_icp_log_muN(sw_icp_logeps))*1[1/(V*m*s)]"
    zero = "0[1/(V*m*s)]"
    tensor = [mobility if index in {0, 4, 8} else zero for index in range(9)]
    lines.append(
        f'    {target}.set("muN", new String[]{{'
        + ", ".join(java_string(value) for value in tensor)
        + "});"
    )
    lines.extend(_function_eedf_lines(mapping))
    for reaction in (
        mapping.model_mapping.reactions.elastic,
        mapping.model_mapping.reactions.excitation,
        mapping.model_mapping.reactions.superelastic,
        mapping.model_mapping.reactions.ionization,
        mapping.model_mapping.reactions.stepwise_ionization,
    ):
        reaction_target = f"{physics_target}.feature({java_string(reaction)})"
        lines.extend(
            [
                f'    {reaction_target}.set("SpecifyReactionUsing", '
                '"UseCrossSectionData");',
                f'    {reaction_target}.set("eedf", "FromPhysicsInterfaceProperty");',
                f'    {reaction_target}.set("UseTownsendTwoTermsBoltzmann", false);',
            ]
        )
    lines.extend(
        [
            "    model.save(" + java_string(java_path(mapping.run.output_mph)) + ");",
            '    ModelUtil.remove("swarmGecIcpApply");',
            "  }",
            "",
            "  private static void replaceInterpolation(Model model, String tag) {",
            "    try { model.func().remove(tag); } catch (Exception ignored) {}",
            '    model.func().create(tag, "Interpolation");',
            "  }",
            "}",
            "",
        ]
    )
    return "\n".join(lines)


def generate_solve_java(mapping: GecIcpRunMapping) -> str:
    """Cold-start the single configured ICP case and export convergence data."""

    model = mapping.model_mapping.model
    parameters = mapping.model_mapping.parameters
    output = mapping.output_directory
    lines = _java_header(SOLVE_CLASS, mapping.run.output_mph, "swarmGecIcpSolve")
    time_list = " ".join(f"{value:.17g}" for value in output_times(mapping))
    lines.extend(
        [
            "    // One physical operating point; no coefficient or power sweep.",
            f"    model.param().set({java_string(parameters.power)}, "
            f"{java_string(f'{mapping.run.power_W:.17g}[W]')});",
            f"    model.param().set({java_string(parameters.gas_temperature)}, "
            f"{java_string(f'{mapping.run.gas_temperature_K:.17g}[K]')});",
            f"    model.param().set({java_string(parameters.pressure)}, "
            f"{java_string(f'{mapping.run.pressure_Pa:.17g}[Pa]')});",
            f"    model.study({java_string(model.study)})"
            f".feature({java_string(model.study_feature)})"
            f'.set("freq", {java_string(f"{mapping.run.frequency_Hz:.17g}")});',
            f"    model.study({java_string(model.study)})"
            f".feature({java_string(model.study_feature)})"
            f'.set("tlist", {java_string(time_list)});',
            "    // Discard the library saved vector: each solver layer starts",
            "    // from the same native initial conditions.",
            f"    model.sol({java_string(model.solution)}).clearSolutionData();",
            f"    model.study({java_string(model.study)}).run();",
            "    model.save(" + java_string(java_path(mapping.run.output_mph)) + ");",
            f"    int[] plasmaDomains = model.component({java_string(model.component)})"
            f".physics({java_string(model.plasma_physics)}).selection().entities(2);",
            "    if (plasmaDomains.length == 0) {",
            '      throw new IllegalStateException("plasma domain selection is empty");',
            "    }",
        ]
    )
    lines.extend(
        _numerical_export_lines(
            "swIcpVolume",
            "IntSurface",
            model.dataset,
            (
                "1",
                "plas.ne",
                "plas.ne*e_const*plas.ebar",
                "plas.n_wArs",
                "plas.n_wAr_1p",
                "mf.Qrh",
            ),
            ("m^3", "1", "eV", "1", "1", "W"),
            (
                "axisymmetric plasma volume",
                "electron inventory",
                "electron mean-energy inventory",
                "argon metastable inventory",
                "argon ion inventory",
                "plasma absorbed RF power",
            ),
            output / "convergence_volume.csv",
            selection="plasmaDomains",
            volume=True,
        )
    )
    lines.extend(
        _numerical_export_lines(
            "swIcpMinEnergy",
            "MinSurface",
            model.dataset,
            ("e_const*plas.ebar",),
            ("eV",),
            ("minimum electron mean energy",),
            output / "convergence_mean_energy_min.csv",
            selection="plasmaDomains",
        )
    )
    lines.extend(
        _numerical_export_lines(
            "swIcpMaxEnergy",
            "MaxSurface",
            model.dataset,
            ("e_const*plas.ebar",),
            ("eV",),
            ("maximum electron mean energy",),
            output / "convergence_mean_energy_max.csv",
            selection="plasmaDomains",
        )
    )
    lines.extend(
        _numerical_export_lines(
            "swIcpCoilPower",
            "EvalGlobal",
            model.dataset,
            ("mf.PCoil_1",),
            ("W",),
            ("configured coil-power readback",),
            output / "convergence_coil_power.csv",
        )
    )
    lines.extend(
        [
            '    ModelUtil.remove("swarmGecIcpSolve");',
            "  }",
            "}",
            "",
        ]
    )
    return "\n".join(lines)


def _mobility_function_lines(
    energies: list[float], mobilities: list[float]
) -> list[str]:
    table = ", ".join(
        "new String[]{"
        + java_string(f"{math.log(energy):.17e}")
        + ", "
        + java_string(f"{math.log(mobility):.17e}")
        + "}"
        for energy, mobility in zip(energies, mobilities, strict=True)
    )
    return [
        '    replaceInterpolation(model, "sw_icp_log_muN");',
        '    model.func("sw_icp_log_muN").set("source", "table");',
        '    model.func("sw_icp_log_muN").set("nargs", 1);',
        '    model.func("sw_icp_log_muN").set("argunit", "1");',
        '    model.func("sw_icp_log_muN").set("fununit", "1");',
        '    model.func("sw_icp_log_muN").set("interp", "piecewisecubic");',
        '    model.func("sw_icp_log_muN").set("extrap", "const");',
        f'    model.func("sw_icp_log_muN").set("table", new String[][]{{{table}}});',
    ]


def _function_eedf_lines(mapping: GecIcpRunMapping) -> list[str]:
    model = mapping.model_mapping.model
    spec = mapping.closure.function_eedf
    component = java_string(model.component)
    physics = (
        f"model.component({component}).physics({java_string(model.plasma_physics)})"
    )
    tag = java_string(spec.function_tag)
    source = java_string(java_path(mapping.bundle.path / spec.table))
    names = "new String[][]{new String[]{" + tag + ', "1"}}'
    common = [
        f"    try {{ model.component({component}).func().remove({tag}); }} "
        "catch (Exception missing) {}",
        f'    model.component({component}).func().create({tag}, "Interpolation");',
        f'    model.component({component}).func({tag}).set("source", "file");',
    ]
    return [
        *common,
        f'    model.component({component}).func({tag}).set("filename", {source});',
        f'    model.component({component}).func({tag}).set("nargs", 2);',
        f'    model.component({component}).func({tag}).set("struct", "spreadsheet");',
        f'    model.component({component}).func({tag}).set("scaledata", "auto");',
        f'    model.component({component}).func({tag}).set("funcnametable", {names});',
        f'    model.component({component}).func({tag}).set("argunit", "eV,eV");',
        f'    model.component({component}).func({tag}).set("fununit", "1");',
        f'    model.component({component}).func({tag}).set("interp", "linear");',
        f'    model.component({component}).func({tag}).set("extrap", "const");',
        f"    model.component({component}).func({tag}).importData();",
        f'    {physics}.prop("EEDFSettings").set("eedf", {tag});',
    ]


def output_times(mapping: GecIcpRunMapping) -> tuple[float, ...]:
    final_time = mapping.run.final_time_s
    first_time = min(1.0e-8, final_time)
    if final_time <= first_time:
        return (0.0, final_time)
    decades = math.log10(final_time) - math.log10(first_time)
    intervals = max(1, math.ceil(decades * mapping.run.output_points_per_decade))
    times = [0.0]
    for index in range(intervals + 1):
        fraction = index / intervals
        times.append(first_time * (final_time / first_time) ** fraction)
    times[-1] = final_time
    return tuple(times)


def _numerical_export_lines(
    tag: str,
    kind: str,
    dataset: str,
    expressions: tuple[str, ...],
    units: tuple[str, ...],
    descriptions: tuple[str, ...],
    output: Path,
    *,
    selection: str | None = None,
    volume: bool = False,
) -> list[str]:
    qtag = java_string(tag)
    table = java_string(f"{tag}Table")
    export = java_string(f"{tag}Export")
    expr = ", ".join(java_string(value) for value in expressions)
    unit = ", ".join(java_string(value) for value in units)
    descr = ", ".join(java_string(value) for value in descriptions)
    if not (len(expressions) == len(units) == len(descriptions)):
        raise GecIcpWorkflowError(
            "COMSOL numerical export expressions, units, and descriptions must align"
        )
    lines = [
        f'    model.result().table().create({table}, "Table");',
        f"    model.result().numerical().create({qtag}, {java_string(kind)});",
        f'    model.result().numerical({qtag}).set("data", {java_string(dataset)});',
        f'    model.result().numerical({qtag}).set("expr", new String[]{{{expr}}});',
        f'    model.result().numerical({qtag}).set("unit", new String[]{{{unit}}});',
        f'    model.result().numerical({qtag}).set("descr", new String[]{{{descr}}});',
    ]
    if selection is not None:
        lines.append(
            f"    model.result().numerical({qtag}).selection().set({selection});"
        )
    if volume:
        lines.append(f'    model.result().numerical({qtag}).set("intvolume", true);')
    lines.extend(
        [
            f'    model.result().numerical({qtag}).set("table", {table});',
            f"    model.result().numerical({qtag}).setResult();",
            f'    model.result().export().create({export}, "Table");',
            f'    model.result().export({export}).set("table", {table});',
            f'    model.result().export({export}).set("header", "on");',
            f'    model.result().export({export}).set("prec", "full");',
            f'    model.result().export({export}).set("ifexists", "overwrite");',
            f'    model.result().export({export}).set("filename", '
            + java_string(java_path(output))
            + ");",
            f"    model.result().export({export}).run();",
        ]
    )
    return lines


def _java_header(class_name: str, input_mph: Path, model_name: str) -> list[str]:
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


__all__ = [
    "APPLY_CLASS",
    "CLOSURE_CLASS",
    "NATIVE_EEDF_AUDIT_CLASS",
    "SOLVE_CLASS",
    "TRANSPORT_TABLE",
    "generate_apply_java",
    "generate_solve_java",
    "output_times",
    "read_mobility_table",
]
