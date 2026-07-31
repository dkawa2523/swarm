"""GEC CCP-specific COMSOL preparation, execution, and result export.

The local model is a time-periodic plasma model, so the electron mean energy
must remain a solved dependent variable.  Swarm data replace the four reduced
transport tables and the three reaction-rate tables as functions of mean
electron energy.  The DC swarm sweep is therefore a local mean-energy closure,
not an RF Boltzmann solve.
"""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
import re
from typing import Any, Iterable
from zipfile import ZipFile

import yaml

from ._io import write_json
from .comsol_adapter import (
    ComsolAdapterError,
    ComsolExecutionSummary,
    execute_generated_comsol_java,
)
from .comsol_java import java_path, java_string


APPLY_CLASS = "SwarmGecCcpApply"
BASELINE_RUN_CLASS = "SwarmGecCcpBaselineRun"
EXTERNAL_RUN_CLASS = "SwarmGecCcpExternalRun"
BASELINE_EXPORT_CLASS = "SwarmGecCcpBaselineExport"
EXTERNAL_EXPORT_CLASS = "SwarmGecCcpExternalExport"
REQUIRED_TABLE_COLUMNS = {
    "transport_vs_mean_energy.csv": (
        "mean_energy_eV",
        "reduced_mobility_m2_V_s_m3",
        "reduced_diffusion_L_m2_s_m3",
        "reduced_electron_energy_mobility_m2_V_s_m3",
        "reduced_electron_energy_diffusion_m2_s_m3",
    ),
    "rates_vs_mean_energy.csv": (
        "mean_energy_eV",
        "process_type",
        "rate_coefficient_m3_s",
    ),
    "eedf_f0.csv": (
        "electron_energy_eV",
        "mean_energy_eV",
        "E_over_N_Td",
        "eepf_eV_m32",
    ),
    "quality.csv": (
        "E_over_N_Td",
        "passed",
        "eedf_normalization_error",
    ),
}
AVOGADRO_PER_MOL = "6.02214076e23[1/mol]"


class GecCcpWorkflowError(RuntimeError):
    """Raised when the GEC CCP workflow cannot be prepared or completed."""


@dataclass(frozen=True, slots=True)
class GecModelSpec:
    input_mph: Path
    baseline_output_mph: Path
    output_mph: Path
    component: str
    physics: str
    plasma_feature: str
    study: str
    conversion_study: str
    expected_original_eedf: str


@dataclass(frozen=True, slots=True)
class GecReactionSpec:
    name: str
    feature: str
    process_type: str


@dataclass(frozen=True, slots=True)
class GecRunSpec:
    power_parameter: str
    power_W: float
    initial_mean_energy_scale: float
    axis_dataset: str
    radial_dataset: str
    period_dataset: str
    phase_dataset: str
    waveform_dataset: str


@dataclass(frozen=True, slots=True)
class GecPathSpec:
    path: Path


@dataclass(frozen=True, slots=True)
class GecResultSpec:
    output_directory: Path


@dataclass(frozen=True, slots=True)
class GecCcpMapping:
    path: Path
    root: Path
    model: GecModelSpec
    bundle: GecPathSpec
    reactions: tuple[GecReactionSpec, ...]
    run: GecRunSpec
    results: GecResultSpec
    logs: GecPathSpec


@dataclass(frozen=True, slots=True)
class MphContract:
    comsol_version: str | None
    physics_operation: str
    physics_tag: str
    original_eedf: str
    plasma_feature: str
    reaction_features: tuple[str, ...]
    studies: tuple[str, ...]
    datasets: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class GecCcpPlan:
    mapping: GecCcpMapping
    contract: MphContract
    output_directory: Path
    plan_json: Path
    baseline_run_java: Path
    baseline_export_java: Path
    apply_java: Path
    external_run_java: Path
    external_export_java: Path
    expected_result_files: tuple[Path, ...]


@dataclass(frozen=True, slots=True)
class GecCcpRunSummary:
    plan: GecCcpPlan
    executions: tuple[ComsolExecutionSummary, ...]
    status_json: Path


def load_gec_ccp_mapping(
    mapping_path: str | Path,
    *,
    bundle_path: str | Path | None = None,
) -> GecCcpMapping:
    path = Path(mapping_path).resolve()
    raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    if not isinstance(raw, dict):
        raise GecCcpWorkflowError("GEC CCP mapping root must be a mapping")
    if int(raw.get("schema_version", 0)) != 2:
        raise GecCcpWorkflowError(
            "GEC CCP mapping requires schema_version: 2; migrate run.powers_W "
            "to the single direct-solve setting run.power_W"
        )
    root = path.parent
    model_raw = _mapping(raw, "model")
    bundle_raw = _mapping(raw, "bundle")
    run_raw = _mapping(raw, "run")
    results_raw = _mapping(raw, "results")
    logs_raw = _mapping(raw, "logs")
    model = GecModelSpec(
        input_mph=_resolve(root, model_raw.get("input_mph"), "model.input_mph"),
        baseline_output_mph=_resolve(
            root,
            model_raw.get("baseline_output_mph"),
            "model.baseline_output_mph",
        ),
        output_mph=_resolve(
            root,
            model_raw.get("external_output_mph"),
            "model.external_output_mph",
        ),
        component=_text(model_raw, "component"),
        physics=_text(model_raw, "physics"),
        plasma_feature=_text(model_raw, "plasma_feature"),
        study=_text(model_raw, "time_periodic_study"),
        conversion_study=_text(model_raw, "conversion_study"),
        expected_original_eedf=_text(model_raw, "expected_original_eedf"),
    )
    raw_reactions = raw.get("reactions")
    if not isinstance(raw_reactions, list) or not raw_reactions:
        raise GecCcpWorkflowError("reactions must be a non-empty list")
    reactions = tuple(
        GecReactionSpec(
            name=_text(item, "name"),
            feature=_text(item, "feature"),
            process_type=_text(item, "process_type"),
        )
        for item in raw_reactions
        if isinstance(item, dict)
    )
    if len(reactions) != len(raw_reactions):
        raise GecCcpWorkflowError("every reactions entry must be a mapping")
    if "powers_W" in run_raw:
        raise GecCcpWorkflowError(
            "run.powers_W is obsolete; use one direct-solve value in run.power_W"
        )
    try:
        power = float(run_raw["power_W"])
    except (KeyError, TypeError, ValueError) as exc:
        raise GecCcpWorkflowError("run.power_W must be a positive number") from exc
    if not math.isfinite(power) or power <= 0.0:
        raise GecCcpWorkflowError("run.power_W must be a positive finite value")
    try:
        energy_scale = float(run_raw["initial_mean_energy_scale"])
    except (KeyError, TypeError, ValueError) as exc:
        raise GecCcpWorkflowError(
            "run.initial_mean_energy_scale must be a positive number"
        ) from exc
    if not math.isfinite(energy_scale) or energy_scale <= 0.0:
        raise GecCcpWorkflowError(
            "run.initial_mean_energy_scale must be a positive finite value"
        )
    run = GecRunSpec(
        power_parameter=_text(run_raw, "power_parameter"),
        power_W=power,
        initial_mean_energy_scale=energy_scale,
        axis_dataset=_text(run_raw, "axis_dataset"),
        radial_dataset=_text(run_raw, "radial_dataset"),
        period_dataset=_text(run_raw, "period_dataset"),
        phase_dataset=_text(run_raw, "phase_dataset"),
        waveform_dataset=_text(run_raw, "waveform_dataset"),
    )
    return GecCcpMapping(
        path=path,
        root=root,
        model=model,
        bundle=GecPathSpec(
            Path(bundle_path).resolve()
            if bundle_path is not None
            else _resolve(root, bundle_raw.get("path"), "bundle.path")
        ),
        reactions=reactions,
        run=run,
        results=GecResultSpec(
            _resolve(
                root,
                results_raw.get("output_directory"),
                "results.output_directory",
            )
        ),
        logs=GecPathSpec(_resolve(root, logs_raw.get("path"), "logs.path")),
    )


def inspect_gec_ccp_mph(path: str | Path) -> MphContract:
    mph = Path(path)
    if not mph.exists():
        raise GecCcpWorkflowError(f"COMSOL model does not exist: {mph}")
    try:
        with ZipFile(mph) as archive:
            xml = archive.read("dmodel.xml").decode("utf-8", errors="replace")
            file_version = (
                archive.read("fileversion").decode("utf-8", errors="replace")
                if "fileversion" in archive.namelist()
                else ""
            )
    except (OSError, KeyError) as exc:
        raise GecCcpWorkflowError(f"cannot read dmodel.xml from {mph}") from exc
    physics = re.search(
        r'<Physics\b[^>]*op="([^"]+)"[^>]*tag="([^"]+)"', xml
    )
    if physics is None:
        raise GecCcpWorkflowError("model has no readable physics declaration")
    eedf_values = _property_values(xml, "eedf")
    eedf = next(
        (
            value
            for value in reversed(eedf_values)
            if value != "FromPhysicsInterfaceProperty"
        ),
        eedf_values[-1],
    )
    plasma = _tag_for_operation(xml, "PhysicsFeature", "PlasmaEsModel")
    reactions = tuple(
        re.findall(
            r'<PhysicsFeature\b[^>]*op="ElectronImpactReaction"[^>]*tag="([^"]+)"',
            xml,
        )
    )
    studies = tuple(re.findall(r'<Study\b[^>]*tag="([^"]+)"', xml))
    datasets = tuple(re.findall(r'<DatasetFeature\b[^>]*tag="([^"]+)"', xml))
    version = re.search(r"COMSOL\s+([0-9.]+)", file_version)
    return MphContract(
        comsol_version=version.group(1) if version else None,
        physics_operation=physics.group(1),
        physics_tag=physics.group(2),
        original_eedf=eedf,
        plasma_feature=plasma,
        reaction_features=reactions,
        studies=studies,
        datasets=datasets,
    )


def prepare_gec_ccp_run(
    mapping_path: str | Path,
    *,
    bundle_path: str | Path | None = None,
    write_java: bool = True,
) -> GecCcpPlan:
    mapping = load_gec_ccp_mapping(mapping_path, bundle_path=bundle_path)
    contract = inspect_gec_ccp_mph(mapping.model.input_mph)
    _validate_contract(mapping, contract)
    bundle_summary = validate_gec_ccp_bundle(mapping)
    output = mapping.results.output_directory
    output.mkdir(parents=True, exist_ok=True)
    mapping.model.baseline_output_mph.parent.mkdir(parents=True, exist_ok=True)
    mapping.model.output_mph.parent.mkdir(parents=True, exist_ok=True)
    baseline_dir = output / "builtin_druyvesteyn"
    external_dir = output / "swarm_tables"
    baseline_dir.mkdir(parents=True, exist_ok=True)
    external_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "baseline_run": output / f"{BASELINE_RUN_CLASS}.java",
        "baseline_export": output / f"{BASELINE_EXPORT_CLASS}.java",
        "apply": output / f"{APPLY_CLASS}.java",
        "external_run": output / f"{EXTERNAL_RUN_CLASS}.java",
        "external_export": output / f"{EXTERNAL_EXPORT_CLASS}.java",
    }
    expected = tuple(
        directory / name
        for directory in (baseline_dir, external_dir)
        for name in (
            "axis_period_average.csv",
            "radial_period_average.csv",
            "axis_phase_resolved.csv",
            "electrode_waveform.csv",
            "closure_phase.csv",
        )
    )
    plan_json = output / "gec_ccp_plan.json"
    plan = GecCcpPlan(
        mapping=mapping,
        contract=contract,
        output_directory=output,
        plan_json=plan_json,
        baseline_run_java=paths["baseline_run"],
        baseline_export_java=paths["baseline_export"],
        apply_java=paths["apply"],
        external_run_java=paths["external_run"],
        external_export_java=paths["external_export"],
        expected_result_files=expected,
    )
    if write_java:
        plan.baseline_run_java.write_text(
            generate_run_java(
                mapping,
                class_name=BASELINE_RUN_CLASS,
                input_mph=mapping.model.input_mph,
                output_mph=mapping.model.baseline_output_mph,
            ),
            encoding="utf-8",
        )
        plan.baseline_export_java.write_text(
            generate_export_java(
                mapping,
                class_name=BASELINE_EXPORT_CLASS,
                input_mph=mapping.model.baseline_output_mph,
                output_dir=baseline_dir,
            ),
            encoding="utf-8",
        )
        plan.apply_java.write_text(generate_apply_java(mapping), encoding="utf-8")
        plan.external_run_java.write_text(
            generate_run_java(
                mapping,
                class_name=EXTERNAL_RUN_CLASS,
                input_mph=mapping.model.output_mph,
                output_mph=mapping.model.output_mph,
                robust_nonlinear=True,
            ),
            encoding="utf-8",
        )
        plan.external_export_java.write_text(
            generate_export_java(
                mapping,
                class_name=EXTERNAL_EXPORT_CLASS,
                input_mph=mapping.model.output_mph,
                output_dir=external_dir,
            ),
            encoding="utf-8",
        )
    write_json(
        plan_json,
        {
            "stage": "prepare-gec-ccp",
            "status": "ready",
            "mapping": str(mapping.path),
            "bundle": bundle_summary,
            "model_contract": asdict(contract),
            "closure": {
                "kind": "direct_atomic_local_mean_energy_closure",
                "swarm_field": "steady_dc",
                "comsol_field": "time_periodic_rf",
                "mean_energy": "solved_by_COMSOL",
                "power_W": mapping.run.power_W,
                "coefficient_sweep": False,
                "initial_guess": {
                    "source": "converged built-in Druyvesteyn periodic solution",
                    "electron_mean_energy_scale": (
                        mapping.run.initial_mean_energy_scale
                    ),
                    "other_dependent_variables": "unchanged",
                    "basis": (
                        "match the characteristic Swarm ionization rate to "
                        "the baseline ionization rate"
                    ),
                },
                "transport_tables": [
                    "muN",
                    "DN",
                    "mueN",
                    "DeN",
                ],
                "reaction_tables": [
                    f"{reaction.feature}:{reaction.process_type}"
                    for reaction in mapping.reactions
                ],
                "reaction_rate_units": {
                    "swarm": "m^3/s per target particle",
                    "comsol": "m^3/(mol*s)",
                    "conversion": f"{AVOGADRO_PER_MOL} * k_swarm",
                },
                "interpolation": {
                    "transport": "shape-preserving piecewise cubic",
                    "reaction_rates": (
                        "shape-preserving piecewise cubic in natural-log(k)"
                    ),
                    "reason": (
                        "positive rates span many decades; log interpolation "
                        "preserves the tabulated points and smooths relative slopes"
                    ),
                },
                "original_eedf": contract.original_eedf,
                "external_eedf_role": "audit_only; rates and transport are imported",
                "newton_iteration_guard": {
                    "mean_energy_eV": [-1.0e6, 1.0e6],
                    "transport": "constant endpoint hold",
                    "inelastic_rates_below_zero": 0.0,
                    "role": "numerical guard outside the validated physical range",
                },
                "direct_solve": {
                    "initial_solution": "converged built-in Druyvesteyn model",
                    "activation": "four transport coefficients and three rates together",
                    "solver_settings": {
                        "baseline": "preserve original model settings",
                        "external": (
                            "COMSOL automatic highly nonlinear Newton; "
                            "no coefficient sweep"
                        ),
                        "maximum_iterations": 200,
                    },
                    "stabilization": {
                        "source": True,
                        "reaction_source": True,
                        "zeta": 1,
                        "iota": 1,
                    },
                },
            },
            "generated_java": {key: str(value) for key, value in paths.items()},
            "expected_results": [str(path) for path in expected],
        },
    )
    return plan


def validate_gec_ccp_bundle(mapping: GecCcpMapping) -> dict[str, Any]:
    manifest_path = mapping.bundle.path / "manifest.json"
    if not manifest_path.exists():
        raise GecCcpWorkflowError(f"missing COMSOL bundle manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "ok":
        raise GecCcpWorkflowError("COMSOL bundle manifest status is not ok")
    tables = manifest.get("tables", {})
    for name, columns in REQUIRED_TABLE_COLUMNS.items():
        entry = tables.get(name)
        listed = entry.get("columns", []) if isinstance(entry, dict) else []
        missing = [column for column in columns if column not in listed]
        if missing:
            raise GecCcpWorkflowError(
                f"bundle table {name} missing columns: {', '.join(missing)}"
            )
        if not (mapping.bundle.path / name).exists():
            raise GecCcpWorkflowError(f"bundle table does not exist: {name}")
    rate_rows = _read_csv(mapping.bundle.path / "rates_vs_mean_energy.csv")
    present = {row.get("process_type") for row in rate_rows}
    missing_processes = [
        reaction.process_type
        for reaction in mapping.reactions
        if reaction.process_type not in present
    ]
    if missing_processes:
        raise GecCcpWorkflowError(
            "bundle lacks mapped reaction process types: "
            + ", ".join(missing_processes)
        )
    monotonic = manifest.get("monotonicity", {}).get(
        "mean_energy_strictly_monotonic"
    )
    if monotonic is not True:
        raise GecCcpWorkflowError(
            "GEC CCP lookup requires strictly monotonic mean energy"
        )
    transport_rows = _read_csv(
        mapping.bundle.path / "transport_vs_mean_energy.csv"
    )
    for column in REQUIRED_TABLE_COLUMNS["transport_vs_mean_energy.csv"]:
        if column == "mean_energy_eV":
            continue
        _positive_numeric_values(transport_rows, column)
    quality_rows = _read_csv(mapping.bundle.path / "quality.csv")
    failed_quality = [
        row.get("E_over_N_Td", "?")
        for row in quality_rows
        if str(row.get("passed", "")).strip().lower() not in {"1", "true"}
    ]
    if failed_quality:
        raise GecCcpWorkflowError(
            "bundle has failed quality points: " + ", ".join(failed_quality)
        )
    normalization_errors = [
        abs(float(row["eedf_normalization_error"]))
        for row in quality_rows
        if row.get("eedf_normalization_error") not in {None, ""}
    ]
    process_stats: dict[str, dict[str, Any]] = {}
    for reaction in mapping.reactions:
        selected = [
            row
            for row in rate_rows
            if row.get("process_type") == reaction.process_type
        ]
        values = _nonnegative_numeric_values(
            selected, "rate_coefficient_m3_s"
        )
        process_stats[reaction.process_type] = {
            "points": len(values),
            "minimum_m3_s": min(values),
            "maximum_m3_s": max(values),
        }
    return {
        "path": str(mapping.bundle.path),
        "source": manifest.get("source"),
        "valid_ranges": manifest.get("valid_ranges", {}),
        "mean_energy_strictly_monotonic": True,
        "transport_points": len(transport_rows),
        "quality_points": len(quality_rows),
        "maximum_eedf_normalization_error": (
            max(normalization_errors) if normalization_errors else None
        ),
        "reaction_rate_ranges": process_stats,
        "tables_checked": sorted(REQUIRED_TABLE_COLUMNS),
    }


def generate_apply_java(
    mapping: GecCcpMapping,
    *,
    input_mph: Path | None = None,
) -> str:
    transport = mapping.bundle.path / "transport_vs_mean_energy.csv"
    rows = _read_csv(transport)
    target = (
        f"model.component({java_string(mapping.model.component)})"
        f".physics({java_string(mapping.model.physics)})"
        f".feature({java_string(mapping.model.plasma_feature)})"
    )
    lines = _java_header(
        APPLY_CLASS,
        input_mph or mapping.model.baseline_output_mph,
        "swarmGecApply",
    )
    lines.extend(
        [
            "    // Preserve COMSOL's solved mean energy and regularize en/ne only",
            "    // during nonphysical Newton iterates at vanishing density.",
            f'    {target}.set("SpecifyElectronDensityAndEnergy", "SpecifyAll");',
            "    // COMSOL recommends both source stabilizations for log-form",
            "    // plasma equations near vanishing ne, en, or mass fractions.",
            f'    model.component({java_string(mapping.model.component)})'
            f'.physics({java_string(mapping.model.physics)})'
            '.prop("Stabilization").set("SourceStabilization", true);',
            f'    model.component({java_string(mapping.model.component)})'
            f'.physics({java_string(mapping.model.physics)})'
            '.prop("Stabilization").set("zeta", "1");',
            f'    model.component({java_string(mapping.model.component)})'
            f'.physics({java_string(mapping.model.physics)})'
            '.prop("Stabilization").set('
            '"ReactionSourceStabilization", true);',
            f'    model.component({java_string(mapping.model.component)})'
            f'.physics({java_string(mapping.model.physics)})'
            '.prop("Stabilization").set("iota", "1");',
        ]
    )
    safe_mean_energy = "ptp.en/max(ptp.ne,1[1/m^3])"
    for tag, property_name, column, unit in (
        (
            "sw_muN_e",
            "mue",
            "reduced_mobility_m2_V_s_m3",
            "1/(V*m*s)",
        ),
        (
            "sw_DeN_e",
            "De",
            "reduced_diffusion_L_m2_s_m3",
            "1/(m*s)",
        ),
        (
            "sw_muenN_e",
            "muen",
            "reduced_electron_energy_mobility_m2_V_s_m3",
            "1/(V*m*s)",
        ),
        (
            "sw_DenN_e",
            "Den",
            "reduced_electron_energy_diffusion_m2_s_m3",
            "1/(m*s)",
        ),
    ):
        x_values, y_values = _lookup(rows, "mean_energy_eV", column)
        x_values, y_values = _guard_mean_energy_lookup(
            x_values,
            y_values,
            low_value=y_values[0],
        )
        lines.extend(
            _inline_interpolation_lines(
                tag,
                x_values,
                y_values,
                argunit="V",
                fununit=unit,
            )
        )
        # Convert the reduced coefficient using COMSOL's solved local neutral
        # density.  Keep this dependence in the Jacobian: it is part of the
        # requested closure, not a continuation or frozen-coefficient device.
        expression = f"{tag}({safe_mean_energy})/max(ptp.Nn,1[1/m^3])"
        tensor = ", ".join(
            java_string(expression if index in {0, 4, 8} else "0")
            for index in range(9)
        )
        lines.append(
            f"    {target}.set({java_string(property_name)}, "
            f"new String[]{{{tensor}}});"
    )
    rate_rows = _read_csv(mapping.bundle.path / "rates_vs_mean_energy.csv")
    for reaction in mapping.reactions:
        selected = [
            row
            for row in rate_rows
            if row.get("process_type") == reaction.process_type
        ]
        x_values, y_values = _lookup(
            selected, "mean_energy_eV", "rate_coefficient_m3_s"
        )
        x_values, y_values = _guard_mean_energy_lookup(
            x_values,
            y_values,
            low_value=y_values[0],
        )
        reaction_target = (
            f"model.component({java_string(mapping.model.component)})"
            f".physics({java_string(mapping.model.physics)})"
            f".feature({java_string(reaction.feature)})"
        )
        function_tag = f"sw_logk_{reaction.process_type}_e"
        log_rate_values = [math.log(max(value, 1.0e-300)) for value in y_values]
        lines.extend(
            _inline_interpolation_lines(
                function_tag,
                x_values,
                log_rate_values,
                argunit="V",
                fununit="1",
            )
        )
        particle_rate = (
            f"exp({function_tag}({safe_mean_energy}))*1[m^3/s]"
        )
        lines.extend(
            [
                f"    // {reaction.name}: preserve the Electron Impact Reaction",
                "    // stoichiometry and energy loss, while replacing its rate.",
                f'    {reaction_target}.set("SpecifyReactionUsing", "RateConstant");',
                f'    {reaction_target}.set("RateConstantForm", "UseRate");',
                f'    {reaction_target}.set("UseTownsendTwoTermsBoltzmann", false);',
                "    // Swarm k is per target particle [m^3/s]; COMSOL kf is",
                "    // molar [m^3/(mol*s)], hence the explicit Avogadro factor.",
                f'    {reaction_target}.set("kf", '
                + java_string(
                    f"{AVOGADRO_PER_MOL}*{particle_rate}"
                )
                + ");",
            ]
        )
    lines.extend(
        [
            "    // Use one physically informed initial state, not a parameter",
            "    // or coefficient sweep. En_per is the log electron-energy",
            "    // density field, so this preserves its RF/spatial shape.",
            "    shiftPeriodicLogEnergySolution(model, "
            f"{mapping.run.initial_mean_energy_scale:.17g});",
        ]
    )
    lines.extend(
        [
            "    model.save("
            + java_string(java_path(mapping.model.output_mph))
            + ");",
            '    ModelUtil.remove("swarmGecApply");',
            "  }",
            "",
            "  private static void ensureInterpolation(Model model, String tag) {",
            "    try {",
            "      model.func(tag);",
            "    } catch (Exception ex) {",
            '      model.func().create(tag, "Interpolation");',
            "    }",
            "  }",
            "",
            "  private static void shiftPeriodicLogEnergySolution(",
            "      Model model, double meanEnergyScale) {",
            '    SolverSequence periodic = model.sol("sol1");',
            "    int solnum = periodic.getDefaultSolnum();",
            "    double[] values = periodic.getU(solnum);",
            '    SolverFeature variables = periodic.feature("v1");',
            "    XmeshInfo mesh = variables.xmeshInfo();",
            "    XmeshInfoDofs dofs = mesh.dofs();",
            "    String[] names = dofs.dofNames();",
            "    int[] nameIndices = dofs.nameInds();",
            "    int[] vectorIndices = dofs.solVectorInds();",
            "    double logScale = Math.log(meanEnergyScale);",
            "    int shifted = 0;",
            "    for (int dof = 0; dof < nameIndices.length; dof++) {",
            "      int nameIndex = nameIndices[dof];",
            "      int vectorIndex = vectorIndices[dof];",
            '      if (nameIndex >= 0 && nameIndex < names.length',
            '          && "comp1.En_per".equals(names[nameIndex])',
            "          && vectorIndex >= 0 && vectorIndex < values.length) {",
            "        values[vectorIndex] += logScale;",
            "        shifted++;",
            "      }",
            "    }",
            "    variables.clearXmesh();",
            "    if (shifted == 0) {",
            '      throw new IllegalStateException("En_per initial DOFs not found");',
            "    }",
            "    periodic.setU(solnum, values);",
            "    periodic.createSolution();",
            "  }",
            "}",
            "",
        ]
    )
    return "\n".join(lines)


def generate_run_java(
    mapping: GecCcpMapping,
    *,
    class_name: str,
    input_mph: Path,
    output_mph: Path,
    robust_nonlinear: bool = False,
) -> str:
    lines = _java_header(class_name, input_mph, class_name)
    lines.extend(
        [
            "    // Match the original model: one solve at its configured 1 W",
            "    // operating point, with no power or coefficient continuation.",
            *(
                [
                    "    // The imported closure is much stiffer in the EEDF tail.",
                    "    // Use COMSOL's cautious automatic Newton method for",
                    "    // strongly nonlinear problems on the refined tables.",
                    '    model.sol("sol1").feature("s1").feature("fc1")'
                    '.set("dtech", "hnlin");',
                    '    model.sol("sol1").feature("s1").feature("fc1")'
                    '.set("minsteph", 1.0e-12);',
                    '    model.sol("sol1").feature("s1").feature("fc1")'
                    '.set("useminsteprecovery", "off");',
                    '    model.sol("sol1").feature("s1").feature("fc1")'
                    '.set("maxiter", 200);',
                ]
                if robust_nonlinear
                else []
            ),
            "      model.param().set("
            + java_string(mapping.run.power_parameter)
            + f', "{mapping.run.power_W:.17g}[W]");',
            f"    model.study({java_string(mapping.model.study)}).run();",
            "    model.save(" + java_string(java_path(output_mph)) + ");",
        ]
    )
    lines.extend(
        [
            "    // Convert the converged periodic solution to physical RF time.",
            f"    model.study({java_string(mapping.model.conversion_study)}).run();",
            "    model.save(" + java_string(java_path(output_mph)) + ");",
            "    ModelUtil.remove(" + java_string(class_name) + ");",
            "  }",
            "}",
            "",
        ]
    )
    return "\n".join(lines)


def generate_export_java(
    mapping: GecCcpMapping,
    *,
    class_name: str,
    input_mph: Path,
    output_dir: Path,
) -> str:
    period_expr = (
        "ptp.neav",
        "ptp.Teav",
        "ptp.Vav",
        "ptp.Re_av",
    )
    period_units = ("1/m^3", "V", "V", "1/(m^3*s)")
    phase_expr = ("ptp.ne", "ptp.Te", "V")
    phase_units = ("1/m^3", "V", "V")
    waveform_expr = ("ptp.mct1.V", "ptp.mct1.I")
    waveform_units = ("V", "A")
    lines = _java_header(class_name, input_mph, class_name)
    lines.extend(
        [
            f"    model.result().dataset({java_string(mapping.run.axis_dataset)})"
            f'.set("data", {java_string(mapping.run.period_dataset)});',
            f"    model.result().dataset({java_string(mapping.run.radial_dataset)})"
            f'.set("data", {java_string(mapping.run.period_dataset)});',
        ]
    )
    lines.extend(
        _data_export_lines(
            "swAxisAvg",
            mapping.run.axis_dataset,
            period_expr,
            period_units,
            output_dir / "axis_period_average.csv",
        )
    )
    lines.extend(
        _data_export_lines(
            "swRadialAvg",
            mapping.run.radial_dataset,
            period_expr,
            period_units,
            output_dir / "radial_period_average.csv",
        )
    )
    lines.extend(
        [
            f"    model.result().dataset({java_string(mapping.run.axis_dataset)})"
            f'.set("data", {java_string(mapping.run.phase_dataset)});',
        ]
    )
    lines.extend(
        _data_export_lines(
            "swAxisPhase",
            mapping.run.axis_dataset,
            phase_expr,
            phase_units,
            output_dir / "axis_phase_resolved.csv",
        )
    )
    lines.extend(
        _data_export_lines(
            "swClosurePhase",
            mapping.run.axis_dataset,
            (
                "ptp.ebar",
                "ptp.Nn",
                "ptp.muerr",
                "ptp.Derr",
                "ptp.muenrr",
                "ptp.Denrr",
                "ptp.kf_1",
                "ptp.kf_2",
                "ptp.kf_3",
            ),
            (
                "V",
                "1/m^3",
                "m^2/(V*s)",
                "m^2/s",
                "m^2/(V*s)",
                "m^2/s",
                "m^3/(s*mol)",
                "m^3/(s*mol)",
                "m^3/(s*mol)",
            ),
            output_dir / "closure_phase.csv",
        )
    )
    lines.extend(
        _data_export_lines(
            "swWaveform",
            mapping.run.waveform_dataset,
            waveform_expr,
            waveform_units,
            output_dir / "electrode_waveform.csv",
        )
    )
    lines.extend(
        [
            "    ModelUtil.remove(" + java_string(class_name) + ");",
            "  }",
            "}",
            "",
        ]
    )
    return "\n".join(lines)


def execute_gec_ccp_run(
    mapping_path: str | Path,
    *,
    bundle_path: str | Path | None = None,
    comsol_executable: str | Path | None = None,
    reuse_baseline: bool = False,
) -> GecCcpRunSummary:
    plan = prepare_gec_ccp_run(
        mapping_path,
        bundle_path=bundle_path,
        write_java=True,
    )
    status_path = plan.output_directory / "gec_ccp_run_status.json"
    if reuse_baseline and not plan.mapping.model.baseline_output_mph.exists():
        raise GecCcpWorkflowError(
            "cannot reuse baseline because its output MPH does not exist: "
            f"{plan.mapping.model.baseline_output_mph}"
        )
    steps = (
        *(()
          if reuse_baseline
          else (("gec_baseline_run", plan.baseline_run_java),)),
        ("gec_baseline_export", plan.baseline_export_java),
        ("apply", plan.apply_java),
        ("gec_external_run", plan.external_run_java),
        ("gec_external_export", plan.external_export_java),
    )
    executions: list[ComsolExecutionSummary] = []
    try:
        for operation, java_file in steps:
            executions.append(
                execute_generated_comsol_java(
                    plan.mapping,
                    java_file,
                    operation=operation,
                    comsol_executable=comsol_executable,
                )
            )
    except ComsolAdapterError as exc:
        blocker_code = _comsol_blocker_code(exc)
        write_json(
            status_path,
            {
                "stage": "run-gec-ccp",
                "status": "blocked",
                "blocker_code": blocker_code,
                "reason": str(exc),
                "completed_operations": [
                    execution.operation for execution in executions
                ],
                "reuse_baseline": reuse_baseline,
                "coefficient_sweep": False,
                "model_contract": asdict(plan.contract),
                "plan": str(plan.plan_json),
            },
        )
        raise
    missing = [str(path) for path in plan.expected_result_files if not path.exists()]
    if missing:
        raise GecCcpWorkflowError(
            "COMSOL completed without all expected result files: "
            + ", ".join(missing)
        )
    write_json(
        status_path,
        {
            "stage": "run-gec-ccp",
            "status": "completed",
            "operations": [execution.operation for execution in executions],
            "reuse_baseline": reuse_baseline,
            "power_W": plan.mapping.run.power_W,
            "coefficient_sweep": False,
            "results": [str(path) for path in plan.expected_result_files],
        },
    )
    return GecCcpRunSummary(plan, tuple(executions), status_path)


def format_gec_ccp_plan(plan: GecCcpPlan) -> str:
    return "\n".join(
        [
            "GEC CCP COMSOL dry-run ready",
            f"model: {plan.mapping.model.input_mph}",
            f"detected EEDF: {plan.contract.original_eedf}",
            f"bundle: {plan.mapping.bundle.path}",
            "closure: COMSOL-solved mean energy + Swarm transport/rate lookups",
            f"plan: {plan.plan_json}",
            f"generated Java: {plan.output_directory}",
            "COMSOL command: not executed",
        ]
    )


def format_gec_ccp_summary(summary: GecCcpRunSummary) -> str:
    return "\n".join(
        [
            "GEC CCP COMSOL workflow completed",
            f"status: {summary.status_json}",
            f"results: {summary.plan.output_directory}",
        ]
    )


def _validate_contract(mapping: GecCcpMapping, contract: MphContract) -> None:
    expected = mapping.model
    failures: list[str] = []
    if contract.physics_operation != "ColdPlasmaTimePeriodic":
        failures.append(f"physics operation={contract.physics_operation}")
    if contract.physics_tag != expected.physics:
        failures.append(f"physics tag={contract.physics_tag}")
    if contract.plasma_feature != expected.plasma_feature:
        failures.append(f"plasma feature={contract.plasma_feature}")
    if contract.original_eedf != expected.expected_original_eedf:
        failures.append(
            f"EEDF={contract.original_eedf}, expected={expected.expected_original_eedf}"
        )
    for reaction in mapping.reactions:
        if reaction.feature not in contract.reaction_features:
            failures.append(f"missing reaction={reaction.feature}")
    for study in (expected.study, expected.conversion_study):
        if study not in contract.studies:
            failures.append(f"missing study={study}")
    for dataset in (
        mapping.run.axis_dataset,
        mapping.run.radial_dataset,
        mapping.run.period_dataset,
        mapping.run.phase_dataset,
        mapping.run.waveform_dataset,
    ):
        if dataset not in contract.datasets:
            failures.append(f"missing dataset={dataset}")
    if failures:
        raise GecCcpWorkflowError(
            "local MPH does not match the GEC CCP mapping: " + "; ".join(failures)
        )
    if expected.input_mph in {
        expected.baseline_output_mph,
        expected.output_mph,
    }:
        raise GecCcpWorkflowError("output MPH paths must not overwrite the input MPH")


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


def _data_export_lines(
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


def _mapping(raw: dict[str, Any], name: str) -> dict[str, Any]:
    value = raw.get(name)
    if not isinstance(value, dict):
        raise GecCcpWorkflowError(f"{name} must be a mapping")
    return value


def _text(raw: dict[str, Any], name: str) -> str:
    value = raw.get(name)
    if not isinstance(value, str) or not value.strip():
        raise GecCcpWorkflowError(f"{name} must be a non-empty string")
    return value.strip()


def _resolve(root: Path, value: object, name: str) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise GecCcpWorkflowError(f"{name} must be a path")
    path = Path(value)
    return (root / path).resolve() if not path.is_absolute() else path.resolve()


def _property_values(xml: str, name: str) -> list[str]:
    matches = re.findall(
        rf'<param\b[^>]*param="{re.escape(name)}"[^>]*value="([^"]+)"',
        xml,
    )
    if not matches:
        raise GecCcpWorkflowError(f"model property not found: {name}")
    values: list[str] = []
    for encoded in matches:
        quoted = re.search(r"'([^']+)'", encoded)
        values.append(quoted.group(1) if quoted else encoded)
    return values


def _tag_for_operation(xml: str, element: str, operation: str) -> str:
    match = re.search(
        rf'<{element}\b[^>]*op="{re.escape(operation)}"[^>]*tag="([^"]+)"',
        xml,
    )
    if match is None:
        raise GecCcpWorkflowError(f"model operation not found: {operation}")
    return match.group(1)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _lookup(
    rows: Iterable[dict[str, str]],
    x_column: str,
    y_column: str,
) -> tuple[list[float], list[float]]:
    values: dict[float, float] = {}
    for row in rows:
        try:
            x = float(row[x_column])
            y = float(row[y_column])
        except (KeyError, ValueError) as exc:
            raise GecCcpWorkflowError(
                f"invalid lookup columns {x_column}, {y_column}"
            ) from exc
        if not math.isfinite(x) or not math.isfinite(y) or y < 0.0:
            raise GecCcpWorkflowError(f"invalid lookup value in {y_column}")
        values[x] = y
    if len(values) < 2:
        raise GecCcpWorkflowError(f"lookup {y_column} needs at least two points")
    x_values = sorted(values)
    return x_values, [values[x] for x in x_values]


def _array(values: Iterable[float]) -> str:
    return ", ".join(f"{value:.17e}" for value in values)


def _inline_interpolation_lines(
    tag: str,
    x_values: list[float],
    y_values: list[float],
    *,
    argunit: str,
    fununit: str,
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
        f"    ensureInterpolation(model, {qtag});",
        f'    model.func({qtag}).set("source", "table");',
        f'    model.func({qtag}).set("nargs", 1);',
        f'    model.func({qtag}).set("argunit", {java_string(argunit)});',
        f'    model.func({qtag}).set("fununit", {java_string(fununit)});',
        f'    model.func({qtag}).set("interp", "piecewisecubic");',
        f'    model.func({qtag}).set("extrap", "const");',
        f'    model.func({qtag}).set("table", new String[][]{{{table}}});',
    ]


def _guard_mean_energy_lookup(
    x_values: list[float],
    y_values: list[float],
    *,
    low_value: float,
) -> tuple[list[float], list[float]]:
    """Keep COMSOL Newton iterates finite outside the physical table range.

    The reported/validated physical range remains the Swarm range.  The two
    low-side points and one high-side point only prevent COMSOL's internal
    interpolators from returning zero during nonphysical Newton iterates.
    """

    guard = max(1.0e6, 10.0 * x_values[-1])
    return (
        [-guard, 0.0, *x_values, guard],
        [low_value, low_value, *y_values, y_values[-1]],
    )


def _comsol_blocker_code(exc: ComsolAdapterError) -> str:
    result = exc.step_result or {}
    text_parts: list[str] = []
    for field in ("stdout", "stderr"):
        value = result.get(field)
        if not isinstance(value, str):
            continue
        path = Path(value)
        if path.exists():
            text_parts.append(
                path.read_text(encoding="utf-8", errors="replace").lower()
            )
    text = "\n".join(text_parts)
    if "license has expired" in text or "product has expired" in text:
        return "comsol_license_expired"
    return "comsol_execution_failed"


def _positive_numeric_values(
    rows: Iterable[dict[str, str]],
    column: str,
) -> list[float]:
    values = _nonnegative_numeric_values(rows, column)
    if any(value <= 0.0 for value in values):
        raise GecCcpWorkflowError(f"bundle column {column} must be positive")
    return values


def _nonnegative_numeric_values(
    rows: Iterable[dict[str, str]],
    column: str,
) -> list[float]:
    values: list[float] = []
    for row in rows:
        try:
            value = float(row[column])
        except (KeyError, TypeError, ValueError) as exc:
            raise GecCcpWorkflowError(
                f"bundle column {column} contains a nonnumeric value"
            ) from exc
        if not math.isfinite(value) or value < 0.0:
            raise GecCcpWorkflowError(
                f"bundle column {column} contains an invalid value"
            )
        values.append(value)
    if not values:
        raise GecCcpWorkflowError(f"bundle column {column} has no values")
    return values
