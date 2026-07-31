"""Parse and validate explicit COMSOL model mapping YAML files."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

import yaml

from ._paths import discover_repo_root, resolve_path


class ComsolMappingError(ValueError):
    """Raised when a COMSOL mapping file is invalid or incomplete."""


@dataclass(frozen=True, slots=True)
class ModelSpec:
    input_mph: Path
    output_mph: Path
    study: str
    component: str
    physics: str


@dataclass(frozen=True, slots=True)
class BundleSpec:
    path: Path


@dataclass(frozen=True, slots=True)
class LogSpec:
    path: Path


@dataclass(frozen=True, slots=True)
class VerifySpec:
    output_path: Path
    relative_tolerance: float
    absolute_tolerance: float


@dataclass(frozen=True, slots=True)
class FunctionMapping:
    name: str
    tag: str
    file: str
    path: Path
    column: str | None
    nargs: int
    argunit: str
    fununit: str
    interp: str
    extrap: str


@dataclass(frozen=True, slots=True)
class MeanEnergyLookup:
    table: str
    path: Path
    x_column: str
    y_column: str


@dataclass(frozen=True, slots=True)
class MeanEnergyFormulation:
    """Canonical mean-energy mode and its exact COMSOL property contract."""

    mode: str
    property_group: str
    property: str
    comsol_value: str


@dataclass(frozen=True, slots=True)
class ClosureMapping:
    physics: str
    feature: str
    mean_energy_formulation: MeanEnergyFormulation
    mean_energy: MeanEnergyLookup


@dataclass(frozen=True, slots=True)
class ReactionLookupMapping:
    name: str
    physics: str
    feature: str
    form: str
    source_model: str
    table: str
    path: Path
    argument: str
    x_column: str
    y_column: str
    process_type: str


@dataclass(frozen=True, slots=True)
class RunMapping:
    voltage_parameter: str
    voltage_feature: str
    voltage_property: str
    voltages_V: tuple[float, ...]


@dataclass(frozen=True, slots=True)
class ResultProbe:
    name: str
    expression: str
    unit: str


@dataclass(frozen=True, slots=True)
class ComsolModelMapping:
    path: Path
    root: Path
    model: ModelSpec
    bundle: BundleSpec
    logs: LogSpec
    verify: VerifySpec
    functions: tuple[FunctionMapping, ...]
    closure: ClosureMapping
    reaction_lookups: tuple[ReactionLookupMapping, ...]
    run: RunMapping
    probes: tuple[ResultProbe, ...]


CLOSURE_QUANTITIES = (
    "mean_energy",
    "reduced_mobility",
    "reduced_longitudinal_diffusion",
    "reduced_energy_mobility",
    "reduced_energy_diffusion",
    "excitation_townsend",
    "ionization_townsend",
)

_MEAN_ENERGY_COMSOL_VALUES = {
    "local_energy": "LocalEnergyApproximationE",
    "local_field": "LocalFieldApproximationE",
}

_MEAN_ENERGY_ACTIVITY = {
    "local_energy": {
        "mean_energy": (
            "inactive",
            "COMSOL solves the mean-electron-energy equation; the E/N-to-mean-energy "
            "table is written for audit but does not prescribe the solved mean energy.",
        ),
        "reduced_mobility": (
            "active",
            "The lookup supplies electron mobility to the solved drift-diffusion system.",
        ),
        "reduced_longitudinal_diffusion": (
            "active",
            "Swarm longitudinal diffusion is mapped to COMSOL's scalar electron "
            "diffusivity lookup used by the 1D drift-diffusion system.",
        ),
        "reduced_energy_mobility": (
            "active",
            "The lookup supplies energy mobility to the solved mean-energy equation.",
        ),
        "reduced_energy_diffusion": (
            "active",
            "The lookup supplies energy diffusivity to the solved mean-energy equation.",
        ),
        "excitation_townsend": (
            "active",
            "The eir2 reaction uses the written Townsend lookup.",
        ),
        "ionization_townsend": (
            "active",
            "The eir4 reaction uses the written Townsend lookup.",
        ),
    },
    "local_field": {
        "mean_energy": (
            "active",
            "COMSOL prescribes mean electron energy from the E/N lookup under the "
            "local-field approximation.",
        ),
        "reduced_mobility": (
            "active",
            "The lookup supplies electron mobility to the solved drift-diffusion system.",
        ),
        "reduced_longitudinal_diffusion": (
            "active",
            "Swarm longitudinal diffusion is mapped to COMSOL's scalar electron "
            "diffusivity lookup used by the 1D drift-diffusion system.",
        ),
        "reduced_energy_mobility": (
            "inactive",
            "The local-field approximation does not solve the mean-energy equation, "
            "so its energy-flux mobility is written for audit but not used.",
        ),
        "reduced_energy_diffusion": (
            "inactive",
            "The local-field approximation does not solve the mean-energy equation, "
            "so its energy-flux diffusivity is written for audit but not used.",
        ),
        "excitation_townsend": (
            "active",
            "The eir2 reaction uses the written Townsend lookup.",
        ),
        "ionization_townsend": (
            "active",
            "The eir4 reaction uses the written Townsend lookup.",
        ),
    },
}


def closure_quantity_activity(
    formulation: MeanEnergyFormulation,
) -> dict[str, dict[str, str]]:
    """Describe which written closure properties enter equations for this mode."""

    contract = _MEAN_ENERGY_ACTIVITY[formulation.mode]
    return {
        quantity: {"status": status, "reason": reason}
        for quantity, (status, reason) in contract.items()
    }


def load_comsol_mapping(
    path: str | Path,
    *,
    validate_files: bool = False,
) -> ComsolModelMapping:
    mapping_path = Path(path).resolve()
    root = discover_repo_root(mapping_path)
    data = _read_yaml_mapping(mapping_path)
    allowed = {
        "model",
        "bundle",
        "logs",
        "verify",
        "functions",
        "closure",
        "reaction_lookups",
        "run",
        "results",
    }
    unknown = sorted(set(data) - allowed)
    if unknown:
        raise ComsolMappingError(
            "unsupported COMSOL mapping field(s): " + ", ".join(unknown)
        )
    model = _parse_model(_required_mapping(data, "model"), root)
    bundle = _parse_bundle(_required_mapping(data, "bundle"), root)
    logs = _parse_logs(data.get("logs", {}), root)
    verify = _parse_verify(data.get("verify", {}), root)
    functions = _parse_functions(_required_mapping(data, "functions"), bundle)
    closure = _parse_closure(_required_mapping(data, "closure"), model, bundle)
    reaction_lookups = _parse_reaction_lookups(
        data.get("reaction_lookups"),
        model,
        bundle,
    )
    run = _parse_run(_required_mapping(data, "run"))
    probes = _parse_results(data.get("results", {}))
    mapping = ComsolModelMapping(
        path=mapping_path,
        root=root,
        model=model,
        bundle=bundle,
        logs=logs,
        verify=verify,
        functions=functions,
        closure=closure,
        reaction_lookups=reaction_lookups,
        run=run,
        probes=probes,
    )
    if validate_files:
        validate_comsol_mapping_files(mapping)
    return mapping


def validate_comsol_mapping_files(
    mapping: ComsolModelMapping,
    *,
    require_unit_metadata: bool = False,
) -> None:
    if not mapping.model.input_mph.exists():
        raise ComsolMappingError(
            f"COMSOL input .mph does not exist: {mapping.model.input_mph}. "
            + _missing_input_mph_hint(mapping.model.input_mph)
        )
    if not mapping.bundle.path.exists() or not mapping.bundle.path.is_dir():
        raise ComsolMappingError(
            f"COMSOL table bundle does not exist: {mapping.bundle.path}. "
            + _missing_bundle_hint(mapping.bundle.path)
        )
    manifest_path = mapping.bundle.path / "manifest.json"
    if not manifest_path.exists():
        raise ComsolMappingError(
            f"COMSOL table bundle manifest.json does not exist: {manifest_path}. "
            + _missing_bundle_hint(mapping.bundle.path)
        )
    manifest = _read_manifest(manifest_path)
    tables = manifest.get("tables", {})
    if not isinstance(tables, dict):
        tables = {}

    for function in mapping.functions:
        if not function.path.exists():
            raise ComsolMappingError(
                f"mapped function CSV does not exist for {function.name}: "
                f"{function.path}"
            )
        header = _csv_header(function.path)
        if function.column is not None and function.column not in header:
            raise ComsolMappingError(
                f"mapped function column {function.column!r} not found in "
                f"{function.file}"
            )
        table_meta = tables.get(function.file)
        if (
            function.column is not None
            and isinstance(table_meta, dict)
            and isinstance(table_meta.get("columns"), list)
            and function.column not in table_meta["columns"]
        ):
            raise ComsolMappingError(
                f"mapped function column {function.column!r} is not listed in "
                f"manifest table {function.file}"
            )
        if require_unit_metadata:
            _validate_function_units(function, table_meta, manifest)

    for lookup in mapping.reaction_lookups:
        if not lookup.path.exists():
            raise ComsolMappingError(
                f"reaction lookup table does not exist for {lookup.name}: "
                f"{lookup.path}"
            )
        header = _csv_header(lookup.path)
        missing = [
            column
            for column in ("process_type", lookup.x_column, lookup.y_column)
            if column not in header
        ]
        if missing:
            raise ComsolMappingError(
                f"reaction lookup {lookup.name} table {lookup.table} is missing "
                + ", ".join(missing)
            )
        if not _csv_has_value(lookup.path, "process_type", lookup.process_type):
            raise ComsolMappingError(
                f"reaction lookup {lookup.name} process_type "
                f"{lookup.process_type!r} not found in {lookup.table}"
            )

    table = mapping.closure.mean_energy
    if not table.path.exists():
        raise ComsolMappingError(
            f"mean-energy lookup table does not exist: {table.path}"
        )
    header = _csv_header(table.path)
    missing = [column for column in (table.x_column, table.y_column) if column not in header]
    if missing:
        raise ComsolMappingError(
            f"mean-energy lookup table {table.table} is missing " + ", ".join(missing)
        )
    transport_path = mapping.bundle.path / "transport_vs_mean_energy.csv"
    transport_columns = {
        "mean_energy_eV",
        "reduced_mobility_m2_V_s_m3",
        "reduced_diffusion_L_m2_s_m3",
        "reduced_electron_energy_mobility_m2_V_s_m3",
        "reduced_electron_energy_diffusion_m2_s_m3",
    }
    if not transport_path.exists():
        raise ComsolMappingError(f"transport lookup table does not exist: {transport_path}")
    missing_transport = sorted(transport_columns - set(_csv_header(transport_path)))
    if missing_transport:
        raise ComsolMappingError(
            "transport lookup table is missing " + ", ".join(missing_transport)
        )


def dry_run_mapping_summary(mapping: ComsolModelMapping) -> dict[str, Any]:
    return {
        "mapping_path": str(mapping.path),
        "root": str(mapping.root),
        "model": {
            "input_mph": str(mapping.model.input_mph),
            "output_mph": str(mapping.model.output_mph),
            "study": mapping.model.study,
            "component": mapping.model.component,
            "physics": mapping.model.physics,
        },
        "bundle": {"path": str(mapping.bundle.path)},
        "logs": {"path": str(mapping.logs.path)},
        "verify": {
            "output_path": str(mapping.verify.output_path),
            "relative_tolerance": mapping.verify.relative_tolerance,
            "absolute_tolerance": mapping.verify.absolute_tolerance,
        },
        "functions": [
            {
                "name": function.name,
                "tag": function.tag,
                "file": function.file,
                "path": str(function.path),
                "column": function.column,
                "nargs": function.nargs,
                "argunit": function.argunit,
                "fununit": function.fununit,
                "interp": function.interp,
                "extrap": function.extrap,
            }
            for function in mapping.functions
        ],
        "closure": {
            "physics": mapping.closure.physics,
            "feature": mapping.closure.feature,
            "mean_energy_formulation": {
                "mode": mapping.closure.mean_energy_formulation.mode,
                "property_group": (
                    mapping.closure.mean_energy_formulation.property_group
                ),
                "property": mapping.closure.mean_energy_formulation.property,
                "comsol_value": (
                    mapping.closure.mean_energy_formulation.comsol_value
                ),
                "quantity_activity": closure_quantity_activity(
                    mapping.closure.mean_energy_formulation
                ),
            },
            "mean_energy": {
                "table": mapping.closure.mean_energy.table,
                "path": str(mapping.closure.mean_energy.path),
                "x_column": mapping.closure.mean_energy.x_column,
                "y_column": mapping.closure.mean_energy.y_column,
            },
        },
        "reaction_lookups": [
            {
                "name": lookup.name,
                "physics": lookup.physics,
                "feature": lookup.feature,
                "form": lookup.form,
                "source_model": lookup.source_model,
                "table": lookup.table,
                "path": str(lookup.path),
                "argument": lookup.argument,
                "x_column": lookup.x_column,
                "y_column": lookup.y_column,
                "process_type": lookup.process_type,
            }
            for lookup in mapping.reaction_lookups
        ],
        "run": {
            "voltage_parameter": mapping.run.voltage_parameter,
            "voltage_feature": mapping.run.voltage_feature,
            "voltage_property": mapping.run.voltage_property,
            "voltages_V": list(mapping.run.voltages_V),
        },
        "probes": [
            {
                "name": probe.name,
                "expression": probe.expression,
                "unit": probe.unit,
            }
            for probe in mapping.probes
        ],
    }


def _missing_input_mph_hint(path: Path) -> str:
    return (
        "Check model.input_mph in the mapping YAML. Use an unchanged source MPH; "
        "the workflow writes all modifications to model.output_mph."
    )


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


def _read_yaml_mapping(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise ComsolMappingError(f"mapping YAML does not exist: {path}")
    with path.open("r", encoding="utf-8") as fp:
        data = yaml.safe_load(fp) or {}
    if not isinstance(data, dict):
        raise ComsolMappingError("mapping YAML root must be a mapping")
    return data


def _parse_model(raw: dict[str, Any], root: Path) -> ModelSpec:
    return ModelSpec(
        input_mph=resolve_path(_required_str(raw, "input_mph"), root),
        output_mph=resolve_path(_required_str(raw, "output_mph"), root),
        study=_required_str(raw, "study"),
        component=_required_str(raw, "component"),
        physics=_required_str(raw, "physics"),
    )


def _parse_bundle(raw: dict[str, Any], root: Path) -> BundleSpec:
    return BundleSpec(path=resolve_path(_required_str(raw, "path"), root))


def _parse_logs(raw: object, root: Path) -> LogSpec:
    if raw is None:
        raw = {}
    if not isinstance(raw, dict):
        raise ComsolMappingError("logs must be a mapping")
    path_value = raw.get("path", "model/logs")
    if not isinstance(path_value, str) or not path_value:
        raise ComsolMappingError("logs.path must be a non-empty string")
    return LogSpec(path=resolve_path(path_value, root))


def _parse_verify(raw: object, root: Path) -> VerifySpec:
    if raw is None:
        raw = {}
    if not isinstance(raw, dict):
        raise ComsolMappingError("verify must be a mapping")
    path_value = raw.get("output_path", "outputs/comsol_verify")
    if not isinstance(path_value, str) or not path_value:
        raise ComsolMappingError("verify.output_path must be a non-empty string")
    return VerifySpec(
        output_path=resolve_path(path_value, root),
        relative_tolerance=_optional_nonnegative_float(
            raw,
            "relative_tolerance",
            1.0e-6,
        ),
        absolute_tolerance=_optional_nonnegative_float(
            raw,
            "absolute_tolerance",
            1.0e-12,
        ),
    )


def _parse_functions(
    raw: dict[str, Any],
    bundle: BundleSpec,
) -> tuple[FunctionMapping, ...]:
    if not raw:
        raise ComsolMappingError("functions must contain at least one mapping")
    functions: list[FunctionMapping] = []
    seen_tags: set[str] = set()
    for name, item in raw.items():
        if not isinstance(name, str) or not name:
            raise ComsolMappingError("function names must be non-empty strings")
        if not isinstance(item, dict):
            raise ComsolMappingError(f"functions.{name} must be a mapping")
        tag = _required_str(item, "tag")
        if tag in seen_tags:
            raise ComsolMappingError(f"duplicate function tag: {tag}")
        seen_tags.add(tag)
        file_name = _required_str(item, "file")
        column = item.get("column")
        if column is not None and not isinstance(column, str):
            raise ComsolMappingError(f"functions.{name}.column must be a string")
        nargs = _required_int(item, "nargs")
        if nargs <= 0:
            raise ComsolMappingError(f"functions.{name}.nargs must be positive")
        path = _resolve_bundle_file(bundle.path, file_name, function_name=name)
        functions.append(
            FunctionMapping(
                name=name,
                tag=tag,
                file=file_name,
                path=path,
                column=column,
                nargs=nargs,
                argunit=_required_str(item, "argunit"),
                fununit=_required_str(item, "fununit"),
                interp=_required_str(item, "interp"),
                extrap=_required_str(item, "extrap"),
            )
        )
    return tuple(functions)


def _parse_closure(
    raw: dict[str, Any],
    model: ModelSpec,
    bundle: BundleSpec,
) -> ClosureMapping:
    formulation_raw = _required_mapping(raw, "mean_energy_formulation")
    mode = _required_str(formulation_raw, "mode")
    if mode not in _MEAN_ENERGY_COMSOL_VALUES:
        raise ComsolMappingError(
            "closure.mean_energy_formulation.mode must be 'local_energy' "
            "or 'local_field'"
        )
    property_group = _required_str(formulation_raw, "property_group")
    property_name = _required_str(formulation_raw, "property")
    comsol_value = _required_str(formulation_raw, "comsol_value")
    if property_group != "ElectronProperties":
        raise ComsolMappingError(
            "closure.mean_energy_formulation.property_group must be "
            "'ElectronProperties'"
        )
    if property_name != "MeanElectronEnergyModel":
        raise ComsolMappingError(
            "closure.mean_energy_formulation.property must be "
            "'MeanElectronEnergyModel'"
        )
    expected_value = _MEAN_ENERGY_COMSOL_VALUES[mode]
    if comsol_value != expected_value:
        raise ComsolMappingError(
            "closure.mean_energy_formulation.comsol_value is inconsistent with "
            f"mode {mode!r}; expected {expected_value!r}"
        )
    mean_energy = _required_mapping(raw, "mean_energy")
    table = _required_str(mean_energy, "table")
    return ClosureMapping(
        physics=_optional_str(raw, "physics", model.physics),
        feature=_required_str(raw, "feature"),
        mean_energy_formulation=MeanEnergyFormulation(
            mode=mode,
            property_group=property_group,
            property=property_name,
            comsol_value=comsol_value,
        ),
        mean_energy=MeanEnergyLookup(
            table=table,
            path=_resolve_bundle_file(bundle.path, table, function_name="mean_energy"),
            x_column=_optional_str(mean_energy, "x_column", "E_over_N_Td"),
            y_column=_optional_str(mean_energy, "y_column", "mean_energy_eV"),
        ),
    )


def _parse_reaction_lookups(
    raw: object,
    model: ModelSpec,
    bundle: BundleSpec,
) -> tuple[ReactionLookupMapping, ...]:
    if raw is None:
        return ()
    if not isinstance(raw, dict):
        raise ComsolMappingError("reaction_lookups must be a mapping")
    lookups = []
    for name, item in raw.items():
        if not isinstance(name, str) or not name:
            raise ComsolMappingError("reaction lookup names must be non-empty strings")
        if not isinstance(item, dict):
            raise ComsolMappingError(f"reaction_lookups.{name} must be a mapping")
        form = _required_str(item, "form")
        if form != "townsend":
            raise ComsolMappingError(
                f"reaction_lookups.{name}.form must be 'townsend'"
            )
        source_model = _optional_str(item, "source_model", "townsend_flux")
        if source_model != "townsend_flux":
            raise ComsolMappingError(
                f"reaction_lookups.{name}.source_model must be 'townsend_flux'"
            )
        argument = _optional_str(item, "argument", "mean_energy_eV")
        table = _optional_str(
            item,
            "table",
            "rates_vs_mean_energy.csv"
            if argument == "mean_energy_eV"
            else "rates_vs_en.csv",
        )
        x_column = _optional_str(item, "x_column", argument)
        y_column = _optional_str(item, "y_column", "reduced_townsend_m2")
        if y_column not in {
            "reduced_townsend_m2",
            "mixture_weighted_reduced_townsend_m2",
        }:
            raise ComsolMappingError(
                f"reaction_lookups.{name}.y_column must be a Townsend column "
                "when form is 'townsend'; use reduced_townsend_m2 or "
                "mixture_weighted_reduced_townsend_m2"
            )
        path = _resolve_bundle_file(bundle.path, table, function_name=name)
        lookups.append(
            ReactionLookupMapping(
                name=name,
                physics=_optional_str(item, "physics", model.physics),
                feature=_required_str(item, "feature"),
                form=form,
                source_model=source_model,
                table=table,
                path=path,
                argument=argument,
                x_column=x_column,
                y_column=y_column,
                process_type=_required_str(item, "process_type"),
            )
        )
    return tuple(lookups)


def _parse_run(raw: dict[str, Any]) -> RunMapping:
    values = raw.get("voltages_V")
    if not isinstance(values, list) or not values:
        raise ComsolMappingError("run.voltages_V must be a non-empty list")
    voltages: list[float] = []
    for value in values:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ComsolMappingError("run.voltages_V values must be numbers")
        number = float(value)
        if number <= 0.0:
            raise ComsolMappingError("run.voltages_V values must be positive")
        voltages.append(number)
    if any(right <= left for left, right in zip(voltages, voltages[1:])):
        raise ComsolMappingError("run.voltages_V must be strictly increasing")
    return RunMapping(
        voltage_parameter=_optional_str(raw, "voltage_parameter", "V0"),
        voltage_feature=_required_str(raw, "voltage_feature"),
        voltage_property=_optional_str(raw, "voltage_property", "V0"),
        voltages_V=tuple(voltages),
    )


def _parse_results(raw: object) -> tuple[ResultProbe, ...]:
    if raw is None:
        return ()
    if not isinstance(raw, dict):
        raise ComsolMappingError("results must be a mapping")
    probes = raw.get("probes", [])
    if probes is None:
        probes = []
    if not isinstance(probes, list):
        raise ComsolMappingError("results.probes must be a list")
    parsed = []
    for index, item in enumerate(probes):
        if not isinstance(item, dict):
            raise ComsolMappingError(f"results.probes[{index}] must be a mapping")
        parsed.append(
            ResultProbe(
                name=_required_str(item, "name"),
                expression=_required_str(item, "expression"),
                unit=_required_str(item, "unit"),
            )
        )
    return tuple(parsed)


def _required_mapping(raw: dict[str, Any], key: str) -> dict[str, Any]:
    value = raw.get(key)
    if not isinstance(value, dict):
        raise ComsolMappingError(f"{key} must be a mapping")
    return value


def _required_str(raw: dict[str, Any], key: str) -> str:
    value = raw.get(key)
    if not isinstance(value, str) or not value:
        raise ComsolMappingError(f"{key} must be a non-empty string")
    return value


def _optional_str(raw: dict[str, Any], key: str, default: str) -> str:
    value = raw.get(key, default)
    if not isinstance(value, str) or not value:
        raise ComsolMappingError(f"{key} must be a non-empty string")
    return value


def _required_int(raw: dict[str, Any], key: str) -> int:
    value = raw.get(key)
    if not isinstance(value, int):
        raise ComsolMappingError(f"{key} must be an integer")
    return value


def _optional_nonnegative_float(
    raw: dict[str, Any],
    key: str,
    default: float,
) -> float:
    value = raw.get(key, default)
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise ComsolMappingError(f"verify.{key} must be a number")
    number = float(value)
    if number < 0.0:
        raise ComsolMappingError(f"verify.{key} must be non-negative")
    return number


def _resolve_bundle_file(
    bundle_path: Path,
    file_name: str,
    *,
    function_name: str,
) -> Path:
    path = Path(file_name)
    if path.is_absolute():
        raise ComsolMappingError(
            f"functions.{function_name}.file must be relative to bundle.path"
        )
    resolved = (bundle_path / path).resolve()
    try:
        resolved.relative_to(bundle_path.resolve())
    except ValueError as exc:
        raise ComsolMappingError(
            f"functions.{function_name}.file must stay under bundle.path"
        ) from exc
    return resolved


def _read_manifest(path: Path) -> dict[str, Any]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ComsolMappingError(f"invalid bundle manifest JSON: {path}") from exc
    if not isinstance(data, dict):
        raise ComsolMappingError(f"bundle manifest must be a mapping: {path}")
    return data


def _csv_header(path: Path) -> tuple[str, ...]:
    with path.open("r", encoding="utf-8", newline="") as fp:
        try:
            header = next(csv.reader(fp))
        except StopIteration as exc:
            raise ComsolMappingError(f"mapped function CSV is empty: {path}") from exc
    return tuple(header)


def _csv_has_value(path: Path, column: str, expected: str) -> bool:
    with path.open("r", encoding="utf-8", newline="") as fp:
        for row in csv.DictReader(fp):
            if row.get(column) == expected:
                return True
    return False


def _validate_function_units(
    function: FunctionMapping,
    table_meta: object,
    manifest: dict[str, Any],
) -> None:
    table_units: object = {}
    if isinstance(table_meta, dict):
        table_units = table_meta.get("units", {})
    manifest_units = manifest.get("units", {})
    unit_columns: set[str] = set()
    if isinstance(table_units, dict):
        unit_columns.update(str(key) for key in table_units)
    if isinstance(manifest_units, dict):
        unit_columns.update(str(key) for key in manifest_units)
    if not unit_columns:
        raise ComsolMappingError(
            f"mapped function {function.name} has no unit metadata in "
            f"bundle manifest for {function.file}"
        )
    if function.column is not None and function.column not in unit_columns:
        raise ComsolMappingError(
            f"mapped function column {function.column!r} has no unit metadata "
            f"in bundle manifest for {function.file}"
        )
