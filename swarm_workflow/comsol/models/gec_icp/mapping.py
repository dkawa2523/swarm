"""Load the strict schema-v2 mapping for the GEC ICP MPH contract."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from .contracts import (
    GecIcpContractError,
    GecIcpMapping,
    GecIcpModelSpec,
    GecIcpParameterTags,
    GecIcpReactionTags,
    GecIcpSpeciesTags,
)


ROOT_FIELDS = {"schema_version", "model", "parameters", "species", "reactions"}
MODEL_FIELDS = {
    "input_mph",
    "output_mph",
    "component",
    "geometry",
    "mesh",
    "geometry_dimension",
    "axisymmetric",
    "plasma_physics",
    "plasma_feature",
    "magnetic_physics",
    "coil_feature",
    "plasma_conductivity_coupling",
    "electron_heat_source_coupling",
    "study",
    "study_feature",
    "solution",
    "dataset",
}
PARAMETER_FIELDS = {"power", "gas_temperature", "pressure"}
SPECIES_FIELDS = {"ground", "excited", "ion"}
REACTION_FIELDS = {
    "elastic",
    "excitation",
    "superelastic",
    "ionization",
    "stepwise_ionization",
}


def load_gec_icp_mapping(
    mapping_path: str | Path,
    *,
    validate_files: bool = False,
) -> GecIcpMapping:
    """Load one model-local mapping; relative paths use its directory."""

    path = Path(mapping_path).resolve()
    raw = _read_yaml_mapping(path)
    _reject_unknown(raw, ROOT_FIELDS, "mapping")
    schema_version = raw.get("schema_version")
    if isinstance(schema_version, bool) or schema_version != 2:
        raise GecIcpContractError("GEC ICP mapping requires schema_version: 2")

    root = path.parent
    model_raw = _required_mapping(raw, "model", "mapping")
    parameters_raw = _required_mapping(raw, "parameters", "mapping")
    species_raw = _required_mapping(raw, "species", "mapping")
    reactions_raw = _required_mapping(raw, "reactions", "mapping")
    _reject_unknown(model_raw, MODEL_FIELDS, "model")
    _reject_unknown(parameters_raw, PARAMETER_FIELDS, "parameters")
    _reject_unknown(species_raw, SPECIES_FIELDS, "species")
    _reject_unknown(reactions_raw, REACTION_FIELDS, "reactions")

    input_mph = _resolved_mph_path(model_raw, "input_mph", root)
    output_mph = _resolved_mph_path(model_raw, "output_mph", root)
    if input_mph == output_mph:
        raise GecIcpContractError("model.output_mph must differ from model.input_mph")
    dimension = _required_integer(model_raw, "geometry_dimension", "model")
    if dimension != 2:
        raise GecIcpContractError("model.geometry_dimension must be 2")
    axisymmetric = _required_boolean(model_raw, "axisymmetric", "model")
    if not axisymmetric:
        raise GecIcpContractError("model.axisymmetric must be true")

    mapping = GecIcpMapping(
        path=path,
        root=root,
        model=GecIcpModelSpec(
            input_mph=input_mph,
            output_mph=output_mph,
            component=_required_string(model_raw, "component", "model"),
            geometry=_required_string(model_raw, "geometry", "model"),
            mesh=_required_string(model_raw, "mesh", "model"),
            geometry_dimension=dimension,
            axisymmetric=axisymmetric,
            plasma_physics=_required_string(model_raw, "plasma_physics", "model"),
            plasma_feature=_required_string(model_raw, "plasma_feature", "model"),
            magnetic_physics=_required_string(model_raw, "magnetic_physics", "model"),
            coil_feature=_required_string(model_raw, "coil_feature", "model"),
            plasma_conductivity_coupling=_required_string(
                model_raw, "plasma_conductivity_coupling", "model"
            ),
            electron_heat_source_coupling=_required_string(
                model_raw, "electron_heat_source_coupling", "model"
            ),
            study=_required_string(model_raw, "study", "model"),
            study_feature=_required_string(model_raw, "study_feature", "model"),
            solution=_required_string(model_raw, "solution", "model"),
            dataset=_required_string(model_raw, "dataset", "model"),
        ),
        parameters=GecIcpParameterTags(
            power=_required_string(parameters_raw, "power", "parameters"),
            gas_temperature=_required_string(
                parameters_raw, "gas_temperature", "parameters"
            ),
            pressure=_required_string(parameters_raw, "pressure", "parameters"),
        ),
        species=GecIcpSpeciesTags(
            ground=_required_string(species_raw, "ground", "species"),
            excited=_required_string(species_raw, "excited", "species"),
            ion=_required_string(species_raw, "ion", "species"),
        ),
        reactions=GecIcpReactionTags(
            elastic=_required_string(reactions_raw, "elastic", "reactions"),
            excitation=_required_string(reactions_raw, "excitation", "reactions"),
            superelastic=_required_string(reactions_raw, "superelastic", "reactions"),
            ionization=_required_string(reactions_raw, "ionization", "reactions"),
            stepwise_ionization=_required_string(
                reactions_raw, "stepwise_ionization", "reactions"
            ),
        ),
    )
    _require_unique_tags(mapping)
    if validate_files and not input_mph.is_file():
        raise GecIcpContractError(f"model.input_mph does not exist: {input_mph}")
    return mapping


def _read_yaml_mapping(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise GecIcpContractError(f"GEC ICP mapping does not exist: {path}")
    try:
        raw = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except (OSError, UnicodeError, yaml.YAMLError) as exc:
        raise GecIcpContractError(f"cannot read GEC ICP mapping: {path}") from exc
    if not isinstance(raw, dict):
        raise GecIcpContractError("GEC ICP mapping root must be a mapping")
    return raw


def _required_mapping(
    raw: dict[str, Any],
    name: str,
    owner: str,
) -> dict[str, Any]:
    value = raw.get(name)
    if not isinstance(value, dict):
        raise GecIcpContractError(f"{owner}.{name} must be a mapping")
    return value


def _required_string(raw: dict[str, Any], name: str, owner: str) -> str:
    value = raw.get(name)
    if not isinstance(value, str) or not value.strip():
        raise GecIcpContractError(f"{owner}.{name} must be a nonempty string")
    return value.strip()


def _required_integer(raw: dict[str, Any], name: str, owner: str) -> int:
    value = raw.get(name)
    if isinstance(value, bool) or not isinstance(value, int):
        raise GecIcpContractError(f"{owner}.{name} must be an integer")
    return value


def _required_boolean(raw: dict[str, Any], name: str, owner: str) -> bool:
    value = raw.get(name)
    if not isinstance(value, bool):
        raise GecIcpContractError(f"{owner}.{name} must be a boolean")
    return value


def _resolved_mph_path(raw: dict[str, Any], name: str, root: Path) -> Path:
    value = _required_string(raw, name, "model")
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = root / path
    path = path.resolve()
    if path.suffix.lower() != ".mph":
        raise GecIcpContractError(f"model.{name} must name an .mph file")
    return path


def _reject_unknown(raw: dict[str, Any], allowed: set[str], owner: str) -> None:
    unknown = sorted(set(raw) - allowed)
    if unknown:
        raise GecIcpContractError(
            f"unsupported {owner} field(s): " + ", ".join(unknown)
        )


def _require_unique_tags(mapping: GecIcpMapping) -> None:
    groups = {
        "species": (
            mapping.species.ground,
            mapping.species.excited,
            mapping.species.ion,
        ),
        "reactions": (
            mapping.reactions.elastic,
            mapping.reactions.excitation,
            mapping.reactions.superelastic,
            mapping.reactions.ionization,
            mapping.reactions.stepwise_ionization,
        ),
        "parameters": (
            mapping.parameters.power,
            mapping.parameters.gas_temperature,
            mapping.parameters.pressure,
        ),
    }
    for owner, tags in groups.items():
        if len(tags) != len(set(tags)):
            raise GecIcpContractError(f"{owner} tags must be unique")
