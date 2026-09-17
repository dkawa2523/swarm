"""Safely inspect and validate the offline GEC ICP MPH contract."""

from __future__ import annotations

import math
from pathlib import Path
import re
from xml.etree import ElementTree
from zipfile import BadZipFile, ZipFile

from .contracts import (
    GecIcpContractError,
    GecIcpMapping,
    GecIcpMphContract,
    MphComponentContract,
    MphCouplingContract,
    MphDatasetContract,
    MphGeometryContract,
    MphParameterContract,
    MphPhysicsContract,
    MphPhysicsFeatureContract,
    MphReactionContract,
    MphSolutionContract,
    MphSpeciesContract,
    MphStudyFeatureContract,
)
from .cross_sections import load_canonical_ground_cross_sections


_MAX_DMODEL_BYTES = 16 * 1024 * 1024
_MAX_MODELINFO_BYTES = 1024 * 1024
_MAX_FILEVERSION_BYTES = 4096

_REACTION_SIGNATURES = {
    "elastic": ("e+Ar=>e+Ar", "Elastic", 0.0),
    "excitation": ("e+Ar=>e+Ars", "Excitation", 11.50),
    "superelastic": ("e+Ars=>e+Ar", "Excitation", -11.50),
    "ionization": ("e+Ar=>2e+Ar+", "Ionization", 15.80),
    "stepwise_ionization": ("e+Ars=>2e+Ar+", "Ionization", 4.427),
}


def inspect_gec_icp_mph(path: str | Path) -> GecIcpMphContract:
    """Read only bounded XML members from a ZIP-format MPH."""

    mph_path = Path(path).resolve()
    if not mph_path.is_file():
        raise GecIcpContractError(f"GEC ICP MPH does not exist: {mph_path}")
    try:
        with ZipFile(mph_path) as archive:
            dmodel = _read_bounded_member(
                archive, "dmodel.xml", _MAX_DMODEL_BYTES, required=True
            )
            modelinfo = _read_bounded_member(
                archive, "modelinfo.xml", _MAX_MODELINFO_BYTES, required=False
            )
            fileversion = _read_bounded_member(
                archive, "fileversion", _MAX_FILEVERSION_BYTES, required=False
            )
    except (BadZipFile, OSError) as exc:
        raise GecIcpContractError(
            f"GEC ICP MPH is not a readable ZIP archive: {mph_path}"
        ) from exc

    root = _parse_xml(dmodel, "dmodel.xml")
    info_root = _parse_xml(modelinfo, "modelinfo.xml") if modelinfo else None
    version_text = fileversion.decode("utf-8", errors="replace")
    version_match = re.search(r"COMSOL\s+([0-9.]+)", version_text)
    comsol_version = version_match.group(1) if version_match else None
    if comsol_version is None and info_root is not None:
        comsol_version = info_root.attrib.get("comsolVersion")

    components = tuple(
        MphComponentContract(
            tag=_required_attribute(element, "tag", "ModelNode"),
            axisymmetric=_component_axisymmetric(root, element),
        )
        for element in _elements(root, "ModelNode")
    )
    geometries = _geometry_contracts(info_root)
    meshes = tuple(
        _required_attribute(element, "tag", "MeshSequence")
        for element in _elements(root, "MeshSequence")
    )
    physics = tuple(
        _physics_contract(element) for element in _elements(root, "Physics")
    )
    physics_features = tuple(
        feature
        for physics_element in _elements(root, "Physics")
        for feature in _physics_features(physics_element)
    )
    couplings = tuple(
        _coupling_contract(element)
        for element in _elements(root, "MultiphysicsCoupling")
    )
    studies = tuple(
        _required_attribute(element, "tag", "Study")
        for element in _elements(root, "Study")
    )
    study_features = tuple(
        feature
        for study in _elements(root, "Study")
        for feature in _study_features(study)
    )
    solutions = tuple(
        MphSolutionContract(
            tag=_required_attribute(element, "tag", "SolverSequence"),
            study=_reference_child_tag(element, "study"),
        )
        for element in _elements(root, "SolverSequence")
    )
    datasets = tuple(
        _dataset_contract(element) for element in _elements(root, "DatasetFeature")
    )
    parameters = tuple(
        MphParameterContract(
            name=_required_attribute(element, "name", "parameter expression"),
            expression=_required_attribute(element, "expr", "parameter expression"),
        )
        for model_parameters in _elements(root, "ModelParam")
        for element in model_parameters
        if _local_name(element.tag) == "expressions"
    )
    species = tuple(
        _species_contract(feature) for feature in _feature_elements(root, "Species")
    )
    reactions = tuple(
        _reaction_contract(feature)
        for feature in _feature_elements(root, "ElectronImpactReaction")
    )
    model_title = None
    if info_root is not None:
        model_title = info_root.attrib.get("title")
    if model_title is None:
        title = next(iter(_elements(root, "modelTitle")), None)
        model_title = title.text.strip() if title is not None and title.text else None

    return GecIcpMphContract(
        path=mph_path,
        comsol_version=comsol_version,
        model_title=model_title,
        components=components,
        geometries=geometries,
        meshes=meshes,
        physics=physics,
        physics_features=physics_features,
        couplings=couplings,
        studies=studies,
        study_features=study_features,
        solutions=solutions,
        datasets=datasets,
        parameters=parameters,
        species=species,
        reactions=reactions,
    )


def validate_gec_icp_mph(
    mapping: GecIcpMapping,
    contract: GecIcpMphContract | None = None,
) -> GecIcpMphContract:
    """Validate the researched GEC ICP topology against a typed mapping."""

    inspected = contract or inspect_gec_icp_mph(mapping.model.input_mph)
    failures: list[str] = []
    model = mapping.model

    component = _single_tagged(
        inspected.components, model.component, "component", failures
    )
    if component is not None and component.axisymmetric != model.axisymmetric:
        failures.append(
            f"component {model.component} axisymmetric={component.axisymmetric}"
        )
    geometry = _single_tagged(
        inspected.geometries, model.geometry, "geometry", failures
    )
    if geometry is not None and geometry.dimension != model.geometry_dimension:
        failures.append(f"geometry {model.geometry} dimension={geometry.dimension}")
    if inspected.meshes.count(model.mesh) != 1:
        failures.append(f"mesh tag {model.mesh!r} must occur exactly once")

    _require_physics(
        inspected,
        model.plasma_physics,
        "ColdPlasma",
        model.geometry,
        failures,
    )
    _require_physics(
        inspected,
        model.magnetic_physics,
        "InductionCurrents",
        model.geometry,
        failures,
    )
    _require_physics_feature(
        inspected,
        model.plasma_physics,
        model.plasma_feature,
        "PlasmaEsModel",
        failures,
    )
    _require_physics_feature(
        inspected,
        model.magnetic_physics,
        model.coil_feature,
        "Coil",
        failures,
    )
    _require_coupling(
        inspected,
        model.plasma_conductivity_coupling,
        "PlasmaConductivityMultiphysicsCoupling",
        model.component,
        failures,
    )
    _require_coupling(
        inspected,
        model.electron_heat_source_coupling,
        "ElectronHeatSourceMultiphysicsCoupling",
        model.component,
        failures,
    )
    if inspected.studies.count(model.study) != 1:
        failures.append(f"study tag {model.study!r} must occur exactly once")
    _require_study_feature(
        inspected,
        model.study,
        model.study_feature,
        "FrequencyTransient",
        failures,
    )
    solution = _single_tagged(inspected.solutions, model.solution, "solution", failures)
    if solution is not None and solution.study != model.study:
        failures.append(f"solution {model.solution} refers to study {solution.study!r}")
    dataset = _single_tagged(inspected.datasets, model.dataset, "dataset", failures)
    if dataset is not None:
        if dataset.operation != "Solution":
            failures.append(f"dataset {model.dataset} operation={dataset.operation!r}")
        if dataset.solution != model.solution:
            failures.append(
                f"dataset {model.dataset} refers to solution {dataset.solution!r}"
            )

    for name in (
        mapping.parameters.power,
        mapping.parameters.gas_temperature,
        mapping.parameters.pressure,
    ):
        count = sum(item.name == name for item in inspected.parameters)
        if count != 1:
            failures.append(f"parameter {name!r} must occur exactly once")

    _validate_species(mapping, inspected, failures)
    _validate_reactions(mapping, inspected, failures)
    if failures:
        raise GecIcpContractError(
            "local MPH does not match the GEC ICP mapping: " + "; ".join(failures)
        )
    return inspected


def _read_bounded_member(
    archive: ZipFile,
    name: str,
    maximum_bytes: int,
    *,
    required: bool,
) -> bytes:
    try:
        info = archive.getinfo(name)
    except KeyError:
        if required:
            raise GecIcpContractError(f"GEC ICP MPH is missing {name}") from None
        return b""
    if info.file_size > maximum_bytes:
        raise GecIcpContractError(
            f"GEC ICP MPH member {name} exceeds {maximum_bytes} bytes"
        )
    data = archive.read(info)
    if len(data) != info.file_size:
        raise GecIcpContractError(f"GEC ICP MPH member {name} is truncated")
    return data


def _parse_xml(data: bytes, name: str) -> ElementTree.Element:
    upper = data.upper()
    if b"<!DOCTYPE" in upper or b"<!ENTITY" in upper:
        raise GecIcpContractError(f"GEC ICP MPH member {name} contains a DTD")
    try:
        return ElementTree.fromstring(data)
    except ElementTree.ParseError as exc:
        raise GecIcpContractError(
            f"GEC ICP MPH member {name} is not valid XML"
        ) from exc


def _local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _elements(root: ElementTree.Element, name: str) -> list[ElementTree.Element]:
    return [element for element in root.iter() if _local_name(element.tag) == name]


def _required_attribute(
    element: ElementTree.Element,
    name: str,
    owner: str,
) -> str:
    value = element.attrib.get(name)
    if not value:
        raise GecIcpContractError(f"{owner} lacks attribute {name}")
    return value


def _direct_child(
    element: ElementTree.Element,
    name: str,
) -> ElementTree.Element | None:
    return next(
        (child for child in element if _local_name(child.tag) == name),
        None,
    )


def _boolean_child(element: ElementTree.Element, name: str) -> bool:
    child = _direct_child(element, name)
    if (
        child is None
        or child.text is None
        or child.text.strip() not in {"true", "false"}
    ):
        raise GecIcpContractError(f"{_local_name(element.tag)} lacks boolean {name}")
    return child.text.strip() == "true"


def _component_axisymmetric(
    root: ElementTree.Element,
    component: ElementTree.Element,
) -> bool:
    direct = _direct_child(component, "axisymmetric")
    if direct is not None:
        return _boolean_child(component, "axisymmetric")

    component_tag = _required_attribute(component, "tag", "ModelNode")
    geometries = [
        geometry
        for geometry in _elements(root, "GeomSequence")
        if _reference_child_tag(geometry, "storedParent") == component_tag
    ]
    if len(geometries) != 1:
        raise GecIcpContractError(
            f"component {component_tag} must own exactly one geometry sequence"
        )
    return _boolean_child(geometries[0], "axisymmetric")


def _reference_child_tag(element: ElementTree.Element, name: str) -> str | None:
    child = _direct_child(element, name)
    if child is None or child.text is None:
        return None
    return _last_path_tag(child.text)


def _last_path_tag(value: str | None) -> str | None:
    if value is None or not value.strip():
        return None
    return value.strip().rstrip("/").rsplit("/", 1)[-1]


def _geometry_contracts(
    info_root: ElementTree.Element | None,
) -> tuple[MphGeometryContract, ...]:
    if info_root is None:
        return ()
    contracts: list[MphGeometryContract] = []
    for element in _elements(info_root, "geom"):
        tag = _required_attribute(element, "tag", "geometry info")
        raw_dimension = _required_attribute(element, "dimension", "geometry info")
        try:
            dimension = int(raw_dimension)
        except ValueError as exc:
            raise GecIcpContractError(
                f"geometry {tag} has noninteger dimension {raw_dimension!r}"
            ) from exc
        contracts.append(MphGeometryContract(tag=tag, dimension=dimension))
    return tuple(contracts)


def _physics_contract(element: ElementTree.Element) -> MphPhysicsContract:
    geometry = _reference_child_tag(element, "geom")
    return MphPhysicsContract(
        tag=_required_attribute(element, "tag", "Physics"),
        operation=_required_attribute(element, "op", "Physics"),
        geometry=geometry,
    )


def _physics_features(
    physics: ElementTree.Element,
) -> tuple[MphPhysicsFeatureContract, ...]:
    physics_tag = _required_attribute(physics, "tag", "Physics")
    return tuple(
        MphPhysicsFeatureContract(
            physics=physics_tag,
            tag=_required_attribute(feature, "tag", "PhysicsFeature"),
            operation=_required_attribute(feature, "op", "PhysicsFeature"),
            enabled=_enabled(feature),
        )
        for feature in physics.iter()
        if _local_name(feature.tag) == "PhysicsFeature"
    )


def _coupling_contract(element: ElementTree.Element) -> MphCouplingContract:
    return MphCouplingContract(
        tag=_required_attribute(element, "tag", "MultiphysicsCoupling"),
        operation=_required_attribute(element, "op", "MultiphysicsCoupling"),
        component=_reference_child_tag(element, "storedParent"),
        enabled=_enabled(element),
    )


def _study_features(study: ElementTree.Element) -> tuple[MphStudyFeatureContract, ...]:
    study_tag = _required_attribute(study, "tag", "Study")
    return tuple(
        MphStudyFeatureContract(
            study=study_tag,
            tag=_required_attribute(feature, "tag", "StudyFeature"),
            operation=_required_attribute(feature, "op", "StudyFeature"),
            enabled=_enabled(feature),
        )
        for feature in study.iter()
        if _local_name(feature.tag) == "StudyFeature"
    )


def _dataset_contract(element: ElementTree.Element) -> MphDatasetContract:
    solution = None
    for child in element:
        if child.attrib.get("name") == "p:solution":
            solution = _last_path_tag(child.attrib.get("Reference"))
            break
    return MphDatasetContract(
        tag=_required_attribute(element, "tag", "DatasetFeature"),
        operation=_required_attribute(element, "op", "DatasetFeature"),
        solution=solution,
    )


def _feature_elements(
    root: ElementTree.Element,
    operation: str,
) -> list[tuple[str, ElementTree.Element]]:
    matches: list[tuple[str, ElementTree.Element]] = []
    for physics in _elements(root, "Physics"):
        physics_tag = _required_attribute(physics, "tag", "Physics")
        for feature in physics.iter():
            if (
                _local_name(feature.tag) == "PhysicsFeature"
                and feature.attrib.get("op") == operation
            ):
                matches.append((physics_tag, feature))
    return matches


def _species_contract(
    item: tuple[str, ElementTree.Element],
) -> MphSpeciesContract:
    physics, feature = item
    return MphSpeciesContract(
        physics=physics,
        tag=_required_attribute(feature, "tag", "Species"),
        species_type=_required_param(feature, "sType"),
        enabled=_enabled(feature),
    )


def _reaction_contract(
    item: tuple[str, ElementTree.Element],
) -> MphReactionContract:
    physics, feature = item
    raw_energy = _required_param(feature, "de")
    try:
        energy_loss = float(raw_energy)
    except ValueError as exc:
        raise GecIcpContractError(
            f"reaction {_required_attribute(feature, 'tag', 'reaction')} "
            f"has nonnumeric energy loss {raw_energy!r}"
        ) from exc
    return MphReactionContract(
        physics=physics,
        tag=_required_attribute(feature, "tag", "ElectronImpactReaction"),
        formula=_required_param(feature, "formula"),
        collision_type=_required_param(feature, "type"),
        energy_loss_eV=energy_loss,
        specification=_required_param(feature, "SpecifyReactionUsing"),
        energy_data_eV=_required_numeric_param_array(feature, "xdata"),
        cross_section_data_m2=_required_numeric_param_array(feature, "ydata"),
        enabled=_enabled(feature),
    )


def _required_param(element: ElementTree.Element, name: str) -> str:
    values = [
        child.attrib["value"]
        for child in element
        if _local_name(child.tag) == "param"
        and child.attrib.get("param") == name
        and "value" in child.attrib
    ]
    if len(values) != 1:
        tag = element.attrib.get("tag", _local_name(element.tag))
        raise GecIcpContractError(
            f"MPH feature {tag} must contain exactly one parameter {name}"
        )
    quoted = re.findall(r"'([^']*)'", values[0])
    return quoted[0] if quoted else values[0]


def _required_numeric_param_array(
    element: ElementTree.Element,
    name: str,
) -> tuple[float, ...]:
    values = [
        child.attrib["value"]
        for child in element
        if _local_name(child.tag) == "param"
        and child.attrib.get("param") == name
        and "value" in child.attrib
    ]
    tag = element.attrib.get("tag", _local_name(element.tag))
    if len(values) != 1:
        raise GecIcpContractError(
            f"MPH feature {tag} must contain exactly one parameter {name}"
        )
    raw = values[0]
    count_match = re.match(r"([0-9]+)\|", raw)
    quoted = re.findall(r"'([^']*)'", raw)
    if count_match is None or int(count_match.group(1)) != len(quoted):
        raise GecIcpContractError(
            f"MPH feature {tag} parameter {name} has inconsistent array encoding"
        )
    try:
        numeric = tuple(float(value) for value in quoted)
    except ValueError as exc:
        raise GecIcpContractError(
            f"MPH feature {tag} parameter {name} contains nonnumeric data"
        ) from exc
    if not numeric or any(not math.isfinite(value) for value in numeric):
        raise GecIcpContractError(
            f"MPH feature {tag} parameter {name} must contain finite numeric data"
        )
    return numeric


def _enabled(element: ElementTree.Element) -> bool:
    flags = _direct_child(element, "entityFlags")
    return flags is None or flags.text is None or "DISABLED" not in flags.text


def _single_tagged(items: tuple, tag: str, owner: str, failures: list[str]):
    matches = [item for item in items if item.tag == tag]
    if len(matches) != 1:
        failures.append(f"{owner} tag {tag!r} must occur exactly once")
        return None
    return matches[0]


def _require_physics(
    contract: GecIcpMphContract,
    tag: str,
    operation: str,
    geometry: str,
    failures: list[str],
) -> None:
    item = _single_tagged(contract.physics, tag, "physics", failures)
    if item is None:
        return
    if item.operation != operation:
        failures.append(f"physics {tag} operation={item.operation!r}")
    if item.geometry != geometry:
        failures.append(f"physics {tag} geometry={item.geometry!r}")


def _require_physics_feature(
    contract: GecIcpMphContract,
    physics: str,
    tag: str,
    operation: str,
    failures: list[str],
) -> None:
    matches = [
        item
        for item in contract.physics_features
        if item.physics == physics and item.tag == tag
    ]
    if len(matches) != 1:
        failures.append(f"physics feature {physics}/{tag} must occur exactly once")
        return
    item = matches[0]
    if item.operation != operation:
        failures.append(f"physics feature {physics}/{tag} operation={item.operation!r}")
    if not item.enabled:
        failures.append(f"physics feature {physics}/{tag} is disabled")


def _require_coupling(
    contract: GecIcpMphContract,
    tag: str,
    operation: str,
    component: str,
    failures: list[str],
) -> None:
    item = _single_tagged(contract.couplings, tag, "coupling", failures)
    if item is None:
        return
    if item.operation != operation:
        failures.append(f"coupling {tag} operation={item.operation!r}")
    if item.component != component:
        failures.append(f"coupling {tag} component={item.component!r}")
    if not item.enabled:
        failures.append(f"coupling {tag} is disabled")


def _require_study_feature(
    contract: GecIcpMphContract,
    study: str,
    tag: str,
    operation: str,
    failures: list[str],
) -> None:
    matches = [
        item
        for item in contract.study_features
        if item.study == study and item.tag == tag
    ]
    if len(matches) != 1:
        failures.append(f"study feature {study}/{tag} must occur exactly once")
        return
    item = matches[0]
    if item.operation != operation:
        failures.append(f"study feature {study}/{tag} operation={item.operation!r}")
    if not item.enabled:
        failures.append(f"study feature {study}/{tag} is disabled")


def _validate_species(
    mapping: GecIcpMapping,
    contract: GecIcpMphContract,
    failures: list[str],
) -> None:
    expected = {
        mapping.species.ground: "neutral",
        mapping.species.excited: "neutral",
        mapping.species.ion: "ion",
    }
    active_heavy = {
        item.tag: item
        for item in contract.species
        if item.enabled and item.species_type != "electron"
    }
    if set(active_heavy) != set(expected):
        failures.append(
            f"active heavy species={sorted(active_heavy)}, expected={sorted(expected)}"
        )
        return
    for tag, species_type in expected.items():
        item = active_heavy[tag]
        if item.physics != mapping.model.plasma_physics:
            failures.append(f"species {tag} physics={item.physics!r}")
        if item.species_type != species_type:
            failures.append(f"species {tag} type={item.species_type!r}")


def _validate_reactions(
    mapping: GecIcpMapping,
    contract: GecIcpMphContract,
    failures: list[str],
) -> None:
    canonical = load_canonical_ground_cross_sections()
    mapped = {name: getattr(mapping.reactions, name) for name in _REACTION_SIGNATURES}
    active = {item.tag: item for item in contract.reactions if item.enabled}
    if set(active) != set(mapped.values()):
        failures.append(
            "active electron-impact reactions="
            f"{sorted(active)}, expected={sorted(mapped.values())}"
        )
        return
    for role, tag in mapped.items():
        item = active[tag]
        formula, collision_type, energy_loss = _REACTION_SIGNATURES[role]
        if item.physics != mapping.model.plasma_physics:
            failures.append(f"reaction {tag} physics={item.physics!r}")
        if item.formula != formula:
            failures.append(f"reaction {tag} formula={item.formula!r}")
        if item.collision_type != collision_type:
            failures.append(f"reaction {tag} type={item.collision_type!r}")
        if not math.isclose(
            item.energy_loss_eV,
            energy_loss,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ):
            failures.append(f"reaction {tag} energy_loss_eV={item.energy_loss_eV}")
        if item.specification != "UseCrossSectionData":
            failures.append(f"reaction {tag} specification={item.specification!r}")
        energies = item.energy_data_eV
        cross_sections = item.cross_section_data_m2
        if len(energies) < 2 or len(cross_sections) < 2:
            failures.append(f"reaction {tag} cross-section arrays must contain data")
            continue
        if len(energies) != len(cross_sections):
            failures.append(
                f"reaction {tag} xdata/ydata lengths differ: "
                f"{len(energies)} != {len(cross_sections)}"
            )
            continue
        if any(right <= left for left, right in zip(energies, energies[1:])):
            failures.append(f"reaction {tag} xdata must be strictly increasing")
        if any(value < 0.0 for value in cross_sections):
            failures.append(f"reaction {tag} ydata must be nonnegative")
        if not any(value > 0.0 for value in cross_sections):
            failures.append(f"reaction {tag} ydata must contain a positive value")

        if role not in canonical:
            continue
        reference = canonical[role]
        if (
            reference.species != "Ar"
            or item.formula != reference.process
            or item.collision_type.lower() != reference.process_type
            or item.energy_loss_eV != reference.threshold_eV
        ):
            failures.append(
                f"reaction {tag} ground-state identity does not match canonical "
                f"Argon {role} species/process/type/threshold"
            )
        if energies != reference.energy_eV:
            mismatch = _first_mismatch(energies, reference.energy_eV)
            failures.append(
                f"reaction {tag} xdata does not match canonical Argon {role} "
                f"at index {mismatch}"
            )
        if cross_sections != reference.cross_section_m2:
            mismatch = _first_mismatch(cross_sections, reference.cross_section_m2)
            failures.append(
                f"reaction {tag} ydata does not match canonical Argon {role} "
                f"at index {mismatch}"
            )


def _first_mismatch(
    actual: tuple[float, ...],
    expected: tuple[float, ...],
) -> int:
    for index, (left, right) in enumerate(zip(actual, expected)):
        if left != right:
            return index
    return min(len(actual), len(expected))
