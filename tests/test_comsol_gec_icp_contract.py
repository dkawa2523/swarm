from __future__ import annotations

import ast
import csv
from dataclasses import replace
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import pytest
import yaml

from swarm_workflow.comsol.models.gec_icp.contracts import (
    GecIcpContractError,
    GecIcpMapping,
    GecIcpMphContract,
)
from swarm_workflow.comsol.models.gec_icp.mapping import load_gec_icp_mapping
from swarm_workflow.comsol.models.gec_icp.mph import (
    inspect_gec_icp_mph,
    validate_gec_icp_mph,
)


ROOT = Path(__file__).resolve().parents[1]
REAL_ICP_MPH = ROOT / "comsol_modes" / "argon_gec_icp.mph"
CANONICAL_ARGON_CROSS_SECTIONS = (
    ROOT / "examples" / "cross_sections" / "argon_application_library.csv"
)


def _canonical_arrays(process: str) -> tuple[tuple[float, ...], tuple[float, ...]]:
    with CANONICAL_ARGON_CROSS_SECTIONS.open(
        "r", encoding="utf-8-sig", newline=""
    ) as stream:
        rows = [row for row in csv.DictReader(stream) if row["process"] == process]
    return (
        tuple(float(row["energy_eV"]) for row in rows),
        tuple(float(row["cross_section_m2"]) for row in rows),
    )


def _mapping_payload(input_mph: str = "model.mph") -> dict[str, object]:
    return {
        "schema_version": 2,
        "model": {
            "input_mph": input_mph,
            "output_mph": "work/result.mph",
            "component": "comp1",
            "geometry": "geom1",
            "mesh": "mesh1",
            "geometry_dimension": 2,
            "axisymmetric": True,
            "plasma_physics": "plas",
            "plasma_feature": "pes1",
            "magnetic_physics": "mf",
            "coil_feature": "coil1",
            "plasma_conductivity_coupling": "pcc1",
            "electron_heat_source_coupling": "ehs1",
            "study": "std1",
            "study_feature": "ftrans",
            "solution": "sol1",
            "dataset": "dset1",
        },
        "parameters": {
            "power": "Psp",
            "gas_temperature": "T0",
            "pressure": "p0",
        },
        "species": {"ground": "Ar", "excited": "Ars", "ion": "Ar_1p"},
        "reactions": {
            "elastic": "eir1",
            "excitation": "eir2",
            "superelastic": "eir3",
            "ionization": "eir4",
            "stepwise_ionization": "eir5",
        },
    }


def _write_mapping(tmp_path: Path, payload: dict[str, object]) -> Path:
    path = tmp_path / "maps" / "gec_icp.yaml"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    return path


def _feature(
    tag: str,
    operation: str,
    parameters: dict[str, str | tuple[float, ...]],
) -> str:
    def encoded(value: str | tuple[float, ...]) -> str:
        if isinstance(value, tuple):
            entries = "|".join(f"1,&apos;{item:.17g}&apos;" for item in value)
            return f"{len(value)}|{entries}"
        return f"1|1,&apos;{value}&apos;"

    values = "".join(
        f'<param param="{name}" value="{encoded(value)}" />'
        for name, value in parameters.items()
    )
    return f'<PhysicsFeature op="{operation}" tag="{tag}">{values}</PhysicsFeature>'


def _dmodel_xml(*, plasma_operation: str = "ColdPlasma") -> str:
    elastic_x, elastic_y = _canonical_arrays("e+Ar=>e+Ar")
    excitation_x, excitation_y = _canonical_arrays("e+Ar=>e+Ars")
    ionization_x, ionization_y = _canonical_arrays("e+Ar=>2e+Ar+")
    species = "".join(
        (
            _feature("e", "Species", {"sType": "electron"}),
            _feature("Ar", "Species", {"sType": "neutral"}),
            _feature("Ars", "Species", {"sType": "neutral"}),
            _feature("Ar_1p", "Species", {"sType": "ion"}),
        )
    )
    reactions = "".join(
        (
            _feature(
                "eir1",
                "ElectronImpactReaction",
                {
                    "formula": "e+Ar=&gt;e+Ar",
                    "type": "Elastic",
                    "de": "0",
                    "SpecifyReactionUsing": "UseCrossSectionData",
                    "xdata": elastic_x,
                    "ydata": elastic_y,
                },
            ),
            _feature(
                "eir2",
                "ElectronImpactReaction",
                {
                    "formula": "e+Ar=&gt;e+Ars",
                    "type": "Excitation",
                    "de": "11.50",
                    "SpecifyReactionUsing": "UseCrossSectionData",
                    "xdata": excitation_x,
                    "ydata": excitation_y,
                },
            ),
            _feature(
                "eir3",
                "ElectronImpactReaction",
                {
                    "formula": "e+Ars=&gt;e+Ar",
                    "type": "Excitation",
                    "de": "-11.50",
                    "SpecifyReactionUsing": "UseCrossSectionData",
                    "xdata": (-11.5, 0.0, 1.2),
                    "ydata": (0.0, 0.0, 6.2e-22),
                },
            ),
            _feature(
                "eir4",
                "ElectronImpactReaction",
                {
                    "formula": "e+Ar=&gt;2e+Ar+",
                    "type": "Ionization",
                    "de": "15.80",
                    "SpecifyReactionUsing": "UseCrossSectionData",
                    "xdata": ionization_x,
                    "ydata": ionization_y,
                },
            ),
            _feature(
                "eir5",
                "ElectronImpactReaction",
                {
                    "formula": "e+Ars=&gt;2e+Ar+",
                    "type": "Ionization",
                    "de": "4.427",
                    "SpecifyReactionUsing": "UseCrossSectionData",
                    "xdata": (0.0, 4.427, 4.628),
                    "ydata": (0.0, 0.0, 1.849e-20),
                },
            ),
        )
    )
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<Model>
  <modelTitle>GEC ICP Reactor, Argon Chemistry</modelTitle>
  <ModelParam tag="param">
    <expressions name="Psp" expr="1500[W]" />
    <expressions name="T0" expr="300[K]" />
    <expressions name="p0" expr="0.02[torr]" />
  </ModelParam>
  <ModelNode tag="comp1"><axisymmetric>true</axisymmetric></ModelNode>
  <GeomSequence tag="geom1" />
  <MeshSequence tag="mesh1" />
  <Physics op="{plasma_operation}" tag="plas">
    <geom>/geom/geom1</geom>
    <PhysicsFeatureList>
      <PhysicsFeature op="PlasmaEsModel" tag="pes1" />
      {species}
      {reactions}
    </PhysicsFeatureList>
  </Physics>
  <Physics op="InductionCurrents" tag="mf">
    <geom>/geom/geom1</geom>
    <PhysicsFeatureList><PhysicsFeature op="Coil" tag="coil1" /></PhysicsFeatureList>
  </Physics>
  <MultiphysicsCoupling op="PlasmaConductivityMultiphysicsCoupling" tag="pcc1">
    <storedParent>/modelNode/comp1</storedParent>
  </MultiphysicsCoupling>
  <MultiphysicsCoupling op="ElectronHeatSourceMultiphysicsCoupling" tag="ehs1">
    <storedParent>/modelNode/comp1</storedParent>
  </MultiphysicsCoupling>
  <Study tag="std1">
    <StudyFeatureList>
      <StudyFeature op="FrequencyTransient" tag="ftrans" />
    </StudyFeatureList>
  </Study>
  <SolverSequence tag="sol1"><study>/study/std1</study></SolverSequence>
  <DatasetFeature op="Solution" tag="dset1">
    <propertyValue name="p:solution" Reference="/sol/sol1" />
  </DatasetFeature>
</Model>
"""


def _write_mph(
    path: Path,
    *,
    dmodel: str | None = None,
    include_dmodel: bool = True,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with ZipFile(path, "w", compression=ZIP_DEFLATED) as archive:
        if include_dmodel:
            archive.writestr("dmodel.xml", dmodel or _dmodel_xml())
        archive.writestr(
            "modelinfo.xml",
            '<modelInfo comsolVersion="6.4.0.257" title="GEC ICP">'
            '<geometryInfo><geom tag="geom1" dimension="2" /></geometryInfo>'
            "</modelInfo>",
        )
        archive.writestr("fileversion", "2091:COMSOL 6.4.0.257")


def test_mapping_loads_strict_schema_v2_with_mapping_directory_paths(
    tmp_path: Path,
) -> None:
    mapping_path = _write_mapping(tmp_path, _mapping_payload())
    input_mph = mapping_path.parent / "model.mph"
    input_mph.touch()

    mapping = load_gec_icp_mapping(mapping_path, validate_files=True)

    assert type(mapping) is GecIcpMapping
    assert mapping.root == mapping_path.parent.resolve()
    assert mapping.model.input_mph == input_mph.resolve()
    assert (
        mapping.model.output_mph == (mapping_path.parent / "work/result.mph").resolve()
    )
    assert mapping.model.component == "comp1"
    assert mapping.parameters.power == "Psp"
    assert mapping.species.excited == "Ars"
    assert mapping.reactions.stepwise_ionization == "eir5"


@pytest.mark.parametrize(
    ("section", "field", "value", "message"),
    (
        (None, "schema_version", 1, "requires schema_version: 2"),
        (None, "schema_version", True, "requires schema_version: 2"),
        ("model", "geometry_dimension", 3, "geometry_dimension must be 2"),
        ("model", "axisymmetric", False, "axisymmetric must be true"),
        ("model", "component", "", "component must be a nonempty string"),
    ),
)
def test_mapping_rejects_invalid_contract_values(
    tmp_path: Path,
    section: str | None,
    field: str,
    value: object,
    message: str,
) -> None:
    payload = _mapping_payload()
    owner = payload if section is None else payload[section]
    assert isinstance(owner, dict)
    owner[field] = value
    path = _write_mapping(tmp_path, payload)

    with pytest.raises(GecIcpContractError, match=message):
        load_gec_icp_mapping(path)


def test_mapping_rejects_unknown_fields_and_duplicate_role_tags(
    tmp_path: Path,
) -> None:
    payload = _mapping_payload()
    model = payload["model"]
    assert isinstance(model, dict)
    model["legacy_study"] = "std0"
    path = _write_mapping(tmp_path, payload)
    with pytest.raises(GecIcpContractError, match="unsupported model field"):
        load_gec_icp_mapping(path)

    payload = _mapping_payload()
    reactions = payload["reactions"]
    assert isinstance(reactions, dict)
    reactions["stepwise_ionization"] = reactions["ionization"]
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    with pytest.raises(GecIcpContractError, match="reactions tags must be unique"):
        load_gec_icp_mapping(path)


def test_inspector_and_validator_cover_the_icp_topology(tmp_path: Path) -> None:
    mph_path = tmp_path / "maps" / "model.mph"
    _write_mph(mph_path)
    mapping_path = _write_mapping(tmp_path, _mapping_payload())
    mapping = load_gec_icp_mapping(mapping_path, validate_files=True)

    contract = inspect_gec_icp_mph(mph_path)
    validated = validate_gec_icp_mph(mapping, contract)

    assert type(validated) is GecIcpMphContract
    assert contract.comsol_version == "6.4.0.257"
    assert contract.model_title == "GEC ICP"
    assert contract.components[0].tag == "comp1"
    assert contract.components[0].axisymmetric is True
    assert contract.geometries[0].dimension == 2
    assert {item.tag: item.operation for item in contract.physics} == {
        "plas": "ColdPlasma",
        "mf": "InductionCurrents",
    }
    assert {item.tag for item in contract.species} == {"e", "Ar", "Ars", "Ar_1p"}
    assert [item.tag for item in contract.reactions] == [
        "eir1",
        "eir2",
        "eir3",
        "eir4",
        "eir5",
    ]
    ground_lengths = {
        item.tag: len(item.energy_data_eV)
        for item in contract.reactions
        if item.tag in {"eir1", "eir2", "eir4"}
    }
    assert ground_lengths == {"eir1": 67, "eir2": 34, "eir4": 31}
    assert contract.datasets[0].solution == "sol1"


def test_validator_rejects_ground_cross_section_point_mismatch(
    tmp_path: Path,
) -> None:
    mph_path = tmp_path / "maps" / "model.mph"
    _write_mph(mph_path)
    mapping = load_gec_icp_mapping(_write_mapping(tmp_path, _mapping_payload()))
    contract = inspect_gec_icp_mph(mph_path)
    reactions = tuple(
        replace(
            item,
            cross_section_data_m2=(
                *item.cross_section_data_m2[:5],
                item.cross_section_data_m2[5] * 1.01,
                *item.cross_section_data_m2[6:],
            ),
        )
        if item.tag == "eir2"
        else item
        for item in contract.reactions
    )

    with pytest.raises(
        GecIcpContractError,
        match="reaction eir2 ydata does not match canonical Argon excitation at index 5",
    ):
        validate_gec_icp_mph(mapping, replace(contract, reactions=reactions))


def test_validator_requires_positive_comsol_owned_cross_section_arrays(
    tmp_path: Path,
) -> None:
    mph_path = tmp_path / "maps" / "model.mph"
    _write_mph(mph_path)
    mapping = load_gec_icp_mapping(_write_mapping(tmp_path, _mapping_payload()))
    contract = inspect_gec_icp_mph(mph_path)
    reactions = tuple(
        replace(
            item,
            cross_section_data_m2=tuple(0.0 for _ in item.cross_section_data_m2),
        )
        if item.tag == "eir3"
        else item
        for item in contract.reactions
    )

    with pytest.raises(
        GecIcpContractError,
        match="reaction eir3 ydata must contain a positive value",
    ):
        validate_gec_icp_mph(mapping, replace(contract, reactions=reactions))


def test_validator_rejects_wrong_physics_operation(tmp_path: Path) -> None:
    mph_path = tmp_path / "maps" / "model.mph"
    _write_mph(mph_path)
    mapping = load_gec_icp_mapping(_write_mapping(tmp_path, _mapping_payload()))
    contract = inspect_gec_icp_mph(mph_path)
    bad_physics = tuple(
        replace(item, operation="ColdPlasmaTimePeriodic")
        if item.tag == "plas"
        else item
        for item in contract.physics
    )

    with pytest.raises(
        GecIcpContractError,
        match="physics plas operation='ColdPlasmaTimePeriodic'",
    ):
        validate_gec_icp_mph(mapping, replace(contract, physics=bad_physics))


def test_inspector_fails_closed_for_nonzip_and_missing_dmodel(tmp_path: Path) -> None:
    nonzip = tmp_path / "not_zip.mph"
    nonzip.write_text("not a ZIP", encoding="utf-8")
    with pytest.raises(GecIcpContractError, match="not a readable ZIP archive"):
        inspect_gec_icp_mph(nonzip)

    missing = tmp_path / "missing.mph"
    _write_mph(missing, include_dmodel=False)
    with pytest.raises(GecIcpContractError, match="missing dmodel.xml"):
        inspect_gec_icp_mph(missing)


def test_icp_contract_modules_do_not_import_gec_ccp() -> None:
    package = ROOT / "swarm_workflow" / "comsol" / "models" / "gec_icp"
    for path in package.glob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        imports = {
            node.module
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom) and node.module is not None
        }
        imports.update(
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.Import)
            for alias in node.names
        )
        assert not any("gec_ccp" in name for name in imports)


@pytest.mark.skipif(not REAL_ICP_MPH.is_file(), reason="local ignored ICP MPH absent")
def test_local_repository_icp_mph_matches_researched_contract(tmp_path: Path) -> None:
    payload = _mapping_payload(str(REAL_ICP_MPH))
    mapping = load_gec_icp_mapping(_write_mapping(tmp_path, payload))

    contract = validate_gec_icp_mph(mapping)

    assert contract.comsol_version == "6.4.0.257"
    assert contract.model_title == "GEC ICP Reactor, Argon Chemistry"
    assert {item.tag for item in contract.reactions} == {
        "eir1",
        "eir2",
        "eir3",
        "eir4",
        "eir5",
    }
