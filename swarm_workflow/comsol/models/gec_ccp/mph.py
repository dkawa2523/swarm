"""Inspect and validate the COMSOL MPH contract for GEC-CCP runs."""

from __future__ import annotations

import html
import math
from pathlib import Path
import re
from zipfile import ZipFile

from swarm_workflow.comsol.models.gec_ccp.contracts import (
    GecCcpMapping,
    GecCcpWorkflowError,
    MphContract,
    MphHeavySpeciesContract,
    MphSurfaceReactionContract,
)


GEC_ARGON_MOLAR_MASS_KG_MOL = 0.04


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
    heavy_species: list[tuple[str, float]] = []
    heavy_species_contracts: list[MphHeavySpeciesContract] = []
    for species_match in re.finditer(
        r'<PhysicsFeature\b[^>]*op="Species"[^>]*tag="([^"]+)"[^>]*>'
        r'(.*?)</PhysicsFeature>',
        xml,
        flags=re.DOTALL,
    ):
        tag, block = species_match.groups()
        species_type = _property_values(block, "sType")[-1]
        if species_type == "electron":
            continue
        molar_mass_expression = _property_values(block, "M")[-1]
        molar_mass_match = re.fullmatch(
            r"([+\-0-9.eE]+)\[kg/mol\]", molar_mass_expression
        )
        if molar_mass_match is None:
            raise GecCcpWorkflowError(
                f"heavy species {tag} has a noncanonical molar mass: "
                f"{molar_mass_expression}"
            )
        molar_mass = float(molar_mass_match.group(1))
        try:
            thermal_diffusion = float(_property_values(block, "DT")[-1])
            charge_number = float(_property_values(block, "z")[-1])
        except ValueError as exc:
            raise GecCcpWorkflowError(
                f"heavy species {tag} has nonnumeric DT or z"
            ) from exc
        mass_constraint = _property_values(block, "FromMassConstraint")[-1]
        if mass_constraint not in {"0", "1"}:
            raise GecCcpWorkflowError(
                f"heavy species {tag} has a nonbinary mass constraint flag: "
                f"{mass_constraint}"
            )
        status_region = block.split("<param", maxsplit=1)[0]
        enabled = not bool(
            re.search(
                r'<entityFlags\b[^>]*>\s*DISABLED\s*</entityFlags>',
                status_region,
            )
        )
        heavy_species.append((tag, molar_mass))
        heavy_species_contracts.append(
            MphHeavySpeciesContract(
                tag=tag,
                species_type=species_type,
                enabled=enabled,
                molar_mass_kg_mol=molar_mass,
                from_mass_constraint=mass_constraint == "1",
                thermal_diffusion=thermal_diffusion,
                charge_number=charge_number,
            )
        )
    surface_reactions: list[MphSurfaceReactionContract] = []
    for surface_match in re.finditer(
        r'<PhysicsFeature\b[^>]*op="SurfaceReaction"[^>]*tag="([^"]+)"'
        r'[^>]*>(.*?)</PhysicsFeature>',
        xml,
        flags=re.DOTALL,
    ):
        tag, block = surface_match.groups()
        status_region = block.split("<param", maxsplit=1)[0]
        surface_reactions.append(
            MphSurfaceReactionContract(
                tag=tag,
                enabled=not bool(
                    re.search(
                        r'<entityFlags\b[^>]*>\s*DISABLED\s*</entityFlags>',
                        status_region,
                    )
                ),
                formula=html.unescape(_property_values(block, "formula")[-1]),
            )
        )
    common_heavy_mass: float | None = None
    if heavy_species:
        candidate = heavy_species[0][1]
        if all(
            math.isclose(value, candidate, rel_tol=0.0, abs_tol=1.0e-15)
            for _, value in heavy_species
        ):
            common_heavy_mass = candidate
    version = re.search(r"COMSOL\s+([0-9.]+)", file_version)
    axisymmetric_values = re.findall(
        r'<axisymmetric\b[^>]*>\s*(true|false)\s*</axisymmetric>',
        xml,
    )
    if len(axisymmetric_values) != 1:
        raise GecCcpWorkflowError(
            "model must contain exactly one finalized axisymmetric setting"
        )
    return MphContract(
        comsol_version=version.group(1) if version else None,
        physics_operation=physics.group(1),
        physics_tag=physics.group(2),
        original_eedf=eedf,
        plasma_feature=plasma,
        reaction_features=reactions,
        studies=studies,
        datasets=datasets,
        heavy_species_molar_masses_kg_mol=tuple(heavy_species),
        common_heavy_species_molar_mass_kg_mol=common_heavy_mass,
        heavy_species=tuple(heavy_species_contracts),
        heavy_species_selection=_property_values(
            xml, "HeavySpeciesSelection"
        )[-1],
        axisymmetric=axisymmetric_values[0] == "true",
        heavy_species_formulation=_single_property_value(
            xml, "Formulation"
        ),
        heavy_species_diffusion_model=_single_property_value(
            xml, "DiffusionModel"
        ),
        heavy_species_migration=_binary_property_value(xml, "Migration"),
        heavy_species_convection=_binary_property_value(xml, "Convection"),
        mixture_diffusion_correction=_binary_property_value(
            xml, "MixtureDiffusionCorrection"
        ),
        ion_tensor_properties=_binary_property_value(xml, "IonTensorProps"),
        ion_electric_field_time_model=_single_property_value(
            xml, "ElectricFieldAppliedToIons"
        ),
        surface_reactions=tuple(surface_reactions),
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
    for study in (
        expected.study,
        expected.external_study,
        expected.conversion_study,
    ):
        if study not in contract.studies:
            failures.append(f"missing study={study}")
    datasets = [
        mapping.run.axis_dataset,
        mapping.run.radial_dataset,
        mapping.run.external_period_dataset,
        mapping.run.phase_dataset,
        mapping.run.external_waveform_dataset,
    ]
    if mapping.run.period_dataset is not None:
        datasets.append(mapping.run.period_dataset)
    if mapping.run.baseline_waveform_dataset is not None:
        datasets.append(mapping.run.baseline_waveform_dataset)
    for dataset in datasets:
        if dataset not in contract.datasets:
            failures.append(f"missing dataset={dataset}")
    if not contract.heavy_species_molar_masses_kg_mol:
        failures.append("heavy-species molar masses are absent")
    elif not math.isclose(
        contract.common_heavy_species_molar_mass_kg_mol or math.nan,
        GEC_ARGON_MOLAR_MASS_KG_MOL,
        rel_tol=0.0,
        abs_tol=1.0e-15,
    ):
        failures.append(
            "heavy-species molar masses are not uniformly "
            f"{GEC_ARGON_MOLAR_MASS_KG_MOL} kg/mol"
        )
    active_heavy = {
        species.tag: species
        for species in contract.heavy_species
        if species.enabled
    }
    if set(active_heavy) != {"Ar", "Ar_1p"}:
        failures.append(
            "active heavy species must be exactly binary Ar and Ar_1p; "
            f"found={sorted(active_heavy)}"
        )
    else:
        neutral = active_heavy["Ar"]
        ion = active_heavy["Ar_1p"]
        if neutral.species_type != "neutral":
            failures.append("Ar must be a neutral species")
        if ion.species_type != "ion":
            failures.append("Ar_1p must be an ion species")
        if not neutral.from_mass_constraint:
            failures.append("Ar must provide the binary mass constraint")
        if ion.from_mass_constraint:
            failures.append("Ar_1p cannot provide the binary mass constraint")
        if not math.isclose(
            neutral.thermal_diffusion, 0.0, rel_tol=0.0, abs_tol=1.0e-15
        ) or not math.isclose(
            ion.thermal_diffusion, 0.0, rel_tol=0.0, abs_tol=1.0e-15
        ):
            failures.append("Ar and Ar_1p thermal diffusion must both be zero")
        if not math.isclose(
            neutral.charge_number, 0.0, rel_tol=0.0, abs_tol=1.0e-15
        ) or not math.isclose(
            ion.charge_number, 1.0, rel_tol=0.0, abs_tol=1.0e-15
        ):
            failures.append("Ar/Ar_1p charge numbers must be 0/1")
    if contract.heavy_species_selection != "BaseGeometry":
        failures.append(
            "heavy species must be solved on BaseGeometry; "
            f"found={contract.heavy_species_selection}"
        )
    expected_transport_contract = {
        "axisymmetric": (contract.axisymmetric, True),
        "heavy-species formulation": (
            contract.heavy_species_formulation,
            "FEMLogLinear",
        ),
        "heavy-species diffusion model": (
            contract.heavy_species_diffusion_model,
            "MixtureAveraged",
        ),
        "heavy-species migration": (contract.heavy_species_migration, True),
        "heavy-species convection": (
            contract.heavy_species_convection,
            False,
        ),
        "mixture diffusion correction": (
            contract.mixture_diffusion_correction,
            False,
        ),
        "ion tensor properties": (contract.ion_tensor_properties, False),
        "ion electric-field time model": (
            contract.ion_electric_field_time_model,
            "Instantaneous",
        ),
    }
    for label, (actual, required) in expected_transport_contract.items():
        if actual != required:
            failures.append(f"{label}={actual}, required={required}")
    active_surface_reactions = {
        reaction.tag: reaction.formula
        for reaction in contract.surface_reactions
        if reaction.enabled
    }
    if active_surface_reactions != {"sr1": "Ar+=>Ar"}:
        failures.append(
            "active surface reaction must be exactly sr1: Ar+=>Ar; "
            f"found={active_surface_reactions}"
        )
    if failures:
        raise GecCcpWorkflowError(
            "local MPH does not match the GEC CCP mapping: " + "; ".join(failures)
        )
    output_paths = {expected.output_mph}
    if expected.baseline_output_mph is not None:
        output_paths.add(expected.baseline_output_mph)
    if expected.input_mph in output_paths:
        raise GecCcpWorkflowError("output MPH paths must not overwrite the input MPH")
    if (
        expected.baseline_output_mph is not None
        and expected.baseline_output_mph == expected.output_mph
    ):
        raise GecCcpWorkflowError(
            "built-in reference and external output MPH paths must differ"
        )


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


def _single_property_value(xml: str, name: str) -> str:
    values = _property_values(xml, name)
    if len(values) != 1:
        raise GecCcpWorkflowError(
            f"model property must occur exactly once: {name}; "
            f"found={len(values)}"
        )
    return values[0]


def _binary_property_value(xml: str, name: str) -> bool:
    value = _single_property_value(xml, name)
    if value not in {"0", "1"}:
        raise GecCcpWorkflowError(
            f"model property must be binary: {name}={value}"
        )
    return value == "1"


def _tag_for_operation(xml: str, element: str, operation: str) -> str:
    match = re.search(
        rf'<{element}\b[^>]*op="{re.escape(operation)}"[^>]*tag="([^"]+)"',
        xml,
    )
    if match is None:
        raise GecCcpWorkflowError(f"model operation not found: {operation}")
    return match.group(1)
