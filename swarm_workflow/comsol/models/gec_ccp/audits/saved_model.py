"""Audit saved GEC-CCP MPH transport, reaction, and elastic bindings."""

from __future__ import annotations

from dataclasses import dataclass
import html
import math
import re
from typing import Any

import numpy as np

from swarm_workflow.comsol.models.gec_ccp.closure import (
    GEC_RESTRICTED_TRANSPORT_FUNCTIONS,
    GEC_TRANSPORT_PROPERTIES,
    GEC_TWO_TERM_HYBRID_TRANSPORT_FUNCTIONS,
    _active_rate_table,
    _reaction_uses_preintegrated_rate,
    _thermal_diffusion_enabled,
    _transport_audit_components,
    _transport_input_contract,
    _uses_direct_mc_rates,
    _uses_external_elastic_energy_loss,
    _uses_hybrid_einstein_transport,
)
from swarm_workflow.comsol.models.gec_ccp.audits.conservation import (
    _read_comsol_numeric_csv,
)
from swarm_workflow.comsol.models.gec_ccp.contracts import (
    AVOGADRO_PER_MOL,
    GEC_RESTRICTED_TRANSPORT_CLOSURE,
    GEC_TWO_TERM_HYBRID_TRANSPORT_CLOSURE,
    GecCcpMapping,
    GecCcpWorkflowError,
)
from swarm_workflow.comsol.models.gec_ccp.data import _read_csv
from swarm_workflow.comsol.models.gec_ccp.closure_arguments import (
    closure_argument_range,
    smooth_log_energy_argument,
)
from swarm_workflow.comsol.models.gec_ccp.validation.joint_consistency import (
    _active_closure_mean_energy_support,
)
from swarm_workflow.comsol.models.gec_ccp.mph import _property_values
from swarm_workflow.comsol.models.gec_ccp.audits.transport_support import (
    _operating_mean_energy_support_audit,
)
from swarm_workflow.comsol.models.gec_ccp.audits.transport_values import (
    _independent_domain_transport_value_audit,
    _independent_elastic_energy_loss_value_audit,
    _independent_external_rate_value_audit,
    _independent_transport_value_audit,
)


GEC_REACTION_PHYSICS_CONTRACT = {
    "elastic": {"formula": "e+Ar=>e+Ar", "energy_loss_eV": 0.0},
    "excitation": {"formula": "e+Ar=>e+Ars", "energy_loss_eV": 11.5},
    "ionization": {"formula": "e+Ar=>2e+Ar+", "energy_loss_eV": 15.8},
}


def _saved_sw_logeps_audit(
    model_xml: str,
    *,
    minimum_eV: float,
    maximum_eV: float,
    variable_name: str = "sw_logeps",
) -> dict[str, Any]:
    """Require one saved closure variable with the exact generated formula."""

    blocks = []
    for attributes, body in re.findall(
        r"<Expr\b([^>]*)>(.*?)</Expr>", model_xml, flags=re.DOTALL
    ):
        tag = re.search(r'\btag="([^"]*)"', attributes)
        if tag is not None and tag.group(1) == "swClosureVars":
            blocks.append(body)
    saved: list[str] = []
    for body in blocks:
        for attributes in re.findall(r"<expressions\b([^>]*)>", body):
            name = re.search(r'\bname="([^"]*)"', attributes)
            expression = re.search(r'\bexpr="([^"]*)"', attributes)
            if name is not None and name.group(1) == variable_name:
                saved.append(
                    html.unescape(expression.group(1)) if expression is not None else ""
                )
    expected = smooth_log_energy_argument(
        "log(ptp.ebar/1[V])",
        minimum_eV=minimum_eV,
        maximum_eV=maximum_eV,
    )
    passed = len(blocks) == 1 and saved == [expected]
    return {
        "passed": passed,
        "variable_name": variable_name,
        "variable_feature_count": len(blocks),
        "sw_logeps_expression_count": len(saved),
        "expected_expression": expected,
        "saved_expression": saved[0] if len(saved) == 1 else None,
        "native_physical_mean_energy_argument": bool(
            saved
            and "log(ptp.ebar/1[V])" in saved[0]
            and "ptp.en/" not in saved[0]
            and "ptp.ne" not in saved[0]
        ),
    }


@dataclass(frozen=True)
class _SavedPhysicsSections:
    passed: bool
    flags: dict[str, dict[str, Any]]
    periodic_restoration: dict[str, Any]
    heavy_species: dict[str, Any]
    reaction_handling: dict[str, Any]


@dataclass(frozen=True)
class _SavedTransportXmlAudit:
    passed: bool
    closure_support: dict[str, Any]
    mode: dict[str, Any]
    bindings: dict[str, dict[str, Any]]
    saved_log_argument: dict[str, Any]
    physics: _SavedPhysicsSections


@dataclass(frozen=True)
class _SavedTransportNumericAudit:
    passed: bool
    operating_range: dict[str, Any]
    orientation_regularization: dict[str, Any]
    comparisons: dict[str, Any]
    transport_values: dict[str, Any]
    elastic_energy_loss_values: dict[str, Any]
    external_rate_values: dict[str, Any]


def _audit_transport_binding(
    mapping: GecCcpMapping,
    model_xml: str,
) -> dict[str, Any]:
    """Audit saved expressions and independent domain values by responsibility."""

    audit_components = _transport_audit_components(
        mapping.closure,
        source=mapping.bundle.expected_source,
    )
    if not audit_components:
        return {"passed": True, "status": "not_applicable"}
    feature_match = re.search(
        rf'<PhysicsFeature\b[^>]*tag="'
        rf'{re.escape(mapping.model.plasma_feature)}"[^>]*>'
        rf"(.*?)</PhysicsFeature>",
        model_xml,
        flags=re.DOTALL,
    )
    if feature_match is None:
        return {"passed": False, "reason": "plasma feature XML is absent"}
    xml_audit = _audit_saved_transport_xml(
        mapping,
        model_xml,
        feature_match.group(1),
    )
    numeric_audit = _audit_saved_transport_numeric(
        mapping,
        xml_audit.closure_support,
    )
    return {
        "passed": xml_audit.passed and numeric_audit.passed,
        "mode": xml_audit.mode,
        "saved_mph_properties": xml_audit.bindings,
        "electron_properties_flags": xml_audit.physics.flags,
        "saved_sw_logeps": xml_audit.saved_log_argument,
        "periodic_solver_equation_view_restoration": (
            xml_audit.physics.periodic_restoration
        ),
        "native_heavy_species_ownership": xml_audit.physics.heavy_species,
        "reaction_and_elastic_energy_ownership": (xml_audit.physics.reaction_handling),
        "tensor_basis": "r_phi_z",
        "numeric_binding_dataset": "phase_resolved_radial_midplane_cut",
        "zero_field_isotropization_Td": (mapping.closure.zero_field_isotropization_Td),
        "operating_mean_energy_range": numeric_audit.operating_range,
        "orientation_regularization": (numeric_audit.orientation_regularization),
        "phase_local_values": numeric_audit.comparisons,
        "independent_bundle_value_audit": numeric_audit.transport_values,
        "independent_elastic_energy_loss_value_audit": (
            numeric_audit.elastic_energy_loss_values
        ),
        "independent_external_rate_value_audit": (numeric_audit.external_rate_values),
    }


def _audit_saved_transport_xml(
    mapping: GecCcpMapping,
    model_xml: str,
    feature_xml: str,
) -> _SavedTransportXmlAudit:
    expected_mode = (
        "SpecifyMueOnly"
        if mapping.closure.electron_transport == "swarm_mobility_einstein"
        else "SpecifyAll"
    )
    try:
        actual_mode = _property_values(
            feature_xml,
            "SpecifyElectronDensityAndEnergy",
        )[-1]
    except GecCcpWorkflowError:
        actual_mode = None
    expected_tokens, forbidden_tokens = _transport_binding_tokens(mapping)
    closure_support = _active_closure_mean_energy_support(mapping)
    saved_log_argument = _saved_closure_argument_audit(
        mapping,
        model_xml,
        closure_support,
    )
    bindings, bindings_passed = _saved_transport_property_bindings(
        mapping,
        feature_xml,
        expected_tokens,
        forbidden_tokens,
    )
    physics = _saved_physics_sections(mapping, model_xml, feature_xml)
    passed = bool(
        actual_mode == expected_mode
        and saved_log_argument["passed"]
        and bindings_passed
        and physics.passed
    )
    return _SavedTransportXmlAudit(
        passed=passed,
        closure_support=closure_support,
        mode={
            "expected": expected_mode,
            "actual": actual_mode,
            "passed": actual_mode == expected_mode,
        },
        bindings=bindings,
        saved_log_argument=saved_log_argument,
        physics=physics,
    )


def _transport_binding_tokens(
    mapping: GecCcpMapping,
) -> tuple[dict[str, tuple[str, ...]], dict[str, tuple[str, ...]]]:
    expected: dict[str, tuple[str, ...]] = {"muN": ("sw_log_muN_e",)}
    forbidden: dict[str, tuple[str, ...]] = {}
    if _uses_hybrid_einstein_transport(mapping.closure):
        expected.update(
            {
                "DeN": ("sw_log_muN_e", "ptp.Te"),
                "muenN": ("sw_log_muenN_e",),
                "DenN": (
                    ("sw_log_DenN_L_e", "sw_log_DenN_T_e")
                    if mapping.bundle.expected_source == "monte_carlo"
                    else ("sw_log_DenN_e",)
                ),
            }
        )
        forbidden["DeN"] = (
            "sw_log_DeN_e",
            "sw_log_DeN_L_e",
            "sw_log_DeN_T_e",
        )
        forbidden["DenN"] = (
            ("sw_log_DenN_e",)
            if mapping.bundle.expected_source == "monte_carlo"
            else ("sw_log_DenN_L_e", "sw_log_DenN_T_e")
        )
    elif (
        mapping.closure.electron_transport == GEC_RESTRICTED_TRANSPORT_CLOSURE
        and mapping.bundle.expected_source == "two_term"
    ):
        expected.update(
            {
                "DeN": ("sw_log_DeN_e",),
                "muenN": ("sw_log_muenN_e",),
                "DenN": ("sw_log_DenN_e",),
            }
        )
        forbidden.update(
            {
                "DeN": (
                    "sw_log_DeN_L_e",
                    "sw_log_DeN_T_e",
                    "ptp.Er",
                    "ptp.Ez",
                ),
                "DenN": (
                    "sw_log_DenN_L_e",
                    "sw_log_DenN_T_e",
                    "ptp.Er",
                    "ptp.Ez",
                ),
            }
        )
    elif mapping.closure.electron_transport == GEC_RESTRICTED_TRANSPORT_CLOSURE:
        expected.update(
            {
                "DeN": ("sw_log_DeN_L_e", "sw_log_DeN_T_e"),
                "muenN": ("sw_log_muenN_e",),
                "DenN": ("sw_log_DenN_L_e", "sw_log_DenN_T_e"),
            }
        )
    return expected, forbidden


def _saved_closure_argument_audit(
    mapping: GecCcpMapping,
    model_xml: str,
    closure_support: dict[str, Any],
) -> dict[str, Any]:
    if not closure_support.get("passed"):
        return {
            "passed": False,
            "reason": "active closure tables have no common support",
        }
    if mapping.bundle.expected_source != "monte_carlo":
        support_minimum, support_maximum = closure_support[
            "common_intersection_mean_energy_eV"
        ]
        return _saved_sw_logeps_audit(
            model_xml,
            minimum_eV=support_minimum,
            maximum_eV=support_maximum,
        )
    expected_arguments: dict[str, list[float]] = {}
    transport_support = closure_support.get("transport_intersection_mean_energy_eV")
    if isinstance(transport_support, list):
        expected_arguments["sw_logeps_transport"] = list(
            closure_argument_range(mapping, transport_support)
        )
    for process_type, support in closure_support.get(
        "rate_supports_mean_energy_eV",
        {},
    ).items():
        expected_arguments[f"sw_logeps_rate_{process_type}"] = list(
            closure_argument_range(mapping, support)
        )
    elastic_support = closure_support.get("elastic_energy_loss_support_mean_energy_eV")
    if isinstance(elastic_support, list):
        expected_arguments["sw_logeps_el"] = list(
            closure_argument_range(mapping, elastic_support)
        )
    argument_audits = {
        name: _saved_sw_logeps_audit(
            model_xml,
            minimum_eV=float(support[0]),
            maximum_eV=float(support[1]),
            variable_name=name,
        )
        for name, support in expected_arguments.items()
    }
    return {
        "passed": bool(argument_audits)
        and all(item["passed"] for item in argument_audits.values()),
        "policy": "separate_transport_rate_and_elastic_energy_loss_guards",
        "arguments": argument_audits,
    }


def _saved_transport_property_bindings(
    mapping: GecCcpMapping,
    feature_xml: str,
    expected_tokens: dict[str, tuple[str, ...]],
    forbidden_tokens: dict[str, tuple[str, ...]],
) -> tuple[dict[str, dict[str, Any]], bool]:
    all_external_tokens = tuple(
        sorted(
            {
                f"sw_log_{item[0].removeprefix('sw_')}"
                for item in (
                    *GEC_RESTRICTED_TRANSPORT_FUNCTIONS,
                    *GEC_TWO_TERM_HYBRID_TRANSPORT_FUNCTIONS,
                )
            }
        )
    )
    input_contract = _transport_input_contract(
        mapping.closure,
        source=mapping.bundle.expected_source,
    )
    bindings: dict[str, dict[str, Any]] = {}
    all_passed = True
    for _, property_name, _, _, _ in GEC_TRANSPORT_PROPERTIES:
        matches = re.findall(
            rf'<param\b[^>]*param="{re.escape(property_name)}"'
            rf'[^>]*value="([^"]+)"',
            feature_xml,
        )
        expression = matches[-1] if matches else None
        tokens = expected_tokens.get(property_name, ())
        forbidden = forbidden_tokens.get(property_name, ())
        if tokens:
            passed = bool(
                expression
                and all(token in expression for token in tokens)
                and not any(token in expression for token in forbidden)
            )
        else:
            passed = not bool(
                expression and any(token in expression for token in all_external_tokens)
            )
        all_passed = all_passed and passed
        bindings[property_name] = {
            "input_source": input_contract[property_name]["source"],
            "active_external": input_contract[property_name]["source"].startswith(
                "external_"
            ),
            "expected_tokens": list(tokens),
            "forbidden_tokens": list(forbidden),
            "passed": passed,
        }
    return bindings, all_passed


def _saved_physics_sections(
    mapping: GecCcpMapping,
    model_xml: str,
    feature_xml: str,
) -> _SavedPhysicsSections:
    expected_flags = {
        "ReducedProps": "1",
        "TensorElectronProps": "0",
        "IncludeThermalDiffusion": (
            "1" if _thermal_diffusion_enabled(mapping.closure) else "0"
        ),
    }
    physics_flags: dict[str, dict[str, Any]] = {}
    for name, expected in expected_flags.items():
        try:
            actual = _property_values(model_xml, name)[-1]
        except GecCcpWorkflowError:
            actual = None
        physics_flags[name] = {
            "expected": expected,
            "actual": actual,
            "passed": actual == expected,
        }
    periodic_lock_counts = {
        identifier: len(
            re.findall(
                rf'<lock\b[^>]*param="{re.escape(identifier)}"',
                feature_xml,
            )
        )
        for identifier in ("ptp.Mn", "ptp.ebar", "ptp.Te")
    }
    periodic_restored = all(count == 0 for count in periodic_lock_counts.values())
    heavy_species, heavy_passed = _saved_heavy_species_ownership(model_xml)
    reaction_handling = (
        _audit_saved_reaction_handling(mapping, model_xml)
        if _uses_external_elastic_energy_loss(mapping.closure)
        else {"passed": True, "status": "not_applicable"}
    )
    passed = bool(
        all(item["passed"] for item in physics_flags.values())
        and periodic_restored
        and heavy_passed
        and reaction_handling["passed"]
    )
    return _SavedPhysicsSections(
        passed=passed,
        flags=physics_flags,
        periodic_restoration={
            "passed": periodic_restored,
            "saved_lock_counts": periodic_lock_counts,
            "expected_saved_lock_count": 0,
            "restored_definitions": "comsol_native",
        },
        heavy_species=heavy_species,
        reaction_handling=reaction_handling,
    )


def _saved_heavy_species_ownership(
    model_xml: str,
) -> tuple[dict[str, Any], bool]:
    ion_feature_match = re.search(
        r'<PhysicsFeature\b[^>]*op="Species"[^>]*tag="Ar_1p"[^>]*>'
        r"(.*?)</PhysicsFeature>",
        model_xml,
        flags=re.DOTALL,
    )
    ion_feature_lock_params: list[str] = []
    if ion_feature_match is not None:
        ion_feature_lock_params = re.findall(
            r'<lock\b[^>]*param="([^"]+)"',
            ion_feature_match.group(1),
        )
    native_ion_feature = bool(
        ion_feature_match is not None and not ion_feature_lock_params
    )
    custom_initialization_count = len(
        re.findall(
            r'<expressions\b[^>]*name="WAr_1p_per"',
            model_xml,
        )
    )
    native_ion_initialization = custom_initialization_count == 0
    obsolete_identifiers = (
        "ptp.dfluxr_wAr_1p",
        "ptp.dfluxz_wAr_1p",
        "ptp.ucr",
        "ptp.ucz",
    )
    obsolete_lock_counts = {
        identifier: len(
            re.findall(
                rf'<lock\b[^>]*param="{re.escape(identifier)}"',
                model_xml,
            )
        )
        for identifier in obsolete_identifiers
    }
    obsolete_locks_absent = all(count == 0 for count in obsolete_lock_counts.values())
    surface_features = re.findall(
        r'<PhysicsFeature\b[^>]*op="SurfaceReaction"[^>]*tag="sr1"[^>]*>'
        r"(.*?)</PhysicsFeature>",
        model_xml,
        flags=re.DOTALL,
    )
    surface_formula = None
    if len(surface_features) == 1:
        try:
            surface_formula = html.unescape(
                _property_values(surface_features[0], "formula")[-1]
            )
        except GecCcpWorkflowError:
            surface_formula = None
    surface_weak_override_count = len(
        re.findall(
            r'<lock\b[^>]*param="root\.comp1\.ptp\.sr1\.weak\$1"',
            model_xml,
        )
    )
    surface_retained = bool(
        len(surface_features) == 1
        and "DISABLED" not in surface_features[0]
        and surface_formula == "Ar+=>Ar"
        and surface_weak_override_count == 0
    )
    passed = bool(
        native_ion_feature
        and native_ion_initialization
        and obsolete_locks_absent
        and surface_retained
    )
    return {
        "passed": passed,
        "owner": "comsol_plasma_interface",
        "species_feature": "Ar_1p",
        "species_feature_present": ion_feature_match is not None,
        "species_feature_lock_params": ion_feature_lock_params,
        "custom_species_lock_count": len(ion_feature_lock_params),
        "custom_initialization_count": custom_initialization_count,
        "obsolete_flux_locks": {
            "passed": obsolete_locks_absent,
            "saved_lock_counts": obsolete_lock_counts,
        },
        "surface_reaction": {
            "tag": "sr1",
            "formula": surface_formula,
            "feature_count": len(surface_features),
            "weak_override_count": surface_weak_override_count,
            "handling": "retained_exactly_once",
            "passed": surface_retained,
        },
        "mass_fraction_mapping_override": False,
        "domain_weak_override": False,
        "boundary_weak_override": False,
        "custom_initialization": False,
    }, passed


def _audit_saved_transport_numeric(
    mapping: GecCcpMapping,
    closure_support: dict[str, Any],
) -> _SavedTransportNumericAudit:
    headers, values = _read_comsol_numeric_csv(
        mapping.results.output_directory / "swarm_tables" / "closure_phase_radial.csv"
    )
    comparisons: dict[str, Any] = {
        "status": "replaced_by_independent_bundle_value_audit",
        "comsol_expected_expression_columns_used": False,
    }
    domain_path = (
        mapping.results.output_directory / "swarm_tables" / "domain_phase_closure.csv"
    )
    domain_headers: list[str] = []
    domain_values = np.empty((0, 0))
    operating_range: dict[str, Any] = {
        "status": "unavailable",
        "passed": False,
        "reason": "domain phase closure export is absent",
    }
    if domain_path.exists():
        domain_headers, domain_values = _read_comsol_numeric_csv(domain_path)
        operating_range = _operating_mean_energy_support_audit(
            mapping,
            domain_headers,
            domain_values,
            closure_support,
        )
    orientation = _orientation_regularization_audit(
        mapping,
        domain_path_exists=domain_path.exists(),
        domain_headers=domain_headers,
        domain_values=domain_values,
    )
    if domain_path.exists():
        if mapping.results.role == "physical_target":
            transport_values = _independent_domain_transport_value_audit(
                mapping,
                domain_headers,
                domain_values,
                closure_support,
                qualified_low_endpoint_continuation=(
                    operating_range.get("support_acceptance_basis")
                    == "qualified_low_energy_guard_sensitivity"
                ),
            )
        else:
            transport_values = _independent_transport_value_audit(
                mapping,
                headers,
                values,
                domain_headers,
                domain_values,
                closure_support,
            )
        elastic_values = _independent_elastic_energy_loss_value_audit(
            mapping,
            domain_headers,
            domain_values,
            closure_support,
        )
        rate_values = _independent_external_rate_value_audit(
            mapping,
            domain_headers,
            domain_values,
            closure_support,
        )
    else:
        transport_values = {
            "passed": False,
            "reason": "domain phase closure export is absent",
        }
        elastic_values = (
            {
                "passed": False,
                "reason": "domain phase closure export is absent",
            }
            if _uses_external_elastic_energy_loss(mapping.closure)
            else {"passed": True, "status": "not_applicable"}
        )
        rate_values = (
            {
                "passed": False,
                "reason": "domain phase closure export is absent",
            }
            if mapping.closure.reaction_model == "external_rates"
            else {"passed": True, "status": "not_applicable"}
        )
    passed = bool(
        operating_range["passed"]
        and transport_values["passed"]
        and elastic_values["passed"]
        and rate_values["passed"]
    )
    return _SavedTransportNumericAudit(
        passed=passed,
        operating_range=operating_range,
        orientation_regularization=orientation,
        comparisons=comparisons,
        transport_values=transport_values,
        elastic_energy_loss_values=elastic_values,
        external_rate_values=rate_values,
    )


def _orientation_regularization_audit(
    mapping: GecCcpMapping,
    *,
    domain_path_exists: bool,
    domain_headers: list[str],
    domain_values: np.ndarray,
) -> dict[str, Any]:
    default = {"status": "not_applicable", "diagnostic_only": True}
    if (
        mapping.closure.electron_transport
        not in {
            GEC_RESTRICTED_TRANSPORT_CLOSURE,
            GEC_TWO_TERM_HYBRID_TRANSPORT_CLOSURE,
        }
        or mapping.bundle.expected_source != "monte_carlo"
    ):
        return default
    if not domain_path_exists:
        return {
            "status": "unavailable",
            "diagnostic_only": True,
            "reason": "domain phase closure export is absent",
        }

    def expression_indices(expression: str) -> list[int]:
        return [
            index
            for index, header in enumerate(domain_headers)
            if header == expression or header.startswith(expression + " ")
        ]

    neutral_indices = expression_indices("ptp.Nn")
    radial_indices = expression_indices("ptp.Er")
    axial_indices = expression_indices("ptp.Ez")
    if (
        not neutral_indices
        or len(neutral_indices) != len(radial_indices)
        or len(neutral_indices) != len(axial_indices)
    ):
        return {
            "status": "unavailable",
            "diagnostic_only": True,
            "reason": "Nn/Er/Ez phase columns are absent or inconsistent",
        }
    neutral = domain_values[:, neutral_indices].reshape(-1)
    radial = domain_values[:, radial_indices].reshape(-1)
    axial = domain_values[:, axial_indices].reshape(-1)
    finite = (
        np.isfinite(neutral)
        & np.isfinite(radial)
        & np.isfinite(axial)
        & (neutral > 0.0)
    )
    if not np.any(finite):
        return {
            "status": "unavailable",
            "diagnostic_only": True,
            "reason": "no finite positive-neutral-density samples",
        }
    reduced_field_Td = (
        np.hypot(radial[finite], axial[finite]) / neutral[finite] / 1.0e-21
    )
    scale_Td = float(mapping.closure.zero_field_isotropization_Td)

    def orientation_weight(scale: float) -> np.ndarray:
        squared = reduced_field_Td * reduced_field_Td
        return squared / (squared + scale * scale)

    nominal = orientation_weight(scale_Td)
    sensitivity: dict[str, Any] = {}
    for label, factor in (("half_scale", 0.5), ("double_scale", 2.0)):
        varied = orientation_weight(factor * scale_Td)
        delta = np.abs(varied - nominal)
        sensitivity[label] = {
            "scale_Td": factor * scale_Td,
            "mean_absolute_weight_change": float(np.mean(delta)),
            "maximum_absolute_weight_change": float(np.max(delta)),
        }
    return {
        "status": "evaluated",
        "diagnostic_only": True,
        "coverage": "unweighted exported plasma-domain nodes and RF phases",
        "scale_Td": scale_Td,
        "sample_count": int(reduced_field_Td.size),
        "reduced_field_Td": {
            "minimum": float(np.min(reduced_field_Td)),
            "median": float(np.median(reduced_field_Td)),
            "maximum": float(np.max(reduced_field_Td)),
        },
        "orientation_weight": {
            "minimum": float(np.min(nominal)),
            "median": float(np.median(nominal)),
            "maximum": float(np.max(nominal)),
            "fraction_below_0_5": float(np.mean(nominal < 0.5)),
            "fraction_below_0_9": float(np.mean(nominal < 0.9)),
            "fraction_below_0_99": float(np.mean(nominal < 0.99)),
        },
        "frozen_field_sensitivity": sensitivity,
    }


def _audit_saved_elastic_energy_loss_binding(
    mapping: GecCcpMapping,
    model_xml: str,
) -> dict[str, Any]:
    if not _uses_external_elastic_energy_loss(mapping.closure):
        return {"passed": True, "status": "not_applicable"}
    features: list[tuple[str, str, str]] = []
    for attributes, body in re.findall(
        r"<PhysicsFeature\b([^>]*)>(.*?)</PhysicsFeature>",
        model_xml,
        flags=re.DOTALL,
    ):
        operation = re.search(r'\bop="([^"]*)"', attributes)
        tag = re.search(r'\btag="([^"]*)"', attributes)
        if operation is None or operation.group(1) != "GeneralPowerDeposition":
            continue
        features.append((tag.group(1) if tag else "", attributes, body))
    selected = [item for item in features if item[0] == "swElLoss"]
    expected_qgen = (
        "-ptp.ne*ptp.n_wAr*exp(sw_logKel("
        + (
            "sw_logeps_el"
            if mapping.bundle.expected_source == "monte_carlo"
            else "sw_logeps"
        )
        + "))*1[eV*m^3/s]"
    )
    saved_qgen: str | None = None
    active = False
    domain_selection_valid = False
    saved_domain_entities: str | None = None
    if len(selected) == 1:
        _, attributes, body = selected[0]
        active = not bool(
            re.search(
                r"<entityFlags\b[^>]*>\s*DISABLED\s*</entityFlags>",
                body.split("<param", maxsplit=1)[0],
            )
        )
        try:
            saved_qgen = html.unescape(_property_values(body, "Qgen")[-1])
        except GecCcpWorkflowError:
            saved_qgen = None
        selection = re.search(
            r'<selection\b[^>]*selType="GEOMDIM"[^>]*>'
            r"(.*?)</selection>",
            body,
            flags=re.DOTALL,
        )
        explicit = (
            re.search(
                r'<explicit\b[^>]*dim="2"[^>]*geom="/geom/geom1"'
                r'[^>]*entities="([^"]+)"',
                selection.group(1),
            )
            if selection is not None
            else None
        )
        saved_domain_entities = explicit.group(1) if explicit else None
        domain_selection_valid = bool(
            explicit is not None and saved_domain_entities == "2,1"
        )
    weak_override_count = len(
        re.findall(
            r'<PhysicsFeature\b[^>]*op="WeakContribution"[^>]*tag="swElLoss"',
            model_xml,
        )
    )
    passed = bool(
        len(features) == 1
        and len(selected) == 1
        and active
        and domain_selection_valid
        and saved_qgen == expected_qgen
        and weak_override_count == 0
    )
    return {
        "passed": passed,
        "status": "active" if passed else "invalid",
        "feature": "swElLoss",
        "operation": "GeneralPowerDeposition",
        "feature_count": len(features),
        "selected_feature_count": len(selected),
        "active": active,
        "expected_Qgen": expected_qgen,
        "saved_Qgen": saved_qgen,
        "generated_domain_selection": [1],
        "saved_domain_entities_encoding": saved_domain_entities,
        "saved_domain_selection_valid": domain_selection_valid,
        "saved_selection_contract": (
            "runtime_COMSOL_API_asserts_exact_domain_1; saved XML also "
            "contains explicit geom1 entities=2,1"
        ),
        "weak_override_count": weak_override_count,
        "sign_convention": "negative_Qgen_is_electron_energy_loss",
        "coefficient_unit": "eV*m^3/s",
    }


def _audit_saved_reaction_handling(
    mapping: GecCcpMapping,
    model_xml: str,
) -> dict[str, Any]:
    """Audit active kinetics without weakening the original reaction nodes."""

    def property_or_none(block: str, name: str) -> str | None:
        try:
            return _property_values(block, name)[-1]
        except GecCcpWorkflowError:
            return None

    expected_counts = {"UseCrossSectionData": 0, "RateConstant": 0}
    actual_counts = {"UseCrossSectionData": 0, "RateConstant": 0}
    active_eedf_bindings = 0
    reactions: list[dict[str, Any]] = []
    external_rate_rows = (
        _read_csv(mapping.bundle.path / "rates_vs_mean_energy.csv")
        if mapping.closure.reaction_model == "external_rates"
        else []
    )
    for reaction in mapping.reactions:
        external_elastic = bool(
            reaction.process_type == "elastic"
            and _uses_external_elastic_energy_loss(mapping.closure)
        )
        preintegrated = _reaction_uses_preintegrated_rate(mapping.closure, reaction)
        expected_binding = (
            "DisabledElasticReaction"
            if external_elastic
            else ("RateConstant" if preintegrated else "UseCrossSectionData")
        )
        if not external_elastic:
            expected_counts[expected_binding] += 1
        feature = re.search(
            rf'<PhysicsFeature\b[^>]*op="ElectronImpactReaction"[^>]*tag="'
            rf'{re.escape(reaction.feature)}"[^>]*>(.*?)</PhysicsFeature>',
            model_xml,
            flags=re.DOTALL,
        )
        block = feature.group(1) if feature is not None else ""
        disabled = bool(
            re.search(
                r"<entityFlags\b[^>]*>\s*DISABLED\s*</entityFlags>",
                block.split("<param", maxsplit=1)[0],
            )
        )
        actual_binding = property_or_none(block, "SpecifyReactionUsing")
        if not disabled and actual_binding in actual_counts:
            actual_counts[actual_binding] += 1
        eedf_binding = property_or_none(block, "eedf")
        if (
            not disabled
            and actual_binding == "UseCrossSectionData"
            and eedf_binding == "FromPhysicsInterfaceProperty"
        ):
            active_eedf_bindings += 1
        rate_form = property_or_none(block, "RateConstantForm")
        rate_expression = property_or_none(block, "kf")
        energy_loss = property_or_none(block, "de")
        formula_value = property_or_none(block, "formula")
        formula_value = (
            html.unescape(formula_value) if formula_value is not None else None
        )
        physics_contract = GEC_REACTION_PHYSICS_CONTRACT[reaction.process_type]
        try:
            energy_loss_value = (
                float(energy_loss) if energy_loss is not None else math.nan
            )
        except ValueError:
            energy_loss_value = math.nan
        formula_valid = formula_value == physics_contract["formula"]
        energy_loss_valid = math.isclose(
            energy_loss_value,
            float(physics_contract["energy_loss_eV"]),
            rel_tol=0.0,
            abs_tol=1.0e-12,
        )
        binding_valid = (
            disabled if external_elastic else actual_binding == expected_binding
        )
        if external_elastic:
            pass
        elif preintegrated:
            expected_rate_token = (
                f"sw_knorm_{reaction.process_type}_e"
                if _uses_direct_mc_rates(mapping)
                else f"sw_logk_{reaction.process_type}_e"
            )
            expected_rate_expression: str | None = None
            if mapping.closure.reaction_model == "external_rates":
                selected_rates, rate_metadata = _active_rate_table(
                    mapping,
                    external_rate_rows,
                    process_type=reaction.process_type,
                )
                if _uses_direct_mc_rates(mapping):
                    rate_scale = float(rate_metadata["normalization_rate_m3_s"])
                    particle_rate_value = (
                        f"{rate_scale:.17e}*{expected_rate_token}("
                        f"sw_logeps_rate_{reaction.process_type})"
                    )
                else:
                    particle_rate_value = f"exp({expected_rate_token}(sw_logeps))"
                if len(selected_rates) < 2:
                    raise GecCcpWorkflowError(
                        "external reaction binding lacks two rate anchors"
                    )
                expected_rate_expression = (
                    f"{AVOGADRO_PER_MOL}*{particle_rate_value}*1[m^3/s]"
                )
            binding_valid = bool(
                binding_valid
                and rate_form == "UseRate"
                and rate_expression is not None
                and expected_rate_token in rate_expression
                and (
                    expected_rate_expression is None
                    or rate_expression == expected_rate_expression
                )
            )
        else:
            binding_valid = bool(
                binding_valid and eedf_binding == "FromPhysicsInterfaceProperty"
            )
        reaction_passed = bool(
            feature is not None
            and formula_valid
            and energy_loss_valid
            and binding_valid
        )
        reactions.append(
            {
                "feature": reaction.feature,
                "process_type": reaction.process_type,
                "node_operation": (
                    "ElectronImpactReaction" if feature is not None else None
                ),
                "expected_binding": expected_binding,
                "actual_binding": actual_binding,
                "active": not disabled,
                "expected_active": not external_elastic,
                "reaction_eedf_binding": eedf_binding,
                "rate_constant_form": rate_form if preintegrated else None,
                "rate_expression": rate_expression if preintegrated else None,
                "expected_formula": physics_contract["formula"],
                "saved_formula": formula_value,
                "formula_preserved": formula_valid,
                "energy_loss_de": energy_loss,
                "expected_energy_loss_eV": physics_contract["energy_loss_eV"],
                "energy_loss_preserved": energy_loss_valid,
                "energy_equation_owner": (
                    "GeneralPowerDeposition:swElLoss"
                    if external_elastic
                    else "ElectronImpactReaction:de"
                ),
                "passed": reaction_passed,
            }
        )
    elastic_energy_loss = _audit_saved_elastic_energy_loss_binding(mapping, model_xml)
    return {
        "passed": (
            actual_counts == expected_counts
            and active_eedf_bindings == expected_counts["UseCrossSectionData"]
            and all(item["passed"] for item in reactions)
            and elastic_energy_loss["passed"]
        ),
        "expected_binding_counts": expected_counts,
        "actual_binding_counts": actual_counts,
        "active_reaction_eedf_bindings": active_eedf_bindings,
        "elastic_energy_loss": elastic_energy_loss,
        "reactions": reactions,
    }
