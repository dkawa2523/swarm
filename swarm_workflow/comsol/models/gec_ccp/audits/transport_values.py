"""Independently reconstruct GEC-CCP transport, rate, and elastic values."""

from __future__ import annotations

import math
from typing import Any, Literal

import numpy as np
from scipy.interpolate import LinearNDInterpolator, PchipInterpolator

from electron_swarm.core.constants import E_CHARGE_C

from swarm_workflow.comsol.input.mean_energy import MeanEnergyArgument
from swarm_workflow.comsol.models.gec_ccp.validation.bundle_guards import (
    GEC_ELASTIC_ENERGY_LOSS_COLUMN,
    GEC_ELASTIC_ENERGY_LOSS_TABLE,
)
from swarm_workflow.comsol.models.gec_ccp.closure import (
    GEC_AXISYMMETRIC_ACTIVE_TRANSPORT_COMPONENTS,
    GEC_AXISYMMETRIC_ACTIVE_TRANSPORT_INDICES,
    GEC_TENSOR_COMPONENTS,
    _active_rate_table,
    _active_transport_functions,
    _reaction_uses_preintegrated_rate,
    _transport_audit_components,
    _uses_direct_mc_rates,
    _uses_external_elastic_energy_loss,
    _uses_hybrid_einstein_transport,
)
from swarm_workflow.comsol.models.gec_ccp.contracts import (
    FUNCTION_EEDF_PREINTEGRATED_INELASTIC,
    GEC_RESTRICTED_TRANSPORT_CLOSURE,
    GecCcpMapping,
    GecCcpWorkflowError,
)
from swarm_workflow.comsol.models.gec_ccp.data import _lookup, _read_csv
from swarm_workflow.comsol.models.gec_ccp.closure_arguments import (
    closure_argument_range,
)
from swarm_workflow.comsol.models.gec_ccp.validation.joint_consistency import (
    _build_dense_function_eedf_rate_closure,
    _independent_log_piecewise_cubic_transport,
)


GEC_INDEPENDENT_TRANSPORT_RELATIVE_TOLERANCE = 2.0e-2


def _phase_expression_columns(
    headers: list[str], expression: str
) -> tuple[list[str], list[int]]:
    labels: list[str] = []
    indices: list[int] = []
    for index, header in enumerate(headers):
        # COMSOL uses semicolons for function arguments in CSV headers so the
        # header remains comma-delimited. Normalize only for expression matching.
        normalized = header.replace(";", ",")
        if normalized != expression and not normalized.startswith(expression + " "):
            continue
        labels.append(
            normalized.split(" @ ", maxsplit=1)[1]
            if " @ " in normalized
            else f"column_{index}"
        )
        indices.append(index)
    return labels, indices


def _independent_elastic_energy_loss_value_audit(
    mapping: GecCcpMapping,
    headers: list[str],
    values: np.ndarray,
    support: dict[str, Any],
) -> dict[str, Any]:
    """Recompute Qgen from the bundle without a saved COMSOL function."""

    if not _uses_external_elastic_energy_loss(mapping.closure):
        return {"passed": True, "status": "not_applicable"}
    phase_labels, energy_indices = _phase_expression_columns(headers, "ptp.ebar")
    matrices: dict[str, np.ndarray] = {}
    elastic_argument = (
        "sw_logeps_el"
        if mapping.bundle.expected_source == "monte_carlo"
        else "sw_logeps"
    )
    exported_qgen_expression = (
        f"-ptp.ne*ptp.n_wAr*exp(sw_logKel({elastic_argument}))*1[eV*m^3/s]"
    )
    for name, expression in (
        ("electron_density", "ptp.ne"),
        ("argon_neutral_density", "ptp.n_wAr"),
        ("Qgen", exported_qgen_expression),
    ):
        labels, indices = _phase_expression_columns(headers, expression)
        if not phase_labels or labels != phase_labels:
            return {
                "passed": False,
                "reason": f"domain export lacks aligned {expression} columns",
            }
        matrices[name] = values[:, indices]
    energy = values[:, energy_indices]
    electron_density = matrices["electron_density"]
    argon_density = matrices["argon_neutral_density"]
    actual = matrices["Qgen"]
    finite = bool(
        np.all(np.isfinite(energy))
        and np.all(np.isfinite(electron_density))
        and np.all(np.isfinite(argon_density))
        and np.all(np.isfinite(actual))
        and np.all(energy > 0.0)
        and np.all(electron_density > 0.0)
        and np.all(argon_density > 0.0)
    )
    elastic_support = support.get("elastic_energy_loss_support_mean_energy_eV")
    if not finite or not isinstance(elastic_support, list):
        return {
            "passed": False,
            "reason": "elastic energy-loss state/support is invalid",
        }
    rows = _read_csv(mapping.bundle.path / GEC_ELASTIC_ENERGY_LOSS_TABLE)
    argument_minimum, argument_maximum = closure_argument_range(
        mapping, elastic_support
    )
    coefficient = _independent_log_piecewise_cubic_transport(
        rows,
        GEC_ELASTIC_ENERGY_LOSS_COLUMN,
        energy,
        support_minimum_eV=argument_minimum,
        support_maximum_eV=argument_maximum,
    )
    expected = -electron_density * argon_density * coefficient * E_CHARGE_C
    error = np.abs(actual - expected)
    allowed = GEC_INDEPENDENT_TRANSPORT_RELATIVE_TOLERANCE * np.abs(
        expected
    ) + 1.0e-12 * max(float(np.max(np.abs(expected))), 1.0e-300)
    normalized = error / np.maximum(allowed, 1.0e-300)
    passed = bool(
        np.all(np.isfinite(coefficient))
        and np.all(coefficient > 0.0)
        and np.all(actual <= 0.0)
        and np.max(normalized) <= 1.0
    )
    return {
        "passed": passed,
        "status": "evaluated",
        "method": "independent_bundle_log_PCHIP",
        "argument_range_eV": [argument_minimum, argument_maximum],
        "formula_SI": "-ne*n_wAr*K_epsilon_el*elementary_charge",
        "sign_convention": "nonpositive_Qgen_is_electron_energy_loss",
        "relative_tolerance": GEC_INDEPENDENT_TRANSPORT_RELATIVE_TOLERANCE,
        "phase_count": len(phase_labels),
        "spatial_rows": int(values.shape[0]),
        "maximum_relative_error": float(
            np.max(error / np.maximum(np.abs(expected), 1.0e-300))
        ),
        "maximum_normalized_error": float(np.max(normalized)),
        "actual_Qgen_W_m3": [float(np.min(actual)), float(np.max(actual))],
        "expected_Qgen_W_m3": [float(np.min(expected)), float(np.max(expected))],
    }


def _independent_rate_coefficient_values(
    mapping: GecCcpMapping,
    process_type: str,
    mean_energy_eV: np.ndarray,
    support: list[float],
) -> np.ndarray:
    if mapping.closure.reaction_model == FUNCTION_EEDF_PREINTEGRATED_INELASTIC:
        rows, _ = _build_dense_function_eedf_rate_closure(
            mapping, input_mph=mapping.model.input_mph
        )
    else:
        rows = _read_csv(mapping.bundle.path / "rates_vs_mean_energy.csv")
    selected, metadata = _active_rate_table(
        mapping,
        rows,
        process_type=process_type,
    )
    x_values, y_values = _lookup(
        selected,
        "mean_energy_eV",
        "rate_coefficient_m3_s",
    )
    smooth = MeanEnergyArgument(*closure_argument_range(mapping, support)).log_values(
        mean_energy_eV
    )
    x = np.log(np.asarray(x_values, dtype=float))
    y = np.asarray(y_values, dtype=float)
    if _uses_direct_mc_rates(mapping):
        scale = float(metadata["normalization_rate_m3_s"])
        return scale * PchipInterpolator(
            x,
            y / scale,
            extrapolate=False,
        )(smooth)
    if np.any(y <= 0.0):
        raise GecCcpWorkflowError(f"positive external rate expected for {process_type}")
    return np.exp(PchipInterpolator(x, np.log(y), extrapolate=False)(smooth))


def _independent_external_rate_value_audit(
    mapping: GecCcpMapping,
    headers: list[str],
    values: np.ndarray,
    support: dict[str, Any],
) -> dict[str, Any]:
    """Compare saved COMSOL molar kf with independent bundle evaluation."""

    if mapping.closure.reaction_model != "external_rates":
        return {"passed": True, "status": "not_applicable"}
    phase_labels, energy_indices = _phase_expression_columns(headers, "ptp.ebar")
    if not phase_labels:
        return {"passed": False, "reason": "rate audit lacks mean energy"}
    energy = values[:, energy_indices]
    if np.any(~np.isfinite(energy)) or np.any(energy <= 0.0):
        return {"passed": False, "reason": "rate audit mean energy is invalid"}
    rate_supports = support.get("rate_supports_mean_energy_eV", {})
    process_results: dict[str, Any] = {}
    passed = True
    for reaction_index, reaction in enumerate(mapping.reactions, start=1):
        if not _reaction_uses_preintegrated_rate(mapping.closure, reaction):
            continue
        labels, rate_indices = _phase_expression_columns(
            headers, f"ptp.kf_{reaction_index}"
        )
        process_support = rate_supports.get(reaction.process_type)
        if labels != phase_labels or not isinstance(process_support, list):
            process_results[reaction.process_type] = {
                "passed": False,
                "reason": "aligned rate columns or support are absent",
            }
            passed = False
            continue
        actual_molar = values[:, rate_indices]
        particle_rate = _independent_rate_coefficient_values(
            mapping,
            reaction.process_type,
            energy,
            process_support,
        )
        expected_molar = 6.02214076e23 * particle_rate
        error = np.abs(actual_molar - expected_molar)
        scale = max(float(np.max(np.abs(expected_molar))), 1.0e-300)
        allowed = (
            GEC_INDEPENDENT_TRANSPORT_RELATIVE_TOLERANCE * np.abs(expected_molar)
            + 1.0e-10 * scale
        )
        normalized = error / np.maximum(allowed, 1.0e-300)
        process_passed = bool(
            np.all(np.isfinite(actual_molar))
            and np.all(np.isfinite(expected_molar))
            and np.all(actual_molar >= 0.0)
            and np.all(expected_molar >= 0.0)
            and np.max(normalized) <= 1.0
        )
        passed = passed and process_passed
        process_results[reaction.process_type] = {
            "passed": process_passed,
            "feature": reaction.feature,
            "exported_expression": f"ptp.kf_{reaction_index}",
            "COMSOL_unit": "m^3/(mol*s)",
            "bundle_unit": "m^3/s",
            "Avogadro_factor_per_mol": 6.02214076e23,
            "phase_count": len(phase_labels),
            "spatial_rows": int(values.shape[0]),
            "maximum_relative_error": float(
                np.max(error / np.maximum(np.abs(expected_molar), 1.0e-300))
            ),
            "maximum_normalized_error": float(np.max(normalized)),
        }
    return {
        "passed": passed,
        "status": "evaluated",
        "method": "independent_bundle_PCHIP_and_Avogadro_conversion",
        "relative_tolerance": GEC_INDEPENDENT_TRANSPORT_RELATIVE_TOLERANCE,
        "processes": process_results,
    }


def _interpolate_domain_state_to_radial_cut(
    domain_coordinates: np.ndarray,
    domain_state: np.ndarray,
    radial_coordinates: np.ndarray,
) -> tuple[np.ndarray, str, float]:
    """Project phase-local domain state independently onto the radial cut."""

    unique_coordinates, unique_indices = np.unique(
        domain_coordinates, axis=0, return_index=True
    )
    unique_state = domain_state[unique_indices]
    projected = np.full((len(radial_coordinates), domain_state.shape[1]), np.nan)
    method = "nearest_exact_only"
    if (
        len(unique_coordinates) >= 3
        and np.linalg.matrix_rank(
            unique_coordinates - np.mean(unique_coordinates, axis=0)
        )
        == 2
    ):
        try:
            projected = np.asarray(
                LinearNDInterpolator(
                    unique_coordinates, unique_state, fill_value=np.nan
                )(radial_coordinates),
                dtype=float,
            )
            method = "linear_nd"
        except (ValueError, RuntimeError):
            pass
    distances = np.linalg.norm(
        radial_coordinates[:, None, :] - unique_coordinates[None, :, :],
        axis=2,
    )
    nearest = np.argmin(distances, axis=1)
    nearest_distance = distances[np.arange(len(radial_coordinates)), nearest]
    coordinate_span = max(
        float(np.ptp(unique_coordinates[:, 0])),
        float(np.ptp(unique_coordinates[:, 1])),
        1.0e-300,
    )
    exact = nearest_distance <= max(1.0e-12, 1.0e-10 * coordinate_span)
    missing = ~np.all(np.isfinite(projected), axis=1)
    projected[missing & exact] = unique_state[nearest[missing & exact]]
    if np.any(missing & exact):
        method += "+exact_node"
    maximum_normalized_distance = (
        float(np.max(nearest_distance[missing & exact] / coordinate_span))
        if np.any(missing & exact)
        else 0.0
    )
    return projected, method, maximum_normalized_distance


def _independent_isotropic_transport_tensor(
    reduced: np.ndarray,
    neutral_density_m3: np.ndarray,
    *,
    include_azimuthal: bool = True,
) -> np.ndarray:
    zero = np.zeros_like(reduced)
    phiphi = reduced if include_azimuthal else zero
    return (
        np.stack(
            [reduced, zero, zero, zero, phiphi, zero, zero, zero, reduced],
            axis=1,
        )
        / neutral_density_m3[:, None]
    )


def _independent_field_aligned_transport_tensor(
    longitudinal: np.ndarray,
    transverse: np.ndarray,
    neutral_density_m3: np.ndarray,
    radial_field_V_m: np.ndarray,
    axial_field_V_m: np.ndarray,
    *,
    zero_field_isotropization_Td: float,
) -> np.ndarray:
    field_squared = radial_field_V_m**2 + axial_field_V_m**2
    denominator = (
        field_squared
        + (zero_field_isotropization_Td * 1.0e-21 * neutral_density_m3) ** 2
    )
    trace_part = field_squared / (3.0 * denominator)
    delta = longitudinal - transverse
    isotropic = (longitudinal + 2.0 * transverse) / 3.0
    rr = isotropic + delta * (radial_field_V_m**2 / denominator - trace_part)
    phiphi = isotropic - delta * trace_part
    zz = isotropic + delta * (axial_field_V_m**2 / denominator - trace_part)
    rz = delta * radial_field_V_m * axial_field_V_m / denominator
    zero = np.zeros_like(longitudinal)
    return (
        np.stack([rr, zero, rz, zero, phiphi, zero, rz, zero, zz], axis=1)
        / neutral_density_m3[:, None]
    )


def _independent_transport_tensors(
    mapping: GecCcpMapping,
    mean_energy_eV: np.ndarray,
    neutral_density_m3: np.ndarray,
    radial_field_V_m: np.ndarray,
    axial_field_V_m: np.ndarray,
    support_mean_energy_eV: list[float],
) -> dict[str, np.ndarray]:
    """Reconstruct local COMSOL transport tensors from the bundle only."""

    energy = np.asarray(mean_energy_eV, dtype=float).reshape(-1)
    neutral = np.asarray(neutral_density_m3, dtype=float).reshape(-1)
    radial_field = np.asarray(radial_field_V_m, dtype=float).reshape(-1)
    axial_field = np.asarray(axial_field_V_m, dtype=float).reshape(-1)
    if not (energy.shape == neutral.shape == radial_field.shape == axial_field.shape):
        raise GecCcpWorkflowError(
            "independent transport reconstruction requires aligned states"
        )
    rows = _read_csv(mapping.bundle.path / "transport_vs_mean_energy.csv")
    reduced: dict[str, np.ndarray] = {}
    argument_minimum, argument_maximum = closure_argument_range(
        mapping, support_mean_energy_eV
    )
    for _, name, column, _ in _active_transport_functions(
        mapping.closure, source=mapping.bundle.expected_source
    ):
        reduced[name] = _independent_log_piecewise_cubic_transport(
            rows,
            column,
            energy,
            support_minimum_eV=argument_minimum,
            support_maximum_eV=argument_maximum,
        )

    tensors = {
        "muN": _independent_isotropic_transport_tensor(
            reduced["muN"],
            neutral,
            include_azimuthal=(
                mapping.closure.electron_transport != "swarm_mobility_einstein"
            ),
        )
    }
    if _uses_hybrid_einstein_transport(mapping.closure):
        tensors["DeN"] = _independent_isotropic_transport_tensor(
            reduced["muN"] * (2.0 / 3.0) * energy, neutral
        )
        tensors["muenN"] = _independent_isotropic_transport_tensor(
            reduced["muenN"], neutral
        )
        if mapping.bundle.expected_source == "monte_carlo":
            tensors["DenN"] = _independent_field_aligned_transport_tensor(
                reduced["DenN_L"],
                reduced["DenN_T"],
                neutral,
                radial_field,
                axial_field,
                zero_field_isotropization_Td=float(
                    mapping.closure.zero_field_isotropization_Td or 0.0
                ),
            )
        else:
            tensors["DenN"] = _independent_isotropic_transport_tensor(
                reduced["DenN"], neutral
            )
    elif (
        mapping.closure.electron_transport == GEC_RESTRICTED_TRANSPORT_CLOSURE
        and mapping.bundle.expected_source == "two_term"
    ):
        for quantity, reduced_name in (
            ("DeN", "DeN"),
            ("muenN", "muenN"),
            ("DenN", "DenN"),
        ):
            tensors[quantity] = _independent_isotropic_transport_tensor(
                reduced[reduced_name], neutral
            )
    elif mapping.closure.electron_transport == GEC_RESTRICTED_TRANSPORT_CLOSURE:
        tensors["muenN"] = _independent_isotropic_transport_tensor(
            reduced["muenN"], neutral
        )
        scale = float(mapping.closure.zero_field_isotropization_Td or 0.0)
        for quantity, longitudinal_name, transverse_name in (
            ("DeN", "DeN_L", "DeN_T"),
            ("DenN", "DenN_L", "DenN_T"),
        ):
            tensors[quantity] = _independent_field_aligned_transport_tensor(
                reduced[longitudinal_name],
                reduced[transverse_name],
                neutral,
                radial_field,
                axial_field,
                zero_field_isotropization_Td=scale,
            )
    return tensors


def _independent_transport_value_audit(
    mapping: GecCcpMapping,
    radial_headers: list[str],
    radial_values: np.ndarray,
    domain_headers: list[str],
    domain_values: np.ndarray,
    support: dict[str, Any],
) -> dict[str, Any]:
    """Check representative saved coefficients against the bundle itself."""

    if not support.get("passed"):
        return {"passed": False, "reason": "common support is unavailable"}
    intersection = support["common_intersection_mean_energy_eV"]
    phase_labels, energy_indices = _phase_expression_columns(radial_headers, "ptp.ebar")
    neutral_labels, neutral_indices = _phase_expression_columns(
        radial_headers, "ptp.Nn"
    )
    if not phase_labels or neutral_labels != phase_labels:
        return {
            "passed": False,
            "reason": "radial export lacks aligned ebar/Nn phase columns",
        }
    energy = radial_values[:, energy_indices]
    neutral = radial_values[:, neutral_indices]

    actual_matrices: dict[tuple[str, str], np.ndarray] = {}
    for quantity, component, expression, _, _ in _transport_audit_components(
        mapping.closure, source=mapping.bundle.expected_source
    ):
        labels, indices = _phase_expression_columns(radial_headers, expression)
        if labels != phase_labels:
            return {
                "passed": False,
                "reason": f"radial export lacks aligned {expression} columns",
            }
        actual_matrices[(quantity, component)] = radial_values[:, indices]

    radial_coordinates = radial_values[:, :2]
    domain_coordinates = domain_values[:, :2]
    state_matrices: list[np.ndarray] = []
    for expression in ("ptp.Nn", "ptp.Er", "ptp.Ez", "ptp.ebar"):
        labels, indices = _phase_expression_columns(domain_headers, expression)
        if labels != phase_labels:
            return {
                "passed": False,
                "reason": f"domain export lacks aligned {expression} columns",
            }
        state_matrices.append(domain_values[:, indices])
    domain_state = np.concatenate(state_matrices, axis=1)
    projected, spatial_method, nearest_distance = (
        _interpolate_domain_state_to_radial_cut(
            domain_coordinates, domain_state, radial_coordinates
        )
    )
    phase_count = len(phase_labels)
    domain_neutral, radial_er, radial_ez, domain_energy = np.split(
        projected, [phase_count, 2 * phase_count, 3 * phase_count], axis=1
    )
    finite_state = (
        np.isfinite(energy)
        & np.isfinite(neutral)
        & np.isfinite(domain_neutral)
        & np.isfinite(radial_er)
        & np.isfinite(radial_ez)
        & np.isfinite(domain_energy)
        & (energy > 0.0)
        & (neutral > 0.0)
        & (domain_neutral > 0.0)
    )
    state_relative_error = np.maximum(
        np.abs(domain_neutral - neutral) / np.maximum(neutral, 1.0e-300),
        np.abs(domain_energy - energy) / np.maximum(energy, 1.0e-300),
    )
    finite_state &= state_relative_error <= (
        GEC_INDEPENDENT_TRANSPORT_RELATIVE_TOLERANCE
    )
    for matrix in actual_matrices.values():
        finite_state &= np.isfinite(matrix)

    source_is_mc = mapping.bundle.expected_source == "monte_carlo"
    if source_is_mc:
        finite_fields = np.concatenate(
            (radial_er[np.isfinite(radial_er)], radial_ez[np.isfinite(radial_ez)])
        )
        if finite_fields.size == 0:
            return {
                "passed": False,
                "reason": "domain projection has no finite electric field",
                "spatial_projection": spatial_method,
            }
        field_scale = max(float(np.max(np.abs(finite_fields))), 1.0e-300)
        mixed_field = (np.abs(radial_er) > 1.0e-10 * field_scale) & (
            np.abs(radial_ez) > 1.0e-10 * field_scale
        )
        candidates = finite_state & mixed_field
        score = np.abs(radial_er * radial_ez)
    else:
        candidates = finite_state
        midpoint = 0.5 * (math.log(intersection[0]) + math.log(intersection[1]))
        score = -np.abs(np.log(np.maximum(energy, 1.0e-300)) - midpoint)
    flat_candidates = np.flatnonzero(candidates)
    if flat_candidates.size == 0:
        return {
            "passed": False,
            "reason": (
                "no phase-local radial point has independently aligned "
                + ("nonzero Er and Ez" if source_is_mc else "domain state")
            ),
            "spatial_projection": spatial_method,
        }
    ranked = flat_candidates[np.argsort(score.reshape(-1)[flat_candidates])[::-1]][:1]
    row_indices, phase_indices = np.unravel_index(ranked, energy.shape)
    sample_energy = energy[row_indices, phase_indices]
    sample_neutral = neutral[row_indices, phase_indices]
    sample_er = radial_er[row_indices, phase_indices]
    sample_ez = radial_ez[row_indices, phase_indices]

    expected_tensors = _independent_transport_tensors(
        mapping,
        sample_energy,
        sample_neutral,
        sample_er,
        sample_ez,
        intersection,
    )

    quantity_results: dict[str, Any] = {}
    passed = True
    for quantity, expected in expected_tensors.items():
        if (
            quantity == "muN"
            and mapping.closure.electron_transport == "swarm_mobility_einstein"
        ):
            quantity_passed = bool(
                np.all(np.isfinite(expected)) and np.all(expected[:, (0, 8)] > 0.0)
            )
            passed = passed and quantity_passed
            quantity_results[quantity] = {
                "passed": quantity_passed,
                "comparison": (
                    "saved_SpecifyMueOnly_binding_plus_independent_bundle_evaluation"
                ),
                "postprocessed_ptp_mue_tensor_compared": False,
                "reason": (
                    "the TimePeriodic phase dataset exposes COMSOL's native "
                    "scalar-property postprocess, not the feature-bound muN "
                    "expression used to assemble the solved residual"
                ),
                "active_components": ["rr", "zz"],
                "minimum_expected_mobility_m2_V_s": float(np.min(expected[:, (0, 8)])),
                "maximum_expected_mobility_m2_V_s": float(np.max(expected[:, (0, 8)])),
            }
            continue
        actual = np.stack(
            [
                actual_matrices[(quantity, component)][row_indices, phase_indices]
                for component in GEC_AXISYMMETRIC_ACTIVE_TRANSPORT_COMPONENTS
            ],
            axis=1,
        )
        equation_active_expected = expected[
            :, GEC_AXISYMMETRIC_ACTIVE_TRANSPORT_INDICES
        ]
        scale = max(float(np.max(np.abs(equation_active_expected))), 1.0e-300)
        error = np.abs(actual - equation_active_expected)
        allowed = (
            GEC_INDEPENDENT_TRANSPORT_RELATIVE_TOLERANCE
            * np.abs(equation_active_expected)
            + 1.0e-10 * scale
        )
        normalized = error / np.maximum(allowed, 1.0e-300)
        quantity_passed = bool(
            np.all(np.isfinite(actual)) and np.max(normalized) <= 1.0
        )
        passed = passed and quantity_passed
        quantity_results[quantity] = {
            "passed": quantity_passed,
            "tensor_components_checked": list(
                GEC_AXISYMMETRIC_ACTIVE_TRANSPORT_COMPONENTS
            ),
            "maximum_relative_error": float(
                np.max(error / np.maximum(np.abs(equation_active_expected), 1.0e-300))
            ),
            "maximum_normalized_error": float(np.max(normalized)),
        }
    return {
        "passed": passed,
        "method": "python_log_PCHIP_piecewise_cubic_bundle_interpolation",
        "shape_preserving_interpolator": "scipy_PchipInterpolator",
        "comsol_expected_expression_columns_used": False,
        "tensor_component_scope": "2d_axisymmetric_equation_active_rz_block",
        "equation_active_tensor_components": list(
            GEC_AXISYMMETRIC_ACTIVE_TRANSPORT_COMPONENTS
        ),
        "stored_but_not_equation_active_tensor_components": [
            component
            for component in GEC_TENSOR_COMPONENTS
            if component not in GEC_AXISYMMETRIC_ACTIVE_TRANSPORT_COMPONENTS
        ],
        "relative_tolerance": GEC_INDEPENDENT_TRANSPORT_RELATIVE_TOLERANCE,
        "spatial_projection": spatial_method,
        "maximum_exact_fallback_distance_over_domain_span": nearest_distance,
        "representative_samples": [
            {
                "R_m": float(radial_coordinates[row, 0]),
                "Z_m": float(radial_coordinates[row, 1]),
                "phase": phase_labels[phase],
                "mean_energy_eV": float(sample_energy[index]),
                "neutral_density_m3": float(sample_neutral[index]),
                "Er_V_m": float(sample_er[index]),
                "Ez_V_m": float(sample_ez[index]),
            }
            for index, (row, phase) in enumerate(
                zip(row_indices, phase_indices, strict=True)
            )
        ],
        "mixed_Er_Ez_sample_required": source_is_mc,
        "quantities": quantity_results,
    }


def _independent_domain_transport_value_audit(
    mapping: GecCcpMapping,
    domain_headers: list[str],
    domain_values: np.ndarray,
    support: dict[str, Any],
    *,
    qualified_low_endpoint_continuation: bool = False,
) -> dict[str, Any]:
    """Fail closed on every exported domain node and RF phase."""

    if not support.get("passed"):
        return {
            "passed": False,
            "status": "unavailable",
            "reason": "common support is unavailable",
            "acceptance_scope": "full_domain_all_exported_rf_phases",
        }
    intersection = support["common_intersection_mean_energy_eV"]
    phase_labels, energy_indices = _phase_expression_columns(domain_headers, "ptp.ebar")
    if not phase_labels:
        return {
            "passed": False,
            "status": "unavailable",
            "reason": "domain export lacks phase-resolved ptp.ebar",
            "acceptance_scope": "full_domain_all_exported_rf_phases",
        }

    state_matrices: dict[str, np.ndarray] = {}
    for expression in ("ptp.Nn", "ptp.Er", "ptp.Ez"):
        labels, indices = _phase_expression_columns(domain_headers, expression)
        if labels != phase_labels:
            return {
                "passed": False,
                "status": "unavailable",
                "reason": f"domain export lacks aligned {expression} columns",
                "acceptance_scope": "full_domain_all_exported_rf_phases",
            }
        state_matrices[expression] = domain_values[:, indices]
    energy = domain_values[:, energy_indices]
    neutral = state_matrices["ptp.Nn"]
    radial_field = state_matrices["ptp.Er"]
    axial_field = state_matrices["ptp.Ez"]

    audit_components = _transport_audit_components(
        mapping.closure, source=mapping.bundle.expected_source
    )
    actual_matrices: dict[tuple[str, str], np.ndarray] = {}
    missing_expressions: list[str] = []
    for quantity, component, expression, _, _ in audit_components:
        labels, indices = _phase_expression_columns(domain_headers, expression)
        if labels != phase_labels:
            missing_expressions.append(expression)
            continue
        actual_matrices[(quantity, component)] = domain_values[:, indices]
    expected_component_count = len(audit_components)
    if missing_expressions:
        return {
            "passed": False,
            "status": "unavailable",
            "reason": "domain export lacks equation-active transport columns",
            "missing_expressions": missing_expressions,
            "acceptance_scope": "full_domain_all_exported_rf_phases",
            "expected_equation_active_component_count": expected_component_count,
            "exported_equation_active_component_count": len(actual_matrices),
        }

    coordinates = np.asarray(domain_values[:, :2], dtype=float)
    coordinate_finite = np.all(np.isfinite(coordinates), axis=1)[:, None]
    state_finite = (
        coordinate_finite
        & np.isfinite(energy)
        & np.isfinite(neutral)
        & np.isfinite(radial_field)
        & np.isfinite(axial_field)
        & (energy > 0.0)
        & (neutral > 0.0)
    )
    support_covered = (energy >= float(intersection[0])) & (
        energy <= float(intersection[1])
    )
    below_support = energy < float(intersection[0])
    above_support = energy > float(intersection[1])
    accepted_support = support_covered | (
        below_support & qualified_low_endpoint_continuation
    )
    comparable = state_finite & accepted_support
    total_sample_count = int(energy.size)
    state_finite_count = int(np.count_nonzero(state_finite))
    support_covered_count = int(np.count_nonzero(state_finite & support_covered))
    evaluated_sample_count = int(np.count_nonzero(comparable))
    full_state_coverage = bool(
        total_sample_count > 0 and evaluated_sample_count == total_sample_count
    )
    coverage = {
        "spatial_node_count": int(domain_values.shape[0]),
        "rf_phase_count": len(phase_labels),
        "required_sample_count": total_sample_count,
        "finite_positive_state_count": state_finite_count,
        "within_bundle_support_count": support_covered_count,
        "qualified_low_endpoint_continuation_count": int(
            np.count_nonzero(
                state_finite & below_support & qualified_low_endpoint_continuation
            )
        ),
        "above_bundle_support_count": int(
            np.count_nonzero(state_finite & above_support)
        ),
        "evaluated_sample_count": evaluated_sample_count,
        "fraction": (
            evaluated_sample_count / total_sample_count if total_sample_count else 0.0
        ),
        "complete": full_state_coverage,
        "qualified_low_endpoint_continuation": (qualified_low_endpoint_continuation),
    }
    if evaluated_sample_count == 0:
        return {
            "passed": False,
            "status": "evaluated",
            "reason": "no finite positive domain state lies within support",
            "acceptance_scope": "full_domain_all_exported_rf_phases",
            "coverage": coverage,
            "expected_equation_active_component_count": expected_component_count,
            "exported_equation_active_component_count": len(audit_components),
        }

    selected_flat = np.flatnonzero(comparable)
    row_indices, phase_indices = np.unravel_index(selected_flat, energy.shape)
    sample_energy = energy[row_indices, phase_indices]
    sample_neutral = neutral[row_indices, phase_indices]
    sample_er = radial_field[row_indices, phase_indices]
    sample_ez = axial_field[row_indices, phase_indices]
    expected_tensors = _independent_transport_tensors(
        mapping,
        sample_energy,
        sample_neutral,
        sample_er,
        sample_ez,
        intersection,
    )

    def finite_number(value: float) -> float | None:
        return float(value) if math.isfinite(float(value)) else None

    def finite_statistic(
        values: np.ndarray, statistic: Literal["maximum", "p95"]
    ) -> float | None:
        finite = values[np.isfinite(values)]
        if finite.size == 0:
            return None
        if statistic == "maximum":
            return float(np.max(finite))
        return float(np.percentile(finite, 95.0))

    quantities: dict[str, Any] = {}
    passed = full_state_coverage
    global_worst: dict[str, Any] | None = None
    global_worst_score = -math.inf
    all_finite_normalized: list[np.ndarray] = []
    for quantity, expected in expected_tensors.items():
        actual = np.stack(
            [
                actual_matrices[(quantity, component)][row_indices, phase_indices]
                for component in GEC_AXISYMMETRIC_ACTIVE_TRANSPORT_COMPONENTS
            ],
            axis=1,
        )
        equation_active_expected = expected[
            :, GEC_AXISYMMETRIC_ACTIVE_TRANSPORT_INDICES
        ]
        finite_pairs = np.isfinite(actual) & np.isfinite(equation_active_expected)
        expected_finite = equation_active_expected[
            np.isfinite(equation_active_expected)
        ]
        scale = max(
            float(np.max(np.abs(expected_finite))) if expected_finite.size else 0.0,
            1.0e-300,
        )
        error = np.abs(actual - equation_active_expected)
        allowed = (
            GEC_INDEPENDENT_TRANSPORT_RELATIVE_TOLERANCE
            * np.abs(equation_active_expected)
            + 1.0e-10 * scale
        )
        normalized = error / np.maximum(allowed, 1.0e-300)
        relative = error / np.maximum(np.abs(equation_active_expected), 1.0e-300)
        finite_value_count = int(np.count_nonzero(finite_pairs))
        all_finite_normalized.append(normalized[np.isfinite(normalized)])
        required_value_count = int(actual.size)
        maximum_normalized = finite_statistic(normalized, "maximum")
        quantity_passed = bool(
            full_state_coverage
            and finite_value_count == required_value_count
            and maximum_normalized is not None
            and maximum_normalized <= 1.0
        )
        passed = passed and quantity_passed

        nonfinite = np.flatnonzero(~finite_pairs.reshape(-1))
        if nonfinite.size:
            worst_flat = int(nonfinite[0])
            worst_score = math.inf
            worst_kind = "nonfinite"
        else:
            worst_flat = int(np.argmax(normalized.reshape(-1)))
            worst_score = float(normalized.reshape(-1)[worst_flat])
            worst_kind = "maximum_normalized_error"
        sample_index, component_index = np.unravel_index(worst_flat, actual.shape)
        row = int(row_indices[sample_index])
        phase = int(phase_indices[sample_index])
        worst = {
            "kind": worst_kind,
            "quantity": quantity,
            "component": GEC_AXISYMMETRIC_ACTIVE_TRANSPORT_COMPONENTS[component_index],
            "R_m": finite_number(coordinates[row, 0]),
            "Z_m": finite_number(coordinates[row, 1]),
            "phase": phase_labels[phase],
            "mean_energy_eV": finite_number(sample_energy[sample_index]),
            "neutral_density_m3": finite_number(sample_neutral[sample_index]),
            "Er_V_m": finite_number(sample_er[sample_index]),
            "Ez_V_m": finite_number(sample_ez[sample_index]),
            "actual": finite_number(actual[sample_index, component_index]),
            "expected": finite_number(
                equation_active_expected[sample_index, component_index]
            ),
            "absolute_error": finite_number(error[sample_index, component_index]),
            "relative_error": finite_number(relative[sample_index, component_index]),
            "normalized_error": finite_number(
                normalized[sample_index, component_index]
            ),
        }
        if worst_score > global_worst_score:
            global_worst = worst
            global_worst_score = worst_score
        quantities[quantity] = {
            "passed": quantity_passed,
            "tensor_components_checked": list(
                GEC_AXISYMMETRIC_ACTIVE_TRANSPORT_COMPONENTS
            ),
            "required_value_count": required_value_count,
            "finite_value_count": finite_value_count,
            "coverage_fraction": (
                finite_value_count / required_value_count
                if required_value_count
                else 0.0
            ),
            "maximum_relative_error": finite_statistic(relative, "maximum"),
            "p95_relative_error": finite_statistic(relative, "p95"),
            "maximum_normalized_error": maximum_normalized,
            "p95_normalized_error": finite_statistic(normalized, "p95"),
            "worst_sample": worst,
        }

    exported_component_count = len(actual_matrices)
    expected_quantities = {item[0] for item in audit_components}
    active_transport_contract = bool(
        exported_component_count == expected_component_count
        and set(expected_tensors) == expected_quantities
    )
    passed = passed and active_transport_contract
    return {
        "passed": passed,
        "status": "evaluated",
        "acceptance_gate": True,
        "acceptance_scope": "full_domain_all_exported_rf_phases",
        "support_treatment": (
            "primary_support_plus_qualified_low_endpoint_continuation"
            if qualified_low_endpoint_continuation
            else "strict_primary_support"
        ),
        "method": "independent_bundle_log_PCHIP_and_tensor_reconstruction",
        "shape_preserving_interpolator": "scipy.PchipInterpolator",
        "comsol_expected_expression_columns_used": False,
        "relative_tolerance": GEC_INDEPENDENT_TRANSPORT_RELATIVE_TOLERANCE,
        "expected_equation_active_component_count": expected_component_count,
        "exported_equation_active_component_count": exported_component_count,
        "active_external_quantities": sorted(expected_quantities),
        "active_transport_contract_complete": active_transport_contract,
        "equation_active_tensor_components": list(
            GEC_AXISYMMETRIC_ACTIVE_TRANSPORT_COMPONENTS
        ),
        "coverage": coverage,
        "maximum_normalized_error": max(
            (
                item["maximum_normalized_error"]
                for item in quantities.values()
                if item["maximum_normalized_error"] is not None
            ),
            default=None,
        ),
        "p95_normalized_error": finite_statistic(
            np.concatenate(all_finite_normalized)
            if any(values.size for values in all_finite_normalized)
            else np.asarray([], dtype=float),
            "p95",
        ),
        "worst_sample": global_worst,
        "quantities": quantities,
    }
