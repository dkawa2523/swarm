"""Audit active closure support and low-energy guard influence."""

from __future__ import annotations

import json
import math
from typing import Any

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.spatial import Delaunay, QhullError

from electron_swarm.core.constants import E_CHARGE_C

from swarm_workflow.comsol.models.gec_ccp.validation.bundle_guards import (
    GEC_ELASTIC_ENERGY_LOSS_COLUMN,
    GEC_ELASTIC_ENERGY_LOSS_TABLE,
    GEC_LOW_ENERGY_GUARD_RELATIVE_INFLUENCE_LIMIT,
)
from swarm_workflow.comsol.models.gec_ccp.closure import _uses_external_elastic_energy_loss
from swarm_workflow.comsol.models.gec_ccp.contracts import GecCcpMapping, GecCcpWorkflowError
from swarm_workflow.comsol.models.gec_ccp.data import _lookup, _read_csv


def _axisymmetric_node_volume_weights(
    radial_m: np.ndarray,
    axial_m: np.ndarray,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Build a linear-triangle 2*pi*r volume lumping on exported mesh nodes."""

    radial = np.asarray(radial_m, dtype=float)
    axial = np.asarray(axial_m, dtype=float)
    if (
        radial.ndim != 1
        or radial.shape != axial.shape
        or radial.size < 3
        or np.any(~np.isfinite(radial))
        or np.any(~np.isfinite(axial))
        or np.any(radial < 0.0)
    ):
        raise GecCcpWorkflowError(
            "axisymmetric support audit requires finite R/Z mesh nodes"
        )
    try:
        triangulation = Delaunay(np.column_stack((radial, axial)))
    except QhullError as exc:
        raise GecCcpWorkflowError(
            "axisymmetric support audit could not reconstruct the R/Z mesh"
        ) from exc
    candidate_simplices = triangulation.simplices
    axis_radius = float(np.min(radial))
    chamber_radius = float(np.max(radial))
    radial_tolerance = max(1.0e-10, 1.0e-6 * chamber_radius)
    on_axis = np.isclose(
        radial, axis_radius, atol=radial_tolerance, rtol=0.0
    )
    if int(np.count_nonzero(on_axis)) < 2:
        raise GecCcpWorkflowError(
            "axisymmetric support audit cannot identify the symmetry-axis boundary"
        )
    gap_low = float(np.min(axial[on_axis]))
    gap_high = float(np.max(axial[on_axis]))
    outer_low = float(np.min(axial))
    outer_high = float(np.max(axial))
    axial_tolerance = max(
        1.0e-10, 1.0e-6 * max(outer_high - outer_low, 1.0e-12)
    )
    outside_gap = (axial < gap_low - axial_tolerance) | (
        axial > gap_high + axial_tolerance
    )
    if np.any(outside_gap):
        geometry_kind = "stepped_gec_gas_domain"
        throat_radius = float(np.min(radial[outside_gap]))
        triangle_points = np.column_stack((radial, axial))[
            candidate_simplices
        ]
        center_r = np.mean(triangle_points[:, :, 0], axis=1)
        center_z = np.mean(triangle_points[:, :, 1], axis=1)
        in_gap = (
            (center_r <= throat_radius + radial_tolerance)
            & (center_z >= gap_low - axial_tolerance)
            & (center_z <= gap_high + axial_tolerance)
        )
        in_outer = (
            (center_r >= throat_radius - radial_tolerance)
            & (center_r <= chamber_radius + radial_tolerance)
            & (center_z >= outer_low - axial_tolerance)
            & (center_z <= outer_high + axial_tolerance)
        )
        simplices = candidate_simplices[in_gap | in_outer]
        exact_volume = math.pi * (
            (throat_radius**2 - axis_radius**2) * (gap_high - gap_low)
            + (chamber_radius**2 - throat_radius**2)
            * (outer_high - outer_low)
        )
    else:
        geometry_kind = "axisymmetric_rectangle"
        throat_radius = chamber_radius
        simplices = candidate_simplices
        exact_volume = (
            math.pi
            * (chamber_radius**2 - axis_radius**2)
            * (outer_high - outer_low)
        )
    weights = np.zeros(radial.size, dtype=float)
    for simplex in simplices:
        points = np.column_stack((radial[simplex], axial[simplex]))
        first = points[1] - points[0]
        second = points[2] - points[0]
        area = 0.5 * abs(
            first[0] * second[1] - first[1] * second[0]
        )
        volume = area * 2.0 * math.pi * float(np.mean(radial[simplex]))
        weights[simplex] += volume / 3.0
    total = float(np.sum(weights))
    if not math.isfinite(total) or total <= 0.0:
        raise GecCcpWorkflowError(
            "axisymmetric support audit reconstructed zero mesh volume"
        )
    volume_relative_error = abs(total - exact_volume) / exact_volume
    if volume_relative_error > 1.0e-10:
        raise GecCcpWorkflowError(
            "axisymmetric support-audit mesh does not reproduce the inferred "
            "GEC gas-domain volume"
        )
    return weights, {
        "method": (
            "domain_masked_Delaunay_linear_triangle_lumped_2piR_volume"
        ),
        "geometry": geometry_kind,
        "node_count": int(radial.size),
        "candidate_triangle_count": int(len(candidate_simplices)),
        "triangle_count": int(len(simplices)),
        "reconstructed_volume_m3": total,
        "inferred_exact_volume_m3": exact_volume,
        "volume_relative_error": volume_relative_error,
        "dimensions_m": {
            "axis_radius": axis_radius,
            "throat_radius": throat_radius,
            "chamber_radius": chamber_radius,
            "gap_z": [gap_low, gap_high],
            "outer_z": [outer_low, outer_high],
        },
    }


def _guard_log_pchip_values(
    rows: list[dict[str, str]],
    column: str,
    mean_energy_eV: np.ndarray,
) -> np.ndarray:
    energies, values = _lookup(rows, "mean_energy_eV", column)
    if any(value <= 0.0 for value in values):
        raise GecCcpWorkflowError(
            f"low-energy guard requires positive {column}"
        )
    query = np.clip(
        np.asarray(mean_energy_eV, dtype=float),
        float(energies[0]),
        float(energies[-1]),
    )
    return np.exp(
        PchipInterpolator(
            np.log(np.asarray(energies, dtype=float)),
            np.log(np.asarray(values, dtype=float)),
            extrapolate=False,
        )(np.log(query))
    )


def _low_energy_guard_influence_audit(
    mapping: GecCcpMapping,
    headers: list[str],
    values: np.ndarray,
    *,
    primary_support_minimum_eV: float,
    primary_support_maximum_eV: float,
) -> dict[str, Any]:
    """Bound endpoint-continuation influence with a qualified two-term guard."""

    if mapping.run.support_policy != "qualified_guard_sensitivity":
        return {
            "status": "not_requested",
            "passed": False,
            "constant_extrapolation_accepted": False,
        }
    guard_root = mapping.run.low_energy_guard_bundle
    if guard_root is None:
        raise GecCcpWorkflowError("low-energy guard bundle is missing")

    def indices(expression: str) -> list[int]:
        return [
            index
            for index, header in enumerate(headers)
            if header == expression or header.startswith(expression + " ")
        ]

    active_mobility = (
        mapping.closure.electron_transport == "swarm_mobility_einstein"
    )
    active_rate_processes = (
        set(mapping.closure.external_rate_processes)
        if mapping.closure.reaction_model == "external_rates"
        else set()
    )
    active_rates = bool(active_rate_processes)
    active_elastic = _uses_external_elastic_energy_loss(mapping.closure)
    energy_indices = indices("ptp.ebar")
    electron_indices = indices("ptp.ne")
    neutral_indices = indices("ptp.Nn")
    aligned: list[list[int]] = [
        energy_indices,
        electron_indices,
        neutral_indices,
    ]
    if active_mobility:
        radial_field_indices = indices("ptp.Er")
        axial_field_indices = indices("ptp.Ez")
        mu_rr_indices = indices("ptp.muerr")
        mu_zr_indices = indices("ptp.muezr")
        mu_rz_indices = indices("ptp.muerz")
        mu_zz_indices = indices("ptp.muezz")
        aligned.extend(
            [
                radial_field_indices,
                axial_field_indices,
                mu_rr_indices,
                mu_zr_indices,
                mu_rz_indices,
                mu_zz_indices,
            ]
        )
    phase_count = len(energy_indices)
    if phase_count == 0 or any(len(item) != phase_count for item in aligned):
        return {
            "status": "failed",
            "passed": False,
            "reason": "domain closure export lacks aligned guard-sensitivity fields",
            "constant_extrapolation_accepted": False,
        }
    radial = values[:, 0]
    axial = values[:, 1]
    node_weights, quadrature = _axisymmetric_node_volume_weights(radial, axial)
    weight = node_weights[:, None]
    energy = values[:, energy_indices]
    electron = values[:, electron_indices]
    neutral = values[:, neutral_indices]
    finite_arrays = [energy, electron, neutral]
    if active_mobility:
        radial_field = values[:, radial_field_indices]
        axial_field = values[:, axial_field_indices]
        mu_rr = values[:, mu_rr_indices]
        mu_zr = values[:, mu_zr_indices]
        mu_rz = values[:, mu_rz_indices]
        mu_zz = values[:, mu_zz_indices]
        finite_arrays.extend(
            [radial_field, axial_field, mu_rr, mu_zr, mu_rz, mu_zz]
        )
    finite = all(np.all(np.isfinite(item)) for item in finite_arrays)
    low = energy < float(primary_support_minimum_eV)
    high = energy > float(primary_support_maximum_eV)
    guard_manifest = json.loads(
        (guard_root / "manifest.json").read_text(encoding="utf-8")
    )
    guard_cross_sections = str(
        guard_manifest.get("hashes", {}).get("cross_sections_sha256", "")
    )
    primary_manifest = json.loads(
        (mapping.bundle.path / "manifest.json").read_text(encoding="utf-8")
    )
    same_cross_sections = guard_cross_sections == str(
        primary_manifest.get("hashes", {}).get("cross_sections_sha256", "")
    )
    guard_transport = _read_csv(
        guard_root / "transport_vs_mean_energy.csv"
    )
    guard_transport_energy, _ = _lookup(
        guard_transport,
        "mean_energy_eV",
        "reduced_mobility_m2_V_s_m3",
    )
    guard_minimum = float(guard_transport_energy[0])
    guard_maximum = float(guard_transport_energy[-1])
    guard_covers = bool(
        finite
        and float(np.min(energy)) >= guard_minimum
        and float(primary_support_minimum_eV) <= guard_maximum
        and not np.any(high)
    )

    def relative_l1(delta: np.ndarray, nominal: np.ndarray) -> float:
        denominator = float(np.sum(np.abs(nominal) * weight))
        numerator = float(np.sum(np.where(low, np.abs(delta), 0.0) * weight))
        return numerator / denominator if denominator > 0.0 else math.inf

    metrics: dict[str, float] = {}
    if active_mobility:
        guard_muN = _guard_log_pchip_values(
            guard_transport,
            "reduced_mobility_m2_V_s_m3",
            energy,
        )
        nominal_muN = mu_rr * neutral
        if np.any(nominal_muN <= 0.0):
            return {
                "status": "failed",
                "passed": False,
                "reason": "nominal COMSOL mobility is nonpositive",
                "constant_extrapolation_accepted": False,
            }
        mobility_ratio = guard_muN / nominal_muN
        electron_flux_r = electron * (
            mu_rr * radial_field + mu_zr * axial_field
        )
        electron_flux_z = electron * (
            mu_rz * radial_field + mu_zz * axial_field
        )
        drift_flux = np.hypot(electron_flux_r, electron_flux_z)
        electric_field = np.hypot(radial_field, axial_field)
        metrics["mobility_drift_flux_relative_L1"] = relative_l1(
            np.abs(mobility_ratio - 1.0) * drift_flux,
            drift_flux,
        )
        metrics["mobility_drift_power_proxy_relative_L1"] = relative_l1(
            np.abs(mobility_ratio - 1.0) * drift_flux * electric_field,
            drift_flux * electric_field,
        )

    if active_rates:
        guard_rates = _read_csv(guard_root / "rates_vs_mean_energy.csv")
        for reaction_index, process_type in (
            (2, "excitation"),
            (3, "ionization"),
        ):
            if process_type not in active_rate_processes:
                continue
            nominal_indices = indices(f"ptp.kf_{reaction_index}")
            if len(nominal_indices) != phase_count:
                return {
                    "status": "failed",
                    "passed": False,
                    "reason": (
                        f"domain export lacks aligned {process_type} rate"
                    ),
                    "constant_extrapolation_accepted": False,
                }
            selected = [
                row
                for row in guard_rates
                if row.get("process_type") == process_type
            ]
            guard_rate = _guard_log_pchip_values(
                selected,
                "rate_coefficient_m3_s",
                energy,
            )
            nominal_rate = values[:, nominal_indices] / 6.02214076e23
            nominal_source = electron * neutral * nominal_rate
            metrics[f"{process_type}_source_relative_L1"] = relative_l1(
                electron * neutral * np.abs(guard_rate - nominal_rate),
                nominal_source,
            )

    if active_elastic:
        elastic_expression = (
            "-ptp.ne*ptp.n_wAr*exp(sw_logKel(sw_logeps_el))*1[eV*m^3/s]"
        )
        elastic_indices = indices(elastic_expression)
        neutral_argon_indices = indices("ptp.n_wAr")
        if (
            len(elastic_indices) != phase_count
            or len(neutral_argon_indices) != phase_count
        ):
            return {
                "status": "failed",
                "passed": False,
                "reason": "domain export lacks aligned elastic-loss fields",
                "constant_extrapolation_accepted": False,
            }
        nominal_elastic_power = -values[:, elastic_indices]
        argon_density = values[:, neutral_argon_indices]
        guard_elastic = _guard_log_pchip_values(
            _read_csv(guard_root / GEC_ELASTIC_ENERGY_LOSS_TABLE),
            GEC_ELASTIC_ENERGY_LOSS_COLUMN,
            energy,
        )
        guard_elastic_power = (
            electron * argon_density * guard_elastic * E_CHARGE_C
        )
        metrics["elastic_power_relative_L1"] = relative_l1(
            guard_elastic_power - nominal_elastic_power,
            nominal_elastic_power,
        )

    if not metrics:
        return {
            "status": "failed",
            "passed": False,
            "reason": "no active external closure quantity requires a guard",
            "constant_extrapolation_accepted": False,
        }
    maximum = max(metrics.values())
    limit = GEC_LOW_ENERGY_GUARD_RELATIVE_INFLUENCE_LIMIT
    passed = bool(
        finite
        and same_cross_sections
        and guard_covers
        and np.any(low)
        and maximum <= limit
    )
    return {
        "status": "passed" if passed else "failed",
        "passed": passed,
        "method": "qualified_two_term_guard_frozen_field_axisymmetric_L1",
        "role": "application_specific_endpoint_continuation_influence_bound",
        "primary_source": "monte_carlo",
        "guard_source": "two_term",
        "active_external_quantities": [
            name
            for name, active in (
                ("mobility", active_mobility),
                ("reaction_rates", active_rates),
                ("elastic_energy_loss", active_elastic),
            )
            if active
        ],
        "guard_bundle": str(guard_root),
        "same_cross_sections": same_cross_sections,
        "guard_support_mean_energy_eV": [guard_minimum, guard_maximum],
        "primary_support_mean_energy_eV": [
            float(primary_support_minimum_eV),
            float(primary_support_maximum_eV),
        ],
        "operating_range_mean_energy_eV": [
            float(np.min(energy)),
            float(np.max(energy)),
        ],
        "guard_covers_all_low_excursions": guard_covers,
        "high_energy_excursion_count": int(np.count_nonzero(high)),
        "low_energy_excursion_count": int(np.count_nonzero(low)),
        "exported_state_count": int(energy.size),
        "relative_influence_limit": limit,
        "maximum_relative_influence": maximum,
        "relative_influence": metrics,
        "quadrature": quadrature,
        "constant_extrapolation_accepted": passed,
        "scope_limit": (
            "frozen_converged_field influence bound for this GEC operating "
            "point; it does not extend the primary MC coefficient support"
        ),
    }


def _operating_mean_energy_support_audit(
    mapping: GecCcpMapping,
    headers: list[str],
    values: np.ndarray,
    closure_support: dict[str, Any],
) -> dict[str, Any]:
    """Qualify every active external table over the converged RF state."""

    energy_indices = [
        index
        for index, header in enumerate(headers)
        if header == "ptp.ebar" or header.startswith("ptp.ebar ")
    ]
    if not energy_indices:
        return {
            "status": "unavailable",
            "passed": False,
            "reason": "domain phase export lacks ptp.ebar",
            **closure_support,
        }
    if not closure_support.get("passed"):
        return {"status": "unavailable", **closure_support}

    operating = values[:, energy_indices]
    finite = bool(np.all(np.isfinite(operating)))
    operating_min = float(np.min(operating)) if finite else math.nan
    operating_max = float(np.max(operating)) if finite else math.nan
    common_min, common_max = closure_support[
        "common_intersection_mean_energy_eV"
    ]
    accepted_min, accepted_max = closure_support[
        "accepted_interior_mean_energy_eV"
    ]
    strictly_inside = bool(
        finite
        and common_min < operating_min
        and operating_max < common_max
    )
    margin_passed = bool(
        finite
        and accepted_min <= operating_min
        and operating_max <= accepted_max
    )
    strict_primary_support_passed = strictly_inside and margin_passed
    low_energy_guard_influence: dict[str, Any] = {
        "status": "not_requested",
        "passed": False,
        "constant_extrapolation_accepted": False,
    }
    if (
        not strict_primary_support_passed
        and mapping.run.support_policy == "qualified_guard_sensitivity"
    ):
        low_energy_guard_influence = _low_energy_guard_influence_audit(
            mapping,
            headers,
            values,
            primary_support_minimum_eV=float(common_min),
            primary_support_maximum_eV=float(common_max),
        )
    guard_support_passed = bool(
        low_energy_guard_influence.get("passed", False)
    )
    effective_support_passed = bool(
        strict_primary_support_passed or guard_support_passed
    )
    return {
        **closure_support,
        "status": "evaluated",
        "passed": effective_support_passed,
        "phase_local_mean_energy_eV": [operating_min, operating_max],
        "strictly_inside_common_support": strictly_inside,
        "endpoint_margin_passed": margin_passed,
        "strict_primary_support_passed": strict_primary_support_passed,
        "support_policy": mapping.run.support_policy,
        "support_acceptance_basis": (
            "strict_primary_support"
            if strict_primary_support_passed
            else (
                "qualified_low_energy_guard_sensitivity"
                if guard_support_passed
                else "rejected"
            )
        ),
        "constant_extrapolation_accepted": bool(
            not strict_primary_support_passed and guard_support_passed
        ),
        "low_energy_guard_influence": low_energy_guard_influence,
        "exported_domain_nodes": int(values.shape[0]),
        "exported_phase_samples_per_node": len(energy_indices),
    }
