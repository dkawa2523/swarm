"""Independent acceptance analysis for COMSOL-sampled Function-EEDF values."""

from __future__ import annotations

import hashlib
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from ..function_eedf import FunctionEedfError, piecewise_linear_weighted_moments
from .contracts import (
    _PROJECTION_FIDELITY_KINDS,
    _PROJECTION_MAXIMUM_RELATIVE_ERROR_LIMIT,
    _PROJECTION_NORMALIZED_RMSE_LIMIT,
    _PROJECTION_SIGNIFICANCE_FRACTION,
    _QUERY_COLUMNS,
    _STRICT_BINDING_KINDS,
    _VALUE_COLUMNS,
    ComsolEedfAuditError,
    ComsolEedfAuditPlan,
)
from .io import (
    _finite_float,
    _read_contract,
    _read_cross_sections,
    _read_rows,
    _validate_audit_inputs,
)


def analyze_comsol_eedf_audit(
    plan: ComsolEedfAuditPlan,
    *,
    cross_section_csv: str | Path | None = None,
    point_rtol: float = 1.0e-8,
    point_atol: float = 1.0e-13,
    moment_tolerance: float = 1.0e-7,
) -> dict[str, Any]:
    """Fail closed on saved settings and independently sampled COMSOL values."""

    _validate_audit_inputs(plan)
    contract = _read_contract(plan.contract_path)
    contract_expected = {
        "function_names": plan.function_tag,
        "source": "file",
        "struct": plan.import_contract.structure,
        "interp": "linear",
        "extrap": "const",
        "fununit": "1",
        "argunit": "eV|eV",
        "problem_count": "0",
        "eedf_selection": plan.function_tag,
        "audit_table_sha256": plan.table_sha256,
        "audit_queries_sha256": plan.query_sha256,
    }
    contract_mismatches = {
        key: {"expected": expected, "actual": contract.get(key)}
        for key, expected in contract_expected.items()
        if contract.get(key) != expected
    }
    if not contract.get("imported_name"):
        contract_mismatches["imported_name"] = {
            "expected": "nonempty imported resource",
            "actual": contract.get("imported_name"),
        }

    queries = _read_rows(plan.query_path, _QUERY_COLUMNS)
    values = _read_rows(plan.values_path, _VALUE_COLUMNS)
    expected_ids = [row["point_id"] for row in queries]
    actual_by_id = {row["point_id"]: row for row in values}
    if (
        len(expected_ids) != plan.point_count
        or len(set(expected_ids)) != plan.point_count
        or len(actual_by_id) != len(values)
        or set(actual_by_id) != set(expected_ids)
    ):
        raise ComsolEedfAuditError(
            "COMSOL EEDF audit values do not cover the query plan"
        )

    joined: list[dict[str, Any]] = []
    minimum_value = math.inf
    for query in queries:
        if any(
            actual_by_id[query["point_id"]][column] != query[column]
            for column in _QUERY_COLUMNS
        ):
            raise ComsolEedfAuditError(
                "COMSOL audit values contain mismatched probe coordinates"
            )
        actual = _finite_float(actual_by_id[query["point_id"]]["comsol_value"])
        expected_a = _finite_float(query["expected_a"])
        expected_b = _finite_float(query["expected_b"])
        scale = max(abs(expected_a), abs(expected_b), point_atol)
        error = min(abs(actual - expected_a), abs(actual - expected_b))
        passed = error <= point_atol + point_rtol * scale and actual >= -point_atol
        minimum_value = min(minimum_value, actual)
        joined.append(
            {
                **query,
                "comsol_value": actual,
                "candidate_error": error,
                "candidate_relative_error": error / scale,
                "strict_candidate_match": passed,
            }
        )

    moment_rows: dict[str, list[dict[str, Any]]] = {}
    for row in joined:
        if row["kind"] == "moment":
            moment_rows.setdefault(row["group"], []).append(row)
    moment_results: list[dict[str, Any]] = []
    observed_rows: dict[float, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for group, rows in sorted(moment_rows.items()):
        rows.sort(key=lambda item: _finite_float(item["electron_energy_eV"]))
        energy = np.asarray([_finite_float(row["electron_energy_eV"]) for row in rows])
        actual = np.asarray([float(row["comsol_value"]) for row in rows])
        expected = np.asarray([_finite_float(row["expected_a"]) for row in rows])
        requested_mean = _finite_float(rows[0]["mean_energy_eV"])
        try:
            normalization, first = piecewise_linear_weighted_moments(energy, actual)
        except FunctionEedfError as exc:
            raise ComsolEedfAuditError(
                "invalid COMSOL-sampled EEDF moment row"
            ) from exc
        mean = first / max(normalization, 1.0e-300)
        moment_results.append(
            {
                "group": group,
                "requested_mean_energy_eV": requested_mean,
                "normalization": normalization,
                "normalization_error": abs(normalization - 1.0),
                "reconstructed_mean_energy_eV": mean,
                "mean_energy_relative_error": abs(mean - requested_mean)
                / requested_mean,
            }
        )
        observed_rows[requested_mean] = (energy, actual, expected)

    moment_error = max(
        (
            max(row["normalization_error"], row["mean_energy_relative_error"])
            for row in moment_results
        ),
        default=math.inf,
    )
    derivative = _derivative_audit(joined)
    rates = _rate_audit(observed_rows, cross_section_csv)
    strict_rows = [row for row in joined if _is_strict_binding_row(row)]
    projection_rows = [row for row in joined if _is_projection_fidelity_row(row)]
    strict_binding_passed = bool(
        strict_rows and all(bool(row["strict_candidate_match"]) for row in strict_rows)
    )
    projection_fidelity = _projection_fidelity_audit(
        projection_rows,
        absolute_floor=point_atol,
    )
    passed = bool(
        not contract_mismatches
        and strict_binding_passed
        and projection_fidelity["passed"]
        and minimum_value >= -point_atol
        and moment_results
        and moment_error <= moment_tolerance
        and derivative["passed"]
        and rates["passed"]
    )
    return {
        "passed": passed,
        "method": "COMSOL parameter-engine evaluation of the saved component Function",
        "provenance": {
            "table_sha256": plan.table_sha256,
            "query_sha256": plan.query_sha256,
            "model_path": str(plan.model_path),
            "model_sha256": hashlib.sha256(plan.model_path.read_bytes()).hexdigest(),
            "values_sha256": hashlib.sha256(plan.values_path.read_bytes()).hexdigest(),
            "contract_sha256": hashlib.sha256(
                plan.contract_path.read_bytes()
            ).hexdigest(),
        },
        "active_table_binding": {
            "passed": strict_binding_passed,
            "table_path": str(plan.table_path),
            "table_sha256": plan.table_sha256,
            "comparison": (
                "strict COMSOL-evaluated saved Function values at table nodes, "
                "one-coordinate midpoints, derivative probes, and constant "
                "extrapolation probes"
            ),
            "strict_kinds": [
                "anchor",
                "energy_midpoint",
                "derivative",
                "outside",
            ],
            "samples": len(strict_rows),
            "failed_samples": sum(
                not bool(row["strict_candidate_match"]) for row in strict_rows
            ),
            "maximum_candidate_relative_error": max(
                (float(row["candidate_relative_error"]) for row in strict_rows),
                default=math.inf,
            ),
        },
        "native_projection_fidelity": projection_fidelity,
        "saved_function_contract": {
            "passed": not contract_mismatches,
            "mismatches": contract_mismatches,
            "settings": contract,
            "expected": contract_expected,
        },
        "point_values": {
            "passed": strict_binding_passed and projection_fidelity["passed"],
            "samples": len(joined),
            "strict_samples": len(strict_rows),
            "projection_fidelity_samples": len(projection_rows),
            "minimum_value": minimum_value,
            "kinds": sorted({str(row["kind"]) for row in joined}),
        },
        "moments": {
            "passed": bool(moment_results and moment_error <= moment_tolerance),
            "tolerance": moment_tolerance,
            "maximum_error": moment_error,
            "rows": moment_results,
        },
        "derivatives": derivative,
        "rates": rates,
    }


def _is_strict_binding_row(row: Mapping[str, Any]) -> bool:
    return str(row["kind"]) in _STRICT_BINDING_KINDS


def _is_projection_fidelity_row(row: Mapping[str, Any]) -> bool:
    kind = str(row["kind"])
    return kind in _PROJECTION_FIDELITY_KINDS or kind == "moment"


def _projection_fidelity_audit(
    rows: list[dict[str, Any]],
    *,
    absolute_floor: float,
) -> dict[str, Any]:
    errors = np.asarray([float(row["candidate_error"]) for row in rows])
    references = np.asarray(
        [
            max(
                abs(float(row["comsol_value"])),
                abs(_finite_float(row["expected_a"])),
                abs(_finite_float(row["expected_b"])),
            )
            for row in rows
        ],
        dtype=float,
    )
    peak = float(np.max(references)) if references.size else 0.0
    significant = references >= peak * _PROJECTION_SIGNIFICANCE_FRACTION
    significant_relative = errors[significant] / np.maximum(
        references[significant], absolute_floor
    )
    maximum_relative = (
        float(np.max(significant_relative)) if significant_relative.size else math.inf
    )
    normalized_rmse = (
        float(
            np.sqrt(np.mean(np.square(errors)))
            / max(np.sqrt(np.mean(np.square(references))), absolute_floor)
        )
        if references.size
        else math.inf
    )
    passed = bool(
        rows
        and np.all(np.isfinite(errors))
        and np.all(np.isfinite(references))
        and maximum_relative <= _PROJECTION_MAXIMUM_RELATIVE_ERROR_LIMIT
        and normalized_rmse <= _PROJECTION_NORMALIZED_RMSE_LIMIT
    )
    return {
        "passed": passed,
        "interpretation": (
            "fidelity of COMSOL's structured tensor-product linear grid at "
            "dense moment rows and two-coordinate interior probes"
        ),
        "kinds": ["cell_center", "moment_anchor", "moment_mid"],
        "samples": len(rows),
        "significant_samples": int(np.count_nonzero(significant)),
        "significance_fraction_of_peak": _PROJECTION_SIGNIFICANCE_FRACTION,
        "maximum_relative_error_significant": maximum_relative,
        "maximum_relative_error_limit": (_PROJECTION_MAXIMUM_RELATIVE_ERROR_LIMIT),
        "all_point_normalized_rmse": normalized_rmse,
        "normalized_rmse_limit": _PROJECTION_NORMALIZED_RMSE_LIMIT,
    }


def _derivative_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    groups: dict[str, dict[str, dict[str, Any]]] = {}
    for row in rows:
        if row["kind"] == "derivative":
            groups.setdefault(row["group"], {})[row["position"]] = row
    maximum_slope_error = 0.0
    maximum_jump = 0.0
    for group in groups.values():
        if set(group) != {"left", "center", "right"}:
            return {"passed": False, "reason": "incomplete finite-difference group"}
        left, center, right = group["left"], group["center"], group["right"]
        x_left = _finite_float(left["mean_energy_eV"])
        x_center = _finite_float(center["mean_energy_eV"])
        x_right = _finite_float(right["mean_energy_eV"])
        actual_left = (float(center["comsol_value"]) - float(left["comsol_value"])) / (
            x_center - x_left
        )
        actual_right = (
            float(right["comsol_value"]) - float(center["comsol_value"])
        ) / (x_right - x_center)
        expected_left = (
            _finite_float(center["expected_a"]) - _finite_float(left["expected_a"])
        ) / (x_center - x_left)
        expected_right = (
            _finite_float(right["expected_a"]) - _finite_float(center["expected_a"])
        ) / (x_right - x_center)
        scale = max(abs(expected_left), abs(expected_right), 1.0e-13)
        maximum_slope_error = max(
            maximum_slope_error,
            abs(actual_left - expected_left) / scale,
            abs(actual_right - expected_right) / scale,
        )
        maximum_jump = max(maximum_jump, abs(actual_right - actual_left))
    return {
        "passed": bool(
            groups
            and math.isfinite(maximum_slope_error)
            and maximum_slope_error <= 1.0e-4
        ),
        "continuity_contract": "C0; finite one-sided slopes may jump at mean-energy anchors",
        "groups": len(groups),
        "maximum_relative_one_sided_slope_error": maximum_slope_error,
        "maximum_observed_slope_jump": maximum_jump,
    }


def _rate_audit(
    observed_rows: Mapping[float, tuple[np.ndarray, np.ndarray, np.ndarray]],
    cross_section_csv: str | Path | None,
) -> dict[str, Any]:
    if cross_section_csv is None:
        return {"passed": True, "status": "not_requested", "processes": {}}
    cross_sections = _read_cross_sections(Path(cross_section_csv))
    processes: dict[str, Any] = {}
    for process, (cross_energy, sigma) in cross_sections.items():
        actual_rates: list[float] = []
        expected_rates: list[float] = []
        for energy, actual, expected in observed_rows.values():
            actual_rates.append(_integrate_rate(energy, actual, cross_energy, sigma))
            expected_rates.append(
                _integrate_rate(energy, expected, cross_energy, sigma)
            )
        actual_array = np.asarray(actual_rates)
        expected_array = np.asarray(expected_rates)
        difference = actual_array - expected_array
        normalized_rmse = float(
            np.sqrt(np.mean(np.square(difference)))
            / max(np.sqrt(np.mean(np.square(expected_array))), 1.0e-300)
        )
        processes[process] = {
            "samples": len(actual_rates),
            "normalized_rmse": normalized_rmse,
            "passed": normalized_rmse <= 1.0e-6,
        }
    return {
        "passed": bool(
            processes and all(item["passed"] for item in processes.values())
        ),
        "status": "evaluated",
        "processes": processes,
    }


def _integrate_rate(
    energy: np.ndarray,
    f0_nodes: np.ndarray,
    cross_energy: np.ndarray,
    sigma_nodes: np.ndarray,
) -> float:
    knots = np.unique(
        np.concatenate(
            (
                energy,
                cross_energy[(cross_energy > energy[0]) & (cross_energy < energy[-1])],
            )
        )
    )
    nodes, weights = np.polynomial.legendre.leggauss(4)
    left = knots[:-1, None]
    right = knots[1:, None]
    half_width = 0.5 * (right - left)
    sample_energy = 0.5 * (right + left) + half_width * nodes
    f0 = np.interp(
        sample_energy, energy, f0_nodes, left=f0_nodes[0], right=f0_nodes[-1]
    )
    sigma = np.interp(
        sample_energy, cross_energy, sigma_nodes, left=0.0, right=sigma_nodes[-1]
    )
    speed = math.sqrt(2.0 * 1.602176634e-19 / 9.1093837139e-31)
    return float(np.sum(half_width * weights * speed * sigma * sample_energy * f0))
