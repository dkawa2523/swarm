"""Validate Function-EEDF artifacts, cross sections, and reconstructed rates."""

from __future__ import annotations

import csv
import hashlib
import math
from pathlib import Path
import re
from typing import Any

import numpy as np

from swarm_workflow._paths import discover_repo_root
from swarm_workflow.comsol.input.function_eedf import (
    ComsolFunctionEedfGrid,
    evaluate_comsol_function_eedf_grid,
    FunctionEedfError,
    piecewise_linear_weighted_moments,
    read_comsol_function_eedf_grid,
)
from swarm_workflow.comsol.models.gec_ccp.contracts import (
    GecCcpMapping,
    GecCcpWorkflowError,
    GecReactionSpec,
)
from swarm_workflow.comsol.models.gec_ccp.data import (
    _float_or_none,
    _read_csv,
    _required_float,
)


FUNCTION_EEDF_RATE_SIGNIFICANCE_FRACTION = 1.0e-5
FUNCTION_EEDF_RATE_RELATIVE_ERROR_LIMIT = 0.10
FUNCTION_EEDF_RATE_P95_RELATIVE_ERROR_LIMIT = 0.05
FUNCTION_EEDF_RATE_NORMALIZED_RMSE_LIMIT = 0.02
UPSTREAM_EEDF_RATE_AUDIT_TABLE = "upstream_eedf_rate_audit.csv"


def _reaction_cross_sections_from_model_xml(
    model_xml: str,
    reactions: tuple[GecReactionSpec, ...],
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    tables: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for reaction in reactions:
        feature = re.search(
            rf'<PhysicsFeature\b[^>]*op="ElectronImpactReaction"[^>]*tag="'
            rf'{re.escape(reaction.feature)}"[^>]*>(.*?)</PhysicsFeature>',
            model_xml,
            flags=re.DOTALL,
        )
        if feature is None:
            raise GecCcpWorkflowError(
                f"saved MPH lacks reaction feature {reaction.feature}"
            )
        block = feature.group(1)
        arrays: dict[str, np.ndarray] = {}
        for name in ("xdata", "ydata"):
            match = re.search(
                rf'<param\b[^>]*param="{name}"[^>]*value="([^"]+)"',
                block,
            )
            if match is None:
                raise GecCcpWorkflowError(
                    f"saved MPH reaction {reaction.feature} lacks {name}"
                )
            arrays[name] = np.asarray(
                [float(value) for value in re.findall(r"'([^']+)'", match.group(1))],
                dtype=float,
            )
        if arrays["xdata"].size < 2 or arrays["xdata"].shape != arrays["ydata"].shape:
            raise GecCcpWorkflowError(
                f"saved MPH reaction {reaction.feature} has invalid cross sections"
            )
        tables[reaction.process_type] = (arrays["xdata"], arrays["ydata"])
    return tables


def _audit_gec_argon_cross_section_identity(
    model_tables: dict[str, tuple[np.ndarray, np.ndarray]],
    reactions: tuple[GecReactionSpec, ...],
) -> dict[str, Any]:
    reference = (
        discover_repo_root(__file__)
        / "examples"
        / "cross_sections"
        / "argon_application_library.csv"
    )
    if not reference.exists():
        return {
            "passed": False,
            "reference": str(reference),
            "reason": "reference cross-section table is missing",
        }
    rows = _read_csv(reference)
    processes: dict[str, Any] = {}
    for reaction in reactions:
        selected = [row for row in rows if row.get("type") == reaction.process_type]
        reference_energy = np.asarray(
            [_required_float(row, "energy_eV") for row in selected], dtype=float
        )
        reference_sigma = np.asarray(
            [_required_float(row, "cross_section_m2") for row in selected],
            dtype=float,
        )
        model_energy, model_sigma = model_tables[reaction.process_type]
        energy_match = bool(
            model_energy.shape == reference_energy.shape
            and np.allclose(model_energy, reference_energy, rtol=0.0, atol=0.0)
        )
        sigma_match = bool(
            model_sigma.shape == reference_sigma.shape
            and np.allclose(model_sigma, reference_sigma, rtol=1.0e-14, atol=0.0)
        )
        processes[reaction.process_type] = {
            "mph_points": int(model_energy.size),
            "reference_points": int(reference_energy.size),
            "energy_pointwise_identical": energy_match,
            "cross_section_pointwise_identical": sigma_match,
            "passed": energy_match and sigma_match,
        }
    return {
        "passed": all(item["passed"] for item in processes.values()),
        "reference": str(reference),
        "comparison": "pointwise energy and cross-section values",
        "processes": processes,
    }


def _rate_significance(
    first: np.ndarray,
    second: np.ndarray,
) -> tuple[float, float, np.ndarray]:
    """Return a shared absolute-significance floor for pointwise rate QA."""
    process_peak = max(
        float(np.max(np.abs(first))),
        float(np.max(np.abs(second))),
        1.0e-300,
    )
    significance_floor = process_peak * FUNCTION_EEDF_RATE_SIGNIFICANCE_FRACTION
    significant = np.maximum(np.abs(first), np.abs(second)) >= significance_floor
    return process_peak, significance_floor, significant


def _rate_significance_metadata(
    *,
    process_peak: float,
    significance_floor: float,
    significant: np.ndarray,
) -> dict[str, Any]:
    return {
        "process_peak_m3_s": process_peak,
        "significance_fraction_of_process_peak": (
            FUNCTION_EEDF_RATE_SIGNIFICANCE_FRACTION
        ),
        "significance_floor_m3_s": significance_floor,
        "significant_samples": int(np.count_nonzero(significant)),
        "record_only_samples": int(significant.size - np.count_nonzero(significant)),
        "below_significance_floor_policy": ("record_only_for_pointwise_relative_error"),
    }


def _audit_function_eedf_source_rates(
    mapping: GecCcpMapping,
    function_grid: ComsolFunctionEedfGrid,
    cross_sections: dict[str, tuple[np.ndarray, np.ndarray]],
    reactions: tuple[GecReactionSpec, ...],
    *,
    rate_table_name: str = "rates_vs_mean_energy.csv",
    audit_csv_name: str = "function_eedf_rate_audit.csv",
) -> dict[str, Any]:
    rate_rows = _read_csv(mapping.bundle.path / rate_table_name)
    solver_source = mapping.bundle.expected_source
    audit_rows: list[dict[str, Any]] = []
    process_results: dict[str, Any] = {}
    for reaction in reactions:
        selected = [
            row for row in rate_rows if row.get("process_type") == reaction.process_type
        ]
        if not selected:
            raise GecCcpWorkflowError(
                f"missing {solver_source} rate evidence for {reaction.process_type}"
            )
        cross_energy, cross_sigma = cross_sections[reaction.process_type]
        derived = np.asarray(
            [
                _integrate_native_function_eedf_rate(
                    function_grid,
                    _required_float(row, "mean_energy_eV"),
                    cross_energy,
                    cross_sigma,
                )
                for row in selected
            ],
            dtype=float,
        )
        source = np.asarray(
            [_required_float(row, "rate_coefficient_m3_s") for row in selected],
            dtype=float,
        )
        process_peak, significance_floor, relevant = _rate_significance(derived, source)
        difference = derived - source
        normalized_rmse = float(
            np.sqrt(np.mean(np.square(difference)))
            / max(np.sqrt(np.mean(np.square(source))), 1.0e-300)
        )
        relative = np.abs(difference[relevant]) / np.maximum(
            np.maximum(derived[relevant], source[relevant]), 1.0e-300
        )
        maximum_relative = float(np.max(relative)) if relative.size else 0.0
        p95_relative = float(np.quantile(relative, 0.95)) if relative.size else 0.0
        process_passed = (
            normalized_rmse <= FUNCTION_EEDF_RATE_NORMALIZED_RMSE_LIMIT
            and maximum_relative <= FUNCTION_EEDF_RATE_RELATIVE_ERROR_LIMIT
            and p95_relative <= FUNCTION_EEDF_RATE_P95_RELATIVE_ERROR_LIMIT
        )
        row_passes: list[bool] = []
        for index, row in enumerate(selected):
            row_relevant = bool(relevant[index])
            row_pass = (
                not row_relevant
                or abs(derived[index] - source[index])
                / max(derived[index], source[index], 1.0e-300)
                <= FUNCTION_EEDF_RATE_RELATIVE_ERROR_LIMIT
            )
            row_passes.append(row_pass)
            audit_rows.append(
                {
                    "process_type": reaction.process_type,
                    "E_over_N_Td": _float_or_none(row.get("E_over_N_Td")),
                    "mean_energy_eV": _required_float(row, "mean_energy_eV"),
                    "f0_integrated_rate_m3_s": derived[index],
                    "solver_rate_m3_s": source[index],
                    "solver_rate_is_zero": bool(source[index] == 0.0),
                    "relevant": row_relevant,
                    "relative_error_if_relevant": (
                        abs(derived[index] - source[index])
                        / max(derived[index], source[index], 1.0e-300)
                        if row_relevant
                        else None
                    ),
                    "criterion": (
                        f"{solver_source}_relative_10_percent"
                        if row_relevant
                        else "below_1e-5_process_peak_record_only"
                    ),
                    "passed": row_pass,
                }
            )
        process_results[reaction.process_type] = {
            "passed": process_passed and all(row_passes),
            "anchors": int(source.size),
            "relevant_anchors": int(np.count_nonzero(relevant)),
            "zero_rate_anchors": int(np.count_nonzero(source == 0.0)),
            "relevant_zero_rate_anchors": int(
                np.count_nonzero((source == 0.0) & relevant)
            ),
            "normalized_rmse": normalized_rmse,
            "normalized_rmse_limit": (FUNCTION_EEDF_RATE_NORMALIZED_RMSE_LIMIT),
            "normalized_rmse_scope": "all_samples",
            "maximum_relative_error_significant": maximum_relative,
            "maximum_relative_error_limit": (FUNCTION_EEDF_RATE_RELATIVE_ERROR_LIMIT),
            "p95_relative_error_significant": p95_relative,
            "p95_relative_error_limit": (FUNCTION_EEDF_RATE_P95_RELATIVE_ERROR_LIMIT),
            **_rate_significance_metadata(
                process_peak=process_peak,
                significance_floor=significance_floor,
                significant=relevant,
            ),
        }

    audit_path = mapping.results.output_directory / audit_csv_name
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    columns = (
        "process_type",
        "E_over_N_Td",
        "mean_energy_eV",
        "f0_integrated_rate_m3_s",
        "solver_rate_m3_s",
        "solver_rate_is_zero",
        "relevant",
        "relative_error_if_relevant",
        "criterion",
        "passed",
    )
    with audit_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(audit_rows)
    return {
        "passed": all(item["passed"] for item in process_results.values()),
        "source": f"{solver_source}_rates",
        "artificial_rate_floor": False,
        "rate_table": rate_table_name,
        "rate_table_sha256": hashlib.sha256(
            (mapping.bundle.path / rate_table_name).read_bytes()
        ).hexdigest(),
        "audit_csv": str(audit_path),
        "audit_csv_sha256": hashlib.sha256(audit_path.read_bytes()).hexdigest(),
        "processes": process_results,
    }


def _integrate_native_function_eedf_rate(
    function_grid: ComsolFunctionEedfGrid,
    mean_energy_eV: float,
    cross_section_energy_eV: np.ndarray,
    cross_section_m2: np.ndarray,
) -> float:
    energy = function_grid.electron_energies_eV
    f0_nodes = _native_function_eedf_row(function_grid, mean_energy_eV)
    knots = np.unique(
        np.concatenate(
            (
                energy,
                cross_section_energy_eV[
                    (cross_section_energy_eV > energy[0])
                    & (cross_section_energy_eV < energy[-1])
                ],
            )
        )
    )
    nodes, weights = np.polynomial.legendre.leggauss(4)
    left = knots[:-1, None]
    right = knots[1:, None]
    half_width = 0.5 * (right - left)
    sample_energy = 0.5 * (right + left) + half_width * nodes
    f0 = np.interp(
        sample_energy,
        energy,
        f0_nodes,
        left=float(f0_nodes[0]),
        right=float(f0_nodes[-1]),
    )
    sigma = np.interp(
        sample_energy,
        cross_section_energy_eV,
        cross_section_m2,
        left=0.0,
        right=float(cross_section_m2[-1]),
    )
    speed_probability_factor = np.sqrt(2.0 * 1.602176634e-19 / 9.1093837139e-31)
    integrand = speed_probability_factor * sigma * sample_energy * f0
    return float(np.sum(half_width * weights * integrand))


def _read_function_eedf(path: Path) -> ComsolFunctionEedfGrid:
    try:
        return read_comsol_function_eedf_grid(path)
    except FunctionEedfError as exc:
        raise GecCcpWorkflowError(f"invalid Function-EEDF: {exc}") from exc


def _audit_upstream_eedf_artifact(
    path: Path,
    entry: dict[str, Any],
    *,
    expected_sha256: str,
) -> dict[str, Any]:
    """Re-read inactive EEDF evidence instead of trusting manifest claims."""

    expected_columns = [
        "electron_energy_eV",
        "mean_energy_eV",
        "eepf_eV_m32",
    ]
    if (
        entry.get("artifact_role") != "canonical_comsol_function_eedf_input"
        or entry.get("canonical_comsol_input") is not True
        or entry.get("format") != "csv"
        or entry.get("representation") != "physical_2d_adaptive_moment_rate_projected"
        or entry.get("argument_order") != ["electron_energy_eV", "mean_energy_eV"]
        or entry.get("columns") != expected_columns
        or entry.get("grid_axis_order") != ["electron_energy_eV", "mean_energy_eV"]
        or entry.get("row_order") != "mean_energy_major_then_electron_energy"
    ):
        raise GecCcpWorkflowError(
            "upstream EEDF evidence lacks the canonical projected-grid contract"
        )
    actual_sha256 = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual_sha256 != expected_sha256:
        raise GecCcpWorkflowError("upstream EEDF SHA-256 mismatch")
    grid = _read_function_eedf(path)
    shape = [
        len(grid.mean_energies_eV),
        len(grid.electron_energies_eV),
    ]
    if (
        entry.get("grid_shape") != shape
        or entry.get("mean_energy_grid_points") != shape[0]
        or entry.get("energy_grid_points") != shape[1]
        or entry.get("electron_energy_range_eV")
        != [
            float(grid.electron_energies_eV[0]),
            float(grid.electron_energies_eV[-1]),
        ]
        or entry.get("mean_energy_range_eV")
        != [
            float(grid.mean_energies_eV[0]),
            float(grid.mean_energies_eV[-1]),
        ]
    ):
        raise GecCcpWorkflowError(
            "upstream EEDF manifest dimensions disagree with the structured grid"
        )
    audit = _native_function_eedf_moment_audit(grid)
    reported = {
        "normalization_error_max": _float_or_none(
            entry.get("projected_normalization_error_max")
        ),
        "mean_energy_relative_error_max": _float_or_none(
            entry.get("projected_mean_energy_relative_error_max")
        ),
        "nonnegative_minimum": _float_or_none(
            entry.get("projected_nonnegative_minimum")
        ),
    }
    if (
        reported["normalization_error_max"] is None
        or reported["mean_energy_relative_error_max"] is None
        or reported["nonnegative_minimum"] is None
        or audit["normalization_error_max"] > 1.0e-8
        or audit["mean_energy_relative_error_max"] > 1.0e-8
        or audit["nonnegative_minimum"] < 0.0
    ):
        raise GecCcpWorkflowError(
            "upstream EEDF normalization, mean-energy, or nonnegativity failed"
        )
    for name, actual in (
        ("normalization_error_max", audit["normalization_error_max"]),
        (
            "mean_energy_relative_error_max",
            audit["mean_energy_relative_error_max"],
        ),
        ("nonnegative_minimum", audit["nonnegative_minimum"]),
    ):
        if not math.isclose(
            float(reported[name]),
            float(actual),
            rel_tol=1.0e-9,
            abs_tol=1.0e-13,
        ):
            raise GecCcpWorkflowError(
                f"upstream EEDF manifest {name} disagrees with the CSV"
            )
    return {
        **audit,
        "table_sha256": actual_sha256,
        "grid_shape": shape,
        "status": "passed",
        "active_in_comsol": False,
        "role": "independent_upstream_swarm_evidence",
    }


def _native_function_eedf_moment_audit(
    grid: ComsolFunctionEedfGrid,
) -> dict[str, Any]:
    normalizations: list[float] = []
    reconstructed_means: list[float] = []
    try:
        for requested_mean, row in zip(
            grid.mean_energies_eV, grid.values_eV_m32, strict=True
        ):
            normalization, first_moment = piecewise_linear_weighted_moments(
                grid.electron_energies_eV, row
            )
            if normalization <= 0.0:
                raise GecCcpWorkflowError(
                    "native Function-EEDF row has zero normalization"
                )
            normalizations.append(normalization)
            reconstructed_means.append(first_moment / normalization)
    except FunctionEedfError as exc:
        raise GecCcpWorkflowError(
            "native Function-EEDF piecewise-linear moment audit failed"
        ) from exc
    normalization_error = float(np.max(np.abs(np.asarray(normalizations) - 1.0)))
    mean_error = float(
        np.max(
            np.abs(np.asarray(reconstructed_means) - grid.mean_energies_eV)
            / grid.mean_energies_eV
        )
    )
    return {
        "integration": "exact_piecewise_linear_energy_segments",
        "rows": len(normalizations),
        "normalization_error_max": normalization_error,
        "mean_energy_relative_error_max": mean_error,
        "nonnegative_minimum": float(np.min(grid.values_eV_m32)),
    }


def _native_function_eedf_row(
    grid: ComsolFunctionEedfGrid, mean_energy_eV: float
) -> np.ndarray:
    """Linearly blend adjacent physical-grid rows at a requested mean energy."""
    try:
        return evaluate_comsol_function_eedf_grid(
            grid, grid.electron_energies_eV, mean_energy_eV
        )
    except ValueError as exc:
        raise GecCcpWorkflowError(
            "Function-EEDF rate audit received an invalid mean energy"
        ) from exc
