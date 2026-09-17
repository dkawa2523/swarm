"""Product comparison rows for solver smoke comparisons."""

from __future__ import annotations

import numpy as np

from electron_swarm.core.config import ComparisonConfig
from electron_swarm.core.results import SwarmCaseResult, SwarmRunResult
SCALAR_METRICS = {
    "mean_energy_eV": "mean_energy_eV_relative_difference",
    "drift_velocity_m_s": "drift_velocity_relative_difference",
    "mobility_m2_V_s": "mobility_relative_difference",
    "diffusion_L_m2_s": "diffusion_L_relative_difference",
    "net_ionization_frequency_s": "net_ionization_frequency_relative_difference",
}
ANGULAR_COMPARE_KEYS = (
    "angular_model",
    "angular_moment_source",
    "exact_dcs_based",
    "ordinary_integral_xs_closure",
)


def _relative_difference(candidate: float, reference: float) -> float:
    denom = max(abs(float(reference)), 1.0e-300)
    return float((float(candidate) - float(reference)) / denom)


def _finite_scalar(value: object) -> float | None:
    if value is None:
        return None
    number = float(value)
    return number if np.isfinite(number) else None


def _scalar_comparison(
    candidate: object,
    reference: object,
) -> tuple[float | None, str]:
    candidate_value = _finite_scalar(candidate)
    reference_value = _finite_scalar(reference)
    if candidate_value is None:
        return None, "candidate_unavailable"
    if reference_value is None:
        return None, "reference_unavailable"
    return (
        _relative_difference(candidate_value, reference_value),
        "available",
    )


def _eedf_l1_error(candidate: SwarmCaseResult, reference: SwarmCaseResult) -> float:
    """Compare the reported cell densities without point interpolation.

    Solver grids can differ substantially near thresholds and in the tail.
    Treating each reported value as a cell-average density makes the L1 norm
    conservative and symmetric with respect to the two grids.
    """

    def cell_edges(case: SwarmCaseResult) -> np.ndarray:
        centers = np.asarray(case.energy_eV, dtype=float)
        widths = np.asarray(case.energy_widths_eV, dtype=float)
        edges = np.concatenate(([0.0], np.cumsum(widths)))
        scale = max(1.0, float(edges[-1]))
        if np.any(centers < edges[:-1] - 1.0e-12 * scale) or np.any(
            centers > edges[1:] + 1.0e-12 * scale
        ):
            raise ValueError(
                f"{case.solver} EEDF centers do not lie in their reported cells"
            )
        return edges

    candidate_edges = cell_edges(candidate)
    reference_edges = cell_edges(reference)
    common_edges = np.unique(np.concatenate((candidate_edges, reference_edges)))
    midpoints = 0.5 * (common_edges[:-1] + common_edges[1:])

    def density_on_intervals(case: SwarmCaseResult, edges: np.ndarray) -> np.ndarray:
        indices = np.searchsorted(edges, midpoints, side="right") - 1
        valid = (indices >= 0) & (indices < len(case.eedf))
        values = np.zeros_like(midpoints)
        values[valid] = np.asarray(case.eedf, dtype=float)[indices[valid]]
        return values

    candidate_density = density_on_intervals(candidate, candidate_edges)
    reference_density = density_on_intervals(reference, reference_edges)
    return float(
        np.sum(np.abs(candidate_density - reference_density) * np.diff(common_edges))
    )


def _angular_value(case: SwarmCaseResult, key: str) -> object | None:
    value = case.metadata.get(key)
    if value in {"", None, "unknown"}:
        return None
    return value


def _angular_model_status(
    candidate: SwarmCaseResult,
    reference: SwarmCaseResult,
) -> str:
    candidate_values = [_angular_value(candidate, key) for key in ANGULAR_COMPARE_KEYS]
    reference_values = [_angular_value(reference, key) for key in ANGULAR_COMPARE_KEYS]
    if any(value is None for value in candidate_values + reference_values):
        return "unknown"
    if candidate_values == reference_values:
        return "match"
    return "mismatch"


def comparison_summary_rows(
    result: SwarmRunResult,
    *,
    reference_solver: str,
    candidate_solvers: list[str] | tuple[str, ...],
    compare_eedf: bool = True,
) -> list[dict[str, object]]:
    by_key: dict[tuple[str, float, str], SwarmCaseResult] = {}
    for case in result.cases:
        key = (case.case_id, float(case.e_over_n_Td), case.solver)
        if key in by_key:
            raise ValueError(f"duplicate solver result in comparison: {key!r}")
        by_key[key] = case

    keys = sorted({(case.case_id, float(case.e_over_n_Td)) for case in result.cases})
    rows: list[dict[str, object]] = []
    for case_id, eover in keys:
        reference = by_key.get((case_id, eover, reference_solver))
        if reference is None:
            continue
        for solver in candidate_solvers:
            candidate = by_key.get((case_id, eover, solver))
            if candidate is None:
                continue
            row: dict[str, object] = {
                "case_id": case_id,
                "E_over_N_Td": eover,
                "reference_solver": reference_solver,
                "candidate_solver": solver,
                "angular_model_status": _angular_model_status(candidate, reference),
            }
            for attr, column in SCALAR_METRICS.items():
                difference, status = _scalar_comparison(
                    getattr(candidate, attr),
                    getattr(reference, attr),
                )
                row[column] = difference
                row[f"{attr}_status"] = status
            if compare_eedf:
                row["eedf_l1_error"] = _eedf_l1_error(candidate, reference)
            rows.append(row)
    return rows


def comparison_reference(
    result: SwarmRunResult, comparison: ComparisonConfig
) -> str | None:
    requested = {case.solver for case in result.cases}
    reference = comparison.reference_solver
    if reference is None and "monte_carlo" in requested:
        reference = "monte_carlo"
    if reference not in requested:
        msg = (
            "comparison reference solver is not present in runnable results: "
            f"{reference!r}"
        )
        result.metadata.setdefault("comparison_warnings", []).append(msg)
        if comparison.required:
            raise ValueError(msg)
        return None
    return reference


def build_comparison_summary(
    result: SwarmRunResult,
    comparison: ComparisonConfig,
) -> list[dict[str, object]]:
    if not comparison.enabled:
        return []
    reference = comparison_reference(result, comparison)
    candidates = comparison.candidate_solvers or [
        solver for solver in sorted(result.by_solver()) if solver != reference
    ]
    if reference in candidates:
        raise ValueError("comparison reference solver cannot also be a candidate")
    if comparison.required and not candidates:
        raise ValueError("required comparison has no candidate solvers")
    summary_rows: list[dict[str, object]] = []
    if reference is not None:
        summary_rows = comparison_summary_rows(
            result,
            reference_solver=reference,
            candidate_solvers=candidates,
            compare_eedf=comparison.compare_eedf,
        )
    if comparison.required:
        reference_keys = {
            (case.case_id, float(case.e_over_n_Td))
            for case in result.cases
            if case.solver == reference
        }
        expected_rows = {
            (case_id, e_over_n_Td, solver)
            for case_id, e_over_n_Td in reference_keys
            for solver in candidates
        }
        actual_rows = {
            (
                str(row["case_id"]),
                float(row["E_over_N_Td"]),
                str(row["candidate_solver"]),
            )
            for row in summary_rows
        }
        missing = sorted(expected_rows - actual_rows)
        if missing:
            raise ValueError(
                "required comparison is missing candidate results: "
                f"{missing}"
            )
        bad_angular = [
            row
            for row in summary_rows
            if "monte_carlo"
            in {row.get("reference_solver"), row.get("candidate_solver")}
            and row.get("angular_model_status") != "match"
        ]
        if bad_angular:
            status = bad_angular[0].get("angular_model_status")
            raise ValueError(
                "comparison required same-angular PN vs MC rows, "
                f"but angular_model_status={status!r}"
            )
    result.metadata["comparison_summary_rows"] = summary_rows
    return summary_rows
