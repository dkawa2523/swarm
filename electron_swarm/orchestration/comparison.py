"""Product comparison rows for solver smoke comparisons."""

from __future__ import annotations

import numpy as np

from electron_swarm.core.config import ComparisonConfig
from electron_swarm.core.results import SwarmCaseResult, SwarmRunResult
from electron_swarm.diagnostics.common import widths_from_centers


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


def _eedf_l1_error(candidate: SwarmCaseResult, reference: SwarmCaseResult) -> float:
    ref_energy = np.asarray(reference.energy_eV, dtype=float)
    ref_eedf = np.asarray(reference.eedf, dtype=float)
    candidate_eedf = np.interp(
        ref_energy,
        np.asarray(candidate.energy_eV, dtype=float),
        np.asarray(candidate.eedf, dtype=float),
        left=0.0,
        right=0.0,
    )
    widths = widths_from_centers(ref_energy)
    return float(np.sum(np.abs(candidate_eedf - ref_eedf) * widths))


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


def _angular_model_warning(status: str) -> str:
    if status == "match":
        return ""
    if status == "unknown":
        return "angular_model_unknown"
    return "angular_model_mismatch"


def _angular_sampler_treatment(
    candidate: SwarmCaseResult,
    reference: SwarmCaseResult,
) -> str:
    for case in (candidate, reference):
        if case.solver == "monte_carlo":
            treatment = case.metadata.get("effective_angular_scattering")
            if treatment not in {"", None}:
                return str(treatment)
            treatment = case.metadata.get("monte_carlo_angular_scattering")
            if treatment not in {"", None}:
                return str(treatment)
    return "not_applicable"


def comparison_summary_rows(
    result: SwarmRunResult,
    *,
    reference_solver: str,
    candidate_solvers: list[str] | tuple[str, ...],
    compare_eedf: bool = True,
) -> list[dict[str, object]]:
    by_key: dict[tuple[str, float, str], SwarmCaseResult] = {}
    for case in result.cases:
        by_key[(case.case_id, float(case.e_over_n_Td), case.solver)] = case

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
                "reference_angular_model": reference.metadata.get("angular_model", ""),
                "candidate_angular_model": candidate.metadata.get("angular_model", ""),
                "reference_angular_moment_source": reference.metadata.get(
                    "angular_moment_source", ""
                ),
                "candidate_angular_moment_source": candidate.metadata.get(
                    "angular_moment_source", ""
                ),
            }
            row["same_angular_model"] = row["angular_model_status"] == "match"
            row["angular_model_warning"] = _angular_model_warning(
                str(row["angular_model_status"])
            )
            row["angular_model_reference"] = row["reference_angular_model"]
            row["angular_model_candidate"] = row["candidate_angular_model"]
            row["angular_sampler_treatment"] = _angular_sampler_treatment(
                candidate,
                reference,
            )
            row["angular_model_mismatch_reason"] = row["angular_model_warning"]
            if row["same_angular_model"]:
                row["angular_model"] = row["reference_angular_model"]
            else:
                row["angular_model"] = (
                    f"{row['reference_angular_model']}->{row['candidate_angular_model']}"
                )
            for attr, column in SCALAR_METRICS.items():
                row[column] = _relative_difference(
                    getattr(candidate, attr),
                    getattr(reference, attr),
                )
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
    summary_rows: list[dict[str, object]] = []
    if reference is not None:
        summary_rows = comparison_summary_rows(
            result,
            reference_solver=reference,
            candidate_solvers=candidates,
            compare_eedf=comparison.compare_eedf,
        )
    if comparison.required:
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
