"""Solver comparisons with explicit numerical and observable limitations."""

from __future__ import annotations

import numpy as np
from dataclasses import fields, is_dataclass
from electron_swarm.core.observables import observable_status

from electron_swarm.core.config import ComparisonConfig
from electron_swarm.core.solver_ids import CANONICAL_SOLVER_IDS
from electron_swarm.core.results import (
    SwarmCaseResult,
    SwarmRunResult,
    case_numerical_status,
)
from electron_swarm.grids.energy import cell_density_l1, reported_cell_edges

SCALAR_METRICS = {
    "mean_energy_eV": "mean_energy_eV_relative_difference",
    "drift_velocity_m_s": "drift_velocity_relative_difference",
    "mobility_m2_V_s": "mobility_relative_difference",
    "diffusion_L_m2_s": "diffusion_L_relative_difference",
    "diffusion_T_m2_s": "diffusion_T_relative_difference",
    "net_ionization_frequency_s": "net_ionization_frequency_relative_difference",
}
ANGULAR_COMPARE_KEYS = (
    "angular_model",
    "angular_moment_source",
    "exact_dcs_based",
    "ordinary_integral_xs_closure",
    "angular_scattering_assumption",
)


def _relative_difference(candidate: float, reference: float) -> float:
    denom = abs(float(reference))
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
    if reference_value == 0.0:
        return (
            (0.0, "available") if candidate_value == 0.0 else (None, "reference_zero")
        )
    difference = _relative_difference(candidate_value, reference_value)
    return (
        (difference, "available")
        if np.isfinite(difference)
        else (None, "relative_difference_overflow")
    )


def _eedf_l1_error(candidate: SwarmCaseResult, reference: SwarmCaseResult) -> float:
    """Compare the reported cell densities without point interpolation.

    Solver grids can differ substantially near thresholds and in the tail.
    Treating each reported value as a cell-average density makes the L1 norm
    conservative and symmetric with respect to the two grids.
    """

    return cell_density_l1(
        reported_cell_edges(candidate.energy_eV, candidate.energy_widths_eV),
        candidate.eedf,
        reported_cell_edges(reference.energy_eV, reference.energy_widths_eV),
        reference.eedf,
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


def _physical_context_status(candidate, reference):
    if candidate.physical_context is None or reference.physical_context is None:
        return "unknown", "physical_context_not_provided"

    def differences_at(left, right, prefix=""):
        differences = []
        for field in fields(left):
            value, expected = getattr(left, field.name), getattr(right, field.name)
            name = prefix + field.name
            if is_dataclass(value) and is_dataclass(expected):
                differences.extend(differences_at(value, expected, name + "."))
            elif value != expected:
                differences.append(name)
        return differences

    differences = differences_at(candidate.physical_context, reference.physical_context)
    return ("mismatch", ";".join(differences)) if differences else ("match", "")


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

    selected = {reference_solver, *candidate_solvers}
    keys = {
        (case.case_id, float(case.e_over_n_Td))
        for case in result.cases
        if case.solver in selected
    }
    run = result.metadata.get("configuration", {}).get("run", {})
    if run:
        keys.update(
            (f"{run['case_prefix']}_{i:04d}", float(value))
            for i, value in enumerate(run["e_over_n_Td"])
        )
    keys = sorted(keys)
    rows: list[dict[str, object]] = []
    for case_id, eover in keys:
        reference = by_key.get((case_id, eover, reference_solver))
        for solver in candidate_solvers:
            candidate = by_key.get((case_id, eover, solver))
            if reference is None or candidate is None:
                row = {
                    "case_id": case_id,
                    "E_over_N_Td": eover,
                    "reference_solver": reference_solver,
                    "candidate_solver": solver,
                    "angular_model_status": "unavailable",
                    "physics_status": "unavailable",
                    "physics_differences": "case_missing",
                    "comparison_status": "reference_missing"
                    if reference is None
                    else "candidate_missing",
                    "reference_numerical_status": "missing"
                    if reference is None
                    else case_numerical_status(reference),
                    "candidate_numerical_status": "missing"
                    if candidate is None
                    else case_numerical_status(candidate),
                    "reference_transport_definition": None
                    if reference is None
                    else reference.transport.definition,
                    "candidate_transport_definition": None
                    if candidate is None
                    else candidate.transport.definition,
                }
                for attr, column in SCALAR_METRICS.items():
                    row[column] = None
                    row[f"{attr}_status"] = (
                        "reference_unavailable"
                        if reference is None
                        else "candidate_unavailable"
                    )
                if compare_eedf:
                    row["eedf_l1_error"] = None
                rows.append(row)
                continue
            row: dict[str, object] = {
                "case_id": case_id,
                "E_over_N_Td": eover,
                "reference_solver": reference_solver,
                "candidate_solver": solver,
                "angular_model_status": _angular_model_status(candidate, reference),
                "reference_numerical_status": case_numerical_status(reference),
                "candidate_numerical_status": case_numerical_status(candidate),
                "reference_transport_definition": reference.transport.definition,
                "candidate_transport_definition": candidate.transport.definition,
            }
            row["physics_status"], row["physics_differences"] = (
                _physical_context_status(candidate, reference)
            )
            numerical = {
                row["reference_numerical_status"],
                row["candidate_numerical_status"],
            }
            row["comparison_status"] = (
                "physics_" + row["physics_status"]
                if row["physics_status"] != "match"
                else "angular_assumptions_" + str(row["angular_model_status"])
                if row["angular_model_status"] != "match"
                else "transport_stale"
                if "postprocessed_transport_stale" in numerical
                else "requires_independent_ensemble"
                if "requires_independent_ensemble" in numerical
                else "numerically_unqualified"
                if numerical != {"numerically_converged"}
                else "numerical_comparison"
            )
            for attr, column in SCALAR_METRICS.items():
                difference, status = _scalar_comparison(
                    getattr(candidate, attr),
                    getattr(reference, attr),
                )
                if row["physics_status"] == "mismatch":
                    difference, status = None, "physics_mismatch"
                elif row["physics_status"] == "unknown":
                    difference, status = None, "physics_unknown"
                elif row["angular_model_status"] != "match":
                    difference, status = None, "angular_" + row["angular_model_status"]
                if attr in {
                    "drift_velocity_m_s",
                    "mobility_m2_V_s",
                    "diffusion_L_m2_s",
                    "diffusion_T_m2_s",
                }:
                    if difference is not None:
                        meaning = observable_status(
                            (
                                candidate.transport_observables
                                or candidate.transport.observable_definitions
                            ).get(attr),
                            (
                                reference.transport_observables
                                or reference.transport.observable_definitions
                            ).get(attr),
                        )
                        if meaning in {"definition_unknown", "definition_mismatch"}:
                            difference, status = None, meaning
                        elif meaning != "available":
                            status = meaning
                    for label, case in (
                        ("reference", reference),
                        ("candidate", candidate),
                    ):
                        if (
                            case.metadata.get("electron_electron_transport_stale")
                            is True
                        ):
                            difference, status = None, f"{label}_transport_stale"
                if (
                    status == "available"
                    and "requires_independent_ensemble" in numerical
                ):
                    status = "available_unqualified_ensemble"
                row[column] = difference
                row[f"{attr}_status"] = status
            if compare_eedf:
                row["eedf_l1_error"] = (
                    _eedf_l1_error(candidate, reference)
                    if row["physics_status"] == "match"
                    and row["angular_model_status"] == "match"
                    else None
                )
            rows.append(row)
    return rows


def build_comparison_summary(
    result: SwarmRunResult,
    comparison: ComparisonConfig,
) -> list[dict[str, object]]:
    if not comparison.enabled:
        result.metadata["comparison"] = {"status": "disabled"}
        result.metadata["comparison_summary_rows"] = []
        return []
    configured = (
        result.metadata.get("configuration", {}).get("run", {}).get("solvers", [])
    )
    requested = {row["id"] for row in configured if row["enabled"]} or set(
        result.by_solver()
    )
    reference = comparison.reference_solver
    if reference is None:
        reference = (
            "monte_carlo"
            if "monte_carlo" in requested
            else next(
                (solver for solver in CANONICAL_SOLVER_IDS if solver in requested),
                None,
            )
        )
    candidates = comparison.candidate_solvers or [
        solver
        for solver in CANONICAL_SOLVER_IDS
        if solver in requested and solver != reference
    ]
    if len(set(candidates)) != len(candidates):
        raise ValueError("comparison.candidate_solvers must be unique")
    if reference in candidates:
        raise ValueError("comparison reference solver cannot also be a candidate")
    rows = (
        comparison_summary_rows(
            result,
            reference_solver=reference,
            candidate_solvers=candidates,
            compare_eedf=comparison.compare_eedf,
        )
        if reference is not None and candidates
        else []
    )
    missing = [
        row
        for row in rows
        if row["comparison_status"] in {"reference_missing", "candidate_missing"}
    ]
    coverage = (
        "unavailable"
        if not rows or len(missing) == len(rows)
        else "partial"
        if missing
        else "complete"
    )

    def input_status(key):
        if not rows:
            return "unavailable"
        values = {row[key] for row in rows}
        return (
            "mismatch"
            if "mismatch" in values
            else "match"
            if values == {"match"}
            else "unknown"
        )

    physics_status = input_status("physics_status")
    angular_status = input_status("angular_model_status")
    status = (
        coverage
        if coverage != "complete"
        else "incompatible"
        if "mismatch" in {physics_status, angular_status}
        else "unverified"
        if {physics_status, angular_status} != {"match"}
        else "available"
    )
    numerical = {
        row[key]
        for row in rows
        for key in ("candidate_numerical_status", "reference_numerical_status")
    }
    qualification = (
        "not_evaluated"
        if status != "available"
        else "transport_stale"
        if "postprocessed_transport_stale" in numerical
        else "requires_independent_ensemble"
        if "requires_independent_ensemble" in numerical
        else "numerically_converged"
        if numerical == {"numerically_converged"}
        else "unqualified"
    )
    report = {
        "status": status,
        "coverage_status": coverage,
        "physics_status": physics_status,
        "angular_model_status": angular_status,
        "qualification_status": qualification,
        "reference_solver": reference,
        "candidate_solvers": candidates,
        "expected_pairs": len(rows),
        "available_pairs": len(rows) - len(missing),
        "missing_pairs": [
            {
                key: row[key]
                for key in (
                    "case_id",
                    "E_over_N_Td",
                    "candidate_solver",
                    "comparison_status",
                )
            }
            for row in missing
        ],
    }
    if not rows:
        report["reason"] = "no_reference_or_candidates"
    result.metadata["comparison"] = report
    result.metadata["comparison_summary_rows"] = rows
    if comparison.required:
        if coverage != "complete":
            raise ValueError(
                f"required comparison is missing candidate results or reference results: {report}"
            )
        if physics_status != "match":
            raise ValueError(
                f"required comparison has {physics_status} physical context: {report}"
            )
        if angular_status != "match":
            raise ValueError(
                f"required comparison needs matching angular assumptions, but angular_model_status={angular_status!r}"
            )
    return rows
