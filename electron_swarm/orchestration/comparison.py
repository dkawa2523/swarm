"""Product comparison rows for solver smoke comparisons."""

from __future__ import annotations

import numpy as np

from electron_swarm.core.results import SwarmCaseResult, SwarmRunResult
from electron_swarm.diagnostics.common import widths_from_centers


SCALAR_METRICS = {
    "mean_energy_eV": "mean_energy_eV_relative_difference",
    "drift_velocity_m_s": "drift_velocity_relative_difference",
    "mobility_m2_V_s": "mobility_relative_difference",
    "diffusion_L_m2_s": "diffusion_L_relative_difference",
    "net_ionization_frequency_s": "net_ionization_frequency_relative_difference",
}


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
            }
            for attr, column in SCALAR_METRICS.items():
                row[column] = _relative_difference(
                    getattr(candidate, attr),
                    getattr(reference, attr),
                )
            if compare_eedf:
                row["eedf_l1_error"] = _eedf_l1_error(candidate, reference)
            rows.append(row)
    return rows
