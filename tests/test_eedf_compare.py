from __future__ import annotations

import numpy as np

from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.diagnostics.eedf_compare import (
    FAILURE_CATEGORIES,
    classify_eedf_failure,
    compare_eedf_cases,
)


def _case(
    solver: str,
    eedf: np.ndarray,
    metadata: dict[str, object] | None = None,
) -> SwarmCaseResult:
    energy = np.array([0.5, 1.5, 2.5, 3.5], dtype=float)
    eedf = np.asarray(eedf, dtype=float)
    eedf = eedf / np.sum(eedf)
    return SwarmCaseResult(
        solver=solver,
        case_id="case_0000",
        e_over_n_Td=50.0,
        mean_energy_eV=float(np.sum(energy * eedf)),
        drift_velocity_m_s=1.0,
        mobility_m2_V_s=1.0,
        reduced_mobility_m2_V_s_m3=1.0,
        diffusion_L_m2_s=1.0,
        diffusion_T_m2_s=1.0,
        reduced_diffusion_L_m2_s_m3=1.0,
        reduced_diffusion_T_m2_s_m3=1.0,
        net_ionization_frequency_s=0.0,
        effective_townsend_m2=0.0,
        energy_eV=energy,
        eedf=eedf,
        eepf=eedf / np.sqrt(energy),
        rates=[],
        metadata=metadata or {},
    )


def test_eedf_compare_classifies_surrogate_mismatch() -> None:
    ref = _case("two_term", np.array([0.7, 0.2, 0.08, 0.02]))
    cand = _case(
        "multi_term",
        np.array([0.2, 0.2, 0.2, 0.4]),
        {"solver_method": "pn_closure_surrogate"},
    )
    comparison = compare_eedf_cases(ref, cand)
    categories = {row["category"] for row in comparison.failures}
    assert "non_independent_or_surrogate_solver_mode" in categories
    assert categories <= set(FAILURE_CATEGORIES)


def test_eedf_compare_classifies_direct_reduction_mismatch_if_present() -> None:
    ref = _case("two_term", np.array([0.7, 0.2, 0.08, 0.02]))
    cand = _case(
        "multi_term",
        np.array([0.2, 0.2, 0.2, 0.4]),
        {"solver_method": "pn_closure_direct"},
    )
    comparison = compare_eedf_cases(ref, cand)
    categories = {row["category"] for row in comparison.failures}
    assert "lmax1_reduction_failure" in categories
    assert categories <= set(FAILURE_CATEGORIES)


def test_eedf_compare_keeps_matching_cases_clean() -> None:
    ref = _case("two_term", np.array([0.7, 0.2, 0.08, 0.02]))
    cand = _case("multi_term", np.array([0.7, 0.2, 0.08, 0.02]))
    comparison = compare_eedf_cases(ref, cand)
    assert comparison.metrics["eedf_relative_l1"] < 1.0e-14
    assert comparison.failures == []


def test_failure_classifier_reports_normalization_evidence() -> None:
    rows = classify_eedf_failure(
        {"normalization_error_candidate": 1.0e-4},
        affected_solver="multi_term",
        affected_e_over_n_Td=50.0,
    )
    assert rows[0]["category"] == "normalization_or_convention_mismatch"
    assert rows[0]["affected_solver"] == "multi_term"
