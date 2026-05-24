"""Compact EEDF comparison metrics for product benchmarks."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from electron_swarm.core.results import SwarmCaseResult


FAILURE_CATEGORIES = (
    "lmax1_reduction_failure",
    "non_independent_or_surrogate_solver_mode",
    "normalization_or_convention_mismatch",
    "energy_grid_or_boundary_mismatch",
    "cross_section_projection_mismatch",
    "field_coupling_error",
    "collision_moment_or_angular_closure_error",
    "source_sink_model_mismatch",
    "nonphysical_solver_mode",
    "MC_statistical_uncertainty",
)


@dataclass(frozen=True, slots=True)
class EedfComparison:
    metrics: dict[str, float]
    failures: list[dict[str, Any]]


def cell_widths_from_centers(energy_eV: np.ndarray) -> np.ndarray:
    energy = np.asarray(energy_eV, dtype=float)
    if energy.ndim != 1 or len(energy) < 2 or np.any(np.diff(energy) <= 0.0):
        raise ValueError("energy grid must be one-dimensional and strictly increasing")
    edges = np.empty(len(energy) + 1, dtype=float)
    edges[1:-1] = 0.5 * (energy[:-1] + energy[1:])
    edges[0] = max(0.0, energy[0] - 0.5 * (energy[1] - energy[0]))
    edges[-1] = energy[-1] + 0.5 * (energy[-1] - energy[-2])
    return np.diff(edges)


def normalize_eedf(energy_eV: np.ndarray, eedf: np.ndarray) -> tuple[np.ndarray, float]:
    widths = cell_widths_from_centers(energy_eV)
    values = np.asarray(eedf, dtype=float)
    if values.shape != widths.shape or not np.all(np.isfinite(values)):
        raise ValueError("EEDF must be finite and match the energy grid")
    norm = float(np.sum(values * widths))
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError("EEDF normalization must be positive and finite")
    return values / norm, norm


def _energy_quantile(energy: np.ndarray, eedf: np.ndarray, widths: np.ndarray, q: float) -> float:
    cdf = np.cumsum(eedf * widths)
    cdf /= max(float(cdf[-1]), 1.0e-300)
    return float(np.interp(q, cdf, energy))


def _major_rate_difference(reference: SwarmCaseResult, candidate: SwarmCaseResult) -> float:
    ref = {
        (rate.species, rate.process, rate.process_type): rate.mixture_weighted_rate_m3_s
        for rate in reference.rates
        if abs(rate.mixture_weighted_rate_m3_s) > 0.0
    }
    if not ref:
        return 0.0
    cand = {
        (rate.species, rate.process, rate.process_type): rate.mixture_weighted_rate_m3_s
        for rate in candidate.rates
    }
    scale = max(abs(value) for value in ref.values())
    important = {key: value for key, value in ref.items() if abs(value) >= 0.01 * scale}
    if not important:
        important = ref
    return float(
        max(
            abs(cand.get(key, 0.0) - value) / max(abs(value), 1.0e-300)
            for key, value in important.items()
        )
    )


def eedf_metrics(reference: SwarmCaseResult, candidate: SwarmCaseResult) -> dict[str, float]:
    ref_energy = np.asarray(reference.energy_eV, dtype=float)
    cand_energy = np.asarray(candidate.energy_eV, dtype=float)
    ref_f, ref_norm = normalize_eedf(ref_energy, reference.eedf)
    cand_f, cand_norm = normalize_eedf(cand_energy, candidate.eedf)
    widths = cell_widths_from_centers(ref_energy)
    cand_on_ref = np.interp(ref_energy, cand_energy, cand_f, left=0.0, right=0.0)
    cand_on_ref, _ = normalize_eedf(ref_energy, cand_on_ref)

    diff = cand_on_ref - ref_f
    l1 = float(np.sum(np.abs(diff) * widths))
    l2 = float(np.sqrt(np.sum(diff * diff * widths)))
    tail_threshold = max(20.0, 0.75 * float(ref_energy[-1]))
    tail = ref_energy >= tail_threshold
    if np.any(tail):
        ref_tail = float(np.sum(ref_f[tail] * widths[tail]))
        cand_tail = float(np.sum(cand_on_ref[tail] * widths[tail]))
        log_tail_error = float(
            np.mean(
                np.abs(
                    np.log(np.maximum(cand_on_ref[tail], 1.0e-300))
                    - np.log(np.maximum(ref_f[tail], 1.0e-300))
                )
            )
        )
    else:
        ref_tail = cand_tail = log_tail_error = 0.0
    q50_ref = _energy_quantile(ref_energy, ref_f, widths, 0.50)
    q90_ref = _energy_quantile(ref_energy, ref_f, widths, 0.90)
    q99_ref = _energy_quantile(ref_energy, ref_f, widths, 0.99)
    q50_cand = _energy_quantile(ref_energy, cand_on_ref, widths, 0.50)
    q90_cand = _energy_quantile(ref_energy, cand_on_ref, widths, 0.90)
    q99_cand = _energy_quantile(ref_energy, cand_on_ref, widths, 0.99)
    mean_ref = float(np.sum(ref_energy * ref_f * widths))
    mean_cand = float(np.sum(ref_energy * cand_on_ref * widths))
    return {
        "eedf_l1": l1,
        "eedf_relative_l1": l1 / max(float(np.sum(np.abs(ref_f) * widths)), 1.0e-300),
        "eedf_l2_weighted": l2,
        "log_tail_error": log_tail_error,
        "mean_energy_relative_difference": abs(candidate.mean_energy_eV - reference.mean_energy_eV)
        / max(abs(reference.mean_energy_eV), 1.0e-300),
        "eedf_mean_energy_relative_difference": abs(mean_cand - mean_ref)
        / max(abs(mean_ref), 1.0e-300),
        "tail_probability_difference": abs(cand_tail - ref_tail),
        "rate_weighted_eedf_error": _major_rate_difference(reference, candidate),
        "drift_velocity_relative_difference": abs(
            candidate.drift_velocity_m_s - reference.drift_velocity_m_s
        )
        / max(abs(reference.drift_velocity_m_s), 1.0e-300),
        "major_rate_relative_difference": _major_rate_difference(reference, candidate),
        "normalization_error_reference": abs(ref_norm - 1.0),
        "normalization_error_candidate": abs(cand_norm - 1.0),
        "E50_difference_eV": abs(q50_cand - q50_ref),
        "E90_difference_eV": abs(q90_cand - q90_ref),
        "E99_difference_eV": abs(q99_cand - q99_ref),
    }


def classify_eedf_failure(
    metrics: dict[str, float],
    *,
    affected_solver: str,
    affected_e_over_n_Td: float,
    candidate_metadata: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    def add(category: str, evidence: str, code_area: str, fix: str, severity: str) -> None:
        rows.append(
            {
                "category": category,
                "evidence": evidence,
                "affected_E_over_N_Td": affected_e_over_n_Td,
                "affected_solver": affected_solver,
                "suspected_code_area": code_area,
                "recommended_fix": fix,
                "severity": severity,
            }
        )

    norm = max(
        metrics.get("normalization_error_reference", 0.0),
        metrics.get("normalization_error_candidate", 0.0),
    )
    if norm > 1.0e-8:
        add(
            "normalization_or_convention_mismatch",
            f"normalization_error={norm:.3e}",
            "EEDF normalization/convention",
            "normalize all outputs as F(E) in 1/eV before comparison",
            "high",
        )
    solver_method = str((candidate_metadata or {}).get("solver_method", ""))
    eedf_l1 = metrics.get("eedf_relative_l1", 0.0)
    if affected_solver == "monte_carlo" and eedf_l1 > 0.05:
        add(
            "MC_statistical_uncertainty",
            f"eedf_relative_l1={eedf_l1:.3e}",
            "Monte Carlo sampling / histogram statistics",
            "increase particles, collisions, and seed ensemble before declaring solver mismatch",
            "medium",
        )
    elif affected_solver == "multi_term" and solver_method == "pn_closure_surrogate" and eedf_l1 > 0.01:
        add(
            "non_independent_or_surrogate_solver_mode",
            f"eedf_relative_l1={eedf_l1:.3e}; solver_method={solver_method}",
            "multi_term surrogate closure",
            "do not treat pn_closure_surrogate as an independent direct PN validation path",
            "medium",
        )
    elif solver_method == "pn_closure_direct" and eedf_l1 > 0.01:
        add(
            "lmax1_reduction_failure",
            f"eedf_relative_l1={eedf_l1:.3e}",
            "multi_term direct PN operator",
            "inspect l=0/l=1 field coupling, source/sink, and boundary rows",
            "high",
        )
    if metrics.get("drift_velocity_relative_difference", 0.0) > 0.01:
        add(
            "field_coupling_error",
            f"drift_velocity_relative_difference={metrics['drift_velocity_relative_difference']:.3e}",
            "field coupling / momentum damping",
            "compare momentum relaxation and E-field scaling against two_term",
            "high",
        )
    if metrics.get("major_rate_relative_difference", 0.0) > 0.02:
        add(
            "source_sink_model_mismatch",
            f"major_rate_relative_difference={metrics['major_rate_relative_difference']:.3e}",
            "rate convolution or inelastic source/sink",
            "recompute rates from solved f0 and compare threshold interpolation",
            "medium",
        )
    return rows


def compare_eedf_cases(reference: SwarmCaseResult, candidate: SwarmCaseResult) -> EedfComparison:
    metrics = eedf_metrics(reference, candidate)
    failures = classify_eedf_failure(
        metrics,
        affected_solver=candidate.solver,
        affected_e_over_n_Td=candidate.e_over_n_Td,
        candidate_metadata=candidate.metadata,
    )
    return EedfComparison(metrics=metrics, failures=failures)
