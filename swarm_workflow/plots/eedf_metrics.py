"""Conservative-grid metrics for normalized EEDF distributions."""

from __future__ import annotations

import math
from typing import Iterable, Sequence

import numpy as np

from swarm_workflow.campaign.statistics import conservative_rebin_probability_mass
from swarm_workflow.plots.eedf_contracts import (
    DEFAULT_TAIL_THRESHOLDS_EV,
    EedfCase,
    EedfComparisonError,
    MAX_SOURCE_NORMALIZATION_ERROR,
)


def union_edges(cases: Iterable[EedfCase]) -> np.ndarray:
    values = np.sort(np.concatenate([case.edges_eV for case in cases]))
    merged = [float(values[0])]
    for raw in values[1:]:
        value = float(raw)
        tolerance = 64.0 * np.finfo(float).eps * max(1.0, abs(value), abs(merged[-1]))
        if value - merged[-1] <= tolerance:
            merged[-1] = 0.5 * (merged[-1] + value)
        else:
            merged.append(value)
    edges = np.asarray(merged, dtype=float)
    if len(edges) < 2 or np.any(np.diff(edges) <= 0.0):
        raise EedfComparisonError("common EEDF grid is not strictly increasing")
    return edges


def mass_on(case: EedfCase, edges: np.ndarray) -> np.ndarray:
    source_mass = case.density_eV_inv * case.widths_eV
    normalization = float(np.sum(source_mass))
    if (
        not math.isfinite(normalization)
        or normalization <= 0.0
        or abs(normalization - 1.0) > MAX_SOURCE_NORMALIZATION_ERROR
    ):
        raise EedfComparisonError(
            f"{case.solver} {case.e_over_n_td:g} Td EEDF is not normalized"
        )
    return conservative_rebin_probability_mass(
        case.edges_eV,
        source_mass / normalization,
        edges,
    )


def half_density_on(case: EedfCase, edges: np.ndarray) -> np.ndarray | None:
    if case.ci95_half_density_eV_inv is None:
        return None
    source_half_mass = case.ci95_half_density_eV_inv * case.widths_eV
    target_half_mass = conservative_rebin_probability_mass(
        case.edges_eV,
        source_half_mass / case.source_normalization,
        edges,
    )
    return target_half_mass / np.diff(edges)


def tail_probability(
    edges: np.ndarray,
    mass: np.ndarray,
    threshold_eV: float,
) -> float:
    left = edges[:-1]
    right = edges[1:]
    fraction = np.clip((right - float(threshold_eV)) / (right - left), 0.0, 1.0)
    return float(np.sum(mass * fraction))


def _jensen_shannon(p: np.ndarray, q: np.ndarray) -> float:
    """Return JS divergence without underflowing a positive tail mass.

    A direct ``p / ((p + q) / 2)`` evaluation can divide a positive
    subnormal value by a midpoint rounded to zero.  Log-add-exp evaluates the
    same expression without clipping, smoothing, or changing either input
    distribution.
    """

    p_values = np.asarray(p, dtype=float)
    q_values = np.asarray(q, dtype=float)
    log_two = math.log(2.0)

    def contribution(left: np.ndarray, right: np.ndarray) -> float:
        active = left > 0.0
        log_left = np.log(left[active])
        log_right = np.full(log_left.shape, -np.inf, dtype=float)
        positive_right = right[active] > 0.0
        log_right[positive_right] = np.log(right[active][positive_right])
        log_midpoint = np.logaddexp(log_left, log_right) - log_two
        return float(np.sum(left[active] * (log_left - log_midpoint)))

    value = 0.5 * contribution(p_values, q_values)
    value += 0.5 * contribution(q_values, p_values)
    if not math.isfinite(value) or value < -1.0e-15 or value > log_two + 1.0e-12:
        raise EedfComparisonError(f"invalid Jensen-Shannon divergence: {value}")
    return max(0.0, value)


def _wasserstein_one(edges: np.ndarray, p: np.ndarray, q: np.ndarray) -> float:
    """Exact W1 for piecewise-uniform densities on common cells."""

    total = 0.0
    cdf_delta = 0.0
    for width, mass_delta in zip(np.diff(edges), p - q, strict=True):
        end = cdf_delta + float(mass_delta)
        if cdf_delta == 0.0 or end == 0.0 or cdf_delta * end >= 0.0:
            total += 0.5 * (abs(cdf_delta) + abs(end)) * float(width)
        else:
            crossing = abs(cdf_delta) / (abs(cdf_delta) + abs(end))
            total += (
                0.5
                * (abs(cdf_delta) * crossing + abs(end) * (1.0 - crossing))
                * float(width)
            )
        cdf_delta = end
    return float(total)


def scaled_case(case: EedfCase) -> EedfCase:
    mean = case.reconstructed_mean_energy_eV
    return EedfCase(
        solver=case.solver,
        e_over_n_td=case.e_over_n_td,
        edges_eV=case.edges_eV / mean,
        density_eV_inv=case.density_eV_inv * mean,
        reported_mean_energy_eV=1.0,
        ci95_half_density_eV_inv=(
            None
            if case.ci95_half_density_eV_inv is None
            else case.ci95_half_density_eV_inv * mean
        ),
        uncertainty_label=case.uncertainty_label,
    )


def compare_eedf_cases(
    reference: EedfCase,
    candidate: EedfCase,
    *,
    tail_thresholds_eV: Sequence[float] = DEFAULT_TAIL_THRESHOLDS_EV,
) -> dict[str, float]:
    """Return bounded distribution distances and moment/tail differences."""

    edges = union_edges((reference, candidate))
    ref_mass = mass_on(reference, edges)
    candidate_mass = mass_on(candidate, edges)
    scaled_reference = scaled_case(reference)
    scaled_candidate = scaled_case(candidate)
    scaled_edges = union_edges((scaled_reference, scaled_candidate))
    scaled_ref_mass = mass_on(scaled_reference, scaled_edges)
    scaled_candidate_mass = mass_on(scaled_candidate, scaled_edges)
    ref_mean = reference.reconstructed_mean_energy_eV
    candidate_mean = candidate.reconstructed_mean_energy_eV
    wasserstein = _wasserstein_one(edges, ref_mass, candidate_mass)
    metrics: dict[str, float] = {
        "total_variation": float(0.5 * np.sum(np.abs(ref_mass - candidate_mass))),
        "shape_total_variation": float(
            0.5 * np.sum(np.abs(scaled_ref_mass - scaled_candidate_mass))
        ),
        "hellinger_distance": float(
            np.sqrt(0.5 * np.sum((np.sqrt(ref_mass) - np.sqrt(candidate_mass)) ** 2))
        ),
        "jensen_shannon_divergence_nats": _jensen_shannon(ref_mass, candidate_mass),
        "wasserstein_1_eV": wasserstein,
        "wasserstein_1_over_reference_mean": wasserstein / ref_mean,
        "reported_mean_energy_relative_difference": abs(
            candidate.reported_mean_energy_eV - reference.reported_mean_energy_eV
        )
        / reference.reported_mean_energy_eV,
        "reconstructed_mean_energy_relative_difference": abs(candidate_mean - ref_mean)
        / ref_mean,
        "reference_source_normalization_error": abs(
            reference.source_normalization - 1.0
        ),
        "candidate_source_normalization_error": abs(
            candidate.source_normalization - 1.0
        ),
        "reference_reconstructed_mean_energy_eV": ref_mean,
        "candidate_reconstructed_mean_energy_eV": candidate_mean,
    }
    for threshold in tail_thresholds_eV:
        label = f"{float(threshold):g}".replace(".", "p")
        reference_tail = tail_probability(edges, ref_mass, float(threshold))
        candidate_tail = tail_probability(edges, candidate_mass, float(threshold))
        metrics[f"reference_tail_probability_gt_{label}_eV"] = reference_tail
        metrics[f"candidate_tail_probability_gt_{label}_eV"] = candidate_tail
        metrics[f"tail_probability_difference_gt_{label}_eV"] = abs(
            candidate_tail - reference_tail
        )
    return metrics
