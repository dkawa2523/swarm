"""Histogram and tail-quality helpers for the internal Monte Carlo backend."""

from __future__ import annotations

import numpy as np

from electron_swarm.core.cross_sections import CrossSectionSet, ProcessType


def _mc_energy_edges(max_energy_eV: float) -> np.ndarray:
    max_energy = float(max_energy_eV)
    if not np.isfinite(max_energy) or max_energy <= 0.0:
        raise ValueError("physics.energy_grid_policy.max_eV_limit must be positive")
    if max_energy <= 80.0:
        return np.linspace(0.0, max_energy, 81)
    parts = [np.linspace(0.0, 80.0, 161)]
    if max_energy > 80.0:
        parts.append(np.arange(90.0, min(max_energy, 260.0) + 1.0e-9, 10.0))
    if max_energy > 260.0:
        parts.append(np.arange(280.0, max_energy + 1.0e-9, 20.0))
    edges = np.unique(np.concatenate(parts))
    if edges[-1] < max_energy:
        edges = np.append(edges, max_energy)
    return edges

def _mc_bin_relative_standard_error(counts: np.ndarray) -> np.ndarray:
    values = np.asarray(counts, dtype=float)
    out = np.full_like(values, np.nan, dtype=float)
    mask = values > 0.0
    out[mask] = 1.0 / np.sqrt(values[mask])
    return out

def _mc_effective_bin_counts(
    weighted_sum: np.ndarray,
    weighted_square_sum: np.ndarray,
) -> np.ndarray:
    values = np.asarray(weighted_sum, dtype=float)
    squares = np.asarray(weighted_square_sum, dtype=float)
    out = np.zeros_like(values, dtype=float)
    mask = squares > 0.0
    out[mask] = values[mask] * values[mask] / squares[mask]
    return out

def _build_eedf_histogram(
    edges_eV: np.ndarray,
    counts: np.ndarray,
    weighted_sum: np.ndarray,
    weighted_square_sum: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    edges = np.asarray(edges_eV, dtype=float)
    if edges.ndim != 1 or edges.size < 2 or not np.all(np.diff(edges) > 0.0):
        raise ValueError("MC EEDF histogram edges must be strictly increasing")
    widths = np.diff(edges)
    energy = 0.5 * (edges[:-1] + edges[1:])
    raw_counts = np.asarray(counts, dtype=int).copy()
    weighted = np.asarray(weighted_sum, dtype=float)
    weighted_square = np.asarray(weighted_square_sum, dtype=float)
    if raw_counts.shape != energy.shape or weighted.shape != energy.shape:
        raise ValueError("MC EEDF histogram arrays must match the edge count")
    if weighted_square.shape != energy.shape:
        raise ValueError("MC EEDF weighted-square array must match the edge count")

    effective_counts = _mc_effective_bin_counts(weighted, weighted_square)
    total_hist_weight = float(np.sum(weighted))
    eedf = weighted / max(total_hist_weight, 1.0e-300) / widths
    if not np.any(eedf > 0.0):
        eedf = np.zeros_like(energy, dtype=float)
        eedf[0] = 1.0 / widths[0]
    eedf = eedf / max(float(np.sum(eedf * widths)), 1.0e-300)
    return energy, widths, eedf, raw_counts, effective_counts

def _mc_tail_uncertainty_metadata(
    energy_eV: np.ndarray,
    counts: np.ndarray,
    threshold_eV: float,
    *,
    min_count: int = 20,
    max_weak_tail_fraction: float = 0.05,
) -> dict[str, float | int | str]:
    energy = np.asarray(energy_eV, dtype=float)
    values = np.asarray(counts, dtype=int)
    tail = energy >= float(threshold_eV)
    nonzero_tail = tail & (values > 0)
    if np.any(values >= min_count):
        max_resolved = float(np.max(energy[values >= min_count]))
    else:
        max_resolved = 0.0
    if not np.any(nonzero_tail):
        min_tail_count = 0
        status = "insufficient"
        weak_fraction = 1.0
    else:
        min_tail_count = int(np.min(values[nonzero_tail]))
        tail_total = float(np.sum(values[tail]))
        weak_total = float(np.sum(values[tail & (values < min_count)]))
        weak_fraction = weak_total / max(tail_total, 1.0e-300)
        status = (
            "ok"
            if min_tail_count >= min_count
            or weak_fraction <= max_weak_tail_fraction
            else "insufficient"
        )
    return {
        "mc_tail_uncertainty_status": status,
        "mc_min_tail_bin_count": min_tail_count,
        "mc_tail_effective_sample_count_min": float(min_tail_count),
        "mc_tail_weak_probability_fraction": float(weak_fraction),
        "mc_max_resolved_energy_eV": max_resolved,
    }

def _mc_tail_uncertainty_metadata_from_effective_counts(
    energy_eV: np.ndarray,
    effective_counts: np.ndarray,
    threshold_eV: float,
    *,
    min_count: int = 20,
    bin_probability: np.ndarray | None = None,
    max_weak_tail_fraction: float = 0.05,
) -> dict[str, float | int | str]:
    energy = np.asarray(energy_eV, dtype=float)
    values = np.asarray(effective_counts, dtype=float)
    tail = energy >= float(threshold_eV)
    nonzero_tail = tail & (values > 0.0)
    if np.any(values >= float(min_count)):
        max_resolved = float(np.max(energy[values >= float(min_count)]))
    else:
        max_resolved = 0.0
    if not np.any(nonzero_tail):
        min_tail_count = 0.0
        status = "insufficient"
        weak_fraction = 1.0
    else:
        min_tail_count = float(np.min(values[nonzero_tail]))
        weak = tail & (values < float(min_count))
        if bin_probability is None:
            mass = values
        else:
            mass = np.asarray(bin_probability, dtype=float)
        tail_total = float(np.sum(mass[tail]))
        weak_total = float(np.sum(mass[weak]))
        weak_fraction = weak_total / max(tail_total, 1.0e-300)
        status = (
            "ok"
            if min_tail_count >= float(min_count)
            or weak_fraction <= max_weak_tail_fraction
            else "insufficient"
        )
    return {
        "mc_tail_uncertainty_status": status,
        "mc_min_tail_bin_count": float(min_tail_count),
        "mc_tail_effective_sample_count_min": float(min_tail_count),
        "mc_tail_weak_probability_fraction": float(weak_fraction),
        "mc_max_resolved_energy_eV": max_resolved,
    }

def _mc_tail_comparison_status(
    *,
    energy_balance_status: str,
    tail_uncertainty_status: str,
) -> str:
    if energy_balance_status == "fail":
        return "energy_balance_fail"
    if energy_balance_status == "warning":
        return "energy_balance_warning"
    if tail_uncertainty_status != "ok":
        return "weak_tail_statistics"
    return "ok"

def _tail_threshold_eV(cross_sections: CrossSectionSet) -> float:
    thresholds = [
        float(proc.threshold_eV)
        for proc in cross_sections.processes
        if proc.process_type
        in {
            ProcessType.EXCITATION,
            ProcessType.IONIZATION,
            ProcessType.ATTACHMENT,
            ProcessType.SUPERELASTIC,
        }
        and proc.threshold_eV is not None
        and np.isfinite(proc.threshold_eV)
        and proc.threshold_eV > 0.0
    ]
    return max(thresholds) if thresholds else 20.0

def _record_histogram_residence(
    *,
    edges: np.ndarray,
    counts: np.ndarray,
    weighted_hist: np.ndarray,
    weighted_square_hist: np.ndarray,
    run_audit,
    energy_eV: float,
    contribution: float,
) -> None:
    bin_index = int(np.searchsorted(edges, energy_eV, side="right") - 1)
    if 0 <= bin_index < counts.size:
        counts[bin_index] += 1
        weighted_hist[bin_index] += contribution
        weighted_square_hist[bin_index] += contribution * contribution
        run_audit.record_histogram_sample(energy_eV)
