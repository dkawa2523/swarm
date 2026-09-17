"""Histogram and tail-quality helpers for the internal Monte Carlo backend."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from electron_swarm.core.cross_sections import (
    CrossSectionProcess,
    ProcessType,
)


_MC_VELOCITY_CORE_LIMIT_EV = 20.0
_MC_VELOCITY_CORE_MAX_WIDTH_EV = 0.025


def _mc_energy_edges(
    max_energy_eV: float,
    thresholds_eV: Iterable[float] | None = None,
    processes: Iterable[CrossSectionProcess] | None = None,
) -> np.ndarray:
    max_energy = float(max_energy_eV)
    if not np.isfinite(max_energy) or max_energy <= 0.0:
        raise ValueError("physics.energy_grid_policy.max_eV_limit must be positive")
    core_limit = min(max_energy, _MC_VELOCITY_CORE_LIMIT_EV)
    # Electron acceleration is linear in velocity, while both residence and
    # reaction kernels contain the speed factor sqrt(E).  A uniform-speed
    # finite-volume grid therefore gives quadratic energy edges and resolves
    # the steep sub-eV distribution without imposing thousands of equally
    # small cells on the high-energy tail.  The cell next to 20 eV is at most
    # 0.025 eV wide; at the O2 vibrational thresholds it is about 0.002--0.005
    # eV wide instead of the former fixed 0.1 eV.
    core_cells = max(
        1,
        int(np.ceil(2.0 * core_limit / _MC_VELOCITY_CORE_MAX_WIDTH_EV)),
    )
    speed_coordinate = np.linspace(0.0, 1.0, core_cells + 1)
    parts = [core_limit * speed_coordinate * speed_coordinate]
    if max_energy > _MC_VELOCITY_CORE_LIMIT_EV:
        parts.append(
            _smooth_tail_edges(_MC_VELOCITY_CORE_LIMIT_EV, max_energy)[1:]
        )
    threshold_values = () if thresholds_eV is None else tuple(thresholds_eV)
    thresholds = np.asarray(threshold_values, dtype=float)
    if thresholds.ndim != 1 or np.any(~np.isfinite(thresholds)):
        raise ValueError("MC EEDF reaction thresholds must be finite")
    thresholds = thresholds[(thresholds >= 0.0) & (thresholds <= max_energy)]
    kernel_edges = _reaction_kernel_landmarks(
        () if processes is None else tuple(processes),
        max_energy_eV=max_energy,
    )
    return _coalesce_energy_edges(
        np.concatenate((*parts, kernel_edges)),
        np.concatenate((thresholds, kernel_edges)),
        max_energy_eV=max_energy,
    )


def _reaction_kernel_landmarks(
    processes: Iterable[CrossSectionProcess],
    *,
    max_energy_eV: float,
) -> np.ndarray:
    """Return physical boundaries of tabulated ``sigma(E) sqrt(E)`` kernels.

    Dense LXCat tables often contain many points on smooth segments, so using
    every cross-section knot as a histogram edge would make sampling quality
    depend on the source file's arbitrary tabulation density.  The finite-
    volume grid instead retains only support transitions and genuine extrema
    of each active rate kernel.  Piecewise-linear knots are still integrated
    exactly by ``cell_integrated_rate_coefficient``; these landmarks prevent a
    narrow peak or onset from being averaged across a cell.
    """

    maximum = float(max_energy_eV)
    landmarks: list[float] = []
    for process in processes:
        energy = np.asarray(process.energy_eV, dtype=float)
        sigma = np.asarray(process.cross_section_m2, dtype=float)
        if (
            energy.ndim != 1
            or sigma.shape != energy.shape
            or energy.size < 2
            or np.any(~np.isfinite(energy))
            or np.any(~np.isfinite(sigma))
            or np.any(np.diff(energy) <= 0.0)
            or np.any(energy < 0.0)
            or np.any(sigma < 0.0)
        ):
            raise ValueError("MC EEDF reaction kernel table is invalid")
        within = (energy >= 0.0) & (energy <= maximum)
        if not np.any(within):
            continue
        kernel = sigma * np.sqrt(energy)
        positive = kernel > 0.0
        transitions = np.flatnonzero(positive[1:] != positive[:-1])
        for index in transitions:
            landmarks.extend((float(energy[index]), float(energy[index + 1])))

        delta = np.diff(kernel)
        scale = max(float(np.max(kernel, initial=0.0)), np.finfo(float).tiny)
        significant = np.maximum(np.abs(delta[:-1]), np.abs(delta[1:])) > (
            64.0 * np.finfo(float).eps * scale
        )
        turning = significant & (delta[:-1] * delta[1:] <= 0.0)
        landmarks.extend(map(float, energy[1:-1][turning]))
        if process.threshold_eV is not None:
            landmarks.append(float(process.threshold_eV))

    if not landmarks:
        return np.empty(0, dtype=float)
    values = np.asarray(landmarks, dtype=float)
    return values[(values >= 0.0) & (values <= maximum)]


def _coalesce_energy_edges(
    grid_edges_eV: np.ndarray,
    thresholds_eV: np.ndarray,
    *,
    max_energy_eV: float,
) -> np.ndarray:
    """Merge numerically identical boundaries before histogram accumulation.

    Decimal reaction thresholds can differ by one ulp from an analytically
    identical regular-grid boundary.  Keeping both creates a cell too narrow
    to survive the centre/width representation used by workflow storage.  The
    cluster representative is chosen by physical role: domain endpoints first,
    then an exact reaction threshold, then an auxiliary grid boundary.
    """

    maximum = float(max_energy_eV)
    candidates: list[tuple[float, int]] = [
        (float(value), 2 if value == 0.0 or value == maximum else 0)
        for value in np.asarray(grid_edges_eV, dtype=float)
    ]
    candidates.extend((float(value), 1) for value in thresholds_eV)
    candidates.extend(((0.0, 2), (maximum, 2)))
    candidates.sort(key=lambda item: item[0])

    clusters: list[list[tuple[float, int]]] = []
    for candidate in candidates:
        if not clusters:
            clusters.append([candidate])
            continue
        previous = clusters[-1][-1][0]
        scale = max(1.0, abs(previous), abs(candidate[0]))
        tolerance = 64.0 * np.finfo(float).eps * scale
        if candidate[0] - previous <= tolerance:
            clusters[-1].append(candidate)
        else:
            clusters.append([candidate])

    edges = np.asarray(
        [max(cluster, key=lambda item: item[1])[0] for cluster in clusters],
        dtype=float,
    )
    if edges[0] != 0.0 or edges[-1] != maximum or np.any(np.diff(edges) <= 0.0):
        raise ValueError("MC EEDF histogram edges could not be represented")
    return edges


def _smooth_tail_edges(start_eV: float, stop_eV: float) -> np.ndarray:
    """Use a deterministic tail grid without abrupt changes in cell width."""

    current = float(start_eV)
    stop = float(stop_eV)
    step = _MC_VELOCITY_CORE_MAX_WIDTH_EV
    edges = [current]
    while current < stop:
        current = min(stop, current + step)
        edges.append(current)
        # Grow continuously from the core width while retaining a bounded
        # relative energy width in the tail.  This resolves broad reaction
        # kernels without inheriting the raw table's point density.
        step = min(20.0, 1.04 * step, 0.025 * max(current, 1.0))
    return np.asarray(edges, dtype=float)


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


def _tail_threshold_eV(processes: Iterable[CrossSectionProcess]) -> float:
    thresholds = [
        float(proc.threshold_eV)
        for proc in processes
        if proc.process_type
        in {
            ProcessType.EXCITATION,
            ProcessType.IONIZATION,
            ProcessType.SUPERELASTIC,
        }
        and proc.threshold_eV is not None
        and np.isfinite(proc.threshold_eV)
        and proc.threshold_eV > 0.0
    ]
    return max(thresholds) if thresholds else 20.0


def _record_histogram_flight(
    *,
    edges: np.ndarray,
    counts: np.ndarray,
    weighted_hist: np.ndarray,
    weighted_square_hist: np.ndarray,
    run_audit,
    energy_eV: np.ndarray,
    contributions: np.ndarray,
) -> None:
    """Accumulate one flight while retaining the flight as the ESS unit."""

    energy = np.asarray(energy_eV, dtype=float)
    weight = np.asarray(contributions, dtype=float)
    if energy.ndim != 1 or weight.shape != energy.shape or energy.size == 0:
        raise ValueError("MC histogram flight samples must be one-dimensional peers")
    bins = np.searchsorted(edges, energy, side="right") - 1
    for bin_index in np.unique(bins):
        if 0 <= int(bin_index) < counts.size:
            selected = bins == bin_index
            contribution = float(np.sum(weight[selected]))
            counts[int(bin_index)] += 1
            weighted_hist[int(bin_index)] += contribution
            # Correlated quadrature nodes from one flight are one sample.
            weighted_square_hist[int(bin_index)] += contribution * contribution
    for value in energy:
        run_audit.record_histogram_sample(float(value))
