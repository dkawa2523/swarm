"""Threshold-aware weighted-ensemble resampling for rare electron energies.

The resampler changes only how a fixed Monte Carlo budget is distributed.  It
preserves the total statistical weight in every occupied energy stratum, so it
does not introduce a tail shape or borrow a closure from another solver.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from electron_swarm.core.constants import ELECTRON_MASS_KG, EV_TO_J
from electron_swarm.core.cross_sections import CrossSectionProcess, ProcessType


_REACTION_TYPES = {
    ProcessType.EXCITATION,
    ProcessType.IONIZATION,
    ProcessType.SUPERELASTIC,
}

# Tail resampling follows the null-collision clock, not the unrelated transport
# correlation lag.  Sixteen trial periods leave descendants time to decorrelate
# while retaining repeated opportunities to advance a rare reaction front.
REACTION_WE_RESAMPLE_TRIAL_PERIODS = 16


@dataclass(frozen=True, slots=True)
class WeightedEnsemblePlan:
    """Particle selection and replacement weights for one resampling step."""

    indices: np.ndarray
    weights: np.ndarray
    occupied_strata: int
    tail_particles_before: int
    tail_particles_after: int


def reaction_kernel_strata_importance(
    processes: Sequence[CrossSectionProcess],
    strata_edges_eV: np.ndarray,
) -> np.ndarray | None:
    """Return one normalized ``sigma(E) v(E)`` profile per reaction.

    Threshold-only strata say where a rare reaction becomes possible, but not
    where its rate integral is accumulated.  These profiles let the resampler
    reserve descendants for every active reaction kernel instead of sending
    the entire exploration budget to the highest occupied energy stratum.
    Each row is normalized independently, so a small but physically requested
    channel is not erased by a larger cross section.
    """

    edges = np.asarray(strata_edges_eV, dtype=float)
    if (
        edges.ndim != 1
        or edges.size < 2
        or np.any(~np.isfinite(edges))
        or np.any(np.diff(edges) <= 0.0)
    ):
        raise ValueError("reaction-kernel stratum edges must strictly increase")

    rows: list[np.ndarray] = []
    for process in processes:
        if process.process_type not in _REACTION_TYPES:
            continue
        table_energy = np.asarray(process.energy_eV, dtype=float)
        # Cross-section knots plus stratum boundaries capture every linear
        # segment and every physical threshold without a dense auxiliary grid.
        sample_energy = np.unique(
            np.concatenate(
                (
                    edges,
                    table_energy[
                        (table_energy >= edges[0]) & (table_energy <= edges[-1])
                    ],
                )
            )
        )
        speed = np.sqrt(
            np.maximum(2.0 * sample_energy * EV_TO_J / ELECTRON_MASS_KG, 0.0)
        )
        kernel = np.asarray(process.sigma(sample_energy), dtype=float) * speed
        peak = float(np.max(kernel))
        if not np.isfinite(peak) or peak <= 0.0:
            continue
        row = np.zeros(edges.size - 1, dtype=float)
        for index, (left, right) in enumerate(zip(edges[:-1], edges[1:], strict=True)):
            in_stratum = (sample_energy >= left) & (sample_energy <= right)
            if np.any(in_stratum):
                row[index] = float(np.max(kernel[in_stratum])) / peak
        rows.append(row)
    if not rows:
        return None
    return np.asarray(rows, dtype=float)


def threshold_strata_edges(
    thresholds_eV: tuple[float, ...],
    *,
    max_energy_eV: float,
) -> np.ndarray:
    """Build generic energy strata leading into the physical thresholds.

    Three guide levels below the first positive threshold let a rare path be
    replicated before it has already crossed the reaction threshold.  All
    later thresholds remain explicit boundaries, which also supports
    state-resolved cross-section sets without model-specific constants.
    """

    upper = float(max_energy_eV)
    if not np.isfinite(upper) or upper <= 0.0:
        raise ValueError("weighted-ensemble maximum energy must be positive")
    positive = sorted(
        {
            float(value)
            for value in thresholds_eV
            if np.isfinite(value) and 0.0 < float(value) < upper
        }
    )
    if not positive:
        return np.array([0.0, upper], dtype=float)
    first = positive[0]
    reference_step = first / 8.0
    guides: list[float] = []
    lower = 0.0
    for threshold in positive:
        gap = threshold - lower
        intervals = min(8, max(1, int(np.ceil(gap / reference_step))))
        guides.extend(lower + gap * index / intervals for index in range(1, intervals))
        lower = threshold
    return np.asarray(sorted({0.0, *guides, *positive, upper}), dtype=float)


def plan_threshold_weighted_resample(
    energies_eV: np.ndarray,
    weights: np.ndarray,
    *,
    target_particles: int,
    strata_edges_eV: np.ndarray,
    rng: np.random.Generator,
    reaction_importance: np.ndarray | None = None,
    exploration_fraction: float = 0.5,
) -> WeightedEnsemblePlan | None:
    """Plan unbiased fixed-population resampling within occupied strata.

    Half of the budget follows probability mass.  The exploration half is
    shared between the advancing rare-event front and independently normalized
    reaction kernels.  This is a multiple-observable allocation: excitation is
    not starved merely because an ionization stratum happens to be the highest
    occupied one.  Exact stratum weights are retained by construction.
    """

    energy = np.asarray(energies_eV, dtype=float)
    particle_weights = np.asarray(weights, dtype=float)
    edges = np.asarray(strata_edges_eV, dtype=float)
    target = int(target_particles)
    explore = float(exploration_fraction)
    if energy.ndim != 1 or particle_weights.shape != energy.shape:
        raise ValueError("weighted-ensemble energies and weights must be 1D peers")
    if target <= 0:
        raise ValueError("weighted-ensemble target population must be positive")
    if (
        edges.ndim != 1
        or edges.size < 2
        or not np.all(np.isfinite(edges))
        or edges[0] != 0.0
        or np.any(np.diff(edges) <= 0.0)
    ):
        raise ValueError("weighted-ensemble stratum edges must increase from zero")
    if not 0.0 <= explore <= 1.0:
        raise ValueError("weighted-ensemble exploration fraction must be in [0, 1]")
    kernels: np.ndarray | None = None
    if reaction_importance is not None:
        kernels = np.asarray(reaction_importance, dtype=float)
        if (
            kernels.ndim != 2
            or kernels.shape[1] != edges.size - 1
            or np.any(~np.isfinite(kernels))
            or np.any(kernels < 0.0)
        ):
            raise ValueError(
                "weighted-ensemble reaction importance must be a finite "
                "nonnegative reaction-by-stratum array"
            )
    if (
        energy.size == 0
        or np.any(~np.isfinite(energy))
        or np.any(energy < 0.0)
        or np.any(~np.isfinite(particle_weights))
        or np.any(particle_weights < 0.0)
        or not np.any(particle_weights > 0.0)
    ):
        raise ValueError("weighted-ensemble particles must have finite positive state")

    # A long sequence of unbiased splits can eventually round a statistically
    # negligible descendant to exactly zero in binary64.  Zero-weight states
    # carry no measure, so exclude them from the next selection instead of
    # aborting an otherwise valid tail trajectory.
    positive_indices = np.flatnonzero(particle_weights > 0.0)
    active_energy = energy[positive_indices]
    active_weights = particle_weights[positive_indices]
    strata = np.searchsorted(edges[1:-1], active_energy, side="right")
    occupied = np.unique(strata)
    # With one occupied stratum this reduces to ordinary resampling and cannot
    # advance a rare front.  If the target cannot retain every occupied
    # stratum, leave the ensemble unchanged instead of deleting support.
    if occupied.size < 2 or occupied.size > target:
        return None

    stratum_weights = np.asarray(
        [float(np.sum(active_weights[strata == item])) for item in occupied],
        dtype=float,
    )
    total_weight = float(np.sum(stratum_weights))
    mass_fraction = stratum_weights / total_weight
    ranks_from_front = np.arange(occupied.size, dtype=float) - (occupied.size - 1)
    front_score = np.exp2(ranks_from_front)
    front_score /= float(np.sum(front_score))
    exploration_score = front_score
    if kernels is not None:
        occupied_kernels = kernels[:, occupied]
        active_rows = np.sum(occupied_kernels, axis=1) > 0.0
        normalized = np.empty((0, occupied.size), dtype=float)
        if np.any(active_rows):
            # Neyman allocation for a weighted stratum scales with both the
            # observable magnitude and the probability mass represented by
            # that stratum.  Without this factor an extremely improbable,
            # high-energy ionization stratum can consume the exploration
            # budget needed by the excitation band that actually contributes
            # to the requested rate integral.
            normalized = occupied_kernels[active_rows] * stratum_weights[None, :]
            active_rows = np.sum(normalized, axis=1) > 0.0
            normalized = normalized[active_rows]
        if normalized.size:
            normalized /= np.sum(normalized, axis=1)[:, None]
            kernel_score = np.mean(normalized, axis=0)
            # Retain a front component for strata below every reaction
            # threshold; their instantaneous kernels are zero but they are
            # the only paths by which a rare reaction can first be reached.
            exploration_score = 0.25 * front_score + 0.75 * kernel_score
            exploration_score /= float(np.sum(exploration_score))
    allocation_score = explore * exploration_score + (1.0 - explore) * mass_fraction
    # Retain every occupied stratum, then distribute the remaining integer
    # budget with the largest-remainder method.
    remaining = target - int(occupied.size)
    raw_extra = remaining * allocation_score
    allocated = np.ones(occupied.size, dtype=int) + np.floor(raw_extra).astype(int)
    remainder = target - int(np.sum(allocated))
    if remainder:
        order = np.argsort(-(raw_extra - np.floor(raw_extra)), kind="stable")
        allocated[order[:remainder]] += 1

    # Do not split a representable subnormal stratum weight into zero-weight
    # descendants.  Retain one carrier for such a stratum and move the freed
    # particle slots to strata whose per-particle weight remains representable.
    freed = 0
    for index, (count, stratum_weight) in enumerate(
        zip(allocated, stratum_weights, strict=True)
    ):
        if stratum_weight / int(count) == 0.0:
            freed += int(count) - 1
            allocated[index] = 1
    if freed:
        order = np.argsort(-allocation_score, kind="stable")
        while freed:
            assigned = False
            for index in order:
                proposed = int(allocated[index]) + 1
                if stratum_weights[index] / proposed == 0.0:
                    continue
                allocated[index] = proposed
                freed -= 1
                assigned = True
                if not freed:
                    break
            if not assigned:
                raise FloatingPointError(
                    "weighted-ensemble allocation cannot retain positive weights"
                )

    selected: list[np.ndarray] = []
    replacement_weights: list[np.ndarray] = []
    for item, count, stratum_weight in zip(
        occupied, allocated, stratum_weights, strict=True
    ):
        active_candidates = np.flatnonzero(strata == item)
        candidates = positive_indices[active_candidates]
        local_weights = active_weights[active_candidates]
        cumulative = np.cumsum(local_weights)
        step = stratum_weight / int(count)
        positions = (float(rng.random()) + np.arange(int(count), dtype=float)) * step
        local_indices = np.searchsorted(cumulative, positions, side="right")
        selected.append(candidates[np.clip(local_indices, 0, candidates.size - 1)])
        replacement_weights.append(np.full(int(count), step, dtype=float))

    indices = np.concatenate(selected)
    new_weights = np.concatenate(replacement_weights)
    permutation = rng.permutation(target)
    indices = indices[permutation]
    new_weights = new_weights[permutation]
    highest = int(occupied[-1])
    return WeightedEnsemblePlan(
        indices=indices,
        weights=new_weights,
        occupied_strata=int(occupied.size),
        tail_particles_before=int(np.count_nonzero(strata == highest)),
        tail_particles_after=int(allocated[-1]),
    )
