from __future__ import annotations

import numpy as np
import pytest

from electron_swarm.core.cross_sections import CrossSectionProcess, ProcessType
from electron_swarm.solvers.monte_carlo.weighted_ensemble import (
    plan_threshold_weighted_resample,
    reaction_kernel_strata_importance,
    threshold_strata_edges,
)


def test_threshold_strata_adds_guides_before_first_reaction() -> None:
    edges = threshold_strata_edges(
        (15.8, 11.6, 13.0, 0.0),
        max_energy_eV=100.0,
    )
    assert edges == pytest.approx(
        [
            0.0,
            1.45,
            2.9,
            4.35,
            5.8,
            7.25,
            8.7,
            10.15,
            11.6,
            13.0,
            14.4,
            15.8,
            100.0,
        ]
    )


def test_threshold_weighted_resample_preserves_each_stratum_weight() -> None:
    energy = np.array([1.0] * 6 + [6.0] * 3 + [9.0], dtype=float)
    weights = np.ones(10, dtype=float)
    edges = np.array([0.0, 5.0, 8.0, 12.0], dtype=float)
    original_strata = np.searchsorted(edges[1:-1], energy, side="right")

    plan = plan_threshold_weighted_resample(
        energy,
        weights,
        target_particles=20,
        strata_edges_eV=edges,
        rng=np.random.default_rng(7),
    )

    assert plan is not None
    selected_energy = energy[plan.indices]
    selected_strata = np.searchsorted(edges[1:-1], selected_energy, side="right")
    for stratum in range(3):
        before = float(np.sum(weights[original_strata == stratum]))
        after = float(np.sum(plan.weights[selected_strata == stratum]))
        assert after == pytest.approx(before, rel=1.0e-14)
    assert plan.tail_particles_before == 1
    assert plan.tail_particles_after > plan.tail_particles_before
    assert plan.indices.size == 20
    assert float(np.sum(plan.weights)) == pytest.approx(float(np.sum(weights)))


def test_threshold_weighted_resample_is_deterministic_for_seed() -> None:
    energy = np.array([1.0, 2.0, 6.0, 9.0], dtype=float)
    weights = np.array([1.0, 2.0, 0.5, 0.1], dtype=float)
    edges = np.array([0.0, 5.0, 8.0, 12.0], dtype=float)

    left = plan_threshold_weighted_resample(
        energy,
        weights,
        target_particles=12,
        strata_edges_eV=edges,
        rng=np.random.default_rng(19),
    )
    right = plan_threshold_weighted_resample(
        energy,
        weights,
        target_particles=12,
        strata_edges_eV=edges,
        rng=np.random.default_rng(19),
    )

    assert left is not None and right is not None
    assert np.array_equal(left.indices, right.indices)
    assert np.array_equal(left.weights, right.weights)


def test_threshold_weighted_resample_keeps_subnormal_stratum_representable() -> None:
    smallest = np.nextafter(0.0, 1.0)
    energy = np.array([1.0, 9.0, 9.0], dtype=float)
    weights = np.array([1.0, smallest, 0.0], dtype=float)
    edges = np.array([0.0, 5.0, 12.0], dtype=float)

    plan = plan_threshold_weighted_resample(
        energy,
        weights,
        target_particles=20,
        strata_edges_eV=edges,
        rng=np.random.default_rng(23),
    )

    assert plan is not None
    assert np.all(plan.weights > 0.0)
    assert 2 not in plan.indices
    selected_energy = energy[plan.indices]
    assert np.count_nonzero(selected_energy >= 5.0) == 1
    assert float(np.sum(plan.weights[selected_energy >= 5.0])) == smallest
    assert plan.indices.size == 20


def test_reaction_kernel_allocation_balances_contributing_strata() -> None:
    energy = np.array([1.0] * 20 + [6.0] * 3 + [9.0], dtype=float)
    weights = np.ones_like(energy)
    edges = np.array([0.0, 5.0, 8.0, 12.0], dtype=float)
    # The first observable is accumulated mainly in the middle band while the
    # second requires the highest band.  Neither should be replaced by a
    # highest-energy-only front allocation.
    importance = np.array(
        [[0.0, 1.0, 0.1], [0.0, 0.0, 1.0]],
        dtype=float,
    )

    front_only = plan_threshold_weighted_resample(
        energy,
        weights,
        target_particles=24,
        strata_edges_eV=edges,
        rng=np.random.default_rng(7),
    )
    reaction_aware = plan_threshold_weighted_resample(
        energy,
        weights,
        target_particles=24,
        strata_edges_eV=edges,
        rng=np.random.default_rng(7),
        reaction_importance=importance,
    )
    zero_kernel = plan_threshold_weighted_resample(
        energy,
        weights,
        target_particles=24,
        strata_edges_eV=edges,
        rng=np.random.default_rng(7),
        reaction_importance=np.zeros((1, 3)),
    )

    assert front_only is not None
    assert reaction_aware is not None
    assert zero_kernel is not None
    front_strata = np.searchsorted(
        edges[1:-1], energy[front_only.indices], side="right"
    )
    aware_strata = np.searchsorted(
        edges[1:-1], energy[reaction_aware.indices], side="right"
    )
    assert np.count_nonzero(aware_strata == 1) > np.count_nonzero(front_strata == 1)
    assert np.count_nonzero(aware_strata == 2) > 1
    assert np.array_equal(zero_kernel.indices, front_only.indices)
    assert np.array_equal(zero_kernel.weights, front_only.weights)


def test_reaction_importance_is_derived_from_sigma_v_not_threshold_rank() -> None:
    excitation = CrossSectionProcess(
        species="Ar",
        process="excitation",
        process_type=ProcessType.EXCITATION,
        threshold_eV=5.0,
        energy_eV=np.array([0.0, 5.0, 10.0]),
        cross_section_m2=np.array([0.0, 0.0, 2.0e-20]),
    )

    importance = reaction_kernel_strata_importance(
        (excitation,),
        np.array([0.0, 5.0, 10.0]),
    )

    assert importance is not None
    assert importance.shape == (1, 2)
    assert importance[0] == pytest.approx([0.0, 1.0])
