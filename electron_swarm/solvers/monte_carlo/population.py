"""Particle population storage and resampling state for internal Monte Carlo."""

from __future__ import annotations

import numpy as np

from electron_swarm.core.constants import ELECTRON_MASS_KG, EV_TO_J
from electron_swarm.solvers.monte_carlo.kinematics import (
    energy_from_velocity_eV,
    random_direction,
)
from electron_swarm.solvers.monte_carlo.weighted_ensemble import (
    WeightedEnsemblePlan,
    plan_threshold_weighted_resample,
)


class _ParticleEnsemble:
    """Structure-of-arrays particle storage with geometric spare capacity.

    Ionization may create several daughters inside one synchronization
    barrier.  Keeping capacity separate from the active size avoids copying
    every particle array for every daughter while preserving the same active
    array interface used by the reference solver.
    """

    __slots__ = (
        "_positions",
        "_velocities",
        "_times",
        "_weights",
        "_lineages",
        "_size",
    )

    def __init__(
        self,
        positions: np.ndarray,
        velocities: np.ndarray,
        times: np.ndarray,
        weights: np.ndarray,
        lineages: np.ndarray | None = None,
    ) -> None:
        active_positions = np.asarray(positions, dtype=float)
        active_velocities = np.asarray(velocities, dtype=float)
        active_times = np.asarray(times, dtype=float)
        active_weights = np.asarray(weights, dtype=float)
        size = int(active_weights.size)
        if (
            active_positions.shape != (size, 3)
            or active_velocities.shape != (size, 3)
            or active_times.shape != (size,)
        ):
            raise ValueError("particle state arrays must have matching active sizes")
        active_lineages = (
            np.arange(size, dtype=np.int64)
            if lineages is None
            else np.asarray(lineages, dtype=np.int64)
        )
        if active_lineages.shape != active_weights.shape:
            raise ValueError("particle lineage ids must match ensemble weights")
        self._replace_state(
            active_positions,
            active_velocities,
            active_times,
            active_weights,
            active_lineages,
        )

    @staticmethod
    def _capacity_for(size: int) -> int:
        return max(1, 2 * int(size))

    def _replace_state(
        self,
        positions: np.ndarray,
        velocities: np.ndarray,
        times: np.ndarray,
        weights: np.ndarray,
        lineages: np.ndarray,
    ) -> None:
        size = int(np.asarray(weights).size)
        capacity = self._capacity_for(size)
        self._positions = np.empty((capacity, 3), dtype=float)
        self._velocities = np.empty((capacity, 3), dtype=float)
        self._times = np.empty(capacity, dtype=float)
        self._weights = np.empty(capacity, dtype=float)
        self._lineages = np.empty(capacity, dtype=np.int64)
        self._positions[:size] = positions
        self._velocities[:size] = velocities
        self._times[:size] = times
        self._weights[:size] = weights
        self._lineages[:size] = lineages
        self._size = size

    @property
    def positions(self) -> np.ndarray:
        return self._positions[: self._size]

    @property
    def velocities(self) -> np.ndarray:
        return self._velocities[: self._size]

    @property
    def times(self) -> np.ndarray:
        return self._times[: self._size]

    @property
    def weights(self) -> np.ndarray:
        return self._weights[: self._size]

    @property
    def lineages(self) -> np.ndarray:
        return self._lineages[: self._size]

    @lineages.setter
    def lineages(self, values: np.ndarray) -> None:
        lineages = np.asarray(values, dtype=np.int64)
        if lineages.shape != (self._size,):
            raise ValueError("particle lineage ids must match ensemble weights")
        self._lineages[: self._size] = lineages

    @classmethod
    def initialize(
        cls,
        particles: int,
        initial_speed: float,
        rng: np.random.Generator,
    ) -> "_ParticleEnsemble":
        return cls(
            positions=np.zeros((particles, 3), dtype=float),
            velocities=np.array(
                [random_direction(rng) * initial_speed for _ in range(particles)],
                dtype=float,
            ),
            times=np.zeros(particles, dtype=float),
            weights=np.ones(particles, dtype=float),
            lineages=np.arange(particles, dtype=np.int64),
        )

    def __len__(self) -> int:
        return self._size

    def copy(self) -> "_ParticleEnsemble":
        return _ParticleEnsemble(
            positions=self.positions.copy(),
            velocities=self.velocities.copy(),
            times=self.times.copy(),
            weights=self.weights.copy(),
            lineages=self.lineages.copy(),
        )

    def total_weighted_energy_eV(self) -> float:
        velocity_square = np.einsum(
            "ij,ij->i",
            self.velocities,
            self.velocities,
        )
        return float(
            0.5 * ELECTRON_MASS_KG / EV_TO_J * np.sum(self.weights * velocity_square)
        )

    def _ensure_capacity(self, required: int) -> None:
        if required <= self._weights.size:
            return
        capacity = max(int(required), 2 * int(self._weights.size))
        positions = np.empty((capacity, 3), dtype=float)
        velocities = np.empty((capacity, 3), dtype=float)
        times = np.empty(capacity, dtype=float)
        weights = np.empty(capacity, dtype=float)
        lineages = np.empty(capacity, dtype=np.int64)
        positions[: self._size] = self.positions
        velocities[: self._size] = self.velocities
        times[: self._size] = self.times
        weights[: self._size] = self.weights
        lineages[: self._size] = self.lineages
        self._positions = positions
        self._velocities = velocities
        self._times = times
        self._weights = weights
        self._lineages = lineages

    def append_particle(
        self,
        *,
        position: np.ndarray,
        velocity: np.ndarray,
        time_s: float,
        weight: float,
        lineage: int,
    ) -> None:
        self._ensure_capacity(self._size + 1)
        index = self._size
        self._positions[index] = np.asarray(position, dtype=float)
        self._velocities[index] = np.asarray(velocity, dtype=float)
        self._times[index] = float(time_s)
        self._weights[index] = float(weight)
        self._lineages[index] = int(lineage)
        self._size += 1

    def normalize_total_weight(self, target_total_weight: float) -> None:
        total = float(np.sum(self.weights))
        if total > 0.0:
            self._weights[: self._size] *= float(target_total_weight) / total

    def systematic_resample(
        self,
        target_particles: int,
        rng: np.random.Generator,
    ) -> bool:
        target = int(target_particles)
        if target <= 0 or len(self) <= target:
            return False
        total = float(np.sum(self.weights))
        if not np.isfinite(total) or total <= 0.0:
            return False
        cumulative = np.cumsum(self.weights)
        positions = (float(rng.random()) + np.arange(target, dtype=float)) * (
            total / target
        )
        indices = np.searchsorted(cumulative, positions, side="right")
        indices = np.clip(indices, 0, len(self) - 1)
        self._replace_state(
            self.positions[indices],
            self.velocities[indices],
            self.times[indices],
            np.full(target, total / target, dtype=float),
            self.lineages[indices],
        )
        return True

    def threshold_weighted_resample(
        self,
        target_particles: int,
        strata_edges_eV: np.ndarray,
        rng: np.random.Generator,
        *,
        reaction_importance: np.ndarray | None = None,
    ) -> WeightedEnsemblePlan | None:
        energies = np.asarray(
            [energy_from_velocity_eV(velocity) for velocity in self.velocities],
            dtype=float,
        )
        plan = plan_threshold_weighted_resample(
            energies,
            self.weights,
            target_particles=target_particles,
            strata_edges_eV=strata_edges_eV,
            rng=rng,
            reaction_importance=reaction_importance,
        )
        if plan is None:
            return None
        self._replace_state(
            self.positions[plan.indices],
            self.velocities[plan.indices],
            self.times[plan.indices],
            plan.weights,
            self.lineages[plan.indices],
        )
        return plan

    def compiled_storage(
        self,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
        return (
            self._positions,
            self._velocities,
            self._times,
            self._weights,
            self._lineages,
            self._size,
        )

    def adopt_compiled_storage(
        self,
        positions: np.ndarray,
        velocities: np.ndarray,
        times: np.ndarray,
        weights: np.ndarray,
        lineages: np.ndarray,
        size: int,
    ) -> None:
        active = int(size)
        capacity = int(weights.size)
        if (
            active < 0
            or active > capacity
            or positions.shape != (capacity, 3)
            or velocities.shape != (capacity, 3)
            or times.shape != (capacity,)
            or lineages.shape != (capacity,)
        ):
            raise ValueError("compiled particle storage has inconsistent dimensions")
        self._positions = positions
        self._velocities = velocities
        self._times = times
        self._weights = weights
        self._lineages = lineages
        self._size = active
