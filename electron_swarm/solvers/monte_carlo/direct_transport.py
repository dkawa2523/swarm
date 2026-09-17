"""Fixed-population direct transport estimators for internal Monte Carlo."""

from __future__ import annotations

from dataclasses import dataclass, fields

import numpy as np

from electron_swarm.core.constants import ELECTRON_MASS_KG, EV_TO_J
import electron_swarm.solvers.monte_carlo.evidence as _mc_evidence
import electron_swarm.solvers.monte_carlo.transport_common as _common


def direct_transport_observation_times(
    *,
    production_trial_steps: int,
    trial_frequency_s_inv: float,
) -> np.ndarray:
    """Return a fixed-cadence transport schedule.

    Cadence and every correlation lag are independent of production length.
    A longer production run adds time origins; it can never lengthen the lag.
    Very short smoke runs retain one deliberately unqualified observation.
    """

    steps = int(production_trial_steps)
    frequency = float(trial_frequency_s_inv)
    if steps <= 0 or not np.isfinite(frequency) or frequency <= 0.0:
        raise ValueError(
            "MC transport schedule requires positive production and frequency"
        )
    interval = _mc_evidence.DIRECT_MC_TRANSPORT_CADENCE_TRIAL_PERIODS / frequency
    count = int(
        np.floor(
            _mc_evidence.DIRECT_MC_TRANSPORT_OBSERVATION_END_FRACTION
            * steps
            / _mc_evidence.DIRECT_MC_TRANSPORT_CADENCE_TRIAL_PERIODS
        )
    )
    if count < 2:
        return np.asarray(
            [0.05 * steps / frequency],
            dtype=float,
        )
    return interval * np.arange(1, count + 1, dtype=float)


def _validate_snapshot_arrays(
    *,
    positions_m: np.ndarray,
    velocities_m_s: np.ndarray,
    ages_s: np.ndarray,
    weights: np.ndarray,
    electric_field_V_m: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    positions = np.asarray(positions_m, dtype=float)
    velocities = np.asarray(velocities_m_s, dtype=float)
    ages = np.asarray(ages_s, dtype=float)
    sample_weights = np.asarray(weights, dtype=float)
    field = np.asarray(electric_field_V_m, dtype=float)
    if (
        positions.ndim != 2
        or positions.shape[1] != 3
        or velocities.shape != positions.shape
        or ages.shape != (positions.shape[0],)
        or sample_weights.shape != ages.shape
        or field.shape != (3,)
        or positions.shape[0] < 2
    ):
        raise ValueError("direct MC transport arrays have incompatible shapes")
    if (
        np.any(~np.isfinite(positions))
        or np.any(~np.isfinite(velocities))
        or np.any(~np.isfinite(ages))
        or np.any(ages <= 0.0)
        or np.any(~np.isfinite(sample_weights))
        or np.any(sample_weights < 0.0)
        or np.any(~np.isfinite(field))
        or float(np.sum(sample_weights)) <= 0.0
        or float(np.dot(field, field)) <= 0.0
    ):
        raise ValueError("direct MC transport requires finite positive sampling data")
    return positions, velocities, ages, sample_weights, field


def direct_flux_transport_snapshot(
    *,
    positions_m: np.ndarray,
    velocities_m_s: np.ndarray,
    ages_s: np.ndarray,
    weights: np.ndarray,
    electric_field_V_m: np.ndarray,
) -> _common.DirectFluxTransport:
    """Evaluate one direct fixed-population trajectory block.

    ``positions_m`` are displacements over ``ages_s``.  Particle diffusion is
    the centered displacement cumulant divided by twice the lag.  The energy
    diffusion moment intentionally remains a terminal energy-current
    correlation over the same fixed short lag.
    """

    positions, velocities, ages, sample_weights, field = _validate_snapshot_arrays(
        positions_m=positions_m,
        velocities_m_s=velocities_m_s,
        ages_s=ages_s,
        weights=weights,
        electric_field_V_m=electric_field_V_m,
    )
    normalized_weights = sample_weights / float(np.sum(sample_weights))
    mean_age = float(np.sum(normalized_weights * ages))
    mean_velocity = np.sum(
        sample_weights[:, None] * positions,
        axis=0,
    ) / float(np.sum(sample_weights * ages))
    drift_frame_positions = positions - ages[:, None] * mean_velocity[None, :]
    centered_positions = drift_frame_positions - _common._weighted_mean(
        drift_frame_positions,
        normalized_weights,
    )

    field_squared = float(np.dot(field, field))
    field_unit = field / np.sqrt(field_squared)
    displacement_covariance = np.einsum(
        "n,ni,nj->ij",
        normalized_weights,
        centered_positions,
        centered_positions,
    )
    particle_diffusion = displacement_covariance / (2.0 * mean_age)

    energy_eV = (
        0.5 * ELECTRON_MASS_KG * np.sum(velocities * velocities, axis=1) / EV_TO_J
    )
    mean_energy = float(np.sum(normalized_weights * energy_eV))
    if not np.isfinite(mean_energy) or mean_energy <= 0.0:
        raise ValueError("direct MC transport requires positive mean energy")
    energy_flux_velocity = (
        _common._weighted_mean(
            energy_eV[:, None] * velocities,
            normalized_weights,
        )
        / mean_energy
    )
    position_energy = _common._weighted_mean(
        centered_positions * energy_eV[:, None],
        normalized_weights,
    )
    position_energy_flux = np.einsum(
        "n,ni,n,nj->ij",
        normalized_weights,
        centered_positions,
        energy_eV,
        velocities,
    )
    energy_diffusion = (
        position_energy_flux - np.outer(position_energy, energy_flux_velocity)
    ) / mean_energy

    particle_longitudinal = float(field_unit @ particle_diffusion @ field_unit)
    particle_transverse = float(
        (np.trace(particle_diffusion) - particle_longitudinal) / 2.0
    )
    energy_longitudinal = float(field_unit @ energy_diffusion @ field_unit)
    energy_transverse = float((np.trace(energy_diffusion) - energy_longitudinal) / 2.0)
    drift_velocity = -float(np.dot(mean_velocity, field_unit))
    mobility = -float(np.dot(mean_velocity, field)) / field_squared
    energy_mobility = -float(np.dot(energy_flux_velocity, field)) / field_squared
    return _common.DirectFluxTransport(
        drift_velocity_m_s=drift_velocity,
        mobility_m2_V_s=mobility,
        diffusion_L_m2_s=particle_longitudinal,
        diffusion_T_m2_s=particle_transverse,
        energy_mobility_m2_V_s=energy_mobility,
        energy_diffusion_L_m2_s=energy_longitudinal,
        energy_diffusion_T_m2_s=energy_transverse,
        mean_energy_eV=mean_energy,
    )


@dataclass(slots=True)
class DirectFluxTransportAccumulator:
    """Streaming mean of correlated origins inside one independent replica."""

    _count: int
    _sums: np.ndarray

    def __init__(self) -> None:
        self._count = 0
        self._sums = np.zeros(len(fields(_common.DirectFluxTransport)), dtype=float)

    def record(
        self, snapshot: _common.DirectFluxTransport | None = None, **arrays: np.ndarray
    ) -> None:
        value = snapshot or direct_flux_transport_snapshot(**arrays)
        self._sums += np.asarray(
            [getattr(value, item.name) for item in fields(_common.DirectFluxTransport)],
            dtype=float,
        )
        self._count += 1

    @property
    def samples(self) -> int:
        return self._count

    def mean(self) -> _common.DirectFluxTransport:
        if self._count <= 0:
            raise RuntimeError("direct MC transport has no production blocks")
        return _common.DirectFluxTransport(
            **{
                item.name: float(self._sums[index] / self._count)
                for index, item in enumerate(fields(_common.DirectFluxTransport))
            }
        )


@dataclass(slots=True)
class _PendingPlane:
    positions: np.ndarray
    velocities: np.ndarray
    weights: np.ndarray
    filled: np.ndarray


@dataclass(frozen=True, slots=True)
class _CompletePlane:
    positions: np.ndarray
    velocities: np.ndarray
    weights: np.ndarray


class SynchronizedFluxTransportObserver:
    """Fixed-cadence rolling time origins with bounded memory."""

    def __init__(
        self,
        *,
        observation_times_s: np.ndarray,
        particles: int,
        electric_field_V_m: np.ndarray,
        cadence_trial_periods: int = (
            _mc_evidence.DIRECT_MC_TRANSPORT_CADENCE_TRIAL_PERIODS
        ),
    ) -> None:
        times = np.asarray(observation_times_s, dtype=float)
        if (
            times.ndim != 1
            or len(times) == 0
            or np.any(~np.isfinite(times))
            or np.any(times <= 0.0)
            or np.any(np.diff(times) <= 0.0)
        ):
            raise ValueError("MC transport observation times must increase")
        if len(times) > 1 and not np.allclose(
            np.diff(times),
            np.diff(times)[0],
            rtol=1.0e-12,
            atol=0.0,
        ):
            raise ValueError("MC transport observations require a fixed cadence")
        particle_count = int(particles)
        cadence = int(cadence_trial_periods)
        if particle_count < 2 or cadence <= 0:
            raise ValueError("MC transport observation settings are invalid")
        field = np.asarray(electric_field_V_m, dtype=float)
        if (
            field.shape != (3,)
            or np.any(~np.isfinite(field))
            or np.dot(field, field) <= 0.0
        ):
            raise ValueError("MC transport observation requires a positive field")

        self.observation_times_s = times
        self.electric_field_V_m = field
        self.cadence_trial_periods = cadence
        self._particles = particle_count
        self._next_index = np.zeros(particle_count, dtype=int)
        self._pending: dict[int, _PendingPlane] = {}
        self._ready: dict[int, _CompletePlane] = {}
        self._planes: dict[int, _CompletePlane] = {}
        self._next_complete_index = 0
        self._lag_accumulators = {
            lag: DirectFluxTransportAccumulator()
            for lag in _mc_evidence.DIRECT_MC_TRANSPORT_LAG_PLANES
        }
        self._selected_early = DirectFluxTransportAccumulator()
        self._selected_late = DirectFluxTransportAccumulator()

        self._weighted_time_s = 0.0
        self._weighted_displacement_m = np.zeros(3, dtype=float)
        self._weighted_energy_time_eV_s = 0.0
        self._weighted_energy_displacement_eV_m = np.zeros(3, dtype=float)

    def cap_step(
        self, particle_index: int, current_time_s: float, dt_s: float
    ) -> float:
        particle = int(particle_index)
        index = int(self._next_index[particle])
        if index >= len(self.observation_times_s):
            return float(dt_s)
        remaining = float(self.observation_times_s[index] - current_time_s)
        tolerance = max(1.0e-13 * self.observation_times_s[index], 1.0e-300)
        if remaining <= tolerance:
            return 0.0
        return min(float(dt_s), remaining)

    def record_residence(
        self,
        *,
        displacement_m: np.ndarray,
        sample_energy_eV: float,
        dt_s: float,
        weight: float,
    ) -> None:
        displacement = np.asarray(displacement_m, dtype=float)
        energy = float(sample_energy_eV)
        duration = float(dt_s)
        sample_weight = float(weight)
        if (
            displacement.shape != (3,)
            or np.any(~np.isfinite(displacement))
            or not np.isfinite(energy)
            or energy <= 0.0
            or not np.isfinite(duration)
            or duration <= 0.0
            or not np.isfinite(sample_weight)
            or sample_weight < 0.0
        ):
            raise ValueError("MC transport residence sample is invalid")
        self._weighted_time_s += sample_weight * duration
        self._weighted_displacement_m += sample_weight * displacement
        self._weighted_energy_time_eV_s += sample_weight * energy * duration
        self._weighted_energy_displacement_eV_m += sample_weight * energy * displacement

    def record_due(
        self,
        *,
        particle_index: int,
        time_s: float,
        position_m: np.ndarray,
        velocity_m_s: np.ndarray,
        weight: float,
    ) -> bool:
        particle = int(particle_index)
        index = int(self._next_index[particle])
        if index >= len(self.observation_times_s):
            return False
        target = float(self.observation_times_s[index])
        tolerance = max(1.0e-11 * target, 1.0e-300)
        if abs(float(time_s) - target) > tolerance:
            return False
        pending = self._pending.get(index)
        if pending is None:
            pending = _PendingPlane(
                positions=np.full((self._particles, 3), np.nan),
                velocities=np.full((self._particles, 3), np.nan),
                weights=np.full(self._particles, np.nan),
                filled=np.zeros(self._particles, dtype=bool),
            )
            self._pending[index] = pending
        pending.positions[particle] = np.asarray(position_m, dtype=float)
        pending.velocities[particle] = np.asarray(velocity_m_s, dtype=float)
        pending.weights[particle] = float(weight)
        pending.filled[particle] = True
        self._next_index[particle] += 1
        if bool(np.all(pending.filled)):
            self._ready[index] = _CompletePlane(
                positions=pending.positions,
                velocities=pending.velocities,
                weights=pending.weights,
            )
            del self._pending[index]
            self._consume_ready_planes()
        return True

    def _consume_ready_planes(self) -> None:
        while self._next_complete_index in self._ready:
            index = self._next_complete_index
            plane = self._ready.pop(index)
            self._planes[index] = plane
            for lag, accumulator in self._lag_accumulators.items():
                origin_index = index - lag
                origin = self._planes.get(origin_index)
                if origin is None:
                    continue
                lag_s = float(
                    self.observation_times_s[index]
                    - self.observation_times_s[origin_index]
                )
                snapshot = direct_flux_transport_snapshot(
                    positions_m=plane.positions - origin.positions,
                    velocities_m_s=plane.velocities,
                    ages_s=np.full(self._particles, lag_s),
                    weights=plane.weights,
                    electric_field_V_m=self.electric_field_V_m,
                )
                accumulator.record(snapshot)
                if lag == _mc_evidence.DIRECT_MC_TRANSPORT_LAG_PLANES[-1]:
                    expected_origins = max(
                        len(self.observation_times_s) - lag,
                        1,
                    )
                    origin_order = origin_index
                    target = (
                        self._selected_early
                        if origin_order < expected_origins // 2
                        else self._selected_late
                    )
                    target.record(snapshot)
            oldest = index - _mc_evidence.DIRECT_MC_TRANSPORT_LAG_PLANES[-1]
            if oldest in self._planes:
                del self._planes[oldest]
            self._next_complete_index += 1

    def _fallback_accumulator(self) -> tuple[int, DirectFluxTransportAccumulator]:
        complete = self._next_complete_index
        if complete <= 0:
            raise RuntimeError("direct MC transport has no complete time plane")
        target_index = complete - 1
        target = self._planes[target_index]
        if complete >= 2:
            origin_index = 0
            origin = self._planes[origin_index]
            displacement = target.positions - origin.positions
            lag_s = float(
                self.observation_times_s[target_index]
                - self.observation_times_s[origin_index]
            )
            lag_planes = target_index - origin_index
        else:
            displacement = target.positions
            lag_s = float(self.observation_times_s[target_index])
            lag_planes = 1
        accumulator = DirectFluxTransportAccumulator()
        accumulator.record(
            positions_m=displacement,
            velocities_m_s=target.velocities,
            ages_s=np.full(self._particles, lag_s),
            weights=target.weights,
            electric_field_V_m=self.electric_field_V_m,
        )
        return lag_planes, accumulator

    def finalize(self) -> tuple[_common.DirectFluxTransport, dict[str, object]]:
        if self._weighted_time_s <= 0.0 or self._weighted_energy_time_eV_s <= 0.0:
            raise RuntimeError("direct MC transport has no residence samples")
        available = [
            lag
            for lag in _mc_evidence.DIRECT_MC_TRANSPORT_LAG_PLANES
            if self._lag_accumulators[lag].samples > 0
        ]
        accumulators = self._lag_accumulators
        if not available:
            fallback_lag, fallback = self._fallback_accumulator()
            available = [fallback_lag]
            accumulators = {fallback_lag: fallback}

        block_means = [accumulators[lag].mean() for lag in available]
        selected_block = block_means[-1]
        field_squared = float(np.dot(self.electric_field_V_m, self.electric_field_V_m))
        field_unit = self.electric_field_V_m / np.sqrt(field_squared)
        mean_velocity = self._weighted_displacement_m / self._weighted_time_s
        mean_energy = self._weighted_energy_time_eV_s / self._weighted_time_s
        energy_flux_velocity = (
            self._weighted_energy_displacement_eV_m / self._weighted_energy_time_eV_s
        )
        drift_velocity = -float(np.dot(mean_velocity, field_unit))
        mobility = (
            -float(np.dot(mean_velocity, self.electric_field_V_m)) / field_squared
        )
        energy_mobility = (
            -float(np.dot(energy_flux_velocity, self.electric_field_V_m))
            / field_squared
        )
        estimate = _common.DirectFluxTransport(
            drift_velocity_m_s=drift_velocity,
            mobility_m2_V_s=mobility,
            diffusion_L_m2_s=selected_block.diffusion_L_m2_s,
            diffusion_T_m2_s=selected_block.diffusion_T_m2_s,
            energy_mobility_m2_V_s=energy_mobility,
            energy_diffusion_L_m2_s=selected_block.energy_diffusion_L_m2_s,
            energy_diffusion_T_m2_s=selected_block.energy_diffusion_T_m2_s,
            mean_energy_eV=mean_energy,
        )

        interval = float(
            np.diff(self.observation_times_s)[0]
            if len(self.observation_times_s) > 1
            else self.observation_times_s[0]
        )
        lag_seconds = [float(lag * interval) for lag in available]
        lag_fields = (
            "diffusion_L_m2_s",
            "diffusion_T_m2_s",
            "energy_diffusion_L_m2_s",
            "energy_diffusion_T_m2_s",
        )
        estimates = {
            name: [float(getattr(value, name)) for value in block_means]
            for name in lag_fields
        }

        selected_mean = selected_block
        early = (
            self._selected_early.mean()
            if self._selected_early.samples
            else selected_mean
        )
        late = (
            self._selected_late.mean() if self._selected_late.samples else selected_mean
        )

        means = [
            getattr(estimate, item.name) for item in fields(_common.DirectFluxTransport)
        ]
        positive_names = (
            "mobility_m2_V_s",
            "diffusion_L_m2_s",
            "diffusion_T_m2_s",
            "energy_mobility_m2_V_s",
            "energy_diffusion_L_m2_s",
            "energy_diffusion_T_m2_s",
            "mean_energy_eV",
        )
        complete_observation = bool(
            len(self.observation_times_s)
            >= _mc_evidence.DIRECT_MC_TRANSPORT_OBSERVATION_PLANES
            and self._next_complete_index == len(self.observation_times_s)
            and available == list(_mc_evidence.DIRECT_MC_TRANSPORT_LAG_PLANES)
        )
        diagnostics: dict[str, object] = {
            "estimator": "direct_mc_fixed_lag_block_helfand_energy_flux",
            "estimator_schema_version": (
                _mc_evidence.DIRECT_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
            ),
            "finite": all(np.isfinite(value) for value in means),
            "strictly_positive": all(
                getattr(estimate, name) > 0.0 for name in positive_names
            ),
            "multiple_time_origins": bool(accumulators[available[-1]].samples >= 2),
            "complete_observation": complete_observation,
            "origin_stationarity_window_mean_early": {
                "mean_energy_eV": early.mean_energy_eV,
            },
            "origin_stationarity_window_mean_late": {
                "mean_energy_eV": late.mean_energy_eV,
            },
            "component_estimators": {
                "drift_mobility": (
                    "production_trajectory_displacement_over_residence_time"
                ),
                "particle_diffusion": "block_helfand_mean_square_displacement",
                "energy_mobility": "production_residence_energy_current",
                "energy_diffusion": ("restricted_density_packet_energy_flux_fixed_lag"),
            },
            "production_estimates": {
                "drift_velocity_m_s": drift_velocity,
                "mobility_m2_V_s": mobility,
                "energy_mobility_m2_V_s": energy_mobility,
                "mean_energy_eV": mean_energy,
            },
            "lag_scan": {
                "lag_planes": available,
                "lag_s": lag_seconds,
                "estimates": estimates,
                "completed_origins": [accumulators[lag].samples for lag in available],
                "selected_lag_index": len(available) - 1,
            },
            "transport_sampling": {
                "lag_grid_source": (
                    "fixed_physical_time_independent_of_production_length"
                ),
                "origin_cadence_s": interval,
                "origin_cadence_trial_periods": self.cadence_trial_periods,
                "rolling_origins": True,
                "production_length_controls_lag": False,
                "warmup_accumulators_reset": True,
                "observation_planes_requested": len(self.observation_times_s),
                "observation_planes_complete": self._next_complete_index,
                "minimum_observation_planes": (
                    _mc_evidence.DIRECT_MC_TRANSPORT_OBSERVATION_PLANES
                ),
            },
            "energy_transport_from_eedf": False,
            "two_term_closure_used": False,
            "energy_diffusion_interpretation": (
                "restricted_direct_density_packet_energy_flux_correlation"
            ),
            "standard_energy_density_gradient_closure_available": True,
            "cross_gradient_response_identified": False,
        }
        diagnostics["energy_transport_status"] = (
            "direct_mc_candidate_for_ensemble_qualification"
            if diagnostics["finite"]
            and diagnostics["strictly_positive"]
            and diagnostics["multiple_time_origins"]
            and diagnostics["complete_observation"]
            else "direct_mc_unqualified"
        )
        return estimate, diagnostics
