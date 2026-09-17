"""Weighted-growth transport estimators for internal Monte Carlo."""

from __future__ import annotations

from dataclasses import dataclass, fields

import numpy as np

from electron_swarm.core.constants import ELECTRON_MASS_KG, EV_TO_J
import electron_swarm.solvers.monte_carlo.evidence as _mc_evidence
import electron_swarm.solvers.monte_carlo.transport_common as _common


# Weighted branching is a different stochastic process from tagged-particle
# tracking.  Its flux coefficients are sampled from finite displacement blocks
# at synchronized physical-time barriers; never mix this schema with the
# fixed-population Helfand estimator.
@dataclass(frozen=True, slots=True)
class WeightedGrowthTransportLagPlan:
    correlation_lag_barriers: int
    lag_barriers: tuple[int, int, int, int]
    block_barriers: int
    production_lag_barriers: int
    hard_convergence_pair_barriers: tuple[int, int]
    lineage_qualification_lag_barriers: tuple[int, int]
    supplemental_pair_barriers: tuple[int, int]


def weighted_growth_transport_lag_plan(
    correlation_lag_barriers: int = 64,
) -> WeightedGrowthTransportLagPlan:
    """Build the complete weighted-growth correlation plan from one lag."""

    if isinstance(correlation_lag_barriers, bool) or not isinstance(
        correlation_lag_barriers, (int, np.integer)
    ):
        raise ValueError("weighted-growth correlation lag must be an integer")
    lag = int(correlation_lag_barriers)
    if lag < 4 or lag & (lag - 1):
        raise ValueError(
            "weighted-growth correlation lag must be a positive power of two >= 4"
        )
    hard_pair = (lag // 2, lag)
    return WeightedGrowthTransportLagPlan(
        correlation_lag_barriers=lag,
        lag_barriers=(lag // 4, lag // 2, lag, 2 * lag),
        block_barriers=2 * lag,
        production_lag_barriers=lag,
        hard_convergence_pair_barriers=hard_pair,
        lineage_qualification_lag_barriers=hard_pair,
        supplemental_pair_barriers=(lag, 2 * lag),
    )


def weighted_growth_flux_transport_snapshot(
    *,
    positions_m: np.ndarray,
    velocities_m_s: np.ndarray,
    weights: np.ndarray,
    electric_field_V_m: np.ndarray,
) -> _common.DirectFluxTransport:
    """Evaluate flux moments of a growing swarm at one physical time.

    For nonconservative swarms the endpoint displacement variance is a bulk
    coefficient.  The flux diffusion tensor instead follows directly from
    the synchronized covariance ``<delta r_i delta v_j>``.  The energy
    coefficient is the analogous one-sided density-packet energy-current
    covariance; it does not identify the independent energy-gradient
    response.
    """

    positions = np.asarray(positions_m, dtype=float)
    velocities = np.asarray(velocities_m_s, dtype=float)
    sample_weights = np.asarray(weights, dtype=float)
    field = np.asarray(electric_field_V_m, dtype=float)
    if (
        positions.ndim != 2
        or positions.shape[1] != 3
        or velocities.shape != positions.shape
        or sample_weights.shape != (positions.shape[0],)
        or positions.shape[0] < 2
        or field.shape != (3,)
        or np.any(~np.isfinite(positions))
        or np.any(~np.isfinite(velocities))
        or np.any(~np.isfinite(sample_weights))
        or np.any(sample_weights < 0.0)
        or not np.isfinite(float(np.sum(sample_weights)))
        or float(np.sum(sample_weights)) <= 0.0
        or np.any(~np.isfinite(field))
        or float(np.dot(field, field)) <= 0.0
    ):
        raise ValueError("weighted-growth MC transport plane is invalid")

    normalized = sample_weights / float(np.sum(sample_weights))
    mean_position = _common._weighted_mean(positions, normalized)
    mean_velocity = _common._weighted_mean(velocities, normalized)
    centered_position = positions - mean_position
    centered_velocity = velocities - mean_velocity
    particle_diffusion = np.einsum(
        "n,ni,nj->ij",
        normalized,
        centered_position,
        centered_velocity,
    )

    energy_eV = (
        0.5 * ELECTRON_MASS_KG * np.sum(velocities * velocities, axis=1) / EV_TO_J
    )
    mean_energy = float(np.sum(normalized * energy_eV))
    if not np.isfinite(mean_energy) or mean_energy <= 0.0:
        raise ValueError("weighted-growth MC transport requires positive energy")
    energy_flux_velocity = (
        _common._weighted_mean(
            energy_eV[:, None] * velocities,
            normalized,
        )
        / mean_energy
    )
    position_energy = _common._weighted_mean(
        centered_position * energy_eV[:, None],
        normalized,
    )
    position_energy_flux = np.einsum(
        "n,ni,n,nj->ij",
        normalized,
        centered_position,
        energy_eV,
        velocities,
    )
    energy_diffusion = (
        position_energy_flux - np.outer(position_energy, energy_flux_velocity)
    ) / mean_energy

    field_squared = float(np.dot(field, field))
    field_unit = field / np.sqrt(field_squared)
    particle_longitudinal = float(field_unit @ particle_diffusion @ field_unit)
    particle_transverse = float(
        (np.trace(particle_diffusion) - particle_longitudinal) / 2.0
    )
    energy_longitudinal = float(field_unit @ energy_diffusion @ field_unit)
    energy_transverse = float((np.trace(energy_diffusion) - energy_longitudinal) / 2.0)
    return _common.DirectFluxTransport(
        drift_velocity_m_s=-float(np.dot(mean_velocity, field_unit)),
        mobility_m2_V_s=-float(np.dot(mean_velocity, field)) / field_squared,
        diffusion_L_m2_s=particle_longitudinal,
        diffusion_T_m2_s=particle_transverse,
        energy_mobility_m2_V_s=(
            -float(np.dot(energy_flux_velocity, field)) / field_squared
        ),
        energy_diffusion_L_m2_s=energy_longitudinal,
        energy_diffusion_T_m2_s=energy_transverse,
        mean_energy_eV=mean_energy,
    )


def _effective_lineage_count(lineages: np.ndarray, weights: np.ndarray) -> float:
    ids = np.asarray(lineages)
    sample_weights = np.asarray(weights, dtype=float)
    if ids.shape != sample_weights.shape or ids.ndim != 1:
        raise ValueError("weighted-growth lineage ids must match weights")
    _, inverse = np.unique(ids, return_inverse=True)
    family_weights = np.bincount(inverse, weights=sample_weights)
    total = float(np.sum(family_weights))
    if not np.isfinite(total) or total <= 0.0:
        raise ValueError("weighted-growth lineage weights are invalid")
    fractions = family_weights / total
    return float(1.0 / np.sum(fractions * fractions))


class SynchronizedWeightedGrowthFluxObserver:
    """Block-lag flux sampler for a population-controlled branching swarm.

    A homogeneous swarm is translation invariant.  The caller therefore
    resets particle positions after every transport block and this observer
    samples several finite displacement lags inside each block.  This avoids
    treating a long absolute-position history as if additional time planes
    were independent transport origins.
    """

    def __init__(
        self,
        *,
        particles: int,
        electric_field_V_m: np.ndarray,
        barrier_cadence_s: float,
        lag_plan: WeightedGrowthTransportLagPlan | None = None,
    ) -> None:
        particle_count = int(particles)
        field = np.asarray(electric_field_V_m, dtype=float)
        cadence = float(barrier_cadence_s)
        if (
            particle_count < 2
            or field.shape != (3,)
            or np.any(~np.isfinite(field))
            or float(np.dot(field, field)) <= 0.0
            or not np.isfinite(cadence)
            or cadence <= 0.0
        ):
            raise ValueError("weighted-growth MC observer settings are invalid")
        self.particles = particle_count
        self.electric_field_V_m = field
        self.barrier_cadence_s = cadence
        self.lag_plan = lag_plan or weighted_growth_transport_lag_plan()
        self._lag_planes: dict[int, list[_common.DirectFluxTransport]] = {}
        self._plane_times_s: list[float] = []
        self._effective_lineages: dict[int, list[float]] = {}
        self._weighted_time_s = 0.0
        self._weighted_displacement_m = np.zeros(3, dtype=float)
        self._weighted_energy_time_eV_s = 0.0
        self._weighted_energy_displacement_eV_m = np.zeros(3, dtype=float)
        self._block_weighted_time_s = 0.0
        self._block_weighted_displacement_m = np.zeros(3, dtype=float)
        self._block_weighted_energy_time_eV_s = 0.0
        self._block_weighted_energy_displacement_eV_m = np.zeros(3, dtype=float)
        self._completed_residence_blocks: list[
            tuple[float, np.ndarray, float, np.ndarray]
        ] = []

    def cap_step(
        self, particle_index: int, current_time_s: float, dt_s: float
    ) -> float:
        del particle_index, current_time_s
        return float(dt_s)

    def record_due(self, **_: object) -> bool:
        return False

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
            or sample_weight <= 0.0
        ):
            raise ValueError("weighted-growth residence sample is invalid")
        self._weighted_time_s += sample_weight * duration
        self._weighted_displacement_m += sample_weight * displacement
        self._weighted_energy_time_eV_s += sample_weight * energy * duration
        self._weighted_energy_displacement_eV_m += sample_weight * energy * displacement
        self._block_weighted_time_s += sample_weight * duration
        self._block_weighted_displacement_m += sample_weight * displacement
        self._block_weighted_energy_time_eV_s += sample_weight * energy * duration
        self._block_weighted_energy_displacement_eV_m += (
            sample_weight * energy * displacement
        )

    def record_aggregate_residence(
        self,
        *,
        weighted_time_s: float,
        weighted_displacement_m: np.ndarray,
        weighted_energy_time_eV_s: float,
        weighted_energy_displacement_eV_m: np.ndarray,
    ) -> None:
        """Commit sufficient statistics accumulated by one numerical barrier."""

        weighted_time = float(weighted_time_s)
        displacement = np.asarray(weighted_displacement_m, dtype=float)
        weighted_energy_time = float(weighted_energy_time_eV_s)
        energy_displacement = np.asarray(
            weighted_energy_displacement_eV_m,
            dtype=float,
        )
        if (
            not np.isfinite(weighted_time)
            or weighted_time <= 0.0
            or displacement.shape != (3,)
            or np.any(~np.isfinite(displacement))
            or not np.isfinite(weighted_energy_time)
            or weighted_energy_time <= 0.0
            or energy_displacement.shape != (3,)
            or np.any(~np.isfinite(energy_displacement))
        ):
            raise ValueError("weighted-growth aggregate residence sample is invalid")
        self._weighted_time_s += weighted_time
        self._weighted_displacement_m += displacement
        self._weighted_energy_time_eV_s += weighted_energy_time
        self._weighted_energy_displacement_eV_m += energy_displacement
        self._block_weighted_time_s += weighted_time
        self._block_weighted_displacement_m += displacement
        self._block_weighted_energy_time_eV_s += weighted_energy_time
        self._block_weighted_energy_displacement_eV_m += energy_displacement

    def record_plane(
        self,
        *,
        time_s: float,
        positions_m: np.ndarray,
        velocities_m_s: np.ndarray,
        weights: np.ndarray,
        lineages: np.ndarray,
        lag_barriers: int,
    ) -> None:
        time_value = float(time_s)
        lag = int(lag_barriers)
        if (
            not np.isfinite(time_value)
            or time_value <= 0.0
            or (self._plane_times_s and time_value <= self._plane_times_s[-1])
            or lag <= 0
            or lag > self.lag_plan.block_barriers
        ):
            raise ValueError("weighted-growth block-lag observation is invalid")
        snapshot = weighted_growth_flux_transport_snapshot(
            positions_m=positions_m,
            velocities_m_s=velocities_m_s,
            weights=weights,
            electric_field_V_m=self.electric_field_V_m,
        )
        self._lag_planes.setdefault(lag, []).append(snapshot)
        self._plane_times_s.append(time_value)
        self._effective_lineages.setdefault(lag, []).append(
            _effective_lineage_count(
                lineages,
                np.asarray(weights, dtype=float),
            )
        )
        if lag == self.lag_plan.block_barriers:
            if (
                self._block_weighted_time_s <= 0.0
                or self._block_weighted_energy_time_eV_s <= 0.0
            ):
                raise RuntimeError(
                    "weighted-growth transport block has no residence evidence"
                )
            self._completed_residence_blocks.append(
                (
                    self._block_weighted_time_s,
                    self._block_weighted_displacement_m.copy(),
                    self._block_weighted_energy_time_eV_s,
                    self._block_weighted_energy_displacement_eV_m.copy(),
                )
            )
            self._block_weighted_time_s = 0.0
            self._block_weighted_displacement_m.fill(0.0)
            self._block_weighted_energy_time_eV_s = 0.0
            self._block_weighted_energy_displacement_eV_m.fill(0.0)

    def finalize(self) -> tuple[_common.DirectFluxTransport, dict[str, object]]:
        if (
            self._weighted_time_s <= 0.0
            or self._weighted_energy_time_eV_s <= 0.0
            or not self._lag_planes
        ):
            raise RuntimeError(
                "weighted-growth MC transport has no production evidence"
            )

        def mean_fields(
            values: list[_common.DirectFluxTransport],
        ) -> _common.DirectFluxTransport:
            return _common.DirectFluxTransport(
                **{
                    item.name: float(
                        np.mean([getattr(value, item.name) for value in values])
                    )
                    for item in fields(_common.DirectFluxTransport)
                }
            )

        complete_blocks = len(self._completed_residence_blocks)
        canonical_lags_complete = all(
            len(self._lag_planes.get(lag, ())) >= complete_blocks
            for lag in self.lag_plan.lag_barriers
        )
        complete = bool(
            complete_blocks >= _mc_evidence.WEIGHTED_MC_TRANSPORT_MIN_BLOCKS
            and canonical_lags_complete
        )
        production_lag = self.lag_plan.production_lag_barriers
        production_lag_available = production_lag in self._lag_planes
        selected_lag = (
            production_lag if production_lag_available else max(self._lag_planes)
        )
        available_lags = (
            list(self.lag_plan.lag_barriers)
            if complete_blocks > 0 and canonical_lags_complete
            else sorted(self._lag_planes)
        )
        if complete_blocks > 0 and production_lag_available:
            selected_planes = self._lag_planes[production_lag][:complete_blocks]
            selected_residence_blocks = self._completed_residence_blocks
        else:
            # Short smoke runs may report a partial-lag estimate, but remain
            # explicitly ineligible for ensemble qualification.
            selected_planes = self._lag_planes[selected_lag]
            selected_residence_blocks = [
                (
                    self._block_weighted_time_s,
                    self._block_weighted_displacement_m.copy(),
                    self._block_weighted_energy_time_eV_s,
                    self._block_weighted_energy_displacement_eV_m.copy(),
                )
            ]
        half = max(len(selected_planes) // 2, 1)
        early_planes = selected_planes[:half]
        late_planes = selected_planes[-half:]
        early_residence_blocks = selected_residence_blocks[:half]
        late_residence_blocks = selected_residence_blocks[-half:]
        plane_mean = mean_fields(selected_planes)
        field_squared = float(np.dot(self.electric_field_V_m, self.electric_field_V_m))
        field_unit = self.electric_field_V_m / np.sqrt(field_squared)

        def residence_ratio_estimate(
            residence_blocks: list[tuple[float, np.ndarray, float, np.ndarray]],
            diffusion_snapshot: _common.DirectFluxTransport,
        ) -> _common.DirectFluxTransport:
            weighted_time = float(sum(item[0] for item in residence_blocks))
            weighted_displacement = np.sum(
                [item[1] for item in residence_blocks], axis=0
            )
            weighted_energy_time = float(sum(item[2] for item in residence_blocks))
            weighted_energy_displacement = np.sum(
                [item[3] for item in residence_blocks], axis=0
            )
            if weighted_time <= 0.0 or weighted_energy_time <= 0.0:
                raise RuntimeError(
                    "weighted-growth stationarity window has no residence evidence"
                )
            mean_velocity = weighted_displacement / weighted_time
            energy_flux_velocity = weighted_energy_displacement / weighted_energy_time
            return _common.DirectFluxTransport(
                drift_velocity_m_s=-float(np.dot(mean_velocity, field_unit)),
                mobility_m2_V_s=(
                    -float(np.dot(mean_velocity, self.electric_field_V_m))
                    / field_squared
                ),
                diffusion_L_m2_s=diffusion_snapshot.diffusion_L_m2_s,
                diffusion_T_m2_s=diffusion_snapshot.diffusion_T_m2_s,
                energy_mobility_m2_V_s=(
                    -float(np.dot(energy_flux_velocity, self.electric_field_V_m))
                    / field_squared
                ),
                energy_diffusion_L_m2_s=(diffusion_snapshot.energy_diffusion_L_m2_s),
                energy_diffusion_T_m2_s=(diffusion_snapshot.energy_diffusion_T_m2_s),
                mean_energy_eV=weighted_energy_time / weighted_time,
            )

        early = residence_ratio_estimate(
            early_residence_blocks,
            mean_fields(early_planes),
        )
        late = residence_ratio_estimate(
            late_residence_blocks,
            mean_fields(late_planes),
        )
        mean_velocity = self._weighted_displacement_m / self._weighted_time_s
        mean_energy = self._weighted_energy_time_eV_s / self._weighted_time_s
        energy_flux_velocity = (
            self._weighted_energy_displacement_eV_m / self._weighted_energy_time_eV_s
        )
        estimate = _common.DirectFluxTransport(
            drift_velocity_m_s=-float(np.dot(mean_velocity, field_unit)),
            mobility_m2_V_s=(
                -float(np.dot(mean_velocity, self.electric_field_V_m)) / field_squared
            ),
            diffusion_L_m2_s=plane_mean.diffusion_L_m2_s,
            diffusion_T_m2_s=plane_mean.diffusion_T_m2_s,
            energy_mobility_m2_V_s=(
                -float(np.dot(energy_flux_velocity, self.electric_field_V_m))
                / field_squared
            ),
            energy_diffusion_L_m2_s=plane_mean.energy_diffusion_L_m2_s,
            energy_diffusion_T_m2_s=plane_mean.energy_diffusion_T_m2_s,
            mean_energy_eV=mean_energy,
        )
        numeric = [
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
        lineage_minimum_counts = []
        for lag in available_lags:
            lineage_values = self._effective_lineages.get(lag, [])
            if complete_blocks > 0:
                lineage_values = lineage_values[:complete_blocks]
            lineage_minimum_counts.append(
                float(min(lineage_values)) if lineage_values else 0.0
            )
        lineage_minimum_fractions = [
            value / self.particles for value in lineage_minimum_counts
        ]
        lineage_count_by_lag = dict(
            zip(
                available_lags,
                lineage_minimum_counts,
                strict=True,
            )
        )
        lineage_ok = all(
            lag in lineage_count_by_lag
            and lineage_count_by_lag[lag]
            >= _mc_evidence.WEIGHTED_MC_TRANSPORT_MIN_EFFECTIVE_LINEAGE_COUNT
            and lineage_count_by_lag[lag] / self.particles
            >= _mc_evidence.WEIGHTED_MC_TRANSPORT_MIN_EFFECTIVE_LINEAGE_FRACTION
            for lag in self.lag_plan.lineage_qualification_lag_barriers
        )
        status_ready = bool(
            all(np.isfinite(value) for value in numeric)
            and all(getattr(estimate, name) > 0.0 for name in positive_names)
            and len(selected_planes) >= 2
            and complete
            and production_lag_available
            and lineage_ok
        )

        def lag_planes_for_evidence(lag: int) -> list[_common.DirectFluxTransport]:
            values = self._lag_planes[lag]
            return values[:complete_blocks] if complete_blocks > 0 else values

        stationarity_fields = (
            "mobility_m2_V_s",
            "diffusion_L_m2_s",
            "diffusion_T_m2_s",
            "energy_mobility_m2_V_s",
            "energy_diffusion_L_m2_s",
            "energy_diffusion_T_m2_s",
            "mean_energy_eV",
        )
        diagnostics: dict[str, object] = {
            "estimator": "direct_mc_weighted_growth_block_lag_flux",
            "estimator_schema_version": (
                _mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
            ),
            "finite": all(np.isfinite(value) for value in numeric),
            "strictly_positive": all(
                getattr(estimate, name) > 0.0 for name in positive_names
            ),
            "multiple_time_blocks": len(selected_planes) >= 2,
            "complete_observation": complete,
            "lineage_qualified": lineage_ok,
            "time_stationarity_window_estimate_early": {
                name: float(getattr(early, name)) for name in stationarity_fields
            },
            "time_stationarity_window_estimate_late": {
                name: float(getattr(late, name)) for name in stationarity_fields
            },
            "component_estimators": {
                "drift_mobility": "normalized_weight_residence_velocity_moment",
                "particle_diffusion": (
                    "synchronized_block_lag_position_velocity_flux_covariance"
                ),
                "energy_mobility": ("normalized_weight_residence_energy_current"),
                "energy_diffusion": (
                    "synchronized_block_lag_restricted_density_packet_"
                    "energy_current_covariance"
                ),
            },
            "production_estimates": {
                item.name: float(getattr(estimate, item.name))
                for item in fields(_common.DirectFluxTransport)
            },
            "block_lag_sampling": {
                "physical_time_barriers": True,
                "resampling_at_common_time_only": True,
                "barrier_cadence_s": self.barrier_cadence_s,
                "configured_correlation_lag_barriers": (
                    self.lag_plan.correlation_lag_barriers
                ),
                "transport_block_barriers": self.lag_plan.block_barriers,
                "lag_barriers": list(self.lag_plan.lag_barriers),
                "complete_blocks": complete_blocks,
                "minimum_blocks": _mc_evidence.WEIGHTED_MC_TRANSPORT_MIN_BLOCKS,
                "positions_reset_each_block": True,
                "warmup_accumulators_reset": True,
            },
            "stationarity_sampling": {
                "window_definition": (
                    "first_half_vs_second_half_complete_transport_blocks"
                    if complete_blocks > 0 and production_lag_available
                    else "single_partial_transport_block_unqualified_smoke"
                ),
                "state_moment_aggregation": (
                    "sum_raw_weighted_residence_moments_before_ratio"
                ),
                "state_estimator": ("same_residence_ratio_formulas_as_production"),
                "diffusion_estimator": (
                    "mean_fixed_production_lag_snapshot_covariance_per_window"
                ),
                "diffusion_lag_barriers": production_lag,
                "early_block_count": len(early_residence_blocks),
                "late_block_count": len(late_residence_blocks),
                "odd_center_block_excluded": bool(
                    len(selected_residence_blocks) > 1
                    and len(selected_residence_blocks) % 2
                ),
            },
            "lag_scan": {
                "lag_barriers": available_lags,
                "lag_s": [
                    float(lag * self.barrier_cadence_s) for lag in available_lags
                ],
                "estimates": {
                    name: [
                        float(
                            getattr(
                                mean_fields(lag_planes_for_evidence(lag)),
                                name,
                            )
                        )
                        for lag in available_lags
                    ]
                    for name in (
                        "diffusion_L_m2_s",
                        "diffusion_T_m2_s",
                        "energy_diffusion_L_m2_s",
                        "energy_diffusion_T_m2_s",
                    )
                },
                "completed_blocks": [
                    len(lag_planes_for_evidence(lag)) for lag in available_lags
                ],
                "production_lag_barriers": production_lag,
                "hard_convergence_pair_barriers": list(
                    self.lag_plan.hard_convergence_pair_barriers
                ),
                "supplemental_pair_barriers": list(
                    self.lag_plan.supplemental_pair_barriers
                ),
            },
            "lineage_sampling": {
                "lineage_horizon": "each_sampled_lag",
                "lag_barriers": available_lags,
                "minimum_effective_lineage_count": lineage_minimum_counts,
                "minimum_effective_lineage_fraction": (lineage_minimum_fractions),
                "production_lag_barriers": production_lag,
                "qualification_lag_barriers": list(
                    self.lag_plan.lineage_qualification_lag_barriers
                ),
                "required_effective_lineage_count": (
                    _mc_evidence.WEIGHTED_MC_TRANSPORT_MIN_EFFECTIVE_LINEAGE_COUNT
                ),
                "required_effective_lineage_fraction": (
                    _mc_evidence.WEIGHTED_MC_TRANSPORT_MIN_EFFECTIVE_LINEAGE_FRACTION
                ),
            },
            "transport_definition": "flux",
            "bulk_transport_identified": False,
            "energy_transport_from_eedf": False,
            "two_term_closure_used": False,
            "energy_diffusion_interpretation": (
                "restricted_direct_density_packet_energy_flux_block_lag_correlation"
            ),
            "standard_energy_density_gradient_closure_available": True,
            "cross_gradient_response_identified": False,
        }
        diagnostics["energy_transport_status"] = (
            "direct_mc_candidate_for_ensemble_qualification"
            if status_ready
            else "direct_mc_unqualified"
        )
        return estimate, diagnostics
