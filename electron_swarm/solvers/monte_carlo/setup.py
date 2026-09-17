"""Validated static preparation for internal Monte Carlo runs."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from electron_swarm.core.config import MagneticFieldConfig, SwarmConfig
from electron_swarm.core.cross_sections import ActiveMixtureInputs, ProcessType
from electron_swarm.core.solver_configs import MonteCarloAdapterConfig
from electron_swarm.physics.angular_scattering import (
    AngularMomentProvider,
    build_angular_model,
)
from electron_swarm.physics.kinetics import gas_number_density
import electron_swarm.solvers.monte_carlo.result_evidence as _result_evidence
import electron_swarm.solvers.monte_carlo.weighted_runtime as _weighted_runtime
from electron_swarm.solvers.monte_carlo.collisions import (
    TRAJECTORY_REACTION_TYPES,
    PreparedCollisionSampler,
    _ProjectedProcess,
    _is_scattering_event_process,
    _max_cross_section_energy,
    _project_processes,
)
from electron_swarm.solvers.monte_carlo.histogram import _tail_threshold_eV
from electron_swarm.solvers.monte_carlo.kinematics import magnetic_field_vector
from electron_swarm.solvers.monte_carlo.source_identity import (
    monte_carlo_source_sha256,
)
from electron_swarm.solvers.monte_carlo.weighted_transport import (
    WeightedGrowthTransportLagPlan,
    weighted_growth_transport_lag_plan,
)
from electron_swarm.solvers.monte_carlo.weighted_ensemble import (
    REACTION_WE_RESAMPLE_TRIAL_PERIODS,
    reaction_kernel_strata_importance,
    threshold_strata_edges,
)


@dataclass(frozen=True, slots=True)
class _MonteCarloRunSetup:
    config: SwarmConfig
    cfg: MonteCarloAdapterConfig
    collect_audit: bool
    particles: int
    population_model: Literal["fixed_particle_single_daughter", "weighted_branching"]
    branching_active: bool
    transport_lag_plan: WeightedGrowthTransportLagPlan
    collisions: int
    tail_collisions: int
    warmup_collisions: int
    density: float
    projected: list[_ProjectedProcess]
    angular_model: AngularMomentProvider
    angular_name: str
    event_sampled: tuple[bool, ...]
    max_energy_limit: float
    collision_sampler: PreparedCollisionSampler
    trial_frequency: float
    barrier_frequency: float
    barrier_dt: float
    max_cross_section_energy: float
    tail_threshold: float
    reaction_thresholds: tuple[float, ...]
    tail_strata_edges: np.ndarray | None
    tail_reaction_importance: np.ndarray | None
    tail_resample_interval: int | None
    zero_high_energy_policy: bool
    magnetic: MagneticFieldConfig
    magnetic_enabled: bool
    effective_B_T: float
    B: np.ndarray
    solver_source_sha256: str
    compiled_plan: _weighted_runtime._CompiledKernelPlan | None
    compiled_fallback_reason: str | None


def _validate_trajectory_processes(projected: list[_ProjectedProcess]) -> None:
    for item in projected:
        process = item.process
        if process.process_type == ProcessType.ATTACHMENT:
            raise NotImplementedError(
                "internal monte_carlo does not yet implement attachment particle "
                "loss in either population model; refusing to report an EEDF or "
                "transport from trajectories that omit the sink"
            )


def _prepare_monte_carlo_run(
    config: SwarmConfig,
    cross_sections: ActiveMixtureInputs,
    solver_config: MonteCarloAdapterConfig,
    *,
    collect_audit: bool,
) -> _MonteCarloRunSetup:
    cfg = solver_config
    particles = int(cfg.particles or 256)
    population_model = cfg.population_model
    branching_active = _result_evidence._weighted_branching_active(
        config,
        population_model,
    )
    transport_lag_plan = weighted_growth_transport_lag_plan(
        cfg.transport_correlation_lag_barriers
    )
    if not branching_active and transport_lag_plan.correlation_lag_barriers != 64:
        raise ValueError(
            "transport_correlation_lag_barriers is only used by an active "
            "weighted-branching Monte Carlo estimator"
        )
    collisions = int(cfg.max_collisions or 100)
    # The rare-tail phase is an explicit MC opt-in.  Keeping the disabled
    # state as a zero executable budget prevents the production length from
    # being misreported (and potentially spent) as tail work.
    tail_collisions = (
        int(cfg.tail_max_collisions)
        if cfg.tail_max_collisions is not None
        else 0
    )
    warmup_collisions = int(cfg.warmup_collisions or 0)
    density = gas_number_density(config)
    projected = _project_processes(config, cross_sections)
    _validate_trajectory_processes(projected)
    if population_model == "weighted_branching":
        if config.physics.ionization.energy_sharing == "loss_only":
            raise NotImplementedError(
                "monte_carlo population_model=weighted_branching requires an "
                "explicit two-daughter ionization model; energy_sharing=loss_only "
                "would silently select fixed-population transport"
            )
    angular_model = build_angular_model(config)
    angular_name = config.physics.angular_scattering.model
    if angular_name not in {"isotropic", "maxent_p1"}:
        raise NotImplementedError(
            f"internal monte_carlo has no product sampler for angular model {angular_name!r}"
        )
    event_sampled = tuple(
        _is_scattering_event_process(item, projected, angular_name)
        or item.process.process_type in TRAJECTORY_REACTION_TYPES
        for item in projected
    )
    max_energy_limit = float(config.physics.energy_grid_policy.max_eV_limit)
    collision_sampler = PreparedCollisionSampler.build(
        config,
        projected,
        angular_model_name=angular_name,
        max_energy_eV=max_energy_limit,
    )
    trial_frequency = collision_sampler.trial_collision_frequency_s_inv(density)
    barrier_frequency = trial_frequency
    barrier_dt = 1.0 / float(barrier_frequency)
    max_cross_section_energy = _max_cross_section_energy(projected)
    tail_threshold = _tail_threshold_eV(item.process for item in projected)
    reaction_thresholds = tuple(
        float(item.process.threshold_eV)
        for item in projected
        if item.process.threshold_eV is not None
        and np.isfinite(item.process.threshold_eV)
        and item.process.threshold_eV >= 0.0
    )
    # Tail sampling is an MC estimator request.  It must not be disabled by
    # ``energy_grid_policy.threshold_refinement``, which controls deterministic
    # solver grids.  An explicit MC tail budget is the sole public opt-in.
    tail_sampling_requested = tail_collisions > 0
    tail_strata_edges = (
        threshold_strata_edges(
            reaction_thresholds,
            max_energy_eV=max_energy_limit,
        )
        if branching_active and tail_sampling_requested
        else None
    )
    tail_resample_interval = (
        REACTION_WE_RESAMPLE_TRIAL_PERIODS if tail_strata_edges is not None else None
    )
    tail_reaction_importance = (
        reaction_kernel_strata_importance(
            tuple(item.process for item in projected),
            tail_strata_edges,
        )
        if tail_strata_edges is not None
        else None
    )
    zero_high_energy_policy = config.cross_sections.high_energy_extrapolation == "zero"

    magnetic = config.physics.field.magnetic_field
    magnetic_enabled = bool(magnetic.enabled and magnetic.B_T > 0.0)
    effective_B_T = float(magnetic.B_T) if magnetic_enabled else 0.0
    B = magnetic_field_vector(effective_B_T, magnetic.angle_EB_deg)
    if cfg.transport_estimator not in {"single_field", "paired_field_parity"}:
        raise ValueError("unsupported Monte Carlo transport estimator")
    if cfg.transport_estimator == "paired_field_parity" and (
        config.physics.field.type != "dc" or magnetic_enabled
    ):
        raise NotImplementedError(
            "monte_carlo transport_estimator=paired_field_parity requires "
            "a DC electric field with magnetic field disabled"
        )
    compiled_plan, compiled_fallback_reason = _weighted_runtime._compiled_kernel_plan(
        config=config,
        solver_config=cfg,
        projected=projected,
        angular_name=angular_name,
        branching_active=branching_active,
        magnetic_enabled=magnetic_enabled,
        collect_audit=collect_audit,
    )
    return _MonteCarloRunSetup(
        config=config,
        cfg=cfg,
        collect_audit=collect_audit,
        particles=particles,
        population_model=population_model,
        branching_active=branching_active,
        transport_lag_plan=transport_lag_plan,
        collisions=collisions,
        tail_collisions=tail_collisions,
        warmup_collisions=warmup_collisions,
        density=density,
        projected=projected,
        angular_model=angular_model,
        angular_name=angular_name,
        event_sampled=event_sampled,
        max_energy_limit=max_energy_limit,
        collision_sampler=collision_sampler,
        trial_frequency=trial_frequency,
        barrier_frequency=barrier_frequency,
        barrier_dt=barrier_dt,
        max_cross_section_energy=max_cross_section_energy,
        tail_threshold=tail_threshold,
        reaction_thresholds=reaction_thresholds,
        tail_strata_edges=tail_strata_edges,
        tail_reaction_importance=tail_reaction_importance,
        tail_resample_interval=tail_resample_interval,
        zero_high_energy_policy=zero_high_energy_policy,
        magnetic=magnetic,
        magnetic_enabled=magnetic_enabled,
        effective_B_T=effective_B_T,
        B=B,
        solver_source_sha256=monte_carlo_source_sha256(),
        compiled_plan=compiled_plan,
        compiled_fallback_reason=compiled_fallback_reason,
    )
