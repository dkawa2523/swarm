"""Internal solver configuration projected from product schema v2."""

from __future__ import annotations

from dataclasses import dataclass, field as dc_field
from typing import Literal, Protocol


@dataclass(slots=True)
class EnergyGridRefinementConfig:
    enabled: bool = False
    threshold_padding_eV: float = 0.15
    points_per_threshold: int = 8
    max_extra_points: int = 120


@dataclass(slots=True)
class EnergyGridConfig:
    min_eV: float = 1.0e-4
    max_eV: float = 100.0
    n: int = 600
    spacing: Literal["linear", "quadratic", "log"] = "quadratic"
    refine: EnergyGridRefinementConfig = dc_field(
        default_factory=EnergyGridRefinementConfig
    )


@dataclass(slots=True)
class MultiTermEnergyGridConfig:
    min_eV: float = 1.0e-3
    max_eV: float = 200.0
    n: int = 220
    spacing: Literal["linear", "log", "log_linear"] = "log_linear"
    linear_until_eV: float = 2.0
    refine: EnergyGridRefinementConfig = dc_field(
        default_factory=EnergyGridRefinementConfig
    )


@dataclass(slots=True)
class AdaptiveGridConfig:
    enabled: bool = True
    # High-field EEDFs can require several 25% grid extensions after the
    # initial mean-energy jump. Four cycles allowed cases to stop far above
    # the requested tail-probability target without reaching max_max_eV.
    max_cycles: int = 16
    mean_energy_multiplier: float = 15.0
    min_max_eV: float = 20.0
    max_max_eV: float = 2000.0
    tail_probability: float = 1.0e-8
    tail_cells_fraction: float = 0.05
    edge_to_peak: float = 1.0e-10


@dataclass(slots=True)
class ConvergenceConfig:
    # Low-field argon cases can require more than 120 fixed-point updates to
    # meet the unchanged shape, eigenvalue, and residual tolerances.  Retain
    # the strict tolerances and permit convergence instead of accepting an
    # iteration-limit result as a formal table point.
    max_iterations: int = 600
    tolerance: float = 1.0e-8
    eigenvalue_tolerance: float = 1.0e-8
    # The residual is normalized with a unit floor, so low-field,
    # near-conservative argon cases report an absolute operator residual.
    # 1e-6 remains strict while avoiding false nonconvergence at the
    # sparse-solve/discretization floor (typically a few 1e-7).
    residual_tolerance: float = 1.0e-6
    # Conservative/weakly growing gases can form a fixed-point two-cycle
    # when the growth eigenvalue is updated too aggressively.
    relaxation: float = 0.1
    clip_negative: bool = True


class MomentumCollisionConfig(Protocol):
    """Solver-neutral momentum-collision regularization settings."""

    min_momentum_cross_section_m2: float


class InelasticCollisionConfig(Protocol):
    """Solver-neutral inelastic and ionization collision settings."""

    nonconservative_model: Literal["growth", "ignore"]
    ionization_energy_sharing: Literal[
        "equal", "primary_secondary", "loss_only"
    ]
    secondary_electron_energy_eV: float


@dataclass(slots=True)
class TwoTermInternalConfig:
    enabled: bool = True
    energy_grid: EnergyGridConfig = dc_field(default_factory=EnergyGridConfig)
    adaptive_grid: AdaptiveGridConfig = dc_field(default_factory=AdaptiveGridConfig)
    convergence: ConvergenceConfig = dc_field(default_factory=ConvergenceConfig)
    initial_electron_temperature_eV: float = 2.0
    nonconservative_model: Literal["growth", "ignore"] = "growth"
    ionization_energy_sharing: Literal["equal", "primary_secondary", "loss_only"] = (
        "equal"
    )
    secondary_electron_energy_eV: float = 0.0
    min_momentum_cross_section_m2: float = 1.0e-24


@dataclass(slots=True)
class MultiTermInternalConfig:
    enabled: bool = True
    lmax: int = 4
    method: Literal["pn_closure_direct", "pn_dcs"] = "pn_closure_direct"
    energy_grid: MultiTermEnergyGridConfig = dc_field(
        default_factory=MultiTermEnergyGridConfig
    )
    convergence: ConvergenceConfig = dc_field(default_factory=ConvergenceConfig)
    nonconservative_model: Literal["growth", "ignore"] = "growth"
    ionization_energy_sharing: Literal["equal", "primary_secondary", "loss_only"] = (
        "equal"
    )
    secondary_electron_energy_eV: float = 0.0
    min_momentum_cross_section_m2: float = 1.0e-24


@dataclass(slots=True)
class MonteCarloAdapterConfig:
    population_model: Literal[
        "fixed_particle_single_daughter", "weighted_branching"
    ] = "fixed_particle_single_daughter"
    seed: int | None = None
    particles: int | None = None
    warmup_collisions: int | None = None
    max_collisions: int | None = None
    tail_max_collisions: int | None = None
    tail_rate_rse_trigger: float = 0.25
    transport_correlation_lag_barriers: int = 64
    transport_estimator: Literal["single_field", "paired_field_parity"] = "single_field"
    numeric_kernel: Literal["auto", "python", "numba"] = "auto"
    collect_audit: bool = False


@dataclass(slots=True)
class PropagatorInternalConfig:
    method: Literal["stationary_response"] = "stationary_response"
    energy_cells: int = 600
    polar_cells: int = 72
    max_iterations: int = 2000
    convergence_tolerance: float = 1.0e-8
    max_memory_mb: int = 1024
    base_max_eV: float = 100.0
    max_eV_limit: float = 2000.0
    adaptive_grid: bool = True
    threshold_refinement: bool = True
    tail_probability_target: float = 1.0e-8


@dataclass(slots=True)
class InternalSolverConfigs:
    two_term: TwoTermInternalConfig = dc_field(default_factory=TwoTermInternalConfig)
    multi_term: MultiTermInternalConfig = dc_field(
        default_factory=MultiTermInternalConfig
    )
    monte_carlo: MonteCarloAdapterConfig = dc_field(
        default_factory=MonteCarloAdapterConfig
    )
    propagator: PropagatorInternalConfig = dc_field(
        default_factory=PropagatorInternalConfig
    )


def build_internal_solver_configs(
    solvers: object, physics: object
) -> InternalSolverConfigs:
    """Project product solver schema into implementation-only configs."""

    def refinement_config() -> EnergyGridRefinementConfig:
        return EnergyGridRefinementConfig(
            enabled=physics.energy_grid_policy.threshold_refinement
        )

    adaptive = AdaptiveGridConfig(
        enabled=physics.energy_grid_policy.adaptive,
        max_max_eV=physics.energy_grid_policy.max_eV_limit,
        tail_probability=physics.energy_grid_policy.tail_probability_target,
    )
    return InternalSolverConfigs(
        two_term=TwoTermInternalConfig(
            energy_grid=EnergyGridConfig(refine=refinement_config()),
            adaptive_grid=adaptive,
            nonconservative_model=solvers.two_term.nonconservative_model,
            ionization_energy_sharing=physics.ionization.energy_sharing,
            secondary_electron_energy_eV=(
                physics.ionization.secondary_electron_energy_eV
            ),
            min_momentum_cross_section_m2=(
                solvers.two_term.min_momentum_cross_section_m2
            ),
        ),
        multi_term=MultiTermInternalConfig(
            lmax=solvers.multi_term.lmax,
            method=solvers.multi_term.method,
            energy_grid=MultiTermEnergyGridConfig(refine=refinement_config()),
            ionization_energy_sharing=physics.ionization.energy_sharing,
            secondary_electron_energy_eV=(
                physics.ionization.secondary_electron_energy_eV
            ),
        ),
        monte_carlo=MonteCarloAdapterConfig(
            population_model=solvers.monte_carlo.population_model,
            seed=solvers.monte_carlo.seed,
            particles=solvers.monte_carlo.particles,
            warmup_collisions=solvers.monte_carlo.warmup_collisions,
            max_collisions=solvers.monte_carlo.max_collisions,
            tail_max_collisions=solvers.monte_carlo.tail_max_collisions,
            tail_rate_rse_trigger=solvers.monte_carlo.tail_rate_rse_trigger,
            transport_correlation_lag_barriers=(
                solvers.monte_carlo.transport_correlation_lag_barriers
            ),
            transport_estimator=solvers.monte_carlo.transport_estimator,
            numeric_kernel=solvers.monte_carlo.numeric_kernel,
        ),
        propagator=PropagatorInternalConfig(
            method=solvers.propagator.method,
            energy_cells=solvers.propagator.energy_cells,
            polar_cells=solvers.propagator.polar_cells,
            max_iterations=solvers.propagator.max_iterations,
            convergence_tolerance=solvers.propagator.convergence_tolerance,
            max_memory_mb=solvers.propagator.max_memory_mb,
            max_eV_limit=physics.energy_grid_policy.max_eV_limit,
            adaptive_grid=physics.energy_grid_policy.adaptive,
            threshold_refinement=physics.energy_grid_policy.threshold_refinement,
            tail_probability_target=(
                physics.energy_grid_policy.tail_probability_target
            ),
        ),
    )
