"""Internal solver configuration projected from product schema v2."""

from __future__ import annotations

from dataclasses import dataclass, field as dc_field
from typing import Literal


InternalBoltzmannBackend = Literal["native_bolsig"]


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


@dataclass(slots=True)
class TwoTermInternalConfig:
    enabled: bool = True
    backend: InternalBoltzmannBackend = "native_bolsig"
    energy_grid: EnergyGridConfig = dc_field(default_factory=EnergyGridConfig)
    adaptive_grid: AdaptiveGridConfig = dc_field(default_factory=AdaptiveGridConfig)
    convergence: ConvergenceConfig = dc_field(default_factory=ConvergenceConfig)
    initial_electron_temperature_eV: float = 2.0
    nonconservative_model: Literal["growth", "ignore"] = "growth"
    ionization_energy_sharing: Literal[
        "equal", "primary_secondary", "loss_only"
    ] = "equal"
    secondary_electron_energy_eV: float = 0.0
    min_momentum_cross_section_m2: float = 1.0e-24


@dataclass(slots=True)
class MultiTermInternalConfig:
    enabled: bool = True
    lmax: int = 4
    method: Literal["moment_closure"] = "moment_closure"
    product_method: Literal["pn_closure_direct", "pn_dcs"] = "pn_closure_direct"
    energy_grid: MultiTermEnergyGridConfig = dc_field(
        default_factory=MultiTermEnergyGridConfig
    )


@dataclass(slots=True)
class MonteCarloAdapterConfig:
    population_model: Literal[
        "fixed_particle_single_daughter", "weighted_branching"
    ] = "fixed_particle_single_daughter"
    seed: int | None = None
    particles: int | None = None
    warmup_collisions: int | None = None
    max_collisions: int | None = None
    collect_audit: bool = False


@dataclass(slots=True)
class InternalSolverConfigs:
    two_term: TwoTermInternalConfig = dc_field(default_factory=TwoTermInternalConfig)
    multi_term: MultiTermInternalConfig = dc_field(
        default_factory=MultiTermInternalConfig
    )
    monte_carlo: MonteCarloAdapterConfig = dc_field(
        default_factory=MonteCarloAdapterConfig
    )


def build_internal_solver_configs(solvers: object, physics: object) -> InternalSolverConfigs:
    """Project product solver schema into implementation-only configs."""

    refine = EnergyGridRefinementConfig(
        enabled=physics.energy_grid_policy.threshold_refinement
    )
    adaptive = AdaptiveGridConfig(
        enabled=physics.energy_grid_policy.adaptive,
        max_max_eV=physics.energy_grid_policy.max_eV_limit,
        tail_probability=physics.energy_grid_policy.tail_probability_target,
    )
    return InternalSolverConfigs(
        two_term=TwoTermInternalConfig(
            backend="native_bolsig",
            energy_grid=EnergyGridConfig(refine=refine),
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
            method="moment_closure",
            product_method=solvers.multi_term.method,
            energy_grid=MultiTermEnergyGridConfig(refine=refine),
        ),
        monte_carlo=MonteCarloAdapterConfig(
            population_model=solvers.monte_carlo.population_model,
            seed=solvers.monte_carlo.seed,
            particles=solvers.monte_carlo.particles,
            warmup_collisions=solvers.monte_carlo.warmup_collisions,
            max_collisions=solvers.monte_carlo.max_collisions,
        ),
    )
