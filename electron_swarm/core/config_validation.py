"""Physical invariants for schema-v2 product configuration.

Parsing checks syntax.  This module checks the resulting dataclasses so the
Python API and YAML entry point enforce the same small set of invariants.
"""

from __future__ import annotations

from math import isclose, isfinite
from numbers import Integral, Real
from typing import Collection

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.solver_ids import CANONICAL_SOLVER_IDS
from electron_swarm.core.legacy_schema import MIGRATION_ERROR


def _finite(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite number")
    result = float(value)
    if not isfinite(result):
        raise ValueError(f"{name} must be a finite number")
    return result


def _integer(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError(f"{name} must be an integer")
    return int(value)


def _boolean(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{name} must be a boolean")
    return value


def _choice(value: object, name: str, allowed: Collection[str]) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ValueError(f"{name} is invalid")
    return value


def _positive(value: object, name: str, *, allow_zero: bool = False) -> float:
    result = _finite(value, name)
    invalid = result < 0.0 if allow_zero else result <= 0.0
    if invalid:
        qualifier = "nonnegative" if allow_zero else "positive"
        raise ValueError(f"{name} must be {qualifier}")
    return result


def _optional_positive(
    value: object | None,
    name: str,
    *,
    allow_zero: bool = False,
) -> None:
    if value is not None:
        _positive(value, name, allow_zero=allow_zero)


def validate_config(config: SwarmConfig) -> None:
    """Reject malformed or physically impossible product configuration."""

    if isinstance(config.schema_version, bool) or config.schema_version != 2:
        raise ValueError(MIGRATION_ERROR)

    if not config.run.solvers:
        raise ValueError("run.solvers must not be empty")
    solver_ids: list[str] = []
    for index, requested in enumerate(config.run.solvers):
        if requested.id not in CANONICAL_SOLVER_IDS:
            raise ValueError(f"run.solvers[{index}].id is not canonical")
        if not isinstance(requested.enabled, bool):
            raise ValueError(f"run.solvers[{index}].enabled must be a boolean")
        solver_ids.append(requested.id)
    if len(set(solver_ids)) != len(solver_ids):
        raise ValueError("run.solvers must not contain duplicate solver ids")
    if not config.run.e_over_n_Td:
        raise ValueError("run.e_over_n_Td must not be empty")
    for index, value in enumerate(config.run.e_over_n_Td):
        _positive(value, f"run.e_over_n_Td[{index}]")
    if not isinstance(config.run.case_prefix, str) or not config.run.case_prefix.strip():
        raise ValueError("run.case_prefix must be a non-empty string")

    conditions = config.conditions
    _positive(conditions.gas_temperature_K, "conditions.gas_temperature_K")
    _optional_positive(conditions.pressure_Pa, "conditions.pressure_Pa")
    _optional_positive(
        conditions.gas_number_density_m3,
        "conditions.gas_number_density_m3",
    )
    _optional_positive(conditions.length_scale_m, "conditions.length_scale_m")
    if conditions.pressure_Pa is None and conditions.gas_number_density_m3 is None:
        raise ValueError(
            "either conditions.pressure_Pa or conditions.gas_number_density_m3 "
            "must be set"
        )
    if conditions.pressure_Pa is not None and conditions.gas_number_density_m3 is not None:
        raise ValueError(
            "set only one of conditions.pressure_Pa and "
            "conditions.gas_number_density_m3"
        )
    if not conditions.gas_mixture:
        raise ValueError("conditions.gas_mixture must not be empty")
    species: list[str] = []
    fraction_sum = 0.0
    for index, component in enumerate(conditions.gas_mixture):
        name = component.species
        if not isinstance(name, str) or not name.strip():
            raise ValueError(
                f"conditions.gas_mixture[{index}].species must be non-empty"
            )
        species.append(name)
        fraction_sum += _positive(
            component.fraction,
            f"conditions.gas_mixture[{index}].fraction",
            allow_zero=True,
        )
        _positive(
            component.mass_amu,
            f"conditions.gas_mixture[{index}].mass_amu",
        )
    if len(set(species)) != len(species):
        raise ValueError("conditions.gas_mixture species must be unique")
    if not isclose(fraction_sum, 1.0, rel_tol=1.0e-12, abs_tol=1.0e-12):
        raise ValueError("conditions.gas_mixture fractions must sum to 1")

    if not config.cross_sections.files:
        raise ValueError("cross_sections.files must not be empty")
    if config.cross_sections.high_energy_extrapolation not in {"zero", "hold", "error"}:
        raise ValueError("cross_sections.high_energy_extrapolation is invalid")

    field = config.physics.field
    _choice(field.type, "physics.field.type", {"dc", "rf", "time_dependent"})
    magnetic = field.magnetic_field
    _boolean(magnetic.enabled, "physics.field.magnetic_field.enabled")
    _positive(
        magnetic.B_T,
        "physics.field.magnetic_field.B_T",
        allow_zero=True,
    )
    angle = _finite(magnetic.angle_EB_deg, "physics.field.magnetic_field.angle_EB_deg")
    if not 0.0 <= angle <= 180.0:
        raise ValueError("physics.field.magnetic_field.angle_EB_deg must be in [0, 180]")

    periodic = field.time_dependent
    _choice(
        periodic.waveform,
        "physics.field.time_dependent.waveform",
        {"sinusoidal"},
    )
    _choice(
        periodic.amplitude_definition,
        "physics.field.time_dependent.amplitude_definition",
        {"rms"},
    )
    _choice(
        periodic.momentum_response,
        "physics.field.time_dependent.momentum_response",
        {"instantaneous"},
    )
    if field.type == "time_dependent":
        _positive(periodic.frequency_Hz, "physics.field.time_dependent.frequency_Hz")
    elif periodic.frequency_Hz is not None:
        _positive(periodic.frequency_Hz, "physics.field.time_dependent.frequency_Hz")
    phase_steps = _integer(periodic.phase_steps, "physics.field.time_dependent.phase_steps")
    if phase_steps < 8 or phase_steps % 2:
        raise ValueError("physics.field.time_dependent.phase_steps must be even and >= 8")
    if _integer(periodic.max_periods, "physics.field.time_dependent.max_periods") <= 0:
        raise ValueError("physics.field.time_dependent.max_periods must be positive")
    tolerance = _positive(
        periodic.periodic_tolerance,
        "physics.field.time_dependent.periodic_tolerance",
    )
    if tolerance >= 1.0:
        raise ValueError("physics.field.time_dependent.periodic_tolerance must be < 1")

    angular = config.physics.angular_scattering
    expected_closure = {
        "isotropic": "zero",
        "momentum_power": "power",
        "maxent_p1": "maxent",
        "moment_table": "table",
    }
    model = _choice(
        angular.model,
        "physics.angular_scattering.model",
        set(expected_closure),
    )
    _choice(
        angular.higher_moment_closure,
        "physics.angular_scattering.higher_moment_closure",
        {"zero", "power", "maxent", "table"},
    )
    if angular.higher_moment_closure != expected_closure[model]:
        raise ValueError(
            "physics.angular_scattering model and higher_moment_closure mismatch"
        )
    if model == "moment_table" and angular.moment_table is None:
        raise ValueError(
            "physics.angular_scattering.moment_table is required for moment_table"
        )
    if model != "moment_table" and angular.moment_table is not None:
        raise ValueError(
            "physics.angular_scattering.moment_table requires model=moment_table"
        )
    if angular.moment_table is not None:
        table = angular.moment_table
        _choice(
            table.format,
            "physics.angular_scattering.moment_table.format",
            {"normalized_legendre_moments"},
        )
        _choice(
            table.provenance,
            "physics.angular_scattering.moment_table.provenance",
            {"dcs_derived", "model_derived", "unknown"},
        )
        _choice(
            table.extrapolation,
            "physics.angular_scattering.moment_table.extrapolation",
            {"error"},
        )

    electron_electron = config.physics.electron_electron
    _boolean(electron_electron.enabled, "physics.electron_electron.enabled")
    _boolean(
        electron_electron.conserve_mean_energy,
        "physics.electron_electron.conserve_mean_energy",
    )
    electron_electron_model = _choice(
        electron_electron.model,
        "physics.electron_electron.model",
        {"none", "relaxation_postprocess", "fp_energy"},
    )
    _choice(
        electron_electron.strength_model,
        "physics.electron_electron.strength_model",
        {"simple_relaxation", "density_based"},
    )
    if electron_electron.enabled and electron_electron_model == "none":
        raise ValueError(
            "physics.electron_electron.model is required when enabled"
        )
    if not electron_electron.enabled and electron_electron_model != "none":
        raise ValueError(
            "physics.electron_electron.model must be none when disabled"
        )
    relaxation = _finite(
        electron_electron.relaxation_fraction,
        "physics.electron_electron.relaxation_fraction",
    )
    if not 0.0 <= relaxation <= 1.0:
        raise ValueError("physics.electron_electron.relaxation_fraction must be in [0, 1]")
    _positive(
        electron_electron.fallback_temperature_eV,
        "physics.electron_electron.fallback_temperature_eV",
    )
    ionization = config.physics.ionization
    energy_sharing = _choice(
        ionization.energy_sharing,
        "physics.ionization.energy_sharing",
        {"equal", "primary_secondary", "loss_only"},
    )
    secondary_energy = _positive(
        ionization.secondary_electron_energy_eV,
        "physics.ionization.secondary_electron_energy_eV",
        allow_zero=True,
    )
    if energy_sharing != "primary_secondary" and secondary_energy != 0.0:
        raise ValueError(
            "physics.ionization.secondary_electron_energy_eV requires "
            "energy_sharing=primary_secondary"
        )

    grid = config.physics.energy_grid_policy
    _boolean(grid.adaptive, "physics.energy_grid_policy.adaptive")
    _boolean(
        grid.threshold_refinement,
        "physics.energy_grid_policy.threshold_refinement",
    )
    _boolean(grid.tail_metrics, "physics.energy_grid_policy.tail_metrics")
    tail_target = _positive(
        grid.tail_probability_target,
        "physics.energy_grid_policy.tail_probability_target",
    )
    if tail_target >= 1.0:
        raise ValueError("physics.energy_grid_policy.tail_probability_target must be < 1")
    _optional_positive(
        grid.tail_threshold_eV,
        "physics.energy_grid_policy.tail_threshold_eV",
        allow_zero=True,
    )
    warning_fraction = _finite(
        grid.tail_rate_warning_fraction,
        "physics.energy_grid_policy.tail_rate_warning_fraction",
    )
    if not 0.0 <= warning_fraction <= 1.0:
        raise ValueError("physics.energy_grid_policy.tail_rate_warning_fraction must be in [0, 1]")
    _positive(grid.max_eV_limit, "physics.energy_grid_policy.max_eV_limit")

    finite_k = config.physics.finite_k
    _boolean(finite_k.enabled, "physics.finite_k.enabled")
    if finite_k.enabled:
        _positive(finite_k.k_m_inv, "physics.finite_k.k_m_inv")
    elif finite_k.k_m_inv is not None:
        _positive(finite_k.k_m_inv, "physics.finite_k.k_m_inv")

    two_term = config.solvers.two_term
    _choice(two_term.backend, "solvers.two_term.backend", {"native_sg"})
    _choice(
        two_term.nonconservative_model,
        "solvers.two_term.nonconservative_model",
        {"growth", "ignore"},
    )
    _positive(
        two_term.min_momentum_cross_section_m2,
        "solvers.two_term.min_momentum_cross_section_m2",
    )

    multi_term = config.solvers.multi_term
    _choice(
        multi_term.method,
        "solvers.multi_term.method",
        {"pn_closure_direct", "pn_dcs"},
    )
    if _integer(multi_term.lmax, "solvers.multi_term.lmax") < 1:
        raise ValueError("solvers.multi_term.lmax must be >= 1")

    monte_carlo = config.solvers.monte_carlo
    population_model = _choice(
        monte_carlo.population_model,
        "solvers.monte_carlo.population_model",
        {"fixed_particle_single_daughter", "weighted_branching"},
    )
    _choice(
        monte_carlo.transport_estimator,
        "solvers.monte_carlo.transport_estimator",
        {"single_field", "paired_field_parity"},
    )
    _choice(
        monte_carlo.numeric_kernel,
        "solvers.monte_carlo.numeric_kernel",
        {"auto", "python", "numba"},
    )
    if monte_carlo.seed is not None:
        seed = _integer(monte_carlo.seed, "solvers.monte_carlo.seed")
        if seed < 0:
            raise ValueError("solvers.monte_carlo.seed must be nonnegative")
    for name in ("particles", "max_collisions", "tail_max_collisions"):
        value = getattr(monte_carlo, name)
        if value is not None and _integer(value, f"solvers.monte_carlo.{name}") <= 0:
            raise ValueError(f"solvers.monte_carlo.{name} must be positive")
    if monte_carlo.warmup_collisions is not None and _integer(
        monte_carlo.warmup_collisions,
        "solvers.monte_carlo.warmup_collisions",
    ) < 0:
        raise ValueError("solvers.monte_carlo.warmup_collisions must be nonnegative")
    tail_rse = _positive(
        monte_carlo.tail_rate_rse_trigger,
        "solvers.monte_carlo.tail_rate_rse_trigger",
    )
    if tail_rse > 1.0:
        raise ValueError("solvers.monte_carlo.tail_rate_rse_trigger must be <= 1")
    lag = _integer(
        monte_carlo.transport_correlation_lag_barriers,
        "solvers.monte_carlo.transport_correlation_lag_barriers",
    )
    if lag < 4 or lag & (lag - 1):
        raise ValueError(
            "solvers.monte_carlo.transport_correlation_lag_barriers must be a "
            "power of two >= 4"
        )
    if population_model == "fixed_particle_single_daughter":
        if monte_carlo.tail_max_collisions is not None:
            raise ValueError(
                "solvers.monte_carlo.tail_max_collisions requires "
                "population_model=weighted_branching"
            )
        if lag != 64:
            raise ValueError(
                "solvers.monte_carlo.transport_correlation_lag_barriers requires "
                "population_model=weighted_branching"
            )

    propagator = config.solvers.propagator
    _choice(
        propagator.method,
        "solvers.propagator.method",
        {"stationary_response"},
    )
    energy_cells = _integer(propagator.energy_cells, "solvers.propagator.energy_cells")
    polar_cells = _integer(propagator.polar_cells, "solvers.propagator.polar_cells")
    iterations = _integer(propagator.max_iterations, "solvers.propagator.max_iterations")
    memory = _integer(propagator.max_memory_mb, "solvers.propagator.max_memory_mb")
    convergence = _positive(
        propagator.convergence_tolerance,
        "solvers.propagator.convergence_tolerance",
    )
    if not 64 <= energy_cells <= 4096:
        raise ValueError("solvers.propagator.energy_cells must be in [64, 4096]")
    if not 8 <= polar_cells <= 360 or polar_cells % 2:
        raise ValueError("solvers.propagator.polar_cells must be even and in [8, 360]")
    if not 10 <= iterations <= 20000:
        raise ValueError("solvers.propagator.max_iterations must be in [10, 20000]")
    if not 1.0e-12 <= convergence <= 1.0e-4:
        raise ValueError("solvers.propagator.convergence_tolerance must be in [1e-12, 1e-4]")
    if memory < 128:
        raise ValueError("solvers.propagator.max_memory_mb must be >= 128")

    comparison = config.comparison
    _boolean(comparison.enabled, "comparison.enabled")
    _boolean(comparison.compare_eedf, "comparison.compare_eedf")
    _boolean(comparison.required, "comparison.required")
    if (
        comparison.reference_solver is not None
        and comparison.reference_solver not in CANONICAL_SOLVER_IDS
    ):
        raise ValueError("comparison.reference_solver is not canonical")
    for index, solver in enumerate(comparison.candidate_solvers):
        if solver not in CANONICAL_SOLVER_IDS:
            raise ValueError(f"comparison.candidate_solvers[{index}] is not canonical")

    _choice(
        config.feature_policy.unsupported,
        "feature_policy.unsupported",
        {"fail", "skip_solver"},
    )
    _choice(
        config.feature_policy.degraded,
        "feature_policy.degraded",
        {"fail", "record"},
    )
