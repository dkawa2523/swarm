"""Canonical solver-section parsers."""

from __future__ import annotations

from typing import Any, Literal, cast

from electron_swarm.core.config import (
    MonteCarloProductConfig,
    MultiTermProductConfig,
    PropagatorProductConfig,
    SolversConfig,
    TwoTermBackend,
    TwoTermProductConfig,
)
from electron_swarm.core.config_sections.common import (
    as_mapping_section,
    canonical_solver,
    integer_field,
    reject_unknown_fields,
    strict_float,
    string_value,
    validate_literal,
)
from electron_swarm.core.legacy_schema import OBSOLETE_PUBLIC_NAMES
from electron_swarm.core.solver_ids import CANONICAL_SOLVER_IDS

TWO_TERM_FIELDS = {
    "backend",
    "nonconservative_model",
    "min_momentum_cross_section_m2",
}
MULTI_TERM_FIELDS = {"method", "lmax"}
MONTE_CARLO_FIELDS = {
    "population_model",
    "seed",
    "particles",
    "warmup_collisions",
    "max_collisions",
    "tail_max_collisions",
    "tail_rate_rse_trigger",
    "transport_correlation_lag_barriers",
    "transport_estimator",
    "numeric_kernel",
}
PROPAGATOR_FIELDS = {
    "method",
    "energy_cells",
    "polar_cells",
    "max_iterations",
    "convergence_tolerance",
    "max_memory_mb",
}


def parse_solvers(raw: dict[str, Any]) -> SolversConfig:
    solvers_raw = as_mapping_section(raw, "solvers", "solvers")
    unknown = set(solvers_raw) - set(CANONICAL_SOLVER_IDS)
    if unknown:
        for key in unknown:
            if key in OBSOLETE_PUBLIC_NAMES:
                canonical_solver(key, f"solvers.{key}")
        raise ValueError(f"Unsupported solver config sections: {sorted(unknown)}")

    tt_raw = as_mapping_section(solvers_raw, "two_term", "solvers.two_term")
    mt_raw = as_mapping_section(solvers_raw, "multi_term", "solvers.multi_term")
    mc_raw = as_mapping_section(solvers_raw, "monte_carlo", "solvers.monte_carlo")
    propagator_raw = as_mapping_section(
        solvers_raw,
        "propagator",
        "solvers.propagator",
    )
    reject_unknown_fields(tt_raw, TWO_TERM_FIELDS, "solvers.two_term")
    reject_unknown_fields(mt_raw, MULTI_TERM_FIELDS, "solvers.multi_term")
    reject_unknown_fields(mc_raw, MONTE_CARLO_FIELDS, "solvers.monte_carlo")
    reject_unknown_fields(
        propagator_raw,
        PROPAGATOR_FIELDS,
        "solvers.propagator",
    )

    two_term = TwoTermProductConfig(
        backend=cast(
            TwoTermBackend,
            validate_literal(
                string_value(
                    tt_raw.get("backend", "native_sg"),
                    "solvers.two_term.backend",
                ),
                {"native_sg"},
                "solvers.two_term.backend",
            ),
        ),
        nonconservative_model=cast(
            Literal["growth", "ignore"],
            validate_literal(
                string_value(
                    tt_raw.get("nonconservative_model", "growth"),
                    "solvers.two_term.nonconservative_model",
                ),
                {"growth", "ignore"},
                "solvers.two_term.nonconservative_model",
            ),
        ),
        min_momentum_cross_section_m2=strict_float(
            tt_raw.get("min_momentum_cross_section_m2", 1.0e-24),
            "solvers.two_term.min_momentum_cross_section_m2",
        ),
    )

    multi_term = MultiTermProductConfig(
        method=cast(
            Literal["pn_closure_direct", "pn_dcs"],
            validate_literal(
                string_value(
                    mt_raw.get("method", "pn_closure_direct"),
                    "solvers.multi_term.method",
                ),
                {"pn_closure_direct", "pn_dcs"},
                "solvers.multi_term.method",
            ),
        ),
        lmax=integer_field(mt_raw, "lmax", 4, "solvers.multi_term.lmax"),
    )

    mc_particles = (
        integer_field(mc_raw, "particles", 0, "solvers.monte_carlo.particles")
        if mc_raw.get("particles") is not None
        else None
    )
    mc_max_collisions = (
        integer_field(
            mc_raw,
            "max_collisions",
            0,
            "solvers.monte_carlo.max_collisions",
        )
        if mc_raw.get("max_collisions") is not None
        else None
    )
    mc_warmup_collisions = (
        integer_field(
            mc_raw,
            "warmup_collisions",
            0,
            "solvers.monte_carlo.warmup_collisions",
        )
        if mc_raw.get("warmup_collisions") is not None
        else None
    )
    mc_tail_max_collisions = (
        integer_field(
            mc_raw,
            "tail_max_collisions",
            0,
            "solvers.monte_carlo.tail_max_collisions",
        )
        if mc_raw.get("tail_max_collisions") is not None
        else None
    )
    mc_tail_rate_rse_trigger = strict_float(
        mc_raw.get("tail_rate_rse_trigger", 0.25),
        "solvers.monte_carlo.tail_rate_rse_trigger",
    )
    mc_transport_lag = integer_field(
        mc_raw,
        "transport_correlation_lag_barriers",
        64,
        "solvers.monte_carlo.transport_correlation_lag_barriers",
    )

    monte_carlo = MonteCarloProductConfig(
        population_model=cast(
            Literal["fixed_particle_single_daughter", "weighted_branching"],
            validate_literal(
                string_value(
                    mc_raw.get(
                        "population_model", "fixed_particle_single_daughter"
                    ),
                    "solvers.monte_carlo.population_model",
                ),
                {"fixed_particle_single_daughter", "weighted_branching"},
                "solvers.monte_carlo.population_model",
            ),
        ),
        seed=(
            integer_field(mc_raw, "seed", 0, "solvers.monte_carlo.seed")
            if mc_raw.get("seed") is not None
            else None
        ),
        particles=mc_particles,
        warmup_collisions=mc_warmup_collisions,
        max_collisions=mc_max_collisions,
        tail_max_collisions=mc_tail_max_collisions,
        tail_rate_rse_trigger=mc_tail_rate_rse_trigger,
        transport_correlation_lag_barriers=mc_transport_lag,
        transport_estimator=cast(
            Literal["single_field", "paired_field_parity"],
            validate_literal(
                string_value(
                    mc_raw.get("transport_estimator", "single_field"),
                    "solvers.monte_carlo.transport_estimator",
                ),
                {"single_field", "paired_field_parity"},
                "solvers.monte_carlo.transport_estimator",
            ),
        ),
        numeric_kernel=cast(
            Literal["auto", "python", "numba"],
            validate_literal(
                string_value(
                    mc_raw.get("numeric_kernel", "auto"),
                    "solvers.monte_carlo.numeric_kernel",
                ),
                {"auto", "python", "numba"},
                "solvers.monte_carlo.numeric_kernel",
            ),
        ),
    )
    propagator = PropagatorProductConfig(
        method=cast(
            Literal["stationary_response"],
            validate_literal(
                string_value(
                    propagator_raw.get("method", "stationary_response"),
                    "solvers.propagator.method",
                ),
                {"stationary_response"},
                "solvers.propagator.method",
            ),
        ),
        energy_cells=integer_field(
            propagator_raw,
            "energy_cells",
            600,
            "solvers.propagator.energy_cells",
        ),
        polar_cells=integer_field(
            propagator_raw,
            "polar_cells",
            72,
            "solvers.propagator.polar_cells",
        ),
        max_iterations=integer_field(
            propagator_raw,
            "max_iterations",
            2000,
            "solvers.propagator.max_iterations",
        ),
        convergence_tolerance=strict_float(
            propagator_raw.get("convergence_tolerance", 1.0e-8),
            "solvers.propagator.convergence_tolerance",
        ),
        max_memory_mb=integer_field(
            propagator_raw,
            "max_memory_mb",
            1024,
            "solvers.propagator.max_memory_mb",
        ),
    )
    return SolversConfig(
        two_term=two_term,
        multi_term=multi_term,
        monte_carlo=monte_carlo,
        propagator=propagator,
    )
