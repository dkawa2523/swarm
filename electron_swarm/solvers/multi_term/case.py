"""Construction of a solver-ready multi-term case."""

from __future__ import annotations

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.constants import TOWNSEND
from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.core.solver_configs import MultiTermInternalConfig
from electron_swarm.physics.kinetics import gas_number_density

from .grid import make_energy_grid
from .models import MultiTermCase


def build_multiterm_case(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    multi_term_config: MultiTermInternalConfig,
    e_over_n_Td: float,
) -> MultiTermCase:
    number_density = gas_number_density(config)
    return MultiTermCase(
        config=config,
        cross_sections=cross_sections,
        multi_term_config=multi_term_config,
        grid=make_energy_grid(multi_term_config, cross_sections),
        e_over_n_Td=float(e_over_n_Td),
        gas_number_density_m3=number_density,
        electric_field_V_m=float(e_over_n_Td) * TOWNSEND * number_density,
    )


__all__ = ["build_multiterm_case"]
