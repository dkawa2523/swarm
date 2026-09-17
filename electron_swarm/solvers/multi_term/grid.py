"""Multi-term projection onto the shared energy-grid builder."""

from __future__ import annotations

from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.core.solver_configs import MultiTermInternalConfig
from electron_swarm.grids.energy import EnergyGrid, build_energy_grid


def make_energy_grid(
    solver_config: MultiTermInternalConfig,
    cross_sections: CrossSectionSet | None = None,
) -> EnergyGrid:
    grid = solver_config.energy_grid
    refinement = grid.refine
    return build_energy_grid(
        min_eV=grid.min_eV,
        max_eV=grid.max_eV,
        n=grid.n,
        spacing=grid.spacing,
        linear_until_eV=grid.linear_until_eV,
        cross_sections=cross_sections,
        refine=refinement.enabled,
        threshold_padding_eV=refinement.threshold_padding_eV,
        points_per_threshold=refinement.points_per_threshold,
        max_extra_points=refinement.max_extra_points,
    )


__all__ = ["make_energy_grid"]
