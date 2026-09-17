"""Energy-angle finite-volume grid construction."""

from __future__ import annotations

import numpy as np

from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.core.solver_configs import PropagatorInternalConfig
from electron_swarm.physics.electron_neutral import speed_from_energy_m_s
from electron_swarm.solvers.propagator.models import PropagatorGrid


_SPEED_CLUSTERING_STRENGTH = 2.5


def _threshold_edges(
    cross_sections: CrossSectionSet,
    upper_eV: float,
) -> np.ndarray:
    return np.asarray(
        sorted(
            {
                float(process.threshold_eV)
                for process in cross_sections.processes
                if process.threshold_eV is not None
                and 0.0 < float(process.threshold_eV) < upper_eV
            }
        ),
        dtype=float,
    )


def _stretched_tail_edges(
    base_edges: np.ndarray,
    upper_eV: float,
) -> np.ndarray:
    if upper_eV <= base_edges[-1] * (1.0 + 1.0e-14):
        return base_edges
    edges = list(map(float, base_edges))
    width = max(edges[-1] - edges[-2], 1.0e-6)
    while edges[-1] < upper_eV:
        width *= 1.08
        edges.append(min(upper_eV, edges[-1] + width))
    return np.asarray(edges, dtype=float)


def build_propagator_grid(
    config: PropagatorInternalConfig,
    cross_sections: CrossSectionSet,
    *,
    max_eV: float | None = None,
) -> PropagatorGrid:
    """Build one core-plus-tail grid for the complete stationary solve."""

    core_upper = min(float(config.base_max_eV), float(config.max_eV_limit))
    requested_upper = (
        core_upper
        if max_eV is None
        else min(float(max_eV), float(config.max_eV_limit))
    )
    if requested_upper <= 0.0:
        raise ValueError("propagator energy upper bound must be positive")
    coordinate = np.linspace(0.0, 1.0, int(config.energy_cells) + 1)
    # A uniform-speed grid spends too many cells near the core ceiling while
    # low-field transport is controlled by the thermal/sub-eV region.  A sinh
    # map is linear at v=0 (so the characteristic origin retains its regular
    # geometry), nested under factor-two refinement, and smoothly stretches
    # toward the high-energy end without an adaptive solve/remap cycle.
    strength = _SPEED_CLUSTERING_STRENGTH
    speed_coordinate = np.sinh(strength * coordinate) / np.sinh(strength)
    edges = core_upper * speed_coordinate * speed_coordinate
    edges = _stretched_tail_edges(edges, requested_upper)
    if config.threshold_refinement:
        edges = np.unique(
            np.concatenate([edges, _threshold_edges(cross_sections, edges[-1])])
        )
    if np.any(~np.isfinite(edges)) or np.any(np.diff(edges) <= 0.0):
        raise ValueError("propagator energy edges must be finite and increasing")
    if len(edges) - 1 > 4 * int(config.energy_cells):
        raise MemoryError(
            "propagator adaptive grid exceeded four times the configured core cells"
        )

    energy_widths = np.diff(edges)
    energy_centers = 0.5 * (edges[:-1] + edges[1:])
    speed_edges = np.asarray(speed_from_energy_m_s(edges), dtype=float)
    speed_centers = np.asarray(
        speed_from_energy_m_s(energy_centers),
        dtype=float,
    )

    theta_edges = np.linspace(0.0, np.pi, int(config.polar_cells) + 1)
    theta_centers = 0.5 * (theta_edges[:-1] + theta_edges[1:])
    cos_edges = np.cos(theta_edges)
    mu_widths = cos_edges[:-1] - cos_edges[1:]
    mu_centers = 0.5 * (cos_edges[:-1] + cos_edges[1:])
    solid_angles = 2.0 * np.pi * mu_widths
    shell_volume = (
        2.0
        * np.pi
        / 3.0
        * (speed_edges[1:] ** 3 - speed_edges[:-1] ** 3)
    )
    cell_volumes = shell_volume[:, None] * mu_widths[None, :]
    if (
        np.any(mu_widths <= 0.0)
        or np.any(solid_angles <= 0.0)
        or np.any(cell_volumes <= 0.0)
    ):
        raise FloatingPointError("propagator grid has a nonpositive cell measure")
    return PropagatorGrid(
        energy_edges_eV=edges,
        energy_centers_eV=energy_centers,
        energy_widths_eV=energy_widths,
        speed_edges_m_s=speed_edges,
        speed_centers_m_s=speed_centers,
        theta_edges_rad=theta_edges,
        theta_centers_rad=theta_centers,
        mu_centers=mu_centers,
        mu_widths=mu_widths,
        solid_angles_sr=solid_angles,
        cell_volumes_v3=cell_volumes,
    )
