"""Memory estimates shared by Propagator setup and result evidence."""

from __future__ import annotations

import numpy as np

from electron_swarm.solvers.propagator.steady import krylov_subspace_dimension


def preflight_memory_bytes(
    *,
    energy_cells: int,
    polar_cells: int,
    angular_model: str,
) -> int:
    cells = int(energy_cells) * int(polar_cells)
    arrays = 16 * cells * np.dtype(float).itemsize
    angular = (
        int(energy_cells)
        * int(polar_cells)
        * int(polar_cells)
        * np.dtype(float).itemsize
        if angular_model == "maxent_p1"
        else 0
    )
    krylov_vectors = (
        krylov_subspace_dimension(cells) + 8
    ) * cells * np.dtype(float).itemsize
    response_blocks = (
        4
        * int(energy_cells)
        * int(polar_cells)
        * int(polar_cells)
        * np.dtype(float).itemsize
    )
    boundary_factor = (
        4
        * int(energy_cells)
        * int(polar_cells)
        * int(polar_cells)
        * (np.dtype(float).itemsize + np.dtype(np.int32).itemsize)
    )
    origin_bound = 12 * cells * np.dtype(float).itemsize
    return int(
        arrays
        + angular
        + krylov_vectors
        + response_blocks
        + boundary_factor
        + origin_bound
    )
