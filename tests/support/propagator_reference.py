"""Independent deterministic references for propagator verification.

These routines are deliberately excluded from the production solve.  They
integrate the continuous discrete-ordinate shell boundary-value problem with
an adaptive high-order ODE method, then eliminate the opposing half-range
boundary data algebraically.  This supplies an independent reference for the
production Cayley/Redheffer response and its embedded error estimates.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from electron_swarm.solvers.propagator.collisions import (
    angular_collision_loss_matrix,
)
from electron_swarm.solvers.propagator.models import (
    CollisionOperatorData,
    PropagatorGrid,
)
from electron_swarm.solvers.propagator.origin_response import (
    AffineOriginResponse,
)
from electron_swarm.solvers.propagator.shell_response import (
    AffineShellResponse,
    half_range_boundary_order,
    spherical_angular_streaming_matrix,
)


@dataclass(frozen=True, slots=True)
class ContinuousShellReference:
    response: AffineShellResponse
    function_evaluations: int
    boundary_block_condition_number: float


def continuous_shell_reference(
    grid: PropagatorGrid,
    collisions: CollisionOperatorData,
    energy_index: int,
    acceleration_m_s2: float,
    local_scalar_loss_s_inv: float,
    *,
    relative_tolerance: float = 2.0e-11,
    absolute_tolerance: float = 1.0e-13,
    maximum_polar_cells: int = 24,
) -> ContinuousShellReference:
    """Solve one regular shell by continuous high-order ODE integration.

    The reference uses the same finite-volume angular and collision
    coefficients as production, but it does not use thin Cayley layers,
    doubling, Redheffer composition, or Richardson extrapolation.  It is
    intentionally bounded to small angular grids because the augmented
    fundamental matrix has ``O(N_mu**2)`` states.
    """

    if energy_index <= 0 or energy_index >= grid.energy_cells:
        raise ValueError("continuous shell reference requires a regular shell")
    acceleration = float(acceleration_m_s2)
    loss = float(local_scalar_loss_s_inv)
    if acceleration <= 0.0 or not np.isfinite(acceleration):
        raise ValueError("continuous shell reference requires positive acceleration")
    if loss < 0.0 or not np.isfinite(loss):
        raise ValueError("continuous shell reference requires nonnegative loss")
    polar_cells = int(grid.polar_cells)
    if polar_cells > int(maximum_polar_cells):
        raise MemoryError(
            "continuous shell reference is limited to "
            f"{maximum_polar_cells} polar cells"
        )

    mu = np.asarray(grid.mu_centers, dtype=float)
    mu_widths = np.asarray(grid.mu_widths, dtype=float)
    order = half_range_boundary_order(mu)
    half = polar_cells // 2
    positive = order[:half]
    negative = order[half:]
    left_speed = float(grid.speed_edges_m_s[energy_index])
    right_speed = float(grid.speed_edges_m_s[energy_index + 1])
    shell_width = right_speed - left_speed
    if left_speed <= 0.0 or shell_width <= 0.0:
        raise ValueError("continuous shell reference received an invalid shell")

    collision_and_loss = angular_collision_loss_matrix(
        collisions,
        energy_index,
        mu_widths,
    ) + loss * np.eye(polar_cells)
    inverse_radial_streaming = np.diag(1.0 / (acceleration * mu))
    source_density = np.diag(
        1.0 / (2.0 * np.pi * mu_widths * shell_width)
    )

    # Columns are first the left-boundary density basis and then a constant
    # shell-source basis.  The last block accumulates the exact radial cell
    # integral for both bases.
    rows = 3 * polar_cells
    columns = 2 * polar_cells
    initial = np.zeros((rows, columns), dtype=float)
    initial[:polar_cells, :polar_cells] = np.eye(polar_cells)
    initial[
        polar_cells : 2 * polar_cells,
        polar_cells :,
    ] = np.eye(polar_cells)

    def derivative(speed_m_s: float, flat: np.ndarray) -> np.ndarray:
        state = flat.reshape(rows, columns)
        generator = -inverse_radial_streaming @ (
            collision_and_loss
            + spherical_angular_streaming_matrix(
                grid,
                speed_m_s,
                acceleration,
            )
        )
        result = np.zeros_like(state)
        result[:polar_cells] = (
            generator @ state[:polar_cells]
            + inverse_radial_streaming
            @ source_density
            @ state[polar_cells : 2 * polar_cells]
        )
        result[2 * polar_cells :] = state[:polar_cells]
        return result.reshape(-1)

    integration = solve_ivp(
        derivative,
        (left_speed, right_speed),
        initial.reshape(-1),
        method="DOP853",
        rtol=float(relative_tolerance),
        atol=float(absolute_tolerance),
        max_step=shell_width / 16.0,
    )
    if not integration.success:
        raise FloatingPointError(
            "continuous shell reference integration failed: "
            + str(integration.message)
        )
    final = integration.y[:, -1].reshape(rows, columns)
    transfer = final[:polar_cells, :polar_cells]
    source_transfer = final[:polar_cells, polar_cells:]
    integral = final[2 * polar_cells :, :polar_cells]
    source_integral = final[2 * polar_cells :, polar_cells:]

    t_pp = transfer[np.ix_(positive, positive)]
    t_pn = transfer[np.ix_(positive, negative)]
    t_np = transfer[np.ix_(negative, positive)]
    t_nn = transfer[np.ix_(negative, negative)]
    condition_number = float(np.linalg.cond(t_nn))
    if not np.isfinite(condition_number) or condition_number > 1.0e12:
        raise FloatingPointError(
            "continuous shell reference half-range elimination is ill-conditioned"
        )
    inverse_nn = np.linalg.inv(t_nn)
    source_positive = source_transfer[positive]
    source_negative = source_transfer[negative]

    left_negative_from_left_positive = -inverse_nn @ t_np
    left_negative_from_right_negative = inverse_nn
    left_negative_from_source = -inverse_nn @ source_negative
    right_positive_from_left_positive = (
        t_pp + t_pn @ left_negative_from_left_positive
    )
    right_positive_from_right_negative = (
        t_pn @ left_negative_from_right_negative
    )
    right_positive_from_source = (
        source_positive + t_pn @ left_negative_from_source
    )
    density_scattering = np.block(
        [
            [
                right_positive_from_left_positive,
                right_positive_from_right_negative,
            ],
            [
                left_negative_from_left_positive,
                left_negative_from_right_negative,
            ],
        ]
    )
    density_source = np.vstack(
        (right_positive_from_source, left_negative_from_source)
    )

    flux_scale = 2.0 * np.pi * acceleration * np.abs(mu) * mu_widths
    ordered_flux = flux_scale[order]
    scattering = (
        ordered_flux[:, None]
        * density_scattering
        / ordered_flux[None, :]
    )
    source_to_outflow = ordered_flux[:, None] * density_source

    left_density_from_incoming_flux = np.zeros(
        (polar_cells, polar_cells), dtype=float
    )
    left_density_from_incoming_flux[
        np.ix_(positive, np.arange(half))
    ] = np.diag(1.0 / ordered_flux[:half])
    left_density_from_incoming_flux[
        np.ix_(negative, np.arange(half))
    ] = left_negative_from_left_positive @ np.diag(
        1.0 / ordered_flux[:half]
    )
    left_density_from_incoming_flux[
        np.ix_(negative, np.arange(half, polar_cells))
    ] = left_negative_from_right_negative @ np.diag(
        1.0 / ordered_flux[half:]
    )
    left_density_from_source = np.zeros(
        (polar_cells, polar_cells), dtype=float
    )
    left_density_from_source[negative] = left_negative_from_source
    population_measure = 2.0 * np.pi * mu_widths[:, None]
    incoming_to_population = population_measure * (
        integral @ left_density_from_incoming_flux
    )
    source_to_population = population_measure * (
        integral @ left_density_from_source + source_integral
    )

    response = AffineShellResponse(
        scattering=scattering,
        source_to_outflow=source_to_outflow,
        incoming_to_population=incoming_to_population,
        source_to_population=source_to_population,
        boundary_order=order,
        conceptual_sublayers=0,
        minimum_entry=float(
            min(
                np.min(scattering),
                np.min(source_to_outflow),
                np.min(incoming_to_population),
                np.min(source_to_population),
            )
        ),
        maximum_column_balance_error=0.0,
    )
    return ContinuousShellReference(
        response=response,
        function_evaluations=int(integration.nfev),
        boundary_block_condition_number=condition_number,
    )


def affine_response_relative_distance(
    left: AffineShellResponse,
    right: AffineShellResponse,
) -> float:
    """Maximum normalized column-L1 distance over all affine maps."""

    errors: list[float] = []
    for left_map, right_map in (
        (left.scattering, right.scattering),
        (left.source_to_outflow, right.source_to_outflow),
        (left.incoming_to_population, right.incoming_to_population),
        (left.source_to_population, right.source_to_population),
    ):
        difference = float(
            np.max(np.sum(np.abs(left_map - right_map), axis=0), initial=0.0)
        )
        scale = max(
            float(np.max(np.sum(np.abs(left_map), axis=0), initial=0.0)),
            float(np.max(np.sum(np.abs(right_map), axis=0), initial=0.0)),
            1.0e-300,
        )
        errors.append(difference / scale)
    return max(errors, default=0.0)


def affine_origin_response_relative_distance(
    left: AffineOriginResponse,
    right: AffineOriginResponse,
) -> float:
    """Maximum normalized column-L1 distance over origin affine maps."""

    errors: list[float] = []
    for left_map, right_map in (
        (left.incoming_to_outflow, right.incoming_to_outflow),
        (left.source_to_outflow, right.source_to_outflow),
        (left.incoming_to_population, right.incoming_to_population),
        (left.source_to_population, right.source_to_population),
    ):
        difference = float(
            np.max(np.sum(np.abs(left_map - right_map), axis=0), initial=0.0)
        )
        scale = max(
            float(np.max(np.sum(np.abs(left_map), axis=0), initial=0.0)),
            float(np.max(np.sum(np.abs(right_map), axis=0), initial=0.0)),
            1.0e-300,
        )
        errors.append(difference / scale)
    return max(errors, default=0.0)


__all__ = [
    "ContinuousShellReference",
    "affine_origin_response_relative_distance",
    "affine_response_relative_distance",
    "continuous_shell_reference",
]
