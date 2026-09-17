"""Positive stationary scattering responses for non-origin speed shells.

Each shell solves radial acceleration, spherical angular streaming, elastic
angular relaxation, and local loss in one boundary-value problem.  A thin
diamond/Cayley response is scaled until its right-hand matrix is nonnegative,
then composed by Redheffer doubling.  Affine source and shell-population
responses are condensed alongside the boundary scattering matrix.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from electron_swarm.solvers.propagator.collisions import (
    angular_collision_loss_matrix,
)
from electron_swarm.solvers.propagator.models import (
    CollisionOperatorData,
    PropagatorGrid,
)


@dataclass(frozen=True, slots=True)
class AffineShellResponse:
    scattering: np.ndarray
    source_to_outflow: np.ndarray
    incoming_to_population: np.ndarray
    source_to_population: np.ndarray
    boundary_order: np.ndarray
    conceptual_sublayers: int
    minimum_entry: float
    maximum_column_balance_error: float
    rational_error_estimate: float = 0.0
    coefficient_error_estimate: float = 0.0
    coefficient_segments: int = 1

    @property
    def memory_bytes(self) -> int:
        return int(
            self.scattering.nbytes
            + self.source_to_outflow.nbytes
            + self.incoming_to_population.nbytes
            + self.source_to_population.nbytes
            + self.boundary_order.nbytes
        )


def half_range_boundary_order(mu: np.ndarray) -> np.ndarray:
    """Order both half ranges by increasing ``rho=1-mu**2``.

    Positive channels therefore run from ``mu=+1`` toward zero and negative
    channels from ``mu=-1`` toward zero.  This is the natural order of the
    characteristic origin ball.
    """

    values = np.asarray(mu, dtype=float)
    positive = np.flatnonzero(values > 0.0)
    negative = np.flatnonzero(values < 0.0)
    positive = positive[np.argsort(1.0 - values[positive] ** 2)]
    negative = negative[np.argsort(1.0 - values[negative] ** 2)]
    order = np.concatenate((positive, negative))
    if len(positive) != len(negative) or order.size != values.size:
        raise ValueError("stationary shell response requires an even polar grid")
    return order


def spherical_angular_streaming_matrix(
    grid: PropagatorGrid,
    speed_m_s: float,
    acceleration_m_s2: float,
) -> np.ndarray:
    """Return the conservative ``mu``-face divergence for ``p=v**2 f``."""

    speed = float(speed_m_s)
    acceleration = float(acceleration_m_s2)
    if speed <= 0.0 or acceleration <= 0.0:
        raise ValueError("regular shell streaming requires positive speed and field")
    widths = np.asarray(grid.mu_widths, dtype=float)
    mu_faces = np.cos(grid.theta_edges_rad)
    cells = grid.polar_cells
    matrix = np.zeros((cells, cells), dtype=float)
    for index in range(cells):
        high_rate = acceleration * (1.0 - mu_faces[index] ** 2) / speed
        low_rate = acceleration * (1.0 - mu_faces[index + 1] ** 2) / speed
        matrix[index, index] += high_rate / widths[index]
        if index + 1 < cells and low_rate > 0.0:
            matrix[index, index + 1] -= low_rate / widths[index]
    scale = max(float(np.linalg.norm(matrix, ord=1)), 1.0)
    if float(np.max(np.abs(widths @ matrix))) > 2.0e-12 * scale:
        raise FloatingPointError("spherical angular streaming is not conservative")
    off_diagonal = matrix.copy()
    np.fill_diagonal(off_diagonal, 0.0)
    if float(np.max(off_diagonal)) > 2.0e-14 * scale:
        raise FloatingPointError("spherical angular streaming is not an M-matrix")
    return matrix


def _base_affine_response(
    *,
    segment_width_m_s: float,
    shell_width_m_s: float,
    acceleration_m_s2: float,
    matrix_s_inv: np.ndarray,
    mu: np.ndarray,
    mu_widths: np.ndarray,
    boundary_order: np.ndarray,
) -> AffineShellResponse:
    width = float(segment_width_m_s)
    shell_width = float(shell_width_m_s)
    acceleration = float(acceleration_m_s2)
    radial_rate = acceleration * np.abs(mu) / width
    cells = radial_rate.size
    diagonal = np.diag_indices(cells)
    left = 0.5 * np.array(matrix_s_inv, dtype=float, copy=True)
    right = -0.5 * np.array(matrix_s_inv, dtype=float, copy=True)
    left[diagonal] += radial_rate
    right[diagonal] += radial_rate
    scale = max(float(np.linalg.norm(right, ord=np.inf)), 1.0)
    if float(np.min(right)) < -2.0e-13 * scale:
        raise FloatingPointError("thin shell response has a negative Cayley source")
    # Both maps have the same left-hand matrix.  Solve all right-hand sides
    # together so LAPACK factors the small dense matrix exactly once.
    right_hand = np.zeros((cells, 2 * cells), dtype=float)
    right_hand[:, :cells] = right
    right_hand[diagonal[0], cells + diagonal[1]] = (
        1.0 / (2.0 * np.pi * mu_widths * shell_width)
    )
    raw = np.linalg.solve(left, right_hand)
    raw_scattering = raw[:, :cells]
    raw_source = raw[:, cells:]
    flux_scale = 2.0 * np.pi * acceleration * np.abs(mu) * mu_widths
    ordered_flux = flux_scale[boundary_order]
    scattering = (
        ordered_flux[:, None]
        * raw_scattering[np.ix_(boundary_order, boundary_order)]
        / ordered_flux[None, :]
    )
    source_to_outflow = (
        ordered_flux[:, None] * raw_source[boundary_order, :]
    )

    inverse_order = np.empty_like(boundary_order)
    inverse_order[boundary_order] = np.arange(boundary_order.size)
    ordered_to_density = np.zeros(
        (boundary_order.size, boundary_order.size), dtype=float
    )
    ordered_to_density[boundary_order, np.arange(boundary_order.size)] = (
        1.0 / ordered_flux
    )
    population_scale = np.pi * width * mu_widths
    incoming_to_population = population_scale[:, None] * (
        ordered_to_density
        + ordered_to_density @ scattering
    )
    source_to_population = (
        population_scale[:, None]
        * ordered_to_density
        @ source_to_outflow
    )
    # ``inverse_order`` is constructed here as an explicit audit of the
    # permutation; retaining only ``boundary_order`` avoids duplicate state.
    if not np.array_equal(boundary_order[inverse_order], np.arange(len(mu))):
        raise FloatingPointError("invalid half-range boundary permutation")
    return AffineShellResponse(
        scattering=scattering,
        source_to_outflow=source_to_outflow,
        incoming_to_population=incoming_to_population,
        source_to_population=source_to_population,
        boundary_order=boundary_order,
        conceptual_sublayers=1,
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


def compose_affine_responses(
    left_response: AffineShellResponse,
    right_response: AffineShellResponse,
    *,
    audit: bool = True,
) -> AffineShellResponse:
    """Redheffer-compose adjacent responses with a common shell source."""

    left = left_response
    right = right_response
    cells = left.scattering.shape[0]
    if (
        cells % 2
        or right.scattering.shape != (cells, cells)
        or (
            left.boundary_order is not right.boundary_order
            and not np.array_equal(
                left.boundary_order,
                right.boundary_order,
            )
        )
    ):
        raise ValueError("incompatible affine shell responses")
    half = cells // 2
    a1 = left.scattering[:half, :half]
    b1 = left.scattering[:half, half:]
    c1 = left.scattering[half:, :half]
    d1 = left.scattering[half:, half:]
    a2 = right.scattering[:half, :half]
    b2 = right.scattering[:half, half:]
    c2 = right.scattering[half:, :half]
    d2 = right.scattering[half:, half:]
    interface = -(b1 @ c2)
    interface[np.diag_indices(half)] += 1.0

    # Boundary and affine-source interface maps share the same matrix.  The
    # literal Redheffer equations use two solves and several selector matrices;
    # concatenate the right-hand sides and express selectors as block updates.
    # This is algebraically identical, but removes a duplicate factorization
    # and the dominant temporary-array traffic in the shell builder.
    interface_rhs = np.empty((half, 2 * cells), dtype=float)
    interface_rhs[:, :half] = a1
    interface_rhs[:, half:cells] = b1 @ d2
    interface_rhs[:, cells:] = (
        left.source_to_outflow[:half]
        + b1 @ right.source_to_outflow[half:]
    )
    interface_maps = np.linalg.solve(interface, interface_rhs)
    interface_positive = interface_maps[:, :cells]
    source_positive = interface_maps[:, cells:]

    interface_negative = c2 @ interface_positive
    interface_negative[:, half:] += d2
    source_negative = (
        c2 @ source_positive + right.source_to_outflow[half:]
    )

    scattering = np.empty((cells, cells), dtype=float)
    scattering[:half] = a2 @ interface_positive
    scattering[:half, half:] += b2
    scattering[half:] = d1 @ interface_negative
    scattering[half:, :half] += c1

    source_to_outflow = np.empty((cells, cells), dtype=float)
    source_to_outflow[:half] = (
        a2 @ source_positive + right.source_to_outflow[:half]
    )
    source_to_outflow[half:] = (
        d1 @ source_negative + left.source_to_outflow[half:]
    )

    incoming_to_population = (
        left.incoming_to_population[:, half:] @ interface_negative
        + right.incoming_to_population[:, :half] @ interface_positive
    )
    incoming_to_population[:, :half] += (
        left.incoming_to_population[:, :half]
    )
    incoming_to_population[:, half:] += (
        right.incoming_to_population[:, half:]
    )
    source_to_population = (
        left.incoming_to_population[:, half:] @ source_negative
        + left.source_to_population
        + right.incoming_to_population[:, :half] @ source_positive
        + right.source_to_population
    )
    return AffineShellResponse(
        scattering=scattering,
        source_to_outflow=source_to_outflow,
        incoming_to_population=incoming_to_population,
        source_to_population=source_to_population,
        boundary_order=left.boundary_order,
        conceptual_sublayers=(
            left.conceptual_sublayers + right.conceptual_sublayers
        ),
        # Internal doubling/segmentation uses only the final composed response.
        # Deferring its four full-array reductions avoids rescanning every
        # discarded intermediate.  The public/default path remains audited.
        minimum_entry=(
            float(
                min(
                    np.min(scattering),
                    np.min(source_to_outflow),
                    np.min(incoming_to_population),
                    np.min(source_to_population),
                )
            )
            if audit
            else float("nan")
        ),
        maximum_column_balance_error=0.0,
        rational_error_estimate=max(
            left.rational_error_estimate,
            right.rational_error_estimate,
        ),
        coefficient_error_estimate=max(
            left.coefficient_error_estimate,
            right.coefficient_error_estimate,
        ),
        coefficient_segments=(
            left.coefficient_segments + right.coefficient_segments
        ),
    )


def _doubled_segment_response(
    *,
    segment_width_m_s: float,
    shell_width_m_s: float,
    acceleration_m_s2: float,
    matrix_s_inv: np.ndarray,
    mu: np.ndarray,
    mu_widths: np.ndarray,
    boundary_order: np.ndarray,
    additional_doublings: int,
) -> AffineShellResponse:
    required = float(
        np.max(
            np.diag(matrix_s_inv)
            * segment_width_m_s
            / (2.0 * acceleration_m_s2 * np.abs(mu))
        )
    )
    positivity_exponent = max(0, int(np.ceil(np.log2(max(required, 1.0)))))
    exponent = positivity_exponent + int(additional_doublings)
    if exponent > 24:
        raise MemoryError("stationary shell response exceeded its subdivision bound")
    response = _base_affine_response(
        segment_width_m_s=segment_width_m_s / (1 << exponent),
        shell_width_m_s=shell_width_m_s,
        acceleration_m_s2=acceleration_m_s2,
        matrix_s_inv=matrix_s_inv,
        mu=mu,
        mu_widths=mu_widths,
        boundary_order=boundary_order,
    )
    for _ in range(exponent):
        response = compose_affine_responses(response, response, audit=False)
    return response


SHELL_RATIONAL_ERROR_TOLERANCE = 2.0e-5
SHELL_COEFFICIENT_ERROR_TOLERANCE = 5.0e-4
SHELL_MAX_ADDITIONAL_DOUBLINGS = 9
SHELL_MAX_COEFFICIENT_SEGMENTS = 128


def _response_distance(
    left: AffineShellResponse,
    right: AffineShellResponse,
) -> float:
    errors: list[float] = []
    for coarse, fine in (
        (left.scattering, right.scattering),
        (left.source_to_outflow, right.source_to_outflow),
        (left.incoming_to_population, right.incoming_to_population),
        (left.source_to_population, right.source_to_population),
    ):
        difference = float(
            np.max(np.sum(np.abs(coarse - fine), axis=0), initial=0.0)
        )
        scale = max(
            float(np.max(np.sum(np.abs(fine), axis=0), initial=0.0)),
            float(np.max(np.sum(np.abs(coarse), axis=0), initial=0.0)),
            1.0e-300,
        )
        errors.append(difference / scale)
    return max(errors, default=0.0)


def _converged_segment_response(
    *,
    segment_width_m_s: float,
    shell_width_m_s: float,
    acceleration_m_s2: float,
    matrix_s_inv: np.ndarray,
    mu: np.ndarray,
    mu_widths: np.ndarray,
    boundary_order: np.ndarray,
) -> AffineShellResponse:
    coarse = _doubled_segment_response(
        segment_width_m_s=segment_width_m_s,
        shell_width_m_s=shell_width_m_s,
        acceleration_m_s2=acceleration_m_s2,
        matrix_s_inv=matrix_s_inv,
        mu=mu,
        mu_widths=mu_widths,
        boundary_order=boundary_order,
        additional_doublings=2,
    )
    for level in range(3, SHELL_MAX_ADDITIONAL_DOUBLINGS + 1):
        fine = _doubled_segment_response(
            segment_width_m_s=segment_width_m_s,
            shell_width_m_s=shell_width_m_s,
            acceleration_m_s2=acceleration_m_s2,
            matrix_s_inv=matrix_s_inv,
            mu=mu,
            mu_widths=mu_widths,
            boundary_order=boundary_order,
            additional_doublings=level,
        )
        estimate = _response_distance(coarse, fine) / 3.0
        if estimate <= SHELL_RATIONAL_ERROR_TOLERANCE:
            return replace(fine, rational_error_estimate=estimate)
        coarse = fine
    raise FloatingPointError(
        "stationary shell Cayley response did not meet its bounded "
        "Richardson tolerance"
    )


def _segmented_shell_response(
    *,
    grid: PropagatorGrid,
    angular_collision: np.ndarray,
    left_speed_m_s: float,
    shell_width_m_s: float,
    acceleration_m_s2: float,
    local_scalar_loss_s_inv: float,
    coefficient_segments: int,
    mu: np.ndarray,
    mu_widths: np.ndarray,
    boundary_order: np.ndarray,
) -> AffineShellResponse:
    responses: list[AffineShellResponse] = []
    for segment in range(coefficient_segments):
        segment_left = (
            left_speed_m_s
            + segment * shell_width_m_s / coefficient_segments
        )
        segment_right = (
            left_speed_m_s
            + (segment + 1) * shell_width_m_s / coefficient_segments
        )
        segment_width = segment_right - segment_left
        # The only radial coefficient in the shell matrix is the spherical
        # streaming factor 1/v.  Use its exact cell average rather than a
        # midpoint sample; the remaining refinement error is therefore the
        # genuine non-commutativity with collisions/local loss.
        inverse_speed_average = (
            np.log(segment_right / segment_left) / segment_width
        )
        speed = 1.0 / inverse_speed_average
        matrix = (
            spherical_angular_streaming_matrix(
                grid,
                speed,
                acceleration_m_s2,
            )
            + angular_collision
            + local_scalar_loss_s_inv * np.eye(grid.polar_cells)
        )
        responses.append(
            _converged_segment_response(
                segment_width_m_s=shell_width_m_s / coefficient_segments,
                shell_width_m_s=shell_width_m_s,
                acceleration_m_s2=acceleration_m_s2,
                matrix_s_inv=matrix,
                mu=mu,
                mu_widths=mu_widths,
                boundary_order=boundary_order,
            )
        )
    response = responses[0]
    for next_response in responses[1:]:
        response = compose_affine_responses(
            response,
            next_response,
            audit=False,
        )
    return replace(
        response,
        rational_error_estimate=max(
            item.rational_error_estimate for item in responses
        ),
        coefficient_segments=coefficient_segments,
    )


def build_shell_response(
    grid: PropagatorGrid,
    collisions: CollisionOperatorData,
    energy_index: int,
    acceleration_m_s2: float,
    local_scalar_loss_s_inv: float,
) -> AffineShellResponse:
    """Build one variable-coefficient shell response.

    Cell-integrated coefficient segmentation and the Cayley rational
    approximation are refined independently until their embedded second-order
    estimates are bounded.  The controls are numerical-core constants rather
    than public solver knobs.
    """

    if energy_index <= 0 or energy_index >= grid.energy_cells:
        raise ValueError("regular shell response excludes the velocity origin")
    acceleration = float(acceleration_m_s2)
    if acceleration <= 0.0 or not np.isfinite(acceleration):
        raise ValueError("stationary shell response requires positive acceleration")
    loss = float(local_scalar_loss_s_inv)
    if loss < 0.0 or not np.isfinite(loss):
        raise ValueError("stationary shell response requires nonnegative local loss")
    left_speed = float(grid.speed_edges_m_s[energy_index])
    right_speed = float(grid.speed_edges_m_s[energy_index + 1])
    shell_width = right_speed - left_speed
    if left_speed <= 0.0 or shell_width <= 0.0:
        raise ValueError("invalid regular speed shell")
    mu = np.asarray(grid.mu_centers, dtype=float)
    widths = np.asarray(grid.mu_widths, dtype=float)
    order = half_range_boundary_order(mu)
    angular_collision = angular_collision_loss_matrix(
        collisions,
        energy_index,
        widths,
    )
    coarsest = _segmented_shell_response(
        grid=grid,
        angular_collision=angular_collision,
        left_speed_m_s=left_speed,
        shell_width_m_s=shell_width,
        acceleration_m_s2=acceleration,
        local_scalar_loss_s_inv=loss,
        coefficient_segments=1,
        mu=mu,
        mu_widths=widths,
        boundary_order=order,
    )
    coarse = _segmented_shell_response(
        grid=grid,
        angular_collision=angular_collision,
        left_speed_m_s=left_speed,
        shell_width_m_s=shell_width,
        acceleration_m_s2=acceleration,
        local_scalar_loss_s_inv=loss,
        coefficient_segments=2,
        mu=mu,
        mu_widths=widths,
        boundary_order=order,
    )
    previous_difference = _response_distance(coarsest, coarse)
    response: AffineShellResponse | None = None
    coefficient_segments = 4
    while coefficient_segments <= SHELL_MAX_COEFFICIENT_SEGMENTS:
        fine = _segmented_shell_response(
            grid=grid,
            angular_collision=angular_collision,
            left_speed_m_s=left_speed,
            shell_width_m_s=shell_width,
            acceleration_m_s2=acceleration,
            local_scalar_loss_s_inv=loss,
            coefficient_segments=coefficient_segments,
            mu=mu,
            mu_widths=widths,
            boundary_order=order,
        )
        difference = _response_distance(coarse, fine)
        if previous_difference <= np.finfo(float).eps:
            contraction = 0.25
        else:
            contraction = difference / previous_difference
        # A single h/2 difference cannot establish that the nominal
        # second-order midpoint segmentation is in its asymptotic regime.
        # Three successive levels expose the observed contraction q.  The
        # remaining geometric error at the accepted fine level is
        # q/(1-q) times the last difference.  Never claim faster than the
        # designed second order (q=1/4), and require a genuinely contracting
        # sequence before accepting the response.
        effective_contraction = max(float(contraction), 0.25)
        if effective_contraction < 0.8:
            estimate = (
                difference
                * effective_contraction
                / (1.0 - effective_contraction)
            )
        else:
            estimate = float("inf")
        if estimate <= SHELL_COEFFICIENT_ERROR_TOLERANCE:
            response = replace(fine, coefficient_error_estimate=estimate)
            break
        previous_difference = difference
        coarse = fine
        coefficient_segments *= 2
    if response is None:
        raise FloatingPointError(
            "stationary shell coefficient response did not meet its bounded "
            "Richardson tolerance"
        )
    scale = max(
        float(np.max(response.scattering)),
        float(np.max(response.source_to_outflow)),
        float(np.max(response.incoming_to_population)),
        float(np.max(response.source_to_population)),
        1.0,
    )
    minimum_entry = float(
        min(
            np.min(response.scattering),
            np.min(response.source_to_outflow),
            np.min(response.incoming_to_population),
            np.min(response.source_to_population),
        )
    )
    response = replace(response, minimum_entry=minimum_entry)
    if response.minimum_entry < -2.0e-11 * scale:
        raise FloatingPointError("stationary shell response is not positive")
    return response
