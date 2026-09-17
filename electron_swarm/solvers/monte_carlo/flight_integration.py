"""Numerical integration along constant-acceleration, zero-field-B flights.

The null-collision clock still determines each flight endpoint.  These helpers
only integrate residence measures between those endpoints.  Histogram mass is
partitioned at the exact times at which the quadratic kinetic energy crosses a
bin edge, while reaction kernels are partitioned at cross-section knots and
physical thresholds before adaptive quadrature is applied.
"""

from __future__ import annotations

import math

import numpy as np

from electron_swarm.core.constants import ELECTRON_MASS_KG, EV_TO_J
from electron_swarm.solvers.monte_carlo.cross_section_table import NUMBA_AVAILABLE

if NUMBA_AVAILABLE:
    from numba import njit
else:  # pragma: no cover - used only when the optional accelerator is absent

    def njit(*_args, **_kwargs):  # type: ignore[no-untyped-def]
        def decorate(function):  # type: ignore[no-untyped-def]
            return function

        return decorate


FLIGHT_RATE_RELATIVE_TOLERANCE = 1.0e-8
FLIGHT_RATE_MAX_REFINEMENT_DEPTH = 12

_ENERGY_FACTOR = 0.5 * ELECTRON_MASS_KG / EV_TO_J
_SPEED_FACTOR = 2.0 * EV_TO_J / ELECTRON_MASS_KG
_FLOAT_TINY = np.finfo(np.float64).tiny
_FLOAT_EPSILON = np.finfo(np.float64).eps

_GAUSS_3_NODES = np.asarray(
    [-math.sqrt(3.0 / 5.0), 0.0, math.sqrt(3.0 / 5.0)],
    dtype=np.float64,
)
_GAUSS_3_WEIGHTS = np.asarray([5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0])
_GAUSS_5_NODES = np.asarray(
    [
        -0.9061798459386640,
        -0.5384693101056831,
        0.0,
        0.5384693101056831,
        0.9061798459386640,
    ],
    dtype=np.float64,
)
_GAUSS_5_WEIGHTS = np.asarray(
    [
        0.2369268850561891,
        0.4786286704993665,
        0.5688888888888889,
        0.4786286704993665,
        0.2369268850561891,
    ],
    dtype=np.float64,
)


@njit(cache=True, fastmath=False)
def _energy_moment_primitives(
    value: float,
    quadratic_energy: float,
    minimum_energy: float,
) -> tuple[float, float, float]:
    sqrt_quadratic = math.sqrt(quadratic_energy)
    if minimum_energy == 0.0:
        absolute = abs(value)
        first = 0.5 * sqrt_quadratic * value * absolute
        third = 0.25 * quadratic_energy * sqrt_quadratic * value * absolute**3
        fifth = (
            quadratic_energy
            * quadratic_energy
            * sqrt_quadratic
            * value
            * absolute**5
            / 6.0
        )
        return first, third, fifth
    energy = minimum_energy + quadratic_energy * value * value
    root_energy = math.sqrt(max(energy, 0.0))
    inverse_root = math.asinh(
        sqrt_quadratic * value / math.sqrt(minimum_energy)
    ) / sqrt_quadratic
    first = 0.5 * (value * root_energy + minimum_energy * inverse_root)
    third = 0.25 * value * energy * root_energy + 0.75 * minimum_energy * first
    fifth = (
        value * energy * energy * root_energy / 6.0
        + (5.0 / 6.0) * minimum_energy * third
    )
    return first, third, fifth


@njit(cache=True, fastmath=False)
def _analytic_energy_moments(
    left_fraction: float,
    right_fraction: float,
    before_velocity_m_s: np.ndarray,
    velocity_delta_m_s: np.ndarray,
) -> tuple[float, float, float, bool]:
    """Return exact half-integer energy moments when subtraction is stable."""

    dx = velocity_delta_m_s[0]
    dy = velocity_delta_m_s[1]
    dz = velocity_delta_m_s[2]
    delta_squared = dx * dx + dy * dy + dz * dz
    width = right_fraction - left_fraction
    if delta_squared == 0.0:
        energy = _energy_at_fraction(
            before_velocity_m_s, velocity_delta_m_s, 0.5 * (left_fraction + right_fraction)
        )
        root_energy = math.sqrt(max(energy, 0.0))
        return (
            width * root_energy,
            width * energy * root_energy,
            width * energy * energy * root_energy,
            True,
        )

    vertex = -float(np.dot(before_velocity_m_s, velocity_delta_m_s)) / delta_squared
    left = left_fraction - vertex
    right = right_fraction - vertex
    # Primitive subtraction loses roughly |x|/Delta-x digits.  Beyond this
    # bound the embedded quadrature path below is both safer and still cheap,
    # because such a distant vertex implies a slowly varying energy.
    if max(abs(left), abs(right)) > 1.0e5 * width:
        return 0.0, 0.0, 0.0, False

    vertex_vx = before_velocity_m_s[0] + vertex * dx
    vertex_vy = before_velocity_m_s[1] + vertex * dy
    vertex_vz = before_velocity_m_s[2] + vertex * dz
    minimum_energy = max(
        _ENERGY_FACTOR
        * (
            vertex_vx * vertex_vx
            + vertex_vy * vertex_vy
            + vertex_vz * vertex_vz
        ),
        0.0,
    )
    quadratic_energy = _ENERGY_FACTOR * delta_squared
    left_first, left_third, left_fifth = _energy_moment_primitives(
        left, quadratic_energy, minimum_energy
    )
    right_first, right_third, right_fifth = _energy_moment_primitives(
        right, quadratic_energy, minimum_energy
    )
    first = right_first - left_first
    third = right_third - left_third
    fifth = right_fifth - left_fifth
    valid = bool(
        math.isfinite(first)
        and math.isfinite(third)
        and math.isfinite(fifth)
        and first >= 0.0
        and third >= 0.0
        and fifth >= 0.0
    )
    return first, third, fifth, valid


@njit(cache=True, fastmath=False)
def _energy_at_fraction(
    before_velocity_m_s: np.ndarray,
    velocity_delta_m_s: np.ndarray,
    fraction: float,
) -> float:
    vx = before_velocity_m_s[0] + fraction * velocity_delta_m_s[0]
    vy = before_velocity_m_s[1] + fraction * velocity_delta_m_s[1]
    vz = before_velocity_m_s[2] + fraction * velocity_delta_m_s[2]
    return _ENERGY_FACTOR * (vx * vx + vy * vy + vz * vz)


@njit(cache=True, fastmath=False)
def _lower_bound(values: np.ndarray, target: float) -> int:
    lower = 0
    upper = values.size
    while lower < upper:
        middle = (lower + upper) // 2
        if values[middle] < target:
            lower = middle + 1
        else:
            upper = middle
    return lower


@njit(cache=True, fastmath=False)
def _upper_bound(values: np.ndarray, target: float) -> int:
    lower = 0
    upper = values.size
    while lower < upper:
        middle = (lower + upper) // 2
        if target < values[middle]:
            upper = middle
        else:
            lower = middle + 1
    return lower


@njit(cache=True, fastmath=False)
def _histogram_bin(edges_eV: np.ndarray, energy_eV: float) -> int:
    if energy_eV == edges_eV[edges_eV.size - 1]:
        return edges_eV.size - 2
    return _upper_bound(edges_eV, energy_eV) - 1


@njit(cache=True, fastmath=False)
def _crossing_fraction(
    before_velocity_m_s: np.ndarray,
    velocity_delta_m_s: np.ndarray,
    vertex_fraction: float,
    branch_midpoint: float,
    energy_eV: float,
) -> float:
    """Return the root of ``E(u) = energy_eV`` on one monotone branch."""

    vx = before_velocity_m_s[0]
    vy = before_velocity_m_s[1]
    vz = before_velocity_m_s[2]
    dx = velocity_delta_m_s[0]
    dy = velocity_delta_m_s[1]
    dz = velocity_delta_m_s[2]
    quadratic = dx * dx + dy * dy + dz * dz
    half_linear = vx * dx + vy * dy + vz * dz
    constant = vx * vx + vy * vy + vz * vz - energy_eV / _ENERGY_FACTOR
    discriminant = half_linear * half_linear - quadratic * constant
    scale = max(
        half_linear * half_linear,
        abs(quadratic * constant),
        _FLOAT_TINY,
    )
    if discriminant < 0.0 and discriminant >= -128.0 * _FLOAT_EPSILON * scale:
        discriminant = 0.0
    root_term = math.sqrt(max(discriminant, 0.0))
    if half_linear >= 0.0:
        numerator = -half_linear - root_term
    else:
        numerator = -half_linear + root_term
    if numerator == 0.0:
        first = -half_linear / quadratic
        second = first
    else:
        first = numerator / quadratic
        second = constant / numerator
    lower_root = min(first, second)
    upper_root = max(first, second)
    if branch_midpoint < vertex_fraction:
        return lower_root
    return upper_root


@njit(cache=True, fastmath=False)
def _touch_histogram_bin(
    bin_index: int,
    contribution: float,
    scratch: np.ndarray,
    touched_bins: np.ndarray,
    touched_count: int,
) -> int:
    if bin_index < 0 or bin_index >= scratch.size or contribution <= 0.0:
        return touched_count
    if scratch[bin_index] == 0.0:
        touched_bins[touched_count] = bin_index
        touched_count += 1
    scratch[bin_index] += contribution
    return touched_count


@njit(cache=True, fastmath=False)
def _accumulate_histogram_branch(
    before_velocity_m_s: np.ndarray,
    velocity_delta_m_s: np.ndarray,
    vertex_fraction: float,
    start_fraction: float,
    end_fraction: float,
    duration_s: float,
    weight: float,
    edges_eV: np.ndarray,
    scratch: np.ndarray,
    touched_bins: np.ndarray,
    touched_count: int,
) -> int:
    if end_fraction <= start_fraction:
        return touched_count
    start_energy = _energy_at_fraction(
        before_velocity_m_s, velocity_delta_m_s, start_fraction
    )
    end_energy = _energy_at_fraction(
        before_velocity_m_s, velocity_delta_m_s, end_fraction
    )
    energy_scale = max(abs(start_energy), abs(end_energy), 1.0)
    current = start_fraction
    branch_midpoint = 0.5 * (start_fraction + end_fraction)

    if abs(end_energy - start_energy) <= (
        64.0 * _FLOAT_EPSILON * energy_scale
    ):
        energy = _energy_at_fraction(
            before_velocity_m_s, velocity_delta_m_s, branch_midpoint
        )
        return _touch_histogram_bin(
            _histogram_bin(edges_eV, energy),
            weight * duration_s * (end_fraction - start_fraction),
            scratch,
            touched_bins,
            touched_count,
        )

    if end_energy > start_energy:
        boundary_index = _upper_bound(edges_eV, start_energy)
        while (
            boundary_index < edges_eV.size
            and edges_eV[boundary_index] < end_energy
        ):
            crossing = _crossing_fraction(
                before_velocity_m_s,
                velocity_delta_m_s,
                vertex_fraction,
                branch_midpoint,
                edges_eV[boundary_index],
            )
            crossing = min(max(crossing, current), end_fraction)
            if crossing > current:
                midpoint = 0.5 * (current + crossing)
                energy = _energy_at_fraction(
                    before_velocity_m_s, velocity_delta_m_s, midpoint
                )
                touched_count = _touch_histogram_bin(
                    _histogram_bin(edges_eV, energy),
                    weight * duration_s * (crossing - current),
                    scratch,
                    touched_bins,
                    touched_count,
                )
            current = crossing
            boundary_index += 1
    else:
        boundary_index = _lower_bound(edges_eV, start_energy) - 1
        while boundary_index >= 0 and edges_eV[boundary_index] > end_energy:
            crossing = _crossing_fraction(
                before_velocity_m_s,
                velocity_delta_m_s,
                vertex_fraction,
                branch_midpoint,
                edges_eV[boundary_index],
            )
            crossing = min(max(crossing, current), end_fraction)
            if crossing > current:
                midpoint = 0.5 * (current + crossing)
                energy = _energy_at_fraction(
                    before_velocity_m_s, velocity_delta_m_s, midpoint
                )
                touched_count = _touch_histogram_bin(
                    _histogram_bin(edges_eV, energy),
                    weight * duration_s * (crossing - current),
                    scratch,
                    touched_bins,
                    touched_count,
                )
            current = crossing
            boundary_index -= 1

    if end_fraction > current:
        midpoint = 0.5 * (current + end_fraction)
        energy = _energy_at_fraction(
            before_velocity_m_s, velocity_delta_m_s, midpoint
        )
        touched_count = _touch_histogram_bin(
            _histogram_bin(edges_eV, energy),
            weight * duration_s * (end_fraction - current),
            scratch,
            touched_bins,
            touched_count,
        )
    return touched_count


@njit(cache=True, fastmath=False)
def accumulate_dc_flight_histogram(
    before_velocity_m_s: np.ndarray,
    after_velocity_m_s: np.ndarray,
    duration_s: float,
    weight: float,
    edges_eV: np.ndarray,
    counts: np.ndarray,
    weighted_hist: np.ndarray,
    weighted_square_hist: np.ndarray,
    scratch: np.ndarray,
    touched_bins: np.ndarray,
) -> int:
    """Accumulate exact B=0 residence durations into energy bins."""

    velocity_delta = after_velocity_m_s - before_velocity_m_s
    quadratic = float(np.dot(velocity_delta, velocity_delta))
    if quadratic > 0.0:
        vertex = -float(np.dot(before_velocity_m_s, velocity_delta)) / quadratic
    else:
        vertex = -1.0
    split = min(max(vertex, 0.0), 1.0)
    touched_count = 0
    if split > 0.0:
        touched_count = _accumulate_histogram_branch(
            before_velocity_m_s,
            velocity_delta,
            vertex,
            0.0,
            split,
            duration_s,
            weight,
            edges_eV,
            scratch,
            touched_bins,
            touched_count,
        )
    if split < 1.0:
        touched_count = _accumulate_histogram_branch(
            before_velocity_m_s,
            velocity_delta,
            vertex,
            split,
            1.0,
            duration_s,
            weight,
            edges_eV,
            scratch,
            touched_bins,
            touched_count,
        )
    if touched_count == 0:
        energy = _energy_at_fraction(before_velocity_m_s, velocity_delta, 0.5)
        touched_count = _touch_histogram_bin(
            _histogram_bin(edges_eV, energy),
            weight * duration_s,
            scratch,
            touched_bins,
            touched_count,
        )

    for touched_index in range(touched_count):
        bin_index = touched_bins[touched_index]
        contribution = scratch[bin_index]
        counts[bin_index] += 1
        weighted_hist[bin_index] += contribution
        weighted_square_hist[bin_index] += contribution * contribution
        scratch[bin_index] = 0.0
    return touched_count


@njit(cache=True, fastmath=False)
def _evaluate_local_rate_rule(
    nodes: np.ndarray,
    weights: np.ndarray,
    left_fraction: float,
    right_fraction: float,
    before_velocity_m_s: np.ndarray,
    velocity_delta_m_s: np.ndarray,
    rate_grid_eV: np.ndarray,
    rate_values_m2: np.ndarray,
    rate_slopes_m2_eV: np.ndarray,
    rate_min_eV: np.ndarray,
    rate_max_eV: np.ndarray,
    rate_right_m2: np.ndarray,
    constant_loss_eV: np.ndarray,
    recoil_mass_ratio: np.ndarray,
    rate_out: np.ndarray,
    loss_out: np.ndarray,
) -> None:
    """Integrate rates in knot-local coordinates on a microinterval.

    Computing every node energy independently and then subtracting a nearby
    cross-section knot loses the small energy offset that defines a newly
    opened channel.  One reference energy plus the factorized quadratic
    difference keeps that offset coherent across the Gauss rules.
    """

    for process_index in range(rate_out.size):
        rate_out[process_index] = 0.0
        loss_out[process_index] = 0.0
    midpoint = 0.5 * (left_fraction + right_fraction)
    half_width = 0.5 * (right_fraction - left_fraction)
    reference_energy = _energy_at_fraction(
        before_velocity_m_s, velocity_delta_m_s, midpoint
    )
    velocity_dot_delta = 0.0
    delta_squared = 0.0
    for component in range(3):
        velocity_dot_delta += (
            before_velocity_m_s[component] * velocity_delta_m_s[component]
        )
        delta_squared += velocity_delta_m_s[component] ** 2

    if reference_energy < rate_grid_eV[0]:
        interval = -1
        anchor_index = 0
        reference_offset = 0.0
    elif reference_energy > rate_grid_eV[rate_grid_eV.size - 1]:
        interval = rate_grid_eV.size - 1
        anchor_index = rate_grid_eV.size - 1
        reference_offset = 0.0
    else:
        interval = min(
            _upper_bound(rate_grid_eV, reference_energy) - 1,
            rate_grid_eV.size - 2,
        )
        interval = max(interval, 0)
        if (
            reference_energy - rate_grid_eV[interval]
            <= rate_grid_eV[interval + 1] - reference_energy
        ):
            anchor_index = interval
        else:
            anchor_index = interval + 1
        reference_offset = reference_energy - rate_grid_eV[anchor_index]

    for node_index in range(nodes.size):
        fraction = midpoint + half_width * nodes[node_index]
        fraction_delta = fraction - midpoint
        energy_offset = _ENERGY_FACTOR * fraction_delta * (
            2.0 * velocity_dot_delta
            + (fraction + midpoint) * delta_squared
        )
        energy = max(reference_energy + energy_offset, 0.0)
        speed = math.sqrt(max(_SPEED_FACTOR * energy, 0.0))
        factor = half_width * weights[node_index]
        for process_index in range(rate_out.size):
            if reference_energy < rate_min_eV[process_index] or interval < 0:
                sigma = 0.0
            elif (
                reference_energy > rate_max_eV[process_index]
                or interval >= rate_grid_eV.size - 1
            ):
                sigma = rate_right_m2[process_index]
            else:
                sigma = (
                    rate_values_m2[process_index, anchor_index]
                    + rate_slopes_m2_eV[process_index, interval]
                    * (reference_offset + energy_offset)
                )
            rate = max(sigma, 0.0) * speed
            rate_out[process_index] += factor * rate
            fixed_loss = constant_loss_eV[process_index]
            recoil_ratio = recoil_mass_ratio[process_index]
            if math.isfinite(fixed_loss):
                loss_out[process_index] += factor * rate * fixed_loss
            elif math.isfinite(recoil_ratio):
                loss_out[process_index] += factor * rate * energy * recoil_ratio


@njit(cache=True, fastmath=False)
def _combine_rate_moments(
    representative_energy_eV: float,
    moment_half: float,
    moment_three_halves: float,
    moment_five_halves: float,
    rate_grid_eV: np.ndarray,
    rate_values_m2: np.ndarray,
    rate_slopes_m2_eV: np.ndarray,
    rate_min_eV: np.ndarray,
    rate_max_eV: np.ndarray,
    rate_right_m2: np.ndarray,
    constant_loss_eV: np.ndarray,
    recoil_mass_ratio: np.ndarray,
    rate_out: np.ndarray,
    loss_out: np.ndarray,
) -> None:
    """Combine common energy moments with each piecewise-linear sigma."""

    energy = representative_energy_eV
    if energy < rate_grid_eV[0]:
        interval = -1
    elif energy > rate_grid_eV[rate_grid_eV.size - 1]:
        interval = rate_grid_eV.size - 1
    else:
        interval = min(_upper_bound(rate_grid_eV, energy) - 1, rate_grid_eV.size - 2)
        interval = max(interval, 0)
    speed_factor = math.sqrt(_SPEED_FACTOR)
    for process_index in range(rate_out.size):
        if energy < rate_min_eV[process_index] or interval < 0:
            intercept = 0.0
            slope = 0.0
        elif energy > rate_max_eV[process_index] or interval >= rate_grid_eV.size - 1:
            intercept = rate_right_m2[process_index]
            slope = 0.0
        else:
            slope = rate_slopes_m2_eV[process_index, interval]
            intercept = (
                rate_values_m2[process_index, interval]
                - slope * rate_grid_eV[interval]
            )
        rate = speed_factor * (
            intercept * moment_half + slope * moment_three_halves
        )
        rate_out[process_index] = max(rate, 0.0)
        fixed_loss = constant_loss_eV[process_index]
        recoil_ratio = recoil_mass_ratio[process_index]
        if math.isfinite(fixed_loss):
            loss_out[process_index] = rate * fixed_loss
        elif math.isfinite(recoil_ratio):
            loss_out[process_index] = max(
                speed_factor
                * recoil_ratio
                * (intercept * moment_three_halves + slope * moment_five_halves),
                0.0,
            )
        else:
            loss_out[process_index] = 0.0


@njit(cache=True, fastmath=False)
def _integrate_rate_interval(
    left_fraction: float,
    right_fraction: float,
    before_velocity_m_s: np.ndarray,
    velocity_delta_m_s: np.ndarray,
    rate_grid_eV: np.ndarray,
    rate_values_m2: np.ndarray,
    rate_slopes_m2_eV: np.ndarray,
    rate_min_eV: np.ndarray,
    rate_max_eV: np.ndarray,
    rate_right_m2: np.ndarray,
    constant_loss_eV: np.ndarray,
    recoil_mass_ratio: np.ndarray,
    rate_total: np.ndarray,
    loss_total: np.ndarray,
    quadrature_workspace: np.ndarray,
    interval_workspace: np.ndarray,
    depth_workspace: np.ndarray,
) -> bool:
    stack_size = 1
    interval_workspace[0, 0] = left_fraction
    interval_workspace[1, 0] = right_fraction
    depth_workspace[0] = 0
    while stack_size > 0:
        stack_size -= 1
        left = interval_workspace[0, stack_size]
        right = interval_workspace[1, stack_size]
        depth = depth_workspace[stack_size]
        analytic_half, analytic_three_halves, analytic_five_halves, analytic = (
            _analytic_energy_moments(
                left,
                right,
                before_velocity_m_s,
                velocity_delta_m_s,
            )
        )
        if analytic:
            representative_energy = _energy_at_fraction(
                before_velocity_m_s,
                velocity_delta_m_s,
                0.5 * (left + right),
            )
            _combine_rate_moments(
                representative_energy,
                analytic_half,
                analytic_three_halves,
                analytic_five_halves,
                rate_grid_eV,
                rate_values_m2,
                rate_slopes_m2_eV,
                rate_min_eV,
                rate_max_eV,
                rate_right_m2,
                constant_loss_eV,
                recoil_mass_ratio,
                quadrature_workspace[1],
                quadrature_workspace[3],
            )
            for process_index in range(rate_total.size):
                rate_total[process_index] += quadrature_workspace[1, process_index]
                loss_total[process_index] += quadrature_workspace[3, process_index]
            continue
        _evaluate_local_rate_rule(
            _GAUSS_3_NODES,
            _GAUSS_3_WEIGHTS,
            left,
            right,
            before_velocity_m_s,
            velocity_delta_m_s,
            rate_grid_eV,
            rate_values_m2,
            rate_slopes_m2_eV,
            rate_min_eV,
            rate_max_eV,
            rate_right_m2,
            constant_loss_eV,
            recoil_mass_ratio,
            quadrature_workspace[0],
            quadrature_workspace[2],
        )
        _evaluate_local_rate_rule(
            _GAUSS_5_NODES,
            _GAUSS_5_WEIGHTS,
            left,
            right,
            before_velocity_m_s,
            velocity_delta_m_s,
            rate_grid_eV,
            rate_values_m2,
            rate_slopes_m2_eV,
            rate_min_eV,
            rate_max_eV,
            rate_right_m2,
            constant_loss_eV,
            recoil_mass_ratio,
            quadrature_workspace[1],
            quadrature_workspace[3],
        )
        converged = True
        for process_index in range(rate_total.size):
            low_rate = quadrature_workspace[0, process_index]
            high_rate = quadrature_workspace[1, process_index]
            rate_scale = max(abs(low_rate), abs(high_rate), _FLOAT_TINY)
            if abs(high_rate - low_rate) > (
                FLIGHT_RATE_RELATIVE_TOLERANCE * rate_scale
            ):
                converged = False
                break
            if math.isfinite(constant_loss_eV[process_index]) or math.isfinite(
                recoil_mass_ratio[process_index]
            ):
                low_loss = quadrature_workspace[2, process_index]
                high_loss = quadrature_workspace[3, process_index]
                loss_scale = max(abs(low_loss), abs(high_loss), _FLOAT_TINY)
                if abs(high_loss - low_loss) > (
                    FLIGHT_RATE_RELATIVE_TOLERANCE * loss_scale
                ):
                    converged = False
                    break
        if converged:
            for process_index in range(rate_total.size):
                rate_total[process_index] += quadrature_workspace[1, process_index]
                loss_total[process_index] += quadrature_workspace[3, process_index]
            continue
        if depth >= FLIGHT_RATE_MAX_REFINEMENT_DEPTH:
            return False
        midpoint = 0.5 * (left + right)
        interval_workspace[0, stack_size] = midpoint
        interval_workspace[1, stack_size] = right
        depth_workspace[stack_size] = depth + 1
        stack_size += 1
        interval_workspace[0, stack_size] = left
        interval_workspace[1, stack_size] = midpoint
        depth_workspace[stack_size] = depth + 1
        stack_size += 1
    return True


@njit(cache=True, fastmath=False)
def _integrate_rate_branch(
    start_fraction: float,
    end_fraction: float,
    before_velocity_m_s: np.ndarray,
    velocity_delta_m_s: np.ndarray,
    vertex_fraction: float,
    breakpoints_eV: np.ndarray,
    rate_grid_eV: np.ndarray,
    rate_values_m2: np.ndarray,
    rate_slopes_m2_eV: np.ndarray,
    rate_min_eV: np.ndarray,
    rate_max_eV: np.ndarray,
    rate_right_m2: np.ndarray,
    constant_loss_eV: np.ndarray,
    recoil_mass_ratio: np.ndarray,
    rate_total: np.ndarray,
    loss_total: np.ndarray,
    quadrature_workspace: np.ndarray,
    interval_workspace: np.ndarray,
    depth_workspace: np.ndarray,
) -> bool:
    if end_fraction <= start_fraction:
        return True
    start_energy = _energy_at_fraction(
        before_velocity_m_s, velocity_delta_m_s, start_fraction
    )
    end_energy = _energy_at_fraction(
        before_velocity_m_s, velocity_delta_m_s, end_fraction
    )
    current = start_fraction
    branch_midpoint = 0.5 * (start_fraction + end_fraction)
    energy_scale = max(abs(start_energy), abs(end_energy), 1.0)

    if abs(end_energy - start_energy) > (
        64.0 * _FLOAT_EPSILON * energy_scale
    ):
        if end_energy > start_energy:
            boundary_index = _upper_bound(breakpoints_eV, start_energy)
            while (
                boundary_index < breakpoints_eV.size
                and breakpoints_eV[boundary_index] < end_energy
            ):
                crossing = _crossing_fraction(
                    before_velocity_m_s,
                    velocity_delta_m_s,
                    vertex_fraction,
                    branch_midpoint,
                    breakpoints_eV[boundary_index],
                )
                crossing = min(max(crossing, current), end_fraction)
                if crossing > current and not _integrate_rate_interval(
                    current,
                    crossing,
                    before_velocity_m_s,
                    velocity_delta_m_s,
                    rate_grid_eV,
                    rate_values_m2,
                    rate_slopes_m2_eV,
                    rate_min_eV,
                    rate_max_eV,
                    rate_right_m2,
                    constant_loss_eV,
                    recoil_mass_ratio,
                    rate_total,
                    loss_total,
                    quadrature_workspace,
                    interval_workspace,
                    depth_workspace,
                ):
                    return False
                current = crossing
                boundary_index += 1
        else:
            boundary_index = _lower_bound(breakpoints_eV, start_energy) - 1
            while (
                boundary_index >= 0
                and breakpoints_eV[boundary_index] > end_energy
            ):
                crossing = _crossing_fraction(
                    before_velocity_m_s,
                    velocity_delta_m_s,
                    vertex_fraction,
                    branch_midpoint,
                    breakpoints_eV[boundary_index],
                )
                crossing = min(max(crossing, current), end_fraction)
                if crossing > current and not _integrate_rate_interval(
                    current,
                    crossing,
                    before_velocity_m_s,
                    velocity_delta_m_s,
                    rate_grid_eV,
                    rate_values_m2,
                    rate_slopes_m2_eV,
                    rate_min_eV,
                    rate_max_eV,
                    rate_right_m2,
                    constant_loss_eV,
                    recoil_mass_ratio,
                    rate_total,
                    loss_total,
                    quadrature_workspace,
                    interval_workspace,
                    depth_workspace,
                ):
                    return False
                current = crossing
                boundary_index -= 1

    if end_fraction > current:
        return _integrate_rate_interval(
            current,
            end_fraction,
            before_velocity_m_s,
            velocity_delta_m_s,
            rate_grid_eV,
            rate_values_m2,
            rate_slopes_m2_eV,
            rate_min_eV,
            rate_max_eV,
            rate_right_m2,
            constant_loss_eV,
            recoil_mass_ratio,
            rate_total,
            loss_total,
            quadrature_workspace,
            interval_workspace,
            depth_workspace,
        )
    return True


@njit(cache=True, fastmath=False)
def integrate_dc_flight_rates(
    before_velocity_m_s: np.ndarray,
    after_velocity_m_s: np.ndarray,
    duration_s: float,
    breakpoints_eV: np.ndarray,
    rate_grid_eV: np.ndarray,
    rate_values_m2: np.ndarray,
    rate_slopes_m2_eV: np.ndarray,
    rate_min_eV: np.ndarray,
    rate_max_eV: np.ndarray,
    rate_right_m2: np.ndarray,
    constant_loss_eV: np.ndarray,
    recoil_mass_ratio: np.ndarray,
    rate_total: np.ndarray,
    loss_total: np.ndarray,
    quadrature_workspace: np.ndarray,
    interval_workspace: np.ndarray,
    depth_workspace: np.ndarray,
) -> bool:
    """Integrate all ``sigma(E) v(E)`` kernels over one B=0 flight."""

    rate_total[:] = 0.0
    loss_total[:] = 0.0
    velocity_delta = after_velocity_m_s - before_velocity_m_s
    quadratic = float(np.dot(velocity_delta, velocity_delta))
    if quadratic > 0.0:
        vertex = -float(np.dot(before_velocity_m_s, velocity_delta)) / quadratic
    else:
        vertex = -1.0
    split = min(max(vertex, 0.0), 1.0)
    converged = True
    if split > 0.0:
        converged = _integrate_rate_branch(
            0.0,
            split,
            before_velocity_m_s,
            velocity_delta,
            vertex,
            breakpoints_eV,
            rate_grid_eV,
            rate_values_m2,
            rate_slopes_m2_eV,
            rate_min_eV,
            rate_max_eV,
            rate_right_m2,
            constant_loss_eV,
            recoil_mass_ratio,
            rate_total,
            loss_total,
            quadrature_workspace,
            interval_workspace,
            depth_workspace,
        )
    if converged and split < 1.0:
        converged = _integrate_rate_branch(
            split,
            1.0,
            before_velocity_m_s,
            velocity_delta,
            vertex,
            breakpoints_eV,
            rate_grid_eV,
            rate_values_m2,
            rate_slopes_m2_eV,
            rate_min_eV,
            rate_max_eV,
            rate_right_m2,
            constant_loss_eV,
            recoil_mass_ratio,
            rate_total,
            loss_total,
            quadrature_workspace,
            interval_workspace,
            depth_workspace,
        )
    if converged:
        rate_total *= duration_s
        loss_total *= duration_s
    return converged


def new_rate_integration_workspace(
    process_count: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Allocate reusable work arrays for :func:`integrate_dc_flight_rates`."""

    count = int(process_count)
    return (
        np.zeros(count, dtype=float),
        np.zeros(count, dtype=float),
        np.empty((4, count), dtype=float),
        np.empty((2, FLIGHT_RATE_MAX_REFINEMENT_DEPTH + 2), dtype=float),
        np.empty(FLIGHT_RATE_MAX_REFINEMENT_DEPTH + 2, dtype=np.int64),
    )
