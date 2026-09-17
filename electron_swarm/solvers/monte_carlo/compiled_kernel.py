"""Compiled numerical barrier kernel for the weighted Monte Carlo backend.

The kernel owns only the dense numerical loop inside one synchronized time
barrier.  Population resampling, transport planes, feature dispatch, metadata,
and all validation remain in the Python orchestration layer.
"""

from __future__ import annotations

import math

import numpy as np

from electron_swarm.core.constants import ELECTRON_MASS_KG, EV_TO_J
from electron_swarm.solvers.monte_carlo.cross_section_table import (
    NUMBA_AVAILABLE,
    evaluate_prepared_cross_sections,
)
from electron_swarm.solvers.monte_carlo.flight_integration import (
    FLIGHT_RATE_MAX_REFINEMENT_DEPTH,
    accumulate_dc_flight_histogram,
    integrate_dc_flight_rates,
)

if NUMBA_AVAILABLE:
    from numba import njit
else:  # pragma: no cover - used only when the optional accelerator is absent

    def njit(*_args, **_kwargs):  # type: ignore[no-untyped-def]
        def decorate(function):  # type: ignore[no-untyped-def]
            return function

        return decorate


COMPILED_MC_KERNEL_SCHEMA_VERSION = "numba_weighted_dc_b0_thermal_isotropic.v4"

PROCESS_NONE = 0
PROCESS_SCATTERING = 1
PROCESS_EXCITATION = 2
PROCESS_IONIZATION = 3
PROCESS_SUPERELASTIC = 4

IONIZATION_EQUAL = 1
IONIZATION_PRIMARY_SECONDARY = 2

KERNEL_OK = 0
KERNEL_MAJORANT_EXCEEDED = 1
KERNEL_CROSS_SECTION_RANGE_EXCEEDED = 2
KERNEL_ENERGY_LIMIT_EXCEEDED = 3
KERNEL_ERROR_EXTRAPOLATION = 4
KERNEL_INVALID_RESIDENCE_SAMPLE = 5
KERNEL_FLIGHT_INTEGRATION_FAILED = 6

_ENERGY_FACTOR = 0.5 * ELECTRON_MASS_KG / EV_TO_J
_SPEED_FACTOR = 2.0 * EV_TO_J / ELECTRON_MASS_KG
_TWO_PI = 2.0 * math.pi
_FLOAT_TINY = 2.2250738585072014e-308


@njit(cache=True, fastmath=False)
def _maxwellian_target_velocity(
    component_std_m_s: float,
    generator,
) -> tuple[float, float, float]:
    return (
        generator.normal(0.0, component_std_m_s),
        generator.normal(0.0, component_std_m_s),
        generator.normal(0.0, component_std_m_s),
    )


@njit(cache=True, fastmath=False)
def _speed_weighted_maxwellian_target_velocity(
    component_std_m_s: float,
    generator,
) -> tuple[float, float, float]:
    first = max(generator.random(), _FLOAT_TINY)
    second = max(generator.random(), _FLOAT_TINY)
    radial_gamma = -math.log(first) - math.log(second)
    speed = component_std_m_s * math.sqrt(2.0 * radial_gamma)
    dx, dy, dz = _random_direction(generator)
    return speed * dx, speed * dy, speed * dz


@njit(cache=True, fastmath=False)
def _energy_eV(vx: float, vy: float, vz: float) -> float:
    return _ENERGY_FACTOR * (vx * vx + vy * vy + vz * vz)


@njit(cache=True, fastmath=False)
def _random_direction(generator) -> tuple[float, float, float]:
    mu = generator.uniform(-1.0, 1.0)
    phi = generator.uniform(0.0, _TWO_PI)
    sint = math.sqrt(max(1.0 - mu * mu, 0.0))
    return sint * math.cos(phi), sint * math.sin(phi), mu


@njit(cache=True, fastmath=False)
def _elastic_binary_scatter(
    vx: float,
    vy: float,
    vz: float,
    target_vx: float,
    target_vy: float,
    target_vz: float,
    mu: float,
    target_mass_fraction: float,
    generator,
) -> tuple[float, float, float]:
    """Exact two-body elastic update in the center-of-mass frame."""

    gx = vx - target_vx
    gy = vy - target_vy
    gz = vz - target_vz
    relative_speed = math.sqrt(gx * gx + gy * gy + gz * gz)
    electron_mass_fraction = 1.0 - target_mass_fraction
    center_x = electron_mass_fraction * vx + target_mass_fraction * target_vx
    center_y = electron_mass_fraction * vy + target_mass_fraction * target_vy
    center_z = electron_mass_fraction * vz + target_mass_fraction * target_vz
    if relative_speed <= 0.0:
        return center_x, center_y, center_z

    ax = gx / relative_speed
    ay = gy / relative_speed
    az = gz / relative_speed
    if abs(ax) < 0.9:
        rx, ry, rz = 1.0, 0.0, 0.0
    else:
        rx, ry, rz = 0.0, 1.0, 0.0
    e1x = ay * rz - az * ry
    e1y = az * rx - ax * rz
    e1z = ax * ry - ay * rx
    e1_norm = max(math.sqrt(e1x * e1x + e1y * e1y + e1z * e1z), 1.0e-300)
    e1x /= e1_norm
    e1y /= e1_norm
    e1z /= e1_norm
    e2x = ay * e1z - az * e1y
    e2y = az * e1x - ax * e1z
    e2z = ax * e1y - ay * e1x
    phi = generator.uniform(0.0, _TWO_PI)
    sine = math.sqrt(max(1.0 - mu * mu, 0.0))
    transverse_1 = sine * math.cos(phi)
    transverse_2 = sine * math.sin(phi)
    rotated_x = relative_speed * (mu * ax + transverse_1 * e1x + transverse_2 * e2x)
    rotated_y = relative_speed * (mu * ay + transverse_1 * e1y + transverse_2 * e2y)
    rotated_z = relative_speed * (mu * az + transverse_1 * e1z + transverse_2 * e2z)
    return (
        center_x + target_mass_fraction * rotated_x,
        center_y + target_mass_fraction * rotated_y,
        center_z + target_mass_fraction * rotated_z,
    )


@njit(cache=True, fastmath=False)
def _grow_particle_storage(
    positions: np.ndarray,
    velocities: np.ndarray,
    times: np.ndarray,
    weights: np.ndarray,
    lineages: np.ndarray,
    size: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    capacity = max(1, 2 * weights.size)
    new_positions = np.empty((capacity, 3), dtype=np.float64)
    new_velocities = np.empty((capacity, 3), dtype=np.float64)
    new_times = np.empty(capacity, dtype=np.float64)
    new_weights = np.empty(capacity, dtype=np.float64)
    new_lineages = np.empty(capacity, dtype=np.int64)
    for index in range(size):
        new_positions[index, 0] = positions[index, 0]
        new_positions[index, 1] = positions[index, 1]
        new_positions[index, 2] = positions[index, 2]
        new_velocities[index, 0] = velocities[index, 0]
        new_velocities[index, 1] = velocities[index, 1]
        new_velocities[index, 2] = velocities[index, 2]
        new_times[index] = times[index]
        new_weights[index] = weights[index]
        new_lineages[index] = lineages[index]
    return new_positions, new_velocities, new_times, new_weights, new_lineages


@njit(cache=True, fastmath=False)
def advance_weighted_barrier(
    positions: np.ndarray,
    velocities: np.ndarray,
    times: np.ndarray,
    weights: np.ndarray,
    lineages: np.ndarray,
    size: int,
    barrier_time_s: float,
    trial_frequency_s_inv: float,
    density_m3: float,
    acceleration_m_s2: np.ndarray,
    collision_grid_eV: np.ndarray,
    collision_values_m2: np.ndarray,
    collision_slopes_m2_eV: np.ndarray,
    collision_min_eV: np.ndarray,
    collision_max_eV: np.ndarray,
    collision_right_m2: np.ndarray,
    collision_error_extrapolation: np.ndarray,
    process_collision_group: np.ndarray,
    cold_bound_m3_s: float,
    thermal_group_bounds_m3_s: np.ndarray,
    thermal_group_component_std_m_s: np.ndarray,
    thermal_group_mean_speed_m_s: np.ndarray,
    thermal_group_proposal_mode: np.ndarray,
    process_kind: np.ndarray,
    process_threshold_eV: np.ndarray,
    process_target_mass_fraction: np.ndarray,
    ionization_model: int,
    secondary_electron_energy_eV: float,
    sampling_enabled: bool,
    transport_enabled: bool,
    reaction_block_index: int,
    histogram_edges_eV: np.ndarray,
    histogram_counts: np.ndarray,
    histogram_weighted: np.ndarray,
    histogram_weighted_square: np.ndarray,
    rate_breakpoints_eV: np.ndarray,
    rate_grid_eV: np.ndarray,
    rate_values_m2: np.ndarray,
    rate_slopes_m2_eV: np.ndarray,
    rate_min_eV: np.ndarray,
    rate_max_eV: np.ndarray,
    rate_right_m2: np.ndarray,
    rate_error_extrapolation: np.ndarray,
    rate_constant_loss_eV: np.ndarray,
    rate_recoil_mass_ratio: np.ndarray,
    rate_weighted_time: np.ndarray,
    rate_integrals: np.ndarray,
    rate_loss_integrals: np.ndarray,
    rate_event_counts: np.ndarray,
    rate_weighted_event_counts: np.ndarray,
    rate_event_loss_sample_counts: np.ndarray,
    rate_weighted_event_loss_eV: np.ndarray,
    rate_weighted_event_counts_by_block: np.ndarray,
    rate_weighted_event_loss_eV_by_block: np.ndarray,
    zero_high_energy_policy: bool,
    max_cross_section_energy_eV: float,
    max_energy_limit_eV: float,
    generator,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    int,
    float,
    np.ndarray,
    int,
    int,
    float,
    int,
]:
    """Advance every active particle to one common time barrier."""

    collision_sigma = np.empty(collision_values_m2.shape[0], dtype=np.float64)
    process_count = rate_values_m2.shape[0]
    flight_before_velocity = np.empty(3, dtype=np.float64)
    flight_after_velocity = np.empty(3, dtype=np.float64)
    histogram_scratch = np.zeros(histogram_counts.size, dtype=np.float64)
    histogram_touched = np.empty(histogram_counts.size, dtype=np.int64)
    flight_rate_integrals = np.zeros(process_count, dtype=np.float64)
    flight_loss_integrals = np.zeros(process_count, dtype=np.float64)
    quadrature_workspace = np.empty((4, process_count), dtype=np.float64)
    interval_workspace = np.empty(
        (2, FLIGHT_RATE_MAX_REFINEMENT_DEPTH + 2), dtype=np.float64
    )
    depth_workspace = np.empty(
        FLIGHT_RATE_MAX_REFINEMENT_DEPTH + 2, dtype=np.int64
    )
    # weighted time, displacement xyz, energy*time, energy*displacement xyz
    residence_totals = np.zeros(8, dtype=np.float64)
    event_observation_time_s = 0.0
    secondary_events = 0
    particle = 0
    while particle < size:
        tolerance = max(1.0e-12 * barrier_time_s, 1.0e-300)
        while times[particle] < barrier_time_s - tolerance:
            remaining = barrier_time_s - times[particle]
            sampled_dt = generator.exponential(1.0 / trial_frequency_s_inv)
            collision_due = sampled_dt < remaining
            dt = sampled_dt if collision_due else remaining
            weight = weights[particle]

            before_vx = velocities[particle, 0]
            before_vy = velocities[particle, 1]
            before_vz = velocities[particle, 2]
            after_vx = before_vx + acceleration_m_s2[0] * dt
            after_vy = before_vy + acceleration_m_s2[1] * dt
            after_vz = before_vz + acceleration_m_s2[2] * dt
            midpoint_vx = 0.5 * (before_vx + after_vx)
            midpoint_vy = 0.5 * (before_vy + after_vy)
            midpoint_vz = 0.5 * (before_vz + after_vz)
            displacement_x = midpoint_vx * dt
            displacement_y = midpoint_vy * dt
            displacement_z = midpoint_vz * dt
            positions[particle, 0] += displacement_x
            positions[particle, 1] += displacement_y
            positions[particle, 2] += displacement_z
            velocities[particle, 0] = after_vx
            velocities[particle, 1] = after_vy
            velocities[particle, 2] = after_vz
            times[particle] += dt
            final_energy = max(_energy_eV(after_vx, after_vy, after_vz), 1.0e-6)
            if zero_high_energy_policy and final_energy > max_cross_section_energy_eV:
                return (
                    positions,
                    velocities,
                    times,
                    weights,
                    lineages,
                    size,
                    event_observation_time_s,
                    residence_totals,
                    secondary_events,
                    KERNEL_CROSS_SECTION_RANGE_EXCEEDED,
                    final_energy,
                    -1,
                )
            if final_energy > max_energy_limit_eV:
                return (
                    positions,
                    velocities,
                    times,
                    weights,
                    lineages,
                    size,
                    event_observation_time_s,
                    residence_totals,
                    secondary_events,
                    KERNEL_ENERGY_LIMIT_EXCEEDED,
                    final_energy,
                    -1,
                )

            if sampling_enabled:
                if transport_enabled:
                    residence_totals[0] += weight * dt
                    residence_totals[1] += weight * displacement_x
                    residence_totals[2] += weight * displacement_y
                    residence_totals[3] += weight * displacement_z
                velocity_delta_x = after_vx - before_vx
                velocity_delta_y = after_vy - before_vy
                velocity_delta_z = after_vz - before_vz
                if transport_enabled:
                    root = math.sqrt(3.0 / 5.0)
                    for quadrature_index in range(3):
                        if quadrature_index == 0:
                            fraction = 0.5 * (1.0 - root)
                            time_fraction = 5.0 / 18.0
                        elif quadrature_index == 1:
                            fraction = 0.5
                            time_fraction = 4.0 / 9.0
                        else:
                            fraction = 0.5 * (1.0 + root)
                            time_fraction = 5.0 / 18.0
                        sample_vx = before_vx + fraction * velocity_delta_x
                        sample_vy = before_vy + fraction * velocity_delta_y
                        sample_vz = before_vz + fraction * velocity_delta_z
                        sample_energy = _energy_eV(sample_vx, sample_vy, sample_vz)
                        if not math.isfinite(sample_energy) or sample_energy <= 0.0:
                            return (
                                positions,
                                velocities,
                                times,
                                weights,
                                lineages,
                                size,
                                event_observation_time_s,
                                residence_totals,
                                secondary_events,
                                KERNEL_INVALID_RESIDENCE_SAMPLE,
                                sample_energy,
                                -1,
                            )
                        sample_dt = time_fraction * dt
                        contribution = weight * sample_dt
                        residence_totals[4] += contribution * sample_energy
                        residence_totals[5] += (
                            weight * sample_energy * sample_vx * sample_dt
                        )
                        residence_totals[6] += (
                            weight * sample_energy * sample_vy * sample_dt
                        )
                        residence_totals[7] += (
                            weight * sample_energy * sample_vz * sample_dt
                        )

                flight_before_velocity[0] = before_vx
                flight_before_velocity[1] = before_vy
                flight_before_velocity[2] = before_vz
                flight_after_velocity[0] = after_vx
                flight_after_velocity[1] = after_vy
                flight_after_velocity[2] = after_vz
                maximum_flight_energy = max(
                    _energy_eV(before_vx, before_vy, before_vz),
                    _energy_eV(after_vx, after_vy, after_vz),
                )
                for process_index in range(rate_error_extrapolation.size):
                    if (
                        rate_error_extrapolation[process_index]
                        and maximum_flight_energy > rate_max_eV[process_index]
                    ):
                        return (
                            positions,
                            velocities,
                            times,
                            weights,
                            lineages,
                            size,
                            event_observation_time_s,
                            residence_totals,
                            secondary_events,
                            KERNEL_ERROR_EXTRAPOLATION,
                            maximum_flight_energy,
                            process_index,
                        )
                if not integrate_dc_flight_rates(
                    flight_before_velocity,
                    flight_after_velocity,
                    dt,
                    rate_breakpoints_eV,
                    rate_grid_eV,
                    rate_values_m2,
                    rate_slopes_m2_eV,
                    rate_min_eV,
                    rate_max_eV,
                    rate_right_m2,
                    rate_constant_loss_eV,
                    rate_recoil_mass_ratio,
                    flight_rate_integrals,
                    flight_loss_integrals,
                    quadrature_workspace,
                    interval_workspace,
                    depth_workspace,
                ):
                    return (
                        positions,
                        velocities,
                        times,
                        weights,
                        lineages,
                        size,
                        event_observation_time_s,
                        residence_totals,
                        secondary_events,
                        KERNEL_FLIGHT_INTEGRATION_FAILED,
                        maximum_flight_energy,
                        -1,
                    )
                rate_weighted_time[reaction_block_index] += weight * dt
                for process_index in range(process_count):
                    rate_integrals[process_index, reaction_block_index] += (
                        weight * flight_rate_integrals[process_index]
                    )
                    rate_loss_integrals[process_index, reaction_block_index] += (
                        weight * flight_loss_integrals[process_index]
                    )
                accumulate_dc_flight_histogram(
                    flight_before_velocity,
                    flight_after_velocity,
                    dt,
                    weight,
                    histogram_edges_eV,
                    histogram_counts,
                    histogram_weighted,
                    histogram_weighted_square,
                    histogram_scratch,
                    histogram_touched,
                )
                event_observation_time_s += dt

            if collision_due:
                trial_rate_coefficient = trial_frequency_s_inv / density_m3
                mark = generator.random() * trial_rate_coefficient
                collision_group = -2
                group_bound = 0.0
                if mark < cold_bound_m3_s:
                    collision_group = -1
                    group_bound = cold_bound_m3_s
                else:
                    mark -= cold_bound_m3_s
                    for group_index in range(thermal_group_bounds_m3_s.size):
                        bound = thermal_group_bounds_m3_s[group_index]
                        if mark < bound:
                            collision_group = group_index
                            group_bound = bound
                            break
                        mark -= bound
                if collision_group == -2:
                    continue

                speed = math.sqrt(
                    after_vx * after_vx + after_vy * after_vy + after_vz * after_vz
                )
                target_vx = 0.0
                target_vy = 0.0
                target_vz = 0.0
                collision_energy = final_energy
                if collision_group == -1:
                    for process_index in range(collision_error_extrapolation.size):
                        if (
                            process_collision_group[process_index] == -1
                            and collision_error_extrapolation[process_index]
                            and final_energy > collision_max_eV[process_index]
                        ):
                            return (
                                positions,
                                velocities,
                                times,
                                weights,
                                lineages,
                                size,
                                event_observation_time_s,
                                residence_totals,
                                secondary_events,
                                KERNEL_ERROR_EXTRAPOLATION,
                                final_energy,
                                process_index,
                            )
                    physical_rate_coefficient = speed
                else:
                    component_std = thermal_group_component_std_m_s[collision_group]
                    mean_speed = thermal_group_mean_speed_m_s[collision_group]
                    proposal_mode = thermal_group_proposal_mode[collision_group]
                    if proposal_mode == 0:
                        target_vx, target_vy, target_vz = _maxwellian_target_velocity(
                            component_std, generator
                        )
                    else:
                        mixture_denominator = speed + mean_speed
                        if generator.random() < speed / mixture_denominator:
                            target_vx, target_vy, target_vz = (
                                _maxwellian_target_velocity(component_std, generator)
                            )
                        else:
                            target_vx, target_vy, target_vz = (
                                _speed_weighted_maxwellian_target_velocity(
                                    component_std,
                                    generator,
                                )
                            )
                    relative_vx = after_vx - target_vx
                    relative_vy = after_vy - target_vy
                    relative_vz = after_vz - target_vz
                    relative_speed = math.sqrt(
                        relative_vx * relative_vx
                        + relative_vy * relative_vy
                        + relative_vz * relative_vz
                    )
                    target_speed = math.sqrt(
                        target_vx * target_vx
                        + target_vy * target_vy
                        + target_vz * target_vz
                    )
                    collision_energy = _energy_eV(
                        relative_vx,
                        relative_vy,
                        relative_vz,
                    )
                    if proposal_mode == 0:
                        physical_rate_coefficient = relative_speed
                    else:
                        physical_rate_coefficient = (
                            relative_speed
                            * (speed + mean_speed)
                            / max(speed + target_speed, _FLOAT_TINY)
                        )

                evaluate_prepared_cross_sections(
                    collision_energy,
                    1.0,
                    collision_grid_eV,
                    collision_values_m2,
                    collision_slopes_m2_eV,
                    collision_min_eV,
                    collision_max_eV,
                    collision_right_m2,
                    collision_sigma,
                )
                collision_total = 0.0
                for process_index in range(collision_sigma.size):
                    if process_collision_group[process_index] == collision_group:
                        collision_total += collision_sigma[process_index]
                if collision_total <= 0.0:
                    continue
                physical_rate_coefficient *= collision_total
                collision_probability = physical_rate_coefficient / group_bound
                if collision_probability > 1.0 + 2.0e-12:
                    return (
                        positions,
                        velocities,
                        times,
                        weights,
                        lineages,
                        size,
                        event_observation_time_s,
                        residence_totals,
                        secondary_events,
                        KERNEL_MAJORANT_EXCEEDED,
                        density_m3 * physical_rate_coefficient,
                        -1,
                    )
                if generator.random() > min(max(collision_probability, 0.0), 1.0):
                    continue

                threshold = generator.random() * collision_total
                cumulative = 0.0
                choice = -1
                for process_index in range(collision_sigma.size):
                    if process_collision_group[process_index] != collision_group:
                        continue
                    cumulative += collision_sigma[process_index]
                    if threshold < cumulative:
                        choice = process_index
                        break
                if choice < 0:
                    continue
                kind = process_kind[choice]
                event_loss_eV = math.nan
                if kind == PROCESS_SCATTERING:
                    mu = generator.uniform(-1.0, 1.0)
                    vx, vy, vz = _elastic_binary_scatter(
                        after_vx,
                        after_vy,
                        after_vz,
                        target_vx,
                        target_vy,
                        target_vz,
                        mu,
                        process_target_mass_fraction[choice],
                        generator,
                    )
                    post_energy = _energy_eV(vx, vy, vz)
                    if (
                        zero_high_energy_policy
                        and post_energy > max_cross_section_energy_eV
                    ):
                        return (
                            positions,
                            velocities,
                            times,
                            weights,
                            lineages,
                            size,
                            event_observation_time_s,
                            residence_totals,
                            secondary_events,
                            KERNEL_CROSS_SECTION_RANGE_EXCEEDED,
                            post_energy,
                            -1,
                        )
                    if post_energy > max_energy_limit_eV:
                        return (
                            positions,
                            velocities,
                            times,
                            weights,
                            lineages,
                            size,
                            event_observation_time_s,
                            residence_totals,
                            secondary_events,
                            KERNEL_ENERGY_LIMIT_EXCEEDED,
                            post_energy,
                            -1,
                        )
                    event_loss_eV = final_energy - post_energy
                    velocities[particle, 0] = vx
                    velocities[particle, 1] = vy
                    velocities[particle, 2] = vz
                elif kind == PROCESS_IONIZATION:
                    reaction_threshold = max(process_threshold_eV[choice], 0.0)
                    available = max(final_energy - reaction_threshold, 0.0)
                    if ionization_model == IONIZATION_EQUAL:
                        primary_energy = 0.5 * available
                        secondary_energy = primary_energy
                    else:
                        secondary_energy = min(
                            max(secondary_electron_energy_eV, 0.0),
                            available,
                        )
                        primary_energy = max(available - secondary_energy, 0.0)
                    event_loss_eV = reaction_threshold
                    dx, dy, dz = _random_direction(generator)
                    primary_speed = math.sqrt(
                        max(_SPEED_FACTOR * max(primary_energy, 1.0e-4), 0.0)
                    )
                    velocities[particle, 0] = primary_speed * dx
                    velocities[particle, 1] = primary_speed * dy
                    velocities[particle, 2] = primary_speed * dz
                    if size >= weights.size:
                        (
                            positions,
                            velocities,
                            times,
                            weights,
                            lineages,
                        ) = _grow_particle_storage(
                            positions,
                            velocities,
                            times,
                            weights,
                            lineages,
                            size,
                        )
                    dx, dy, dz = _random_direction(generator)
                    secondary_speed = math.sqrt(
                        max(_SPEED_FACTOR * max(secondary_energy, 1.0e-4), 0.0)
                    )
                    positions[size, 0] = positions[particle, 0]
                    positions[size, 1] = positions[particle, 1]
                    positions[size, 2] = positions[particle, 2]
                    velocities[size, 0] = secondary_speed * dx
                    velocities[size, 1] = secondary_speed * dy
                    velocities[size, 2] = secondary_speed * dz
                    times[size] = times[particle]
                    weights[size] = weight
                    lineages[size] = lineages[particle]
                    size += 1
                    secondary_events += 1
                elif kind == PROCESS_EXCITATION or kind == PROCESS_SUPERELASTIC:
                    reaction_loss = process_threshold_eV[choice]
                    if kind == PROCESS_SUPERELASTIC:
                        reaction_loss = -abs(reaction_loss)
                    post_energy = max(final_energy - reaction_loss, 1.0e-4)
                    event_loss_eV = reaction_loss
                    dx, dy, dz = _random_direction(generator)
                    post_speed = math.sqrt(max(_SPEED_FACTOR * post_energy, 0.0))
                    velocities[particle, 0] = post_speed * dx
                    velocities[particle, 1] = post_speed * dy
                    velocities[particle, 2] = post_speed * dz

                if sampling_enabled and math.isfinite(event_loss_eV):
                    rate_event_counts[choice] += 1
                    rate_weighted_event_counts[choice] += weight
                    rate_event_loss_sample_counts[choice] += 1
                    rate_weighted_event_loss_eV[choice] += weight * event_loss_eV
                    rate_weighted_event_counts_by_block[
                        choice, reaction_block_index
                    ] += weight
                    rate_weighted_event_loss_eV_by_block[
                        choice, reaction_block_index
                    ] += weight * event_loss_eV

        times[particle] = barrier_time_s
        particle += 1

    return (
        positions,
        velocities,
        times,
        weights,
        lineages,
        size,
        event_observation_time_s,
        residence_totals,
        secondary_events,
        KERNEL_OK,
        0.0,
        -1,
    )
