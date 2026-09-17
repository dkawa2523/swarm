"""Observable evaluation and canonical result assembly for Propagator cases."""

from __future__ import annotations

from time import perf_counter
from typing import Mapping

import numpy as np

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.constants import TOWNSEND
from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.core.solver_configs import PropagatorInternalConfig
from electron_swarm.core.transport import ElectronTransport
from electron_swarm.solvers.propagator.memory import preflight_memory_bytes
from electron_swarm.solvers.propagator.models import PropagatorSteadySolution
from electron_swarm.solvers.propagator.observables import (
    cell_integrated_rates,
    eedf_from_population,
    elastic_energy_loss_rate_coefficient_eV_m3_s,
    energy_angle_distribution,
    flux_drift_velocity_m_s,
    mean_energy_eV,
)
from electron_swarm.solvers.propagator.operator import PropagatorOperator
from electron_swarm.solvers.propagator.shell_response import (
    SHELL_COEFFICIENT_ERROR_TOLERANCE,
    SHELL_MAX_COEFFICIENT_SEGMENTS,
    SHELL_RATIONAL_ERROR_TOLERANCE,
)


def build_case_result(
    *,
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    solver_config: PropagatorInternalConfig,
    qualification: Mapping[str, object],
    operator: PropagatorOperator,
    solution: PropagatorSteadySolution,
    gas_number_density_m3: float,
    e_over_n_Td: float,
    case_id: str,
    solver_name: str,
    timing: dict[str, float],
    total_started: float,
) -> SwarmCaseResult:
    """Project one converged population into the public result schema."""

    observable_started = perf_counter()
    population = solution.population
    grid = operator.grid
    eedf = eedf_from_population(grid, population)
    mean_energy = mean_energy_eV(grid, population)
    drift = flux_drift_velocity_m_s(grid, population)
    reduced_mobility = drift / (float(e_over_n_Td) * TOWNSEND)
    transport = ElectronTransport(
        definition="propagator_flux_velocity_moment",
        gas_number_density_m3=gas_number_density_m3,
        drift_velocity_m_s=drift,
        reduced_mobility_m2_V_s_m3=reduced_mobility,
        reduced_diffusion_L_m2_s_m3=None,
        reduced_diffusion_T_m2_s_m3=None,
    )
    rate_data = cell_integrated_rates(
        config,
        cross_sections,
        grid,
        population,
        case_id=case_id,
        e_over_n_Td=e_over_n_Td,
        solver_name=solver_name,
    )
    net_frequency = rate_data.net_ionization_frequency_s
    discrete_growth_frequency = operator.reaction_number_rate_s_inv(population)
    effective_townsend = net_frequency / max(
        abs(drift) * gas_number_density_m3,
        1.0e-300,
    )
    angle_distribution = energy_angle_distribution(grid, population)
    elastic_loss = elastic_energy_loss_rate_coefficient_eV_m3_s(
        grid,
        operator.collisions,
        population,
        gas_number_density_m3,
    )
    timing["observables"] += perf_counter() - observable_started
    timing["total"] = perf_counter() - total_started
    diagnostics = _case_diagnostics(
        config=config,
        solver_config=solver_config,
        qualification=qualification,
        operator=operator,
        solution=solution,
        eedf=eedf,
        elastic_loss=elastic_loss,
        net_frequency=net_frequency,
        discrete_growth_frequency=discrete_growth_frequency,
        timing=timing,
    )
    metadata = _case_metadata(config, operator, transport)
    return SwarmCaseResult(
        solver=solver_name,
        case_id=case_id,
        e_over_n_Td=e_over_n_Td,
        mean_energy_eV=mean_energy,
        net_ionization_frequency_s=net_frequency,
        effective_townsend_m2=effective_townsend,
        transport=transport,
        energy_eV=grid.energy_centers_eV.copy(),
        eedf=eedf,
        energy_widths_eV=grid.energy_widths_eV.copy(),
        rates=rate_data.rates,
        metadata=metadata,
        diagnostics={"propagator": diagnostics},
        energy_angle_distribution=angle_distribution,
    )


def _case_diagnostics(
    *,
    config: SwarmConfig,
    solver_config: PropagatorInternalConfig,
    qualification: Mapping[str, object],
    operator: PropagatorOperator,
    solution: PropagatorSteadySolution,
    eedf: np.ndarray,
    elastic_loss: float,
    net_frequency: float,
    discrete_growth_frequency: float,
    timing: dict[str, float],
) -> dict[str, object]:
    grid = operator.grid
    algorithm_diagnostics = dict(solution.diagnostics.timings_s)
    preflight_bytes = preflight_memory_bytes(
        energy_cells=grid.energy_cells,
        polar_cells=grid.polar_cells,
        angular_model=config.physics.angular_scattering.model,
    )
    response_memory_bytes = int(
        algorithm_diagnostics.get("response_operator_memory_bytes", 0.0)
    )
    measured_core_memory_bytes = operator.memory_bytes + response_memory_bytes
    estimated_peak_memory_bytes = max(preflight_bytes, measured_core_memory_bytes)
    memory_limit_bytes = int(solver_config.max_memory_mb) * 1024 * 1024
    if measured_core_memory_bytes > memory_limit_bytes:
        raise MemoryError(
            "propagator stationary response uses "
            f"{measured_core_memory_bytes} bytes, exceeding "
            f"max_memory_mb={solver_config.max_memory_mb}"
        )
    solution.diagnostics.timings_s = {**timing, **algorithm_diagnostics}
    diagnostics = solution.diagnostics.as_dict()
    diagnostics.update(
        {
            "schema": "swarm.propagator_diagnostics.v2",
            "energy_cells": grid.energy_cells,
            "polar_cells": grid.polar_cells,
            "energy_max_eV": float(grid.energy_edges_eV[-1]),
            "adaptive_cycles": 0,
            "final_grid_iterations": solution.diagnostics.iterations,
            "residual_L1": solution.diagnostics.operator_residual_L1,
            "residual_tolerance": max(
                solver_config.convergence_tolerance,
                1.0e-8,
            ),
            "tail_probability_target": solver_config.tail_probability_target,
            "edge_to_peak": float(eedf[-1])
            / max(float(np.max(eedf)), 1.0e-300),
            "edge_to_peak_target": max(
                solver_config.tail_probability_target,
                1.0e-10,
            ),
            "grid_max_eV": float(grid.energy_edges_eV[-1]),
            "grid_max_limit_eV": solver_config.max_eV_limit,
            "growth_frequency_s_inv": solution.growth_frequency_s_inv,
            "reported_net_ionization_frequency_s_inv": net_frequency,
            "discrete_population_growth_frequency_s_inv": (
                discrete_growth_frequency
            ),
            "growth_number_balance_relative_difference": abs(
                solution.growth_frequency_s_inv - discrete_growth_frequency
            )
            / max(
                abs(solution.growth_frequency_s_inv),
                abs(discrete_growth_frequency),
                1.0,
            ),
            "estimated_peak_memory_bytes": estimated_peak_memory_bytes,
            "measured_core_memory_bytes": measured_core_memory_bytes,
            **_algorithm_metadata(config),
            "elastic_energy_loss": {
                "schema": "swarm.elastic_energy_loss.v1",
                "status": "available",
                "symbol": "K_epsilon_el",
                "rate_coefficient_eV_m3_s": elastic_loss,
                "estimator": (
                    "same_discrete_finite_temperature_elastic_sg_operator"
                ),
                "operator": (
                    "finite_volume_scharfetter_gummel_elastic_energy_generator"
                ),
                "source_eedf": "same_solved_energy_angle_distribution",
                "neutral_thermal_motion_model": (
                    "finite_temperature_fokker_planck"
                ),
                "gas_temperature_terms_included": True,
                "gas_temperature_K": float(config.conditions.gas_temperature_K),
                "sign_convention": "positive_is_net_electron_energy_loss",
                "uncertainty_status": "deterministic_kinetic_solver",
            },
            "qualification": dict(qualification),
        }
    )
    return diagnostics


def _algorithm_metadata(config: SwarmConfig) -> dict[str, object]:
    return {
        "acceleration_scheme": "stationary_collision_coupled_shell_response",
        "steady_algorithm": (
            "positive_perron_cone_certified_safeguarded_"
            "secant_brent_number_balance"
        ),
        "perron_branch_selection": (
            "staged_largest_real_cone_certificate_with_multiritz_fallback"
        ),
        "angular_kernel_source": (
            "isotropic_reciprocal_cell_kernel"
            if config.physics.angular_scattering.model == "isotropic"
            else "maxent_p1_reversible_cell_joint_exact_p1"
        ),
        "elastic_recoil_model": (
            "momentum_transfer_driven_finite_temperature_"
            "first_mass_ratio_reversible_sg"
        ),
        "elastic_frequency_model": "separate_cell_integrated_total_and_momentum",
        "elastic_equilibrium_measure": "cell_integrated_maxwell_energy_mass",
        "neutral_thermal_motion_model": "finite_temperature_fokker_planck",
        "origin_treatment": (
            "characteristic_aligned_cut_cell_common_refinement_"
            "with_radial_rate_quadrature"
        ),
        "shell_coefficient_model": (
            "exact_inverse_speed_average_with_bounded_refinement"
        ),
        "shell_response_error_control": {
            "rational_richardson_limit": SHELL_RATIONAL_ERROR_TOLERANCE,
            "coefficient_richardson_limit": SHELL_COEFFICIENT_ERROR_TOLERANCE,
            "maximum_coefficient_segments": SHELL_MAX_COEFFICIENT_SEGMENTS,
        },
        "outer_energy_boundary": (
            "zero_incoming_with_measured_characteristic_escape"
        ),
        "energy_domain_strategy": (
            "one_shot_core_plus_stretched_tail_to_configured_ceiling"
        ),
        "energy_core_coordinate": "sinh_stretched_speed",
        "energy_core_coordinate_strength": 2.5,
        "response_shift_policy": "minimum_nonnegative_local_outflow_shift",
        "reaction_rate_estimator": (
            "xs_knot_threshold_partitioned_cell_quadrature"
        ),
        "inelastic_energy_transfer": (
            "positive_weak_projection_with_incident_and_daughter_energy_identity"
        ),
        "transport_moment_quadrature": (
            "cell_integrated_piecewise_constant_energy_density"
        ),
    }


def _case_metadata(
    config: SwarmConfig,
    operator: PropagatorOperator,
    transport: ElectronTransport,
) -> dict[str, str]:
    return {
        "velocity_space_representation": "axisymmetric_energy_theta_cells",
        "transport_definition": transport.definition,
        "transport_components": "drift_velocity;mobility",
        "inelastic_angular_model": "isotropic_integral_xs_closure",
        "inelastic_radial_transfer": (
            "xs_knot_threshold_partitioned_positive_weak_projection"
        ),
        "elastic_recoil_model": (
            "momentum_transfer_driven_finite_temperature_"
            "first_mass_ratio_reversible_sg"
        ),
        "neutral_thermal_motion_model": "finite_temperature_fokker_planck",
        "elastic_collision_xs_role": ";".join(
            operator.collisions.elastic_xs_roles
        ),
        "elastic_total_xs_process_ids": ";".join(
            operator.collisions.elastic_total_process_ids
        ),
        "elastic_momentum_xs_process_ids": ";".join(
            operator.collisions.elastic_momentum_process_ids
        ),
    }
