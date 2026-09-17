from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse
from scipy.special import gammainc

from electron_swarm import load_config
from electron_swarm.core.constants import (
    AMU_KG,
    BOLTZMANN_J_K,
    E_CHARGE_C,
    ELECTRON_MASS_KG,
    EV_TO_J,
    TOWNSEND,
)
from electron_swarm.core.cross_sections import (
    CrossSectionProcess,
    CrossSectionSet,
    ProcessType,
    ScatteringRole,
    load_cross_sections,
)
from electron_swarm.core.solver_configs import build_internal_solver_configs
from electron_swarm.physics.kinetics import gas_number_density
from electron_swarm.physics.electron_neutral import speed_from_energy_m_s
from electron_swarm.solvers.propagator.angular_kernel import (
    isotropic_destination_weights,
    maxent_p1_transition_matrices,
)
from electron_swarm.solvers.propagator.collisions import (
    build_collision_operator,
    collision_generator_action,
    reaction_number_rate_s_inv,
)
from electron_swarm.solvers.propagator.grid import build_propagator_grid
from electron_swarm.solvers.propagator.observables import (
    elastic_energy_loss_rate_coefficient_eV_m3_s,
    flux_drift_velocity_m_s,
)
from electron_swarm.solvers.propagator.operator import PropagatorOperator
from electron_swarm.solvers.propagator.origin_response import (
    _radial_band_average_rates,
    build_origin_operator,
    build_origin_response,
)
from electron_swarm.solvers.propagator.shell_response import (
    AffineShellResponse,
    SHELL_COEFFICIENT_ERROR_TOLERANCE,
    SHELL_RATIONAL_ERROR_TOLERANCE,
    build_shell_response,
    compose_affine_responses,
)
from electron_swarm.solvers.propagator.steady import (
    SteadyControls,
    maxwellian_initial_population,
    solve_steady_state,
)
from tests.support.propagator_reference import (
    affine_origin_response_relative_distance,
    affine_response_relative_distance,
    continuous_shell_reference,
)

from product_helpers import base_product_config, write_config


def _build_operator(
    tmp_path: Path,
    *,
    energy_cells: int = 64,
    polar_cells: int = 8,
) -> tuple[object, PropagatorOperator]:
    data = base_product_config(tmp_path, ["propagator"])
    data["physics"]["energy_grid_policy"]["threshold_refinement"] = False
    data["solvers"]["propagator"] = {
        "energy_cells": max(64, energy_cells),
        "polar_cells": polar_cells,
        "max_iterations": 1000,
        "max_memory_mb": 128,
    }
    config = load_config(write_config(tmp_path, data))
    cross_sections = load_cross_sections(config.cross_sections, config.conditions)
    internal = replace(
        build_internal_solver_configs(
            config.solvers,
            config.physics,
        ).propagator,
        energy_cells=energy_cells,
        polar_cells=polar_cells,
    )
    grid = build_propagator_grid(internal, cross_sections)
    collisions = build_collision_operator(config, cross_sections, grid)
    return config, PropagatorOperator.build(grid, collisions)


def _acceleration(config: object, e_over_n_Td: float) -> float:
    density = gas_number_density(config)
    return (
        E_CHARGE_C
        * float(e_over_n_Td)
        * TOWNSEND
        * density
        / ELECTRON_MASS_KG
    )


def _build_maxent_validation_operator(
    tmp_path: Path,
    *,
    energy_cells: int = 64,
    polar_cells: int = 8,
) -> tuple[object, PropagatorOperator]:
    data = base_product_config(tmp_path, ["propagator"])
    data["cross_sections"] = {
        "format": "csv",
        "high_energy_extrapolation": "error",
        "files": [
            {
                "path": (
                    Path(__file__).resolve().parents[1]
                    / "examples"
                    / "cross_sections"
                    / "argon_magboltz_11_17_gas2_elastic.csv"
                ).as_posix(),
                "species": "Ar",
                "format": "csv",
            }
        ],
    }
    data["physics"]["angular_scattering"] = {
        "model": "maxent_p1",
        "higher_moment_closure": "maxent",
    }
    data["physics"]["energy_grid_policy"]["threshold_refinement"] = False
    data["solvers"]["propagator"] = {
        "energy_cells": energy_cells,
        "polar_cells": polar_cells,
        "max_iterations": 1000,
        "max_memory_mb": 128,
    }
    config = load_config(write_config(tmp_path, data, "maxent_operator.yaml"))
    cross_sections = load_cross_sections(
        config.cross_sections,
        config.conditions,
    )
    internal = build_internal_solver_configs(
        config.solvers,
        config.physics,
    ).propagator
    grid = build_propagator_grid(internal, cross_sections)
    collisions = build_collision_operator(config, cross_sections, grid)
    return config, PropagatorOperator.build(grid, collisions)


def test_one_shot_grid_keeps_the_resolved_core_and_adds_a_stretched_tail(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["propagator"])
    data["physics"]["energy_grid_policy"]["threshold_refinement"] = True
    data["solvers"]["propagator"] = {
        "energy_cells": 64,
        "polar_cells": 8,
        "max_memory_mb": 128,
    }
    config = load_config(write_config(tmp_path, data))
    cross_sections = load_cross_sections(config.cross_sections, config.conditions)
    internal = build_internal_solver_configs(config.solvers, config.physics).propagator
    core = build_propagator_grid(internal, cross_sections, max_eV=100.0)
    extended = build_propagator_grid(internal, cross_sections, max_eV=160.0)

    assert np.array_equal(
        core.energy_edges_eV,
        extended.energy_edges_eV[: len(core.energy_edges_eV)],
    )
    assert extended.energy_edges_eV[-1] == pytest.approx(160.0)
    assert extended.energy_cells > core.energy_cells
    assert 11.55 in core.energy_edges_eV
    assert 15.76 in core.energy_edges_eV


def test_sinh_speed_core_is_nested_and_resolves_sub_electronvolt_energy(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["propagator"])
    data["physics"]["energy_grid_policy"]["threshold_refinement"] = False
    data["solvers"]["propagator"] = {
        "energy_cells": 64,
        "polar_cells": 8,
        "max_memory_mb": 128,
    }
    config = load_config(write_config(tmp_path, data))
    cross_sections = load_cross_sections(
        config.cross_sections,
        config.conditions,
    )
    internal = build_internal_solver_configs(
        config.solvers,
        config.physics,
    ).propagator
    coarse = build_propagator_grid(internal, cross_sections)
    fine = build_propagator_grid(
        replace(internal, energy_cells=128),
        cross_sections,
    )

    assert fine.energy_edges_eV[::2] == pytest.approx(
        coarse.energy_edges_eV,
        abs=2.0e-14,
    )
    assert int(np.count_nonzero(coarse.energy_edges_eV <= 1.0)) >= 14
    speed_widths = np.diff(coarse.speed_edges_m_s)
    assert speed_widths[-1] > 4.0 * speed_widths[0]


def test_collision_operator_is_conservative_positive_and_thermalizing(
    tmp_path: Path,
) -> None:
    config, operator = _build_operator(tmp_path)
    grid = operator.grid
    collisions = operator.collisions
    energy_generator = collisions.elastic_energy_generator_s_inv
    scale = max(float(np.max(np.abs(energy_generator.data))), 1.0)
    column_error = float(
        np.max(np.abs(np.asarray(energy_generator.sum(axis=0)).reshape(-1)))
    )
    off_diagonal = energy_generator - sparse.diags(energy_generator.diagonal())

    assert column_error <= 2.0e-12 * scale
    assert float(np.min(off_diagonal.data, initial=0.0)) >= -2.0e-14 * scale
    equilibrium = collisions.elastic_equilibrium_mass
    assert np.sum(equilibrium) == pytest.approx(1.0, abs=2.0e-14)
    assert np.linalg.norm(energy_generator @ equilibrium, ord=1) <= (
        2.0e-12 * scale
    )

    population = np.zeros((grid.energy_cells, grid.polar_cells), dtype=float)
    energy_index = int(np.searchsorted(grid.energy_centers_eV, 40.0))
    population[energy_index, 3] = 1.0
    action = collision_generator_action(collisions, population)
    assert float(np.sum(action)) == pytest.approx(
        reaction_number_rate_s_inv(collisions, population),
        rel=2.0e-13,
        abs=2.0e-6,
    )
    for transfer in collisions.inelastic:
        gain_columns = np.asarray(
            transfer.gain_matrix_s_inv.sum(axis=0)
        ).reshape(-1)
        assert gain_columns == pytest.approx(
            transfer.daughter_count * transfer.frequency_s_inv,
            rel=2.0e-13,
            abs=2.0e-7,
        )

    isotropic_population = (
        equilibrium[:, None] * collisions.isotropic_weights[None, :]
    )
    density = gas_number_density(config)
    reported_equilibrium_loss = (
        elastic_energy_loss_rate_coefficient_eV_m3_s(
            grid,
            collisions,
            isotropic_population,
            density,
        )
        * density
    )
    assert reported_equilibrium_loss == pytest.approx(0.0, abs=2.0e-10)

    shell_mass = np.linspace(1.0, 2.0, grid.energy_cells)
    shell_mass /= float(np.sum(shell_mass))
    test_population = (
        shell_mass[:, None] * collisions.isotropic_weights[None, :]
    )
    elastic_action = energy_generator @ test_population
    expected_loss = -float(
        np.sum(elastic_action * grid.energy_centers_eV[:, None])
    )
    reported_loss = (
        elastic_energy_loss_rate_coefficient_eV_m3_s(
            grid,
            collisions,
            test_population,
            density,
        )
        * density
    )
    assert reported_loss == pytest.approx(expected_loss, rel=2.0e-13, abs=2.0e-9)

    cell_probe = np.zeros_like(test_population)
    cell_probe[0, 0] = 1.0
    lower, upper = grid.energy_edges_eV[:2]
    exact_cell_speed = (
        float(speed_from_energy_m_s(np.asarray([1.0]))[0])
        * (2.0 / 3.0)
        * (upper**1.5 - lower**1.5)
        / (upper - lower)
    )
    assert flux_drift_velocity_m_s(grid, cell_probe) == pytest.approx(
        exact_cell_speed * grid.mu_centers[0],
        rel=2.0e-15,
    )


def test_inelastic_cell_quadrature_splits_thresholds_and_preserves_kinematics(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["propagator"])
    data["physics"]["ionization"] = {
        "energy_sharing": "primary_secondary",
        "secondary_electron_energy_eV": 0.37,
    }
    data["physics"]["energy_grid_policy"]["threshold_refinement"] = False
    data["solvers"]["propagator"] = {
        "energy_cells": 64,
        "polar_cells": 8,
        "max_memory_mb": 128,
    }
    config = load_config(write_config(tmp_path, data, "weak_transfer.yaml"))
    energy = np.asarray((0.0, 1000.0))
    sigma = 2.0e-20
    thresholds = {
        ProcessType.ATTACHMENT: 3.21,
        ProcessType.EXCITATION: 12.345,
        ProcessType.IONIZATION: 17.891,
    }
    processes = [
        CrossSectionProcess(
            species="Ar",
            process="constant elastic",
            process_type=ProcessType.ELASTIC,
            energy_eV=energy,
            cross_section_m2=np.full(2, 3.0e-20),
            mass_amu=39.948,
        )
    ]
    for process_type, threshold in thresholds.items():
        processes.append(
            CrossSectionProcess(
                species="Ar",
                process=f"constant {process_type.value}",
                process_type=process_type,
                energy_eV=energy,
                cross_section_m2=np.full(2, sigma),
                threshold_eV=threshold,
                mass_amu=39.948,
            )
        )
    cross_sections = CrossSectionSet(processes)
    internal = build_internal_solver_configs(
        config.solvers,
        config.physics,
    ).propagator
    grid = build_propagator_grid(internal, cross_sections)
    collisions = build_collision_operator(config, cross_sections, grid)
    density = gas_number_density(config)
    speed_scale = float(speed_from_energy_m_s(np.asarray([1.0]))[0])

    for transfer in collisions.inelastic:
        process_type = ProcessType(transfer.process_type)
        threshold = thresholds[process_type]
        source = int(np.searchsorted(grid.energy_edges_eV, threshold) - 1)
        lower = float(grid.energy_edges_eV[source])
        upper = float(grid.energy_edges_eV[source + 1])
        assert lower < threshold < upper
        expected_frequency = (
            density
            * sigma
            * speed_scale
            * (2.0 / 3.0)
            * (upper**1.5 - threshold**1.5)
            / (upper - lower)
        )
        expected_incident_energy_rate = (
            density
            * sigma
            * speed_scale
            * (2.0 / 5.0)
            * (upper**2.5 - threshold**2.5)
            / (upper - lower)
        )
        assert transfer.frequency_s_inv[source] == pytest.approx(
            expected_frequency,
            rel=3.0e-14,
        )
        assert transfer.incident_energy_rate_eV_s_inv[source] == pytest.approx(
            expected_incident_energy_rate,
            rel=3.0e-14,
        )
        if transfer.daughter_count == 0:
            assert transfer.gain_matrix_s_inv.nnz == 0
            assert transfer.daughter_energy_rate_eV_s_inv[source] == 0.0
            continue
        expected_daughter_energy_rate = (
            expected_incident_energy_rate
            - threshold * expected_frequency
        )
        assert transfer.daughter_energy_rate_eV_s_inv[source] == pytest.approx(
            expected_daughter_energy_rate,
            rel=3.0e-13,
            abs=2.0e-8,
        )
        gain_columns = np.asarray(
            transfer.gain_matrix_s_inv.sum(axis=0)
        ).reshape(-1)
        assert gain_columns[source] == pytest.approx(
            transfer.daughter_count * expected_frequency,
            rel=3.0e-14,
        )


def test_inelastic_discrete_energy_balance_converges_quadratically(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["propagator"])
    data["physics"]["energy_grid_policy"]["threshold_refinement"] = True
    data["solvers"]["propagator"] = {
        "energy_cells": 64,
        "polar_cells": 8,
        "max_memory_mb": 128,
    }
    config = load_config(write_config(tmp_path, data, "inelastic_refinement.yaml"))
    cross_sections = load_cross_sections(
        config.cross_sections,
        config.conditions,
    )
    base = build_internal_solver_configs(
        config.solvers,
        config.physics,
    ).propagator
    errors: dict[str, list[float]] = {}

    for energy_cells in (64, 128, 256):
        grid = build_propagator_grid(
            replace(base, energy_cells=energy_cells),
            cross_sections,
        )
        collisions = build_collision_operator(config, cross_sections, grid)
        for transfer in collisions.inelastic:
            if transfer.daughter_count == 0:
                continue
            deposited_energy_rate = np.asarray(
                grid.energy_centers_eV @ transfer.gain_matrix_s_inv
            ).reshape(-1)
            discrete_change = (
                deposited_energy_rate
                - grid.energy_centers_eV * transfer.frequency_s_inv
            )
            exact_threshold_change = (
                -transfer.threshold_eV * transfer.frequency_s_inv
            )
            scale = float(np.sum(np.abs(exact_threshold_change)))
            error = float(
                np.sum(np.abs(discrete_change - exact_threshold_change))
                / max(scale, 1.0e-300)
            )
            errors.setdefault(transfer.process_type, []).append(error)

    assert set(errors) == {"excitation", "ionization"}
    for process_errors in errors.values():
        assert process_errors[0] < 3.0e-3
        assert process_errors[1] < 0.27 * process_errors[0]
        assert process_errors[2] < 0.27 * process_errors[1]


def test_angular_kernel_is_positive_reciprocal_and_exact_p1(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["propagator"])
    data["solvers"]["propagator"] = {
        "energy_cells": 64,
        "polar_cells": 8,
        "max_memory_mb": 128,
    }
    config = load_config(write_config(tmp_path, data))
    cross_sections = load_cross_sections(config.cross_sections, config.conditions)
    internal = build_internal_solver_configs(config.solvers, config.physics).propagator
    requested_moments = np.asarray((-0.95, 0.0, 0.5, 0.8, 0.95, 0.999))

    for polar_cells in (36, 72):
        grid = build_propagator_grid(
            replace(internal, polar_cells=polar_cells),
            cross_sections,
        )
        weights = isotropic_destination_weights(grid)
        kernels = maxent_p1_transition_matrices(grid, requested_moments)
        for moment, kernel in zip(requested_moments, kernels, strict=True):
            joint = kernel * weights[None, :]
            assert float(np.min(kernel)) >= 0.0
            assert float(
                np.max(np.abs(np.sum(kernel, axis=0) - 1.0))
            ) <= 5.0e-10
            assert float(np.sum(np.abs(kernel @ weights - weights))) <= 5.0e-10
            assert float(np.max(np.abs(joint - joint.T))) <= 5.0e-10
            assert float(
                np.max(
                    np.abs(
                        grid.mu_centers @ kernel
                        - moment * grid.mu_centers
                    )
                )
            ) <= 5.0e-9
        # Closely spaced radial-band moments exercise continuation without
        # leaving the centrosymmetric dual subspace.
        band_moments = np.concatenate(
            (
                np.asarray((4.0e-5, 1.8e-4, 4.6e-4)),
                np.linspace(0.01, 0.8, 13),
            )
        )
        band_kernels = maxent_p1_transition_matrices(grid, band_moments)
        assert float(
            np.max(
                np.abs(
                    np.einsum(
                        "i,eij->ej",
                        grid.mu_centers,
                        band_kernels,
                        optimize=True,
                    )
                    - band_moments[:, None] * grid.mu_centers[None, :]
                )
            )
        ) <= 5.0e-9


def test_elastic_total_and_momentum_integrals_have_distinct_core_roles(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["propagator"])
    data["physics"]["angular_scattering"] = {
        "model": "maxent_p1",
        "higher_moment_closure": "maxent",
    }
    data["solvers"]["propagator"] = {
        "energy_cells": 64,
        "polar_cells": 8,
        "max_memory_mb": 128,
    }
    config = load_config(write_config(tmp_path, data, "anisotropic.yaml"))
    energy = np.asarray((0.0, 1000.0))
    total = CrossSectionProcess(
        species="Ar",
        process="synthetic total",
        process_type=ProcessType.ELASTIC,
        scattering_role=ScatteringRole.ELASTIC_TOTAL,
        energy_eV=energy,
        cross_section_m2=np.full(2, 5.0e-20),
        mass_amu=39.948,
    )
    momentum = CrossSectionProcess(
        species="Ar",
        process="synthetic momentum",
        process_type=ProcessType.MOMENTUM,
        scattering_role=ScatteringRole.ELASTIC_MOMENTUM_TRANSFER,
        energy_eV=energy,
        cross_section_m2=np.full(2, 1.0e-20),
        mass_amu=39.948,
    )
    cross_sections = CrossSectionSet([total, momentum])
    internal = build_internal_solver_configs(
        config.solvers, config.physics
    ).propagator
    grid = build_propagator_grid(internal, cross_sections)
    collisions = build_collision_operator(config, cross_sections, grid)
    [transfer] = collisions.elastic

    assert transfer.collision_frequency_s_inv == pytest.approx(
        5.0 * transfer.momentum_frequency_s_inv,
        rel=3.0e-14,
    )
    assert transfer.mean_cosine == pytest.approx(
        np.full(grid.energy_cells, 0.8),
        abs=3.0e-14,
    )
    kT_eV = config.conditions.gas_temperature_K * BOLTZMANN_J_K / EV_TO_J
    mass_ratio = ELECTRON_MASS_KG / (39.948 * AMU_KG)
    assert collisions.elastic_A_eV_s == pytest.approx(
        2.0
        * mass_ratio
        * transfer.momentum_frequency_s_inv
        * (0.5 * kT_eV - grid.energy_centers_eV),
        rel=3.0e-14,
        abs=1.0e-12,
    )
    assert collisions.elastic_D_eV2_s == pytest.approx(
        2.0
        * mass_ratio
        * transfer.momentum_frequency_s_inv
        * grid.energy_centers_eV
        * kT_eV,
        rel=3.0e-14,
        abs=1.0e-12,
    )

    expected_equilibrium = np.diff(
        gammainc(1.5, grid.energy_edges_eV / kT_eV)
    )
    expected_equilibrium /= float(np.sum(expected_equilibrium))
    assert collisions.elastic_equilibrium_mass == pytest.approx(
        expected_equilibrium,
        abs=3.0e-14,
    )
    generator_scale = max(
        float(
            np.max(
                np.abs(collisions.elastic_energy_generator_s_inv.data),
                initial=0.0,
            )
        ),
        1.0,
    )
    assert float(
        np.linalg.norm(
            collisions.elastic_energy_generator_s_inv
            @ collisions.elastic_equilibrium_mass,
            ord=1,
        )
    ) <= 2.0e-12 * generator_scale

    origin_total, origin_momentum = _radial_band_average_rates(
        transfer,
        float(grid.energy_edges_eV[1]),
        16,
    )
    assert np.all(np.diff(origin_total) > 0.0)
    assert origin_total == pytest.approx(5.0 * origin_momentum, rel=3.0e-14)


def test_isotropic_closure_rejects_distinct_total_and_momentum_integrals(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["propagator"])
    data["solvers"]["propagator"] = {
        "energy_cells": 64,
        "polar_cells": 8,
        "max_memory_mb": 128,
    }
    config = load_config(write_config(tmp_path, data, "isotropic.yaml"))
    energy = np.asarray((0.0, 1000.0))
    cross_sections = CrossSectionSet(
        [
            CrossSectionProcess(
                species="Ar",
                process="synthetic total",
                process_type=ProcessType.ELASTIC,
                scattering_role=ScatteringRole.ELASTIC_TOTAL,
                energy_eV=energy,
                cross_section_m2=np.full(2, 2.0e-20),
                mass_amu=39.948,
            ),
            CrossSectionProcess(
                species="Ar",
                process="synthetic momentum",
                process_type=ProcessType.MOMENTUM,
                scattering_role=ScatteringRole.ELASTIC_MOMENTUM_TRANSFER,
                energy_eV=energy,
                cross_section_m2=np.full(2, 1.0e-20),
                mass_amu=39.948,
            ),
        ]
    )
    internal = build_internal_solver_configs(
        config.solvers, config.physics
    ).propagator
    grid = build_propagator_grid(internal, cross_sections)
    with pytest.raises(ValueError, match="use maxent_p1"):
        build_collision_operator(config, cross_sections, grid)


def test_shell_and_characteristic_origin_responses_are_positive_and_balanced(
    tmp_path: Path,
) -> None:
    config, operator = _build_operator(tmp_path)
    grid = operator.grid
    acceleration = _acceleration(config, 10.0)
    rng = np.random.default_rng(23)

    energy_index = 10
    shell_loss = float(
        operator.collisions.nonlocal_outflow_s_inv[energy_index] + 1.0
    )
    shell = build_shell_response(
        grid,
        operator.collisions,
        energy_index,
        acceleration,
        shell_loss,
    )
    shell_incoming = rng.random(grid.polar_cells)
    shell_source = rng.random(grid.polar_cells)
    shell_outgoing = (
        shell.scattering @ shell_incoming
        + shell.source_to_outflow @ shell_source
    )
    shell_population = (
        shell.incoming_to_population @ shell_incoming
        + shell.source_to_population @ shell_source
    )
    assert shell.minimum_entry >= -1.0e-12
    assert shell.rational_error_estimate <= 2.0e-5
    assert shell.coefficient_error_estimate <= 5.0e-4
    assert shell.coefficient_segments in {2, 4, 8, 16, 32, 64, 128}
    assert float(np.min(shell_outgoing)) >= 0.0
    assert float(np.min(shell_population)) >= 0.0
    assert (
        float(np.sum(shell_outgoing) - np.sum(shell_incoming))
        + shell_loss * float(np.sum(shell_population))
        - float(np.sum(shell_source))
    ) == pytest.approx(0.0, abs=2.0e-11)

    origin_loss = float(operator.collisions.nonlocal_outflow_s_inv[0] + 1.0)
    origin = build_origin_response(
        operator.origin,
        float(grid.speed_edges_m_s[1]),
        acceleration,
        origin_loss,
    )
    origin_incoming = rng.random(grid.polar_cells // 2)
    origin_source = rng.random(grid.polar_cells)
    origin_outgoing = (
        origin.incoming_to_outflow @ origin_incoming
        + origin.source_to_outflow @ origin_source
    )
    origin_population = (
        origin.incoming_to_population @ origin_incoming
        + origin.source_to_population @ origin_source
    )
    generator_scale = max(
        float(np.max(np.abs(operator.origin.collision_generator_s_inv.data))),
        1.0,
    )
    assert operator.origin.column_balance_error <= 2.0e-11 * generator_scale
    assert operator.origin.detailed_balance_error <= 2.0e-11 * generator_scale
    assert origin.minimum_entry >= -1.0e-12
    assert float(np.min(origin_outgoing)) >= 0.0
    assert float(np.min(origin_population)) >= 0.0
    assert (
        float(np.sum(origin_outgoing) - np.sum(origin_incoming))
        + origin_loss * float(np.sum(origin_population))
        - float(np.sum(origin_source))
    ) == pytest.approx(0.0, abs=2.0e-11)


def test_fused_redheffer_composition_matches_literal_block_equations() -> None:
    rng = np.random.default_rng(914)
    cells = 8
    half = cells // 2
    order = np.arange(cells, dtype=int)

    def response() -> AffineShellResponse:
        return AffineShellResponse(
            scattering=0.025 * rng.random((cells, cells)),
            source_to_outflow=rng.random((cells, cells)),
            incoming_to_population=rng.random((cells, cells)),
            source_to_population=rng.random((cells, cells)),
            boundary_order=order,
            conceptual_sublayers=3,
            minimum_entry=0.0,
            maximum_column_balance_error=0.0,
            rational_error_estimate=2.0e-7,
            coefficient_error_estimate=3.0e-6,
            coefficient_segments=4,
        )

    left = response()
    right = response()
    actual = compose_affine_responses(left, right)

    a1 = left.scattering[:half, :half]
    b1 = left.scattering[:half, half:]
    c1 = left.scattering[half:, :half]
    d1 = left.scattering[half:, half:]
    a2 = right.scattering[:half, :half]
    b2 = right.scattering[:half, half:]
    c2 = right.scattering[half:, :half]
    d2 = right.scattering[half:, half:]
    identity = np.eye(half)
    left_selector = np.zeros((half, cells), dtype=float)
    left_selector[:, :half] = identity
    right_selector = np.zeros((half, cells), dtype=float)
    right_selector[:, half:] = identity
    interface = identity - b1 @ c2
    positive = np.linalg.solve(
        interface,
        a1 @ left_selector + b1 @ d2 @ right_selector,
    )
    negative = c2 @ positive + d2 @ right_selector
    source_positive = np.linalg.solve(
        interface,
        left.source_to_outflow[:half]
        + b1 @ right.source_to_outflow[half:],
    )
    source_negative = (
        c2 @ source_positive + right.source_to_outflow[half:]
    )
    expected_scattering = np.vstack(
        (
            a2 @ positive + b2 @ right_selector,
            c1 @ left_selector + d1 @ negative,
        )
    )
    expected_outflow = np.vstack(
        (
            a2 @ source_positive + right.source_to_outflow[:half],
            d1 @ source_negative + left.source_to_outflow[half:],
        )
    )
    expected_incoming_population = (
        left.incoming_to_population @ np.vstack((left_selector, negative))
        + right.incoming_to_population @ np.vstack((positive, right_selector))
    )
    expected_source_population = (
        left.incoming_to_population
        @ np.vstack((np.zeros((half, cells)), source_negative))
        + left.source_to_population
        + right.incoming_to_population
        @ np.vstack((source_positive, np.zeros((half, cells))))
        + right.source_to_population
    )

    np.testing.assert_allclose(
        actual.scattering, expected_scattering, rtol=2.0e-14, atol=2.0e-14
    )
    np.testing.assert_allclose(
        actual.source_to_outflow, expected_outflow, rtol=2.0e-14, atol=2.0e-14
    )
    np.testing.assert_allclose(
        actual.incoming_to_population,
        expected_incoming_population,
        rtol=2.0e-14,
        atol=2.0e-14,
    )
    np.testing.assert_allclose(
        actual.source_to_population,
        expected_source_population,
        rtol=2.0e-14,
        atol=2.0e-14,
    )
    assert actual.conceptual_sublayers == 6
    assert actual.coefficient_segments == 8

def test_shell_error_control_bounds_independent_continuous_references(
    tmp_path: Path,
) -> None:
    config, operator = _build_maxent_validation_operator(tmp_path)
    samples = (
        (1.0, 1),
        (1.0, 3),
        (1.0, 5),
        (1.0, 10),
        (10.0, 20),
        (100.0, 40),
        (100.0, 62),
    )
    for field, energy_index in samples:
        acceleration = _acceleration(config, field)
        loss = float(
            operator.collisions.nonlocal_outflow_s_inv[energy_index] + 1.0
        )
        response = build_shell_response(
            operator.grid,
            operator.collisions,
            energy_index,
            acceleration,
            loss,
        )
        reference = continuous_shell_reference(
            operator.grid,
            operator.collisions,
            energy_index,
            acceleration,
            loss,
            relative_tolerance=2.0e-12,
            absolute_tolerance=1.0e-15,
        )
        actual_error = affine_response_relative_distance(
            response,
            reference.response,
        )
        estimated_error = (
            response.coefficient_error_estimate
            + response.rational_error_estimate
        )
        assert response.coefficient_segments >= 4
        assert reference.boundary_block_condition_number <= 1.0e6
        assert actual_error <= estimated_error
        assert actual_error <= (
            SHELL_COEFFICIENT_ERROR_TOLERANCE
            + SHELL_RATIONAL_ERROR_TOLERANCE
        )


def test_characteristic_origin_is_stable_under_refinement_for_maxent_argon(
    tmp_path: Path,
) -> None:
    config, operator = _build_maxent_validation_operator(
        tmp_path,
        energy_cells=300,
        polar_cells=36,
    )
    refined = build_origin_operator(
        operator.grid,
        operator.collisions,
        u_cells=288,
        radial_bands=144,
    )
    default_scale = max(
        float(
            np.max(
                np.abs(operator.origin.collision_generator_s_inv.data),
                initial=0.0,
            )
        ),
        1.0,
    )
    refined_scale = max(
        float(
            np.max(
                np.abs(refined.collision_generator_s_inv.data),
                initial=0.0,
            )
        ),
        1.0,
    )
    assert operator.origin.column_balance_error <= 2.0e-11 * default_scale
    assert operator.origin.detailed_balance_error <= 2.0e-11 * default_scale
    assert refined.column_balance_error <= 2.0e-11 * refined_scale
    assert refined.detailed_balance_error <= 2.0e-11 * refined_scale

    radius = float(operator.grid.speed_edges_m_s[1])
    loss = float(operator.collisions.nonlocal_outflow_s_inv[0] + 1.0)
    for field in (0.05, 0.1, 1.0):
        acceleration = _acceleration(config, field)
        coarse_response = build_origin_response(
            operator.origin,
            radius,
            acceleration,
            loss,
        )
        refined_response = build_origin_response(
            refined,
            radius,
            acceleration,
            loss,
        )
        assert coarse_response.minimum_entry >= 0.0
        assert refined_response.minimum_entry >= 0.0
        assert affine_origin_response_relative_distance(
            coarse_response,
            refined_response,
        ) <= 3.0e-3


def test_characteristic_origin_obeys_collisionless_and_rate_scaling_limits(
    tmp_path: Path,
) -> None:
    config, operator = _build_maxent_validation_operator(tmp_path)
    radius = float(operator.grid.speed_edges_m_s[1])
    acceleration = _acceleration(config, 1.0)
    zero_collision = replace(
        operator.origin,
        collision_generator_s_inv=sparse.csr_matrix(
            operator.origin.collision_generator_s_inv.shape,
            dtype=float,
        ),
        detailed_balance_error=0.0,
        column_balance_error=0.0,
    )
    ballistic = build_origin_response(
        zero_collision,
        radius,
        acceleration,
        1.0e-12 * acceleration / radius,
    )
    assert ballistic.incoming_to_outflow == pytest.approx(
        np.eye(operator.grid.polar_cells // 2),
        abs=3.0e-12,
    )

    loss = float(operator.collisions.nonlocal_outflow_s_inv[0] + 1.0)
    baseline = build_origin_response(
        operator.origin,
        radius,
        acceleration,
        loss,
    )
    factor = 7.0
    scaled_operator = replace(
        operator.origin,
        collision_generator_s_inv=(
            factor * operator.origin.collision_generator_s_inv
        ),
        detailed_balance_error=factor * operator.origin.detailed_balance_error,
        column_balance_error=factor * operator.origin.column_balance_error,
    )
    scaled = build_origin_response(
        scaled_operator,
        radius,
        factor * acceleration,
        factor * loss,
    )
    assert scaled.incoming_to_outflow == pytest.approx(
        baseline.incoming_to_outflow,
        rel=2.0e-12,
        abs=2.0e-13,
    )
    assert scaled.source_to_outflow == pytest.approx(
        baseline.source_to_outflow,
        rel=2.0e-12,
        abs=2.0e-13,
    )
    assert factor * scaled.incoming_to_population == pytest.approx(
        baseline.incoming_to_population,
        rel=2.0e-12,
        abs=2.0e-13,
    )
    assert factor * scaled.source_to_population == pytest.approx(
        baseline.source_to_population,
        rel=2.0e-12,
        abs=2.0e-13,
    )


def test_global_stationary_response_matches_its_explicit_positive_matrix(
    tmp_path: Path,
) -> None:
    config, operator = _build_operator(tmp_path)
    acceleration = _acceleration(config, 10.0)
    growth = 0.0
    system = operator.build_response_system(
        acceleration,
        growth,
        operator.minimum_response_shift(growth),
    )
    matrix = operator.explicit_response_matrix(system)
    rng = np.random.default_rng(29)
    population = rng.random(
        (operator.grid.energy_cells, operator.grid.polar_cells)
    )
    application = system.apply(population)

    assert float(np.min(matrix.data, initial=0.0)) >= -1.0e-13
    assert matrix @ population.reshape(-1) == pytest.approx(
        application.population.reshape(-1),
        rel=2.0e-12,
        abs=5.0e-11,
    )
    assert application.outer_escape_s_inv >= 0.0
    assert application.boundary_residual_L1 <= 1.0e-12


def test_perron_branch_matches_dense_reference_and_is_seed_independent(
    tmp_path: Path,
) -> None:
    config, operator = _build_operator(
        tmp_path,
        energy_cells=16,
        polar_cells=8,
    )
    acceleration = _acceleration(config, 10.0)
    controls = SteadyControls(
        max_iterations=1000,
        tolerance=1.0e-9,
        tail_probability_target=1.0e-8,
    )
    shape = (operator.grid.energy_cells, operator.grid.polar_cells)
    rng = np.random.default_rng(917)
    random_seed = rng.random(shape) + 1.0e-6
    beam_seed = np.zeros(shape, dtype=float)
    beam_seed[-3, 0] = 1.0
    solutions = [
        solve_steady_state(
            operator,
            acceleration,
            controls,
            initial_population=seed,
        )
        for seed in (
            maxwellian_initial_population(operator),
            random_seed,
            beam_seed,
        )
    ]
    reference = solutions[0]
    for solution in solutions[1:]:
        assert solution.growth_frequency_s_inv == pytest.approx(
            reference.growth_frequency_s_inv,
            rel=2.0e-10,
            abs=2.0e-7,
        )
        assert float(
            np.sum(np.abs(solution.population - reference.population))
        ) <= 2.0e-8

    growth = reference.growth_frequency_s_inv
    system = operator.build_response_system(
        acceleration,
        growth,
        operator.minimum_response_shift(growth),
    )
    matrix = operator.explicit_response_matrix(system).toarray()
    eigenvalues, eigenvectors = np.linalg.eig(matrix)
    principal_index = int(np.argmax(np.abs(eigenvalues)))
    principal_value = complex(eigenvalues[principal_index])
    principal = np.real(eigenvectors[:, principal_index])
    if float(np.sum(principal)) < 0.0:
        principal = -principal
    principal /= float(np.sum(principal))
    moduli = np.sort(np.abs(eigenvalues))

    assert abs(principal_value.imag) <= 2.0e-11
    assert principal_value.real == pytest.approx(1.0, abs=2.0e-9)
    assert moduli[-1] - moduli[-2] > 1.0e-7
    assert float(np.min(principal)) >= -2.0e-12
    assert principal == pytest.approx(
        reference.population.reshape(-1),
        rel=2.0e-7,
        abs=2.0e-10,
    )

    bounds = operator.growth_bounds(acceleration)
    displacement = min(
        0.05 * (bounds.upper_s_inv - bounds.lower_s_inv),
        max(1.0, 0.01 * abs(growth)),
    )
    sampled_radii: list[float] = []
    for sampled_growth in (growth - displacement, growth + displacement):
        sampled_growth = float(
            np.clip(
                sampled_growth,
                bounds.lower_s_inv,
                bounds.upper_s_inv,
            )
        )
        sampled_system = operator.build_response_system(
            acceleration,
            sampled_growth,
            operator.minimum_response_shift(sampled_growth),
        )
        sampled_matrix = operator.explicit_response_matrix(
            sampled_system
        ).toarray()
        sampled_radii.append(
            float(np.max(np.abs(np.linalg.eigvals(sampled_matrix))))
        )
    assert sampled_radii[0] > 1.0
    assert sampled_radii[1] < 1.0
