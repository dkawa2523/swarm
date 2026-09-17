"""Physical moments of a normalized propagator cell population."""

from __future__ import annotations

import numpy as np

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.cross_sections import (
    CrossSectionSet,
    ProcessType,
    mixture_fraction,
)
from electron_swarm.core.results import EnergyAngleDistributionResult
from electron_swarm.core.results import RateResult
from electron_swarm.physics.electron_neutral import speed_from_energy_m_s
from electron_swarm.physics.kinetics import RateConvolution, gas_number_density
from electron_swarm.solvers.propagator.collisions import (
    cell_average_collision_rate,
)
from electron_swarm.solvers.propagator.models import (
    CollisionOperatorData,
    PropagatorGrid,
)


def normalized_population(population: np.ndarray) -> np.ndarray:
    values = np.asarray(population, dtype=float)
    if np.any(~np.isfinite(values)):
        raise FloatingPointError("propagator population is non-finite")
    total = float(np.sum(values))
    if total <= 0.0:
        raise FloatingPointError("propagator population has no positive mass")
    return values / total


def eedf_from_population(
    grid: PropagatorGrid,
    population: np.ndarray,
) -> np.ndarray:
    values = normalized_population(population)
    return np.sum(values, axis=1) / grid.energy_widths_eV


def mean_energy_eV(
    grid: PropagatorGrid,
    population: np.ndarray,
) -> float:
    values = normalized_population(population)
    return float(
        np.sum(values * grid.energy_centers_eV[:, None])
    )


def flux_drift_velocity_m_s(
    grid: PropagatorGrid,
    population: np.ndarray,
) -> float:
    values = normalized_population(population)
    edges = np.asarray(grid.energy_edges_eV, dtype=float)
    widths = np.asarray(grid.energy_widths_eV, dtype=float)
    speed_per_sqrt_eV = float(speed_from_energy_m_s(np.asarray([1.0]))[0])
    # The stored energy-angle population represents a piecewise-constant
    # density in energy.  Integrate v(E) exactly in that same basis instead of
    # applying a center sample of sqrt(E), which is biased in the origin cells.
    cell_average_speed = (
        speed_per_sqrt_eV
        * (2.0 / 3.0)
        * (edges[1:] ** 1.5 - edges[:-1] ** 1.5)
        / widths
    )
    return float(
        np.sum(
            values
            * cell_average_speed[:, None]
            * grid.mu_centers[None, :]
        )
    )


def energy_angle_distribution(
    grid: PropagatorGrid,
    population: np.ndarray,
) -> EnergyAngleDistributionResult:
    values = normalized_population(population)
    density = values / (
        grid.energy_widths_eV[:, None] * grid.mu_widths[None, :]
    )
    normalization = float(
        np.sum(
            density
            * grid.energy_widths_eV[:, None]
            * grid.mu_widths[None, :]
        )
    )
    if not np.isclose(normalization, 1.0, rtol=0.0, atol=2.0e-13):
        raise FloatingPointError(
            "energy-angle distribution failed normalization"
        )
    return EnergyAngleDistributionResult(
        energy_eV=grid.energy_centers_eV.copy(),
        energy_widths_eV=grid.energy_widths_eV.copy(),
        mu=grid.mu_centers.copy(),
        mu_widths=grid.mu_widths.copy(),
        density_eV_inv=density,
    )


def elastic_energy_loss_rate_coefficient_eV_m3_s(
    grid: PropagatorGrid,
    collisions: CollisionOperatorData,
    population: np.ndarray,
    gas_number_density_m3: float,
) -> float:
    values = normalized_population(population)
    elastic_action = np.asarray(
        collisions.elastic_energy_generator_s_inv @ values,
        dtype=float,
    )
    energy_change_rate = float(
        np.sum(elastic_action * grid.energy_centers_eV[:, None])
    )
    density = float(gas_number_density_m3)
    if density <= 0.0 or not np.isfinite(density):
        raise ValueError("gas number density must be finite and positive")
    return -energy_change_rate / density


def cell_integrated_rates(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    grid: PropagatorGrid,
    population: np.ndarray,
    *,
    case_id: str,
    e_over_n_Td: float,
    solver_name: str,
) -> RateConvolution:
    """Integrate every reported process with the collision-cell quadrature."""

    shell_mass = np.sum(normalized_population(population), axis=1)
    density = gas_number_density(config)
    rates: list[RateResult] = []
    ionization_rate = 0.0
    attachment_rate = 0.0
    net_frequency = 0.0
    for process in cross_sections.processes:
        fraction = mixture_fraction(config.conditions, process.species)
        if fraction <= 0.0:
            continue
        process_rate_cells = cell_average_collision_rate(
            (process,),
            1.0,
            1.0,
            grid,
        )
        rate_coefficient = float(np.dot(process_rate_cells, shell_mass))
        mixture_rate = fraction * rate_coefficient
        if process.process_type == ProcessType.IONIZATION:
            ionization_rate += mixture_rate
            net_frequency += density * mixture_rate
        elif process.process_type == ProcessType.ATTACHMENT:
            attachment_rate += mixture_rate
            net_frequency -= density * mixture_rate
        energy_loss = 0.0
        if process.process_type in {
            ProcessType.EXCITATION,
            ProcessType.IONIZATION,
        }:
            energy_loss = float(process.threshold_eV)
        elif process.process_type == ProcessType.SUPERELASTIC:
            energy_loss = -abs(float(process.threshold_eV))
        rates.append(
            RateResult(
                solver=solver_name,
                case_id=case_id,
                e_over_n_Td=e_over_n_Td,
                species=process.species,
                process=process.process,
                process_type=process.process_type.value,
                threshold_eV=process.threshold_eV,
                rate_coefficient_m3_s=rate_coefficient,
                target_species_fraction=fraction,
                energy_loss_eV=energy_loss,
                gas_number_density_m3=density,
            )
        )
    return RateConvolution(
        rates,
        ionization_rate,
        attachment_rate,
        net_frequency,
    )
