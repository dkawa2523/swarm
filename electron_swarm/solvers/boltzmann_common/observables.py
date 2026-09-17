"""Distribution moments, transport primitives, and rate projection."""

from __future__ import annotations

import numpy as np

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.constants import E_CHARGE_C, ELECTRON_MASS_KG
from electron_swarm.core.cross_sections import (
    CrossSectionSet,
    ProcessType,
    mixture_fraction,
)
from electron_swarm.core.numerics import f0_from_eedf
from electron_swarm.core.results import RateResult
from electron_swarm.physics.kinetics import (
    RateConvolution,
    cell_integrated_rate_coefficient,
    gas_number_density,
)
from electron_swarm.solvers.boltzmann_common.grid import cell_edges_from_centers


def weighted_integral(values: np.ndarray, widths: np.ndarray) -> float:
    return float(np.sum(np.asarray(values, dtype=float) * widths))


def normalize_eedf(eedf: np.ndarray, widths: np.ndarray) -> np.ndarray:
    values = np.asarray(eedf, dtype=float)
    widths = np.asarray(widths, dtype=float)
    if (
        values.shape != widths.shape
        or not np.all(np.isfinite(values))
        or not np.all(np.isfinite(widths))
        or np.any(widths <= 0.0)
    ):
        raise FloatingPointError("Cannot normalize non-finite EEDF")
    total = weighted_integral(values, widths)
    if not np.isfinite(total) or abs(total) < 1.0e-300:
        raise FloatingPointError("Cannot normalize EEDF with zero/non-finite integral")
    if total < 0.0:
        values = -values
        total = -total
    return values / total


def mean_energy_from_eedf(
    energy_eV: np.ndarray,
    widths_eV: np.ndarray,
    eedf: np.ndarray,
) -> float:
    energy = np.asarray(energy_eV, dtype=float)
    widths = np.asarray(widths_eV, dtype=float)
    values = np.asarray(eedf, dtype=float)
    if (
        energy.shape != widths.shape
        or values.shape != widths.shape
        or not np.all(np.isfinite(energy))
        or not np.all(np.isfinite(widths))
        or not np.all(np.isfinite(values))
        or np.any(widths <= 0.0)
    ):
        raise ValueError("Energy, widths, and EEDF arrays must have matching shape")
    return weighted_integral(energy * values, widths)


def nonuniform_center_gradient(values: np.ndarray, energy_eV: np.ndarray) -> np.ndarray:
    """Second-order center gradient for nonuniform energy centers."""

    values = np.asarray(values, dtype=float)
    energy = np.asarray(energy_eV, dtype=float)
    if values.shape != energy.shape or values.ndim != 1:
        raise ValueError("Gradient values and energy grid must be matching 1-D arrays")
    n = len(energy)
    if n < 2:
        return np.zeros_like(values)
    if np.any(np.diff(energy) <= 0.0):
        raise ValueError("Energy grid must be strictly increasing")
    if n == 2:
        slope = (values[1] - values[0]) / (energy[1] - energy[0])
        return np.array([slope, slope], dtype=float)

    gradient = np.empty_like(values)
    for i in range(1, n - 1):
        h0 = energy[i] - energy[i - 1]
        h1 = energy[i + 1] - energy[i]
        gradient[i] = (
            h0 * h0 * values[i + 1]
            + (h1 * h1 - h0 * h0) * values[i]
            - h1 * h1 * values[i - 1]
        ) / (h0 * h1 * (h0 + h1))

    h0 = energy[1] - energy[0]
    h1 = energy[2] - energy[1]
    gradient[0] = (
        -(2.0 * h0 + h1) * values[0] / (h0 * (h0 + h1))
        + (h0 + h1) * values[1] / (h0 * h1)
        - h0 * values[2] / (h1 * (h0 + h1))
    )

    h0 = energy[-2] - energy[-3]
    h1 = energy[-1] - energy[-2]
    gradient[-1] = (
        h1 * values[-3] / (h0 * (h0 + h1))
        - (h0 + h1) * values[-2] / (h0 * h1)
        + (h0 + 2.0 * h1) * values[-1] / (h1 * (h0 + h1))
    )
    return gradient


def negative_mass_fraction(eedf: np.ndarray, widths: np.ndarray) -> float:
    values = np.asarray(eedf, dtype=float)
    widths = np.asarray(widths, dtype=float)
    if (
        values.shape != widths.shape
        or not np.all(np.isfinite(values))
        or not np.all(np.isfinite(widths))
        or np.any(widths <= 0.0)
    ):
        raise FloatingPointError("Cannot evaluate negative mass for non-finite EEDF")
    negative = float(np.sum(np.maximum(-values, 0.0) * widths))
    mass = float(np.sum(np.abs(values) * widths))
    return negative / max(mass, 1.0e-300)


def f0_reduced_transport_from_eedf(
    energy: np.ndarray,
    widths: np.ndarray,
    eedf: np.ndarray,
    sigma_m_m2: np.ndarray,
) -> tuple[float, float, float, float]:
    energy = np.asarray(energy, dtype=float)
    widths = np.asarray(widths, dtype=float)
    eedf = np.asarray(eedf, dtype=float)
    sigma_m = np.asarray(sigma_m_m2, dtype=float)
    if energy.shape != widths.shape or eedf.shape != energy.shape:
        raise ValueError("Energy, widths, and EEDF arrays must have matching shape")
    if sigma_m.shape != energy.shape:
        raise ValueError("Momentum cross section must match the energy grid")
    if not np.all(np.isfinite(sigma_m)) or np.any(sigma_m <= 0.0):
        raise FloatingPointError(
            "Momentum cross section must be finite and positive for flux transport"
        )
    return f0_reduced_transport_from_inverse_sigma(
        energy,
        widths,
        eedf,
        1.0 / sigma_m,
    )


def f0_reduced_transport_from_inverse_sigma(
    energy: np.ndarray,
    widths: np.ndarray,
    eedf: np.ndarray,
    inverse_sigma_m2_inv: np.ndarray,
) -> tuple[float, float, float, float]:
    inverse_sigma = np.asarray(inverse_sigma_m2_inv, dtype=float)
    if inverse_sigma.shape != energy.shape:
        raise ValueError("Inverse momentum cross section must match the energy grid")
    if not np.all(np.isfinite(inverse_sigma)) or np.any(inverse_sigma < 0.0):
        raise FloatingPointError(
            "Inverse momentum cross section must be finite and nonnegative"
        )
    mean_energy = mean_energy_from_eedf(energy, widths, eedf)
    if not np.isfinite(mean_energy) or mean_energy <= 0.0:
        raise ValueError("Mean energy must be positive for energy transport")

    f0 = f0_from_eedf(energy, eedf)
    dF0 = nonuniform_center_gradient(f0, energy)
    gamma = float(np.sqrt(2.0 * E_CHARGE_C / ELECTRON_MASS_KG))
    mobility = -gamma / 3.0 * weighted_integral(
        energy * inverse_sigma * dF0, widths
    )
    diffusion = gamma / 3.0 * weighted_integral(
        energy * inverse_sigma * f0, widths
    )
    energy_mobility = (
        -gamma
        / (3.0 * mean_energy)
        * weighted_integral(energy * energy * inverse_sigma * dF0, widths)
    )
    energy_diffusion = (
        gamma
        / (3.0 * mean_energy)
        * weighted_integral(energy * energy * inverse_sigma * f0, widths)
    )
    return (
        float(mobility),
        float(diffusion),
        float(energy_mobility),
        float(energy_diffusion),
    )


def compute_rates_from_eedf(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    energy: np.ndarray,
    widths: np.ndarray,
    eedf: np.ndarray,
    *,
    case_id: str,
    e_over_n_Td: float,
    solver_name: str,
) -> RateConvolution:
    N = gas_number_density(config)
    edges, reconstructed_widths = cell_edges_from_centers(energy)
    if not np.allclose(reconstructed_widths, widths, rtol=1.0e-12, atol=0.0):
        raise ValueError("Rate-convolution widths do not match the energy grid")
    cell_midpoints = 0.5 * (edges[:-1] + edges[1:])
    rates: list[RateResult] = []
    ion_rate = 0.0
    attach_rate = 0.0
    net_freq = 0.0
    for proc in cross_sections.processes:
        frac = mixture_fraction(config.conditions, proc.species)
        if frac <= 0.0:
            continue
        values = np.asarray(eedf, dtype=float)
        if np.any(values < 0.0):
            # PN solves can retain roundoff-scale signed tail values.  The
            # convolution is linear, so preserve them exactly instead of
            # clipping the solved distribution or weakening the histogram
            # validation in the solver-neutral integral.
            k = cell_integrated_rate_coefficient(
                cell_midpoints, widths, np.maximum(values, 0.0), proc
            ) - cell_integrated_rate_coefficient(
                cell_midpoints, widths, np.maximum(-values, 0.0), proc
            )
        else:
            k = cell_integrated_rate_coefficient(
                cell_midpoints, widths, values, proc
            )
        kmix = frac * k
        if proc.process_type == ProcessType.IONIZATION:
            ion_rate += kmix
            net_freq += N * kmix
        elif proc.process_type == ProcessType.ATTACHMENT:
            attach_rate += kmix
            net_freq -= N * kmix
        loss = 0.0
        if proc.process_type in {ProcessType.EXCITATION, ProcessType.IONIZATION}:
            loss = float(proc.threshold_eV)
        elif proc.process_type == ProcessType.SUPERELASTIC:
            loss = -abs(float(proc.threshold_eV))
        rates.append(
            RateResult(
                solver=solver_name,
                case_id=case_id,
                e_over_n_Td=e_over_n_Td,
                species=proc.species,
                process=proc.process,
                process_type=proc.process_type.value,
                threshold_eV=proc.threshold_eV,
                rate_coefficient_m3_s=k,
                target_species_fraction=frac,
                energy_loss_eV=loss,
                gas_number_density_m3=N,
            )
        )
    return RateConvolution(rates, ion_rate, attach_rate, net_freq)


__all__ = [
    "compute_rates_from_eedf",
    "f0_reduced_transport_from_eedf",
    "f0_reduced_transport_from_inverse_sigma",
    "mean_energy_from_eedf",
    "negative_mass_fraction",
    "nonuniform_center_gradient",
    "normalize_eedf",
    "weighted_integral",
]
