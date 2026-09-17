"""Small solver-independent kinetic quantities shared by swarm methods."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.constants import BOLTZMANN_J_K, ELECTRON_MASS_KG, EV_TO_J
from electron_swarm.core.cross_sections import CrossSectionProcess, ProcessType
from electron_swarm.core.results import RateResult


_RATE_GAUSS_NODES, _RATE_GAUSS_WEIGHTS = np.polynomial.legendre.leggauss(3)
_INCIDENT_THRESHOLD_TYPES = frozenset(
    {
        ProcessType.ATTACHMENT,
        ProcessType.EXCITATION,
        ProcessType.IONIZATION,
    }
)


@dataclass(frozen=True, slots=True)
class RateConvolution:
    rates: list[RateResult]
    ionization_rate_m3_s: float
    attachment_rate_m3_s: float
    net_ionization_frequency_s: float


def gas_number_density(config: SwarmConfig) -> float:
    conditions = config.conditions
    if conditions.gas_number_density_m3 is not None:
        return float(conditions.gas_number_density_m3)
    assert conditions.pressure_Pa is not None
    return float(
        conditions.pressure_Pa
        / (BOLTZMANN_J_K * conditions.gas_temperature_K)
    )


def cell_integrated_rate_coefficient(
    energy_eV: np.ndarray,
    energy_width_eV: np.ndarray,
    density_eV_inv: np.ndarray,
    process: CrossSectionProcess,
) -> float:
    """Integrate ``sigma(E) v(E) F(E)`` over finite-volume EEDF cells.

    A histogram represents a cell-average density, not a point value at the
    cell centre.  Evaluating a thresholded cross section only at that centre
    creates an order-one bias when a reaction turns on inside a coarse cell.
    Every cell is therefore split at the process threshold and at all
    piecewise-linear cross-section knots before applying the cached three-point
    Gauss rule.  Quadrature nodes never straddle a discontinuity or a change in
    interpolation slope.  Incident-threshold reactions integrate only the
    physical support at and above that threshold, so narrow near-edge support
    is retained without modifying the histogram or correcting the rate later.
    """

    centers = np.asarray(energy_eV, dtype=float)
    widths = np.asarray(energy_width_eV, dtype=float)
    density = np.asarray(density_eV_inv, dtype=float)
    if (
        centers.ndim != 1
        or widths.shape != centers.shape
        or density.shape != centers.shape
        or np.any(~np.isfinite(centers))
        or np.any(~np.isfinite(widths))
        or np.any(~np.isfinite(density))
        or np.any(centers < 0.0)
        or np.any(widths <= 0.0)
        or np.any(density < 0.0)
    ):
        raise ValueError("rate-convolution EEDF cells are invalid")
    half_widths = 0.5 * widths
    left_edges = centers - half_widths
    scale = np.maximum(1.0, np.abs(centers))
    if np.any(left_edges < -64.0 * np.finfo(float).eps * scale):
        raise ValueError("rate-convolution EEDF cells extend below zero energy")
    left_edges = np.maximum(left_edges, 0.0)
    right_edges = centers + half_widths

    breakpoints = np.asarray(process.energy_eV, dtype=float)
    if process.threshold_eV is not None:
        breakpoints = np.concatenate(
            (breakpoints, np.asarray([float(process.threshold_eV)]))
        )
    breakpoints = np.unique(breakpoints)

    segment_left_parts: list[np.ndarray] = []
    segment_right_parts: list[np.ndarray] = []
    segment_cell_parts: list[np.ndarray] = []
    for cell_index, (left, right) in enumerate(
        zip(left_edges, right_edges, strict=True)
    ):
        first = int(np.searchsorted(breakpoints, left, side="right"))
        stop = int(np.searchsorted(breakpoints, right, side="left"))
        bounds = np.concatenate(
            (
                np.asarray([left]),
                breakpoints[first:stop],
                np.asarray([right]),
            )
        )
        segment_left_parts.append(bounds[:-1])
        segment_right_parts.append(bounds[1:])
        segment_cell_parts.append(
            np.full(bounds.size - 1, cell_index, dtype=np.intp)
        )

    segment_left = np.concatenate(segment_left_parts)
    segment_right = np.concatenate(segment_right_parts)
    segment_cells = np.concatenate(segment_cell_parts)
    if (
        process.threshold_eV is not None
        and process.process_type in _INCIDENT_THRESHOLD_TYPES
    ):
        active = segment_right > float(process.threshold_eV)
        segment_left = segment_left[active]
        segment_right = segment_right[active]
        segment_cells = segment_cells[active]
        if segment_left.size == 0:
            return 0.0
    segment_half_width = 0.5 * (segment_right - segment_left)
    segment_center = segment_left + segment_half_width
    nodes = (
        segment_center[:, None]
        + segment_half_width[:, None] * _RATE_GAUSS_NODES[None, :]
    )
    sigma = np.asarray(process.sigma(nodes.reshape(-1)), dtype=float).reshape(
        nodes.shape
    )
    if sigma.shape != nodes.shape or np.any(~np.isfinite(sigma)) or np.any(sigma < 0.0):
        raise ValueError("rate-convolution cross section is invalid")
    speed = np.sqrt(2.0 * EV_TO_J * nodes / ELECTRON_MASS_KG)
    segment_integrals = segment_half_width * np.sum(
        _RATE_GAUSS_WEIGHTS[None, :] * sigma * speed,
        axis=1,
    )
    return float(np.sum(density[segment_cells] * segment_integrals))


__all__ = (
    "RateConvolution",
    "cell_integrated_rate_coefficient",
    "gas_number_density",
)
