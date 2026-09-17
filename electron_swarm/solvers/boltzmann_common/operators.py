"""Solver-neutral finite-volume energy-transport operators."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse

from electron_swarm.core.constants import E_CHARGE_C, ELECTRON_MASS_KG, EV_TO_J
from electron_swarm.solvers.boltzmann_common.collisions import (
    EffectiveCollisionData,
)
from electron_swarm.solvers.boltzmann_common.grid import KineticGrid


@dataclass(slots=True)
class KineticOperatorBlock:
    grid: KineticGrid
    gas_number_density_m3: float
    electric_field_V_m: float
    collisions: EffectiveCollisionData
    energy_flux_matrix: sparse.csr_matrix
    collision_matrix: sparse.csr_matrix
    matrix: sparse.csr_matrix
    discretization: str = "finite_volume_scharfetter_gummel"

    @property
    def energy_eV(self) -> np.ndarray:
        return self.grid.energy_eV

    @property
    def edges_eV(self) -> np.ndarray:
        return self.grid.edges_eV

    @property
    def widths_eV(self) -> np.ndarray:
        return self.grid.widths_eV


def _bernoulli(x: np.ndarray | float) -> np.ndarray | float:
    x_arr = np.asarray(x, dtype=float)
    out = np.empty_like(x_arr, dtype=float)
    small = np.abs(x_arr) < 1.0e-7
    pos_large = x_arr > 80.0
    neg_large = x_arr < -80.0
    mid = ~(small | pos_large | neg_large)
    out[small] = (
        1.0
        - x_arr[small] / 2.0
        + x_arr[small] ** 2 / 12.0
        - x_arr[small] ** 4 / 720.0
    )
    out[pos_large] = 0.0
    out[neg_large] = -x_arr[neg_large]
    out[mid] = x_arr[mid] / np.expm1(x_arr[mid])
    if np.isscalar(x):
        return float(out)
    return out


def field_diffusion_eV2_s(
    energy: np.ndarray,
    electric_field_V_m: float,
    momentum_frequency_s_inv: np.ndarray,
) -> np.ndarray:
    eps_J = np.maximum(energy * EV_TO_J, 0.0)
    D_J2_s = (
        (2.0 / 3.0)
        * (E_CHARGE_C * electric_field_V_m) ** 2
        / ELECTRON_MASS_KG
        * eps_J
        / np.maximum(momentum_frequency_s_inv, 1.0e-60)
    )
    return D_J2_s / (EV_TO_J * EV_TO_J)


def assemble_energy_flux_operator(
    energy: np.ndarray,
    widths: np.ndarray,
    electric_field_V_m: float,
    collisions: EffectiveCollisionData,
    *,
    momentum_frequency_s_inv: np.ndarray | None = None,
) -> sparse.csr_matrix:
    n = len(energy)
    mat = sparse.lil_matrix((n, n), dtype=float)
    momentum_frequency = (
        collisions.nu_m
        if momentum_frequency_s_inv is None
        else np.asarray(momentum_frequency_s_inv, dtype=float)
    )
    if momentum_frequency.shape != np.asarray(energy).shape:
        raise ValueError("Momentum frequency must match the energy grid")
    D_field = field_diffusion_eV2_s(
        energy, electric_field_V_m, momentum_frequency
    )
    eps_safe = np.maximum(energy, max(energy[1] - energy[0], 1.0e-12) * 0.5)
    A_field = D_field / (2.0 * eps_safe)
    D_center = np.maximum(collisions.elastic_D_eV2_s + D_field, 1.0e-80)
    A_center = collisions.elastic_A_eV_s + A_field
    for i in range(n - 1):
        h = energy[i + 1] - energy[i]
        if h <= 0.0:
            continue
        D = max(0.5 * (D_center[i] + D_center[i + 1]), 1.0e-80)
        A = 0.5 * (A_center[i] + A_center[i + 1])
        peclet = A * h / D
        cL = D / h * _bernoulli(-peclet)
        cR = -D / h * _bernoulli(peclet)
        mat[i, i] += -cL / widths[i]
        mat[i, i + 1] += -cR / widths[i]
        mat[i + 1, i] += cL / widths[i + 1]
        mat[i + 1, i + 1] += cR / widths[i + 1]
    return mat.tocsr()


__all__ = [
    "KineticOperatorBlock",
    "assemble_energy_flux_operator",
    "field_diffusion_eV2_s",
]
