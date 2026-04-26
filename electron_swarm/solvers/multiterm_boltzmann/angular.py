"""Legendre basis helpers for the multi-term Boltzmann solver."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse


@dataclass(frozen=True, slots=True)
class LegendreBasis:
    lmax: int

    def __post_init__(self) -> None:
        if self.lmax < 1:
            raise ValueError("lmax must be >= 1")

    @property
    def n_terms(self) -> int:
        return self.lmax + 1

    def flat_index(self, ell: int, i_energy: int, n_energy: int) -> int:
        if ell < 0 or ell > self.lmax:
            raise IndexError("ell out of range")
        return ell * n_energy + i_energy

    def flatten(self, coeff: np.ndarray) -> np.ndarray:
        coeff = np.asarray(coeff)
        if coeff.shape[0] != self.n_terms:
            raise ValueError("coefficient array first dimension must be lmax+1")
        return coeff.reshape(-1)

    def unflatten(self, flat: np.ndarray, n_energy: int) -> np.ndarray:
        return np.asarray(flat).reshape(self.n_terms, n_energy)

    def density_weights(self, energy_weights: np.ndarray) -> np.ndarray:
        energy_weights = np.asarray(energy_weights, dtype=float)
        out = np.zeros(self.n_terms * len(energy_weights), dtype=float)
        out[: len(energy_weights)] = energy_weights
        return out


def mu_coupling_coefficients(ell: int):
    """Coefficients of mu P_l = c_minus P_{l-1} + c_plus P_{l+1}."""

    c_minus = ell / (2 * ell + 1) if ell > 0 else 0.0
    c_plus = (ell + 1) / (2 * ell + 1)
    return (ell - 1 if ell > 0 else None, c_minus), (ell + 1, c_plus)


def build_mu_multiplication_operator(
    basis: LegendreBasis, n_energy: int
) -> sparse.csr_matrix:
    rows = []
    cols = []
    vals = []
    for ell in range(basis.n_terms):
        (lm, cm), (lp, cp) = mu_coupling_coefficients(ell)
        for i in range(n_energy):
            col = basis.flat_index(ell, i, n_energy)
            if lm is not None and lm <= basis.lmax:
                rows.append(basis.flat_index(lm, i, n_energy))
                cols.append(col)
                vals.append(cm)
            if lp <= basis.lmax:
                rows.append(basis.flat_index(lp, i, n_energy))
                cols.append(col)
                vals.append(cp)
    n = basis.n_terms * n_energy
    return sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))
