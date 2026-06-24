"""Energy-space electron-electron Fokker-Planck relaxation helper."""

from __future__ import annotations

import numpy as np

from electron_swarm.collisions.electron_electron import (
    relaxation_target,
)


def _normalize(energy_eV: np.ndarray, widths_eV: np.ndarray, values: np.ndarray) -> np.ndarray:
    if not np.all(np.isfinite(values)):
        raise ValueError("electron_electron fp_energy produced non-finite EEDF")
    clipped = np.clip(np.asarray(values, dtype=float), 0.0, None)
    norm = float(np.sum(clipped * widths_eV))
    if norm <= 0.0 or not np.isfinite(norm):
        raise ValueError("electron_electron fp_energy produced non-normalizable EEDF")
    return clipped / norm


def _solve_tridiagonal(
    lower: np.ndarray,
    diag: np.ndarray,
    upper: np.ndarray,
    rhs: np.ndarray,
) -> np.ndarray:
    n = len(diag)
    a = np.asarray(lower, dtype=float).copy()
    b = np.asarray(diag, dtype=float).copy()
    c = np.asarray(upper, dtype=float).copy()
    d = np.asarray(rhs, dtype=float).copy()
    for i in range(1, n):
        pivot = b[i - 1]
        if abs(pivot) <= 1.0e-300:
            raise ValueError("electron_electron fp_energy operator is singular")
        factor = a[i - 1] / pivot
        b[i] -= factor * c[i - 1]
        d[i] -= factor * d[i - 1]
    out = np.empty(n, dtype=float)
    pivot = b[-1]
    if abs(pivot) <= 1.0e-300:
        raise ValueError("electron_electron fp_energy operator is singular")
    out[-1] = d[-1] / pivot
    for i in range(n - 2, -1, -1):
        pivot = b[i]
        if abs(pivot) <= 1.0e-300:
            raise ValueError("electron_electron fp_energy operator is singular")
        out[i] = (d[i] - c[i] * out[i + 1]) / pivot
    return out


def apply_fp_energy_operator(
    energy_eV: np.ndarray,
    widths_eV: np.ndarray,
    eedf_eV_inv: np.ndarray,
    *,
    relaxation_fraction: float,
    conserve_mean_energy: bool = True,
    fallback_temperature_eV: float = 2.0,
) -> tuple[np.ndarray, dict[str, object]]:
    """Apply a conservative f0-only energy-space FP relaxation step.

    The operator is a zero-flux finite-volume diffusion in ``f / f_eq`` where
    ``f_eq`` is the same Maxwellian-like target used by the product
    relaxation postprocess. It is intentionally limited to the isotropic energy
    distribution; angular e-e damping is not included here.
    """

    energy = np.asarray(energy_eV, dtype=float)
    widths = np.asarray(widths_eV, dtype=float)
    old = np.asarray(eedf_eV_inv, dtype=float)
    if energy.ndim != 1 or widths.shape != energy.shape or old.shape != energy.shape:
        raise ValueError("electron_electron fp_energy arrays must share one shape")
    if len(energy) < 2 or np.any(np.diff(energy) <= 0.0):
        raise ValueError("electron_electron fp_energy requires increasing energy grid")
    if np.any(widths <= 0.0) or not np.all(np.isfinite(energy + widths + old)):
        raise ValueError("electron_electron fp_energy received invalid EEDF inputs")
    alpha = float(relaxation_fraction)
    if not (0.0 <= alpha <= 1.0):
        raise ValueError("electron_electron relaxation_fraction must be in [0, 1]")

    old = _normalize(energy, widths, old)
    target = relaxation_target(
        energy,
        widths,
        old,
        conserve_mean_energy=conserve_mean_energy,
        fallback_temperature_eV=fallback_temperature_eV,
    )
    if alpha == 0.0:
        return old, {}

    floor = max(float(np.max(target)) * 1.0e-14, 1.0e-300)
    target = _normalize(energy, widths, np.maximum(target, floor))
    de = np.diff(energy)
    face_target = 0.5 * (target[:-1] + target[1:])
    conductance = face_target / np.maximum(de, 1.0e-300)
    step = min(alpha / max(1.0 - alpha, 1.0e-6), 1.0e6)
    step *= float(np.median(widths) ** 2)

    n = len(energy)
    lower = np.zeros(n - 1, dtype=float)
    diag = np.ones(n, dtype=float)
    upper = np.zeros(n - 1, dtype=float)
    for i in range(n):
        if i > 0:
            c_left = conductance[i - 1] / widths[i]
            lower[i - 1] = -step * c_left / target[i - 1]
            diag[i] += step * c_left / target[i]
        if i < n - 1:
            c_right = conductance[i] / widths[i]
            upper[i] = -step * c_right / target[i + 1]
            diag[i] += step * c_right / target[i]

    relaxed = _solve_tridiagonal(lower, diag, upper, old)
    new = _normalize(energy, widths, relaxed)
    return new, {}
