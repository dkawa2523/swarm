"""Finite-k hydrodynamic transport extraction for the multi-term operator."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as spla

from electron_swarm.core.transport import (
    BulkTransport,
    FluxTransport,
    TransportMetadata,
    TransportSet,
)


@dataclass(frozen=True, slots=True)
class HydroConfig:
    component: Literal["longitudinal", "transverse"] = "longitudinal"
    k_values: tuple[float, ...] = (-2e-4, -1e-4, 0.0, 1e-4, 2e-4)
    fit_order: int = 2
    dense_threshold: int = 256


@dataclass(frozen=True, slots=True)
class HydroFit:
    k_values: np.ndarray
    eigenvalues: np.ndarray
    velocities: np.ndarray
    nu_eff_s_inv: float
    W_bulk_m_s: float
    D_bulk_m2_s: float | None
    W_flux_m_s: float
    D_flux_m2_s: float | None
    fit_residual: float
    symmetry_error: float
    mode_continuity_error: float


def _dominant_eigenpair(A, dense_threshold: int = 256, sigma: complex | None = None):
    n = A.shape[0]
    target = 0.0 + 0.0j if sigma is None else sigma
    if sparse.issparse(A):
        if n <= dense_threshold:
            dense = A.toarray()
        else:
            try:
                vals, vecs = spla.eigs(
                    A, k=1, sigma=target, which="LM"
                )
                return vals[0], vecs[:, 0]
            except Exception:
                vals, vecs = spla.eigs(A, k=1, which="LR")
                return vals[0], vecs[:, 0]
    else:
        dense = np.asarray(A, dtype=complex)
    vals, vecs = np.linalg.eig(dense)
    idx = int(np.argmin(np.abs(vals - target)))
    return vals[idx], vecs[:, idx]


def fit_lambda_dispersion(k: np.ndarray, lam: np.ndarray, fit_order: int = 2):
    k = np.asarray(k, dtype=float)
    lam = np.asarray(lam, dtype=complex)
    deg = max(1, min(int(fit_order), len(k) - 1))
    coef_real = np.polyfit(k, lam.real, deg=deg)
    coef_imag = np.polyfit(k, lam.imag, deg=deg)
    cr = coef_real[::-1]
    ci = coef_imag[::-1]
    nu = float(cr[0])
    W = float(-ci[1]) if len(ci) > 1 else 0.0
    D = float(-cr[2]) if len(cr) > 2 else None
    pred = np.polyval(coef_real, k) + 1j * np.polyval(coef_imag, k)
    residual = float(np.linalg.norm(pred - lam) / max(np.linalg.norm(lam), 1e-300))
    return nu, W, D, cr + 1j * ci, residual


def fit_flux_velocity(k: np.ndarray, velocity: np.ndarray, fit_order: int = 1):
    k = np.asarray(k, dtype=float)
    v = np.asarray(velocity, dtype=complex)
    deg = max(1, min(int(fit_order), len(k) - 1))
    coef_real = np.polyfit(k, v.real, deg=deg)
    coef_imag = np.polyfit(k, v.imag, deg=deg)
    cr = coef_real[::-1]
    ci = coef_imag[::-1]
    W = float(cr[0])
    D = float(-ci[1]) if len(ci) > 1 else None
    pred = np.polyval(coef_real, k) + 1j * np.polyval(coef_imag, k)
    residual = float(np.linalg.norm(pred - v) / max(np.linalg.norm(v), 1e-300))
    return W, D, cr + 1j * ci, residual


def _symmetric_pair_error(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=complex)
    if len(values) < 3:
        return 0.0
    errs = []
    half = len(values) // 2
    for i in range(half):
        a = values[i]
        b = np.conjugate(values[-(i + 1)])
        errs.append(abs(a - b) / max(abs(a), abs(b), 1.0))
    return float(max(errs, default=0.0))


def _mode_continuity_error(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=complex)
    if len(values) < 3:
        return 0.0
    jumps = np.abs(np.diff(values))
    scale = max(float(np.max(np.abs(values))), 1.0)
    return float(np.max(jumps) / scale)


def _tracking_order(k_values: np.ndarray) -> list[int]:
    return sorted(range(len(k_values)), key=lambda idx: (abs(k_values[idx]), k_values[idx]))


class HydrodynamicModeSolver:
    """Finite-k hydrodynamic mode solver for flux/bulk transport."""

    def solve(
        self,
        base_operator,
        streaming_operator,
        velocity_moment,
        electric_field_V_m: float,
        normalization_weights,
        ionization_frequency_s_inv: float = 0.0,
        attachment_frequency_s_inv: float = 0.0,
        config: HydroConfig | None = None,
    ):
        config = config or HydroConfig()
        L0 = (
            sparse.csr_matrix(base_operator)
            if sparse.issparse(base_operator)
            else np.asarray(base_operator, dtype=complex)
        )
        S = (
            sparse.csr_matrix(streaming_operator)
            if sparse.issparse(streaming_operator)
            else np.asarray(streaming_operator, dtype=complex)
        )
        v_moment = np.asarray(velocity_moment, dtype=complex)
        norm_w = np.asarray(normalization_weights, dtype=complex)
        k_values = np.asarray(config.k_values, dtype=float)
        eigvals = np.empty(len(k_values), dtype=complex)
        velocities = np.empty(len(k_values), dtype=complex)
        solved: list[int] = []
        for idx in _tracking_order(k_values):
            k = k_values[idx]
            if solved:
                nearest = min(solved, key=lambda j: abs(k_values[j] - k))
                target = eigvals[nearest]
            else:
                target = None
            A = L0 - 1j * float(k) * S
            val, vec = _dominant_eigenpair(
                A, dense_threshold=config.dense_threshold, sigma=target
            )
            den = np.vdot(norm_w, vec)
            if abs(den) < 1e-300:
                den = 1.0
            vec = vec / den
            eigvals[idx] = val
            velocities[idx] = np.vdot(v_moment, vec)
            solved.append(idx)
        if not np.all(np.isfinite(eigvals)) or not np.all(np.isfinite(velocities)):
            raise RuntimeError("hydrodynamic finite-k solve produced non-finite modes")
        symmetry_error = max(
            _symmetric_pair_error(eigvals), _symmetric_pair_error(velocities)
        )
        continuity_error = max(
            _mode_continuity_error(eigvals), _mode_continuity_error(velocities)
        )
        nu, Wb, Db, _, res_lam = fit_lambda_dispersion(
            k_values, eigvals, config.fit_order
        )
        Wf, Df, _, res_v = fit_flux_velocity(
            k_values, velocities, min(1, config.fit_order)
        )
        residual = max(res_lam, res_v)
        if not np.isfinite(residual) or residual > 5.0e-1:
            raise RuntimeError(
                f"hydrodynamic finite-k fit residual is too large ({residual:.3e})"
            )
        if symmetry_error > 5.0e-1:
            raise RuntimeError(
                "hydrodynamic finite-k mode symmetry check failed "
                f"({symmetry_error:.3e})"
            )
        if continuity_error > 5.0e-1:
            raise RuntimeError(
                "hydrodynamic finite-k mode continuity check failed "
                f"({continuity_error:.3e})"
            )
        scale = max(abs(Wf), abs(Wb), 1.0)
        if Df is not None and Df < -1.0e-8 * scale:
            raise RuntimeError(f"flux diffusion from finite-k fit is negative ({Df:.3e})")
        if Db is not None and Db < -1.0e-8 * scale:
            raise RuntimeError(f"bulk diffusion from finite-k fit is negative ({Db:.3e})")
        if Df is not None and Df < 0.0:
            Df = 0.0
        if Db is not None and Db < 0.0:
            Db = 0.0
        fit = HydroFit(
            k_values,
            eigvals,
            velocities,
            nu,
            Wb,
            Db,
            Wf,
            Df,
            residual,
            symmetry_error,
            continuity_error,
        )
        flux = FluxTransport.from_drift_and_field(
            Wf,
            electric_field_V_m,
            diffusion_longitudinal_m2_s=Df
            if config.component == "longitudinal"
            else None,
            diffusion_transverse_m2_s=Df if config.component == "transverse" else None,
        )
        bulk = BulkTransport.from_drift_and_field(
            Wb,
            electric_field_V_m,
            diffusion_longitudinal_m2_s=Db
            if config.component == "longitudinal"
            else None,
            diffusion_transverse_m2_s=Db if config.component == "transverse" else None,
        )
        transport = TransportSet.from_flux_bulk(
            flux,
            bulk,
            ionization_frequency_s_inv=ionization_frequency_s_inv,
            attachment_frequency_s_inv=attachment_frequency_s_inv,
            metadata=TransportMetadata(
                solver="multiterm_boltzmann", swarm_condition="hydrodynamic"
            ),
        )
        return fit, transport
