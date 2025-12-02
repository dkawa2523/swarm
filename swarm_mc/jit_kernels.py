# Copyright (c) 2020-2021 ETH Zurich
"""
Optional numba JIT kernels for hot paths. Falls back to pure NumPy when numba
is unavailable or disabled via Config.use_jit.
"""

from __future__ import annotations

import math
import numpy as np
import scipy.constants as csts

try:
    import numba as nb  # type: ignore
    JIT_AVAILABLE = True
except ImportError:  # pragma: no cover - optional dependency
    nb = None
    JIT_AVAILABLE = False


def _maybe_jit(**kwargs):
    """
    Decorator factory that returns numba.njit when available, else identity.
    """

    def _decorator(fn):
        if JIT_AVAILABLE:
            return nb.njit(**kwargs)(fn)  # pragma: no cover - exercised when numba present
        return fn

    return _decorator


_ENERGY_FACTOR = 0.5 * csts.electron_mass / csts.elementary_charge
_VEL_FACTOR = math.sqrt(2.0 * csts.elementary_charge / csts.electron_mass)


@_maybe_jit(fastmath=True, cache=True)
def determine_timestep_jit(rand: float,
                           max_velocity: float,
                           max_acceleration: float,
                           max_coll_freq: np.ndarray,
                           max_coll_period: np.ndarray,
                           max_coll_period_squared: np.ndarray,
                           max_energy_grid: np.ndarray,
                           max_cross_section_energy: float) -> (float, float):
    """
    JIT-friendly timestep determination. Mirrors MonteCarlo.determine_timestep.

    Args:
        rand: uniform(0,1) random number
        max_velocity: max electron velocity (m/s)
        max_acceleration: max electron acceleration (m/s^2)
        max_coll_freq: cumulative maximum collision frequency (array)
        max_coll_period: helper array (same length as max_coll_freq)
        max_coll_period_squared: helper array (same length as max_coll_freq)
        max_energy_grid: energy grid corresponding to max_coll_freq
        max_cross_section_energy: clip ceiling for energy
    Returns:
        (dt, trial_coll_freq)
    """
    s = -math.log(rand)
    max_energy = _ENERGY_FACTOR * max_velocity * max_velocity

    de = 2.0 * max_velocity * max_acceleration * s * max_coll_period \
        + max_acceleration * max_acceleration * s * s * max_coll_period_squared
    e_end = max_energy + de

    # clip to available energy range
    if max_cross_section_energy > 0.0:
        for i in range(e_end.size):
            if e_end[i] > max_cross_section_energy:
                e_end[i] = max_cross_section_energy

    freq_needed = 0.0
    for i in range(e_end.size):
        val = e_end[i]
        idx = np.searchsorted(max_energy_grid, val, side="right")
        if idx <= 0:
            idx = 1
        elif idx >= max_energy_grid.size:
            idx = max_energy_grid.size - 1
        e_low = max_energy_grid[idx - 1]
        e_high = max_energy_grid[idx]
        denom = e_high - e_low
        w = 0.0
        if denom != 0.0:
            w = (val - e_low) / denom
        freq_val = max_coll_freq[idx - 1] + (max_coll_freq[idx] - max_coll_freq[idx - 1]) * w
        if freq_val > freq_needed:
            freq_needed = freq_val

    idx = np.searchsorted(max_coll_freq, freq_needed, side="left")
    if idx < 0:
        idx = 0
    if idx >= max_coll_freq.size:
        idx = max_coll_freq.size - 1

    trial_coll_freq = 1.01 * max_coll_freq[idx]
    return s / float(trial_coll_freq), float(trial_coll_freq)


@_maybe_jit(fastmath=True, cache=True, parallel=True)
def _scattering_angles_jit(rand_phi: np.ndarray, rand_chi: np.ndarray,
                           energy: np.ndarray, iso: bool,
                           e_low: float, e_high: float):
    """
    Generate scattering angles; inputs are pre-sampled random arrays.
    """
    n = energy.size
    phi = np.empty(n, dtype=np.float64)
    sin_phi = np.empty(n, dtype=np.float64)
    cos_phi = np.empty(n, dtype=np.float64)
    cos_chi = np.empty(n, dtype=np.float64)
    sin_chi = np.empty(n, dtype=np.float64)

    for i in range(n):
        ph = 2.0 * math.pi * rand_phi[i]
        phi[i] = ph
        sin_phi[i] = math.sin(ph)
        cos_phi[i] = math.cos(ph)

        if iso:
            cos_val = 1.0 - 2.0 * rand_chi[i]
        else:
            e_safe = energy[i]
            if e_safe < 1e-9:
                e_safe = 1e-9
            anisotropic = (2.0 + e_safe - 2.0 * math.pow(1.0 + e_safe, rand_chi[i])) / e_safe
            isotropic = 1.0 - 2.0 * rand_chi[i]
            blend = (energy[i] - e_low) / (e_high - e_low)
            if blend < 0.0:
                blend = 0.0
            elif blend > 1.0:
                blend = 1.0
            cos_val = (1.0 - blend) * isotropic + blend * anisotropic
        if cos_val < -1.0:
            cos_val = -1.0
        elif cos_val > 1.0:
            cos_val = 1.0
        cos_chi[i] = cos_val
        sin_chi[i] = math.sqrt(max(0.0, 1.0 - cos_val * cos_val))
    return cos_chi, sin_chi, cos_phi, sin_phi


@_maybe_jit(fastmath=True, cache=True, parallel=True)
def unit_scattered_velocity_jit(energy: np.ndarray, velocity: np.ndarray,
                                iso: bool, rand_phi: np.ndarray,
                                rand_chi: np.ndarray,
                                e_low: float, e_high: float):
    """
    JIT variant of unit_scattered_velocity using pre-sampled random arrays.
    """
    cos_chi, sin_chi, cos_phi, sin_phi = _scattering_angles_jit(
        rand_phi, rand_chi, energy, iso, e_low, e_high
    )
    n = energy.size
    v_new_dir = np.empty_like(velocity)

    for i in range(n):
        vx = velocity[0, i]
        vy = velocity[1, i]
        vz = velocity[2, i]
        norm = math.sqrt(vx * vx + vy * vy + vz * vz)
        if norm == 0.0:
            norm = 1e-30
        v_hat_x = vx / norm
        v_hat_y = vy / norm
        v_hat_z = vz / norm

        sin_theta = math.sqrt(max(0.0, 1.0 - v_hat_x * v_hat_x))
        use_alt = sin_theta < 1e-12
        if use_alt:
            sin_theta = 1.0
            ex_x, ex_y, ex_z = 0.0, 1.0, 0.0
            cross1_x, cross1_y, cross1_z = -v_hat_z, 0.0, v_hat_x
            cross2_x = -v_hat_y * v_hat_x
            cross2_y = v_hat_z * v_hat_z + v_hat_x * v_hat_x
            cross2_z = -v_hat_y * v_hat_z
        else:
            ex_x, ex_y, ex_z = 1.0, 0.0, 0.0
            cross1_x, cross1_y, cross1_z = 0.0, v_hat_z, -v_hat_y
            cross2_x = v_hat_y * v_hat_y + v_hat_z * v_hat_z
            cross2_y = -v_hat_x * v_hat_y
            cross2_z = -v_hat_x * v_hat_z

        scale = sin_chi[i] / sin_theta
        v_new_dir[0, i] = v_hat_x * cos_chi[i] \
            + cross1_x * scale * sin_phi[i] + cross2_x * scale * cos_phi[i]
        v_new_dir[1, i] = v_hat_y * cos_chi[i] \
            + cross1_y * scale * sin_phi[i] + cross2_y * scale * cos_phi[i]
        v_new_dir[2, i] = v_hat_z * cos_chi[i] \
            + cross1_z * scale * sin_phi[i] + cross2_z * scale * cos_phi[i]

    return v_new_dir, cos_chi


@_maybe_jit(fastmath=True, cache=True)
def cos_chi_inverse_cdf_lut(n: int = 1024) -> np.ndarray:
    """
    Precompute inverse CDF for isotropic cos_chi distribution.
    """
    # For isotropic scattering: cos_chi = 1 - 2*u, u~U(0,1)
    u = np.linspace(0.0, 1.0, n, dtype=np.float64)
    return 1.0 - 2.0 * u


@_maybe_jit(fastmath=True, cache=True)
def sample_cos_chi_from_lut(lut: np.ndarray, rand: np.ndarray) -> np.ndarray:
    """
    Sample cos_chi using precomputed LUT and random array in [0,1).
    """
    # lut is linear in u, but we keep interp for generality/consistency
    idx = rand * (lut.size - 1)
    left = np.floor(idx).astype(np.int64)
    right = np.minimum(left + 1, lut.size - 1)
    w = idx - left
    return lut[left] * (1.0 - w) + lut[right] * w


@_maybe_jit(fastmath=True, cache=True)
def build_aniso_cos_chi_lut(e_grid: np.ndarray, n_bins: int = 1024) -> np.ndarray:
    """
    Build a LUT of inverse CDF for anisotropic scattering per energy in e_grid.
    Shape: (len(e_grid), n_bins)
    """
    # Vahedi anisotropic inverse: cos_chi = (2+E) - 2*(1+E)^u / E
    # We tabulate for u in [0,1]
    u = np.linspace(0.0, 1.0, n_bins, dtype=np.float64)
    lut = np.empty((e_grid.size, n_bins), dtype=np.float64)
    for i in range(e_grid.size):
        e = e_grid[i]
        if e < 1e-9:
            e = 1e-9
        lut[i, :] = (2.0 + e - 2.0 * np.power(1.0 + e, u)) / e
        lut[i, :] = np.clip(lut[i, :], -1.0, 1.0)
    return lut


@_maybe_jit(fastmath=True, cache=True, parallel=True)
def sample_aniso_cos_chi_from_lut(lut: np.ndarray, e_grid: np.ndarray,
                                  energy: np.ndarray, rand: np.ndarray) -> np.ndarray:
    """
    Sample anisotropic cos_chi using precomputed LUT over energy grid.
    """
    out = np.empty_like(energy)
    n_bins = lut.shape[1]
    for i in nb.prange(energy.size):
        e = energy[i]
        idx = np.searchsorted(e_grid, e, side="right")
        if idx <= 0:
            idx = 1
        elif idx >= e_grid.size:
            idx = e_grid.size - 1
        e0 = e_grid[idx - 1]
        e1 = e_grid[idx]
        t = 0.0
        denom = e1 - e0
        if denom != 0.0:
            t = (e - e0) / denom
        u_idx = rand[i] * (n_bins - 1)
        u_l = int(math.floor(u_idx))
        u_r = min(u_l + 1, n_bins - 1)
        w = u_idx - u_l
        c0 = lut[idx - 1, u_l] * (1.0 - w) + lut[idx - 1, u_r] * w
        c1 = lut[idx, u_l] * (1.0 - w) + lut[idx, u_r] * w
        out[i] = c0 * (1.0 - t) + c1 * t
        if out[i] < -1.0:
            out[i] = -1.0
        elif out[i] > 1.0:
            out[i] = 1.0
    return out


@_maybe_jit(cache=True, parallel=True)
def classify_collisions(collisions: np.ndarray,
                        is_attachment: np.ndarray,
                        is_ionization: np.ndarray) -> (np.ndarray, np.ndarray):
    """
    Classify collisions into attachment and ionization masks.
    """
    n = collisions.size
    attach = np.empty(n, dtype=np.bool_)
    ion = np.empty(n, dtype=np.bool_)
    for i in nb.prange(n):
        idx = int(collisions[i])
        attach[i] = is_attachment[idx]
        ion[i] = is_ionization[idx]
    return attach, ion


@_maybe_jit(fastmath=True, cache=True, parallel=True)
def energy_after_collision_jit(energy: np.ndarray,
                               thresholds: np.ndarray,
                               mass_ratios: np.ndarray,
                               cos_chi: np.ndarray) -> np.ndarray:
    """
    Compute post-collision energies with losses, clipped to zero.
    """
    n = energy.size
    out = np.empty_like(energy)
    for i in nb.prange(n):
        loss = thresholds[i] + energy[i] * mass_ratios[i] * (1.0 - cos_chi[i])
        val = energy[i] - loss
        if val < 0.0:
            val = 0.0
        out[i] = val
    return out


@_maybe_jit(cache=True, parallel=True)
def assemble_particles(null_p: np.ndarray, null_v: np.ndarray,
                       surv_p: np.ndarray, surv_v: np.ndarray,
                       new_p: np.ndarray, new_v: np.ndarray,
                       out_p: np.ndarray, out_v: np.ndarray) -> None:
    """
    Copy particles into preallocated output buffers.

    Args:
        null_p/null_v: particles with null collisions, shape (3, n_null)
        surv_p/surv_v: surviving/scattered particles, shape (3, n_survive)
        new_p/new_v: newly created electrons, shape (3, n_new)
        out_p/out_v: preallocated (3, total) buffers
    """
    cursor = 0
    # null collisions
    n_null = null_p.shape[1]
    if n_null:
        out_p[:, cursor:cursor + n_null] = null_p
        out_v[:, cursor:cursor + n_null] = null_v
        cursor += n_null

    # survived/scattered
    n_surv = surv_p.shape[1]
    if n_surv:
        out_p[:, cursor:cursor + n_surv] = surv_p
        out_v[:, cursor:cursor + n_surv] = surv_v
        cursor += n_surv

    # new electrons
    n_new = new_p.shape[1]
    if n_new:
        out_p[:, cursor:cursor + n_new] = new_p
        out_v[:, cursor:cursor + n_new] = new_v


@_maybe_jit(fastmath=True, cache=True, parallel=True)
def histogram_increment_jit(values: np.ndarray, bins: np.ndarray) -> np.ndarray:
    """
    Increment histogram counts for given bins. Returns counts array len(bins)-1.
    """
    n_bins = bins.size - 1
    counts = np.zeros(n_bins, dtype=np.int64)
    for i in nb.prange(values.size):
        v = values[i]
        if v < bins[0] or v >= bins[-1]:
            continue
        idx = np.searchsorted(bins, v, side="right") - 1
        if idx >= 0 and idx < n_bins:
            counts[idx] += 1
    return counts


@_maybe_jit(fastmath=True, cache=True)
def velocity_from_energy_jit(energy: np.ndarray) -> np.ndarray:
    """
    Vectorized velocity_from_energy (m/s) using fastmath.
    """
    return np.sqrt(2.0 * energy * csts.elementary_charge / csts.electron_mass)


@_maybe_jit(fastmath=True, cache=True)
def energy_from_velocity_jit(velocity: np.ndarray) -> np.ndarray:
    """
    Vectorized energy_from_velocity (eV) using fastmath.
    """
    return _ENERGY_FACTOR * velocity * velocity


def jit_supported() -> bool:
    """
    Returns True if numba is importable.
    """
    return JIT_AVAILABLE


@_maybe_jit(cache=True)
def rng_uniform_scalar(state: np.ndarray) -> float:
    """
    Xorshift64*-like generator returning one float in [0,1).
    Mutates state in place. state must be shape (1,) uint64.
    """
    s = state[0]
    if s == 0:
        s = np.uint64(0x9e3779b97f4a7c15)
    s ^= (s << np.uint64(13)) & np.uint64(0xFFFFFFFFFFFFFFFF)
    s ^= (s >> np.uint64(7))
    s ^= (s << np.uint64(17)) & np.uint64(0xFFFFFFFFFFFFFFFF)
    state[0] = s
    return float((s >> np.uint64(11)) & np.uint64((1 << 53) - 1)) / float(1 << 53)


@_maybe_jit(cache=True)
def rng_uniform(state: np.ndarray, n: int) -> np.ndarray:
    """
    Generate n floats in [0,1) using the same xorshift64* core.
    """
    out = np.empty(n, dtype=np.float64)
    s = state[0]
    if s == 0:
        s = np.uint64(0x9e3779b97f4a7c15)
    for i in range(n):
        s ^= (s << np.uint64(13)) & np.uint64(0xFFFFFFFFFFFFFFFF)
        s ^= (s >> np.uint64(7))
        s ^= (s << np.uint64(17)) & np.uint64(0xFFFFFFFFFFFFFFFF)
        out[i] = float((s >> np.uint64(11)) & np.uint64((1 << 53) - 1)) / float(1 << 53)
    state[0] = s
    return out
