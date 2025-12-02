#!/usr/bin/env python3
"""
Fit generalized (and optionally shifted) EEDF models to the aggregated `eedf_table.csv` outputs.

For each unique condition (distinguished by the `mean energy (eV)` column), the script:
  - Fits a mixture of 1–N components (default N=3) of the form
        f(E) = sum_i A_i * (E - E_shift)^alpha_i * exp(-((E - E_shift) / E_c_i)^beta_i)
    evaluated on a log-scale loss with quadrature-style weighting to respect the sampled grid.
  - Selects the component count via AICc (default) to avoid overfitting tails while
    still capturing multi-hump / shifted EEDF shapes.
  - Saves comparison plots (data vs. mixture and per-component curves) with loss/coeffs annotated
    and writes a `fit_results.csv` alongside the plots.

By default all `outputs/*/eedf_table.csv` files are processed and plots land in
`<experiment>/plots/eedf_fits/`. Use `--table` to target a single file and `--save-dir`
to place results elsewhere.
"""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import matplotlib

# Use a non-interactive backend to work in headless environments.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.signal import savgol_filter
import yaml


EEDF_COLUMNS = ("mean energy (eV)", "energy (eV)", "eedf (eV-1)")

# Default configuration used when EEDF_optimize.yaml is absent or partial.
DEFAULT_CONFIG = {
    "preprocess": {
        # Keep only energies within this window (set to null to disable each).
        "energy_min": None,
        "energy_max": None,
        # Robust outlier rejection on log10(EEDF) using MAD.
        "log_mad_sigma": 4.0,
        # Floor applied before log to avoid -inf.
        "log_floor": 1e-30,
        # Dynamic floor = min_positive_eedf * factor (applied before log); useful to drop the lowest decade.
        "log_floor_factor": 10.0,
        # Optional smoothing in log space.
        "smooth": {
            "enabled": True,
            "window": 11,  # must be odd; will be clipped to data length
            "polyorder": 2,
        },
        # Resample to a uniform grid to stabilize optimization (0 to disable).
        "resample_points": 800,
        # High-energy tail handling: zero-out negligible tail area to avoid artificial plateaus.
        "tail_area_rel_threshold": 1e-4,  # fraction of total area; tail beyond this is floored
        "tail_min_energy_quantile": 0.9,  # only zero tails above this energy quantile
        # Additional trimming: drop trailing bins once EEDF falls below this fraction of its max.
        "tail_rel_keep": 1e-4,
        "tail_keep_min_quantile": 0.75,  # never trim before this quantile
        "tail_keep_min_points": 100,      # minimum points to keep before trimming
    },
    "model": {
        "max_components": 3,
        "num_starts": 6,
        "random_seed": 12345,
        # Bounds for parameters.
        "alpha_bounds": [-2.0, 4.0],
        "beta_bounds": [0.2, 4.0],
        "amplitude_min": 1e-30,
        # Characteristic energy upper bound factors.
        "e_c_max_factor": 10.0,          # times max energy
        "e_c_mean_energy_factor": 5.0,   # times mean_energy_guess
        # Shift bounds relative to max energy (fractional).
        "shift_bounds_relative": [-0.3, 0.3],
        # Optional absolute shift bounds; if set, they override relative.
        "shift_bounds_absolute": None,
    },
    "selection": {
        # Currently fixed to AICc; tie-breaker by log10 MSE.
        "criterion": "AICc",
    },
    "plot": {
        "y_min": None,
        "y_max": None,
        "y_margin_factor": 0.1,  # expand ylim by this fraction of the range
        "y_clip_percentiles": [0.1, 99.9],  # clip extreme points before setting ylim
    },
    "weights": {
        # Optional boost for low-energy region to improve fit there.
        "low_energy_boost": {
            "amplitude": 1.0,         # weight multiplier peak (0 to disable)
            "scale_rel_to_max": 0.2,  # decay scale as fraction of max energy
        },
    },
}


def _merge_dicts(base: dict, override: dict) -> dict:
    """Recursively merge two dicts (override wins)."""
    out = dict(base)
    for k, v in (override or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge_dicts(out[k], v)
        else:
            out[k] = v
    return out


def load_config(path: Path) -> dict:
    """Load YAML config if present; otherwise return defaults."""
    cfg = dict(DEFAULT_CONFIG)
    if path and path.exists():
        try:
            with path.open("r") as fh:
                user_cfg = yaml.safe_load(fh) or {}
            cfg = _merge_dicts(cfg, user_cfg)
        except Exception:
            # Fall back to defaults on parse errors.
            cfg = dict(DEFAULT_CONFIG)
    return cfg


@dataclass
class FitComponent:
    amplitude: float
    alpha: float
    characteristic_energy_eV: float
    beta: float


@dataclass
class FitResult:
    label: str
    mean_energy_eV: float
    n_components: int
    shift_eV: float
    components: list[FitComponent]
    log10_mse: float
    aicc: float
    log10_mae: float
    area_rel_error: float
    plot_path: Path


def generalized_eedf(energy: np.ndarray, amplitude: float, alpha: float, e_char: float, beta: float) -> np.ndarray:
    """Generalized EEDF model: A * E^alpha * exp(-(E/Ec)^beta)."""
    energy = np.asarray(energy, dtype=float)
    return amplitude * np.power(energy, alpha) * np.exp(-np.power(energy / e_char, beta))


def log_mse(observed: np.ndarray, predicted: np.ndarray, base: float = 10.0, weights: Optional[np.ndarray] = None) -> float:
    """
    Mean squared error in log space (base 10 by default).
    If weights are given, they must sum to 1 and the weighted mean is returned.
    """
    eps = 1e-300
    obs_safe = np.clip(observed, eps, None)
    pred_safe = np.clip(predicted, eps, None)
    if base == 10:
        obs_log = np.log10(obs_safe)
        pred_log = np.log10(pred_safe)
    else:
        obs_log = np.log(obs_safe)
        pred_log = np.log(pred_safe)
    diff2 = (pred_log - obs_log) ** 2
    if weights is not None:
        w = np.asarray(weights, dtype=float)
        w = np.clip(w, 0.0, np.inf)
        if w.sum() == 0:
            return float(np.mean(diff2))
        w = w / w.sum()
        return float(np.sum(w * diff2))
    return float(np.mean(diff2))


def compute_metrics(
    energy: np.ndarray,
    observed: np.ndarray,
    predicted: np.ndarray,
    weights: Optional[np.ndarray] = None,
) -> dict:
    """Compute multiple accuracy indicators in log space and area difference."""
    obs = np.clip(observed, 1e-300, np.inf)
    pred = np.clip(predicted, 1e-300, np.inf)
    if weights is not None:
        w = np.asarray(weights, dtype=float)
        w = np.clip(w, 0.0, np.inf)
        if w.sum() == 0:
            w = None
        else:
            w = w / w.sum()
    else:
        w = None

    log_obs = np.log10(obs)
    log_pred = np.log10(pred)
    diff = log_pred - log_obs
    if w is None:
        log10_mse = float(np.mean(diff ** 2))
        log10_mae = float(np.mean(np.abs(diff)))
    else:
        log10_mse = float(np.sum(w * diff ** 2))
        log10_mae = float(np.sum(w * np.abs(diff)))

    area_obs = float(np.trapezoid(obs, energy))
    area_pred = float(np.trapezoid(pred, energy))
    if area_obs > 0:
        area_rel_error = abs(area_pred - area_obs) / area_obs
    else:
        area_rel_error = math.inf

    return {
        "log10_mse": log10_mse,
        "log10_mae": log10_mae,
        "area_rel_error": float(area_rel_error),
    }


def spacing_weights(energy: np.ndarray) -> np.ndarray:
    """Quadrature-like weights from grid spacing to avoid over-emphasizing dense bins."""
    grad = np.abs(np.gradient(energy))
    grad = np.clip(grad, 1e-12, None)
    grad /= np.sum(grad)
    return grad


def build_weights(energy: np.ndarray, base: np.ndarray, weight_cfg: dict) -> np.ndarray:
    """Compose base weights with optional low-energy boost."""
    w = np.array(base, dtype=float)
    boost_cfg = (weight_cfg or {}).get("low_energy_boost", {}) or {}
    amp = float(boost_cfg.get("amplitude", 0.0) or 0.0)
    if amp > 0 and energy.size:
        scale_rel = float(boost_cfg.get("scale_rel_to_max", 0.0) or 0.0)
        if scale_rel > 0:
            scale = max(scale_rel * float(np.max(energy)), 1e-9)
            boost = 1.0 + amp * np.exp(-energy / scale)
            w = w * boost
    total = np.sum(w)
    if total > 0:
        w = w / total
    return w


def preprocess_eedf(energy: np.ndarray, eedf: np.ndarray, cfg: dict) -> tuple[np.ndarray, np.ndarray]:
    """
    Preprocess raw EEDF samples:
      - apply energy windowing
      - floor + log-domain outlier rejection via MAD
      - optional smoothing (Savitzky-Golay in log-space)
      - optional uniform resampling for stable optimization
      - zero-out negligible high-energy tail (avoid artificial plateaus from single particles)
    """
    energy = np.asarray(energy, dtype=float)
    eedf = np.asarray(eedf, dtype=float)

    mask = (energy > 0) & (eedf > 0)
    emin = cfg.get("energy_min")
    emax = cfg.get("energy_max")
    if emin is not None:
        mask &= energy >= float(emin)
    if emax is not None:
        mask &= energy <= float(emax)

    energy = energy[mask]
    eedf = eedf[mask]
    if energy.size == 0:
        return energy, eedf

    order = np.argsort(energy)
    energy = energy[order]
    eedf = eedf[order]

    floor_base = max(cfg.get("log_floor", 1e-30), 1e-300)
    floor_factor = cfg.get("log_floor_factor", None)
    min_pos = float(np.min(eedf[eedf > 0])) if np.any(eedf > 0) else 0.0
    floor_dyn = min_pos * float(floor_factor) if floor_factor is not None and floor_factor > 0 and min_pos > 0 else 0.0
    floor = max(floor_base, floor_dyn)
    # Guard against degenerate case where floor exceeds most data
    max_val = float(np.max(eedf)) if eedf.size else floor
    if floor >= max_val:
        floor = min(floor_base, max_val * 0.5) if max_val > 0 else floor_base
    keep_floor = eedf >= floor
    energy = energy[keep_floor]
    eedf = eedf[keep_floor]
    if energy.size == 0:
        return energy, eedf
    eedf = np.clip(eedf, floor, None)
    logy = np.log10(eedf)

    sigma = cfg.get("log_mad_sigma", None)
    if sigma is not None and sigma > 0 and np.isfinite(sigma):
        med = np.median(logy)
        mad = np.median(np.abs(logy - med))
        if mad > 0:
            z = 0.6744897501960817 * (logy - med) / mad
            keep = np.abs(z) <= sigma
            energy = energy[keep]
            logy = logy[keep]
            if energy.size == 0:
                return energy, np.array([], dtype=float)

    smooth_cfg = cfg.get("smooth", {}) or {}
    if smooth_cfg.get("enabled", False) and logy.size >= 5:
        win = int(smooth_cfg.get("window", 11))
        if win % 2 == 0:
            win += 1
        win = max(5, min(win, logy.size - (1 - logy.size % 2)))
        poly = int(min(smooth_cfg.get("polyorder", 2), max(1, win - 2)))
        try:
            logy = savgol_filter(logy, window_length=win, polyorder=poly, mode="interp")
        except Exception:
            pass  # fallback to unsmoothed

    resample_points = int(cfg.get("resample_points", 0) or 0)
    if resample_points > 0 and energy.size >= 2:
        e_new = np.linspace(energy.min(), energy.max(), resample_points)
        log_interp = np.interp(e_new, energy, logy)
        energy = e_new
        logy = log_interp

    # High-energy tail zeroing (treat negligible area as zero instead of a flat plateau).
    tail_area_rel = cfg.get("tail_area_rel_threshold", None)
    tail_q = cfg.get("tail_min_energy_quantile", 0.9)
    if tail_area_rel is not None and tail_area_rel > 0 and energy.size >= 3:
        energy_q = np.quantile(energy, min(max(tail_q, 0.0), 1.0))
        # Compute cumulative area from high energy side.
        area_total = float(np.trapezoid(np.power(10.0, logy), energy))
        if area_total > 0:
            dE = np.diff(energy)
            y_lin = np.power(10.0, logy)
            area_segments = 0.5 * (y_lin[:-1] + y_lin[1:]) * dE
            cum_rev = np.cumsum(area_segments[::-1])[::-1]
            # Align length to energy[:-1]; add trailing zero for last point.
            remaining = np.concatenate([cum_rev, [0.0]])
            rel_remaining = remaining / area_total
            mask_tail = (rel_remaining <= tail_area_rel) & (energy >= energy_q)
            if np.any(mask_tail):
                logy = np.where(mask_tail, math.log10(floor), logy)

    # Drop trailing bins after EEDF falls below relative threshold (avoid single-particle plateaus).
    tail_rel_keep = cfg.get("tail_rel_keep", None)
    if tail_rel_keep is not None and tail_rel_keep > 0 and energy.size >= 10:
        y_lin = np.power(10.0, logy)
        max_val = float(np.max(y_lin))
        thresh = max_val * float(tail_rel_keep)
        idx_keep = np.where(y_lin >= thresh)[0]
        if idx_keep.size:
            last_sig = int(idx_keep.max())
            min_q = float(cfg.get("tail_keep_min_quantile", 0.75))
            min_q_idx = int(min_q * (energy.size - 1))
            min_pts = int(cfg.get("tail_keep_min_points", 0) or 0)
            min_idx = max(last_sig, min_q_idx, min_pts)
            cut = min(energy.size, max(min_idx + 1, 1))
            energy = energy[:cut]
            logy = logy[:cut]

    eedf = np.power(10.0, logy)
    eedf = np.clip(eedf, floor, None)
    return energy, eedf


def mixture_eedf(energy: np.ndarray, shift_eV: float, params: np.ndarray, n_components: int) -> np.ndarray:
    """Evaluate mixture of n_components using packed params [A, alpha, E_c, beta]*k after shift."""
    shifted = np.clip(energy - shift_eV, 1e-6, None)
    pred = np.zeros_like(shifted)
    for i in range(n_components):
        a, alpha, e_c, beta = params[i * 4 : (i + 1) * 4]
        pred += generalized_eedf(shifted, a, alpha, e_c, beta)
    return pred


def aicc_from_rss(rss: float, n_samples: int, n_params: int) -> float:
    """Compute corrected AIC; return +inf when undefined (n too small)."""
    if n_samples <= n_params + 1 or rss <= 0:
        return math.inf
    return (
        2 * n_params
        + n_samples * math.log(rss / n_samples)
        + (2 * n_params * (n_params + 1)) / (n_samples - n_params - 1)
    )


def decode_components(params: np.ndarray, n_components: int) -> list[FitComponent]:
    comps: list[FitComponent] = []
    for i in range(n_components):
        a, alpha, e_c, beta = params[i * 4 : (i + 1) * 4]
        comps.append(FitComponent(float(a), float(alpha), float(e_c), float(beta)))
    return comps


def initial_guesses(
    energy: np.ndarray,
    eedf: np.ndarray,
    k: int,
    shift_bounds: tuple[float, float],
    weighted_mean_energy: float,
    num_starts: int,
    rng: np.random.Generator,
    amplitude_floor: float,
) -> list[np.ndarray]:
    """Build a set of initial guesses (with jitter) for robustness."""
    shift0 = 0.0
    quantiles = np.linspace(0.25, 0.9, k)
    e_chars = np.quantile(energy, quantiles)
    base_amp = max(np.max(eedf) / max(k, 1), amplitude_floor)
    guesses: list[np.ndarray] = []

    base = [shift0]
    for i in range(k):
        base.extend([base_amp, 0.5, max(e_chars[i], weighted_mean_energy / 2), 1.0])
    guesses.append(np.array(base, dtype=float))

    for _ in range(max(1, num_starts - 1)):
        g = [float(np.clip(rng.normal(shift0, 0.05 * energy.max()), *shift_bounds))]
        for i in range(k):
            amp_jitter = max(base_amp * rng.uniform(0.5, 1.6), amplitude_floor)
            alpha_jitter = rng.normal(0.5, 0.2)
            e_c_jitter = np.clip(e_chars[i] * rng.uniform(0.6, 1.6), 1e-5, energy.max() * 10)
            beta_jitter = np.clip(rng.normal(1.0, 0.2), 0.2, 4.0)
            g.extend([amp_jitter, alpha_jitter, e_c_jitter, beta_jitter])
        guesses.append(np.array(g, dtype=float))

    return guesses


def fit_mixture_model(
    energy: np.ndarray,
    eedf: np.ndarray,
    mean_energy_guess: float,
    model_cfg: dict,
    weight_cfg: dict,
) -> tuple[np.ndarray, dict, float, int]:
    """
    Fit 1..max_components mixtures with a global shift. Returns best params.

    Output:
        params: [shift, A1, alpha1, E_c1, beta1, A2, alpha2, ...]
        best_metrics: dict with log10_mse/log10_mae/area_rel_error
        best_aicc: AICc of best model
        best_k: number of components selected
    """
    mask = (energy > 0) & (eedf > 0)
    energy = energy[mask]
    eedf = eedf[mask]
    if energy.size == 0:
        raise ValueError("No positive EEDF samples to fit.")

    order = np.argsort(energy)
    energy = energy[order]
    eedf = eedf[order]

    weights = spacing_weights(energy)
    weights = build_weights(energy, weights, weight_cfg)
    weighted_mean_energy = np.average(energy, weights=eedf)
    energy_max = float(energy.max())

    max_components = int(model_cfg.get("max_components", 3))
    num_starts = int(model_cfg.get("num_starts", 6))
    rng = np.random.default_rng(model_cfg.get("random_seed", 12345))
    alpha_bounds = model_cfg.get("alpha_bounds", [-2.0, 4.0])
    beta_bounds = model_cfg.get("beta_bounds", [0.2, 4.0])
    amplitude_min = float(model_cfg.get("amplitude_min", 1e-30))
    e_c_max_factor = float(model_cfg.get("e_c_max_factor", 10.0))
    e_c_mean_energy_factor = float(model_cfg.get("e_c_mean_energy_factor", 5.0))
    shift_abs = model_cfg.get("shift_bounds_absolute")
    if shift_abs is not None and isinstance(shift_abs, (list, tuple)) and len(shift_abs) == 2:
        shift_bounds = (float(shift_abs[0]), float(shift_abs[1]))
    else:
        shift_rel = model_cfg.get("shift_bounds_relative", [-0.3, 0.3])
        shift_bounds = (
            float(shift_rel[0]) * energy_max,
            float(shift_rel[1]) * energy_max,
        )

    best_params: Optional[np.ndarray] = None
    best_aicc = math.inf
    best_metrics: dict = {}
    best_k = 1

    for k in range(1, max_components + 1):
        lb = [shift_bounds[0]]
        ub = [shift_bounds[1]]
        for _ in range(k):
            lb.extend([amplitude_min, alpha_bounds[0], 1e-6, beta_bounds[0]])
            ub.extend([np.inf, alpha_bounds[1], max(energy_max * e_c_max_factor, mean_energy_guess * e_c_mean_energy_factor), beta_bounds[1]])
        bounds = (np.array(lb), np.array(ub))

        for guess in initial_guesses(
            energy,
            eedf,
            k,
            shift_bounds,
            weighted_mean_energy,
            num_starts,
            rng,
            amplitude_min,
        ):
            def residual(params: np.ndarray) -> np.ndarray:
                shift = params[0]
                comps = params[1:]
                model = mixture_eedf(energy, shift, comps, k)
                obs = np.clip(eedf, 1e-300, np.inf)
                pred = np.clip(model, 1e-300, np.inf)
                return np.sqrt(weights) * (np.log(pred) - np.log(obs))

            try:
                result = least_squares(
                    residual,
                    x0=guess,
                    bounds=bounds,
                    method="trf",
                    ftol=1e-12,
                    xtol=1e-12,
                    gtol=1e-12,
                    max_nfev=30000,
                )
            except Exception:
                continue

            rss = float(np.sum(result.fun ** 2))
            n_samples = result.fun.size
            n_params = result.x.size
            aicc = aicc_from_rss(rss, n_samples, n_params)

            model_pred = np.clip(mixture_eedf(energy, result.x[0], result.x[1:], k), 1e-300, np.inf)
            metrics = compute_metrics(energy, eedf, model_pred, weights=weights)

            # Prefer lower AICc; tie-break with log10 MSE.
            if (aicc < best_aicc - 1e-9) or (
                math.isclose(aicc, best_aicc) and metrics["log10_mse"] < best_metrics.get("log10_mse", math.inf)
            ):
                best_params = result.x
                best_aicc = aicc
                best_metrics = metrics
                best_k = k

    if best_params is None:
        raise RuntimeError("Fitting failed for all component counts.")

    return best_params, best_metrics, best_aicc, best_k


def slugify(text: str) -> str:
    """Create a filesystem-friendly slug."""
    return re.sub(r"[^0-9A-Za-z._-]+", "_", text).strip("_")


def find_tables(root: Path, explicit_table: Optional[Path]) -> list[Path]:
    if explicit_table:
        return [explicit_table.resolve()]
    return sorted(root.glob("*/eedf_table.csv"))


def build_label_lookup(summary_path: Path) -> Optional[pd.DataFrame]:
    if not summary_path.exists():
        return None
    try:
        df = pd.read_csv(summary_path)
    except Exception:
        return None
    needed = {"mean energy (eV)", "run_label"}
    if not needed.issubset(df.columns):
        return None
    return df


def label_for_mean_energy(mean_energy: float, summary_df: Optional[pd.DataFrame]) -> str:
    if summary_df is None:
        return f"meanE_{mean_energy:.3f}eV"
    diffs = np.abs(summary_df["mean energy (eV)"] - mean_energy)
    idx = int(diffs.idxmin())
    if diffs.iloc[idx] > max(1e-5, 1e-6 * max(mean_energy, 1.0)):
        return f"meanE_{mean_energy:.3f}eV"

    row = summary_df.loc[idx]
    label = str(row.get("run_label", f"meanE_{mean_energy:.3f}eV"))
    sweep_val = row.get("sweep_value")
    try:
        if not (isinstance(sweep_val, float) and math.isnan(sweep_val)) and sweep_val is not None:
            label = f"{label}_sweep{float(sweep_val):g}"
    except Exception:
        pass
    return label


def plot_fit(
    energy: np.ndarray,
    eedf: np.ndarray,
    params: np.ndarray,
    metrics: dict,
    aicc: float,
    n_components: int,
    title: str,
    out_path: Path,
    plot_cfg: dict,
) -> None:
    energies_sorted = np.linspace(np.min(energy), np.max(energy), 600)
    shift = params[0]
    comps = params[1:]
    fitted = mixture_eedf(energies_sorted, shift, comps, n_components)
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(energy, eedf, label="EEDF (data)", color="C0", linewidth=1.2)
    ax.plot(energies_sorted, fitted, label=f"Mixture fit (k={n_components})", color="C1", linewidth=1.4)

    for i in range(n_components):
        comp = generalized_eedf(np.clip(energies_sorted - shift, 1e-6, None), *comps[i * 4 : (i + 1) * 4])
        ax.plot(
            energies_sorted,
            comp,
            linestyle="--",
            linewidth=1.0,
            label=f"Component {i+1}",
        )

    ax.set_yscale("log")
    # Y-axis tightening: use percentiles of data+fit and an expansion margin.
    y_all = [eedf[eedf > 0], fitted[fitted > 0]]
    for i in range(n_components):
        comp_vals = generalized_eedf(np.clip(energies_sorted - shift, 1e-6, None), *comps[i * 4 : (i + 1) * 4])
        y_all.append(comp_vals[comp_vals > 0])
    y_all = np.concatenate([arr.ravel() for arr in y_all if arr.size])
    if y_all.size:
        clip_p = plot_cfg.get("y_clip_percentiles", [0.1, 99.9]) or [0.1, 99.9]
        try:
            y_lo, y_hi = np.percentile(y_all, clip_p)
        except Exception:
            y_lo, y_hi = y_all.min(), y_all.max()
        margin = plot_cfg.get("y_margin_factor", 0.1) or 0.0
        y_min_auto = max(y_lo * (1 - margin), 1e-12)
        y_max_auto = y_hi * (1 + margin)
        y_min = plot_cfg.get("y_min", None)
        y_max = plot_cfg.get("y_max", None)
        ax.set_ylim(y_min if y_min is not None else y_min_auto, y_max if y_max is not None else y_max_auto)

    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("EEDF (eV$^{-1}$)")
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.3)
    annotation_lines = [
        f"shift={shift:.3e} eV",
        f"log10-MSE={metrics.get('log10_mse', math.nan):.3e}",
        f"log10-MAE={metrics.get('log10_mae', math.nan):.3e}",
        f"area rel err={metrics.get('area_rel_error', math.nan):.3e}",
        f"AICc={aicc:.3f}",
    ]
    for i in range(n_components):
        a, alpha, e_c, beta = comps[i * 4 : (i + 1) * 4]
        annotation_lines.append(f"[{i+1}] A={a:.2e}, a={alpha:.2f}, Ec={e_c:.2f}, b={beta:.2f}")
    annotation = "\n".join(annotation_lines)
    ax.text(
        0.98,
        0.97,
        annotation,
        transform=ax.transAxes,
        fontsize=8,
        ha="right",
        va="top",
        bbox=dict(facecolor="white", alpha=0.8, edgecolor="0.7", boxstyle="round,pad=0.35"),
    )
    ax.legend()
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def process_table(
    table_path: Path,
    save_root: Optional[Path] = None,
    config: Optional[dict] = None,
) -> list[FitResult]:
    df = pd.read_csv(table_path)
    if not set(EEDF_COLUMNS).issubset(df.columns):
        raise ValueError(f"{table_path} is missing required columns {EEDF_COLUMNS}")

    cfg = config or DEFAULT_CONFIG
    summary_df = build_label_lookup(table_path.parent / "summary.csv")
    dest_dir = (save_root / table_path.parent.name) if save_root else table_path.parent / "plots" / "eedf_fits"
    results: list[FitResult] = []

    mean_energies = sorted(df["mean energy (eV)"].unique())
    for mean_energy in mean_energies:
        mask = np.isclose(df["mean energy (eV)"], mean_energy, rtol=0, atol=1e-12)
        subset = df.loc[mask]
        energy = subset["energy (eV)"].to_numpy(dtype=float)
        eedf = subset["eedf (eV-1)"].to_numpy(dtype=float)

        energy_pp, eedf_pp = preprocess_eedf(energy, eedf, cfg["preprocess"])
        if energy_pp.size == 0:
            continue

        params, metrics, aicc, k = fit_mixture_model(
            energy_pp,
            eedf_pp,
            mean_energy_guess=float(mean_energy),
            model_cfg=cfg["model"],
            weight_cfg=cfg.get("weights", {}),
        )
        run_label = label_for_mean_energy(float(mean_energy), summary_df)
        file_label = slugify(run_label)
        out_path = dest_dir / f"eedf_fit_{file_label}.png"
        title = f"{table_path.parent.name} | {run_label}"
        plot_fit(energy_pp, eedf_pp, params, metrics, aicc, k, title, out_path, cfg["plot"])

        results.append(
            FitResult(
                label=run_label,
                mean_energy_eV=float(mean_energy),
                n_components=k,
                shift_eV=float(params[0]),
                components=decode_components(params[1:], k),
                log10_mse=float(metrics.get("log10_mse", math.nan)),
                log10_mae=float(metrics.get("log10_mae", math.nan)),
                area_rel_error=float(metrics.get("area_rel_error", math.nan)),
                aicc=float(aicc),
                plot_path=out_path,
            )
        )

    if results:
        dest_dir.mkdir(parents=True, exist_ok=True)
        max_comps = int(cfg["model"].get("max_components", max(r.n_components for r in results)))
        results_df = pd.DataFrame(
            [
                {
                    "label": r.label,
                    "mean_energy (eV)": r.mean_energy_eV,
                    "n_components": r.n_components,
                    "shift (eV)": r.shift_eV,
                    "log10_mse": r.log10_mse,
                    "log10_mae": r.log10_mae,
                    "area_rel_error": r.area_rel_error,
                    "AICc": r.aicc,
                    "plot_path": str(r.plot_path),
                    "source_table": str(table_path),
                    **{
                        f"component_{i+1}_amplitude": r.components[i].amplitude if i < len(r.components) else None
                        for i in range(max_comps)
                    },
                    **{
                        f"component_{i+1}_alpha": r.components[i].alpha if i < len(r.components) else None
                        for i in range(max_comps)
                    },
                    **{
                        f"component_{i+1}_E_c (eV)": r.components[i].characteristic_energy_eV if i < len(r.components) else None
                        for i in range(max_comps)
                    },
                    **{
                        f"component_{i+1}_beta": r.components[i].beta if i < len(r.components) else None
                        for i in range(max_comps)
                    },
                }
                for r in results
            ]
        )
        results_df.to_csv(dest_dir / "fit_results.csv", index=False)

    return results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fit generalized EEDF models to eedf_table.csv outputs.")
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("outputs"),
        help="Root directory containing experiment folders with eedf_table.csv (ignored if --table is given).",
    )
    parser.add_argument(
        "--table",
        type=Path,
        help="Path to a specific eedf_table.csv to process. When set, only this table is used.",
    )
    parser.add_argument(
        "--save-dir",
        type=Path,
        help="Optional directory to save plots/results. A subfolder named after each experiment is created inside.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("EEDF_optimize.yaml"),
        help="Path to YAML config for preprocessing/fitting. Defaults to EEDF_optimize.yaml.",
    )
    parser.add_argument(
        "--max-components",
        type=int,
        default=None,
        choices=range(1, 5),
        metavar="{1-4}",
        help="Override: maximum number of mixture components (otherwise use YAML).",
    )
    parser.add_argument(
        "--num-starts",
        type=int,
        default=None,
        help="Override: number of initial guesses per component count (otherwise use YAML).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cfg = load_config(args.config)
    if args.max_components is not None:
        cfg.setdefault("model", {})["max_components"] = args.max_components
    if args.num_starts is not None:
        cfg.setdefault("model", {})["num_starts"] = args.num_starts
    tables = find_tables(args.root, args.table)
    if not tables:
        raise SystemExit(f"No eedf_table.csv found under {args.root.resolve()}")

    for table_path in tables:
        print(f"[INFO] Processing {table_path}")
        results = process_table(
            table_path,
            save_root=args.save_dir,
            config=cfg,
        )
        for res in results:
            print(
                f"  - {res.label}: k={res.n_components}, shift={res.shift_eV:.3e} eV, "
                f"log10-MSE={res.log10_mse:.3e}, log10-MAE={res.log10_mae:.3e}, "
                f"area_err={res.area_rel_error:.3e}, AICc={res.aicc:.3f} -> {res.plot_path}"
            )


if __name__ == "__main__":
    main()
