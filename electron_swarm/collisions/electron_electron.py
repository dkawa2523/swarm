"""Lightweight electron-electron collision helpers.

The first mergeable step is intentionally deterministic and solver-neutral: it
provides relaxation targets and metadata helpers that two-term or multi-term
operators can call later.  It does not introduce MC Coulomb-collision sampling or
new public result contracts.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


_EPS = 1.0e-300


@dataclass(frozen=True, slots=True)
class ElectronElectronSettings:
    enabled: bool = False
    model: str = "relaxation"
    electron_density_m3: float | None = None
    coulomb_log: float | str = "auto"
    conserve_mean_energy: bool = True
    strength_scale: float = 1.0


def coulomb_log_auto(
    electron_density_m3: float,
    electron_temperature_eV: float,
) -> float:
    """Return a bounded engineering estimate for the Coulomb logarithm."""

    ne = max(float(electron_density_m3), 1.0)
    te = max(float(electron_temperature_eV), 1.0e-3)
    estimate = 23.0 - 0.5 * np.log(ne * 1.0e-6) + 1.5 * np.log(te)
    return float(np.clip(estimate, 2.0, 20.0))


def _normalize(values: np.ndarray, widths: np.ndarray) -> np.ndarray:
    total = float(np.sum(values * widths))
    if total <= 0.0 or not np.isfinite(total):
        raise ValueError("Cannot normalize non-positive EEDF target")
    return values / total


def mean_energy_eV(energy_eV: np.ndarray, widths_eV: np.ndarray, eedf_eV_inv: np.ndarray) -> float:
    energy = np.asarray(energy_eV, dtype=float)
    widths = np.asarray(widths_eV, dtype=float)
    eedf = np.asarray(eedf_eV_inv, dtype=float)
    norm = max(float(np.sum(eedf * widths)), _EPS)
    return float(np.sum(energy * eedf * widths) / norm)


def maxwellian_eedf_like(
    energy_eV: np.ndarray,
    widths_eV: np.ndarray,
    mean_energy_target_eV: float,
) -> np.ndarray:
    """Return a normalized Maxwellian-like energy PDF matching mean energy.

    For a Maxwellian energy distribution, <eps> = 3/2 kT, so kT = 2/3 <eps>.
    """

    energy = np.asarray(energy_eV, dtype=float)
    widths = np.asarray(widths_eV, dtype=float)
    kT = max(2.0 * float(mean_energy_target_eV) / 3.0, 1.0e-8)
    target = np.sqrt(np.maximum(energy, 0.0)) * np.exp(-np.maximum(energy, 0.0) / kT)
    return _normalize(np.clip(target, 0.0, None), widths)


def relaxation_target(
    energy_eV: np.ndarray,
    widths_eV: np.ndarray,
    eedf_eV_inv: np.ndarray,
    *,
    conserve_mean_energy: bool = True,
    fallback_temperature_eV: float = 2.0,
) -> np.ndarray:
    """Return the EEDF target for a relaxation-model e-e operator."""

    if conserve_mean_energy:
        target_mean = mean_energy_eV(energy_eV, widths_eV, eedf_eV_inv)
    else:
        target_mean = 1.5 * float(fallback_temperature_eV)
    return maxwellian_eedf_like(energy_eV, widths_eV, target_mean)


def metadata_from_settings(settings: ElectronElectronSettings) -> dict[str, bool | float | str]:
    """Flatten settings into existing ``SwarmCaseResult.metadata`` style."""

    data: dict[str, bool | float | str] = {
        "ee_enabled": bool(settings.enabled),
        "ee_model": str(settings.model),
        "ee_conserve_mean_energy": bool(settings.conserve_mean_energy),
        "ee_strength_scale": float(settings.strength_scale),
    }
    if settings.electron_density_m3 is not None:
        data["ee_electron_density_m3"] = float(settings.electron_density_m3)
    if isinstance(settings.coulomb_log, (int, float)):
        data["ee_coulomb_log"] = float(settings.coulomb_log)
    else:
        data["ee_coulomb_log"] = str(settings.coulomb_log)
    return data
