"""Shared value types and formulas for Monte Carlo transport estimators."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, slots=True)
class DirectFluxTransport:
    """Shared value contract for direct particle-flux transport estimates.

    Concrete fixed-population and weighted-growth owners define their distinct
    sampling formulas. Energy diffusion remains a restricted one-sided
    density-packet response, not a complete two-gradient response matrix.
    """

    drift_velocity_m_s: float
    mobility_m2_V_s: float
    diffusion_L_m2_s: float
    diffusion_T_m2_s: float
    energy_mobility_m2_V_s: float
    energy_diffusion_L_m2_s: float
    energy_diffusion_T_m2_s: float
    mean_energy_eV: float


def _weighted_mean(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    shape = (-1,) + (1,) * (values.ndim - 1)
    return np.sum(values * weights.reshape(shape), axis=0)
