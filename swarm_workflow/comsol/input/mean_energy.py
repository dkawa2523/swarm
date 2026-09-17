"""Mean-energy support handling shared by COMSOL input consumers."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True, slots=True)
class MeanEnergyArgument:
    """An unmodified mean-energy argument with finite qualified support."""

    minimum_eV: float
    maximum_eV: float

    def __post_init__(self) -> None:
        if not (
            math.isfinite(self.minimum_eV)
            and math.isfinite(self.maximum_eV)
            and 0.0 < self.minimum_eV < self.maximum_eV
        ):
            raise ValueError(
                "mean-energy argument requires a positive increasing finite range"
            )

    def log_values(self, mean_energy_eV: np.ndarray) -> np.ndarray:
        energy = np.asarray(mean_energy_eV, dtype=float)
        if np.any(~np.isfinite(energy)) or np.any(energy <= 0.0):
            raise ValueError("mean-energy arguments must be positive and finite")
        if np.any(energy < self.minimum_eV) or np.any(energy > self.maximum_eV):
            requested = [float(np.min(energy)), float(np.max(energy))]
            raise ValueError(
                "mean energy is outside qualified support; new solver anchors "
                f"are required for {requested} eV"
            )
        return np.log(energy)

    def values(self, mean_energy_eV: np.ndarray) -> np.ndarray:
        return np.exp(self.log_values(mean_energy_eV))

    def expression(self, log_mean_energy: str) -> str:
        """Use the physical argument directly; support is an acceptance gate."""

        return f"({log_mean_energy})"
