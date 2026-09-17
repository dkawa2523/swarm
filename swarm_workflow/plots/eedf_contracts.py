"""Typed contracts for model-independent EEDF comparisons."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import numpy as np


DEFAULT_REPRESENTATIVE_FIELDS_TD = (1.0, 5.0, 20.0, 100.0, 1000.0, 2500.0)
DEFAULT_TAIL_THRESHOLDS_EV = (11.55, 15.76)
MAX_SOURCE_NORMALIZATION_ERROR = 1.0e-6

COLORS = {
    "monte_carlo": "#2563EB",
    "mc_vs_propagator": "#9B3C7D",
    "two_term": "#D97706",
    "propagator": "#4D7C0F",
}


class EedfComparisonError(RuntimeError):
    """Raised when supplied EEDF evidence is incomplete or inconsistent."""


@dataclass(frozen=True, slots=True)
class EedfCase:
    solver: str
    e_over_n_td: float
    edges_eV: np.ndarray
    density_eV_inv: np.ndarray
    reported_mean_energy_eV: float
    ci95_half_density_eV_inv: np.ndarray | None = None
    uncertainty_label: str | None = None

    @property
    def widths_eV(self) -> np.ndarray:
        return np.diff(self.edges_eV)

    @property
    def centers_eV(self) -> np.ndarray:
        return 0.5 * (self.edges_eV[:-1] + self.edges_eV[1:])

    @property
    def source_normalization(self) -> float:
        return float(np.sum(self.density_eV_inv * self.widths_eV))

    @property
    def reconstructed_mean_energy_eV(self) -> float:
        mass = self.density_eV_inv * self.widths_eV
        return float(np.sum(self.centers_eV * mass) / np.sum(mass))


@dataclass(frozen=True, slots=True)
class EedfDataset:
    solver: str
    cases: Mapping[float, EedfCase]
    provenance: Mapping[str, object]


@dataclass(frozen=True, slots=True)
class EedfComparisonSummary:
    output_directory: Path
    figures: tuple[Path, ...]
    metrics_csv: Path
    manifest: Path
    mc_status: Path
    mc_status_not_passed_fields_td: tuple[float, ...]
