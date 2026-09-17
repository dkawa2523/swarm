"""Contracts for saved-solution GEC-ICP comparison analysis."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping


MAX_COMPARISON_CASES = 4
MIN_COMPARISON_CASES = 2

VOLUME_COLUMNS = (
    "Time (s)",
    "axisymmetric plasma volume (m^3)",
    "electron inventory (1)",
    "electron mean-energy inventory (eV)",
    "argon metastable inventory (1)",
    "argon ion inventory (1)",
    "plasma RF resistive deposition (W)",
)
COIL_COLUMNS = ("Time (s)", "configured coil-power readback (W)")
SOLUTION_TIME_COLUMNS = ("solution_index", "time_s")

METRIC_KEYS = (
    "electron_inventory",
    "electron_weighted_mean_energy_eV",
    "metastable_inventory",
    "ion_inventory",
    "absorbed_power_W",
    "coil_power_W",
)
METRIC_LABELS = {
    "electron_inventory": "electron inventory",
    "electron_weighted_mean_energy_eV": "mean electron energy (eV)",
    "metastable_inventory": "Ar* inventory",
    "ion_inventory": "Ar+ inventory",
    "absorbed_power_W": "absorbed RF power (W)",
    "coil_power_W": "coil power (W)",
}
DEFAULT_STATIONARITY_LIMITS = {
    "electron_inventory": 0.02,
    "electron_weighted_mean_energy_eV": 0.01,
    "metastable_inventory": 0.03,
    "ion_inventory": 0.02,
    "absorbed_power_W": 0.01,
    "coil_power_W": 0.005,
}


class GecIcpComparisonAnalysisError(RuntimeError):
    """Raised when comparison evidence is incomplete or inconsistent."""


@dataclass(frozen=True, slots=True)
class GecIcpIntegratedState:
    """One integrated plasma state at an explicitly identified time."""

    time_s: float
    axisymmetric_volume_m3: float
    electron_inventory: float
    electron_weighted_mean_energy_eV: float
    metastable_inventory: float
    ion_inventory: float
    absorbed_power_W: float
    coil_power_W: float

    def metrics(self) -> dict[str, float]:
        return {key: float(getattr(self, key)) for key in METRIC_KEYS}


@dataclass(frozen=True, slots=True)
class GecIcpComparisonCaseData:
    """Validated evidence loaded from one read-only comparison export."""

    case_id: str
    label: str
    export_directory: Path
    input_mph: Path
    input_mph_sha256: str
    source_manifest: Path
    source_manifest_sha256: str
    model_mapping: Path
    model_mapping_sha256: str
    common: GecIcpIntegratedState
    terminal: GecIcpIntegratedState
    time_series: tuple[GecIcpIntegratedState, ...]
    solution_times_s: tuple[float, ...]
    source_hashes: Mapping[str, str]
    stationarity_changes: Mapping[str, float]
    stationarity_limits: Mapping[str, float]
    stationarity_passed: Mapping[str, bool]

    @property
    def converged(self) -> bool:
        return all(self.stationarity_passed.values())


@dataclass(frozen=True, slots=True)
class GecIcpComparisonAnalysisSummary:
    """Paths and case state emitted by a completed comparison analysis."""

    output_directory: Path
    comparison_csv: Path
    stationarity_csv: Path
    figures: tuple[Path, ...]
    manifest: Path
    cases: tuple[GecIcpComparisonCaseData, ...]


__all__ = (
    "COIL_COLUMNS",
    "DEFAULT_STATIONARITY_LIMITS",
    "GecIcpComparisonAnalysisError",
    "GecIcpComparisonAnalysisSummary",
    "GecIcpComparisonCaseData",
    "GecIcpIntegratedState",
    "MAX_COMPARISON_CASES",
    "METRIC_KEYS",
    "METRIC_LABELS",
    "MIN_COMPARISON_CASES",
    "SOLUTION_TIME_COLUMNS",
    "VOLUME_COLUMNS",
)
