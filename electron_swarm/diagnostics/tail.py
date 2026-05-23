"""Product tail metrics for high-threshold reaction-rate quality."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.constants import ELECTRON_MASS_KG, EV_TO_J
from electron_swarm.core.cross_sections import (
    CrossSectionProcess,
    CrossSectionSet,
    ProcessType,
)
from electron_swarm.core.results import RateResult, SwarmCaseResult

from .common import widths_from_centers

REACTION_TYPES = {
    ProcessType.EXCITATION.value,
    ProcessType.IONIZATION.value,
    ProcessType.ATTACHMENT.value,
    ProcessType.SUPERELASTIC.value,
}
DEFAULT_TAIL_THRESHOLD_EV = 20.0
_EPS = 1.0e-300


@dataclass(frozen=True, slots=True)
class TailMetrics:
    tail_threshold_eV: float
    tail_probability: float
    tail_rate_fraction_max: float
    dominant_tail_process: str
    high_energy_cutoff_rate_fraction: float
    energy_grid_tail_status: str

    def as_metadata(self) -> dict[str, float | str | bool]:
        return {
            "tail_metrics_enabled": True,
            "tail_threshold_eV": self.tail_threshold_eV,
            "tail_probability": self.tail_probability,
            "tail_rate_fraction_max": self.tail_rate_fraction_max,
            "dominant_tail_process": self.dominant_tail_process,
            "high_energy_cutoff_rate_fraction": self.high_energy_cutoff_rate_fraction,
            "energy_grid_tail_status": self.energy_grid_tail_status,
        }


def resolve_tail_threshold_eV(
    config: SwarmConfig, cross_sections: CrossSectionSet
) -> float:
    explicit = config.physics.energy_grid_policy.tail_threshold_eV
    if explicit is not None:
        return float(explicit)
    thresholds = [
        float(proc.threshold_eV)
        for proc in cross_sections.processes
        if proc.process_type.value in REACTION_TYPES
        and proc.threshold_eV is not None
        and np.isfinite(proc.threshold_eV)
        and proc.threshold_eV > 0.0
    ]
    return max(thresholds) if thresholds else DEFAULT_TAIL_THRESHOLD_EV


def tail_probability(
    energy_eV: np.ndarray,
    eedf_eV_inv: np.ndarray,
    threshold_eV: float,
    widths_eV: np.ndarray | None = None,
) -> float:
    energy = np.asarray(energy_eV, dtype=float)
    eedf = np.asarray(eedf_eV_inv, dtype=float)
    widths = widths_from_centers(energy) if widths_eV is None else np.asarray(widths_eV, dtype=float)
    if len(energy) == 0 or len(eedf) != len(energy) or len(widths) != len(energy):
        return float("nan")
    mask = energy >= float(threshold_eV)
    if not np.any(mask):
        return 0.0
    return float(np.sum(np.clip(eedf[mask], 0.0, None) * widths[mask]))


def rate_tail_fraction(
    energy_eV: np.ndarray,
    eedf_eV_inv: np.ndarray,
    process: CrossSectionProcess,
    threshold_eV: float,
    widths_eV: np.ndarray | None = None,
) -> float:
    energy = np.asarray(energy_eV, dtype=float)
    eedf = np.asarray(eedf_eV_inv, dtype=float)
    widths = widths_from_centers(energy) if widths_eV is None else np.asarray(widths_eV, dtype=float)
    if len(energy) == 0 or len(eedf) != len(energy) or len(widths) != len(energy):
        return float("nan")
    speed = np.sqrt(np.maximum(2.0 * energy * EV_TO_J / ELECTRON_MASS_KG, 0.0))
    integrand = process.sigma(energy) * speed * np.clip(eedf, 0.0, None) * widths
    total = float(np.sum(integrand))
    if total <= _EPS or not np.isfinite(total):
        return 0.0
    tail = float(np.sum(integrand[energy >= float(threshold_eV)]))
    return float(np.clip(tail / total, 0.0, 1.0))


def _matching_process(
    rate: RateResult, cross_sections: CrossSectionSet
) -> CrossSectionProcess | None:
    for proc in cross_sections.processes:
        if (
            proc.species == rate.species
            and proc.process == rate.process
            and proc.process_type.value == rate.process_type
        ):
            return proc
    return None


def _status(max_tail: float, max_cutoff: float, warning_fraction: float) -> str:
    if max_cutoff > warning_fraction:
        return "insufficient"
    if max_tail > warning_fraction:
        return "warning"
    return "ok"


def attach_tail_metrics(
    case: SwarmCaseResult,
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
) -> SwarmCaseResult:
    policy = config.physics.energy_grid_policy
    if not policy.tail_metrics:
        case.metadata.update(
            {
                "tail_metrics_enabled": False,
                "tail_threshold_eV": None,
                "tail_probability": None,
                "tail_rate_fraction_max": None,
                "dominant_tail_process": "",
                "high_energy_cutoff_rate_fraction": None,
                "energy_grid_tail_status": "disabled",
            }
        )
        for rate in case.rates:
            rate.tail_fraction = None
        return case

    threshold = resolve_tail_threshold_eV(config, cross_sections)
    widths = widths_from_centers(case.energy_eV)
    probability = tail_probability(case.energy_eV, case.eedf, threshold, widths)
    cutoff_threshold = (
        0.9 * float(np.max(case.energy_eV)) if len(case.energy_eV) else threshold
    )

    max_tail = 0.0
    max_cutoff = 0.0
    dominant = ""
    for rate in case.rates:
        if rate.process_type not in REACTION_TYPES:
            rate.tail_fraction = None
            continue
        proc = _matching_process(rate, cross_sections)
        if proc is None:
            rate.tail_fraction = None
            continue
        fraction = rate_tail_fraction(case.energy_eV, case.eedf, proc, threshold, widths)
        cutoff_fraction = rate_tail_fraction(
            case.energy_eV, case.eedf, proc, cutoff_threshold, widths
        )
        rate.tail_fraction = fraction
        if np.isfinite(fraction) and fraction > max_tail:
            max_tail = float(fraction)
            dominant = f"{rate.species}:{rate.process}"
        if np.isfinite(cutoff_fraction):
            max_cutoff = max(max_cutoff, float(cutoff_fraction))

    metrics = TailMetrics(
        tail_threshold_eV=float(threshold),
        tail_probability=float(probability),
        tail_rate_fraction_max=float(max_tail),
        dominant_tail_process=dominant,
        high_energy_cutoff_rate_fraction=float(max_cutoff),
        energy_grid_tail_status=_status(
            max_tail, max_cutoff, policy.tail_rate_warning_fraction
        ),
    )
    case.metadata.update(metrics.as_metadata())
    return case


def enrich_tail_metrics(
    cases: list[SwarmCaseResult],
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
) -> list[SwarmCaseResult]:
    for case in cases:
        attach_tail_metrics(case, config, cross_sections)
    return cases
