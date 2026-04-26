"""Compact, solver-neutral diagnostics for swarm case results.

The diagnostics here are intentionally lightweight: they add scalar metadata to
existing :class:`SwarmCaseResult` objects instead of introducing another public
result contract.  This keeps CSV compatibility while making two-term,
multi-term, and MC outputs easier to audit.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

from electron_swarm.core.results import SwarmCaseResult, SwarmRunResult


_EPS = 1.0e-300


@dataclass(frozen=True, slots=True)
class EedfQuality:
    normalization_error: float
    negative_fraction: float
    tail_probability: float
    edge_to_peak: float

    def as_metadata(self, prefix: str = "quality") -> dict[str, float]:
        return {
            f"{prefix}_normalization_error": self.normalization_error,
            f"{prefix}_negative_fraction": self.negative_fraction,
            f"{prefix}_tail_probability": self.tail_probability,
            f"{prefix}_edge_to_peak": self.edge_to_peak,
        }


def widths_from_centers(energy_eV: np.ndarray) -> np.ndarray:
    """Return finite-volume cell widths from strictly increasing centers."""

    energy = np.asarray(energy_eV, dtype=float)
    if energy.ndim != 1 or len(energy) == 0:
        return np.asarray([], dtype=float)
    if len(energy) == 1:
        return np.ones_like(energy)
    edges = np.empty(len(energy) + 1, dtype=float)
    edges[1:-1] = 0.5 * (energy[:-1] + energy[1:])
    edges[0] = max(0.0, energy[0] - 0.5 * (energy[1] - energy[0]))
    edges[-1] = energy[-1] + 0.5 * (energy[-1] - energy[-2])
    return np.diff(edges)


def eedf_quality_metrics(
    energy_eV: np.ndarray,
    eedf_eV_inv: np.ndarray,
    widths_eV: np.ndarray | None = None,
    *,
    tail_fraction: float = 0.10,
) -> EedfQuality:
    """Compute minimal EEDF quality diagnostics.

    The tail region is the upper ``tail_fraction`` of the active energy domain.
    This metric is deliberately simple and solver-independent; process-specific
    tail-rate sensitivity can be added by solvers that already own collision
    frequencies.
    """

    energy = np.asarray(energy_eV, dtype=float)
    eedf = np.asarray(eedf_eV_inv, dtype=float)
    if widths_eV is None:
        widths = widths_from_centers(energy)
    else:
        widths = np.asarray(widths_eV, dtype=float)
    if len(energy) == 0 or len(eedf) != len(energy) or len(widths) != len(energy):
        return EedfQuality(float("nan"), float("nan"), float("nan"), float("nan"))

    integral = float(np.sum(eedf * widths))
    negative_weight = float(np.sum(np.abs(np.minimum(eedf, 0.0)) * widths))
    total_weight = float(np.sum(np.abs(eedf) * widths))
    negative_fraction = negative_weight / max(total_weight, _EPS)
    cutoff = float(np.max(energy)) * max(0.0, 1.0 - tail_fraction)
    tail_mask = energy >= cutoff
    tail_probability = float(np.sum(np.clip(eedf[tail_mask], 0.0, None) * widths[tail_mask]))
    edge_to_peak = float(eedf[-1] / max(float(np.max(np.abs(eedf))), _EPS))
    return EedfQuality(
        normalization_error=abs(integral - 1.0),
        negative_fraction=negative_fraction,
        tail_probability=tail_probability,
        edge_to_peak=edge_to_peak,
    )


def _transport_metadata(case: SwarmCaseResult) -> dict[str, str | bool | float]:
    if case.transport is None:
        return {
            "transport_has_transport_set": False,
            "transport_has_bulk": False,
        }
    meta = case.transport.metadata
    data: dict[str, str | bool | float] = {
        "transport_has_transport_set": True,
        "transport_has_bulk": case.transport.bulk is not None,
    }
    if meta is not None:
        data.update(
            {
                "transport_coefficient_definition": meta.coefficient_definition,
                "transport_swarm_condition": meta.swarm_condition,
                "transport_notes": "; ".join(meta.notes),
            }
        )
    return data


def enrich_case_diagnostics(case: SwarmCaseResult) -> SwarmCaseResult:
    """Attach compact scalar diagnostics to a case in-place and return it."""

    widths = widths_from_centers(case.energy_eV)
    quality = eedf_quality_metrics(case.energy_eV, case.eedf, widths)
    for key, value in quality.as_metadata().items():
        case.metadata.setdefault(key, value)
    for key, value in _transport_metadata(case).items():
        case.metadata.setdefault(key, value)
    case.metadata.setdefault("grid_n_cells", int(len(case.energy_eV)))
    if len(case.energy_eV):
        case.metadata.setdefault("grid_min_eV", float(np.min(case.energy_eV)))
        case.metadata.setdefault("grid_max_eV", float(np.max(case.energy_eV)))
    return case


def enrich_run_diagnostics(result: SwarmRunResult) -> SwarmRunResult:
    """Attach diagnostics to all cases in a run result."""

    for case in result.cases:
        enrich_case_diagnostics(case)
    result.metadata.setdefault("diagnostics_enriched", True)
    return result
