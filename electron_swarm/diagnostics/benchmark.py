"""Small benchmark comparison helpers.

These helpers compare already-computed case results.  They do not execute any
external solver; optional tools such as BOLOS/BOLSIG+ can write reference CSVs
that are converted to ``SwarmCaseResult`` elsewhere.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

from electron_swarm.core.results import SwarmCaseResult
from .common import widths_from_centers


@dataclass(frozen=True, slots=True)
class MetricTolerance:
    name: str
    rtol: float
    atol: float = 0.0


@dataclass(frozen=True, slots=True)
class MetricComparison:
    name: str
    candidate: float
    reference: float
    absolute_difference: float
    relative_difference: float
    rtol: float
    atol: float
    passed: bool


DEFAULT_TOLERANCES: tuple[MetricTolerance, ...] = (
    MetricTolerance("mean_energy_eV", 0.03),
    MetricTolerance("mobility_m2_V_s", 0.05),
    MetricTolerance("diffusion_L_m2_s", 0.10),
    MetricTolerance("net_ionization_frequency_s", 0.15),
)


def compare_scalar(
    name: str,
    candidate: float,
    reference: float,
    *,
    rtol: float,
    atol: float = 0.0,
) -> MetricComparison:
    cand = float(candidate)
    ref = float(reference)
    abs_diff = abs(cand - ref)
    rel_diff = abs_diff / max(abs(ref), 1.0e-300)
    passed = bool(np.isfinite(rel_diff) and abs_diff <= (atol + rtol * abs(ref)))
    return MetricComparison(
        name=name,
        candidate=cand,
        reference=ref,
        absolute_difference=float(abs_diff),
        relative_difference=float(rel_diff),
        rtol=float(rtol),
        atol=float(atol),
        passed=passed,
    )


def eedf_l1_error(candidate: SwarmCaseResult, reference: SwarmCaseResult) -> float:
    """Return L1 distance between EEDFs on the reference energy grid."""

    ref_e = np.asarray(reference.energy_eV, dtype=float)
    ref_f = np.asarray(reference.eedf, dtype=float)
    cand_f = np.interp(ref_e, candidate.energy_eV, candidate.eedf, left=0.0, right=0.0)
    widths = widths_from_centers(ref_e)
    return float(np.sum(np.abs(cand_f - ref_f) * widths))


def compare_cases(
    candidate: SwarmCaseResult,
    reference: SwarmCaseResult,
    tolerances: Iterable[MetricTolerance] | None = None,
) -> dict[str, object]:
    """Compare a case result against a reference using a compact schema."""

    active = tuple(tolerances or DEFAULT_TOLERANCES)
    comparisons = [
        compare_scalar(
            tol.name,
            getattr(candidate, tol.name),
            getattr(reference, tol.name),
            rtol=tol.rtol,
            atol=tol.atol,
        )
        for tol in active
    ]
    return {
        "candidate_solver": candidate.solver,
        "reference_solver": reference.solver,
        "case_id": candidate.case_id,
        "passed": all(item.passed for item in comparisons),
        "metrics": comparisons,
        "eedf_l1_error": eedf_l1_error(candidate, reference),
    }
