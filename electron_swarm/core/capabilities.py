"""Solver capability matrix for product schema v2."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class SupportLevel(str, Enum):
    """How faithfully an implementation executes a requested feature.

    ``EXACT`` means that the selected numerical model is executed without an
    implementation fallback.  It does not claim that the selected physical
    model is exact; model fidelity is recorded separately.
    """

    EXACT = "exact"
    APPROXIMATE = "approximate"
    UNSUPPORTED = "unsupported"


class ModelFidelity(str, Enum):
    """Evidence behind a selected angular model, independent of support."""

    INTEGRAL_XS_CLOSURE = "integral_xs_closure"
    MODEL_DERIVED_MOMENTS = "model_derived_moments"
    DCS_DERIVED_MOMENTS = "dcs_derived_moments"
    UNKNOWN_MOMENT_PROVENANCE = "unknown_moment_provenance"


@dataclass(frozen=True, slots=True)
class SolverCapabilities:
    solver: str
    angular_scattering: SupportLevel
    ionization_source: SupportLevel
    electron_electron: SupportLevel
    magnetic_field: SupportLevel
    rf_field: SupportLevel
    tail_refinement: SupportLevel
