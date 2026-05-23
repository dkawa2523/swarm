"""Shared diagnostics helpers for electron swarm solver results."""

from .common import enrich_case_diagnostics, enrich_run_diagnostics, eedf_quality_metrics
from .tail import (
    attach_tail_metrics,
    enrich_tail_metrics,
    rate_tail_fraction,
    resolve_tail_threshold_eV,
    tail_probability,
)

__all__ = [
    "attach_tail_metrics",
    "eedf_quality_metrics",
    "enrich_case_diagnostics",
    "enrich_run_diagnostics",
    "enrich_tail_metrics",
    "rate_tail_fraction",
    "resolve_tail_threshold_eV",
    "tail_probability",
]
