"""Shared diagnostics helpers for electron swarm solver results."""

from .common import enrich_case_diagnostics, enrich_run_diagnostics, eedf_quality_metrics
from .benchmark import MetricTolerance, compare_cases

__all__ = [
    "MetricTolerance",
    "compare_cases",
    "eedf_quality_metrics",
    "enrich_case_diagnostics",
    "enrich_run_diagnostics",
]
