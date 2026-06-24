"""Shared diagnostics helpers for electron swarm solver results."""

from .tail import (
    attach_tail_metrics,
    enrich_tail_metrics,
    rate_tail_fraction,
    resolve_tail_threshold_eV,
    tail_probability,
)

__all__ = [
    "attach_tail_metrics",
    "enrich_tail_metrics",
    "rate_tail_fraction",
    "resolve_tail_threshold_eV",
    "tail_probability",
]
