"""Canonical workflow-table construction API."""

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from .anchor_fallback import (
        CompositeTableSummary as CompositeTableSummary,
        compose_mc_anchor_fallback_tables as compose_mc_anchor_fallback_tables,
    )
    from .builder import build_tables as build_tables

__all__ = [
    "CompositeTableSummary",
    "build_tables",
    "compose_mc_anchor_fallback_tables",
]


def __getattr__(name: str) -> Any:
    if name in {"CompositeTableSummary", "compose_mc_anchor_fallback_tables"}:
        from . import anchor_fallback

        return getattr(anchor_fallback, name)
    if name == "build_tables":
        from .builder import build_tables

        return build_tables
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
