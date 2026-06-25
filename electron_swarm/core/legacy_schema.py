"""Schema v1 migration rejection helpers."""

from __future__ import annotations

from typing import Any


MIGRATION_ERROR = (
    "schema v2 is required; use schema_version: 2 and run.solvers with "
    "two_term, multi_term, or monte_carlo"
)

OBSOLETE_PUBLIC_NAMES = {
    "both",
    "all",
    "boltzmann_two_term",
    "multiterm_boltzmann",
}

def reject_removed_schema(raw: dict[str, Any]) -> None:
    if raw.get("schema_version") != 2:
        raise ValueError(MIGRATION_ERROR)
    if "run" in raw and isinstance(raw["run"], dict) and "mode" in raw["run"]:
        raise ValueError(MIGRATION_ERROR)
