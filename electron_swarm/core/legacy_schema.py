"""Schema v1 and removed public YAML rejection helpers."""

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

OBSOLETE_TOP_LEVEL_SOLVER_KEYS = {
    "boltzmann_two_term",
    "multiterm_boltzmann",
    "monte_carlo",
}

OBSOLETE_OUTPUT_FIELDS = {
    "compatibility",
    "write_summary",
    "write_eedf",
    "write_rates",
    "write_comparison",
    "write_plots",
}


def reject_removed_schema(raw: dict[str, Any], *, output_fields: set[str]) -> None:
    if raw.get("schema_version") != 2:
        raise ValueError(MIGRATION_ERROR)
    if "references" in raw:
        raise ValueError(
            f"{MIGRATION_ERROR}; references.external is benchmark-only and "
            "must not appear in product YAML"
        )
    if "run" in raw and isinstance(raw["run"], dict) and "mode" in raw["run"]:
        raise ValueError(MIGRATION_ERROR)
    present = OBSOLETE_TOP_LEVEL_SOLVER_KEYS.intersection(raw)
    if present:
        raise ValueError(
            f"{MIGRATION_ERROR}; move {sorted(present)} under solvers.*"
        )
    out_raw = raw.get("output", {}) or {}
    if isinstance(out_raw, dict):
        unknown = set(out_raw) - output_fields
        obsolete = sorted(unknown & OBSOLETE_OUTPUT_FIELDS)
        if obsolete:
            raise ValueError(
                "schema v2 writes fixed canonical outputs; remove output fields "
                f"{obsolete}"
            )
