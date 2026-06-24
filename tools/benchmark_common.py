"""Shared helpers for benchmark and audit scripts."""

from __future__ import annotations

import copy
import csv
from pathlib import Path
from typing import Any

from electron_swarm.core.config import RequestedSolverConfig, SwarmConfig


def write_csv(
    path: Path,
    rows: list[dict[str, object]],
    fields: list[str] | None = None,
    *,
    extrasaction: str = "ignore",
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = fields or (list(rows[0]) if rows else ["status"])
    with path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(
            fp,
            fieldnames=fieldnames,
            extrasaction=extrasaction,
        )
        writer.writeheader()
        writer.writerows(rows)


def variant_config(
    cfg: SwarmConfig,
    solver: str,
    *,
    method: str | None = None,
    lmax: int | None = None,
    e_over_n_Td: float | None = None,
) -> SwarmConfig:
    variant = copy.deepcopy(cfg)
    variant.run.solvers = [RequestedSolverConfig(id=solver)]  # type: ignore[arg-type]
    if e_over_n_Td is not None:
        variant.run.e_over_n_Td = [float(e_over_n_Td)]
    if solver == "multi_term":
        variant.solvers.multi_term.method = method or "pn_closure_direct"  # type: ignore[assignment]
        if lmax is not None:
            variant.solvers.multi_term.lmax = int(lmax)
    return variant


def metadata_value(metadata: dict[str, object], key: str, default: Any = "") -> Any:
    if key in metadata:
        return metadata[key]
    prefixed = f"meta_{key}"
    if prefixed in metadata:
        return metadata[prefixed]
    for section in ("mc_run", "mc_energy_audit", "mc_tail_audit"):
        values = metadata.get(section)
        if isinstance(values, dict) and key in values:
            return values[key]
    return default


def case_audit_metadata(case: object) -> dict[str, object]:
    metadata = dict(getattr(case, "metadata", {}))
    audit = getattr(case, "diagnostics", {}).get("internal_monte_carlo_audit")
    if isinstance(audit, dict):
        metadata.update(audit)
    return metadata
