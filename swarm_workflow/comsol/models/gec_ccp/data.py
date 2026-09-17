"""Shared parsing helpers for GEC-CCP tabular inputs."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any, Iterable

from swarm_workflow.comsol.models.gec_ccp.contracts import GecCcpWorkflowError


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _float_or_none(value: Any) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _required_float(row: dict[str, Any], column: str) -> float:
    value = _float_or_none(row.get(column))
    if value is None:
        raise GecCcpWorkflowError(
            f"bundle column {column} contains a nonnumeric value"
        )
    return value


def _lookup(
    rows: Iterable[dict[str, str]],
    x_column: str,
    y_column: str,
) -> tuple[list[float], list[float]]:
    values: list[tuple[float, float]] = []
    seen: set[float] = set()
    for row in rows:
        try:
            x = float(row[x_column])
            y = float(row[y_column])
        except (KeyError, TypeError, ValueError) as exc:
            raise GecCcpWorkflowError(
                f"invalid lookup columns {x_column}, {y_column}"
            ) from exc
        if not math.isfinite(x) or not math.isfinite(y) or y < 0.0:
            raise GecCcpWorkflowError(f"invalid lookup value in {y_column}")
        if x in seen:
            raise GecCcpWorkflowError(
                f"lookup {y_column} contains duplicate {x_column}: {x}"
            )
        seen.add(x)
        values.append((x, y))
    if len(values) < 2:
        raise GecCcpWorkflowError(f"lookup {y_column} needs at least two points")
    values.sort(key=lambda item: item[0])
    return [item[0] for item in values], [item[1] for item in values]
