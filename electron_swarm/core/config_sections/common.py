"""Shared syntax and scalar conversion helpers for config section parsers."""

from __future__ import annotations

from pathlib import Path
from typing import Any, cast

from electron_swarm.core.config import SolverId
from electron_swarm.core.legacy_schema import OBSOLETE_PUBLIC_NAMES
from electron_swarm.core.solver_ids import CANONICAL_SOLVER_IDS


def as_path(value: str | Path | None, base: Path | None) -> Path | None:
    if value is None:
        return None
    if not isinstance(value, (str, Path)):
        raise ValueError("path value must be a string or Path")
    path = Path(value)
    if not path.is_absolute() and base is not None:
        path = (base / path).resolve()
    return path


def strict_float(value: Any, field_name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{field_name} must be numeric")
    return float(value)


def list_float(value: Any) -> list[float]:
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return [strict_float(value, "run.e_over_n_Td")]
    if isinstance(value, list):
        return [
            strict_float(item, f"run.e_over_n_Td[{index}]")
            for index, item in enumerate(value)
        ]
    raise TypeError("run.e_over_n_Td must be a number or list of numbers")


def validate_literal(value: str, allowed: set[str], field_name: str) -> str:
    if value not in allowed:
        raise ValueError(f"{field_name} must be one of {sorted(allowed)}")
    return value


def reject_unknown_fields(
    raw: dict[str, Any], allowed: set[str], field_name: str
) -> None:
    unknown = set(raw) - allowed
    if unknown:
        raise ValueError(f"Unsupported {field_name} fields: {sorted(unknown)}")


def as_mapping_section(
    raw: dict[str, Any], key: str, field_name: str
) -> dict[str, Any]:
    value = raw.get(key, {})
    if value is None:
        value = {}
    if not isinstance(value, dict):
        raise ValueError(f"{field_name} must be a mapping")
    return value


def string_value(value: Any, field_name: str) -> str:
    if not isinstance(value, str):
        raise ValueError(f"{field_name} must be a string")
    return value


def bool_field(
    raw: dict[str, Any], key: str, default: bool, field_name: str
) -> bool:
    value = raw[key] if key in raw else default
    if not isinstance(value, bool):
        raise ValueError(f"{field_name} must be a boolean")
    return value


def float_field(
    raw: dict[str, Any], key: str, default: float, field_name: str
) -> float:
    return strict_float(raw.get(key, default), field_name)


def integer_field(
    raw: dict[str, Any],
    key: str,
    default: int,
    field_name: str,
) -> int:
    value = raw.get(key, default)
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{field_name} must be an integer")
    return int(value)


def canonical_solver(value: Any, field_name: str) -> SolverId:
    solver = string_value(value, field_name)
    if solver in OBSOLETE_PUBLIC_NAMES:
        raise ValueError(
            f"{field_name} uses obsolete solver id {solver!r}; use two_term, "
            "multi_term, monte_carlo, or propagator"
        )
    return cast(
        SolverId,
        validate_literal(solver, set(CANONICAL_SOLVER_IDS), field_name),
    )
