"""Primitive value parsing shared by GEC-CCP mapping sections."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

from .contracts import GecCcpWorkflowError


def mapping_value(raw: dict[str, Any], name: str) -> dict[str, Any]:
    value = raw.get(name)
    if not isinstance(value, dict):
        raise GecCcpWorkflowError(f"{name} must be a mapping")
    return value


def reject_unknown_keys(
    raw: dict[str, Any], allowed: set[str], context: str
) -> None:
    unknown = sorted(set(raw) - allowed)
    if unknown:
        raise GecCcpWorkflowError(
            f"unknown {context} key(s): {', '.join(unknown)}"
        )


def required_text(raw: dict[str, Any], name: str) -> str:
    value = raw.get(name)
    if not isinstance(value, str) or not value.strip():
        raise GecCcpWorkflowError(f"{name} must be a non-empty string")
    return value.strip()


def expected_dc_field_type(raw: dict[str, Any]) -> Literal["dc"]:
    value = required_text(raw, "expected_field_type")
    if value != "dc":
        raise GecCcpWorkflowError(
            "bundle.expected_field_type must be dc for the current GEC "
            "mean-energy closure"
        )
    return "dc"


def boolean_value(
    raw: dict[str, Any], name: str, *, default: bool
) -> bool:
    value = raw.get(name, default)
    if not isinstance(value, bool):
        raise GecCcpWorkflowError(f"run.{name} must be true or false")
    return value


def resolved_path(root: Path, value: object, name: str) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise GecCcpWorkflowError(f"{name} must be a path")
    path = Path(value)
    return (root / path).resolve() if not path.is_absolute() else path.resolve()
