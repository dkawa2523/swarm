"""Typed contracts for the 1D positive-column adapter."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

from ...runtime import ComsolExecutionSummary
from .config import ComsolModelMapping

if TYPE_CHECKING:
    from .verify import VerifyComsolSummary


class PositiveColumnWorkflowError(RuntimeError):
    """Raised when the positive-column E2E workflow cannot proceed."""


@dataclass(frozen=True, slots=True)
class BundleValidationSummary:
    bundle_path: Path
    valid_e_over_n_Td: tuple[float, float]
    quality_passed: bool
    tables_checked: tuple[str, ...]
    manifest_sha256: str
    table_sha256: tuple[tuple[str, str], ...]


@dataclass(frozen=True, slots=True)
class PositiveColumnModelConditions:
    pressure_Pa: float | None
    gas_temperature_K: float | None
    mesh_elements: int | None


@dataclass(frozen=True, slots=True)
class PositiveColumnPlan:
    mapping: ComsolModelMapping
    original_mapping_path: Path
    effective_mapping_path: Path
    bundle: BundleValidationSummary
    model_conditions: PositiveColumnModelConditions
    output_dir: Path
    apply_java_path: Path
    verify_java_path: Path
    run_java_path: Path
    result_java_path: Path
    result_csv_path: Path
    steps: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class PositiveColumnRunSummary:
    plan: PositiveColumnPlan
    apply: ComsolExecutionSummary
    verify: VerifyComsolSummary
    run: ComsolExecutionSummary
    export: ComsolExecutionSummary
    run_summary_json: Path
