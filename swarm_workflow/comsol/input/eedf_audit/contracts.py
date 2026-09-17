"""Stable contracts and artifact columns for COMSOL Function-EEDF audits."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..function_eedf import ComsolEedfImportContract


class ComsolEedfAuditError(ValueError):
    """Raised when a saved COMSOL Function-EEDF fails its import contract."""


@dataclass(frozen=True, slots=True)
class ComsolEedfAuditPlan:
    model_path: Path
    model_sha256: str | None
    table_path: Path
    query_path: Path
    values_path: Path
    contract_path: Path
    java_path: Path
    function_tag: str
    table_sha256: str
    query_sha256: str
    import_contract: ComsolEedfImportContract
    point_count: int
    moment_mean_count: int


_QUERY_COLUMNS = (
    "point_id",
    "kind",
    "group",
    "position",
    "electron_energy_eV",
    "mean_energy_eV",
    "expected_a",
    "expected_b",
)
_VALUE_COLUMNS = (*_QUERY_COLUMNS, "comsol_value")
_STRICT_BINDING_KINDS = frozenset(
    {"anchor", "energy_midpoint", "derivative", "outside"}
)
_PROJECTION_FIDELITY_KINDS = frozenset({"cell_center"})
_PROJECTION_SIGNIFICANCE_FRACTION = 1.0e-5
_PROJECTION_MAXIMUM_RELATIVE_ERROR_LIMIT = 0.10
_PROJECTION_NORMALIZED_RMSE_LIMIT = 0.02
