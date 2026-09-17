"""Preparation of provenance-bound COMSOL Function-EEDF audit plans."""

from __future__ import annotations

import hashlib
from pathlib import Path

from ...._io import write_csv
from ..function_eedf import ComsolEedfImportContract, FunctionEedfError
from .contracts import _QUERY_COLUMNS, ComsolEedfAuditError, ComsolEedfAuditPlan
from .io import read_native_eedf_grid
from .java import render_comsol_eedf_audit_java
from .sampling import _audit_queries


def prepare_comsol_eedf_audit(
    *,
    model_path: str | Path,
    table_path: str | Path,
    output_directory: str | Path,
    component: str,
    physics: str,
    function_tag: str,
    class_name: str = "SwarmComsolEedfAudit",
    moment_mean_count: int = 9,
    require_model_exists: bool = True,
    write_java: bool = True,
) -> ComsolEedfAuditPlan:
    """Create deterministic probe points and a read-only COMSOL Java audit."""

    model = Path(model_path).resolve()
    table = Path(table_path).resolve()
    output = Path(output_directory).resolve()
    if require_model_exists and not model.exists():
        raise ComsolEedfAuditError(f"saved MPH does not exist: {model}")
    try:
        import_contract = ComsolEedfImportContract.from_path(table)
    except FunctionEedfError as exc:
        raise ComsolEedfAuditError(str(exc)) from exc
    grid = read_native_eedf_grid(table)
    rows = _audit_queries(grid, moment_mean_count=moment_mean_count)
    output.mkdir(parents=True, exist_ok=True)
    query_path = output / "comsol_eedf_audit_queries.csv"
    values_path = output / "comsol_eedf_audit_values.csv"
    contract_path = output / "comsol_eedf_audit_contract.tsv"
    java_path = output / f"{class_name}.java"
    write_csv(query_path, _QUERY_COLUMNS, rows)
    table_sha256 = hashlib.sha256(table.read_bytes()).hexdigest()
    query_sha256 = hashlib.sha256(query_path.read_bytes()).hexdigest()
    if write_java:
        java_path.write_text(
            render_comsol_eedf_audit_java(
                class_name=class_name,
                model_path=model,
                component=component,
                physics=physics,
                function_tag=function_tag,
                audit_rows=rows,
                table_sha256=table_sha256,
                query_sha256=query_sha256,
            ),
            encoding="utf-8",
        )
    return ComsolEedfAuditPlan(
        model_path=model,
        model_sha256=(
            hashlib.sha256(model.read_bytes()).hexdigest()
            if require_model_exists
            else None
        ),
        table_path=table,
        query_path=query_path,
        values_path=values_path,
        contract_path=contract_path,
        java_path=java_path,
        function_tag=function_tag,
        table_sha256=table_sha256,
        query_sha256=query_sha256,
        import_contract=import_contract,
        point_count=len(rows),
        moment_mean_count=len(
            {row["group"] for row in rows if row["kind"] == "moment"}
        ),
    )
