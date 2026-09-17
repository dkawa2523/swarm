"""Materialize the shared Function-EEDF table family."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .c1_tables import _write_c1_function_eedf_tables
from .function_source import _read_function_eedf_source


def _write_eedf_mean_energy_function_table(
    output_dir: Path,
    *,
    source: Path | None = None,
    source_kind: str = "",
    source_composition: object = None,
    rate_kernel_path: Path | None = None,
) -> dict[str, dict[str, Any]]:
    """Write one solver-independent canonical Function-EEDF input."""

    source = source or output_dir / "eedf.csv"
    if not source.exists():
        return {}
    prepared = _read_function_eedf_source(
        source,
        source_kind=source_kind,
        source_composition=source_composition,
    )
    if prepared is None:
        return {}
    grouped, source_support = prepared
    return _write_c1_function_eedf_tables(
        output_dir,
        source=source,
        grouped=grouped,
        source_support=source_support,
        rate_kernel_path=rate_kernel_path,
    )
