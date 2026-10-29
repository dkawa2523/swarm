"""Publication and comparison parsers."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from electron_swarm.core.config import (
    ComparisonConfig,
    OutputConfig,
)
from electron_swarm.core.config_sections.common import (
    as_mapping_section,
    as_path,
    bool_field,
    canonical_solver,
    reject_unknown_fields,
    string_value,
)
from electron_swarm.core.output import RateExportConfig


def parse_comparison(raw: dict[str, Any]) -> ComparisonConfig:
    cmp_raw = as_mapping_section(raw, "comparison", "comparison")
    reject_unknown_fields(
        cmp_raw,
        {
            "enabled",
            "reference_solver",
            "candidate_solvers",
            "compare_eedf",
            "required",
        },
        "comparison",
    )
    reference = cmp_raw.get("reference_solver")
    candidates = cmp_raw.get("candidate_solvers", [])
    if candidates is None:
        candidates = []
    if not isinstance(candidates, list):
        raise ValueError("comparison.candidate_solvers must be a list")
    return ComparisonConfig(
        enabled=bool_field(cmp_raw, "enabled", False, "comparison.enabled"),
        reference_solver=(
            canonical_solver(reference, "comparison.reference_solver")
            if reference is not None
            else None
        ),
        candidate_solvers=[
            canonical_solver(item, "comparison.candidate_solvers")
            for item in candidates
        ],
        compare_eedf=bool_field(
            cmp_raw, "compare_eedf", True, "comparison.compare_eedf"
        ),
        required=bool_field(cmp_raw, "required", False, "comparison.required"),
    )


def parse_output(raw: dict[str, Any], base: Path) -> OutputConfig:
    out_raw = as_mapping_section(raw, "output", "output")
    reject_unknown_fields(
        out_raw, {"directory", "base_name", "export", "figures"}, "output"
    )
    figures = out_raw.get("figures", [])
    if not isinstance(figures, list):
        raise ValueError("output.figures must be a list")
    export = None
    if "export" in out_raw:
        selected = out_raw["export"]
        if not isinstance(selected, dict):
            raise ValueError("output.export must be a mapping")
        reject_unknown_fields(selected, {"solver", "quantities"}, "output.export")
        quantities = selected.get("quantities")
        if not isinstance(quantities, list):
            raise ValueError("output.export.quantities must be [reaction_rates]")
        export = RateExportConfig(
            canonical_solver(selected.get("solver"), "output.export.solver"),
            tuple(quantities),
        )
    output = OutputConfig(
        directory=out_raw.get("directory", "outputs"),
        base_name=string_value(out_raw.get("base_name", "swarm"), "output.base_name"),
        export=export,
        figures=tuple(figures),
    )
    output.directory = as_path(output.directory, base)
    return output
