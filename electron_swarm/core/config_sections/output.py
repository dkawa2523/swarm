"""Publication and comparison parsers."""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
from typing import Any

from electron_swarm.chemistry.catalog import PathwayScreen
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
        reject_unknown_fields(
            selected, {"solver", "quantities", "pathways"}, "output.export"
        )
        quantities = selected.get("quantities")
        if not isinstance(quantities, list):
            raise ValueError(
                "output.export.quantities must be a list of reaction_rates or reaction_pathways"
            )
        pathways = None
        if "pathways" in selected:
            settings = selected["pathways"]
            label = "output.export.pathways"
            if not isinstance(settings, dict):
                raise ValueError(f"{label} must be a mapping")
            reject_unknown_fields(
                settings, {"catalog", "focus_species", "minimum_activity_s_inv"}, label
            )
            catalog = string_value(settings.get("catalog"), f"{label}.catalog")
            if not catalog.strip():
                raise ValueError(f"{label}.catalog must be a file path")
            pathways = PathwayScreen(
                base / catalog,
                settings.get("focus_species"),
                settings.get("minimum_activity_s_inv", 0.0),
            )
            # Validate exact dataset bytes at configuration load. The publisher
            # independently checks the files recorded by the saved calculation.
            from .conditions import parse_cross_sections

            files = parse_cross_sections(raw, base).files
            pathways.verify_collision_files(
                sha256(path.read_bytes()).hexdigest()
                for path in {file.path.resolve() for file in files}
            )
        export = RateExportConfig(
            canonical_solver(selected.get("solver"), "output.export.solver"),
            tuple(quantities),
            pathways,
        )
    output = OutputConfig(
        directory=out_raw.get("directory", "outputs"),
        base_name=string_value(out_raw.get("base_name", "swarm"), "output.base_name"),
        export=export,
        figures=tuple(figures),
    )
    output.directory = as_path(output.directory, base)
    return output
