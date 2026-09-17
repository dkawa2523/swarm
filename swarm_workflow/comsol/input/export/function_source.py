"""Read and qualify source EEDF rows before table materialization."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

from swarm_workflow.selection import (
    ANCHOR_FALLBACK_SCHEMA,
)

from .contracts import ComsolExportError
from .writers import _read_csv, _required_float


def _read_function_eedf_source(
    source: Path,
    *,
    source_kind: str,
    source_composition: object = None,
) -> tuple[dict[float, list[tuple[float, float, float]]], dict[str, Any]] | None:
    source_rows = _read_csv(source)

    grouped: dict[float, list[tuple[float, float, float]]] = {}

    grouped_rows: dict[float, list[dict[str, float]]] = {}
    group_sources: dict[float, str] = {}
    composite_sources = (
        _composite_anchor_sources(source_composition)
        if source_kind == "composite"
        else None
    )
    observed_composite_fields: set[float] = set()

    for row in source_rows:
        mean_energy = _required_float(row, "mean_energy_eV")
        energy = _required_float(row, "electron_energy_eV")
        width = _required_float(row, "energy_width_eV")
        eedf = _required_float(row, "eedf")
        effective_source = source_kind
        if composite_sources is not None:
            field = _required_float(row, "E_over_N_Td")
            matches = [
                anchor
                for anchor in composite_sources
                if math.isclose(field, anchor, rel_tol=1.0e-12, abs_tol=1.0e-12)
            ]
            if len(matches) != 1:
                raise ComsolExportError(
                    "composite EEDF row does not match one selected E/N anchor"
                )
            selected_field = matches[0]
            observed_composite_fields.add(selected_field)
            effective_source = composite_sources[selected_field]
        if energy < 0.0:
            raise ComsolExportError("eedf.csv contains a negative electron energy")
        if eedf < 0.0:
            raise ComsolExportError("eedf.csv contains a negative EEDF value")
        if width <= 0.0:
            raise ComsolExportError("eedf.csv contains a nonpositive energy width")
        previous_source = group_sources.setdefault(mean_energy, effective_source)
        if previous_source != effective_source:
            raise ComsolExportError(
                "one mean-energy EEDF anchor contains multiple effective sources"
            )
        grouped_rows.setdefault(mean_energy, []).append(
            {
                "energy_eV": energy,
                "width_eV": width,
                "eedf": eedf,
            }
        )

    if composite_sources is not None and observed_composite_fields != set(
        composite_sources
    ):
        raise ComsolExportError(
            "composite EEDF anchors disagree with source_composition"
        )

    if not grouped_rows:
        return None

    if len(grouped_rows) < 2:
        return None

    grouped = {
        mean_energy: [
            (row["energy_eV"], row["width_eV"], row["eedf"])
            for row in sorted(rows, key=lambda item: item["energy_eV"])
        ]
        for mean_energy, rows in grouped_rows.items()
    }

    return grouped, {
        "policy": "full_source_support",
        "source_values_modified": False,
    }


def _composite_anchor_sources(source_composition: object) -> dict[float, str]:
    if (
        not isinstance(source_composition, dict)
        or source_composition.get("schema") != ANCHOR_FALLBACK_SCHEMA
        or source_composition.get("primary_solver") != "monte_carlo"
        or source_composition.get("fallback_solver") != "two_term"
        or source_composition.get("scope") != "low_e_over_n_anchor_fallback"
        or source_composition.get("component_mixing_within_anchor") is not False
        or source_composition.get("whole_closure_replacement") is not False
        or source_composition.get("high_e_over_n_monte_carlo_preserved") is not True
        or source_composition.get("postprocess_repair") is not False
    ):
        raise ComsolExportError("invalid composite source_composition contract")
    anchors = source_composition.get("anchors")
    if not isinstance(anchors, list) or not anchors:
        raise ComsolExportError("composite source_composition has no anchors")
    result: dict[float, str] = {}
    order: list[float] = []
    for entry in anchors:
        if not isinstance(entry, dict):
            raise ComsolExportError("composite source anchor is malformed")
        field = _required_composite_field(entry.get("E_over_N_Td"))
        effective_source = entry.get("effective_source")
        if effective_source not in {"monte_carlo", "two_term"} or field in result:
            raise ComsolExportError("composite source anchor selection is invalid")
        result[field] = str(effective_source)
        order.append(field)
    if order != sorted(order):
        raise ComsolExportError("composite source anchors must be strictly increasing")
    return result


def _required_composite_field(value: object) -> float:
    if isinstance(value, bool):
        raise ComsolExportError("composite source anchor E/N is invalid")
    try:
        field = float(value)
    except (TypeError, ValueError) as exc:
        raise ComsolExportError("composite source anchor E/N is invalid") from exc
    if not math.isfinite(field) or field <= 0.0:
        raise ComsolExportError("composite source anchor E/N is invalid")
    return field
