"""Reindex native coefficient points for COMSOL, without fitting or repair."""

from __future__ import annotations

import csv
from hashlib import sha256
from io import StringIO
import math
from pathlib import Path
from typing import Any

from swarm_workflow._io import write_csv
from swarm_workflow.comsol.input.validation import ValidatedTableArtifact
from swarm_workflow.tables.math import strictly_monotonic

from .contracts import MEAN_ENERGY_TABLES, TABLE_ROLES, ComsolExportError


def write_mean_energy_tables(
    directory: Path,
    artifacts: dict[str, ValidatedTableArtifact],
    source_manifest: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    """Keep each quantity's own E/N support, units, values and missing cells.

    Validate every projection before writing. Mean energy is an observed column,
    never joined from a narrower transport grid or inferred by interpolation.
    """
    mean_by_field: dict[float, float] = {}
    for row in _read_snapshot(artifacts["mean_energy_vs_en.csv"]):
        field, energy = _anchor(row)
        if field in mean_by_field:
            raise ComsolExportError("duplicate E/N in mean_energy_vs_en.csv")
        mean_by_field[field] = energy
    projections = []
    for source, target in MEAN_ENERGY_TABLES.items():
        if source not in artifacts:
            continue
        artifact = artifacts[source]
        metadata = source_manifest["tables"][source]
        if metadata.get("argument") != "E_over_N_Td":
            raise ComsolExportError(
                f"native coefficient argument must be E_over_N_Td: {source}"
            )
        if not {"E_over_N_Td", "mean_energy_eV"} <= set(artifact.columns):
            raise ComsolExportError(
                f"native coefficient lacks E/N or mean energy: {source}"
            )
        rows = _read_snapshot(artifact)
        axis: dict[float, float] = {}
        for row in rows:
            field, energy = _anchor(row)
            if field in axis and source != "rates_vs_en.csv":
                raise ComsolExportError(f"duplicate E/N in {source}")
            if mean_by_field.setdefault(field, energy) != energy:
                raise ComsolExportError(
                    f"mean energy differs between native quantities at {field:g} Td: {source}"
                )
            axis[field] = energy
        if not strictly_monotonic(axis[field] for field in sorted(axis)):
            raise ComsolExportError(
                f"{source} cannot define a single-valued mean-energy function; "
                "native E/N points remain available, but COMSOL conversion requires a monotonic axis"
            )
        columns = (
            "mean_energy_eV",
            *(name for name in artifact.columns if name != "mean_energy_eV"),
        )
        rows.sort(key=lambda row: float(row["mean_energy_eV"]))
        entry = {
            "argument": "mean_energy_eV",
            "artifact_role": TABLE_ROLES[target],
            "columns": list(columns),
            "units": dict(metadata["units"]),
            "source_table": {"file": source, "sha256": artifact.sha256},
            "projection": "reindex_native_points_without_interpolation",
        }
        if "physics_contract" in metadata:
            entry["physics_contract"] = metadata["physics_contract"]
        projections.append((target, columns, rows, entry))
    tables = {}
    for target, columns, rows, entry in projections:
        path = directory / target
        write_csv(path, columns, rows)
        tables[target] = {**entry, "sha256": sha256(path.read_bytes()).hexdigest()}
    return tables


def _read_snapshot(artifact: ValidatedTableArtifact) -> list[dict[str, str]]:
    content = artifact.path.read_bytes()
    if sha256(content).hexdigest() != artifact.sha256:
        raise ComsolExportError(
            f"source table changed during conversion: {artifact.name}"
        )
    return list(csv.DictReader(StringIO(content.decode("utf-8-sig"), newline="")))


def _anchor(row: dict[str, str]) -> tuple[float, float]:
    try:
        field, energy = float(row["E_over_N_Td"]), float(row["mean_energy_eV"])
    except (KeyError, ValueError) as error:
        raise ComsolExportError(
            "native coefficient requires numeric E/N and mean energy"
        ) from error
    if not math.isfinite(field) or not math.isfinite(energy) or energy <= 0:
        raise ComsolExportError(
            "native coefficient requires finite E/N and positive finite mean energy"
        )
    return field, energy
