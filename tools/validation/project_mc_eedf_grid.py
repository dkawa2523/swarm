"""Create an analysis copy of an MC database with representable EEDF cells.

LXCat reaction thresholds are included as Monte Carlo histogram boundaries.  A
decimal threshold can differ from an otherwise identical regular-grid boundary
by one floating-point ulp.  Serialising cell centres and widths to SQLite can
then make the intervening zero-measure cell impossible to reconstruct.

This external validation utility leaves the solver database untouched.  It
conservatively removes only cells whose width is at floating-point resolution,
preserving probability mass and recording every projection in a sidecar
manifest.  It is intentionally not part of the solver or workflow runtime.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import shutil
import sqlite3
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from swarm_workflow.campaign.statistics import energy_edges_from_cells  # noqa: E402


SCHEMA = "swarm.validation.mc_eedf_grid_projection.v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True)
class Cell:
    energy_eV: float
    width_eV: float
    eedf: float
    sample_count: int | None
    effective_count: float | None

    @property
    def mass(self) -> float:
        return self.eedf * self.width_eV


def _nominal_edges(cells: list[Cell]) -> np.ndarray:
    centers = np.asarray([cell.energy_eV for cell in cells], dtype=float)
    widths = np.asarray([cell.width_eV for cell in cells], dtype=float)
    left = centers - 0.5 * widths
    right = centers + 0.5 * widths
    shared = 0.5 * (right[:-1] + left[1:])
    return np.concatenate(([max(0.0, left[0])], shared, [right[-1]]))


def _coarsen(cells: list[Cell]) -> tuple[list[Cell], dict[str, float | int]]:
    edges = _nominal_edges(cells)
    scale = np.maximum(1.0, np.abs(edges))
    tolerance = 256.0 * np.finfo(float).eps * scale
    keep = np.ones(len(edges), dtype=bool)
    keep[0] = True
    last = 0
    for index in range(1, len(edges) - 1):
        if edges[index] - edges[last] <= max(tolerance[index], tolerance[last]):
            keep[index] = False
        else:
            last = index
    if edges[-1] - edges[last] <= max(tolerance[-1], tolerance[last]):
        keep[last] = False
    keep[-1] = True
    projected_edges = edges[keep]
    if np.any(np.diff(projected_edges) <= 0.0):
        raise ValueError("projected EEDF edges are not strictly increasing")

    masses = np.zeros(len(projected_edges) - 1, dtype=float)
    sample_counts: list[int | None] = [0] * len(masses)
    effective_counts: list[float | None] = [0.0] * len(masses)
    removed_widths: list[float] = []
    for index, cell in enumerate(cells):
        left = edges[index]
        right = edges[index + 1]
        midpoint = 0.5 * (left + right)
        target = int(np.searchsorted(projected_edges, midpoint, side="right") - 1)
        target = max(0, min(target, len(masses) - 1))
        masses[target] += cell.mass
        if cell.sample_count is None:
            sample_counts[target] = None
        elif sample_counts[target] is not None:
            sample_counts[target] += int(cell.sample_count)
        if cell.effective_count is None or not math.isfinite(cell.effective_count):
            effective_counts[target] = None
        elif effective_counts[target] is not None:
            effective_counts[target] += float(cell.effective_count)
        if right - left <= max(tolerance[index], tolerance[index + 1]):
            removed_widths.append(float(cell.width_eV))

    projected: list[Cell] = []
    for index, mass in enumerate(masses):
        width = float(projected_edges[index + 1] - projected_edges[index])
        energy = float(0.5 * (projected_edges[index + 1] + projected_edges[index]))
        projected.append(
            Cell(
                energy_eV=energy,
                width_eV=width,
                eedf=float(mass / width),
                sample_count=sample_counts[index],
                effective_count=effective_counts[index],
            )
        )

    check_edges = energy_edges_from_cells(
        [cell.energy_eV for cell in projected],
        [cell.width_eV for cell in projected],
    )
    if not np.allclose(check_edges, projected_edges, rtol=0.0, atol=1.0e-13):
        raise ValueError("projected EEDF cells do not reconstruct their edges")
    before_mass = float(sum(cell.mass for cell in cells))
    after_mass = float(sum(cell.mass for cell in projected))
    before_moment = float(sum(cell.mass * cell.energy_eV for cell in cells))
    after_moment = float(sum(cell.mass * cell.energy_eV for cell in projected))
    return projected, {
        "input_cells": len(cells),
        "output_cells": len(projected),
        "removed_cells": len(cells) - len(projected),
        "maximum_removed_width_eV": max(removed_widths, default=0.0),
        "probability_mass_change": after_mass - before_mass,
        "first_moment_change_eV": after_moment - before_moment,
    }


def _project_database(source: Path, output: Path) -> list[dict[str, object]]:
    shutil.copy2(source, output)
    connection = sqlite3.connect(output)
    connection.row_factory = sqlite3.Row
    events: list[dict[str, object]] = []
    try:
        groups = connection.execute(
            """
            SELECT DISTINCT mixture_id, solver, e_over_n_Td, replicate
            FROM eedf_bins
            ORDER BY mixture_id, solver, e_over_n_Td, replicate
            """
        ).fetchall()
        with connection:
            for group in groups:
                key = (
                    int(group["mixture_id"]),
                    str(group["solver"]),
                    float(group["e_over_n_Td"]),
                    int(group["replicate"]),
                )
                rows = connection.execute(
                    """
                    SELECT energy_eV, energy_width_eV, eedf, sample_count,
                           effective_sample_count
                    FROM eedf_bins
                    WHERE mixture_id=? AND solver=? AND e_over_n_Td=? AND replicate=?
                    ORDER BY bin_index
                    """,
                    key,
                ).fetchall()
                cells = [
                    Cell(
                        energy_eV=float(row["energy_eV"]),
                        width_eV=float(row["energy_width_eV"]),
                        eedf=float(row["eedf"]),
                        sample_count=(
                            None if row["sample_count"] is None else int(row["sample_count"])
                        ),
                        effective_count=(
                            None
                            if row["effective_sample_count"] is None
                            else float(row["effective_sample_count"])
                        ),
                    )
                    for row in rows
                ]
                try:
                    energy_edges_from_cells(
                        [cell.energy_eV for cell in cells],
                        [cell.width_eV for cell in cells],
                    )
                    continue
                except ValueError:
                    projected, metrics = _coarsen(cells)

                connection.execute(
                    """
                    DELETE FROM eedf_bins
                    WHERE mixture_id=? AND solver=? AND e_over_n_Td=? AND replicate=?
                    """,
                    key,
                )
                connection.executemany(
                    """
                    INSERT INTO eedf_bins(
                        mixture_id, solver, e_over_n_Td, replicate, bin_index,
                        energy_eV, energy_width_eV, eedf, eepf, sample_count,
                        effective_sample_count, relative_standard_error
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    [
                        (
                            *key,
                            index,
                            cell.energy_eV,
                            cell.width_eV,
                            cell.eedf,
                            cell.eedf / math.sqrt(max(cell.energy_eV, 1.0e-300)),
                            cell.sample_count,
                            cell.effective_count,
                            (
                                1.0 / math.sqrt(cell.effective_count)
                                if cell.effective_count is not None
                                and cell.effective_count > 0.0
                                else None
                            ),
                        )
                        for index, cell in enumerate(projected)
                    ],
                )
                events.append(
                    {
                        "mixture_id": key[0],
                        "solver": key[1],
                        "e_over_n_Td": key[2],
                        "replicate": key[3],
                        **metrics,
                    }
                )
        integrity = connection.execute("PRAGMA integrity_check").fetchone()
        if integrity is None or integrity[0] != "ok":
            raise ValueError(f"projected SQLite integrity check failed: {integrity}")
    finally:
        connection.close()
    return events


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--finalize-analysis",
        action="store_true",
        help=(
            "bind an existing projection manifest to the current analysis database "
            "after downstream aggregate tables have been written"
        ),
    )
    args = parser.parse_args()

    source = args.source.resolve()
    output = args.output.resolve()
    manifest = (args.manifest or output.with_suffix(output.suffix + ".projection.json")).resolve()
    if not source.is_file():
        parser.error(f"source database does not exist: {source}")
    if source == output:
        parser.error("output must differ from source; the raw database is immutable")
    if args.finalize_analysis:
        if not output.is_file():
            parser.error(f"analysis database does not exist: {output}")
        if not manifest.is_file():
            parser.error(f"projection manifest does not exist: {manifest}")
        payload = json.loads(manifest.read_text(encoding="utf-8"))
        if payload.get("schema") != SCHEMA:
            parser.error(f"unexpected projection manifest schema: {payload.get('schema')}")
        if Path(str(payload.get("source_database"))).resolve() != source:
            parser.error("projection manifest source database does not match")
        if Path(str(payload.get("output_database"))).resolve() != output:
            parser.error("projection manifest output database does not match")
        if payload.get("source_sha256") != _sha256(source):
            parser.error("raw source database changed after projection")
        projection_hash = payload.pop(
            "output_sha256", payload.get("projection_output_sha256")
        )
        if not projection_hash:
            parser.error("projection manifest has no projection-stage output hash")
        payload["projection_output_sha256"] = projection_hash
        payload["analysis_output_sha256"] = _sha256(output)
        payload["analysis_hash_scope"] = "after downstream aggregate table materialization"
        payload["analysis_finalized"] = True
        manifest.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(json.dumps(payload, sort_keys=True))
        return
    for path in (output, manifest):
        if path.exists() and not args.force:
            parser.error(f"refusing to overwrite {path}; pass --force")
    output.parent.mkdir(parents=True, exist_ok=True)
    manifest.parent.mkdir(parents=True, exist_ok=True)

    source_hash = _sha256(source)
    events = _project_database(source, output)
    payload = {
        "schema": SCHEMA,
        "status": "projected" if events else "unchanged_copy",
        "source_database": str(source),
        "source_sha256": source_hash,
        "output_database": str(output),
        "projection_output_sha256": _sha256(output),
        "projection_rule": (
            "Conservatively merge only finite-volume EEDF cells whose width is "
            "at or below 256 machine eps times the local energy scale."
        ),
        "probability_mass_preserved": all(
            abs(float(event["probability_mass_change"])) <= 1.0e-14 for event in events
        ),
        "events": events,
    }
    manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise
