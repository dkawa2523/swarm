"""Create the explicit low-mean-energy COMSOL benchmark bundle.

The Swarm solve is bounded below by 0.0444 Td, while the COMSOL plasma
solution can contain boundary points whose mean energy is below the first
computed Swarm value.  COMSOL's lookup-table mode therefore needs an explicit,
auditable low-energy policy instead of implicit extrapolation.

For the added mean-energy rows this script:

* holds the computed transport closure at its lowest Swarm value;
* sets the diagnostic drift velocity to zero because the derived rows carry
  ``E/N = 0`` (the velocity column is not injected into COMSOL);
* holds the elastic collision row at its lowest Swarm value; and
* sets excitation and ionization rates/Townsend coefficients to zero.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil


LOW_MEAN_ENERGY_EV = (0.0, 0.05, 0.1, 0.2, 0.35)


def _read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")
        return list(reader.fieldnames), list(reader)


def _write_rows(
    path: Path,
    fieldnames: list[str],
    rows: list[dict[str, str]],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _prepend_transport_floor(path: Path) -> None:
    fieldnames, rows = _read_rows(path)
    if not rows:
        raise ValueError(f"transport table is empty: {path}")
    first = min(rows, key=lambda row: float(row["mean_energy_eV"]))
    added: list[dict[str, str]] = []
    for energy in LOW_MEAN_ENERGY_EV:
        row = dict(first)
        row["mean_energy_eV"] = format(energy, ".17g")
        row["E_over_N_Td"] = "0"
        row["E_over_N_V_m2"] = "0"
        if "drift_velocity_m_s" in row:
            row["drift_velocity_m_s"] = "0"
        added.append(row)
    retained = [
        row
        for row in rows
        if float(row["mean_energy_eV"]) > LOW_MEAN_ENERGY_EV[-1]
    ]
    _write_rows(path, fieldnames, added + retained)


def _prepend_rate_floor(path: Path) -> None:
    fieldnames, rows = _read_rows(path)
    if not rows:
        raise ValueError(f"rate table is empty: {path}")
    lowest_energy = min(float(row["mean_energy_eV"]) for row in rows)
    lowest_rows = [
        row for row in rows if float(row["mean_energy_eV"]) == lowest_energy
    ]
    rate_columns = (
        "rate_coefficient_m3_s",
        "mixture_weighted_rate_m3_s",
        "reduced_townsend_m2",
        "mixture_weighted_reduced_townsend_m2",
    )
    added: list[dict[str, str]] = []
    for energy in LOW_MEAN_ENERGY_EV:
        for source in sorted(lowest_rows, key=lambda row: row["process_type"]):
            row = dict(source)
            row["mean_energy_eV"] = format(energy, ".17g")
            row["E_over_N_Td"] = "0"
            row["E_over_N_V_m2"] = "0"
            if row["process_type"] != "elastic":
                for column in rate_columns:
                    row[column] = "0"
            added.append(row)
    retained = [
        row
        for row in rows
        if float(row["mean_energy_eV"]) > LOW_MEAN_ENERGY_EV[-1]
    ]
    retained.sort(key=lambda row: (float(row["mean_energy_eV"]), row["process_type"]))
    _write_rows(path, fieldnames, added + retained)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def create_bundle(source: Path, output: Path) -> Path:
    source = source.resolve()
    output = output.resolve()
    if not (source / "manifest.json").is_file():
        raise FileNotFoundError(f"source manifest is missing: {source}")
    if output.exists():
        raise FileExistsError(
            f"output already exists; remove it deliberately before rebuilding: {output}"
        )
    shutil.copytree(source, output)
    transport_path = output / "transport_vs_mean_energy.csv"
    rates_path = output / "rates_vs_mean_energy.csv"
    _prepend_transport_floor(transport_path)
    _prepend_rate_floor(rates_path)

    manifest_path = output / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest.setdefault("valid_ranges", {})["mean_energy_eV"] = [
        0.0,
        manifest["valid_ranges"]["mean_energy_eV"][1],
    ]
    manifest["postprocess"] = {
        "low_mean_energy_floor": {
            "enabled": True,
            "added_mean_energy_eV": list(LOW_MEAN_ENERGY_EV),
            "transport_policy": "constant floor from lowest computed Swarm point",
            "drift_velocity_policy": (
                "zero on derived E/N=0 rows; diagnostic column is not mapped "
                "into COMSOL"
            ),
            "elastic_rate_policy": "constant floor from lowest computed Swarm point",
            "inelastic_reaction_policy": "zero below lowest computed mean-energy point",
            "reason": (
                "avoid implicit COMSOL extrapolation below the computed Swarm "
                "mean-energy domain"
            ),
            "source_bundle_relative": Path(
                os.path.relpath(source, start=output)
            ).as_posix(),
            "source_manifest_sha256": _sha256(source / "manifest.json"),
        }
    }
    manifest["derived_table_hashes_sha256"] = {
        "transport_vs_mean_energy.csv": _sha256(transport_path),
        "rates_vs_mean_energy.csv": _sha256(rates_path),
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    print(create_bundle(args.source, args.output))


if __name__ == "__main__":
    main()
