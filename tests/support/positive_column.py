"""Canonical provenance and channel identities for synthetic adapter fixtures."""

import csv
from hashlib import sha256
import json
import io
from pathlib import Path
import re

from swarm_workflow.comsol.input.export.contracts import TABLE_ROLES


def result_csv(text: str, *, dataset: str = "dset1", solution_index: int = 1) -> str:
    """Attach explicit identity to stationary result fixtures used by acceptance tests."""
    reader = csv.DictReader(io.StringIO(text))
    columns = [
        "dataset",
        "solution_index",
        "point_index",
        "time_s",
        "parameter_name",
        "parameter_unit",
        "parameter_value",
        *reader.fieldnames,
    ]
    output = io.StringIO()
    writer = csv.DictWriter(output, columns, lineterminator="\n")
    writer.writeheader()
    for point, row in enumerate(reader):
        writer.writerow(
            {
                "dataset": dataset,
                "solution_index": solution_index,
                "point_index": point,
                **row,
            }
        )
    return output.getvalue()


def complete_bundle_fixture(bundle: Path) -> None:
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    rates = bundle / "rates_vs_mean_energy.csv"
    with rates.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        columns = list(reader.fieldnames)
        rows = list(reader)
    for column in ("species", "process", "threshold_eV"):
        if column not in columns:
            columns.append(column)
    for row in rows:
        row.setdefault("species", "Ar")
        row.setdefault(
            "process",
            "e+Ar=>e+Ars" if row["process_type"] == "excitation" else "e+Ar=>2e+Ar+",
        )
        row.setdefault(
            "threshold_eV", "11.5" if row["process_type"] == "excitation" else "15.8"
        )
    with rates.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    with (bundle / "mean_energy_vs_en.csv").open(
        encoding="utf-8", newline=""
    ) as stream:
        mean_rows = list(csv.DictReader(stream))
    quality = bundle / "quality.csv"
    with quality.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["E_over_N_Td", "E_over_N_V_m2", "passed"])
        writer.writerows(
            (row["E_over_N_Td"], row["E_over_N_V_m2"], "1") for row in mean_rows
        )
    for table in bundle.glob("*.csv"):
        with table.open(encoding="utf-8", newline="") as stream:
            header = next(csv.reader(stream))
        metadata = manifest.setdefault("tables", {}).setdefault(table.name, {})
        metadata.update(
            columns=header,
            sha256=sha256(table.read_bytes()).hexdigest(),
            artifact_role=TABLE_ROLES[table.name],
        )
    manifest.update(
        format_version=2,
        stage="export-comsol",
        status="ok",
        source="two_term",
        source_policy={"source": "two_term", "postprocess": "none"},
        hashes={
            key: "1" * 64
            for key in (
                "workflow_config_sha256",
                "calculation_sha256",
                "cross_sections_sha256",
                "source_database_sha256",
            )
        },
        physical_context={
            "schema": "swarm_physical_context.v2",
            "gas_temperature_K": 293.15,
            "gas_number_density_m3": 1e22,
            "field": {
                "type": "dc",
                "magnetic_field": {"enabled": False, "B_T": 0.0, "angle_EB_deg": 0.0},
            },
            "angular_scattering": {
                "model": "isotropic",
                "higher_moment_closure": "zero",
            },
            "electron_electron": {"enabled": False, "model": "none"},
            "ionization": {
                "energy_sharing": "equal",
                "secondary_electron_energy_eV": 0.0,
            },
            "finite_k": {"enabled": False, "k_m_inv": None},
            "high_energy_extrapolation": "zero",
        },
        mixture={
            "mixture_id": 0,
            "species": [{"species": "Ar", "fraction": 1.0, "mass_amu": 39.948}],
        },
        valid_ranges={
            "E_over_N_Td": [
                float(mean_rows[0]["E_over_N_Td"]),
                float(mean_rows[-1]["E_over_N_Td"]),
            ]
        },
    )
    passed = manifest.get("quality_summary", {}).get("passed", True)
    manifest["quality_summary"] = {
        "passed": passed,
        "failed_points": 0 if passed else 1,
        "total_points": len(mean_rows),
    }
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")


def complete_mapping_fixture(path: Path) -> None:
    text = path.read_text(encoding="utf-8")
    text = text.replace(
        "bundle:\n",
        "physical_context:\n  gas_temperature_K: 293.15\n  gas_mixture:\n    - {species: Ar, fraction: 1.0, mass_amu: 39.948}\nbundle:\n",
        1,
    )
    text = re.sub(
        r"(bundle:\n  path: [^\n]+)", r"\1\n  source: two_term", text, count=1
    )
    text = text.replace(
        "    process_type: excitation\n",
        "    process_type: excitation\n    species: Ar\n    process: e+Ar=>e+Ars\n    threshold_eV: 11.5\n",
    )
    text = text.replace(
        "    process_type: ionization\n",
        "    process_type: ionization\n    species: Ar\n    process: e+Ar=>2e+Ar+\n    threshold_eV: 15.8\n",
    )
    path.write_text(text, encoding="utf-8")
