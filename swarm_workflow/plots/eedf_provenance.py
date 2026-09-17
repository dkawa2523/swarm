"""Source-manifest and Monte Carlo database provenance validation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from swarm_workflow.plots.eedf_contracts import EedfComparisonError


_MC_DATABASE_IDENTITY_KEYS = (
    "workflow_config_sha256",
    "base_config_sha256",
    "mc_sampling_plan_json",
    "mc_transport_estimator_schema_version",
    "mc_eedf_estimator_schema_version",
    "mc_seed_derivation_schema_version",
    "mc_solver_source_sha256",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_source_manifest(path: Path, *, solver: str) -> dict[str, Any]:
    if not path.is_file():
        raise EedfComparisonError(f"missing {solver} source manifest: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise EedfComparisonError(f"invalid {solver} source manifest: {path}") from exc
    if not isinstance(payload, dict) or payload.get("source") != solver:
        raise EedfComparisonError(
            f"{solver} source manifest does not identify the expected solver"
        )
    if not isinstance(payload.get("physical_context"), dict):
        raise EedfComparisonError(f"{solver} source manifest lacks physical context")
    if not isinstance(payload.get("mixture"), dict):
        raise EedfComparisonError(f"{solver} source manifest lacks mixture identity")
    hashes = payload.get("hashes")
    if not isinstance(hashes, dict) or not isinstance(
        hashes.get("cross_sections_sha256"), str
    ):
        raise EedfComparisonError(f"{solver} source manifest lacks cross-section hash")
    return payload


def verify_manifest_file(
    manifest: dict[str, Any],
    *,
    root: Path,
    name: str,
    solver: str,
) -> None:
    tables = manifest.get("tables")
    entry = tables.get(name) if isinstance(tables, dict) else None
    expected = entry.get("sha256") if isinstance(entry, dict) else None
    path = root / name
    if not isinstance(expected, str) or not path.is_file() or sha256(path) != expected:
        raise EedfComparisonError(
            f"{solver} {name} does not match its source manifest"
        )


def manifest_provenance(path: Path, manifest: dict[str, Any]) -> dict[str, object]:
    return {
        "manifest_path": str(path),
        "manifest_sha256": sha256(path),
        "physical_context": manifest["physical_context"],
        "mixture": manifest["mixture"],
        "cross_sections_sha256": manifest["hashes"]["cross_sections_sha256"],
    }


def validate_manifest_artifact(
    manifest_path: str | Path,
    artifact_path: str | Path,
    *,
    solver: str,
) -> None:
    """Verify that an artifact is the file recorded by its solver manifest."""

    resolved_manifest = Path(manifest_path).resolve()
    resolved_artifact = Path(artifact_path).resolve()
    if resolved_artifact.parent != resolved_manifest.parent:
        raise EedfComparisonError(
            f"{solver} artifact is outside its source-manifest directory"
        )
    manifest = read_source_manifest(resolved_manifest, solver=solver)
    verify_manifest_file(
        manifest,
        root=resolved_manifest.parent,
        name=resolved_artifact.name,
        solver=solver,
    )


def validate_mc_database_provenance(
    *,
    manifest_path: Path,
    manifest: dict[str, Any],
    metadata: Mapping[str, str],
    mixture_id: int,
    mixture_row: Any,
    mixture_species_rows: Sequence[Any],
) -> dict[str, object]:
    """Verify that an MC manifest identifies the selected database mixture."""

    try:
        database_context = json.loads(metadata["physical_context_json"])
        database_fractions = json.loads(str(mixture_row["fractions_json"]))
    except (KeyError, TypeError, json.JSONDecodeError) as exc:
        raise EedfComparisonError(
            "Monte Carlo database lacks canonical physical provenance"
        ) from exc
    manifest_mixture_id = manifest["mixture"].get("mixture_id")
    if isinstance(manifest_mixture_id, bool) or manifest_mixture_id != mixture_id:
        raise EedfComparisonError(
            "Monte Carlo manifest mixture_id differs from the requested mixture"
        )
    try:
        manifest_species = sorted(
            (
                str(item["species"]),
                float(item["fraction"]),
                float(item["mass_amu"]),
            )
            for item in manifest["mixture"].get("species", [])
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise EedfComparisonError(
            "Monte Carlo manifest has invalid mixture species"
        ) from exc
    database_species = [
        (
            str(row["species"]),
            float(row["fraction"]),
            float(row["mass_amu"]),
        )
        for row in mixture_species_rows
    ]
    manifest_fractions = {
        species: fraction for species, fraction, _mass in manifest_species
    }
    hashes = manifest["hashes"]
    for key in _MC_DATABASE_IDENTITY_KEYS:
        if not isinstance(hashes.get(key), str) or hashes[key] != metadata.get(key):
            raise EedfComparisonError(
                f"Monte Carlo database and source manifest differ in {key}"
            )
    tail_key = "mc_tail_estimator_schema_version"
    if tail_key in hashes or tail_key in metadata:
        if hashes.get(tail_key) != metadata.get(tail_key):
            raise EedfComparisonError(
                "Monte Carlo database and source manifest differ in " + tail_key
            )
    if (
        database_context != manifest["physical_context"]
        or database_fractions != manifest_fractions
        or database_species != manifest_species
        or metadata.get("cross_sections_sha256")
        != hashes["cross_sections_sha256"]
    ):
        raise EedfComparisonError(
            "Monte Carlo database and source manifest describe different physics"
        )
    return manifest_provenance(manifest_path, manifest)
