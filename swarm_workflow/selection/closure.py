"""Physical context and immutable evidence for whole-closure solver selection."""

from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping

from electron_swarm import SwarmConfig
from electron_swarm.physics.kinetics import gas_number_density


PHYSICAL_CONTEXT_KEY = "physical_context_json"
MC_QUALIFICATION_TABLE = "mc_qualification.csv"


class ClosureSelectionError(ValueError):
    """A selection must refer to unchanged, compatible physical inputs."""


def file_sha256(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def physical_context(config: SwarmConfig) -> dict[str, Any]:
    """Exclude numerical controls and solver ids from the shared physics."""
    field = asdict(config.physics.field)
    if field["type"] == "dc":
        field.pop("time_dependent")
    else:
        for key in ("phase_steps", "max_periods", "periodic_tolerance"):
            field["time_dependent"].pop(key)
    angular = asdict(config.physics.angular_scattering)
    if angular["moment_table"] is not None:
        table = angular["moment_table"]
        table["sha256"] = file_sha256(table.pop("path"))
    return {
        "schema": "swarm_physical_context.v1",
        "gas_temperature_K": config.conditions.gas_temperature_K,
        "gas_number_density_m3": gas_number_density(config),
        "field": field,
        "angular_scattering": angular,
        "electron_electron": asdict(config.physics.electron_electron),
        "ionization": asdict(config.physics.ionization),
        "finite_k": asdict(config.physics.finite_k),
        "high_energy_extrapolation": config.cross_sections.high_energy_extrapolation,
    }


def context_from_metadata(metadata: Mapping[str, str]) -> dict[str, Any] | None:
    if PHYSICAL_CONTEXT_KEY not in metadata:
        return None
    result = json.loads(metadata[PHYSICAL_CONTEXT_KEY])
    if (
        not isinstance(result, dict)
        or result.get("schema") != "swarm_physical_context.v1"
    ):
        raise ClosureSelectionError("invalid physical-context provenance")
    return result


def read_json_object(path: Path) -> dict[str, Any]:
    try:
        result = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, ValueError) as exc:
        raise ClosureSelectionError(f"cannot read evidence: {path}") from exc
    if not isinstance(result, dict):
        raise ClosureSelectionError(f"evidence must be a JSON object: {path}")
    return result


def contained_file(directory: Path, name: str) -> Path:
    candidate = (directory / name).resolve()
    if not candidate.is_relative_to(directory.resolve()) or not candidate.is_file():
        raise ClosureSelectionError(f"invalid or missing evidence artifact: {name}")
    return candidate


def table_evidence(
    directory: Path, manifest: dict[str, Any], quality: Path
) -> dict[str, Any]:
    listed = manifest.get("tables")
    if not isinstance(listed, dict):
        raise ClosureSelectionError("selection requires manifest-listed table evidence")
    tables: dict[str, str] = {}
    for name, entry in sorted(listed.items()):
        path = contained_file(directory, name)
        actual = file_sha256(path)
        if not isinstance(entry, dict) or entry.get("sha256") != actual:
            raise ClosureSelectionError(f"table evidence hash mismatch: {name}")
        tables[name] = actual
    if quality.name not in tables:
        raise ClosureSelectionError(
            "selection quality evidence is not in the table manifest"
        )
    return {
        "manifest": str(directory / "manifest.json"),
        "manifest_sha256": file_sha256(directory / "manifest.json"),
        "quality": str(quality),
        "quality_sha256": file_sha256(quality),
        "tables": tables,
    }


def compatible_physics(mc: dict[str, Any], fallback: dict[str, Any]) -> None:
    for label, manifest in (("monte_carlo", mc), ("two_term", fallback)):
        context = manifest.get("physical_context")
        if (
            not isinstance(context, dict)
            or context.get("schema") != "swarm_physical_context.v1"
        ):
            raise ClosureSelectionError(
                f"{label} lacks verified physical context; rebuild from recorded inputs"
            )
        if not manifest.get("hashes", {}).get("cross_sections_sha256"):
            raise ClosureSelectionError(f"{label} lacks cross-section provenance")
    comparisons = {
        "physical_context": (mc["physical_context"], fallback["physical_context"]),
        "cross_sections": (
            mc["hashes"]["cross_sections_sha256"],
            fallback["hashes"]["cross_sections_sha256"],
        ),
        "mixture": (
            mc.get("mixture", {}).get("species"),
            fallback.get("mixture", {}).get("species"),
        ),
    }
    for key, (left, right) in comparisons.items():
        if left is None or left != right:
            raise ClosureSelectionError(
                f"MC and two_term have different or missing {key}"
            )


def read_selection(path: str | Path) -> dict[str, Any]:
    result = read_json_object(Path(path))
    if (
        result.get("format_version") != 2
        or result.get("stage") != "mc-solver-selection"
        or result.get("selection_scope") != "whole_comsol_closure"
        or result.get("action")
        not in {
            "accept_monte_carlo",
            "select_two_term",
            "extend_time",
            "add_replicas",
            "increase_particles",
            "extend_tail",
            "blocked",
        }
    ):
        raise ClosureSelectionError(
            "unsupported solver-selection artifact; run decide-mc with current evidence"
        )
    return result


def validate_selected_tables(
    selection_path: str | Path, directory: Path
) -> dict[str, Any]:
    selection = read_selection(selection_path)
    source = selection.get("selected_solver")
    if (
        selection.get("status") != "selected"
        or {"accept_monte_carlo": "monte_carlo", "select_two_term": "two_term"}.get(
            selection["action"]
        )
        != source
        or source not in {"monte_carlo", "two_term"}
    ):
        raise ClosureSelectionError(
            "COMSOL export requires a completed whole-closure solver selection"
        )
    manifest = read_json_object(directory / "manifest.json")
    evidence = selection.get("inputs", {}).get(source, {})
    if manifest.get("source") != source or evidence.get(
        "manifest_sha256"
    ) != file_sha256(directory / "manifest.json"):
        raise ClosureSelectionError(
            "selected source or table generation differs from the export input"
        )
    for name, expected in evidence.get("tables", {}).items():
        if file_sha256(contained_file(directory, name)) != expected:
            raise ClosureSelectionError(
                f"selected table changed after decision: {name}"
            )
    if not evidence.get("tables"):
        raise ClosureSelectionError("selection lacks coefficient artifact hashes")
    return selection
