"""Verify portable solver-selection evidence embedded in a COMSOL bundle."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from swarm_workflow.selection import (
    ClosureSelectionError,
    contained_file,
    file_sha256,
    read_json_object,
    read_selection,
)


SELECTION_FILE = "solver_selection.json"
SOURCE_MANIFEST_FILE = "solver_source_manifest.json"


def validate_bundle_selection(
    directory: str | Path,
    *,
    required: bool = False,
) -> dict[str, Any] | None:
    """Verify a selection from portable bundle artifacts, never source paths."""

    directory = Path(directory).resolve()
    manifest = read_json_object(directory / "manifest.json")
    binding = manifest.get("solver_selection")
    if binding is None:
        if required:
            raise ClosureSelectionError(
                "this COMSOL target requires a solver selection; export with --selection"
            )
        return None
    if not isinstance(binding, dict):
        raise ClosureSelectionError("invalid bundle solver-selection binding")
    for name, key in (
        (SELECTION_FILE, "sha256"),
        (SOURCE_MANIFEST_FILE, "source_manifest_sha256"),
    ):
        if file_sha256(contained_file(directory, name)) != binding.get(key):
            raise ClosureSelectionError(f"bundle selection provenance changed: {name}")
    selection = read_selection(directory / SELECTION_FILE)
    source = selection.get("selected_solver")
    if (
        selection.get("status") != "selected"
        or source != manifest.get("source")
        or {"accept_monte_carlo": "monte_carlo", "select_two_term": "two_term"}.get(
            selection["action"]
        )
        != source
    ):
        raise ClosureSelectionError(
            "bundle solver differs from the completed selection"
        )
    evidence = selection.get("inputs", {}).get(source, {})
    if evidence.get("manifest_sha256") != binding["source_manifest_sha256"]:
        raise ClosureSelectionError(
            "bundle source manifest differs from the selected generation"
        )
    source_manifest = read_json_object(directory / SOURCE_MANIFEST_FILE)
    for key in ("physical_context", "mixture", "hashes", "source_policy"):
        if manifest.get(key) != source_manifest.get(key):
            raise ClosureSelectionError(f"bundle changed selected source {key}")
    expected_tables = evidence.get("tables")
    if not isinstance(expected_tables, dict) or not expected_tables:
        raise ClosureSelectionError("selection lacks coefficient artifact hashes")
    for name, expected in expected_tables.items():
        if file_sha256(contained_file(directory, name)) != expected:
            raise ClosureSelectionError(
                f"bundle coefficient differs from solver selection: {name}"
            )
    return {
        "sha256": binding["sha256"],
        "selected_solver": source,
        "quality_scope": selection["quality_scope"],
        "attempt": selection["attempt"],
        "selection_scope": selection["selection_scope"],
    }
