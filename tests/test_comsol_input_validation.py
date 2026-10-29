from __future__ import annotations

from hashlib import sha256
from pathlib import Path

import pytest

from swarm_workflow.comsol.input.validation import (
    ComsolInputArtifactError,
    validate_manifest_tables,
)


def _table(tmp_path: Path) -> tuple[Path, dict[str, object]]:
    path = tmp_path / "tables" / "transport.csv"
    path.parent.mkdir()
    path.write_text("energy,mobility\n1,2\n", encoding="utf-8")
    return path, {
        "tables": {
            "tables/transport.csv": {
                "columns": ["energy", "mobility"],
                "sha256": sha256(path.read_bytes()).hexdigest(),
                "artifact_role": "canonical_comsol_coefficient_input",
            }
        }
    }


def test_input_table_contract_verifies_nested_csv_and_returns_bound_hash(
    tmp_path: Path,
) -> None:
    path, manifest = _table(tmp_path)
    artifacts = validate_manifest_tables(tmp_path, manifest, require_roles=True)
    artifact = artifacts["tables/transport.csv"]
    assert artifact.path == path
    assert artifact.columns == ("energy", "mobility")
    assert artifact.sha256 == sha256(path.read_bytes()).hexdigest()


@pytest.mark.parametrize(
    "name",
    ("../transport.csv", "/transport.csv", "C:transport.csv", "tables\\transport.csv", "tables/./transport.csv"),
)
def test_input_table_contract_rejects_noncanonical_paths(tmp_path: Path, name: str) -> None:
    _path, manifest = _table(tmp_path)
    manifest["tables"][name] = manifest["tables"].pop("tables/transport.csv")
    with pytest.raises(ComsolInputArtifactError, match="artifact path"):
        validate_manifest_tables(tmp_path, manifest)


@pytest.mark.parametrize(
    ("mutation", "error"),
    (
        ("hash_missing", "valid SHA256"),
        ("hash_mismatch", "hash mismatch"),
        ("columns_differ", "columns differ"),
        ("columns_duplicate", "manifest columns"),
        ("columns_missing", "manifest columns"),
        ("row_width", "row width"),
        ("role_missing", "artifact role"),
        ("missing_file", "missing artifact"),
    ),
)
def test_input_table_contract_fails_closed(tmp_path: Path, mutation: str, error: str) -> None:
    path, manifest = _table(tmp_path)
    metadata = manifest["tables"]["tables/transport.csv"]
    if mutation == "hash_missing":
        del metadata["sha256"]
    elif mutation == "hash_mismatch":
        path.write_text("energy,mobility\n1,200\n", encoding="utf-8")
    elif mutation == "columns_differ":
        metadata["columns"] = ["mobility", "energy"]
    elif mutation == "columns_duplicate":
        metadata["columns"] = ["energy", "energy"]
    elif mutation == "columns_missing":
        del metadata["columns"]
    elif mutation == "row_width":
        path.write_text("energy,mobility\n1,2,3\n", encoding="utf-8")
        metadata["sha256"] = sha256(path.read_bytes()).hexdigest()
    elif mutation == "role_missing":
        del metadata["artifact_role"]
    elif mutation == "missing_file":
        path.unlink()
    with pytest.raises(ComsolInputArtifactError, match=error):
        validate_manifest_tables(tmp_path, manifest, require_roles=True)


def test_input_table_contract_rejects_missing_and_unknown_inventory(tmp_path: Path) -> None:
    _path, manifest = _table(tmp_path)
    with pytest.raises(ComsolInputArtifactError, match="omits required"):
        validate_manifest_tables(tmp_path, manifest, required_names=("quality.csv",))
    with pytest.raises(ComsolInputArtifactError, match="unsupported tables"):
        validate_manifest_tables(tmp_path, manifest, allowed_names=("quality.csv",))
