"""Check downloadable models, solver inputs, and physical acceptance records."""

import hashlib
import json
from pathlib import Path
import sqlite3
from zipfile import ZipFile

import pytest
import yaml


ROOT = Path(__file__).resolve().parents[1]
MODELS = ROOT / "comsol_models"


def _json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def _check_file(base, record):
    path = (base / record["file"]).resolve()
    assert path.is_relative_to(base.resolve())
    with path.open("rb") as stream:
        prefix = stream.read(128)
        stream.seek(0)
        if prefix.startswith(b"version https://git-lfs.github.com/spec/v1\n"):
            # CI can inspect the committed pointer without downloading MPH data.
            fields = dict(line.split(" ", 1) for line in path.read_text().splitlines())
            assert fields["oid"] == "sha256:" + record["sha256"]
            assert int(fields["size"]) == record["size_bytes"]
            return False
        digest = hashlib.sha256()
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    assert digest.hexdigest() == record["sha256"], path
    assert path.stat().st_size == record["size_bytes"], path
    return True


def test_distributed_source_models_retain_the_original_bytes():
    sources = _json(MODELS / "sources.json")
    assert len(sources["models"]) == 5
    for record in sources["models"]:
        _check_file(ROOT, record)


@pytest.mark.parametrize("target", ["positive_column", "gec_ccp", "gec_icp"])
def test_configured_model_has_accepted_results_and_qualified_input(target):
    directory = MODELS / target
    release = _json(directory / "release.json")
    assert release["source_solver"] == "two_term"
    assert release["physical_target_accepted"] is True
    assert release["gui_recompute"]["source_csv_paths_unavailable"] is True
    materialized = {
        record["file"]: _check_file(directory, record)
        for record in release["files"]
    }
    _check_file(ROOT, release["source_model"])
    model = release["model"]
    has_mph = _check_file(directory, model)
    acceptance = _json(directory / release["acceptance_record"])
    if target == "positive_column":
        assert acceptance["status"] == "completed" and acceptance["quality_passed"]
    else:
        assert acceptance["physical_target_accepted"] is True
        assert acceptance["quality_accepted_for_declared_closure"] is True
        if target == "gec_ccp":
            assert acceptance["output_mph"]["sha256"] == model["sha256"]
        else:
            export = _json(directory / "results/comparison_export_manifest.json")
            assert export["input_mph"]["sha256"] == model["sha256"]
    configuration = yaml.safe_load((directory / "swarm.yaml").read_text())
    assert configuration["schema_version"] == 2
    assert configuration["run"]["solvers"] == [{"id": "two_term"}]
    assert len(configuration["run"]["e_over_n_Td"]) == release["case_count"]
    bundle = _json(directory / "bundle/manifest.json")
    assert bundle["format_version"] == 2
    assert bundle["valid_ranges"] == release["valid_ranges"]
    database = next(r for r in release["files"] if r["file"] == "swarm_results.sqlite")
    assert bundle["hashes"]["source_database_sha256"] == database["sha256"]
    if has_mph and target != "positive_column":
        with ZipFile(directory / model["file"]) as mph:
            embedded = mph.read("resources/fileResource1.csv")
        assert hashlib.sha256(embedded).hexdigest() == (
            bundle["tables"]["eedf_f0_comsol_2d.csv"]["sha256"]
        )
    if materialized["swarm_results.sqlite"]:
        uri = (directory / "swarm_results.sqlite").as_uri() + "?mode=ro"
        with sqlite3.connect(uri, uri=True) as connection:
            fields = connection.execute(
                "SELECT solver, e_over_n_Td FROM cases ORDER BY e_over_n_Td"
            ).fetchall()
        assert len(fields) == release["case_count"]
        assert {solver for solver, _ in fields} == {"two_term"}
        assert [field for _, field in fields] == sorted(configuration["run"]["e_over_n_Td"])


@pytest.mark.parametrize("source", ["two_term", "monte_carlo"])
def test_historical_mph_matches_the_final_solve_and_embedded_eedf(source):
    directory = MODELS / "saved" / f"gec_icp_{source}"
    archive = _json(directory / "archive.json")
    assert archive["status"] == "historical_solved_model"
    assert archive["current_package_compatibility_verified"] is False
    assert archive["current_physical_target_accepted"] is False
    for record in archive["files"]:
        _check_file(directory, record)
    solve = _json(directory / "evidence/solve_output_provenance.json")
    assert archive["model"]["sha256"] == solve["output_mph"]["sha256"]
    if _check_file(directory, archive["model"]):
        with ZipFile(directory / archive["model"]["file"]) as mph:
            embedded = mph.read("resources/fileResource1.csv")
        table = directory / "bundle/eedf_f0_comsol_2d.csv"
        assert hashlib.sha256(embedded).hexdigest() == hashlib.sha256(table.read_bytes()).hexdigest()
