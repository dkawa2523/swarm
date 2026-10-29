"""Check downloadable model bytes and their distinction from prepared inputs."""

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
def test_prepared_inputs_do_not_claim_a_materialized_or_accepted_model(target):
    directory = MODELS / target
    preparation = _json(directory / "preparation.json")
    assert preparation["source_solver"] == "two_term"
    assert preparation["plasma_solve_executed"] is False
    assert preparation["physical_target_accepted"] is False
    assert not list(directory.glob("*.mph"))
    assert preparation["mph_materialization"] in {
        "blocked_license_expired", "blocked_preflight"
    }
    for record in preparation["files"]:
        _check_file(directory, record)
    _check_file(ROOT, preparation["source_model"])
    configuration = yaml.safe_load((directory / "swarm.yaml").read_text())
    assert configuration["schema_version"] == 2
    assert configuration["run"]["solvers"] == [{"id": "two_term"}]
    assert len(configuration["run"]["e_over_n_Td"]) == preparation["case_count"]
    bundle = _json(directory / "bundle/manifest.json")
    assert bundle["format_version"] == 2
    assert bundle["valid_ranges"] == preparation["valid_ranges"]
    database = next(r for r in preparation["files"] if r["file"] == "swarm_results.sqlite")
    if _check_file(directory, database):
        uri = (directory / "swarm_results.sqlite").as_uri() + "?mode=ro"
        with sqlite3.connect(uri, uri=True) as connection:
            fields = connection.execute(
                "SELECT solver, e_over_n_Td FROM cases ORDER BY e_over_n_Td"
            ).fetchall()
        assert len(fields) == preparation["case_count"]
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
