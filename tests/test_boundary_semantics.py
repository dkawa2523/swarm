from product_helpers import campaign_config

"""Regressions for successful calls that previously changed data meaning."""

from copy import deepcopy
import json
import shutil
import sqlite3

import numpy as np
import pandas as pd
import pytest

from swarm_workflow.campaign.document import load_case_template
from electron_swarm import run
from electron_swarm.core.config import AngularScatteringConfig, MomentTableConfig
from electron_swarm.core.output import ComparisonConfig, OutputConfig
from electron_swarm.core.physical_context import (
    physical_context,
    workflow_physical_context,
)
from electron_swarm.core.results import SwarmRunResult, execution_outcome
from electron_swarm.io import read_product_run, write_outputs
from electron_swarm.io.artifacts import read_scalar_table, write_scalar_table
from electron_swarm.orchestration.comparison import build_comparison_summary
from product_helpers import ROOT, base_product_config, write_config, write_moment_table
from swarm_workflow.campaign.identity import (
    WORKFLOW_IMPLEMENTATION_KEY,
    workflow_execution_identity,
)
from swarm_workflow.campaign.sweep import run_sweep


@pytest.mark.parametrize("body", ["1,2,3", "1", ""])
def test_ragged_csv_is_rejected_instead_of_inferring_an_index(body):
    schema = {
        "columns": [{"name": "a", "kind": "number"}, {"name": "b", "kind": "number"}],
        "rows": 1,
    }
    with pytest.raises(ValueError, match="row width"):
        read_scalar_table(f"a,b\n{body}\n".encode(), schema)


def test_typed_csv_uses_header_order_and_preserves_text_indices_and_integers(tmp_path):
    frame = pd.DataFrame(
        {
            "text": ["NA", "001", 'a,"b"\r\nc'],
            "integer": [2**63 + 1, None, 7],
            "number": [1.2345678901234567, 1e-300, None],
        },
        dtype=object,
        index=[90, 20, 70],
    )
    path = tmp_path / "table.csv"
    schema = write_scalar_table(
        frame, path, {"number": "number", "integer": "integer", "text": "string"}
    )
    decoded = read_scalar_table(path.read_bytes(), schema)
    assert list(decoded) == list(frame)
    assert decoded.text.tolist() == frame.text.tolist()
    assert decoded.integer.tolist() == frame.integer.tolist()
    assert decoded.number.iloc[:2].tolist() == frame.number.iloc[:2].tolist()


@pytest.mark.parametrize(
    "value,kind",
    [
        (True, "number"),
        ("1", "integer"),
        ("True", "boolean"),
        (42, "string"),
        (2**63 + 1, "number"),
        (float(2**53), "integer"),
        (float(2**53 + 2), "integer"),
    ],
)
def test_typed_csv_does_not_coerce_invalid_cells_or_replace_an_existing_file(
    tmp_path, value, kind
):
    path = tmp_path / "table.csv"
    path.write_bytes(b"existing")
    with pytest.raises(ValueError, match="column 'x', row 1"):
        write_scalar_table(
            pd.DataFrame({"x": [value]}, dtype=object), path, {"x": kind}
        )
    assert path.read_bytes() == b"existing"


@pytest.mark.parametrize("directory", [None, ""])
def test_empty_output_path_is_not_replaced_by_a_default(tmp_path, directory):
    raw = base_product_config(tmp_path)
    raw["output"]["directory"] = directory
    with pytest.raises(ValueError, match="output.directory"):
        load_case_template(write_config(tmp_path, raw))


def test_python_output_api_validates_mutated_settings_before_creating_files(tmp_path):
    output = OutputConfig(tmp_path / "new")
    output.base_name = "../escape"
    with pytest.raises(ValueError, match="filename component"):
        write_outputs(SwarmRunResult(cases=[]), output)
    assert not (tmp_path / "new").exists()


@pytest.fixture
def calculated(tmp_path):
    config = load_case_template(write_config(tmp_path, base_product_config(tmp_path)))
    return config, run(config, write=False)


def test_angular_table_contents_are_part_of_direct_and_workflow_physics(
    calculated, tmp_path
):
    config, result = calculated
    a = write_moment_table(tmp_path, "a.csv", m1=0.1)
    b = write_moment_table(tmp_path, "b.csv", m1=0.7)
    config.physics.angular_scattering = AngularScatteringConfig(
        model="moment_table",
        higher_moment_closure="table",
        moment_table=MomentTableConfig(a, provenance="model_derived"),
    )
    first = physical_context(config, "same_collisions", 50.0)
    common_first = workflow_physical_context(config)
    config.physics.angular_scattering.moment_table.path = b
    second = physical_context(config, "same_collisions", 50.0)
    assert second != first and workflow_physical_context(config) != common_first
    copy = tmp_path / "same_data_elsewhere.csv"
    shutil.copy2(a, copy)
    config.physics.angular_scattering.moment_table.path = copy
    assert physical_context(config, "same_collisions", 50.0) == first
    reference = result.cases[0]
    candidate = deepcopy(reference)
    candidate.solver = "multi_term"
    reference.physical_context, candidate.physical_context = first, second
    combined = SwarmRunResult(cases=[reference, candidate])
    [row] = build_comparison_summary(combined, ComparisonConfig(enabled=True))
    assert row["physics_differences"] == "angular_scattering.moment_table_sha256"
    assert (
        row["mean_energy_eV_relative_difference"] is None
        and row["eedf_l1_error"] is None
    )
    assert combined.metadata["comparison"]["status"] == "incompatible"


def test_optional_and_required_comparisons_use_the_same_input_eligibility(calculated):
    _, result = calculated
    reference = result.cases[0]
    candidate = deepcopy(reference)
    candidate.solver = "multi_term"
    candidate.metadata["angular_model"] = "maxent_p1"
    combined = SwarmRunResult(cases=[reference, candidate])
    comparison = ComparisonConfig(enabled=True)
    [row] = build_comparison_summary(combined, comparison)
    assert combined.metadata["comparison"]["coverage_status"] == "complete"
    assert combined.metadata["comparison"]["status"] == "incompatible"
    assert row["mean_energy_eV_relative_difference"] is None
    assert execution_outcome(combined)["status"] == "partial"
    comparison.required = True
    with pytest.raises(ValueError, match="matching angular assumptions"):
        build_comparison_summary(combined, comparison)


def test_comparison_availability_does_not_claim_ensemble_qualification(calculated):
    _, result = calculated
    candidate = deepcopy(result.cases[0])
    candidate.solver = "monte_carlo"
    combined = SwarmRunResult(cases=[result.cases[0], candidate])
    build_comparison_summary(combined, ComparisonConfig(enabled=True))
    assert combined.metadata["comparison"]["status"] == "available"
    assert (
        combined.metadata["comparison"]["qualification_status"]
        == "requires_independent_ensemble"
    )


def test_product_tables_preserve_integer_counts_when_other_cases_have_no_counts(
    calculated,
):
    config, result = calculated
    sampled = deepcopy(result.cases[0])
    sampled.solver = "monte_carlo"
    for rate in sampled.rates:
        rate.solver = sampled.solver
    count = 2**63 + 1
    sampled.eedf_counts = np.full(len(sampled.eedf), count, dtype=np.uint64)
    supplied = SwarmRunResult(cases=[result.cases[0], sampled])
    write_outputs(supplied, config.output)
    artifact = read_product_run(config.output.directory / "prod_manifest.json")
    assert artifact.metadata["evidence_scope"] == "provided_results"
    table = artifact.tables["eedf_csv"]
    assert set(table.loc[table.solver == "monte_carlo", "sample_count"]) == {count}
    assert table.loc[table.solver == "two_term", "sample_count"].isna().all()


def test_publication_uses_one_result_snapshot_if_the_caller_changes_data(
    calculated, monkeypatch
):
    import electron_swarm.io.writers as writers

    config, result = calculated
    original_mean = result.cases[0].mean_energy_eV
    original_eedf = result.cases[0].eedf.copy()
    record_id = result.execution_record.identity
    original_write = writers._write_generation

    def change_callers_data(snapshot, output):
        result.cases[0].mean_energy_eV *= 2
        result.cases[0].eedf[:] = result.cases[0].eedf[::-1]
        result.metadata["configuration"]["conditions"]["gas_temperature_K"] = 999
        return original_write(snapshot, output)

    monkeypatch.setattr(writers, "_write_generation", change_callers_data)
    writers.write_outputs(result, config.output)
    artifact = read_product_run(config.output.directory / "prod_manifest.json")
    assert artifact.metadata["evidence_scope"] == "solver_execution"
    assert artifact.metadata["execution"]["id"] == record_id
    assert artifact.metadata["configuration"]["conditions"]["gas_temperature_K"] == 300
    assert artifact.tables["summary_csv"].mean_energy_eV.iloc[0] == original_mean
    np.testing.assert_array_equal(
        artifact.tables["eedf_csv"].eedf.to_numpy(), original_eedf
    )
    assert not result.execution_record.matches(result)


def test_unused_ee_temperature_does_not_change_the_physical_problem(calculated):
    config, result = calculated
    ee = config.physics.electron_electron
    ee.enabled, ee.model = True, "fp_energy"
    first = physical_context(config, "same_collisions", 50.0)
    ee.fallback_temperature_eV = 99.0
    assert physical_context(config, "same_collisions", 50.0) == first
    ee.conserve_mean_energy = False
    assert physical_context(config, "same_collisions", 50.0) != first
    write_outputs(result, config.output)
    case = read_product_run(config.output.directory / "prod_manifest.json").metadata[
        "cases"
    ][0]
    assert "moment_table_sha256" in case["physical_context"]["angular_scattering"]
    assert "electron_electron_parameters" not in case["physical_context"]


def test_workflow_source_identity_is_portable_and_tracks_seeding_not_cli(tmp_path):
    shutil.copytree(
        ROOT / "swarm_workflow/campaign",
        tmp_path / "swarm_workflow/campaign",
        ignore=shutil.ignore_patterns("__pycache__"),
    )
    statistics = tmp_path / "electron_swarm/core/statistics.py"
    statistics.parent.mkdir(parents=True)
    shutil.copy2(ROOT / "electron_swarm/core/statistics.py", statistics)
    original = workflow_execution_identity(tmp_path)
    assert original == workflow_execution_identity()
    (tmp_path / "swarm_workflow/cli.py").write_text("# output only\n")
    assert workflow_execution_identity(tmp_path) == original
    statistics.write_text(
        statistics.read_text(encoding="utf-8") + "\n# changed statistics\n",
        encoding="utf-8",
    )
    assert workflow_execution_identity(tmp_path) != original
    shutil.copy2(ROOT / "electron_swarm/core/statistics.py", statistics)
    with (tmp_path / "swarm_workflow/campaign/aggregate.py").open("a") as stream:
        stream.write("\n# changed seed derivation\n")
    assert workflow_execution_identity(tmp_path) != original


def test_workflow_cannot_resume_after_its_execution_identity_changed(tmp_path):
    import yaml

    base = write_config(tmp_path, base_product_config(tmp_path))
    path = tmp_path / "workflow.yaml"
    path.write_text(
        yaml.safe_dump(
            campaign_config(
                base,
                database="workflow.sqlite",
                e_over_n_Td=[50.0],
                mixtures=[{"Ar": 1.0}],
            )
        )
    )
    assert run_sweep(path).cases_written == 1
    with sqlite3.connect(tmp_path / "workflow.sqlite") as connection:
        recorded = json.loads(
            connection.execute(
                "SELECT value FROM metadata WHERE key=?", (WORKFLOW_IMPLEMENTATION_KEY,)
            ).fetchone()[0]
        )
        assert recorded == workflow_execution_identity()
        recorded["sha256"] = "0" * 64
        connection.execute(
            "UPDATE metadata SET value=? WHERE key=?",
            (json.dumps(recorded, sort_keys=True), WORKFLOW_IMPLEMENTATION_KEY),
        )
    with pytest.raises(ValueError, match="workflow_implementation|workflow execution"):
        run_sweep(path)
