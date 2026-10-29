from product_helpers import campaign_config
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
import yaml

from swarm_workflow.campaign.document import load_case_template
from electron_swarm import run
from electron_swarm.core.transport import ElectronTransport
from electron_swarm.collisions.electron_electron import (
    maxwellian_eedf_like,
    relaxation_target,
)
from electron_swarm.io import read_product_run, write_outputs
from electron_swarm.io.artifacts import read_scalar_table, write_scalar_table
from electron_swarm.orchestration.preparation import prepare_run
from swarm_workflow.campaign.config import load_workflow
from product_helpers import ROOT, base_product_config, write_config


def test_duplicate_yaml_cannot_override_a_requested_output(tmp_path):
    path = write_config(tmp_path, base_product_config(tmp_path))
    with path.open("a", encoding="utf-8") as stream:
        stream.write("\noutput:\n  base_name: replaced\n")
    with pytest.raises(ValueError, match="duplicate YAML key 'output'"):
        load_case_template(path)


def test_benchmark_loader_cannot_override_a_reference_selection(tmp_path):
    from tools.benchmarks.benchmark_ar_external_references import load_benchmark_config
    from tools.benchmarks.references.common import load_external_reference_configs

    path = tmp_path / "reference.yaml"
    path.write_text("references: {}\nreferences: {external: []}\n", encoding="utf-8")
    for loader in (load_benchmark_config, load_external_reference_configs):
        with pytest.raises(ValueError, match="duplicate YAML key 'references'"):
            loader(path)


@pytest.mark.parametrize(
    "old", ["max_collisions", "warmup_collisions", "tail_max_collisions"]
)
def test_old_collision_count_names_require_explicit_trial_units(tmp_path, old):
    data = base_product_config(tmp_path)
    data["solvers"]["monte_carlo"][old] = 2
    with pytest.raises(ValueError, match="MC budgets count trial periods"):
        load_case_template(write_config(tmp_path, data))


@pytest.mark.parametrize(
    "change", ["frequency", "length", "rounding", "disabled_comparison"]
)
def test_unused_or_lossy_product_requests_fail_before_computation(tmp_path, change):
    data = base_product_config(tmp_path)
    if change == "frequency":
        data["physics"]["field"]["time_dependent"] = {"frequency_Hz": 13.56e6}
    elif change == "length":
        data["conditions"]["length_scale_m"] = 0.01
    elif change == "rounding":
        data["output"]["float_format"] = "%.0f"
    else:
        data["comparison"]["required"] = True
    with pytest.raises(ValueError):
        load_case_template(write_config(tmp_path, data))


@pytest.mark.parametrize(
    "change",
    ["boolean_field", "string_field", "weights", "unused_mc", "duplicate_field"],
)
def test_workflow_inputs_have_the_same_numeric_and_selection_meaning(tmp_path, change):
    base = write_config(tmp_path, base_product_config(tmp_path), "base.yaml")
    data = campaign_config(base, database='swarm.sqlite', e_over_n_Td=[50.0], mixtures=[{'Ar': 1.0}])
    if change == "boolean_field":
        data["run"]["e_over_n_Td"] = [True]
    elif change == "string_field":
        data["run"]["e_over_n_Td"] = ["50"]
    elif change == "weights":
        data["conditions"]["gas_mixture"][0]["fraction"] = 100
    elif change == "unused_mc":
        data["campaign"]["mc"] = {"replicas": 3}
    else:
        data["run"]["e_over_n_Td"] = [50.0, 50.0]
    path = tmp_path / "workflow.yaml"
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    with pytest.raises(ValueError):
        load_workflow(path)
    assert not (tmp_path / "swarm.sqlite").exists()


def test_impossible_conservation_does_not_return_a_different_mean():
    with pytest.raises(ValueError, match="mean energy cannot be conserved"):
        maxwellian_eedf_like(np.array([0.5, 1.5, 2.5, 3.5]), np.ones(4), 3.5)


def test_specified_temperature_is_a_temperature_on_a_bounded_domain():
    energy = np.array([0.5, 1.5, 2.5, 3.5])
    target = relaxation_target(
        energy,
        np.ones(4),
        np.ones(4) / 4,
        conserve_mean_energy=False,
        fallback_temperature_eV=2.0,
    )
    expected = np.sqrt(energy) * np.exp(-energy / 2.0)
    assert target == pytest.approx(expected / expected.sum())
    assert np.dot(energy, target) != pytest.approx(3.0)


def test_conflicting_scalar_and_tensor_diffusion_is_not_accepted():
    with pytest.raises(ValueError, match="contradicts its tensor"):
        ElectronTransport.from_actual(
            definition="flux",
            gas_number_density_m3=1.0,
            drift_velocity_m_s=1.0,
            mobility_m2_V_s=1.0,
            diffusion_L_m2_s=1.0,
            diffusion_T_m2_s=1.0,
            electron_energy_diffusion_m2_s=6.0,
            electron_energy_diffusion_L_m2_s=7.0,
            electron_energy_diffusion_T_m2_s=8.0,
        )


def test_scalar_table_preserves_text_large_integers_and_binary64(tmp_path):
    frame = pd.DataFrame(
        {
            "text": ["NA", "001", ""],
            "number": [1.2345678901234567, 1e-300, None],
            "integer": [2**63 + 1, None, 7],
        },
        dtype=object,
    )
    path = tmp_path / "values.csv"
    schema = write_scalar_table(
        frame, path, {"text": "string", "number": "number", "integer": "integer"}
    )
    decoded = read_scalar_table(path.read_bytes(), schema)
    assert decoded.text.tolist() == ["NA", "001", ""]
    assert decoded.integer.tolist() == [2**63 + 1, None, 7]
    assert decoded.number.iloc[:2].tolist() == frame.number.iloc[:2].tolist()


def test_run_artifact_is_lossless_and_exposes_resolved_sampling(tmp_path):
    data = base_product_config(tmp_path, ["two_term", "monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "seed": 2**63 + 1,
        "particles": 8,
        "production_trial_periods": 4,
    }
    result = run(load_case_template(write_config(tmp_path, data)))
    artifact = read_product_run(tmp_path / "prod_manifest.json")
    assert artifact.metadata["outcome"]["status"] == "complete"
    sampling = artifact.metadata["resolved_solvers"]["monte_carlo"]
    assert sampling["sampling_unit"] == "null_collision_trial_periods"
    assert sampling["warmup_trial_periods"] == 0 and sampling["tail_trial_periods"] == 0
    row = artifact.tables["summary_csv"].query("solver == 'monte_carlo'").iloc[0]
    assert row.meta_monte_carlo_base_seed == 2**63 + 1
    eedf = artifact.tables["eedf_csv"]
    for case in result.cases:
        part = eedf[eedf.solver == case.solver]
        assert np.array_equal(part.eedf.to_numpy(), case.eedf)
        assert np.array_equal(part.energy_width_eV.to_numpy(), case.energy_widths_eV)


def test_solver_cannot_change_the_configuration_that_is_recorded(tmp_path, monkeypatch):
    import electron_swarm.runner as runner

    original = runner.execute_solve_plan

    def mutate(config, *args):
        cases = original(config, *args)
        config.run.case_prefix = "mutated"
        return cases

    monkeypatch.setattr(runner, "execute_solve_plan", mutate)
    config = load_case_template(write_config(tmp_path, base_product_config(tmp_path)))
    saved = deepcopy(config)
    with pytest.raises(RuntimeError, match="configuration changed"):
        run(config)
    assert config == saved
    assert not (tmp_path / "prod_manifest.json").exists()


def test_skipped_cli_has_distinct_status_and_no_simulation_tables(tmp_path):
    data = base_product_config(tmp_path)
    data["run"]["solvers"][0]["enabled"] = False
    path = write_config(tmp_path, data)
    completed = subprocess.run(
        [sys.executable, "-m", "electron_swarm", str(path)],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 3
    artifact = read_product_run(tmp_path / "prod_manifest.json")
    assert artifact.metadata["outcome"]["status"] == "skipped"
    assert set(artifact.tables) == {"solver_plan_csv"}


def test_comparison_is_published_from_the_completed_result(tmp_path):
    data = base_product_config(tmp_path, ["two_term", "multi_term"])
    data["comparison"]["enabled"] = True
    config = load_case_template(write_config(tmp_path, data))
    result = run(config)
    artifact = read_product_run(tmp_path / "prod_manifest.json")
    assert artifact.metadata["comparison"]["status"] == "available"
    assert len(artifact.tables["comparison_summary_csv"]) == 1

    config.output.directory = tmp_path / "republished"
    write_outputs(result, config.output)
    artifact = read_product_run(config.output.directory / "prod_manifest.json")
    assert len(artifact.tables["comparison_summary_csv"]) == 1
    assert artifact.metadata["publication"]["created_directory"] == str(
        config.output.directory
    )


@pytest.mark.parametrize("status", ["declared", "unattributed", "hash_mismatch"])
def test_input_attribution_is_separate_from_execution_qualification(tmp_path, status):
    data = base_product_config(tmp_path)
    xs = tmp_path / "argon.csv"
    xs.write_bytes(Path(data["cross_sections"]["files"][0]["path"]).read_bytes())
    data["cross_sections"]["files"][0]["path"] = str(xs)
    source = (
        {"url": "https://example.invalid/synthetic-test-source"}
        if status == "declared"
        else {}
    )
    digest = (
        "0" * 64 if status == "hash_mismatch" else sha256(xs.read_bytes()).hexdigest()
    )
    sidecar = xs.with_suffix(".provenance.json")
    sidecar.write_text(
        json.dumps({"source": source, "output": {"sha256": digest}}), encoding="utf-8"
    )
    prepared = prepare_run(load_case_template(write_config(tmp_path, data)))
    assert prepared.input_evidence[str(xs)]["attribution_status"] == status
    assert str(sidecar) in prepared.input_files


def test_changed_collision_input_is_rejected_before_publication(tmp_path, monkeypatch):
    import electron_swarm.runner as runner

    data = base_product_config(tmp_path)
    xs = tmp_path / "argon.csv"
    xs.write_bytes(Path(data["cross_sections"]["files"][0]["path"]).read_bytes())
    data["cross_sections"]["files"][0]["path"] = str(xs)
    original = runner.execute_solve_plan

    def mutate(config, *args):
        cases = original(config, *args)
        xs.write_bytes(xs.read_bytes() + b"\n")
        return cases

    monkeypatch.setattr(runner, "execute_solve_plan", mutate)
    with pytest.raises(RuntimeError, match="input files changed"):
        run(load_case_template(write_config(tmp_path, data)))
    assert not (tmp_path / "prod_manifest.json").exists()
