from copy import deepcopy
from dataclasses import replace
import json
import shutil
import sqlite3

import pytest

from electron_swarm import load_config, run
from electron_swarm.cli import main
from electron_swarm.core.results import SwarmRunResult
from electron_swarm.core.source_identity import (
    solver_source_fingerprint,
    execution_source_fingerprint,
)
from electron_swarm.io import read_product_run, write_outputs
from electron_swarm.orchestration.comparison import build_comparison_summary
from electron_swarm.orchestration.preparation import prepare_run
from electron_swarm.core.config import ComparisonConfig
from product_helpers import ROOT, base_product_config, write_config
from swarm_workflow.campaign.store import WorkflowStore
from swarm_workflow.tables.repository import _field_source_policy
from swarm_workflow.tables.contracts import TableBuildError


@pytest.fixture
def executed(tmp_path):
    config = load_config(write_config(tmp_path, base_product_config(tmp_path)))
    return config, run(config, write=False)


def test_mass_is_required_and_conflicts_fail_before_execution(tmp_path):
    data = base_product_config(tmp_path)
    del data["conditions"]["gas_mixture"][0]["mass_amu"]
    with pytest.raises(ValueError, match="mass_amu is required"):
        load_config(write_config(tmp_path, data))
    data["conditions"]["gas_mixture"][0]["mass_amu"] = 28.0134
    with pytest.raises(ValueError, match="conflicting mass"):
        prepare_run(load_config(write_config(tmp_path, data)))


def test_rounded_source_masses_are_canonical_for_every_operator(tmp_path):
    data = base_product_config(tmp_path)
    data["conditions"]["gas_mixture"][0]["mass_amu"] = 39.95
    prepared = prepare_run(load_config(write_config(tmp_path, data)))
    assert {p.mass_amu for p in prepared.inputs.processes} == {39.95}
    assert prepared.physical_cases[0].species == (("Ar", 1.0, 39.95),)
    prepared.inputs.processes[0].cross_section_m2[0] += 1e-20
    with pytest.raises(RuntimeError, match="loaded collision data changed"):
        prepared.verify_unchanged()


def test_population_is_physics_and_disabled_solvers_do_not_make_partial_runs(
    tmp_path, capsys
):
    data = base_product_config(tmp_path, ["two_term", "monte_carlo"])
    data["run"]["solvers"][1]["enabled"] = False
    path = write_config(tmp_path, data)
    prepared = prepare_run(load_config(path))
    assert prepared.physical_cases[0].growth_model == "temporal_growth"
    assert main([str(path), "--no-write"]) == 0
    assert "complete:" in capsys.readouterr().out
    data["solvers"]["monte_carlo"]["population_model"] = "weighted_branching"
    with pytest.raises(ValueError, match="physics.ionization.population_model"):
        load_config(write_config(tmp_path, data))
    del data["solvers"]["monte_carlo"]["population_model"]
    data["physics"]["ionization"]["population_model"] = "fixed_particle_single_daughter"
    with pytest.raises(ValueError, match="fixed_population"):
        prepare_run(load_config(write_config(tmp_path, data)))
    data["feature_policy"]["unsupported"] = "skip_solver"
    prepared = prepare_run(load_config(write_config(tmp_path, data)))
    assert prepared.plan[0].skipped and not prepared.plan[1].requested
    data["run"]["solvers"][1]["enabled"] = True
    data["solvers"]["monte_carlo"] = {
        "seed": 7,
        "particles": 8,
        "production_trial_periods": 2,
    }
    assert main([str(write_config(tmp_path, data)), "--no-write"]) == 4
    assert "partial:" in capsys.readouterr().out


def test_comparison_matches_physics_and_each_observable(executed):
    _, result = executed
    reference = result.cases[0]
    candidate = deepcopy(reference)
    candidate.solver = "monte_carlo"
    candidate.transport = replace(
        candidate.transport, definition="mc_flux_particle_tracking_fixed_population"
    )
    candidate.transport_observables = candidate.transport.observable_definitions
    comparison = ComparisonConfig(
        enabled=True, reference_solver="two_term", candidate_solvers=["monte_carlo"]
    )
    pair = SwarmRunResult(cases=[reference, candidate])
    [row] = build_comparison_summary(pair, comparison)
    assert row["diffusion_L_m2_s_status"] == "definition_mismatch"
    assert row["diffusion_L_relative_difference"] is None
    assert row["diffusion_T_m2_s_status"] == "available_model_approximation"
    assert reference.transport_observables["diffusion_L_m2_s"].approximate
    candidate.physical_context = replace(
        candidate.physical_context, growth_model="fixed_population"
    )
    [row] = build_comparison_summary(pair, comparison)
    assert row["physics_status"] == "mismatch"
    assert row["physics_differences"] == "growth_model"
    assert (
        row["mean_energy_eV_relative_difference"] is None
        and row["eedf_l1_error"] is None
    )
    comparison.required = True
    with pytest.raises(ValueError, match="physical context"):
        build_comparison_summary(pair, comparison)


def test_modified_results_have_immutable_parent_and_no_stale_comparison(
    executed, tmp_path
):
    config, result = executed
    record = result.execution_record
    result.cases[0].mean_energy_eV *= 2
    result.metadata["configuration"]["conditions"]["gas_temperature_K"] = 999
    result.metadata["comparison_summary_rows"] = [{"fake": "old comparison"}]
    write_outputs(result, config.output)
    artifact = read_product_run(tmp_path / "prod_manifest.json")
    assert artifact.metadata["evidence_scope"] == "derived_results"
    assert (
        artifact.metadata["configuration"] is None
        and artifact.metadata["execution"] is None
    )
    parent = artifact.metadata["source_execution"]
    assert parent["id"] == record.identity
    assert parent["metadata"]["configuration"]["conditions"]["gas_temperature_K"] == 300
    assert artifact.metadata["comparison"]["status"] == "not_evaluated"
    assert "comparison_summary_csv" not in artifact.tables
    assert artifact.tables["summary_csv"].meta_numerical_status.tolist() == ["derived"]
    assert result.cases[0].metadata["numerical_status"] != "derived"
    assert (
        json.loads(record.metadata_json)["configuration"]["conditions"][
            "gas_temperature_K"
        ]
        == 300
    )


def test_provided_results_cannot_claim_execution_from_configuration(executed, tmp_path):
    config, result = executed
    supplied = SwarmRunResult(cases=result.cases, metadata=deepcopy(result.metadata))
    write_outputs(supplied, config.output)
    artifact = read_product_run(tmp_path / "prod_manifest.json")
    assert artifact.metadata["evidence_scope"] == "provided_results"
    assert artifact.metadata["configuration"] is None
    assert artifact.metadata["outcome"]["status"] == "provided"
    assert artifact.metadata["source_execution"] is None
    assert artifact.metadata["comparison"]["status"] == "not_evaluated"


def test_workflow_and_table_policy_preserve_context_and_observable_meanings(
    executed, tmp_path
):
    _, result = executed
    case = result.cases[0]
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        store.write_mixture(0, [("Ar", 1.0, 39.948)])
        store.write_case(mixture_id=0, replicate=0, case=case)
        store.connection.row_factory = sqlite3.Row
        policy = _field_source_policy(store.connection, solver="two_term", mixture_id=0)
        assert policy["physical_context"]["species"] == [["Ar", 1.0, 39.948]]
        assert "e_over_n_Td" not in policy["physical_context"]
        assert (
            policy["transport_observables"]["diffusion_L_m2_s"]["approximate"] is True
        )
        changed = deepcopy(case)
        changed.physical_context = replace(case.physical_context, gas_temperature_K=999)
        store.write_case(mixture_id=0, replicate=1, case=changed)
        with pytest.raises(TableBuildError, match="incompatible physical contexts"):
            _field_source_policy(store.connection, solver="two_term", mixture_id=0)


def test_numerical_identity_excludes_cli_output_but_tracks_model_and_execution(
    tmp_path,
):
    shutil.copytree(
        ROOT / "electron_swarm",
        tmp_path / "electron_swarm",
        ignore=shutil.ignore_patterns("__pycache__"),
    )

    def identities():
        return solver_source_fingerprint("propagator", tmp_path)[
            "sha256"
        ], execution_source_fingerprint(tmp_path)["sha256"]

    original = identities()
    for name in (
        "cli.py",
        "io/writers.py",
        "core/config_sections/output.py",
        "core/output.py",
    ):
        with (tmp_path / "electron_swarm" / name).open("a") as stream:
            stream.write("\n# presentation-only change\n")
    assert identities() == original
    with (tmp_path / "electron_swarm/runner.py").open("a") as stream:
        stream.write("\n# execution change\n")
    assert identities()[0] == original[0] and identities()[1] != original[1]
    prior_execution = identities()[1]
    with (tmp_path / "electron_swarm/core/config.py").open("a") as stream:
        stream.write("\n# input-default change\n")
    assert identities()[0] == original[0] and identities()[1] != prior_execution
    with (tmp_path / "electron_swarm/solvers/propagator/collisions.py").open(
        "a"
    ) as stream:
        stream.write("\n# model change\n")
    assert identities()[0] != original[0]
