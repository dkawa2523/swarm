"""Neutral coefficients preserve meaning from config through saved data and CSV."""

from copy import deepcopy
from dataclasses import replace
import json
import sqlite3

import numpy as np
import pandas as pd
import pytest
import yaml

from electron_swarm import load_config, run
from electron_swarm.core.output import OutputConfig, RateExportConfig
from electron_swarm.io import (
    read_product_run,
    read_reaction_rates,
    write_outputs,
    write_reaction_rates,
)
from electron_swarm.io.reaction_rates import build_reaction_rates, case_rate_sample
from electron_swarm.orchestration.preparation import prepare_run
from swarm_workflow.campaign.sweep import plan_sweep, run_sweep
from swarm_workflow.cli import main
from swarm_workflow.reaction_rates import export_reaction_rates
from product_helpers import base_product_config, write_config


def _input(tmp_path, solver="two_term"):
    data = base_product_config(tmp_path, [solver])
    data["run"]["e_over_n_Td"] = [10.0, 30.0]
    data["output"]["export"] = {"solver": solver, "quantities": ["reaction_rates"]}
    if solver == "monte_carlo":
        data["solvers"][solver] = {
            "seed": 17,
            "particles": 2,
            "production_trial_periods": 8,
        }
        data["physics"]["ionization"]["population_model"] = (
            "fixed_particle_single_daughter"
        )
    if solver == "propagator":
        data["solvers"][solver] = {
            "energy_cells": 64,
            "polar_cells": 8,
            "max_memory_mb": 128,
        }
    return data


@pytest.mark.parametrize(
    "export_request",
    [
        None,
        {},
        {"solver": "multiterm_boltzmann", "quantities": ["reaction_rates"]},
        {"solver": "two_term", "quantities": ["transport"]},
        {"solver": "two_term", "quantities": "reaction_rates"},
        {"solver": "two_term", "quantities": ["reaction_rates"], "scenario": "default"},
    ],
)
def test_invalid_exports_fail_at_loading(tmp_path, export_request):
    data = _input(tmp_path)
    data["output"]["export"] = export_request
    with pytest.raises(ValueError):
        load_config(write_config(tmp_path, data))


def test_export_plan_rejects_unavailable_solver_and_rf_before_creating_files(tmp_path):
    data = _input(tmp_path)
    data["output"]["export"]["solver"] = "multi_term"
    with pytest.raises(ValueError, match="enabled"):
        prepare_run(load_config(write_config(tmp_path, data)))
    data = _input(tmp_path)
    data["physics"]["field"]["type"] = "time_dependent"
    data["physics"]["field"]["time_dependent"] = {"frequency_Hz": 1e6}
    with pytest.raises(NotImplementedError, match="DC E/N"):
        prepare_run(load_config(write_config(tmp_path, data)))
    assert not (tmp_path / "prod_reaction_rates.csv").exists()


@pytest.mark.parametrize(
    "solver", ("two_term", "multi_term", "monte_carlo", "propagator")
)
def test_each_solver_publishes_native_k_points_and_context(tmp_path, solver):
    path = write_config(tmp_path, _input(tmp_path, solver))
    config = load_config(path)
    plan = prepare_run(config).description()
    assert plan["export"]["solver"] == solver
    assert plan["export"]["data_role"] == "coefficient_points"
    assert not (tmp_path / "prod_manifest.json").exists()
    result = run(config)
    metadata_path = result.metadata["output_paths"]["manifest_json"]
    product = read_product_run(metadata_path)
    rates = read_reaction_rates(metadata_path)
    assert list(rates.table.E_over_N_Td.unique()) == [10.0, 30.0]
    assert rates.metadata["units"]["rate_coefficient_m3_s"] == "m^3/s"
    assert rates.metadata["interpolation"] == "none"
    expected = product.tables["rates_csv"].sort_values(
        ["E_over_N_Td", "species", "process"]
    )
    actual = rates.table.sort_values(["E_over_N_Td", "species", "process"])
    np.testing.assert_array_equal(
        actual.rate_coefficient_m3_s, expected.rate_coefficient_m3_s
    )
    np.testing.assert_array_equal(
        actual.mixture_weighted_rate_m3_s, expected.mixture_weighted_rate_m3_s
    )
    for series in rates.metadata["series"].values():
        assert series["physical_context"]["gas_temperature_K"] == 300
        assert "e_over_n_Td" not in series["physical_context"]
        assert series["range_Td"] == [10.0, 30.0]
    assert rates.table.standard_error_m3_s.isna().all()
    if solver == "monte_carlo":
        assert set(rates.table.estimate_status) <= {
            "single_replica",
            "censored_all_zero",
        }


def test_saved_product_export_needs_no_transport_or_monotonic_mean_energy(tmp_path):
    config = load_config(write_config(tmp_path, _input(tmp_path)))
    result = run(config, write=False)
    for index, case in enumerate(result.cases):
        case.mean_energy_eV = [7, 3][index]
        case.transport = replace(
            case.transport,
            drift_velocity_m_s=0,
            reduced_diffusion_L_m2_s_m3=None,
            reduced_diffusion_T_m2_s_m3=None,
        )
        case.metadata["electron_electron_transport_stale"] = True
    paths = write_outputs(result, config.output)
    before = paths["manifest_json"].read_bytes()
    exported = export_reaction_rates(
        paths["manifest_json"], OutputConfig(tmp_path / "neutral"), solver="two_term"
    )
    restored = read_reaction_rates(exported["manifest_json"])
    assert restored.metadata["source"]["evidence_scope"] == "derived_results"
    assert set(restored.table.numerical_status) == {"derived"}
    assert paths["manifest_json"].read_bytes() == before
    with pytest.raises(ValueError, match="overwrite its input"):
        export_reaction_rates(paths["manifest_json"], config.output, solver="two_term")
    config.output.export = None
    write_outputs(result, config.output)
    assert not (tmp_path / "prod_reaction_rates.csv").exists()
    assert "reaction_rates_csv" not in read_product_run(paths["manifest_json"]).tables


@pytest.mark.parametrize(
    "solver", ("two_term", "multi_term", "monte_carlo", "propagator")
)
def test_campaign_exports_and_resumes_the_selected_solver(tmp_path, solver):
    data = _input(tmp_path, solver)
    data["campaign"] = {"database": "cases.sqlite"}
    if solver == "monte_carlo":
        data["campaign"]["mc"] = {"replicas": 2}
    path = write_config(tmp_path, data)
    assert plan_sweep(path)["export"]["solver"] == solver
    first = run_sweep(path)
    artifact_path = tmp_path / "prod_manifest.json"
    artifact = read_reaction_rates(artifact_path)
    expected = artifact.table.copy(deep=True)
    assert artifact.metadata["source"]["kind"] == "workflow_database"
    assert set(artifact.table.replicates) == ({2} if solver == "monte_carlo" else {1})
    assert run_sweep(path).cases_written == 0
    pd.testing.assert_frame_equal(read_reaction_rates(artifact_path).table, expected)
    with sqlite3.connect(first.database_path) as connection:
        assert connection.execute("SELECT COUNT(*) FROM cases").fetchone()[0] == (
            4 if solver == "monte_carlo" else 2
        )
    if solver == "monte_carlo":
        with sqlite3.connect(first.database_path) as connection:
            connection.execute(
                "DELETE FROM cases WHERE solver='monte_carlo' AND e_over_n_Td=10 AND replicate=1"
            )
        with pytest.raises(ValueError, match="planned MC replicas"):
            export_reaction_rates(
                first.database_path,
                OutputConfig(tmp_path / "incomplete"),
                solver=solver,
            )
        assert not (tmp_path / "incomplete").exists()


@pytest.fixture
def sample(tmp_path):
    config = load_config(write_config(tmp_path, _input(tmp_path)))
    return case_rate_sample(run(config, write=False).cases[0])


def test_mc_independent_replica_statistics_and_zero_censoring(sample):
    replicas = []
    for index in range(3):
        rates = [
            dict(row, rate_coefficient_m3_s=(index + 1) * 1e-16) for row in sample.rates
        ]
        replicas.append(
            replace(
                sample,
                solver="monte_carlo",
                rates=rates,
                replicate=index,
                seed=index + 10,
            )
        )
    artifact = build_reaction_rates(replicas, solver="monte_carlo", source={})
    row = artifact.table.iloc[0]
    assert row.rate_coefficient_m3_s == pytest.approx(2e-16)
    assert row.standard_error_m3_s == pytest.approx(1e-16 / np.sqrt(3))
    assert row.ci95_high_m3_s == pytest.approx(
        2e-16 + 4.3026527299 * row.standard_error_m3_s
    )
    assert set(artifact.table.estimate_status) == {"independent_replica_estimate"}
    for replica in replicas:
        for channel in replica.rates:
            channel["rate_coefficient_m3_s"] = 0.0
    zero = build_reaction_rates(replicas, solver="monte_carlo", source={})
    assert set(zero.table.estimate_status) == {"censored_all_zero"}
    assert zero.table.ci95_high_m3_s.isna().all()
    with pytest.raises(ValueError, match="independent replica seeds"):
        build_reaction_rates(
            [replicas[0], replicas[0]], solver="monte_carlo", source={}
        )


def test_channel_coverage_and_target_basis_are_never_guessed(sample):
    missing = replace(
        sample,
        e_over_n_Td=30,
        physical_context={**sample.physical_context, "e_over_n_Td": 30},
        rates=sample.rates[:-1],
    )
    with pytest.raises(ValueError, match="channel coverage"):
        build_reaction_rates([sample, missing], solver=sample.solver, source={})
    bad = deepcopy(sample)
    bad.rates[0]["target_species_fraction"] = 0.5
    with pytest.raises(ValueError, match="target fraction"):
        build_reaction_rates([bad], solver=sample.solver, source={})


def test_neutral_artifact_is_hash_bound_and_changed_data_loses_numerical_claim(
    sample, tmp_path
):
    artifact = build_reaction_rates(
        [sample], solver=sample.solver, source={"kind": "source"}
    )
    artifact.table.loc[:, "rate_coefficient_m3_s"] *= 2
    artifact.table.loc[:, "mixture_weighted_rate_m3_s"] *= 2
    paths = write_reaction_rates(artifact, OutputConfig(tmp_path))
    restored = read_reaction_rates(paths["manifest_json"])
    assert restored.metadata["source"]["kind"] == "derived"
    assert set(restored.table.numerical_status) == {"derived"}
    paths["reaction_rates_csv"].write_bytes(b"altered")
    with pytest.raises(ValueError, match="altered"):
        read_reaction_rates(paths["manifest_json"])


def test_saved_data_cli_and_unknown_solver_fail_clearly(tmp_path, capsys):
    data = _input(tmp_path)
    del data["output"]["export"]
    run(load_config(write_config(tmp_path, data)))
    main(
        [
            "export-reaction-rates",
            str(tmp_path / "prod_manifest.json"),
            "--solver",
            "two_term",
            "--output",
            str(tmp_path / "export"),
        ]
    )
    assert "swarm_manifest.json" in capsys.readouterr().out
    with pytest.raises(SystemExit):
        main(
            [
                "export-reaction-rates",
                str(tmp_path / "prod_manifest.json"),
                "--solver",
                "multi_term",
                "--output",
                str(tmp_path / "absent"),
            ]
        )
    assert not (tmp_path / "absent").exists()


def test_yaml_export_does_not_change_the_calculation_identity(tmp_path):
    from swarm_workflow.campaign.document import calculation_sha256

    data = _input(tmp_path)
    path = write_config(tmp_path, data)
    expected = calculation_sha256(path)
    del data["output"]["export"]
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    assert calculation_sha256(path) == expected
    with pytest.raises(ValueError):
        RateExportConfig("both")


def test_updated_ee_rates_remain_exportable_with_stale_transport(tmp_path):
    data = _input(tmp_path)
    data["physics"]["electron_electron"] = {
        "enabled": True,
        "model": "relaxation_postprocess",
        "relaxation_fraction": 0.05,
    }
    config = load_config(write_config(tmp_path, data))
    run(config)
    artifact = read_reaction_rates(tmp_path / "prod_manifest.json")
    assert set(artifact.table.numerical_status) == {"postprocessed_rates"}
    for series in artifact.metadata["series"].values():
        assert all(
            case["solver_metadata"]["electron_electron_transport_stale"]
            for case in series["cases"]
        )


def test_campaign_export_selects_after_all_solver_workers(tmp_path):
    data = _input(tmp_path, "multi_term")
    data["run"]["solvers"].append({"id": "two_term"})
    data["campaign"] = {"database": "cases.sqlite"}
    path = write_config(tmp_path, data)
    assert run_sweep(path).cases_written == 4
    artifact = read_reaction_rates(tmp_path / "prod_manifest.json")
    assert set(artifact.table.solver) == {"multi_term"}


def test_product_reader_does_not_accept_an_unknown_manifest_without_a_role(tmp_path):
    path = tmp_path / "bad_manifest.json"
    path.write_text('{"files": {}}')
    with pytest.raises(ValueError, match="unsupported"):
        read_product_run(path)


def test_adding_campaign_replaces_the_old_product_generation(tmp_path):
    data = _input(tmp_path)
    path = write_config(tmp_path, data)
    run(load_config(path))
    assert (tmp_path / "prod_manifest.json").is_file()
    data["campaign"] = {"database": "cases.sqlite"}
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    run_sweep(path)
    assert (
        json.loads((tmp_path / "prod_manifest.json").read_text(encoding="utf-8"))[
            "artifact_role"
        ]
        == "reaction_rates"
    )
    assert not (tmp_path / "prod_summary.csv").exists()
    assert not (tmp_path / "prod_rates.csv").exists()
    assert not (tmp_path / "prod_eedf.csv").exists()
    assert (
        read_reaction_rates(tmp_path / "prod_manifest.json").metadata["source"]["kind"]
        == "workflow_database"
    )


def test_rate_units_cannot_be_relabelled_and_do_not_mutate_future_exports(
    sample, tmp_path
):
    data = build_reaction_rates([sample], solver=sample.solver, source={})
    data.metadata["units"]["rate_coefficient_m3_s"] = "cm^3/s"
    with pytest.raises(ValueError, match="columns and units"):
        write_reaction_rates(data, OutputConfig(tmp_path / "invalid"))
    original = build_reaction_rates([sample], solver=sample.solver, source={})
    assert original.metadata["units"]["rate_coefficient_m3_s"] == "m^3/s"


def test_saved_export_binds_the_manifest_snapshot_it_actually_read(
    tmp_path, monkeypatch
):
    from electron_swarm.core.identity import content_digest
    import swarm_workflow.reaction_rates as module

    data = _input(tmp_path)
    run(load_config(write_config(tmp_path, data)))
    source = tmp_path / "prod_manifest.json"
    original = read_product_run(source)

    def replace_after_read(path, **_options):
        changed = deepcopy(original.metadata)
        changed["publication"]["created_directory"] = "another-generation"
        path.write_text(json.dumps(changed), encoding="utf-8")
        return original

    monkeypatch.setattr(module, "read_product_run", replace_after_read)
    paths = module.export_reaction_rates(
        source, OutputConfig(tmp_path / "snapshot"), solver="two_term"
    )
    artifact = read_reaction_rates(paths["manifest_json"])
    assert artifact.metadata["source"]["manifest_content_sha256"] == content_digest(
        original.metadata
    )
    assert artifact.metadata["source"]["manifest_content_sha256"] != content_digest(
        json.loads(source.read_text(encoding="utf-8"))
    )


def test_manifest_uses_the_same_base_name_lock_for_all_valid_names(sample, tmp_path):
    artifact = build_reaction_rates([sample], solver=sample.solver, source={})
    paths = write_reaction_rates(
        artifact, OutputConfig(tmp_path, "demo_reaction_rates")
    )
    assert paths["manifest_json"].name == "demo_reaction_rates_manifest.json"
    assert not read_reaction_rates(paths["manifest_json"]).table.empty
    assert (tmp_path / ".demo_reaction_rates.lock").is_file()
    assert not (tmp_path / ".demo.lock").exists()
