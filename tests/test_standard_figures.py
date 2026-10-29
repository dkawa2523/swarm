"""Figures retain physical meaning, numerical provenance and native source values."""

from dataclasses import replace
from pathlib import Path
import json
import sqlite3

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("matplotlib")

from electron_swarm import load_config, run
from electron_swarm.core.output import OutputConfig
from electron_swarm.io import (
    read_figures,
    read_product_run,
    read_reaction_rates,
    write_outputs,
)
from electron_swarm.io.figure_data import FigureCase, eedf_curves, scalar_points
from electron_swarm.orchestration.preparation import prepare_run
from swarm_workflow.campaign.document import calculation_sha256
from swarm_workflow.campaign.sweep import plan_sweep, run_sweep
from swarm_workflow.cli import main
from swarm_workflow.figures import plot_results
from swarm_workflow.reaction_rates import export_reaction_rates
from product_helpers import base_product_config, write_config


def _input(tmp_path, solver="two_term"):
    data = base_product_config(tmp_path, [solver])
    data["run"]["e_over_n_Td"] = [10.0]
    data["output"]["figures"] = ["eedf", "transport", "reaction_rates"]
    if solver == "propagator":
        data["solvers"][solver] = {
            "energy_cells": 64,
            "polar_cells": 8,
            "max_memory_mb": 128,
        }
        data["output"]["figures"].append("energy_angle")
    if solver == "monte_carlo":
        data["solvers"][solver] = {
            "seed": 17,
            "particles": 2,
            "production_trial_periods": 8,
        }
        data["physics"]["ionization"]["population_model"] = (
            "fixed_particle_single_daughter"
        )
    return data


@pytest.mark.parametrize(
    "figures", [None, "eedf", ["eedf", "eedf"], ["composition_time"], [{}]]
)
def test_invalid_figure_requests_fail_at_loading(tmp_path, figures):
    data = _input(tmp_path)
    data["output"]["figures"] = figures
    with pytest.raises(ValueError, match="output.figures"):
        load_config(write_config(tmp_path, data))


@pytest.mark.parametrize("kind", ["energy_angle", "rf_phase"])
def test_unsupported_figure_plan_fails_before_execution(tmp_path, kind):
    data = _input(tmp_path)
    data["output"].update(directory=str(tmp_path / "absent"), figures=[kind])
    path = write_config(tmp_path, data)
    with pytest.raises(ValueError):
        prepare_run(load_config(path))
    data["campaign"] = {"database": "absent/cases.sqlite"}
    with pytest.raises(ValueError):
        plan_sweep(write_config(tmp_path, data))
    assert not (tmp_path / "absent").exists()


@pytest.mark.parametrize(
    "solver", ["two_term", "multi_term", "monte_carlo", "propagator"]
)
def test_direct_figures_use_native_values_and_explicit_units(tmp_path, solver):
    data = _input(tmp_path, solver)
    path = write_config(tmp_path, data)
    identity = calculation_sha256(path)
    config = load_config(path)
    assert (
        prepare_run(config).description()["figures"]["requested"]
        == data["output"]["figures"]
    )
    assert not (tmp_path / "prod_manifest.json").exists()
    result = run(config)
    artifact = read_figures(result.metadata["output_paths"]["manifest_json"])
    assert artifact.metadata["evidence_scope"] == "solver_execution"
    points = artifact.tables["figure_data_csv"]
    selected = points.loc[(points.figure_id == "eedf_0000") & (points.panel == "eedf")]
    np.testing.assert_array_equal(selected.x, result.cases[0].energy_eV)
    np.testing.assert_array_equal(selected.y, result.cases[0].eedf)
    np.testing.assert_allclose(
        selected.x_upper - selected.x_lower, result.cases[0].energy_widths_eV
    )
    assert selected.ci95_low.isna().all()
    assert (
        artifact.metadata["figures"]["eedf_0000"]["panels"]["eedf"]["y_unit"] == "1/eV"
    )
    rates = points.loc[points.figure_id == "reaction_rates_0000"]
    np.testing.assert_array_equal(
        sorted(rates.y),
        sorted(rate.rate_coefficient_m3_s for rate in result.cases[0].rates),
    )
    for spec in artifact.metadata["figures"].values():
        assert set(spec["files"]) == {"png", "svg"}
        assert (
            (tmp_path / spec["files"]["png"]["name"])
            .read_bytes()
            .startswith(b"\x89PNG")
        )
    data["output"]["figures"] = []
    assert calculation_sha256(write_config(tmp_path, data)) == identity
    if solver == "propagator":
        spec = artifact.metadata["figures"]["energy_angle_0000"]
        assert (
            spec["panels"]["density"]["data_table"] == "energy_angle_distribution_csv"
        )
        assert spec["panels"]["density"]["numerical_status"] == "numerically_converged"
        assert any(
            "diffusion" in value.lower()
            for value in artifact.metadata["figures"]["transport_0000"]["unavailable"]
        )


def test_campaign_figures_and_export_share_snapshot_and_resume(tmp_path):
    data = _input(tmp_path, "monte_carlo")
    data["physics"]["ionization"]["population_model"] = "weighted_branching"
    data["output"]["export"] = {
        "solver": "monte_carlo",
        "quantities": ["reaction_rates"],
    }
    data["run"]["solvers"].append({"id": "two_term"})
    data["campaign"] = {"database": "cases.sqlite", "mc": {"replicas": 2}}
    path = write_config(tmp_path, data)
    assert plan_sweep(path)["figures"]["sources"]["eedf"] == ["monte_carlo", "two_term"]
    assert run_sweep(path).cases_written == 3
    artifact = read_figures(tmp_path / "prod_manifest.json")
    assert artifact.metadata["evidence_scope"] == "saved_results_analysis"
    assert set(read_reaction_rates(tmp_path / "prod_manifest.json").table.solver) == {
        "monte_carlo"
    }
    spec = artifact.metadata["figures"]["reaction_rates_0000"]
    assert {curve["solver"] for curve in spec["curves"].values()} == {
        "monte_carlo",
        "two_term",
    }
    spec = artifact.metadata["figures"]["eedf_0000"]
    [mc_id] = [
        key
        for key, curve in spec["curves"].items()
        if curve["solver"] == "monte_carlo"
        and key
        in set(
            artifact.tables["figure_data_csv"].loc[lambda x: x.panel == "eedf", "curve"]
        )
    ]
    points = artifact.tables["figure_data_csv"]
    mc = points.loc[(points.figure_id == "eedf_0000") & (points.curve == mc_id)]
    assert set(mc.replicates) == {2}
    assert mc.ci95_high.notna().any()
    source_before = (tmp_path / "cases.sqlite").read_bytes()
    plot_results(
        tmp_path / "cases.sqlite", OutputConfig(tmp_path / "plots", figures=("eedf",))
    )
    assert (tmp_path / "cases.sqlite").read_bytes() == source_before
    assert run_sweep(path).cases_written == 0


def test_replot_repairs_images_preserves_numeric_execution_and_cleans_stale_files(
    tmp_path, monkeypatch
):
    import electron_swarm.io.publication as publication

    data = _input(tmp_path)
    data["output"]["base_name"] = "prod[1]"
    config = load_config(write_config(tmp_path, data))
    result = run(config)
    manifest = Path(result.metadata["output_paths"]["manifest_json"])
    old = read_product_run(manifest)
    unrelated = tmp_path / "prod[1]_figure_eedf_figure_eedf_0000.png"
    unrelated.write_bytes(b"another base name's figure")
    image_path = Path(result.metadata["output_paths"]["figure_eedf_0000_png"])
    image_path.unlink()
    with pytest.raises(OSError):
        read_figures(manifest)
    plot_results(manifest, OutputConfig(tmp_path, "prod[1]", figures=("eedf",)))
    repaired = read_figures(manifest)
    assert repaired.metadata["execution"] == old.metadata["execution"]
    pd.testing.assert_frame_equal(repaired.tables["eedf_csv"], old.tables["eedf_csv"])
    assert image_path.exists()
    assert unrelated.read_bytes() == b"another base name's figure"
    assert not list(tmp_path.glob("prod[[]1]_figure_transport*"))
    before = {
        path: path.read_bytes()
        for path in tmp_path.iterdir()
        if path.name.startswith("prod[1]_")
    }
    original_replace = publication.os.replace
    calls = 0

    def fail_second(source, target):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise PermissionError("simulated locked output")
        return original_replace(source, target)

    with monkeypatch.context() as patch:
        patch.setattr(publication.os, "replace", fail_second)
        with pytest.raises(PermissionError):
            plot_results(
                manifest, OutputConfig(tmp_path, "prod[1]", figures=("transport",))
            )
    assert all(path.read_bytes() == content for path, content in before.items())
    config.output.figures = ()
    config.output.directory = str(tmp_path)
    write_outputs(result, config.output)
    assert not any(
        "_figure_" in path.name and path != unrelated for path in tmp_path.iterdir()
    )
    assert read_product_run(manifest).metadata.get("figures") is None


def test_neutral_coefficients_plot_and_missing_eedf_is_explicit(tmp_path, capsys):
    data = _input(tmp_path)
    data["output"]["figures"] = []
    run(load_config(write_config(tmp_path, data)))
    source = export_reaction_rates(
        tmp_path / "prod_manifest.json",
        OutputConfig(tmp_path / "neutral"),
        solver="two_term",
    )["manifest_json"]
    main(
        [
            "plot",
            str(source),
            "--figures",
            "reaction_rates",
            "--output",
            str(tmp_path / "plots"),
        ]
    )
    assert "swarm_manifest.json" in capsys.readouterr().out
    assert (
        read_reaction_rates(tmp_path / "plots/swarm_manifest.json").metadata["solver"]
        == "two_term"
    )
    before = source.read_bytes()
    with pytest.raises(ValueError, match="nonempty eedf_csv"):
        plot_results(source, OutputConfig(tmp_path / "bad", figures=("eedf",)))
    assert source.read_bytes() == before
    assert not (tmp_path / "bad/swarm_manifest.json").exists()


def test_transport_meanings_and_stale_metadata_are_visible(tmp_path):
    data = _input(tmp_path)
    data["run"]["solvers"].append({"id": "multi_term"})
    data["physics"]["electron_electron"] = {
        "enabled": True,
        "model": "relaxation_postprocess",
        "relaxation_fraction": 0.05,
    }
    config = load_config(write_config(tmp_path, data))
    result = run(config, write=False)
    # Distinct transport definitions must never become a single comparison curve.
    case = result.cases[1]
    definition = case.transport_observables["diffusion_L_m2_s"]
    case.transport_observables["diffusion_L_m2_s"] = replace(
        definition, quantity="restricted_density_packet_diffusion_L"
    )
    paths = write_outputs(result, config.output)
    artifact = read_figures(paths["manifest_json"])
    specs = [
        spec
        for spec in artifact.metadata["figures"].values()
        if spec["kind"] == "transport"
    ]
    meanings = {
        panel["quantity"] for spec in specs for panel in spec["panels"].values()
    }
    assert "restricted_density_packet_diffusion_L" in meanings
    assert "isotropic_particle_diffusion_approximation" in meanings
    assert any(
        "transport stale" in curve["label"]
        for spec in specs
        for curve in spec["curves"].values()
    )
    assert any(
        "stale" in value
        for value in artifact.tables["figure_data_csv"].numerical_status
    )


def test_replica_distribution_uses_conservative_overlap_without_renormalization():
    context = {"angular_scattering": {"model": "isotropic"}}
    cases = [
        FigureCase(
            "monte_carlo", str(i), 10, context, {}, {"monte_carlo_case_seed": 11 + i}, i
        )
        for i in range(2)
    ]
    frame = pd.DataFrame(
        [
            {
                "solver": "monte_carlo",
                "case_id": "0",
                "E_over_N_Td": 10,
                "energy_eV": e,
                "energy_width_eV": 1,
                "eedf": v,
            }
            for e, v in [(0.5, 0.5), (1.5, 0.5)]
        ]
        + [
            {
                "solver": "monte_carlo",
                "case_id": "1",
                "E_over_N_Td": 10,
                "energy_eV": 1.5,
                "energy_width_eV": 3,
                "eedf": 1 / 3,
            }
        ]
    )
    [(group, _, widths, values, low, high, info)] = eedf_curves(cases, frame)
    assert len(group) == 2
    assert np.sum(widths * values) == pytest.approx(5 / 6)
    assert info["source_mass_outside_common_domain"] == pytest.approx([0, 1 / 3])
    assert all(value is not None for value in low + high)
    unproved = [replace(case, replicate=None) for case in cases]
    assert len(eedf_curves(unproved, frame)) == 2
    points = scalar_points(
        unproved, {tuple(case.key): float(i) for i, case in enumerate(unproved)}
    )
    assert len(points) == 2
    assert all(point["ci95_high"] is None for point in points)


def test_rf_phase_is_electron_response_and_dc_rate_plot_is_rejected(tmp_path):
    data = _input(tmp_path)
    data["physics"]["energy_grid_policy"]["adaptive"] = False
    data["physics"]["field"].update(
        type="time_dependent",
        time_dependent={
            "frequency_Hz": 13.56e6,
            "phase_steps": 8,
            "max_periods": 4,
            "periodic_tolerance": 0.5,
        },
    )
    with pytest.raises(NotImplementedError, match="DC"):
        prepare_run(load_config(write_config(tmp_path, data)))
    data["output"]["figures"] = ["rf_phase", "eedf", "transport"]
    result = run(load_config(write_config(tmp_path, data)))
    artifact = read_figures(result.metadata["output_paths"]["manifest_json"])
    spec = artifact.metadata["figures"]["rf_phase_0000"]
    assert spec["physical_context"]["field_frequency_Hz"] == 13.56e6
    assert any("not chemical-species time evolution" in text for text in spec["notes"])
    points = artifact.tables["figure_data_csv"].loc[
        lambda x: (x.figure_id == "rf_phase_0000") & (x.panel == "mean_energy_eV")
    ]
    np.testing.assert_array_equal(points.y, result.cases[0].rf_phase.mean_energy_eV)
    assert len(points) == 8


def test_historical_database_plot_retains_original_identity_and_checks_numeric_hashes(
    tmp_path,
):
    data = _input(tmp_path)
    data["output"]["figures"] = []
    data["campaign"] = {"database": "cases.sqlite"}
    database = run_sweep(write_config(tmp_path, data)).database_path
    with sqlite3.connect(database) as connection:
        original = json.loads(
            connection.execute(
                "SELECT value FROM metadata WHERE key='solver_implementations_json'"
            ).fetchone()[0]
        )
        original["two_term"]["execution_sha256"] = "0" * 64
        connection.execute(
            "UPDATE metadata SET value=? WHERE key='solver_implementations_json'",
            (json.dumps(original),),
        )
    before = database.read_bytes()
    output = OutputConfig(tmp_path / "analysis", figures=("eedf",))
    paths = plot_results(database, output)
    artifact = read_figures(paths["manifest_json"])
    assert artifact.metadata["source"]["implementations"] == original
    assert "execution" not in artifact.metadata
    assert database.read_bytes() == before
    altered = json.loads(paths["manifest_json"].read_text())
    altered["files"]["../../escape_csv"] = altered["files"]["eedf_csv"]
    unsafe = tmp_path / "analysis/unsafe_manifest.json"
    unsafe.write_text(json.dumps(altered))
    with pytest.raises(ValueError, match="canonical output roles"):
        plot_results(unsafe, OutputConfig(tmp_path / "guarded", figures=("eedf",)))
    assert not (tmp_path / "guarded").exists()
    paths["eedf_csv"].write_bytes(b"altered")
    with pytest.raises(ValueError, match="altered"):
        plot_results(paths["manifest_json"], output)


def test_signed_native_tail_is_not_clipped_and_image_loss_does_not_block_rate_export(
    tmp_path,
):
    data = _input(tmp_path)
    data["output"]["figures"] = ["eedf"]
    config = load_config(write_config(tmp_path, data))
    result = run(config, write=False)
    case = result.cases[0]
    previous = case.eedf[-1] * case.energy_widths_eV[-1]
    case.eedf[-1] = -1e-12 / case.energy_widths_eV[-1]
    case.eedf[0] += (previous + 1e-12) / case.energy_widths_eV[0]
    paths = write_outputs(result, config.output)
    artifact = read_figures(paths["manifest_json"])
    points = artifact.tables["figure_data_csv"].loc[lambda x: x.panel == "eedf"]
    np.testing.assert_array_equal(points.y, case.eedf)
    assert artifact.metadata["figures"]["eedf_0000"]["panels"]["eedf"][
        "scale_y"
    ].startswith("symlog")
    paths["figure_eedf_0000_png"].unlink()
    exported = export_reaction_rates(
        paths["manifest_json"], OutputConfig(tmp_path / "rates"), solver="two_term"
    )
    assert set(
        read_reaction_rates(exported["manifest_json"]).table.numerical_status
    ) == {"derived"}
