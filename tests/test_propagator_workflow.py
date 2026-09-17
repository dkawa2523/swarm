from __future__ import annotations

import json
from pathlib import Path
import sqlite3

import pandas as pd
import pytest
import yaml

from swarm_workflow.comsol.input import export_comsol_bundle
from swarm_workflow.comsol.models.gec_ccp import (
    GecCcpWorkflowError,
    prepare_gec_ccp_run,
)
from swarm_workflow.campaign.sweep import run_sweep
from swarm_workflow.quality.propagator_source import (
    PROPAGATOR_SOLVER_SOURCE_METADATA_KEY,
)
from swarm_workflow.tables.contracts import TableBuildError
from swarm_workflow.tables import build_tables

from product_helpers import ROOT, base_product_config, write_config


def test_propagator_workflow_persists_null_transport_and_exports_bundle(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["propagator"])
    data["run"]["e_over_n_Td"] = [10.0]
    data["solvers"]["propagator"] = {
        "energy_cells": 64,
        "polar_cells": 8,
        "max_iterations": 2000,
        "convergence_tolerance": 1.0e-8,
        "max_memory_mb": 128,
    }
    base = write_config(tmp_path, data, "propagator_base.yaml")
    workflow = tmp_path / "propagator_workflow.yaml"
    workflow.write_text(
        yaml.safe_dump(
            {
                "base_config": base.name,
                "database": "propagator.sqlite",
                "e_over_n_Td": [10.0, 30.0, 100.0],
                "mixtures": [{"Ar": 1.0}],
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    first = run_sweep(workflow)
    second = run_sweep(workflow)
    assert first.cases_written == 3
    assert second.cases_written == 0
    with sqlite3.connect(first.database_path) as connection:
        rows = connection.execute(
            "SELECT solver, diffusion_L_m2_s, diffusion_T_m2_s "
            "FROM cases ORDER BY e_over_n_Td"
        ).fetchall()
        assert rows == [("propagator", None, None)] * 3
        angle_counts = connection.execute(
            "SELECT e_over_n_Td, COUNT(*) FROM energy_angle_bins "
            "GROUP BY e_over_n_Td ORDER BY e_over_n_Td"
        ).fetchall()
        assert [item[0] for item in angle_counts] == [10.0, 30.0, 100.0]
        assert len({item[1] for item in angle_counts}) == 1
        assert angle_counts[0][1] > 64 * 8
        assert angle_counts[0][1] % 8 == 0
        assert connection.execute(
            "SELECT COUNT(*) FROM aggregate_quality WHERE passed = 1"
        ).fetchone() == (3,)
        propagator_source = connection.execute(
            "SELECT value FROM metadata WHERE key = ?",
            (PROPAGATOR_SOLVER_SOURCE_METADATA_KEY,),
        ).fetchone()[0]
        assert len(propagator_source) == 64

    tables_dir = tmp_path / "tables"
    with pytest.raises(
        TableBuildError,
        match="require --solver-qualification",
    ):
        build_tables(
            first.database_path,
            tables_dir,
            source="propagator",
        )
    table_summary = build_tables(
        first.database_path,
        tables_dir,
        source="propagator",
        solver_qualification_path=(
            ROOT
            / "docs"
            / "dev"
            / "results"
            / "propagator_p1_deterministic_qualification_20260908.json"
        ),
    )
    assert table_summary.unqualified_mixtures == 0
    mixture_tables = tables_dir / "mixture_0000"
    table_manifest = json.loads(
        (mixture_tables / "manifest.json").read_text(encoding="utf-8")
    )
    assert table_manifest["source"] == "propagator"
    assert table_manifest["hashes"][PROPAGATOR_SOLVER_SOURCE_METADATA_KEY] == (
        propagator_source
    )
    transport = pd.read_csv(mixture_tables / "transport_vs_mean_energy.csv")
    assert transport["reduced_mobility_m2_V_s_m3"].notna().all()
    assert transport["reduced_diffusion_L_m2_s_m3"].isna().all()
    assert transport["reduced_diffusion_T_m2_s_m3"].isna().all()

    with sqlite3.connect(first.database_path) as connection:
        connection.execute(
            "UPDATE metadata SET value = ? WHERE key = ?",
            ("0" * 64, PROPAGATOR_SOLVER_SOURCE_METADATA_KEY),
        )
        connection.commit()
    with pytest.raises(TableBuildError, match="solver source does not match"):
        build_tables(
            first.database_path,
            tmp_path / "stale_tables",
            source="propagator",
            solver_qualification_path=(
                ROOT
                / "docs"
                / "dev"
                / "results"
                / "propagator_p1_deterministic_qualification_20260908.json"
            ),
        )

    bundle_dir = tmp_path / "bundle"
    export_summary = export_comsol_bundle(tables_dir, bundle_dir)
    assert export_summary.bundles_written == 1
    mixture_bundle = bundle_dir / "mixture_0000"
    bundle_manifest = json.loads(
        (mixture_bundle / "manifest.json").read_text(encoding="utf-8")
    )
    assert bundle_manifest["source"] == "propagator"
    assert bundle_manifest["hashes"][PROPAGATOR_SOLVER_SOURCE_METADATA_KEY] == (
        propagator_source
    )
    assert (mixture_bundle / "eedf_f0_comsol_2d.csv").is_file()
    assert not (mixture_bundle / "eedf_f0_comsol_grid.txt").exists()
    active_eedf = bundle_manifest["tables"]["eedf_f0_comsol_2d.csv"]
    assert active_eedf["format"] == "csv"
    assert active_eedf["mean_energy_grid_points"] >= 3
    assert active_eedf["solver_anchor_mean_energy_count"] == 3
    assert active_eedf["mean_energy_axis"] == (
        "error_controlled_c1_materialization_with_solver_anchors"
    )
    assert active_eedf["representation"] == (
        "physical_2d_adaptive_moment_rate_projected"
    )
    assert (mixture_bundle / "elastic_energy_loss_vs_mean_energy.csv").is_file()

    mapping_data = yaml.safe_load(
        (
            ROOT
            / "comsol_modes"
            / "maps"
            / "argon_gec_ccp_propagator_function_eedf.yaml"
        ).read_text(encoding="utf-8")
    )
    mapping_data["model"]["input_mph"] = (
        ROOT / "comsol_modes" / "argon_gec_ccp.mph"
    ).as_posix()
    mapping_data["model"]["external_output_mph"] = (
        tmp_path / "work" / "propagator.mph"
    ).as_posix()
    mapping_data["results"]["output_directory"] = (
        tmp_path / "comsol_results"
    ).as_posix()
    mapping_data["logs"]["path"] = (tmp_path / "comsol_logs").as_posix()
    mapping_path = tmp_path / "propagator_map.yaml"
    mapping_path.write_text(
        yaml.safe_dump(mapping_data, sort_keys=False), encoding="utf-8"
    )
    with pytest.raises(
        GecCcpWorkflowError,
        match="lacks bound target-range qualification",
    ):
        prepare_gec_ccp_run(
            mapping_path,
            bundle_path=mixture_bundle,
            write_java=False,
        )
