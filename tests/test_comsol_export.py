from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest
import yaml

from swarm_workflow.cli import main as workflow_cli_main
from swarm_workflow.comsol_export import ComsolExportError, export_comsol_bundle
from swarm_workflow.store import WorkflowStore
from swarm_workflow.sweep import run_sweep
from swarm_workflow.tables import build_tables
from product_helpers import (
    base_product_config,
    insert_workflow_case,
    insert_workflow_eedf,
    insert_workflow_rate,
    write_config,
)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fp:
        return list(csv.DictReader(fp))


def _write_two_mixture_workflow(tmp_path: Path) -> Path:
    data = base_product_config(tmp_path, ["two_term"])
    base_path = write_config(tmp_path, data, "base.yaml")
    workflow_path = tmp_path / "workflow.yaml"
    workflow_path.write_text(
        yaml.safe_dump(
            {
                "base_config": base_path.name,
                "database": "outputs/swarm.sqlite",
                "e_over_n_Td": [10.0, 20.0],
                "mixtures": [{"Ar": 1.0}, {"Ar": 1.0}],
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    return workflow_path


def test_two_mixture_two_term_sweep_build_tables_and_export_comsol(
    tmp_path: Path,
) -> None:
    workflow_path = _write_two_mixture_workflow(tmp_path)
    run_sweep(workflow_path)

    workflow_cli_main(
        [
            "build-tables",
            str(tmp_path / "outputs" / "swarm.sqlite"),
            "--output",
            str(tmp_path / "tables"),
            "--source",
            "two_term",
        ]
    )
    workflow_cli_main(
        [
            "export-comsol",
            str(tmp_path / "tables"),
            "--output",
            str(tmp_path / "bundle"),
        ]
    )

    root_manifest = json.loads((tmp_path / "bundle" / "manifest.json").read_text())
    assert len(root_manifest["bundles"]) == 2
    for mixture_id in [0, 1]:
        bundle = tmp_path / "bundle" / f"mixture_{mixture_id:04d}"
        expected = {
            "manifest.json",
            "mean_energy_vs_en.csv",
            "transport_vs_en.csv",
            "rates_vs_en.csv",
            "townsend_vs_en.csv",
            "quality.csv",
            "eedf.csv",
            "eedf_f0.csv",
        }
        assert expected.issubset({path.name for path in bundle.iterdir()})
        function_dir = bundle / "comsol_functions"
        assert {
            "sw_meanE.csv",
            "sw_muN.csv",
        } == {path.name for path in function_dir.iterdir()}
        manifest = json.loads((bundle / "manifest.json").read_text())
        assert manifest["status"] == "ok"
        assert manifest["tables"]["comsol_functions/sw_meanE.csv"][
            "derived_for_comsol_interpolation"
        ]
        assert manifest["tables"]["eedf_f0.csv"]["definition"] == (
            "eepf_eV_m32 = eedf / sqrt(electron_energy_eV)"
        )
        transport = _read_csv(bundle / "transport_vs_en.csv")
        assert len(transport) == 2
        assert float(transport[0]["E_over_N_V_m2"]) == pytest.approx(
            float(transport[0]["E_over_N_Td"]) * 1.0e-21
        )


def test_export_rates_do_not_double_apply_mixture_fraction(tmp_path: Path) -> None:
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        insert_workflow_case(store.connection, solver="two_term", drift=2.0)
        insert_workflow_rate(
            store.connection,
            solver="two_term",
            value=4.0e-15,
            target_fraction=0.25,
            process="excitation",
            process_type="excitation",
            threshold_eV=11.5,
        )
        insert_workflow_eedf(store.connection, solver="two_term")
        store.connection.commit()

    build_tables(tmp_path / "workflow.sqlite", tmp_path / "tables", source="two_term")
    export_comsol_bundle(
        tmp_path / "tables" / "mixture_0000",
        tmp_path / "bundle",
    )
    row = _read_csv(tmp_path / "bundle" / "rates_vs_en.csv")[0]

    assert float(row["rate_coefficient_m3_s"]) == pytest.approx(4.0e-15)
    assert float(row["mixture_weighted_rate_m3_s"]) == pytest.approx(1.0e-15)
    assert float(row["reduced_townsend_m2"]) == pytest.approx(2.0e-15)
    assert float(row["mixture_weighted_reduced_townsend_m2"]) == pytest.approx(0.5e-15)
    assert (tmp_path / "bundle" / "rates_vs_mean_energy.csv").exists()


def test_export_omits_unused_rate_and_energy_loss_functions(tmp_path: Path) -> None:
    table_dir = tmp_path / "tables" / "mixture_0000"
    table_dir.mkdir(parents=True)
    _write_minimal_manifest(table_dir)
    _write_required_base_tables(table_dir)

    export_comsol_bundle(table_dir, tmp_path / "bundle")

    function_dir = tmp_path / "bundle" / "comsol_functions"
    assert {path.name for path in function_dir.iterdir()} == {
        "sw_meanE.csv",
        "sw_muN.csv",
    }
    assert not (tmp_path / "bundle" / "energy_loss.csv").exists()


def test_export_rejects_missing_mobility_with_failure_manifest(
    tmp_path: Path,
) -> None:
    table_dir = tmp_path / "tables" / "mixture_0000"
    table_dir.mkdir(parents=True)
    _write_minimal_manifest(table_dir)
    _write_required_base_tables(table_dir)
    transport = table_dir / "transport_vs_en.csv"
    transport.write_text(
        transport.read_text(encoding="utf-8").replace(
            "10.0,1e-20,2.0,5.0,10.0,20.0",
            "10.0,1e-20,2.0,5.0,,20.0",
        ),
        encoding="utf-8",
    )

    with pytest.raises(ComsolExportError, match="missing required coefficients"):
        export_comsol_bundle(
            table_dir,
            tmp_path / "bundle",
        )

    manifest = json.loads((tmp_path / "bundle" / "manifest.json").read_text())
    assert manifest["status"] == "failed"
    assert any(
        "reduced_mobility_m2_V_s_m3" in item
        for item in manifest["missing_required_coefficients"]
    )


def _write_minimal_manifest(table_dir: Path) -> None:
    (table_dir / "manifest.json").write_text(
        json.dumps(
            {
                "format_version": 1,
                "stage": "build-tables",
                "source": "two_term",
                "hashes": {},
                "mixture": {"mixture_id": 0, "species": []},
                "source_policy": {"source": "two_term"},
                "valid_ranges": {"E_over_N_Td": [10.0, 10.0]},
                "units": {},
                "table_argument": {"primary": "E_over_N_Td", "secondary": None},
                "monotonicity": {"mean_energy_strictly_monotonic": True},
                "quality_summary": {"passed": True, "failed_points": 0, "total_points": 1},
                "tables": {},
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def _write_required_base_tables(table_dir: Path) -> None:
    (table_dir / "mean_energy_vs_en.csv").write_text(
        "E_over_N_Td,E_over_N_V_m2,mean_energy_eV\n10.0,1e-20,2.0\n",
        encoding="utf-8",
    )
    (table_dir / "transport_vs_en.csv").write_text(
        (
            "E_over_N_Td,E_over_N_V_m2,mean_energy_eV,drift_velocity_m_s,"
            "reduced_mobility_m2_V_s_m3,reduced_diffusion_L_m2_s_m3,"
            "reduced_diffusion_T_m2_s_m3,"
            "reduced_electron_energy_mobility_m2_V_s_m3,"
            "reduced_electron_energy_diffusion_m2_s_m3\n"
            "10.0,1e-20,2.0,5.0,10.0,20.0,30.0,7.0,8.0\n"
        ),
        encoding="utf-8",
    )
    (table_dir / "rates_vs_en.csv").write_text(
        (
            "E_over_N_Td,E_over_N_V_m2,species,process,process_type,threshold_eV,"
            "target_species_fraction,rate_coefficient_m3_s,"
            "mixture_weighted_rate_m3_s,reduced_townsend_m2,"
            "mixture_weighted_reduced_townsend_m2\n"
            "10.0,1e-20,Ar,ionization,ionization,15.0,1.0,1e-15,1e-15,2e-16,2e-16\n"
        ),
        encoding="utf-8",
    )
    (table_dir / "townsend_vs_en.csv").write_text(
        "E_over_N_Td,E_over_N_V_m2,effective_townsend_m2\n10.0,1e-20,2e-16\n",
        encoding="utf-8",
    )
    (table_dir / "energy_loss.csv").write_text(
        (
            "E_over_N_Td,E_over_N_V_m2,species,process,process_type,threshold_eV,"
            "target_species_fraction,energy_loss_eV,energy_loss_rate_coefficient_eV_m3_s\n"
            "10.0,1e-20,Ar,ionization,ionization,15.0,1.0,15.0,1.5e-14\n"
        ),
        encoding="utf-8",
    )
    (table_dir / "quality.csv").write_text(
        (
            "E_over_N_Td,E_over_N_V_m2,passed,failure_reasons_json,mobility_rse,"
            "diffusion_L_rse,diffusion_T_rse,max_major_rate_rse,"
            "eedf_normalization_error,valid_replicates,uncertainty_available,"
            "quality_source\n"
            "10.0,1e-20,1,[],0,0,0,0,0,1,0,two_term\n"
        ),
        encoding="utf-8",
    )
