from __future__ import annotations

import csv
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from swarm_workflow.cli import main as workflow_cli_main
from swarm_workflow.comsol_positive_column import (
    PositiveColumnWorkflowError,
    execute_positive_column_run,
    generate_result_export_java,
    prepare_positive_column_run,
    validate_positive_column_result,
    validate_positive_column_bundle,
)


def test_prepare_positive_column_dry_run_writes_effective_mapping_and_java(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_positive_column_repo(tmp_path)

    plan = prepare_positive_column_run(mapping_path, bundle_path=bundle)

    assert plan.effective_mapping_path.exists()
    assert plan.apply_java_path.exists()
    assert plan.verify_java_path.exists()
    assert plan.run_java_path.exists()
    assert plan.result_java_path.exists()
    assert plan.result_csv_path.name == "positive_column_results.csv"
    assert plan.mapping.model.input_mph != plan.mapping.model.output_mph
    assert plan.bundle.valid_e_over_n_Td == (10.0, 100.0)
    assert "apply COMSOL interpolation tables" in plan.steps
    assert "verify imported COMSOL interpolation functions" in plan.steps
    run_source = plan.run_java_path.read_text(encoding="utf-8")
    assert '.feature("mct1").set("V0", "V0")' in run_source
    assert run_source.count('.study("std1").run()') == 1
    assert 'runAtVoltage(model, "200[V]")' in run_source
    assert "Direct target-voltage solve failed" in run_source
    assert (
        "private static void runAtVoltage(Model model, String voltage) "
        "throws Exception"
    ) in run_source
    assert 'System.out.println("SWARM_COMSOL_VOLTAGE_V=" + voltage)' in run_source


def test_prepare_positive_column_applies_explicit_conditions_with_readback(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_positive_column_repo(tmp_path)

    plan = prepare_positive_column_run(
        mapping_path,
        bundle_path=bundle,
        pressure_Pa=13.3322,
        gas_temperature_K=293.15,
        mesh_elements=400,
    )

    source = plan.apply_java_path.read_text(encoding="utf-8")
    assert 'model.param().set("p0", "13.3322[Pa]")' in source
    assert '.feature("pes1").set("T_src", "userdef")' in source
    assert '.feature("pes1").set("T", "293.14999999999998[K]")' in source
    assert '.feature("dis1").set("elemcount", 400)' in source
    assert '.mesh("mesh1").run()' in source
    assert '.feature("pes1").getString("T")' in source
    assert '.feature("dis1").getString("elemcount")' in source
    assert "SWARM_COMSOL_PRESSURE_PA=13.3322" in source
    assert "SWARM_COMSOL_GAS_TEMPERATURE_K=" in source
    assert "SWARM_COMSOL_MESH_ELEMENTS=" in source
    assert plan.model_conditions.pressure_Pa == pytest.approx(13.3322)
    assert plan.model_conditions.gas_temperature_K == pytest.approx(293.15)
    assert plan.model_conditions.mesh_elements == 400

    with pytest.raises(PositiveColumnWorkflowError, match="200 or 400"):
        prepare_positive_column_run(
            mapping_path,
            bundle_path=bundle,
            mesh_elements=300,
        )


def test_positive_column_result_rejects_table_range_escape(tmp_path: Path) -> None:
    mapping_path, bundle = _write_positive_column_repo(tmp_path)
    plan = prepare_positive_column_run(mapping_path, bundle_path=bundle)
    plan.result_csv_path.write_text(
        "E_over_N,mean_electron_energy,reduced_mobility,reduced_diffusion_L,"
        "excitation_townsend,ionization_townsend\n"
        "2e-19,4,10.6666666666667,20.6666666666667,"
        "1.333333333333e-21,1.333333333333e-22\n",
        encoding="utf-8",
    )

    with pytest.raises(PositiveColumnWorkflowError, match="outside the Swarm table range"):
        validate_positive_column_result(plan)


def test_positive_column_result_rejects_inactive_transport(tmp_path: Path) -> None:
    mapping_path, bundle = _write_positive_column_repo(tmp_path)
    plan = prepare_positive_column_run(mapping_path, bundle_path=bundle)
    plan.result_csv_path.write_text(
        "E_over_N,mean_electron_energy,reduced_mobility,reduced_diffusion_L,"
        "excitation_townsend,ionization_townsend\n"
        "4e-20,4,1,20.6666666666667,"
        "1.333333333333e-21,1.333333333333e-22\n",
        encoding="utf-8",
    )

    with pytest.raises(PositiveColumnWorkflowError, match="did not follow"):
        validate_positive_column_result(plan)


def test_positive_column_result_export_java_uses_mapping_result_expressions(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_positive_column_repo(tmp_path)
    plan = prepare_positive_column_run(mapping_path, bundle_path=bundle)

    source = generate_result_export_java(plan)

    assert "SwarmPositiveColumnResultsExport" in source
    assert "ModelUtil.loadCopy" in source
    assert '"Eval"' in source
    assert "getData()" in source
    assert "getCoordinates()" in source
    assert "positive_column_results.csv" in source
    assert "plas.ne" in source
    assert "plas.EN" in source
    assert "plas.Dexx*plas.Nn" in source
    assert "electron_density" in source
    assert "total_current_density" in source


def test_positive_column_bundle_validation_rejects_bad_quality(tmp_path: Path) -> None:
    _mapping_path, bundle = _write_positive_column_repo(tmp_path, quality_passed=False)

    with pytest.raises(PositiveColumnWorkflowError, match="quality_summary"):
        validate_positive_column_bundle(bundle)


def test_positive_column_rejects_missing_result_probe(tmp_path: Path) -> None:
    mapping_path, bundle = _write_positive_column_repo(tmp_path)
    text = mapping_path.read_text(encoding="utf-8")
    text = text.replace(
        "    - name: E_over_N\n      expression: plas.EN\n      unit: V*m^2\n",
        "",
    )
    mapping_path.write_text(text, encoding="utf-8")

    with pytest.raises(PositiveColumnWorkflowError, match="E_over_N"):
        prepare_positive_column_run(mapping_path, bundle_path=bundle)


def test_positive_column_rejects_template_overwrite(tmp_path: Path) -> None:
    mapping_path, bundle = _write_positive_column_repo(tmp_path)
    text = mapping_path.read_text(encoding="utf-8").replace(
        "output_mph: model/work/positive_column_1d_swarm_tables.mph",
        "output_mph: model/positive_column_1d.mph",
    )
    mapping_path.write_text(text, encoding="utf-8")

    with pytest.raises(PositiveColumnWorkflowError, match="overwrite template"):
        prepare_positive_column_run(mapping_path, bundle_path=bundle)


def test_execute_positive_column_run_invokes_steps_in_order(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_path, bundle = _write_positive_column_repo(tmp_path)
    calls: list[str] = []

    def fake_apply(
        mapping,
        *,
        comsol_executable=None,
        pressure_Pa=None,
        gas_temperature_K=None,
        mesh_elements=None,
    ):
        calls.append("apply")
        assert pressure_Pa == pytest.approx(100.0)
        assert gas_temperature_K == pytest.approx(293.15)
        assert mesh_elements == 200
        stdout = tmp_path / "logs" / "apply" / "apply_stdout.txt"
        stdout.parent.mkdir(parents=True)
        stdout.write_text(
            "SWARM_COMSOL_PRESSURE_PA=100\n"
            "SWARM_COMSOL_GAS_TEMPERATURE_K=293.15[K]\n"
            "SWARM_COMSOL_MESH_ELEMENTS=200\n"
            "SWARM_COMSOL_MEAN_ENERGY_MODEL=LocalEnergyApproximationE\n",
            encoding="utf-8",
        )
        return SimpleNamespace(
            log_dir=stdout.parent,
            stdout_paths=(stdout,),
        )

    def fake_verify(mapping, *, comsol_executable=None):
        calls.append("verify")
        summary_json = (
            tmp_path / "outputs" / "comsol_verify" / "function_verify_summary.json"
        )
        summary_json.parent.mkdir(parents=True, exist_ok=True)
        summary_json.write_text(
            json.dumps(
                {
                    "status": "passed",
                    "closure_items": [
                        {
                            "quantity": quantity,
                            "status": "passed",
                            "control_points": 2,
                            "max_relative_error": 0.0,
                        }
                        for quantity in (
                            "mean_energy",
                            "reduced_mobility",
                            "reduced_longitudinal_diffusion",
                            "reduced_energy_mobility",
                            "reduced_energy_diffusion",
                            "excitation_townsend",
                            "ionization_townsend",
                        )
                    ],
                }
            ),
            encoding="utf-8",
        )
        return SimpleNamespace(
            summary_csv=tmp_path / "outputs" / "comsol_verify" / "summary.csv",
            summary_json=summary_json,
        )

    def fake_export(mapping, java_path, *, operation, comsol_executable=None):
        calls.append(operation)
        if operation == "positive_column_export":
            result = (
                tmp_path
                / "outputs"
                / "comsol_positive_column"
                / "positive_column_results.csv"
            )
            result.parent.mkdir(parents=True, exist_ok=True)
            result.write_text(
                "x,electron_density,E_over_N,mean_electron_energy,"
                "reduced_mobility,reduced_diffusion_L,"
                "excitation_townsend,ionization_townsend,"
                "applied_voltage,gas_pressure\n"
                "0,1,4e-20,4,10.6666666666667,20.6666666666667,"
                "1.333333333333e-21,"
                "1.333333333333e-22,200,100\n",
                encoding="utf-8",
            )
        return SimpleNamespace(log_dir=tmp_path / "logs" / "export")

    monkeypatch.setattr("swarm_workflow.comsol_positive_column.execute_apply_comsol", fake_apply)
    monkeypatch.setattr(
        "swarm_workflow.comsol_positive_column.execute_verify_comsol_functions",
        fake_verify,
    )
    monkeypatch.setattr(
        "swarm_workflow.comsol_positive_column.execute_generated_comsol_java",
        fake_export,
    )

    summary = execute_positive_column_run(
        mapping_path,
        bundle_path=bundle,
        comsol_executable=tmp_path / "bin" / "comsol.exe",
        pressure_Pa=100.0,
        gas_temperature_K=293.15,
        mesh_elements=200,
    )

    assert calls == ["apply", "verify", "positive_column_run", "positive_column_export"]
    assert summary.plan.result_csv_path.exists()
    assert summary.run_summary_json.exists()
    run_summary = json.loads(summary.run_summary_json.read_text(encoding="utf-8"))
    assert run_summary["status"] == "completed"
    assert run_summary["result_stats"]["rows"] == 1
    assert run_summary["result_stats"]["electron_density"]["max"] == 1.0
    assert run_summary["final_voltage_V"] == 200.0
    assert run_summary["executed_voltage_sequence_V"] is None
    assert run_summary["provenance"] == {
        "apply": None,
        "verify": None,
        "run": None,
        "export": None,
    }
    assert run_summary["closure_checks"]["mobility_max_relative_error"] == pytest.approx(0)
    assert run_summary["closure_verification"]["status"] == "passed"
    assert len(run_summary["closure_verification"]["items"]) == 7
    assert run_summary["model_conditions"]["status"] == "passed"
    assert run_summary["model_conditions"]["requested"] == {
        "pressure_Pa": 100.0,
        "gas_temperature_K": 293.15,
        "mesh_elements": 200,
    }
    assert (
        run_summary["model_conditions"]["checks"]["pressure_spatial_profile"][
            "control_points"
        ]
        == 1
    )
    assert run_summary["result_stats"]["gas_pressure"] == {
        "min": 100.0,
        "max": 100.0,
    }
    assert run_summary["note"] == (
        "External closure passed; COMSOL agreement is reported separately."
    )
    assert summary.plan.mapping.model.input_mph.read_bytes() == b"base mph"


def test_cli_run_positive_column_dry_run_prints_all_paths(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    mapping_path, bundle = _write_positive_column_repo(tmp_path)

    workflow_cli_main(
        [
            "run-comsol",
            str(mapping_path),
            "--bundle",
            str(bundle),
            "--dry-run",
        ]
    )

    output = capsys.readouterr().out
    assert "COMSOL positive-column dry-run" in output
    assert "effective_mapping:" in output
    assert "apply_java:" in output
    assert "verify_java:" in output
    assert "run_java:" in output
    assert "result_java:" in output
    assert "result_csv:" in output
    assert "COMSOL command: not executed" in output


def _write_positive_column_repo(
    root: Path,
    *,
    quality_passed: bool = True,
) -> tuple[Path, Path]:
    (root / "pyproject.toml").write_text("[project]\nname='tmp'\n", encoding="utf-8")
    model_dir = root / "model"
    work_dir = model_dir / "work"
    maps_dir = model_dir / "maps"
    bundle = root / "outputs" / "bundle" / "mixture_0000"
    work_dir.mkdir(parents=True)
    maps_dir.mkdir(parents=True)
    bundle.mkdir(parents=True)
    (model_dir / "positive_column_1d.mph").write_bytes(b"base mph")
    _write_bundle(bundle, quality_passed=quality_passed)
    mapping_path = maps_dir / "positive_column.yaml"
    mapping_path.write_text(
        """
model:
  input_mph: model/positive_column_1d.mph
  output_mph: model/work/positive_column_1d_swarm_tables.mph
  study: std1
  component: comp1
  physics: plas
bundle:
  path: outputs/placeholder_bundle
logs:
  path: model/logs
verify:
  output_path: outputs/comsol_verify
  relative_tolerance: 1.0e-6
  absolute_tolerance: 1.0e-12
functions:
  mean_energy_vs_en:
    tag: sw_meanE
    file: mean_energy_vs_en.csv
    column: mean_energy_eV
    nargs: 1
    argunit: V*m^2
    fununit: eV
    interp: linear
    extrap: const
  muN:
    tag: sw_muN
    file: transport_vs_en.csv
    column: reduced_mobility_m2_V_s_m3
    nargs: 1
    argunit: V*m^2
    fununit: 1/(V*m*s)
    interp: linear
    extrap: const
closure:
  feature: pes1
  mean_energy_formulation:
    mode: local_energy
    property_group: ElectronProperties
    property: MeanElectronEnergyModel
    comsol_value: LocalEnergyApproximationE
  mean_energy:
    table: mean_energy_vs_en.csv
reaction_lookups:
  excitation:
    feature: eir2
    form: townsend
    process_type: excitation
  ionization:
    feature: eir4
    form: townsend
    process_type: ionization
run:
  voltage_feature: mct1
  voltages_V: [20, 50, 100, 200]
results:
  probes:
    - name: x
      expression: x
      unit: m
    - name: electron_density
      expression: plas.ne
      unit: 1/m^3
    - name: mean_electron_energy
      expression: plas.meanE
      unit: eV
    - name: electric_potential
      expression: V
      unit: V
    - name: electron_current_density
      expression: plas.Jex
      unit: A/m^2
    - name: ion_current_density
      expression: plas.Jix
      unit: A/m^2
    - name: total_current_density
      expression: plas.Jex+plas.Jix
      unit: A/m^2
    - name: excitation_source
      expression: plas.Rexc
      unit: 1/(m^3*s)
    - name: ionization_source
      expression: plas.Ri
      unit: 1/(m^3*s)
    - name: E_over_N
      expression: plas.EN
      unit: V*m^2
    - name: reduced_mobility
      expression: plas.muexx*plas.Nn
      unit: 1/(V*m*s)
    - name: reduced_diffusion_L
      expression: plas.Dexx*plas.Nn
      unit: 1/(m*s)
    - name: excitation_townsend
      expression: plas.eir2.alpha(plas.ebar)
      unit: m^2
    - name: ionization_townsend
      expression: plas.eir4.alpha(plas.ebar)
      unit: m^2
    - name: applied_voltage
      expression: V0
      unit: V
    - name: gas_pressure
      expression: p0
      unit: Pa
""".lstrip(),
        encoding="utf-8",
    )
    return mapping_path, bundle


def _write_bundle(bundle: Path, *, quality_passed: bool) -> None:
    rows = {
        "mean_energy_vs_en.csv": (
            ["E_over_N_Td", "E_over_N_V_m2", "mean_energy_eV"],
            [["10", "1e-20", "2"], ["100", "1e-19", "8"]],
        ),
        "transport_vs_en.csv": (
            [
                "E_over_N_Td",
                "E_over_N_V_m2",
                "mean_energy_eV",
                "reduced_mobility_m2_V_s_m3",
                "reduced_diffusion_L_m2_s_m3",
                "reduced_diffusion_T_m2_s_m3",
                "reduced_electron_energy_mobility_m2_V_s_m3",
                "reduced_electron_energy_diffusion_m2_s_m3",
            ],
            [
                ["10", "1e-20", "2", "10", "20", "30", "15", "25"],
                ["100", "1e-19", "8", "12", "22", "32", "17", "27"],
            ],
        ),
        "rates_vs_mean_energy.csv": (
            [
                "mean_energy_eV",
                "process_type",
                "reduced_townsend_m2",
            ],
            [
                ["2", "excitation", "1e-21"],
                ["8", "excitation", "2e-21"],
                ["2", "ionization", "1e-22"],
                ["8", "ionization", "2e-22"],
            ],
        ),
        "quality.csv": (
            ["E_over_N_Td", "E_over_N_V_m2", "passed"],
            [["10", "1e-20", "1"]],
        ),
    }
    rows["transport_vs_mean_energy.csv"] = rows["transport_vs_en.csv"]
    tables = {}
    units = {
        "E_over_N_Td": "Td",
        "E_over_N_V_m2": "V m^2",
        "mean_energy_eV": "eV",
        "reduced_mobility_m2_V_s_m3": "m^2/(V s) m^3",
        "reduced_diffusion_L_m2_s_m3": "m^2/s m^3",
        "reduced_diffusion_T_m2_s_m3": "m^2/s m^3",
        "reduced_electron_energy_mobility_m2_V_s_m3": "1/(V m s)",
        "reduced_electron_energy_diffusion_m2_s_m3": "1/(m s)",
        "process_type": "",
        "reduced_townsend_m2": "m^2",
        "passed": "",
    }
    for name, (columns, data) in rows.items():
        with (bundle / name).open("w", encoding="utf-8", newline="") as fp:
            writer = csv.writer(fp)
            writer.writerow(columns)
            writer.writerows(data)
        tables[name] = {
            "columns": columns,
            "units": {column: units[column] for column in columns},
        }
    (bundle / "manifest.json").write_text(
        json.dumps(
            {
                "format_version": 1,
                "stage": "export-comsol",
                "status": "ok",
                "valid_ranges": {"E_over_N_Td": [10.0, 100.0]},
                "quality_summary": {"passed": quality_passed},
                "tables": tables,
                "units": units,
            }
        ),
        encoding="utf-8",
    )
