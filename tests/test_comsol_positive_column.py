from __future__ import annotations
from support.positive_column import result_csv

from swarm_workflow.comsol.models.positive_column.bundle import (
    validate_positive_column_bundle,
)

from swarm_workflow.comsol.models.positive_column.audit import (
    validate_positive_column_result,
)

import csv
import json
from hashlib import sha256
from pathlib import Path
from swarm_workflow.comsol.models.positive_column.results import read_mapped_results
from types import SimpleNamespace

import pytest
import yaml

from support.positive_column import complete_bundle_fixture, complete_mapping_fixture

from swarm_workflow.cli import main as workflow_cli_main
from swarm_workflow.comsol.models.positive_column.workflow import (
    PositiveColumnWorkflowError,
    execute_positive_column_run,
    prepare_positive_column_run,
)
from swarm_workflow.comsol.models.positive_column.config import (
    ComsolMappingError,
    load_comsol_mapping,
)
from swarm_workflow.comsol.models.positive_column.lookup import read_lookup
from swarm_workflow.comsol.models.positive_column.java import generate_apply_java_source
from swarm_workflow.comsol.models.positive_column.verify import (
    plan_verify_comsol_functions,
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
        "private static void runAtVoltage(Model model, String voltage) throws Exception"
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
        result_csv(
            "E_over_N,mean_electron_energy,reduced_mobility,reduced_diffusion_L,"
            "excitation_townsend,ionization_townsend\n"
            "2e-19,4,10.6666666666667,20.6666666666667,"
            "1.333333333333e-21,1.333333333333e-22\n"
        ),
        encoding="utf-8",
    )

    with pytest.raises(
        PositiveColumnWorkflowError, match="outside the Swarm table range"
    ):
        validate_positive_column_result(
            plan, read_mapped_results(plan.result_csv_path, plan.mapping.results)
        )


def test_positive_column_result_rejects_inactive_transport(tmp_path: Path) -> None:
    mapping_path, bundle = _write_positive_column_repo(tmp_path)
    plan = prepare_positive_column_run(mapping_path, bundle_path=bundle)
    plan.result_csv_path.write_text(
        result_csv(
            "E_over_N,mean_electron_energy,reduced_mobility,reduced_diffusion_L,"
            "excitation_townsend,ionization_townsend\n"
            "4e-20,4,1,20.6666666666667,"
            "1.333333333333e-21,1.333333333333e-22\n"
        ),
        encoding="utf-8",
    )

    with pytest.raises(PositiveColumnWorkflowError, match="did not follow"):
        validate_positive_column_result(
            plan, read_mapped_results(plan.result_csv_path, plan.mapping.results)
        )


def test_positive_column_bundle_validation_rejects_bad_quality(tmp_path: Path) -> None:
    _mapping_path, bundle = _write_positive_column_repo(tmp_path, quality_passed=False)

    with pytest.raises(PositiveColumnWorkflowError, match="quality_summary"):
        validate_positive_column_bundle(bundle)


@pytest.mark.parametrize(
    "field,value,reason",
    [
        ("format_version", 1, "format_version"),
        ("stage", "build-tables", "stage"),
        ("status", "unqualified", "status"),
        ("source", "monte_carlo", "source"),
        ("hashes", {}, "provenance hash"),
        ("physical_context", None, "physical_context"),
        (
            "quality_summary",
            {"passed": "false", "failed_points": 0, "total_points": 2},
            "quality_summary",
        ),
    ],
)
def test_positive_column_requires_complete_qualified_manifest(
    tmp_path: Path,
    field: str,
    value: object,
    reason: str,
) -> None:
    mapping, bundle = _write_positive_column_repo(tmp_path)
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest[field] = value
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(PositiveColumnWorkflowError, match=reason):
        prepare_positive_column_run(mapping, bundle_path=bundle)


@pytest.mark.parametrize("table", ["transport_vs_en.csv", "quality.csv"])
def test_positive_column_rejects_changed_table_before_licensed_execution(
    tmp_path: Path,
    table: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping, bundle = _write_positive_column_repo(tmp_path)
    path = bundle / table
    path.write_text(path.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    monkeypatch.setattr(
        "swarm_workflow.comsol.models.positive_column.workflow.execute_apply_comsol",
        lambda *args, **kwargs: pytest.fail("COMSOL must not be invoked"),
    )
    with pytest.raises(PositiveColumnWorkflowError, match="artifact hash mismatch"):
        execute_positive_column_run(mapping, bundle_path=bundle)


def test_positive_column_checks_quality_rows_and_model_mixture(tmp_path: Path) -> None:
    mapping, bundle = _write_positive_column_repo(tmp_path)
    quality = bundle / "quality.csv"
    quality.write_text(
        "E_over_N_Td,E_over_N_V_m2,passed\n10,1e-20,0\n100,1e-19,1\n", encoding="utf-8"
    )
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["tables"]["quality.csv"]["sha256"] = sha256(
        quality.read_bytes()
    ).hexdigest()
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(PositiveColumnWorkflowError, match="quality.csv"):
        prepare_positive_column_run(mapping, bundle_path=bundle)
    manifest["mixture"]["species"][0]["species"] = "O2"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(PositiveColumnWorkflowError, match="gas mixture"):
        prepare_positive_column_run(mapping, bundle_path=bundle)


def test_positive_column_temperature_override_is_used_by_every_preparation_stage(
    tmp_path: Path,
) -> None:
    mapping, bundle = _write_positive_column_repo(tmp_path)
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["physical_context"]["gas_temperature_K"] = 300.0
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(PositiveColumnWorkflowError, match="gas temperature"):
        prepare_positive_column_run(mapping, bundle_path=bundle)
    plan = prepare_positive_column_run(
        mapping, bundle_path=bundle, gas_temperature_K=300.0
    )
    assert plan.mapping.physical_context.gas_temperature_K == 300.0
    assert plan.model_conditions.gas_temperature_K == 300.0
    assert plan_verify_comsol_functions(plan.effective_mapping_path).feature_points
    assert '.set("T", "300[K]")' in plan.apply_java_path.read_text(encoding="utf-8")


def test_positive_column_uses_exact_reaction_identity_for_apply_and_readback(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_positive_column_repo(tmp_path)
    table = bundle / "rates_vs_mean_energy.csv"
    with table.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        columns, rows = list(reader.fieldnames), list(reader)
    alternative = [
        dict(row, process="other Ar excitation", reduced_townsend_m2="9e-99")
        for row in rows
        if row["process_type"] == "excitation"
    ]
    with table.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows + alternative)
    _update_table_hash(bundle, table.name)
    plan = prepare_positive_column_run(mapping_path, bundle_path=bundle)
    excitation = next(
        item
        for item in plan.mapping.reaction_lookups
        if item.process_type == "excitation"
    )
    data = read_lookup(
        excitation.path,
        excitation.x_column,
        excitation.y_column,
        reaction=excitation.identity,
    )
    assert data.y == (1e-21, 2e-21)
    verify = plan_verify_comsol_functions(plan.effective_mapping_path)
    assert [
        point.expected_value
        for point in verify.feature_points
        if point.quantity == "excitation_townsend"
    ] == list(data.y)
    assert "e-99" not in generate_apply_java_source(plan.mapping)
    # Duplicating a selected channel at a coordinate is an error, never a last-row selection.
    with table.open("a", encoding="utf-8", newline="") as stream:
        csv.DictWriter(stream, fieldnames=columns).writerow(rows[0])
    _update_table_hash(bundle, table.name)
    with pytest.raises(ComsolMappingError, match="duplicate coordinates"):
        prepare_positive_column_run(mapping_path, bundle_path=bundle)


def test_positive_column_legacy_reaction_mapping_requires_migration(
    tmp_path: Path,
) -> None:
    mapping_path, _bundle = _write_positive_column_repo(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace("    species: Ar\n", ""),
        encoding="utf-8",
    )
    with pytest.raises(ComsolMappingError, match="schema v2 migration required"):
        load_comsol_mapping(mapping_path)


def test_failed_new_positive_column_run_replaces_previous_acceptance(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping, bundle = _write_positive_column_repo(tmp_path)
    plan = prepare_positive_column_run(mapping, bundle_path=bundle)
    summary = plan.output_dir / "positive_column_run_summary.json"
    summary.write_text(
        json.dumps({"status": "completed", "quality_passed": True}), encoding="utf-8"
    )
    path = bundle / "transport_vs_en.csv"
    original = path.read_bytes()
    path.write_bytes(original + b"\n")
    with pytest.raises(PositiveColumnWorkflowError, match="artifact hash mismatch"):
        execute_positive_column_run(mapping, bundle_path=bundle)
    failed = json.loads(summary.read_text(encoding="utf-8"))
    assert failed["status"] == "blocked"
    assert failed["execution_stage"] == "preflight"
    assert failed["quality_passed"] is False
    path.write_bytes(original)

    def fail_apply(*args, **kwargs):
        assert json.loads(summary.read_text(encoding="utf-8"))["status"] == "running"
        raise RuntimeError("test COMSOL compile failure")

    monkeypatch.setattr(
        "swarm_workflow.comsol.models.positive_column.workflow.execute_apply_comsol",
        fail_apply,
    )
    with pytest.raises(RuntimeError, match="compile failure"):
        execute_positive_column_run(mapping, bundle_path=bundle)
    failed = json.loads(summary.read_text(encoding="utf-8"))
    assert failed["status"] == "failed"
    assert failed["execution_stage"] == "apply"
    assert failed["quality_passed"] is False


def _update_table_hash(bundle: Path, name: str) -> None:
    path = bundle / "manifest.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    manifest["tables"][name]["sha256"] = sha256(
        (bundle / name).read_bytes()
    ).hexdigest()
    path.write_text(json.dumps(manifest), encoding="utf-8")


def test_positive_column_rejects_missing_result_field(tmp_path: Path) -> None:
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
        "output_mph: ../work/positive_column_1d_swarm_tables.mph",
        "output_mph: ../positive_column_1d.mph",
    )
    mapping_path.write_text(text, encoding="utf-8")

    with pytest.raises(PositiveColumnWorkflowError, match="overwrite template"):
        prepare_positive_column_run(mapping_path, bundle_path=bundle)


@pytest.mark.parametrize(
    "produce_csv,scaled_units", [(True, False), (False, False), (True, True)]
)
def test_execute_positive_column_run_invokes_steps_in_order(
    tmp_path: Path,
    produce_csv: bool,
    scaled_units: bool,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_path, bundle = _write_positive_column_repo(tmp_path)
    scales = {
        "x": ("cm", 100.0),
        "electron_density": ("1/cm^3", 1e-6),
        "mean_electron_energy": ("J", 1.602176634e-19),
        "E_over_N": ("Td", 1e21),
        "reduced_mobility": ("1/(V*cm*s)", 0.01),
        "reduced_diffusion_L": ("1/(cm*s)", 0.01),
        "excitation_townsend": ("cm^2", 1e4),
        "ionization_townsend": ("cm^2", 1e4),
        "gas_pressure": ("kPa", 0.001),
        "applied_voltage": ("kV", 0.001),
    }
    if scaled_units:
        raw = yaml.safe_load(mapping_path.read_text())
        for field in raw["results"]["fields"]:
            if field["name"] in scales:
                field["unit"] = scales[field["name"]][0]
        mapping_path.write_text(yaml.safe_dump(raw))
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

    def fake_export(
        context,
        java_path,
        *,
        operation,
        comsol_executable=None,
        study,
    ):
        assert study == "std1"
        calls.append(operation)
        if operation == "positive_column_export" and produce_csv:
            result = (
                tmp_path
                / "outputs"
                / "comsol_positive_column"
                / "positive_column_results.csv"
            )
            result.parent.mkdir(parents=True, exist_ok=True)
            result.write_text(
                result_csv(
                    "x,electron_density,E_over_N,mean_electron_energy,"
                    "reduced_mobility,reduced_diffusion_L,"
                    "excitation_townsend,ionization_townsend,"
                    "applied_voltage,gas_pressure,electric_potential,electron_current_density,ion_current_density,total_current_density,excitation_source,ionization_source\n"
                    "0,1,4e-20,4,10.6666666666667,20.6666666666667,"
                    "1.333333333333e-21,"
                    "1.333333333333e-22,200,100,0,-1,1,0,1,1\n"
                ),
                encoding="utf-8",
            )
        if operation == "positive_column_export" and produce_csv and scaled_units:
            with result.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            for row in rows:
                for name, (_, factor) in scales.items():
                    row[name] = float(row[name]) * factor
            with result.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
        return SimpleNamespace(log_dir=tmp_path / "logs" / "export")

    monkeypatch.setattr(
        "swarm_workflow.comsol.models.positive_column.workflow.execute_apply_comsol",
        fake_apply,
    )
    monkeypatch.setattr(
        "swarm_workflow.comsol.models.positive_column.workflow."
        "execute_verify_comsol_functions",
        fake_verify,
    )
    monkeypatch.setattr(
        "swarm_workflow.comsol.models.positive_column.workflow."
        "execute_generated_comsol_java",
        fake_export,
    )

    if not produce_csv:
        previous = prepare_positive_column_run(
            mapping_path, bundle_path=bundle
        ).result_csv_path
        previous.write_text("stale successful result", encoding="utf-8")
        with pytest.raises(PositiveColumnWorkflowError, match="was not produced"):
            execute_positive_column_run(
                mapping_path,
                bundle_path=bundle,
                pressure_Pa=100.0,
                gas_temperature_K=293.15,
                mesh_elements=200,
            )
        assert not previous.exists()
        return

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
    from swarm_workflow.comsol.models.positive_column.results import (
        read_positive_column_results,
    )

    restored = read_positive_column_results(summary.run_summary_json)
    assert restored.units["mean_electron_energy"] == "eV"
    assert restored.rows[0]["mean_electron_energy"] == pytest.approx(4.0)
    assert restored.rows[0]["E_over_N"] == pytest.approx(4e-20)
    assert restored.metadata["source_units"]["mean_electron_energy"] == (
        "J" if scaled_units else "eV"
    )
    assert run_summary["final_voltage_V"] == 200.0
    assert run_summary["executed_voltage_sequence_V"] is None
    assert run_summary["provenance"] == {
        "apply": None,
        "verify": None,
        "run": None,
        "export": None,
    }
    assert run_summary["closure_checks"][
        "mobility_max_relative_error"
    ] == pytest.approx(0)
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
            "run-positive-column",
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
schema_version: 2
model:
  input_mph: ../positive_column_1d.mph
  output_mph: ../work/positive_column_1d_swarm_tables.mph
  study: std1
  component: comp1
  physics: plas
bundle:
  path: ../../outputs/placeholder_bundle
logs:
  path: ../logs
verify:
  output_path: ../../outputs/comsol_verify
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
  output_directory: ../../outputs/comsol_positive_column
  dataset: dset1
  solutions: last
  fields:
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
    complete_mapping_fixture(mapping_path)
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
                "format_version": 2,
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
    complete_bundle_fixture(bundle)
