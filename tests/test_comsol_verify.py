from __future__ import annotations

import csv
import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from swarm_workflow.comsol_adapter import find_comsol_executable
from swarm_workflow.comsol_verify import (
    ComsolVerifyError,
    compare_function_values,
    execute_verify_comsol_functions,
    format_verify_plan,
    generate_verify_java_source,
    plan_verify_comsol_functions,
    write_expected_values_csv,
)


def test_verify_plan_generates_control_points_from_csv_values(tmp_path: Path) -> None:
    mapping_path = _write_verify_repo(tmp_path)

    plan = plan_verify_comsol_functions(mapping_path)

    assert plan.output_dir == (tmp_path / "outputs" / "comsol_verify").resolve()
    assert plan.relative_tolerance == pytest.approx(1.0e-5)
    assert plan.absolute_tolerance == pytest.approx(1.0e-11)
    assert {point.function_tag for point in plan.points} == {"sw_meanE", "sw_muN"}
    mu_points = [point for point in plan.points if point.function_tag == "sw_muN"]
    assert [point.reason for point in mu_points] == ["min", "mid", "max"]
    assert [point.argument for point in mu_points] == [1.0e-20, 2.0e-20, 4.0e-20]
    assert [point.expected_value for point in mu_points] == [10.0, 20.0, 40.0]
    assert {point.quantity for point in plan.feature_points} == {
        "mean_energy",
        "reduced_mobility",
        "reduced_longitudinal_diffusion",
        "reduced_energy_mobility",
        "reduced_energy_diffusion",
        "excitation_townsend",
        "ionization_townsend",
    }
    assert all(point.reason == "feature_table_row" for point in plan.feature_points)
    assert len(plan.feature_points) == 21
    energy_mobility = [
        point
        for point in plan.feature_points
        if point.quantity == "reduced_energy_mobility"
    ]
    energy_diffusion = [
        point
        for point in plan.feature_points
        if point.quantity == "reduced_energy_diffusion"
    ]
    assert {point.column for point in energy_mobility} == {
        "reduced_electron_energy_mobility_m2_V_s_m3"
    }
    assert {point.fununit for point in energy_mobility} == {"1/(V*m*s)"}
    assert {point.column for point in energy_diffusion} == {
        "reduced_electron_energy_diffusion_m2_s_m3"
    }
    assert {point.fununit for point in energy_diffusion} == {"1/(m*s)"}


def test_verify_dry_run_format_reports_planned_points(tmp_path: Path) -> None:
    plan = plan_verify_comsol_functions(_write_verify_repo(tmp_path))
    text = format_verify_plan(plan)

    assert "COMSOL function verify dry-run" in text
    assert "sw_meanE" in text
    assert "sw_muN" in text
    assert "function_values_expected.csv" in text
    assert "verify_manifest.json" in text
    assert "COMSOL command: not executed" in text
    assert "relative_tolerance: 1e-05" in text


def test_verify_java_source_evaluates_mapped_function_tags(tmp_path: Path) -> None:
    plan = plan_verify_comsol_functions(_write_verify_repo(tmp_path))

    source = generate_verify_java_source(plan)

    assert "ModelUtil.loadCopy" in source
    assert "SWARM_COMSOL_VERIFY_CSV_BEGIN" in source
    assert 'getString("table")' in source
    assert "tableValue(model" in source
    assert 'model.func("sw_meanE")' in source
    assert 'model.func("sw_muN")' in source
    assert 'tableValue(model, "sw_muN", 1e-20)' in source
    assert "featureTableValue(model" in source
    assert '.getDoubleArray(xProperty)' in source
    assert '"deNXdata"' in source
    assert '"mueNYdata"' in source
    assert '"denNYdata"' in source
    assert '"xtownratedata"' in source


def test_compare_function_values_writes_pass_fail_summary(tmp_path: Path) -> None:
    plan = plan_verify_comsol_functions(_write_verify_repo(tmp_path))
    write_expected_values_csv(plan)
    _write_comsol_values_from_expected(plan, offset_by_tag={"sw_muN": 1.0e-3})

    summary = compare_function_values(plan)

    rows = _read_csv(plan.summary_csv)
    assert summary.passed > 0
    assert summary.failed > 0
    assert any(row["passed"] == "0" and row["function_tag"] == "sw_muN" for row in rows)
    assert {
        "absolute_error",
        "relative_error",
        "expected_unit",
        "verification_kind",
        "quantity",
    } <= set(rows[0])
    payload = json.loads(plan.summary_json.read_text(encoding="utf-8"))
    assert len(payload["closure_items"]) == 7
    assert all(item["status"] == "passed" for item in payload["closure_items"])


def test_townsend_readback_uses_relative_not_dimensionally_large_absolute_tolerance(
    tmp_path: Path,
) -> None:
    plan = plan_verify_comsol_functions(_write_verify_repo(tmp_path))
    write_expected_values_csv(plan)
    _write_comsol_values_from_expected(
        plan,
        offset_by_tag={"eir4.ytownratedata": 1.0e-22},
    )

    summary = compare_function_values(plan)

    assert summary.failed > 0
    payload = json.loads(plan.summary_json.read_text(encoding="utf-8"))
    ionization = next(
        item
        for item in payload["closure_items"]
        if item["quantity"] == "ionization_townsend"
    )
    assert ionization["status"] == "failed"


def test_execute_verify_writes_expected_java_and_summary_without_comsol(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_path = _write_verify_repo(tmp_path, output_mph_exists=True)

    def fake_execute(mapping, java_path, *, operation, comsol_executable=None):
        assert operation == "verify"
        assert Path(java_path).exists()
        plan = plan_verify_comsol_functions(mapping.path)
        _write_comsol_values_from_expected(plan)
        return SimpleNamespace(
            log_dir=tmp_path / "model" / "logs" / "verify_fake",
            command_json=tmp_path / "command.json",
            result_json=tmp_path / "result.json",
        )

    monkeypatch.setattr(
        "swarm_workflow.comsol_verify.execute_generated_comsol_java",
        fake_execute,
    )

    summary = execute_verify_comsol_functions(
        mapping_path,
        comsol_executable=tmp_path / "bin" / "comsol.exe",
    )

    assert summary.failed == 0
    assert summary.expected_csv.exists()
    assert summary.comsol_csv.exists()
    assert summary.summary_csv.exists()
    assert summary.plan.manifest_path.exists()
    assert summary.java_path.exists()
    manifest = json.loads(summary.plan.manifest_path.read_text(encoding="utf-8"))
    assert manifest["status"] == "passed"
    assert manifest["control_points"] == (
        len(summary.plan.points) + len(summary.plan.feature_points)
    )
    assert manifest["feature_table_control_points"] == len(
        summary.plan.feature_points
    )
    assert len(manifest["closure_items"]) == 7


def test_execute_verify_clears_stale_comsol_values_before_execution(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_path = _write_verify_repo(tmp_path, output_mph_exists=True)
    stale = tmp_path / "outputs" / "comsol_verify" / "function_values_comsol.csv"
    stale.parent.mkdir(parents=True)
    stale.write_text("old,row\n1,2\n", encoding="utf-8")

    def fake_execute(mapping, java_path, *, operation, comsol_executable=None):
        plan = plan_verify_comsol_functions(mapping.path)
        assert not plan.comsol_csv.exists()
        _write_comsol_values_from_expected(plan)
        return SimpleNamespace(
            log_dir=tmp_path / "model" / "logs" / "verify_fake",
            command_json=tmp_path / "command.json",
            result_json=tmp_path / "result.json",
        )

    monkeypatch.setattr(
        "swarm_workflow.comsol_verify.execute_generated_comsol_java",
        fake_execute,
    )

    summary = execute_verify_comsol_functions(mapping_path)

    assert summary.failed == 0
    assert _read_csv(summary.comsol_csv)[0]["function_name"] == "mean_energy_vs_en"


def test_verify_rejects_missing_output_mph_for_execution(tmp_path: Path) -> None:
    mapping_path = _write_verify_repo(tmp_path, output_mph_exists=False)

    with pytest.raises(ComsolVerifyError, match="output .mph does not exist"):
        execute_verify_comsol_functions(mapping_path, comsol_executable="comsol")


def test_verify_rejects_missing_mapped_csv(tmp_path: Path) -> None:
    mapping_path = _write_verify_repo(tmp_path)
    (tmp_path / "outputs" / "comsol_bundle" / "mixture_0000" / "transport_vs_en.csv").unlink()

    with pytest.raises(ComsolVerifyError, match="mapped function CSV does not exist"):
        plan_verify_comsol_functions(mapping_path)


def test_verify_rejects_missing_mapped_column(tmp_path: Path) -> None:
    mapping_path = _write_verify_repo(tmp_path)
    transport = tmp_path / "outputs" / "comsol_bundle" / "mixture_0000" / "transport_vs_en.csv"
    transport.write_text(
        "E_over_N_Td,E_over_N_V_m2,mean_energy_eV\n10,1e-20,2\n",
        encoding="utf-8",
    )

    with pytest.raises(ComsolVerifyError, match="reduced_mobility"):
        plan_verify_comsol_functions(mapping_path)


def test_verify_rejects_unit_mismatch(tmp_path: Path) -> None:
    mapping_path = _write_verify_repo(tmp_path)
    manifest_path = tmp_path / "outputs" / "comsol_bundle" / "mixture_0000" / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["tables"]["transport_vs_en.csv"]["units"][
        "reduced_mobility_m2_V_s_m3"
    ] = "kg"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(ComsolVerifyError, match="unit mismatch"):
        plan_verify_comsol_functions(mapping_path)


def test_verify_rejects_negative_expected_value(tmp_path: Path) -> None:
    mapping_path = _write_verify_repo(tmp_path)
    transport = tmp_path / "outputs" / "comsol_bundle" / "mixture_0000" / "transport_vs_en.csv"
    transport.write_text(
        (
            "E_over_N_Td,E_over_N_V_m2,mean_energy_eV,"
            "reduced_mobility_m2_V_s_m3\n"
            "10,1e-20,2,-10\n"
        ),
        encoding="utf-8",
    )

    with pytest.raises(ComsolVerifyError, match="negative expected value"):
        plan_verify_comsol_functions(mapping_path)


def test_verify_rejects_unsupported_nargs(tmp_path: Path) -> None:
    mapping_path = _write_verify_repo(tmp_path)
    text = mapping_path.read_text(encoding="utf-8").replace(
        "    nargs: 1\n    argunit: V*m^2\n",
        "    nargs: 2\n    argunit: V*m^2\n",
        1,
    )
    mapping_path.write_text(text, encoding="utf-8")

    with pytest.raises(ComsolVerifyError, match="supports nargs=1 only"):
        plan_verify_comsol_functions(mapping_path)


@pytest.mark.comsol
@pytest.mark.slow
def test_optional_comsol_verify_integration_requires_local_prerequisites() -> None:
    resolved = find_comsol_executable()
    if resolved is None:
        pytest.skip(
            "COMSOL executable not found; set COMSOL_BATCH, COMSOL_EXECUTABLE, "
            "or put comsol/comsolbatch on PATH"
        )
    mapping_value = os.environ.get("COMSOL_VERIFY_TEST_MAPPING")
    if not mapping_value:
        pytest.skip(
            "COMSOL_VERIFY_TEST_MAPPING is not set to a valid local mapping YAML"
        )
    mapping_path = Path(mapping_value)
    if not mapping_path.exists():
        pytest.skip(f"COMSOL_VERIFY_TEST_MAPPING does not exist: {mapping_path}")

    summary = execute_verify_comsol_functions(
        mapping_path,
        comsol_executable=resolved.path,
    )

    assert summary.summary_csv.exists()
    assert summary.failed == 0


def _write_verify_repo(root: Path, *, output_mph_exists: bool = False) -> Path:
    (root / "pyproject.toml").write_text("[project]\nname='tmp'\n", encoding="utf-8")
    model_dir = root / "model"
    maps_dir = model_dir / "maps"
    bundle_dir = root / "outputs" / "comsol_bundle" / "mixture_0000"
    maps_dir.mkdir(parents=True)
    bundle_dir.mkdir(parents=True)
    if output_mph_exists:
        output = model_dir / "work" / "fake_out.mph"
        output.parent.mkdir(parents=True)
        output.write_bytes(b"fake output mph")
    _write_bundle(bundle_dir)
    mapping_path = maps_dir / "verify.yaml"
    mapping_path.write_text(
        """
model:
  input_mph: model/fake_input.mph
  output_mph: model/work/fake_out.mph
  study: std1
  component: comp1
  physics: plas
bundle:
  path: outputs/comsol_bundle/mixture_0000
logs:
  path: model/logs
verify:
  output_path: outputs/comsol_verify
  relative_tolerance: 1.0e-5
  absolute_tolerance: 1.0e-11
functions:
  mean_energy_vs_en:
    tag: sw_meanE
    file: mean_energy_vs_en.csv
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
""".lstrip(),
        encoding="utf-8",
    )
    return mapping_path


def _write_bundle(bundle_dir: Path) -> None:
    (bundle_dir / "mean_energy_vs_en.csv").write_text(
        (
            "E_over_N_Td,E_over_N_V_m2,mean_energy_eV\n"
            "10,1e-20,2\n"
            "20,2e-20,4\n"
            "40,4e-20,8\n"
        ),
        encoding="utf-8",
    )
    (bundle_dir / "transport_vs_en.csv").write_text(
        (
            "E_over_N_Td,E_over_N_V_m2,mean_energy_eV,"
            "reduced_mobility_m2_V_s_m3,reduced_diffusion_L_m2_s_m3,"
            "reduced_electron_energy_mobility_m2_V_s_m3,"
            "reduced_electron_energy_diffusion_m2_s_m3\n"
            "10,1e-20,2,10,11,12,13\n"
            "20,2e-20,4,20,21,22,23\n"
            "40,4e-20,8,40,41,42,43\n"
        ),
        encoding="utf-8",
    )
    (bundle_dir / "transport_vs_mean_energy.csv").write_text(
        (
            "mean_energy_eV,E_over_N_Td,E_over_N_V_m2,"
            "reduced_mobility_m2_V_s_m3,reduced_diffusion_L_m2_s_m3,"
            "reduced_electron_energy_mobility_m2_V_s_m3,"
            "reduced_electron_energy_diffusion_m2_s_m3\n"
            "2,10,1e-20,10,11,12,13\n"
            "4,20,2e-20,20,21,22,23\n"
            "8,40,4e-20,40,41,42,43\n"
        ),
        encoding="utf-8",
    )
    (bundle_dir / "rates_vs_mean_energy.csv").write_text(
        (
            "mean_energy_eV,process_type,reduced_townsend_m2\n"
            "2,excitation,1e-21\n"
            "4,excitation,2e-21\n"
            "8,excitation,4e-21\n"
            "2,ionization,1e-22\n"
            "4,ionization,2e-22\n"
            "8,ionization,4e-22\n"
        ),
        encoding="utf-8",
    )
    units = {
        "E_over_N_Td": "Td",
        "E_over_N_V_m2": "V m^2",
        "mean_energy_eV": "eV",
        "reduced_mobility_m2_V_s_m3": "m^2/(V s) m^3",
        "reduced_diffusion_L_m2_s_m3": "m^2/s m^3",
        "reduced_electron_energy_mobility_m2_V_s_m3": (
            "1/(V m s)"
        ),
        "reduced_electron_energy_diffusion_m2_s_m3": "1/(m s)",
        "process_type": "",
        "reduced_townsend_m2": "m^2",
    }
    manifest = {
        "format_version": 1,
        "tables": {
            "mean_energy_vs_en.csv": {
                "columns": ["E_over_N_Td", "E_over_N_V_m2", "mean_energy_eV"],
                "units": {
                    "E_over_N_V_m2": units["E_over_N_V_m2"],
                    "mean_energy_eV": units["mean_energy_eV"],
                },
            },
            "transport_vs_en.csv": {
                "columns": [
                    "E_over_N_Td",
                    "E_over_N_V_m2",
                    "mean_energy_eV",
                    "reduced_mobility_m2_V_s_m3",
                    "reduced_diffusion_L_m2_s_m3",
                    "reduced_electron_energy_mobility_m2_V_s_m3",
                    "reduced_electron_energy_diffusion_m2_s_m3",
                ],
                "units": {
                    "E_over_N_V_m2": units["E_over_N_V_m2"],
                    "reduced_mobility_m2_V_s_m3": units[
                        "reduced_mobility_m2_V_s_m3"
                    ],
                },
            },
            "transport_vs_mean_energy.csv": {
                "columns": [
                    "mean_energy_eV",
                    "E_over_N_Td",
                    "E_over_N_V_m2",
                    "reduced_mobility_m2_V_s_m3",
                    "reduced_diffusion_L_m2_s_m3",
                    "reduced_electron_energy_mobility_m2_V_s_m3",
                    "reduced_electron_energy_diffusion_m2_s_m3",
                ],
                "units": units,
            },
            "rates_vs_mean_energy.csv": {
                "columns": [
                    "mean_energy_eV",
                    "process_type",
                    "reduced_townsend_m2",
                ],
                "units": units,
            },
        },
        "units": units,
    }
    (bundle_dir / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")


def _write_comsol_values_from_expected(
    plan,
    *,
    offset_by_tag: dict[str, float] | None = None,
) -> None:
    offset_by_tag = offset_by_tag or {}
    expected = _read_csv(plan.expected_csv)
    with plan.comsol_csv.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(
            fp,
            fieldnames=["function_name", "function_tag", "argument", "comsol_value"],
        )
        writer.writeheader()
        for row in expected:
            value = float(row["expected_value"]) + offset_by_tag.get(
                row["function_tag"],
                0.0,
            )
            writer.writerow(
                {
                    "function_name": row["function_name"],
                    "function_tag": row["function_tag"],
                    "argument": row["argument"],
                    "comsol_value": f"{value:.17g}",
                }
            )


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fp:
        return list(csv.DictReader(fp))
