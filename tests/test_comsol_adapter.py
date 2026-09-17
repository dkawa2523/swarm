from __future__ import annotations

import json
import os
from pathlib import Path
import re
import subprocess

import pytest

import swarm_workflow.comsol.runtime as comsol_runtime
from swarm_workflow.comsol.java import compose_java_mains
from swarm_workflow.comsol.models.positive_column.config import (
    ComsolMappingError,
    comsol_run_context,
    load_comsol_mapping,
    validate_comsol_mapping_files,
)
from swarm_workflow.comsol.models.positive_column.java import (
    generate_apply_java_source,
)
from swarm_workflow.comsol.models.positive_column.workflow import execute_apply_comsol
from swarm_workflow.comsol.runtime import (
    ComsolAdapterError,
    build_java_batch_command,
    execute_generated_comsol_java,
    resolve_comsol_executable,
)


ROOT = Path(__file__).resolve().parents[1]


def test_generated_java_stages_share_one_ordered_entry_point() -> None:
    source = compose_java_mains(
        "Combined",
        (
            (
                "Apply",
                "import com.comsol.model.Model;\npublic class Apply {\n  public static void main(String[] args) throws Exception {}\n}\n",
            ),
            (
                "Verify",
                "import com.comsol.model.Model;\npublic class Verify {\n  public static void main(String[] args) throws Exception {}\n}\n",
            ),
        ),
    )

    assert source.count("import com.comsol.model.Model;") == 1
    assert "final class Apply" in source
    assert "final class Verify" in source
    assert source.count("public class ") == 1
    assert source.index("Apply.main(args);") < source.index("Verify.main(args);")


def test_adapter_exposes_execution_api_without_legacy_planning_facade() -> None:
    obsolete = (
        "ApplyComsolPlan",
        "JavaWriteSummary",
        "find_comsol_executable",
        "format_apply_plan",
        "plan_apply_comsol",
        "write_apply_comsol_java",
    )
    assert all(not hasattr(comsol_runtime, name) for name in obsolete)


def test_java_source_contains_mapping_functions_and_standard_properties(
    tmp_path: Path,
) -> None:
    mapping = load_comsol_mapping(_write_valid_apply_repo(tmp_path))

    source = generate_apply_java_source(mapping)

    assert "ModelUtil.loadCopy" in source
    assert "model.save" in source
    assert 'ensureInterpolation(model, "sw_meanE")' in source
    assert 'ensureInterpolation(model, "sw_muN")' in source
    assert source.index('"sw_meanE"') < source.index('"sw_muN"')
    assert "mean_energy_vs_en.csv" in source
    assert "transport_vs_en.csv" in source
    assert '.set("source", "file")' in source
    assert '.set("nargs", 1)' in source
    assert '.set("argunit", "V*m^2")' in source
    assert '.set("fununit", "1/(V*m*s)")' in source
    assert '.set("interp", "linear")' in source
    assert '.set("extrap", "const")' in source
    assert ".importData()" in source


def test_java_source_contains_mapping_driven_closure_and_reaction_lookup(
    tmp_path: Path,
) -> None:
    mapping_path = _write_valid_apply_repo(tmp_path)
    mapping = load_comsol_mapping(mapping_path)

    source = generate_apply_java_source(mapping)

    assert '"SpecifyElectronDensityAndEnergy", "UseLookupTables"' in source
    assert '"SpecifyMeanElectronEnergy", "MeanEnergyTable"' in source
    assert '.feature("pes1").set("enrgXdata", new double[]{' in source
    assert '.feature("pes1").set("enrgYdata", new double[]{' in source
    assert '.feature("eir2").set("SpecifyReactionUsing", "UseLookupTable")' in source
    assert '.feature("eir2").set("RateConstantForm", "UseTownsend")' in source
    assert '.feature("eir2").set("UseTownsendTwoTermsBoltzmann", false)' in source
    assert '.feature("eir2").set("xtownratedata", new double[]{' in source
    assert '.feature("eir4").set("ytownratedata", new double[]{' in source
    assert "2.00000000000000000e+00" in source
    assert re.search(r"ytownratedata\", new double\[\]\{[^}]*e-2[12]", source)

    assert mapping.closure.mean_energy.table == "mean_energy_vs_en.csv"
    assert mapping.closure.mean_energy_formulation.mode == "local_energy"
    assert mapping.reaction_lookups[0].feature == "eir2"
    assert mapping.reaction_lookups[0].source_model == "townsend_flux"


def test_apply_comsol_rejects_missing_unit_metadata(tmp_path: Path) -> None:
    mapping_path = _write_valid_apply_repo(tmp_path, include_units=False)
    mapping = load_comsol_mapping(mapping_path)

    with pytest.raises(ComsolMappingError, match="no unit metadata"):
        validate_comsol_mapping_files(mapping, require_unit_metadata=True)


def test_apply_comsol_rejects_missing_required_mapping_field(tmp_path: Path) -> None:
    mapping_path = _write_valid_apply_repo(tmp_path)
    text = mapping_path.read_text(encoding="utf-8").replace("    tag: sw_muN\n", "")
    mapping_path.write_text(text, encoding="utf-8")

    with pytest.raises(ComsolMappingError, match="tag must be a non-empty string"):
        load_comsol_mapping(mapping_path)


def test_resolve_comsol_executable_order(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cli = _fake_exe(tmp_path / "cli" / "comsol.exe")
    batch = _fake_exe(tmp_path / "env" / "comsolbatch.exe")
    executable = _fake_exe(tmp_path / "env" / "comsol.exe")
    path_comsol = _fake_exe(tmp_path / "path" / "comsol.exe")
    path_batch = _fake_exe(tmp_path / "path" / "comsolbatch.exe")

    monkeypatch.setenv("COMSOL_BATCH", str(batch))
    monkeypatch.setenv("COMSOL_EXECUTABLE", str(executable))
    monkeypatch.setattr(
        "swarm_workflow.comsol.runtime.executable.shutil.which",
        lambda name: str(path_comsol if name == "comsol" else path_batch),
    )

    resolved = resolve_comsol_executable(cli)
    assert resolved.path == str(cli.resolve())
    assert resolved.source == "cli"

    resolved = resolve_comsol_executable()
    assert resolved.path == str(batch.resolve())
    assert resolved.source == "env:COMSOL_BATCH"

    monkeypatch.delenv("COMSOL_BATCH")
    resolved = resolve_comsol_executable()
    assert resolved.path == str(executable.resolve())
    assert resolved.source == "env:COMSOL_EXECUTABLE"

    monkeypatch.delenv("COMSOL_EXECUTABLE")
    resolved = resolve_comsol_executable()
    assert resolved.path == str(path_comsol.resolve())
    assert resolved.source == "path:comsol"

    monkeypatch.setattr(
        "swarm_workflow.comsol.runtime.executable.shutil.which", lambda _: None
    )
    with pytest.raises(ComsolAdapterError, match="COMSOL executable not found"):
        resolve_comsol_executable()


def test_installed_license_is_explicitly_passed_to_comsol_batch(
    tmp_path: Path,
) -> None:
    installation = tmp_path / "Multiphysics"
    batch = _fake_exe(installation / "bin" / "win64" / "comsolbatch.exe")
    license_file = installation / "license" / "license.dat"
    license_file.parent.mkdir(parents=True)
    license_file.write_text("FEATURE test\n", encoding="utf-8")
    class_file = tmp_path / "Apply.class"

    resolved = resolve_comsol_executable(batch)
    command = build_java_batch_command(
        resolved.path,
        class_file,
        license_path=resolved.license_path,
    )

    assert resolved.license_path == str(license_file.resolve())
    assert resolved.license_source == "installation:license/license.dat"
    assert command == [
        str(batch.resolve()),
        "-c",
        str(license_file.resolve()),
        "-inputfile",
        str(class_file),
    ]


def test_execute_apply_comsol_runs_subprocess_with_safe_arguments_and_logs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_path = _write_valid_apply_repo(tmp_path)
    mapping = load_comsol_mapping(mapping_path)
    comsol = _fake_exe(tmp_path / "bin" / "comsol.exe")
    calls: list[tuple[list[str], dict[str, object]]] = []

    def fake_run(args: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append((args, kwargs))
        if Path(args[0]).name == "comsolcompile.exe":
            Path(args[1]).with_suffix(".class").write_bytes(b"class")
        if "-inputfile" in args and str(args[args.index("-inputfile") + 1]).endswith(
            ".class"
        ):
            mapping.model.output_mph.write_bytes(b"output mph")
        return subprocess.CompletedProcess(args, 0, stdout="ok out", stderr="ok err")

    monkeypatch.setattr(
        "swarm_workflow.comsol.runtime.process.subprocess.run", fake_run
    )

    summary = execute_apply_comsol(mapping_path, comsol_executable=comsol)

    assert len(calls) == 2
    for args, kwargs in calls:
        assert isinstance(args, list)
        assert kwargs["shell"] is False
        assert kwargs["check"] is True
        assert "capture_output" not in kwargs
        assert "text" not in kwargs
        assert kwargs["stdout"] is not None
        assert kwargs["stderr"] is not None
    assert calls[0][0][0] == str(comsol.with_name("comsolcompile.exe").resolve())
    assert calls[1][0][0] == str(comsol.with_name("comsolbatch.exe").resolve())

    command = json.loads(summary.command_json.read_text(encoding="utf-8"))
    result = json.loads(summary.result_json.read_text(encoding="utf-8"))
    assert command["generated_java_path"].endswith("SwarmComsolApply.java")
    assert command["output_mph"] == str(summary.output_mph)
    assert command["steps"][0]["args"] == calls[0][0]
    assert result["status"] == "ok"
    assert result["steps"][0]["return_code"] == 0
    assert summary.stdout_paths[0].read_text(encoding="utf-8") == "ok out"
    assert summary.stderr_paths[0].read_text(encoding="utf-8") == "ok err"
    provenance = json.loads(summary.provenance_json.read_text(encoding="utf-8"))
    assert provenance["status"] == "ok"
    assert provenance["inputs"]["mapping"]["sha256"]
    assert provenance["inputs"]["input_mph"]["sha256"]
    assert provenance["inputs"]["bundle_manifest"]["sha256"]
    assert provenance["generated_artifacts"]["staged_java"]["sha256"]
    assert provenance["generated_artifacts"]["compiled_class"]["sha256"]
    assert provenance["generated_artifacts"]["compiled_class"]["modified_at_utc"]
    assert "configured_voltage_sequence_V" not in provenance


def test_zero_exit_compile_error_does_not_execute_stale_class(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_path = _write_valid_apply_repo(tmp_path)
    mapping = load_comsol_mapping(mapping_path)
    comsol = _fake_exe(tmp_path / "bin" / "comsol.exe")
    java_path = tmp_path / "generated" / "StaleRun.java"
    java_path.parent.mkdir(parents=True)
    java_path.write_text("public class StaleRun {}\n", encoding="utf-8")
    stale_class = java_path.with_suffix(".class")
    stale_class.write_bytes(b"old class that must never execute")
    calls: list[list[str]] = []

    def fake_run(args: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(args)
        return subprocess.CompletedProcess(
            args,
            0,
            stdout=(
                "Failed to compile java file.\n"
                "Compilation error on line 29:\n"
                "ERROR: Compilation failed. For details, see the log file."
            ),
            stderr="",
        )

    monkeypatch.setattr(
        "swarm_workflow.comsol.runtime.process.subprocess.run", fake_run
    )

    with pytest.raises(ComsolAdapterError, match="batch step failed") as exc_info:
        execute_generated_comsol_java(
            comsol_run_context(mapping),
            java_path,
            operation="positive_column_run",
            comsol_executable=comsol,
            study=mapping.model.study,
        )

    assert len(calls) == 1
    assert stale_class.read_bytes() == b"old class that must never execute"
    step_result = exc_info.value.step_result
    assert step_result is not None
    assert step_result["return_code"] == 0
    assert step_result["detected_comsol_error"] == "compiler_failed"
    log_dir = Path(step_result["stdout"]).parent
    command = json.loads((log_dir / "command.json").read_text(encoding="utf-8"))
    assert command["source_java_path"] == str(java_path.resolve())
    assert Path(command["generated_java_path"]).parent == log_dir
    assert Path(command["generated_java_path"]) != java_path.resolve()
    result = json.loads((log_dir / "result.json").read_text(encoding="utf-8"))
    assert result["status"] == "failed"
    assert len(result["steps"]) == 1
    provenance = json.loads((log_dir / "provenance.json").read_text(encoding="utf-8"))
    assert provenance["status"] == "failed"
    assert provenance["generated_artifacts"]["source_origin"]["sha256"]
    assert provenance["generated_artifacts"]["compiled_class"]["exists"] is False


def test_generated_java_uses_fresh_per_run_class_and_records_comsol_build(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_path = _write_valid_apply_repo(tmp_path)
    mapping = load_comsol_mapping(mapping_path)
    comsol = _fake_exe(tmp_path / "bin" / "comsol.exe")
    java_path = tmp_path / "generated" / "FreshRun.java"
    java_path.parent.mkdir(parents=True)
    java_path.write_text("public class FreshRun {}\n", encoding="utf-8")
    calls: list[list[str]] = []

    def fake_run(args: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(args)
        if Path(args[0]).name == "comsolcompile.exe":
            Path(args[1]).with_suffix(".class").write_bytes(b"fresh class")
            stdout = "Compilation completed."
        else:
            stdout = (
                "***COMSOL 6.4.0.429 progress output file***\n"
                "COMSOL Multiphysics 6.4 (Build: 429) starting in batch mode\n"
                "SWARM_COMSOL_VOLTAGE_V=200[V]\n"
                "Direct target-voltage solve failed; using continuation.\n"
                "SWARM_COMSOL_VOLTAGE_V=20[V]\n"
                "SWARM_COMSOL_VOLTAGE_V=50[V]\n"
                "SWARM_COMSOL_VOLTAGE_V=100[V]\n"
                "SWARM_COMSOL_VOLTAGE_V=200[V]\n"
            )
        return subprocess.CompletedProcess(args, 0, stdout=stdout, stderr="")

    monkeypatch.setattr(
        "swarm_workflow.comsol.runtime.process.subprocess.run", fake_run
    )

    summary = execute_generated_comsol_java(
        comsol_run_context(mapping),
        java_path,
        operation="positive_column_run",
        comsol_executable=comsol,
        study=mapping.model.study,
    )

    assert len(calls) == 2
    executed_class = Path(calls[1][calls[1].index("-inputfile") + 1])
    assert executed_class.parent == summary.log_dir
    assert executed_class.read_bytes() == b"fresh class"
    assert summary.java_path == summary.log_dir / java_path.name
    assert summary.comsol_version == "6.4"
    assert summary.comsol_build == "429"
    assert not hasattr(summary, "executed_voltage_sequence_V")
    provenance = json.loads(summary.provenance_json.read_text(encoding="utf-8"))
    assert provenance["comsol"]["progress_version"] == "6.4.0.429"
    assert "voltage_execution_mode" not in provenance
    assert provenance["generated_artifacts"]["source_origin"]["path"] == str(
        java_path.resolve()
    )
    assert provenance["generated_artifacts"]["staged_java"]["path"] == str(
        summary.java_path
    )
    assert (
        provenance["generated_artifacts"]["source_origin"]["sha256"]
        == provenance["generated_artifacts"]["staged_java"]["sha256"]
    )


def test_generated_java_compiles_support_classes_for_one_batch_process(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_path = _write_valid_apply_repo(tmp_path)
    mapping = load_comsol_mapping(mapping_path)
    comsol = _fake_exe(tmp_path / "bin" / "comsol.exe")
    source_dir = tmp_path / "generated"
    source_dir.mkdir(parents=True)
    entry = source_dir / "Combined.java"
    support = source_dir / "Apply.java"
    entry.write_text(
        "public class Combined { public static void main(String[] args) { "
        "Apply.main(args); } }\n",
        encoding="utf-8",
    )
    support.write_text(
        "public class Apply { public static void main(String[] args) {} }\n",
        encoding="utf-8",
    )
    calls: list[list[str]] = []

    def fake_run(args: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        del kwargs
        calls.append(args)
        if Path(args[0]).name == "comsolcompile.exe":
            Path(args[1]).with_suffix(".class").write_bytes(b"fresh class")
            stdout = "Compilation completed."
        else:
            stdout = "COMSOL Multiphysics 6.4 (Build: 429) starting in batch mode"
        return subprocess.CompletedProcess(args, 0, stdout=stdout, stderr="")

    monkeypatch.setattr(
        "swarm_workflow.comsol.runtime.process.subprocess.run", fake_run
    )

    summary = execute_generated_comsol_java(
        comsol_run_context(mapping),
        entry,
        operation="combined",
        comsol_executable=comsol,
        study="not_applicable",
        support_java_paths=(support,),
    )

    assert [Path(call[0]).name for call in calls] == [
        "comsolcompile.exe",
        "comsolcompile.exe",
        "comsolbatch.exe",
    ]
    assert Path(calls[0][1]).name == "Combined.java"
    assert Path(calls[1][1]).name == "Apply.java"
    provenance = json.loads(summary.provenance_json.read_text(encoding="utf-8"))
    artifacts = provenance["generated_artifacts"]["support_classes"]
    assert len(artifacts) == 1
    assert (
        artifacts[0]["source_origin"]["sha256"] == artifacts[0]["staged_java"]["sha256"]
    )
    assert artifacts[0]["compiled_class"]["sha256"]


def test_compile_success_without_new_class_is_rejected(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_path = _write_valid_apply_repo(tmp_path)
    comsol = _fake_exe(tmp_path / "bin" / "comsol.exe")
    calls: list[list[str]] = []

    def fake_run(args: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(args)
        return subprocess.CompletedProcess(args, 0, stdout="ok", stderr="")

    monkeypatch.setattr(
        "swarm_workflow.comsol.runtime.process.subprocess.run", fake_run
    )

    with pytest.raises(ComsolAdapterError, match="did not produce a class file") as exc:
        execute_apply_comsol(mapping_path, comsol_executable=comsol)

    assert len(calls) == 1
    assert exc.value.step_result is not None
    assert exc.value.step_result["missing_artifact"].endswith(".class")


def test_apply_rejects_unchanged_preexisting_output_mph(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_path = _write_valid_apply_repo(tmp_path)
    mapping = load_comsol_mapping(mapping_path)
    mapping.model.output_mph.parent.mkdir(parents=True)
    mapping.model.output_mph.write_bytes(b"old output that must not be reused")
    comsol = _fake_exe(tmp_path / "bin" / "comsol.exe")
    calls: list[list[str]] = []

    def fake_run(args: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(args)
        if Path(args[0]).name == "comsolcompile.exe":
            Path(args[1]).with_suffix(".class").write_bytes(b"fresh class")
        return subprocess.CompletedProcess(args, 0, stdout="ok", stderr="")

    monkeypatch.setattr(
        "swarm_workflow.comsol.runtime.process.subprocess.run", fake_run
    )

    with pytest.raises(ComsolAdapterError, match="stale result") as exc:
        execute_apply_comsol(mapping_path, comsol_executable=comsol)

    assert len(calls) == 2
    assert exc.value.step_result is not None
    assert exc.value.step_result["stale_artifact"] == str(mapping.model.output_mph)
    log_dir = Path(exc.value.step_result["stdout"]).parent
    provenance = json.loads((log_dir / "provenance.json").read_text("utf-8"))
    assert provenance["status"] == "failed"
    assert provenance["inputs"]["preexisting_output_mph"]["sha256"]
    assert (
        provenance["inputs"]["preexisting_output_mph"]["sha256"]
        == provenance["generated_artifacts"]["output_mph"]["sha256"]
    )


def test_execute_failure_reports_stdout_and_stderr_log_paths(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_path = _write_valid_apply_repo(tmp_path)
    comsol = _fake_exe(tmp_path / "bin" / "comsol.exe")

    def fake_run(args: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        raise subprocess.CalledProcessError(
            7,
            args,
            output="bad stdout",
            stderr="bad stderr",
        )

    monkeypatch.setattr(
        "swarm_workflow.comsol.runtime.process.subprocess.run", fake_run
    )

    with pytest.raises(ComsolAdapterError) as exc_info:
        execute_apply_comsol(mapping_path, comsol_executable=comsol)

    message = str(exc_info.value)
    assert "log_dir:" in message
    assert "stdout:" in message
    assert "stderr:" in message
    step_result = exc_info.value.step_result
    assert step_result is not None
    stdout_path = Path(step_result["stdout"])
    stderr_path = Path(step_result["stderr"])
    assert stdout_path.read_text(encoding="utf-8") == "bad stdout"
    assert stderr_path.read_text(encoding="utf-8") == "bad stderr"
    aggregate = json.loads((stdout_path.parent / "result.json").read_text("utf-8"))
    assert aggregate["status"] == "failed"


def test_electron_swarm_still_has_no_workflow_comsol_or_mph_dependency() -> None:
    forbidden = ["swarm_workflow", "sqlite3", "subprocess", "comsol", "mph"]
    for path in (ROOT / "electron_swarm").rglob("*.py"):
        text = path.read_text(encoding="utf-8").lower()
        for item in forbidden:
            assert item not in text, f"{path} contains forbidden dependency {item}"


@pytest.mark.comsol
@pytest.mark.slow
def test_optional_comsol_apply_integration_requires_local_prerequisites() -> None:
    try:
        resolved = resolve_comsol_executable()
    except ComsolAdapterError:
        pytest.skip(
            "COMSOL executable not found; set COMSOL_BATCH, COMSOL_EXECUTABLE, "
            "or put comsol/comsolbatch on PATH"
        )
    mapping_value = os.environ.get("COMSOL_TEST_MAPPING")
    if not mapping_value:
        pytest.skip("COMSOL_TEST_MAPPING is not set to a valid local mapping YAML")
    mapping_path = Path(mapping_value)
    if not mapping_path.exists():
        pytest.skip(f"COMSOL_TEST_MAPPING does not exist: {mapping_path}")

    summary = execute_apply_comsol(mapping_path, comsol_executable=resolved.path)

    assert summary.result_json.exists()
    assert summary.output_mph.exists()


def _fake_exe(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("fake executable", encoding="utf-8")
    return path


def _write_valid_apply_repo(root: Path, *, include_units: bool = True) -> Path:
    (root / "pyproject.toml").write_text("[project]\nname='tmp'\n", encoding="utf-8")
    model_dir = root / "model"
    maps_dir = model_dir / "maps"
    bundle_dir = root / "outputs" / "comsol_bundle" / "mixture_0000"
    maps_dir.mkdir(parents=True)
    bundle_dir.mkdir(parents=True)
    (model_dir / "fake.mph").write_bytes(b"not a real mph")
    (bundle_dir / "mean_energy_vs_en.csv").write_text(
        "E_over_N_Td,E_over_N_V_m2,mean_energy_eV\n10,1e-20,2\n20,2e-20,3\n",
        encoding="utf-8",
    )
    (bundle_dir / "transport_vs_en.csv").write_text(
        (
            "E_over_N_Td,E_over_N_V_m2,mean_energy_eV,"
            "reduced_mobility_m2_V_s_m3,reduced_diffusion_L_m2_s_m3,"
            "reduced_electron_energy_mobility_m2_V_s_m3,"
            "reduced_electron_energy_diffusion_m2_s_m3\n"
            "10,1e-20,2,12,20,18,22\n"
            "20,2e-20,3,14,24,21,26\n"
        ),
        encoding="utf-8",
    )
    (bundle_dir / "transport_vs_mean_energy.csv").write_text(
        (bundle_dir / "transport_vs_en.csv").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    (bundle_dir / "rates_vs_mean_energy.csv").write_text(
        (
            "mean_energy_eV,E_over_N_Td,E_over_N_V_m2,species,process,"
            "process_type,threshold_eV,target_species_fraction,"
            "rate_coefficient_m3_s,mixture_weighted_rate_m3_s,"
            "reduced_townsend_m2,mixture_weighted_reduced_townsend_m2\n"
            "2,10,1e-20,Ar,e+Ar=>e+Ars,excitation,11.5,1,"
            "1e-15,1e-15,2e-21,2e-21\n"
            "3,20,2e-20,Ar,e+Ar=>e+Ars,excitation,11.5,1,"
            "2e-15,2e-15,3e-21,3e-21\n"
            "2,10,1e-20,Ar,e+Ar=>2e+Ar+,ionization,15.8,1,"
            "1e-16,1e-16,4e-22,4e-22\n"
            "3,20,2e-20,Ar,e+Ar=>2e+Ar+,ionization,15.8,1,"
            "2e-16,2e-16,5e-22,5e-22\n"
        ),
        encoding="utf-8",
    )
    manifest = {
        "format_version": 1,
        "tables": {
            "mean_energy_vs_en.csv": {
                "columns": ["E_over_N_Td", "E_over_N_V_m2", "mean_energy_eV"],
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
            },
            "transport_vs_mean_energy.csv": {
                "columns": [
                    "E_over_N_Td",
                    "E_over_N_V_m2",
                    "mean_energy_eV",
                    "reduced_mobility_m2_V_s_m3",
                    "reduced_diffusion_L_m2_s_m3",
                    "reduced_electron_energy_mobility_m2_V_s_m3",
                    "reduced_electron_energy_diffusion_m2_s_m3",
                ],
            },
            "rates_vs_mean_energy.csv": {
                "columns": [
                    "mean_energy_eV",
                    "process_type",
                    "reduced_townsend_m2",
                ],
            },
        },
    }
    if include_units:
        manifest["tables"]["mean_energy_vs_en.csv"]["units"] = {
            "E_over_N_V_m2": "V m^2",
            "mean_energy_eV": "eV",
        }
        manifest["tables"]["transport_vs_en.csv"]["units"] = {
            "E_over_N_V_m2": "V m^2",
            "reduced_mobility_m2_V_s_m3": "m^2/(V s) m^3",
            "reduced_diffusion_L_m2_s_m3": "m^2/s m^3",
            "reduced_electron_energy_mobility_m2_V_s_m3": "1/(V m s)",
            "reduced_electron_energy_diffusion_m2_s_m3": "1/(m s)",
        }
        manifest["tables"]["transport_vs_mean_energy.csv"]["units"] = dict(
            manifest["tables"]["transport_vs_en.csv"]["units"]
        )
    (bundle_dir / "manifest.json").write_text(
        json.dumps(manifest),
        encoding="utf-8",
    )
    mapping_path = maps_dir / "valid.yaml"
    mapping_path.write_text(
        """
schema_version: 2
model:
  input_mph: model/fake.mph
  output_mph: model/work/fake_out.mph
  study: std1
  component: comp1
  physics: plas
bundle:
  path: outputs/comsol_bundle/mixture_0000
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
results:
  probes:
    - name: electron_density_center
      expression: plas.ne
      unit: 1/m^3
""".lstrip(),
        encoding="utf-8",
    )
    return mapping_path
