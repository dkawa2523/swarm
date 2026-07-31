from __future__ import annotations

import csv
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
RUNNER_PATH = (
    ROOT
    / "reports"
    / "comsol_swarm_benchmark_2026"
    / "repro"
    / "run_builtin_benchmark.py"
)
SPEC = importlib.util.spec_from_file_location("builtin_benchmark_runner", RUNNER_PATH)
assert SPEC is not None and SPEC.loader is not None
runner = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = runner
SPEC.loader.exec_module(runner)


def test_java_source_uses_sol3_exact_profile_and_mesh_contract(
    tmp_path: Path,
) -> None:
    source = runner.render_java_source(
        input_mph=tmp_path / "input.mph",
        output_mph=tmp_path / "output.mph",
        profile_csv=tmp_path / "profile.csv",
        voltage_V=200.0,
        pressure_Pa=13.3322,
        gas_temperature_K=293.15,
        mesh_elements=400,
    )

    assert 'model.sol("sol3").runAll();' in source
    assert '.feature("dis1").set("elemcount", 400);' in source
    assert 'model.component("comp1").mesh("mesh1").run();' in source
    assert "public static void main(String[] args) throws Exception" in source
    assert ",".join(runner.PROFILE_HEADER) in source
    assert "plas.Jix_wAr_1p+plas.Jelx" in source
    assert "SWARM_BUILTIN_SOLVE_TIME_S=" in source


def test_exit_zero_compile_error_refuses_batch_and_records_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_mph = tmp_path / "input.mph"
    input_mph.write_bytes(b"input model")
    comsol = _fake_exe(tmp_path / "bin" / "comsolbatch.exe")
    calls: list[list[str]] = []

    def fake_run(
        args: list[str], **kwargs: object
    ) -> subprocess.CompletedProcess[str]:
        calls.append(args)
        return subprocess.CompletedProcess(
            args,
            0,
            stdout=(
                "Failed to compile java file.\n"
                "ERROR: Compilation failed. For details, see the log file."
            ),
            stderr="",
        )

    monkeypatch.setattr(runner.subprocess, "run", fake_run)

    with pytest.raises(runner.BuiltinBenchmarkError, match="compile failed"):
        runner.run_repetitions(
            input_mph=input_mph,
            output_root=tmp_path / "runs",
            comsol_executable=comsol,
            repetitions=1,
        )

    assert len(calls) == 1
    provenance = json.loads(
        (tmp_path / "runs" / "run_01" / "provenance.json").read_text(
            encoding="utf-8"
        )
    )
    assert provenance["status"] == "failed"
    assert provenance["steps"][0]["return_code"] == 0
    assert (
        provenance["steps"][0]["detected_comsol_error"] == "compiler_failed"
    )
    assert provenance["artifacts"]["class"]["exists"] is False


def test_success_uses_fresh_class_validates_profile_and_writes_provenance(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    input_mph = tmp_path / "input.mph"
    input_mph.write_bytes(b"input model")
    comsol = _fake_exe(tmp_path / "bin" / "comsolbatch.exe")
    calls: list[list[str]] = []

    def fake_run(
        args: list[str], **kwargs: object
    ) -> subprocess.CompletedProcess[str]:
        calls.append(args)
        if Path(args[0]).name.lower() == "comsolcompile.exe":
            java = Path(args[1])
            java.with_suffix(".class").write_bytes(b"fresh class")
            os.utime(java.with_suffix(".class"), ns=(java.stat().st_mtime_ns,) * 2)
            stdout = "Compilation completed."
        else:
            run_dir = Path(str(kwargs["cwd"]))
            (run_dir / "positive_column_builtin_solved.mph").write_bytes(
                b"solved model"
            )
            _write_profile(run_dir / "profile.csv")
            stdout = (
                "***COMSOL 6.4.0.429 progress output file***\n"
                "COMSOL Multiphysics 6.4 (Build: 429) starting in batch mode\n"
                "SWARM_BUILTIN_VOLTAGE_V=200\n"
                "SWARM_BUILTIN_GAS_TEMPERATURE_K=293.14999999999998\n"
                "SWARM_BUILTIN_MESH_ELEMENTS=200\n"
                "SWARM_BUILTIN_SOLVE_TIME_S=370.25\n"
                "SWARM_BUILTIN_PROFILE_ROWS=2\n"
            )
        return subprocess.CompletedProcess(args, 0, stdout=stdout, stderr="")

    monkeypatch.setattr(runner.subprocess, "run", fake_run)

    results = runner.run_repetitions(
        input_mph=input_mph,
        output_root=tmp_path / "runs",
        comsol_executable=comsol,
        repetitions=1,
    )

    assert len(calls) == 2
    assert results[0].solve_time_s == pytest.approx(370.25)
    assert results[0].profile_rows == 2
    provenance = json.loads(
        results[0].paths.provenance_json.read_text(encoding="utf-8")
    )
    assert provenance["status"] == "ok"
    assert provenance["comsol"]["version"] == "6.4"
    assert provenance["comsol"]["build"] == "429"
    assert provenance["conditions"]["voltage_V"] == 200.0
    assert provenance["executed_conditions"]["mesh_elements"] == 200
    assert provenance["artifacts"]["java"]["sha256"]
    assert provenance["artifacts"]["class"]["sha256"]
    assert provenance["artifacts"]["input_mph"]["sha256"]
    assert provenance["artifacts"]["output_mph"]["sha256"]
    assert provenance["artifacts"]["profile_csv"]["sha256"]
    assert provenance["profile_checks"]["all_finite"] is True
    assert (tmp_path / "runs" / "timing_runs.csv").is_file()
    assert (tmp_path / "runs" / "timing_summary.json").is_file()

    with pytest.raises(runner.BuiltinBenchmarkError, match="already exists"):
        runner.run_repetitions(
            input_mph=input_mph,
            output_root=tmp_path / "runs",
            comsol_executable=comsol,
            repetitions=1,
        )


def _write_profile(path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(runner.PROFILE_HEADER)
        writer.writerow(
            [
                0.016,
                1.0e14,
                6.0,
                0.0,
                -0.2,
                -0.8,
                -1.0,
                1.0e19,
                1.0e18,
                4.0e-18,
                200.0,
                13.3322,
            ]
        )
        writer.writerow(
            [
                0.384,
                2.0e14,
                7.0,
                200.0,
                -0.3,
                -0.7,
                -1.0,
                2.0e19,
                2.0e18,
                5.0e-18,
                200.0,
                13.3322,
            ]
        )


def _fake_exe(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"")
    return path
