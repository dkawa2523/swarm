from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import pytest


ROOT = Path(__file__).resolve().parents[1]
RUNNER_PATH = (
    ROOT
    / "reports"
    / "comsol_swarm_benchmark_2026"
    / "repro"
    / "run_external_benchmark.py"
)
SPEC = importlib.util.spec_from_file_location("external_benchmark_runner", RUNNER_PATH)
assert SPEC is not None and SPEC.loader is not None
runner = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = runner
SPEC.loader.exec_module(runner)


def test_external_runner_uses_independent_cli_and_archives_all_stages(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fixture = _write_repo_fixture(tmp_path)
    calls: list[tuple[list[str], dict[str, object]]] = []

    def fake_run(
        args: list[str], **kwargs: object
    ) -> subprocess.CompletedProcess[str]:
        calls.append((args, kwargs))
        run_number = len(calls)
        _write_successful_product_outputs(
            fixture,
            run_number,
            pressure_Pa=float(args[10]),
            gas_temperature_K=float(args[12]),
            mesh_elements=int(args[14]),
        )
        return subprocess.CompletedProcess(
            args,
            0,
            stdout=f"completed external run {run_number}",
            stderr="",
        )

    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    archive_root = fixture["repo"] / "outputs" / "archive"

    results = runner.run_repetitions(
        mapping=fixture["mapping"],
        bundle=fixture["bundle"],
        comsol_executable=fixture["comsol"],
        output_root=archive_root,
        repetitions=2,
        pressure_Pa=15.0,
        gas_temperature_K=310.0,
        mesh_elements=400,
    )

    assert len(calls) == 2
    assert len(results) == 2
    for args, kwargs in calls:
        assert args[:4] == [
            sys.executable,
            "-m",
            "swarm_workflow.cli",
            "run-comsol",
        ]
        assert args[4] == str(fixture["mapping"].resolve())
        assert args[5:7] == ["--bundle", str(fixture["bundle"].resolve())]
        assert args[7:9] == ["--comsol", str(fixture["comsol"].resolve())]
        assert args[9:11] == ["--pressure-Pa", "15"]
        assert args[11:13] == ["--gas-temperature-K", "310"]
        assert args[13:15] == ["--mesh-elements", "400"]
        assert kwargs["cwd"] == str(fixture["repo"].resolve())
        assert kwargs["shell"] is False
        assert kwargs["check"] is False
        assert kwargs["capture_output"] is True
        assert kwargs["text"] is True

    for index, result in enumerate(results, start=1):
        workflow = result.run_dir / "workflow"
        assert (workflow / runner.RESULT_CSV).is_file()
        assert (workflow / runner.RUN_SUMMARY).is_file()
        assert (workflow / runner.EFFECTIVE_MAPPING).is_file()
        for java_name in runner.GENERATED_JAVA:
            assert (workflow / "java" / java_name).is_file()
        for stage in runner.STAGES:
            stage_dir = workflow / "stage_logs" / stage
            assert (stage_dir / "provenance.json").is_file()
            assert (stage_dir / "result.json").is_file()
        assert (
            workflow / "verification" / "function_verify_summary.csv"
        ).is_file()
        provenance = json.loads(
            result.provenance_json.read_text(encoding="utf-8")
        )
        assert provenance["status"] == "ok"
        assert provenance["return_code"] == 0
        assert provenance["inputs"]["mapping"]["sha256"]
        assert provenance["inputs"]["bundle"]["tree_sha256"]
        assert provenance["inputs"]["input_mph"]["sha256"]
        assert provenance["requested_conditions"] == {
            "pressure_Pa": 15.0,
            "gas_temperature_K": 310.0,
            "mesh_elements": 400,
        }
        assert provenance["model_conditions"]["status"] == "passed"
        assert provenance["timing_s"]["apply_s"] == pytest.approx(10.0 + index)
        assert provenance["timing_s"]["solve_s"] == pytest.approx(40.0 + index)
        manifest = json.loads(
            (result.run_dir / "archive_manifest.json").read_text(
                encoding="utf-8"
            )
        )
        assert manifest["copy_policy"].startswith("copy-only")
        assert len(manifest["records"]) == 12

    timing_rows = (archive_root / "timing_runs.csv").read_text(encoding="utf-8")
    assert "apply_s,verify_s,run_s,solve_s,export_s,total_s" in timing_rows
    timing_summary = json.loads(
        (archive_root / "timing_summary.json").read_text(encoding="utf-8")
    )
    assert timing_summary["run_count"] == 2
    assert timing_summary["statistics"]["solve_s"]["median"] == pytest.approx(
        41.5
    )
    aggregate = json.loads(
        (archive_root / "runner_provenance.json").read_text(encoding="utf-8")
    )
    assert aggregate["status"] == "ok"
    assert aggregate["repetitions_completed"] == 2

    # Copy-only archive: source product outputs remain in place.
    assert (
        fixture["repo"]
        / "outputs"
        / "comsol_positive_column"
        / runner.RESULT_CSV
    ).is_file()


def test_external_runner_refuses_existing_archive_root(tmp_path: Path) -> None:
    fixture = _write_repo_fixture(tmp_path)
    archive_root = fixture["repo"] / "outputs" / "archive"
    archive_root.mkdir(parents=True)

    with pytest.raises(runner.ExternalBenchmarkError, match="already exists"):
        runner.run_repetitions(
            mapping=fixture["mapping"],
            bundle=fixture["bundle"],
            comsol_executable=fixture["comsol"],
            output_root=archive_root,
        )


def test_external_runner_records_nonzero_subprocess_failure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fixture = _write_repo_fixture(tmp_path)

    def fake_run(
        args: list[str], **kwargs: object
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(
            args,
            2,
            stdout="partial output",
            stderr="license unavailable",
        )

    monkeypatch.setattr(runner.subprocess, "run", fake_run)
    archive_root = fixture["repo"] / "outputs" / "archive"

    with pytest.raises(runner.ExternalBenchmarkError, match="return code 2"):
        runner.run_repetitions(
            mapping=fixture["mapping"],
            bundle=fixture["bundle"],
            comsol_executable=fixture["comsol"],
            output_root=archive_root,
            repetitions=1,
        )

    run_provenance = json.loads(
        (archive_root / "run_01" / "provenance.json").read_text(
            encoding="utf-8"
        )
    )
    assert run_provenance["status"] == "failed"
    assert run_provenance["return_code"] == 2
    assert (
        archive_root / "run_01" / "workflow_stderr.txt"
    ).read_text(encoding="utf-8") == "license unavailable"
    aggregate = json.loads(
        (archive_root / "runner_provenance.json").read_text(encoding="utf-8")
    )
    assert aggregate["status"] == "failed"
    assert aggregate["repetitions_completed"] == 0


def _write_repo_fixture(tmp_path: Path) -> dict[str, Path]:
    repo = tmp_path / "repo"
    (repo / "Model" / "maps").mkdir(parents=True)
    (repo / "outputs").mkdir()
    (repo / "pyproject.toml").write_text("[project]\nname='fixture'\n", encoding="utf-8")
    input_mph = repo / "Model" / "positive_column_1d.mph"
    input_mph.write_bytes(b"input mph")
    mapping = repo / "Model" / "maps" / "positive_column_external.yaml"
    mapping.write_text(
        "model:\n"
        "  input_mph: Model/positive_column_1d.mph\n"
        "  output_mph: Model/work/output.mph\n",
        encoding="utf-8",
    )
    bundle = repo / "outputs" / "bundle"
    bundle.mkdir()
    (bundle / "manifest.json").write_text('{"status":"ok"}\n', encoding="utf-8")
    (bundle / "table.csv").write_text("x,y\n1,2\n", encoding="utf-8")
    comsol = repo / "comsolbatch.exe"
    comsol.write_bytes(b"fake executable")
    return {
        "repo": repo,
        "input_mph": input_mph,
        "mapping": mapping,
        "bundle": bundle,
        "comsol": comsol,
    }


def _write_successful_product_outputs(
    fixture: dict[str, Path],
    run_number: int,
    *,
    pressure_Pa: float = 13.3322,
    gas_temperature_K: float = 293.15,
    mesh_elements: int = 200,
) -> None:
    repo = fixture["repo"]
    source_output = repo / "outputs" / "comsol_positive_column"
    source_output.mkdir(exist_ok=True)
    result_csv = source_output / runner.RESULT_CSV
    result_csv.write_text("x,electron_density\n0,1\n1,2\n", encoding="utf-8")
    effective_mapping = source_output / runner.EFFECTIVE_MAPPING
    effective_mapping.write_text("effective: true\n", encoding="utf-8")
    for name in runner.GENERATED_JAVA:
        (source_output / name).write_text(
            f"// generated run {run_number}\n",
            encoding="utf-8",
        )

    stage_paths: dict[str, Path] = {}
    provenance_paths: dict[str, str] = {}
    for stage in runner.STAGES:
        stage_dir = repo / "Model" / "logs" / f"{stage}_fixture_{run_number}"
        stage_dir.mkdir(parents=True, exist_ok=True)
        steps = [
            {"name": "compile", "duration_s": 1.0},
            {
                "name": (
                    "positive_column_run"
                    if stage == "run"
                    else f"{stage}_runtime"
                ),
                "duration_s": 40.0 + run_number,
            },
        ]
        (stage_dir / "result.json").write_text(
            json.dumps({"status": "ok", "steps": steps}) + "\n",
            encoding="utf-8",
        )
        provenance = stage_dir / "provenance.json"
        provenance.write_text('{"status":"ok"}\n', encoding="utf-8")
        (stage_dir / f"{stage}_stdout.txt").write_text(
            f"{stage} output\n",
            encoding="utf-8",
        )
        stage_paths[stage] = stage_dir
        provenance_paths[stage] = str(provenance.resolve())

    verify_dir = repo / "outputs" / "comsol_verify"
    verify_dir.mkdir(exist_ok=True)
    verify_csv = verify_dir / "function_verify_summary.csv"
    verify_csv.write_text("quantity,status\nmean_energy,pass\n", encoding="utf-8")

    timing = {
        "apply": 10.0 + run_number,
        "verify": 11.0 + run_number,
        "run": 42.0 + run_number,
        "export": 12.0 + run_number,
    }
    timing["total"] = sum(timing.values())
    summary = {
        "stage": "run-comsol",
        "status": "completed",
        "mapping": str(fixture["mapping"].resolve()),
        "effective_mapping": str(effective_mapping.resolve()),
        "bundle": str(fixture["bundle"].resolve()),
        "input_mph": str(fixture["input_mph"].resolve()),
        "result_csv": str(result_csv.resolve()),
        "logs": {
            "apply": str(stage_paths["apply"].resolve()),
            "verify_summary_csv": str(verify_csv.resolve()),
            "run": str(stage_paths["run"].resolve()),
            "export": str(stage_paths["export"].resolve()),
        },
        "provenance": provenance_paths,
        "timing_s": timing,
        "model_conditions": {
            "status": "passed",
            "requested": {
                "pressure_Pa": pressure_Pa,
                "gas_temperature_K": gas_temperature_K,
                "mesh_elements": mesh_elements,
            },
            "applied": {
                "pressure_Pa": pressure_Pa,
                "gas_temperature_K": gas_temperature_K,
                "mesh_elements": mesh_elements,
            },
            "checks": {
                "apply_feature_readback": {"status": "passed"},
                "pressure_spatial_profile": {"status": "passed"},
                "gas_temperature": {"status": "passed"},
                "mesh_elements": {"status": "passed"},
            },
        },
    }
    (source_output / runner.RUN_SUMMARY).write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )

    # Ensure all mocked product outputs are newer than the runner freshness mark.
    now_ns = time.time_ns() + 10_000_000
    for path in repo.rglob("*"):
        if path.is_file():
            os.utime(path, ns=(now_ns, now_ns))
