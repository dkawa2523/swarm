"""Run and archive independent external-Swarm COMSOL benchmark repetitions.

Each repetition launches a new Python process for the complete product CLI:

``python -m swarm_workflow.cli run-comsol MAPPING --bundle BUNDLE --comsol EXE``

The mutable product output directory is never treated as the record of a run.
After a successful, fresh workflow completion, this runner copies the result,
summary, effective mapping, generated Java, and all four referenced stage-log
directories into a new immutable ``run_NN`` archive.  Existing archive roots
and run directories are refused; nothing is deleted or moved.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
import time
from time import perf_counter
from typing import Any, Sequence

import yaml


RESULT_CSV = "positive_column_results.csv"
RUN_SUMMARY = "positive_column_run_summary.json"
EFFECTIVE_MAPPING = "positive_column_effective_mapping.yaml"
GENERATED_JAVA = (
    "SwarmComsolApply.java",
    "SwarmComsolVerify.java",
    "SwarmPositiveColumnRun.java",
    "SwarmPositiveColumnResultsExport.java",
)
STAGES = ("apply", "verify", "run", "export")
TIMING_KEYS = (
    "apply_s",
    "verify_s",
    "run_s",
    "solve_s",
    "export_s",
    "total_s",
    "workflow_wall_end_to_end_s",
    "runner_end_to_end_s",
)
TIMING_DEFINITIONS = {
    "apply_s": "Product apply-stage total, including Java compile and COMSOL batch.",
    "verify_s": "Product verification-stage total, including Java compile and COMSOL batch.",
    "run_s": "Product positive-column run-stage total, including Java compile and COMSOL batch.",
    "solve_s": (
        "Duration of the positive_column_run COMSOL batch step from the run-stage "
        "result.json; it includes model load, voltage solve/continuation, and "
        "model save, and is not a PDE-kernel-only timer."
    ),
    "export_s": "Product result-export stage total, including Java compile and COMSOL batch.",
    "total_s": "Sum of apply, verify, run, and export product-stage totals.",
    "workflow_wall_end_to_end_s": (
        "Parent wall clock around the independent Python run-comsol subprocess."
    ),
    "runner_end_to_end_s": (
        "Parent wall clock including subprocess execution, validation, hashing, "
        "and copy-only archival."
    ),
}


class ExternalBenchmarkError(RuntimeError):
    """Raised when a clean external-COMSOL repetition cannot be certified."""


@dataclass(frozen=True, slots=True)
class ExternalRunResult:
    index: int
    run_dir: Path
    provenance_json: Path
    return_code: int
    timing_s: dict[str, float]


def run_repetitions(
    *,
    mapping: Path,
    bundle: Path,
    comsol_executable: Path,
    output_root: Path,
    repetitions: int = 3,
    pressure_Pa: float = 13.3322,
    gas_temperature_K: float = 293.15,
    mesh_elements: int = 200,
) -> list[ExternalRunResult]:
    """Run the external path in independent processes and archive each result."""

    mapping = mapping.resolve()
    bundle = bundle.resolve()
    comsol_executable = comsol_executable.resolve()
    output_root = output_root.resolve()
    if repetitions < 1:
        raise ExternalBenchmarkError("repetitions must be at least 1")
    for name, value in (
        ("pressure_Pa", pressure_Pa),
        ("gas_temperature_K", gas_temperature_K),
    ):
        if not math.isfinite(float(value)) or float(value) <= 0.0:
            raise ExternalBenchmarkError(f"{name} must be finite and positive")
    if mesh_elements not in (200, 400):
        raise ExternalBenchmarkError("mesh_elements must be 200 or 400")
    if not mapping.is_file():
        raise ExternalBenchmarkError(f"mapping does not exist: {mapping}")
    if not bundle.is_dir():
        raise ExternalBenchmarkError(f"bundle does not exist: {bundle}")
    if not (bundle / "manifest.json").is_file():
        raise ExternalBenchmarkError(
            f"bundle manifest.json does not exist: {bundle / 'manifest.json'}"
        )
    if not comsol_executable.is_file():
        raise ExternalBenchmarkError(
            f"COMSOL executable does not exist: {comsol_executable}"
        )
    if output_root.exists():
        raise ExternalBenchmarkError(
            "archive output root already exists; refusing reuse: "
            f"{output_root}"
        )

    repo_root = _discover_repo_root(mapping)
    _require_within(repo_root, mapping, "mapping")
    _require_within(repo_root, bundle, "bundle")
    _require_within(repo_root, output_root.parent, "archive output parent")
    input_mph = _mapping_input_mph(mapping, repo_root)
    if not input_mph.is_file():
        raise ExternalBenchmarkError(f"input MPH does not exist: {input_mph}")
    _require_within(repo_root, input_mph, "input MPH")
    source_output = (repo_root / "outputs" / "comsol_positive_column").resolve()

    output_root.mkdir(parents=True)
    static_inputs = {
        "runner": _file_artifact(Path(__file__)),
        "mapping": _file_artifact(mapping),
        "bundle": _tree_artifact(bundle),
        "bundle_manifest": _file_artifact(bundle / "manifest.json"),
        "input_mph": _file_artifact(input_mph),
        "comsol_executable": _file_artifact(comsol_executable),
    }
    requested_conditions = {
        "pressure_Pa": float(pressure_Pa),
        "gas_temperature_K": float(gas_temperature_K),
        "mesh_elements": int(mesh_elements),
    }
    _write_json(
        output_root / "runner_provenance.json",
        {
            "format_version": 1,
            "status": "running",
            "recorded_at_utc": _utc_now(),
            "python_executable": sys.executable,
            "repetitions_requested": repetitions,
            "source_output": str(source_output),
            "inputs": static_inputs,
            "requested_conditions": requested_conditions,
            "timing_definitions": TIMING_DEFINITIONS,
        },
    )

    results: list[ExternalRunResult] = []
    try:
        for index in range(1, repetitions + 1):
            results.append(
                _run_one(
                    index=index,
                    repo_root=repo_root,
                    source_output=source_output,
                    output_root=output_root,
                    mapping=mapping,
                    bundle=bundle,
                    comsol_executable=comsol_executable,
                    input_mph=input_mph,
                    static_inputs=static_inputs,
                    requested_conditions=requested_conditions,
                )
            )
    except Exception as exc:
        _write_json(
            output_root / "runner_provenance.json",
            {
                "format_version": 1,
                "status": "failed",
                "recorded_at_utc": _utc_now(),
                "python_executable": sys.executable,
                "repetitions_requested": repetitions,
                "repetitions_completed": len(results),
                "source_output": str(source_output),
                "inputs": static_inputs,
                "requested_conditions": requested_conditions,
                "timing_definitions": TIMING_DEFINITIONS,
                "failure": f"{type(exc).__name__}: {exc}",
                "runs": [str(result.provenance_json) for result in results],
            },
        )
        if isinstance(exc, ExternalBenchmarkError):
            raise
        raise ExternalBenchmarkError(str(exc)) from exc

    timing_csv, timing_summary = _write_aggregate(output_root, results)
    _write_json(
        output_root / "runner_provenance.json",
        {
            "format_version": 1,
            "status": "ok",
            "recorded_at_utc": _utc_now(),
            "python_executable": sys.executable,
            "repetitions_requested": repetitions,
            "repetitions_completed": len(results),
            "source_output": str(source_output),
            "inputs": static_inputs,
            "requested_conditions": requested_conditions,
            "timing_definitions": TIMING_DEFINITIONS,
            "runs": [str(result.provenance_json) for result in results],
            "timing_runs_csv": str(timing_csv),
            "timing_summary_json": str(timing_summary),
        },
    )
    return results


def _run_one(
    *,
    index: int,
    repo_root: Path,
    source_output: Path,
    output_root: Path,
    mapping: Path,
    bundle: Path,
    comsol_executable: Path,
    input_mph: Path,
    static_inputs: dict[str, Any],
    requested_conditions: dict[str, float | int],
) -> ExternalRunResult:
    runner_started = perf_counter()
    run_dir = output_root / f"run_{index:02d}"
    if run_dir.exists():
        raise ExternalBenchmarkError(
            f"run archive already exists; refusing reuse: {run_dir}"
        )
    run_dir.mkdir()
    command = [
        sys.executable,
        "-m",
        "swarm_workflow.cli",
        "run-comsol",
        str(mapping),
        "--bundle",
        str(bundle),
        "--comsol",
        str(comsol_executable),
        "--pressure-Pa",
        format(float(requested_conditions["pressure_Pa"]), ".17g"),
        "--gas-temperature-K",
        format(float(requested_conditions["gas_temperature_K"]), ".17g"),
        "--mesh-elements",
        str(int(requested_conditions["mesh_elements"])),
    ]
    stdout_path = run_dir / "workflow_stdout.txt"
    stderr_path = run_dir / "workflow_stderr.txt"
    command_json = run_dir / "command.json"
    provenance_json = run_dir / "provenance.json"
    _write_json(
        command_json,
        {
            "format_version": 1,
            "cwd": str(repo_root),
            "args": command,
        },
    )
    provenance: dict[str, Any] = {
        "format_version": 1,
        "status": "running",
        "started_at_utc": _utc_now(),
        "index": index,
        "command": str(command_json),
        "inputs": static_inputs,
        "requested_conditions": requested_conditions,
        "timing_definitions": TIMING_DEFINITIONS,
    }
    _write_json(provenance_json, provenance)

    freshness_epoch_ns = time.time_ns()
    workflow_started = perf_counter()
    try:
        completed = subprocess.run(
            command,
            cwd=str(repo_root),
            shell=False,
            check=False,
            capture_output=True,
            text=True,
        )
        workflow_wall_s = perf_counter() - workflow_started
        stdout = _completed_text(completed.stdout)
        stderr = _completed_text(completed.stderr)
        return_code = int(completed.returncode)
    except OSError as exc:
        workflow_wall_s = perf_counter() - workflow_started
        stdout = ""
        stderr = str(exc)
        return_code = -1
    stdout_path.write_text(stdout, encoding="utf-8")
    stderr_path.write_text(stderr, encoding="utf-8")

    try:
        if return_code != 0:
            raise ExternalBenchmarkError(
                f"external workflow process failed with return code {return_code}; "
                f"stdout={stdout_path}; stderr={stderr_path}"
            )
        summary_path = source_output / RUN_SUMMARY
        _require_fresh_file(summary_path, freshness_epoch_ns, RUN_SUMMARY)
        summary = _read_json(summary_path)
        _validate_summary_contract(
            summary,
            mapping=mapping,
            bundle=bundle,
            input_mph=input_mph,
            requested_conditions=requested_conditions,
        )
        source_files = _fresh_source_files(source_output, freshness_epoch_ns)
        stage_sources, verify_summary = _stage_sources(
            summary,
            repo_root=repo_root,
            freshness_epoch_ns=freshness_epoch_ns,
        )
        stage_timing = _extract_stage_timing(summary, stage_sources["run"])
        stage_timing["workflow_wall_end_to_end_s"] = workflow_wall_s
        archive_manifest = _archive_run(
            run_dir=run_dir,
            source_files=source_files,
            stage_sources=stage_sources,
            verify_summary=verify_summary,
        )
        stage_timing["runner_end_to_end_s"] = perf_counter() - runner_started
        summary_conditions = summary["model_conditions"]
    except Exception as exc:
        provenance.update(
            {
                "status": "failed",
                "recorded_at_utc": _utc_now(),
                "return_code": return_code,
                "workflow_wall_end_to_end_s": workflow_wall_s,
                "failure": f"{type(exc).__name__}: {exc}",
                "artifacts": {
                    "stdout": _file_artifact(stdout_path),
                    "stderr": _file_artifact(stderr_path),
                },
            }
        )
        _write_json(provenance_json, provenance)
        if isinstance(exc, ExternalBenchmarkError):
            raise
        raise ExternalBenchmarkError(str(exc)) from exc

    provenance.update(
        {
            "status": "ok",
            "recorded_at_utc": _utc_now(),
            "return_code": return_code,
            "timing_s": stage_timing,
            "model_conditions": summary_conditions,
            "archive_manifest": str(archive_manifest),
            "artifacts": {
                "stdout": _file_artifact(stdout_path),
                "stderr": _file_artifact(stderr_path),
                "archived_summary": _file_artifact(
                    run_dir / "workflow" / RUN_SUMMARY
                ),
                "archived_result_csv": _file_artifact(
                    run_dir / "workflow" / RESULT_CSV
                ),
                "archive_manifest": _file_artifact(archive_manifest),
            },
        }
    )
    _write_json(provenance_json, provenance)
    return ExternalRunResult(
        index=index,
        run_dir=run_dir,
        provenance_json=provenance_json,
        return_code=return_code,
        timing_s=stage_timing,
    )


def _validate_summary_contract(
    summary: dict[str, Any],
    *,
    mapping: Path,
    bundle: Path,
    input_mph: Path,
    requested_conditions: dict[str, float | int],
) -> None:
    if summary.get("stage") != "run-comsol" or summary.get("status") != "completed":
        raise ExternalBenchmarkError(
            "positive-column summary is not a completed run-comsol result"
        )
    expected = {
        "mapping": mapping,
        "bundle": bundle,
        "input_mph": input_mph,
    }
    for key, requested in expected.items():
        value = summary.get(key)
        if not isinstance(value, str) or Path(value).resolve() != requested:
            raise ExternalBenchmarkError(
                f"summary {key} does not match requested input: {value!r}"
            )
    conditions = summary.get("model_conditions")
    if not isinstance(conditions, dict) or conditions.get("status") != "passed":
        raise ExternalBenchmarkError(
            "summary model_conditions is missing or did not pass"
        )
    recorded = conditions.get("requested")
    if not isinstance(recorded, dict):
        raise ExternalBenchmarkError(
            "summary model_conditions.requested is not a mapping"
        )
    for key, requested in requested_conditions.items():
        actual = recorded.get(key)
        if key == "mesh_elements":
            matches = actual == requested
        else:
            try:
                matches = math.isclose(
                    float(actual),
                    float(requested),
                    rel_tol=1.0e-12,
                    abs_tol=1.0e-12,
                )
            except (TypeError, ValueError):
                matches = False
        if not matches:
            raise ExternalBenchmarkError(
                f"summary model condition {key}={actual!r} does not match "
                f"requested {requested!r}"
            )
    checks = conditions.get("checks")
    if not isinstance(checks, dict):
        raise ExternalBenchmarkError(
            "summary model_conditions.checks is not a mapping"
        )
    required_checks = (
        "apply_feature_readback",
        "pressure_spatial_profile",
        "gas_temperature",
        "mesh_elements",
    )
    failed_checks = [
        key
        for key in required_checks
        if not isinstance(checks.get(key), dict)
        or checks[key].get("status") != "passed"
    ]
    if failed_checks:
        raise ExternalBenchmarkError(
            "summary model-condition verification failed or is missing: "
            + ", ".join(failed_checks)
        )


def _fresh_source_files(
    source_output: Path,
    freshness_epoch_ns: int,
) -> list[Path]:
    required = [
        source_output / RESULT_CSV,
        source_output / RUN_SUMMARY,
        source_output / EFFECTIVE_MAPPING,
        *(source_output / name for name in GENERATED_JAVA),
    ]
    for path in required:
        _require_fresh_file(path, freshness_epoch_ns, path.name)
    return required


def _stage_sources(
    summary: dict[str, Any],
    *,
    repo_root: Path,
    freshness_epoch_ns: int,
) -> tuple[dict[str, Path], Path]:
    provenance = summary.get("provenance")
    logs = summary.get("logs")
    if not isinstance(provenance, dict) or not isinstance(logs, dict):
        raise ExternalBenchmarkError(
            "summary must contain logs and provenance mappings"
        )
    stage_sources: dict[str, Path] = {}
    for stage in STAGES:
        value = provenance.get(stage)
        if not isinstance(value, str):
            raise ExternalBenchmarkError(
                f"summary provenance is missing stage {stage}"
            )
        provenance_path = Path(value).resolve()
        _require_within(repo_root, provenance_path, f"{stage} provenance")
        _require_fresh_file(
            provenance_path,
            freshness_epoch_ns,
            f"{stage} provenance",
        )
        provenance_payload = _read_json(provenance_path)
        if provenance_payload.get("status") != "ok":
            raise ExternalBenchmarkError(
                f"{stage} provenance status is not ok: {provenance_path}"
            )
        stage_dir = provenance_path.parent
        if not stage_dir.is_dir():
            raise ExternalBenchmarkError(
                f"{stage} log directory does not exist: {stage_dir}"
            )
        stage_sources[stage] = stage_dir
        log_value = logs.get(stage)
        if stage != "verify":
            if not isinstance(log_value, str):
                raise ExternalBenchmarkError(
                    f"summary logs is missing stage {stage}"
                )
            if Path(log_value).resolve() != stage_dir:
                raise ExternalBenchmarkError(
                    f"{stage} log directory and provenance parent disagree"
                )

    verify_value = logs.get("verify_summary_csv")
    if not isinstance(verify_value, str):
        raise ExternalBenchmarkError("summary is missing verify_summary_csv")
    verify_summary = Path(verify_value).resolve()
    _require_within(repo_root, verify_summary, "verify summary")
    _require_fresh_file(
        verify_summary,
        freshness_epoch_ns,
        "verify summary CSV",
    )
    return stage_sources, verify_summary


def _extract_stage_timing(
    summary: dict[str, Any],
    run_log_dir: Path,
) -> dict[str, float]:
    timing = summary.get("timing_s")
    if not isinstance(timing, dict):
        raise ExternalBenchmarkError("summary timing_s is not a mapping")
    result = {
        f"{stage}_s": _finite_nonnegative(timing.get(stage), stage)
        for stage in ("apply", "verify", "run")
    }
    solve_value = timing.get("solve")
    if solve_value is None:
        run_result_path = run_log_dir / "result.json"
        if not run_result_path.is_file():
            raise ExternalBenchmarkError(
                f"run stage result.json is missing: {run_result_path}"
            )
        run_result = _read_json(run_result_path)
        steps = run_result.get("steps")
        if not isinstance(steps, list):
            raise ExternalBenchmarkError(
                f"run stage result has no steps: {run_result_path}"
            )
        runtime_step = next(
            (
                step
                for step in steps
                if isinstance(step, dict)
                and step.get("name") in ("positive_column_run", "solve")
            ),
            None,
        )
        if runtime_step is None:
            raise ExternalBenchmarkError(
                "run stage result has no positive_column_run timing"
            )
        solve_value = runtime_step.get("duration_s")
    result["solve_s"] = _finite_nonnegative(solve_value, "solve")
    result["export_s"] = _finite_nonnegative(timing.get("export"), "export")
    result["total_s"] = _finite_nonnegative(timing.get("total"), "total")
    return result


def _archive_run(
    *,
    run_dir: Path,
    source_files: Sequence[Path],
    stage_sources: dict[str, Path],
    verify_summary: Path,
) -> Path:
    workflow_dir = run_dir / "workflow"
    java_dir = workflow_dir / "java"
    stage_root = workflow_dir / "stage_logs"
    verify_dir = workflow_dir / "verification"
    workflow_dir.mkdir()
    java_dir.mkdir()
    stage_root.mkdir()
    verify_dir.mkdir()
    records: list[dict[str, Any]] = []
    for source in source_files:
        destination = (
            java_dir / source.name
            if source.suffix.lower() == ".java"
            else workflow_dir / source.name
        )
        if destination.exists():
            raise ExternalBenchmarkError(
                f"archive destination already exists: {destination}"
            )
        shutil.copy2(source, destination)
        records.append(
            {
                "kind": "file",
                "source": str(source),
                "archive": str(destination.relative_to(run_dir)),
                "artifact": _file_artifact(destination),
            }
        )
    for stage in STAGES:
        source = stage_sources[stage]
        destination = stage_root / stage
        if destination.exists():
            raise ExternalBenchmarkError(
                f"archive stage destination already exists: {destination}"
            )
        shutil.copytree(source, destination)
        records.append(
            {
                "kind": "tree",
                "stage": stage,
                "source": str(source),
                "archive": str(destination.relative_to(run_dir)),
                "artifact": _tree_artifact(destination),
            }
        )
    verify_destination = verify_dir / verify_summary.name
    shutil.copy2(verify_summary, verify_destination)
    records.append(
        {
            "kind": "file",
            "source": str(verify_summary),
            "archive": str(verify_destination.relative_to(run_dir)),
            "artifact": _file_artifact(verify_destination),
        }
    )
    manifest_path = run_dir / "archive_manifest.json"
    _write_json(
        manifest_path,
        {
            "format_version": 1,
            "status": "ok",
            "copy_policy": "copy-only; source files were not moved or deleted",
            "recorded_at_utc": _utc_now(),
            "records": records,
        },
    )
    return manifest_path


def _write_aggregate(
    output_root: Path,
    results: Sequence[ExternalRunResult],
) -> tuple[Path, Path]:
    rows: list[dict[str, Any]] = []
    for result in results:
        row: dict[str, Any] = {
            "run": result.index,
            **result.timing_s,
            "return_code": result.return_code,
            "provenance": str(result.provenance_json),
        }
        rows.append(row)
    csv_path = output_root / "timing_runs.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary_path = output_root / "timing_summary.json"
    _write_json(
        summary_path,
        {
            "format_version": 1,
            "status": "ok",
            "run_count": len(rows),
            "runs_csv": str(csv_path),
            "statistics": {
                key: {
                    "median": statistics.median(
                        float(row[key]) for row in rows
                    ),
                    "min": min(float(row[key]) for row in rows),
                    "max": max(float(row[key]) for row in rows),
                }
                for key in TIMING_KEYS
            },
            "timing_definitions": TIMING_DEFINITIONS,
        },
    )
    return csv_path, summary_path


def _mapping_input_mph(mapping: Path, repo_root: Path) -> Path:
    raw = yaml.safe_load(mapping.read_text(encoding="utf-8")) or {}
    if not isinstance(raw, dict):
        raise ExternalBenchmarkError("mapping YAML root is not a mapping")
    model = raw.get("model")
    if not isinstance(model, dict):
        raise ExternalBenchmarkError("mapping YAML has no model mapping")
    value = model.get("input_mph")
    if not isinstance(value, str) or not value.strip():
        raise ExternalBenchmarkError("mapping model.input_mph is missing")
    path = Path(value)
    return (path if path.is_absolute() else repo_root / path).resolve()


def _discover_repo_root(path: Path) -> Path:
    for candidate in (path.parent, *path.parents):
        if (candidate / "pyproject.toml").is_file():
            return candidate.resolve()
    raise ExternalBenchmarkError(
        f"repository root with pyproject.toml not found above {path}"
    )


def _require_within(root: Path, path: Path, label: str) -> None:
    try:
        path.resolve().relative_to(root.resolve())
    except ValueError as exc:
        raise ExternalBenchmarkError(
            f"{label} must stay inside repository root {root}: {path}"
        ) from exc


def _require_fresh_file(path: Path, epoch_ns: int, label: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise ExternalBenchmarkError(f"{label} is missing or empty: {path}")
    if path.stat().st_mtime_ns < epoch_ns:
        raise ExternalBenchmarkError(
            f"{label} predates this repetition and may be stale: {path}"
        )


def _finite_nonnegative(value: Any, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ExternalBenchmarkError(
            f"timing {label} is not numeric: {value!r}"
        ) from exc
    if not math.isfinite(number) or number < 0.0:
        raise ExternalBenchmarkError(
            f"timing {label} must be finite and nonnegative: {number}"
        )
    return number


def _tree_artifact(root: Path) -> dict[str, Any]:
    resolved = root.resolve()
    if not resolved.is_dir():
        return {"path": str(resolved), "exists": False}
    digest = hashlib.sha256()
    files = sorted(path for path in resolved.rglob("*") if path.is_file())
    total_bytes = 0
    for path in files:
        relative = path.relative_to(resolved).as_posix().encode("utf-8")
        digest.update(relative)
        digest.update(b"\0")
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
                total_bytes += len(chunk)
        digest.update(b"\0")
    return {
        "path": str(resolved),
        "exists": True,
        "file_count": len(files),
        "size_bytes": total_bytes,
        "tree_sha256": digest.hexdigest(),
    }


def _file_artifact(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    artifact: dict[str, Any] = {
        "path": str(resolved),
        "exists": resolved.is_file(),
    }
    if not resolved.is_file():
        return artifact
    stat = resolved.stat()
    artifact.update(
        {
            "size_bytes": stat.st_size,
            "modified_at_utc": datetime.fromtimestamp(
                stat.st_mtime, timezone.utc
            ).isoformat(),
            "sha256": _sha256_file(resolved),
        }
    )
    return artifact


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ExternalBenchmarkError(f"could not read JSON {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise ExternalBenchmarkError(f"JSON root is not a mapping: {path}")
    return value


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _completed_text(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run the complete external-Swarm COMSOL path in independent "
            "processes and archive each repetition."
        )
    )
    parser.add_argument(
        "--mapping",
        type=Path,
        default=Path("Model/maps/positive_column_external.yaml"),
    )
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--comsol", type=Path, required=True)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path(
            "outputs/comsol_swarm_benchmark_2026/external_200elem_runs"
        ),
    )
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--pressure-Pa", type=float, default=13.3322)
    parser.add_argument("--gas-temperature-K", type=float, default=293.15)
    parser.add_argument(
        "--mesh-elements",
        type=int,
        choices=(200, 400),
        default=200,
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    results = run_repetitions(
        mapping=args.mapping,
        bundle=args.bundle,
        comsol_executable=args.comsol,
        output_root=args.output_root,
        repetitions=args.repetitions,
        pressure_Pa=args.pressure_Pa,
        gas_temperature_K=args.gas_temperature_K,
        mesh_elements=args.mesh_elements,
    )
    for result in results:
        print(
            f"run_{result.index:02d}: "
            f"solve={result.timing_s['solve_s']:.6g} s, "
            "wall_end_to_end="
            f"{result.timing_s['workflow_wall_end_to_end_s']:.6g} s, "
            f"provenance={result.provenance_json}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
