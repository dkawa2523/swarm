"""Reproducible COMSOL built-in Boltzmann benchmark runner.

This runner is intentionally separate from the external-table mapping path:
the built-in reference has no Swarm bundle.  It reuses the public COMSOL
command builders and executable resolver from ``swarm_workflow.comsol_adapter``
while recording a report-specific provenance document.

No run directory is ever reused or deleted.  In particular, a pre-existing
class file is a hard error and a compile step that exits with code zero but
prints a COMSOL compilation-error marker is treated as failed.
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
import re
import statistics
import subprocess
from time import perf_counter
from typing import Any, Sequence

from swarm_workflow.comsol_adapter import (
    build_apply_java_batch_command,
    build_java_compile_command,
    resolve_comsol_executable,
)


CLASS_NAME = "RunBuiltinPositiveColumn"
PROFILE_HEADER = (
    "x",
    "electron_density",
    "mean_electron_energy",
    "electric_potential",
    "electron_current_density",
    "ion_current_density",
    "total_current_density",
    "excitation_source",
    "ionization_source",
    "E_over_N",
    "applied_voltage",
    "gas_pressure",
)
ERROR_SIGNATURES = (
    ("compiler_failed", "failed to compile java file"),
    ("compilation_failed", "compilation failed"),
    ("compilation_error", "compilation error on line"),
    ("java_class_error", "error running java class"),
    ("comsol_error_block", "/*****error"),
    ("undefined_variable", "undefined variable"),
    ("unknown_function", "unknown function or operator"),
)
VOLTAGE_MARKER = "SWARM_BUILTIN_VOLTAGE_V="
GAS_TEMPERATURE_MARKER = "SWARM_BUILTIN_GAS_TEMPERATURE_K="
MESH_MARKER = "SWARM_BUILTIN_MESH_ELEMENTS="
SOLVE_TIME_MARKER = "SWARM_BUILTIN_SOLVE_TIME_S="
PROFILE_ROWS_MARKER = "SWARM_BUILTIN_PROFILE_ROWS="


class BuiltinBenchmarkError(RuntimeError):
    """Raised when a built-in COMSOL benchmark run is not reproducible."""


@dataclass(frozen=True, slots=True)
class RunPaths:
    run_dir: Path
    java: Path
    class_file: Path
    output_mph: Path
    profile_csv: Path
    command_json: Path
    provenance_json: Path


@dataclass(frozen=True, slots=True)
class BuiltinRunResult:
    index: int
    paths: RunPaths
    compile_time_s: float
    batch_time_s: float
    solve_time_s: float
    end_to_end_time_s: float
    profile_rows: int
    comsol_version: str | None
    comsol_build: str | None


def render_java_source(
    *,
    input_mph: Path,
    output_mph: Path,
    profile_csv: Path,
    voltage_V: float,
    pressure_Pa: float,
    gas_temperature_K: float,
    mesh_elements: int,
) -> str:
    """Render the fixed 1D built-in-Boltzmann solve and 12-column export."""

    if mesh_elements not in (200, 400):
        raise BuiltinBenchmarkError("mesh_elements must be 200 or 400")
    input_java = _java_path(input_mph)
    output_java = _java_path(output_mph)
    profile_java = _java_path(profile_csv)
    voltage = _java_number(voltage_V)
    pressure = _java_number(pressure_Pa)
    gas_temperature = _java_number(gas_temperature_K)
    return f"""import com.comsol.model.*;
import com.comsol.model.util.*;
import java.io.PrintWriter;

public class {CLASS_NAME} {{
  public static void main(String[] args) throws Exception {{
    Model model = ModelUtil.loadCopy(
      "builtinPositiveColumnBenchmark",
      "{input_java}"
    );
    model.param().set("V0", "{voltage}[V]");
    model.param().set("p0", "{pressure}[Pa]");
    model.component("comp1").physics("plas").feature("pes1")
      .set("T_src", "userdef");
    model.component("comp1").physics("plas").feature("pes1")
      .set("T", "{gas_temperature}[K]");
    model.component("comp1").physics("plas").feature("pes1")
      .set("pA_src", "userdef");
    model.component("comp1").physics("plas").feature("pes1")
      .set("pA", "p0");

    // The source MPH records mesh1/edg1/dis1 with property elemcount=200.
    // Rebuilding the same feature makes both the 200- and 400-element runs
    // explicit and invalidates any saved solution before sol3 is rerun.
    model.component("comp1").mesh("mesh1").feature("edg1")
      .feature("dis1").set("elemcount", {mesh_elements});
    model.component("comp1").mesh("mesh1").run();

    System.out.println("{VOLTAGE_MARKER}{voltage}");
    System.out.println("{GAS_TEMPERATURE_MARKER}{gas_temperature}");
    System.out.println("{MESH_MARKER}{mesh_elements}");
    long solveStarted = System.nanoTime();
    model.sol("sol3").runAll();
    double solveTimeSeconds = (System.nanoTime() - solveStarted) / 1.0e9;
    System.out.println("{SOLVE_TIME_MARKER}" + Double.toString(solveTimeSeconds));
    model.save("{output_java}");

    String evalTag = "swarmBuiltinBenchmarkEval";
    try {{
      model.result().numerical().remove(evalTag);
    }} catch (Exception ignored) {{
    }}
    model.result().numerical().create(evalTag, "Eval");
    model.result().numerical(evalTag).set("expr", new String[]{{
      "x",
      "plas.ne",
      "plas.ebar/1[V/eV]",
      "V",
      "plas.Jelx",
      "plas.Jix_wAr_1p",
      "plas.Jix_wAr_1p+plas.Jelx",
      "root.comp1.plas.eir2.alpha(plas.ebar)*(p0/(k_B_const*{gas_temperature}[K]))*abs(plas.Jelx)/e_const",
      "root.comp1.plas.eir4.alpha(plas.ebar)*(p0/(k_B_const*{gas_temperature}[K]))*abs(plas.Jelx)/e_const",
      "plas.Erd",
      "V0",
      "p0"
    }});
    model.result().numerical(evalTag).set("unit", new String[]{{
      "m", "1/m^3", "eV", "V", "A/m^2", "A/m^2", "A/m^2",
      "1/(m^3*s)", "1/(m^3*s)", "V*m^2", "V", "Pa"
    }});
    double[][][] values = model.result().numerical(evalTag).getData();
    int n = model.result().numerical(evalTag).getNData();
    try (PrintWriter writer = new PrintWriter("{profile_java}", "UTF-8")) {{
      writer.println("{','.join(PROFILE_HEADER)}");
      for (int i = 0; i < n; i++) {{
        StringBuilder row = new StringBuilder();
        for (int j = 0; j < {len(PROFILE_HEADER)}; j++) {{
          if (j > 0) row.append(',');
          double[][] expressionValues = values[j];
          double value = expressionValues[expressionValues.length - 1][i];
          row.append(Double.toString(value));
        }}
        writer.println(row.toString());
      }}
    }}
    System.out.println("{PROFILE_ROWS_MARKER}" + Integer.toString(n));
    ModelUtil.remove(model.tag());
  }}
}}
"""


def run_repetitions(
    *,
    input_mph: Path,
    output_root: Path,
    comsol_executable: str | Path | None,
    repetitions: int = 3,
    start_index: int = 1,
    voltage_V: float = 200.0,
    pressure_Pa: float = 13.3322,
    gas_temperature_K: float = 293.15,
    mesh_elements: int = 200,
) -> list[BuiltinRunResult]:
    """Execute independent compile/batch processes and write timing summaries."""

    input_mph = input_mph.resolve()
    output_root = output_root.resolve()
    if not input_mph.is_file():
        raise BuiltinBenchmarkError(f"input MPH does not exist: {input_mph}")
    if repetitions < 1:
        raise BuiltinBenchmarkError("repetitions must be at least 1")
    if start_index < 1:
        raise BuiltinBenchmarkError("start_index must be at least 1")
    if mesh_elements not in (200, 400):
        raise BuiltinBenchmarkError("mesh_elements must be 200 or 400")
    resolved = resolve_comsol_executable(comsol_executable)
    target_indices = range(start_index, start_index + repetitions)
    existing = [
        output_root / f"run_{index:02d}"
        for index in target_indices
        if (output_root / f"run_{index:02d}").exists()
    ]
    if existing:
        raise BuiltinBenchmarkError(
            "run directory already exists; refusing reuse so a stale class "
            f"cannot be executed: {existing[0]}"
        )
    results: list[BuiltinRunResult] = []
    for index in target_indices:
        results.append(
            _run_one(
                index=index,
                input_mph=input_mph,
                output_root=output_root,
                comsol_executable=resolved.path,
                executable_source=resolved.source,
                voltage_V=voltage_V,
                pressure_Pa=pressure_Pa,
                gas_temperature_K=gas_temperature_K,
                mesh_elements=mesh_elements,
            )
        )
    _write_aggregate(output_root, results)
    return results


def _run_one(
    *,
    index: int,
    input_mph: Path,
    output_root: Path,
    comsol_executable: str,
    executable_source: str,
    voltage_V: float,
    pressure_Pa: float,
    gas_temperature_K: float,
    mesh_elements: int,
) -> BuiltinRunResult:
    end_to_end_started = perf_counter()
    paths = _prepare_run_paths(output_root, index)
    paths.run_dir.mkdir(parents=True)
    source = render_java_source(
        input_mph=input_mph,
        output_mph=paths.output_mph,
        profile_csv=paths.profile_csv,
        voltage_V=voltage_V,
        pressure_Pa=pressure_Pa,
        gas_temperature_K=gas_temperature_K,
        mesh_elements=mesh_elements,
    )
    paths.java.write_text(source, encoding="utf-8", newline="\n")
    compile_args = build_java_compile_command(comsol_executable, paths.java)
    batch_args = build_apply_java_batch_command(comsol_executable, paths.class_file)
    command_payload = {
        "format_version": 1,
        "operation": "builtin_positive_column_benchmark",
        "cwd": str(paths.run_dir),
        "input_mph": str(input_mph),
        "output_mph": str(paths.output_mph),
        "profile_csv": str(paths.profile_csv),
        "conditions": {
            "voltage_V": voltage_V,
            "pressure_Pa": pressure_Pa,
            "gas_temperature_K": gas_temperature_K,
            "mesh_elements": mesh_elements,
        },
        "mesh_contract": {
            "mesh": "comp1/mesh1",
            "edge_feature": "edg1",
            "distribution_feature": "dis1",
            "element_property": "elemcount",
        },
        "solver": "sol3",
        "commands": {
            "compile": compile_args,
            "batch": batch_args,
        },
    }
    _write_json(paths.command_json, command_payload)
    provenance: dict[str, Any] = {
        **command_payload,
        "status": "running",
        "started_at_utc": _utc_now(),
        "comsol": {
            "executable": comsol_executable,
            "executable_source": executable_source,
            "version": None,
            "build": None,
        },
        "artifacts": {
            "runner": _file_artifact(Path(__file__)),
            "input_mph": _file_artifact(input_mph),
            "java": _file_artifact(paths.java),
            "class": _file_artifact(paths.class_file),
            "output_mph": _file_artifact(paths.output_mph),
            "profile_csv": _file_artifact(paths.profile_csv),
        },
    }
    _write_json(paths.provenance_json, provenance)
    compile_result: dict[str, Any] | None = None
    batch_result: dict[str, Any] | None = None
    try:
        if paths.class_file.exists():
            raise BuiltinBenchmarkError(
                f"refusing pre-existing class file: {paths.class_file}"
            )
        compile_result = _run_logged(
            name="compile",
            args=compile_args,
            cwd=paths.run_dir,
        )
        _require_success(compile_result)
        _require_fresh_class(paths.java, paths.class_file)
        batch_result = _run_logged(
            name="batch",
            args=batch_args,
            cwd=paths.run_dir,
        )
        _require_success(batch_result)
        if not paths.output_mph.is_file() or paths.output_mph.stat().st_size == 0:
            raise BuiltinBenchmarkError(
                f"batch did not create a nonempty output MPH: {paths.output_mph}"
            )
        profile_checks = _validate_profile(
            paths.profile_csv,
            voltage_V=voltage_V,
            pressure_Pa=pressure_Pa,
        )
        runtime_text = _read_step_text(batch_result)
        executed_voltage = _parse_float_marker(runtime_text, VOLTAGE_MARKER)
        executed_temperature = _parse_float_marker(
            runtime_text, GAS_TEMPERATURE_MARKER
        )
        executed_mesh = _parse_int_marker(runtime_text, MESH_MARKER)
        solve_time_s = _parse_float_marker(runtime_text, SOLVE_TIME_MARKER)
        reported_rows = _parse_int_marker(runtime_text, PROFILE_ROWS_MARKER)
        if not math.isclose(executed_voltage, voltage_V, rel_tol=0.0, abs_tol=1e-12):
            raise BuiltinBenchmarkError(
                f"executed voltage marker {executed_voltage} != {voltage_V}"
            )
        if executed_mesh != mesh_elements:
            raise BuiltinBenchmarkError(
                f"executed mesh marker {executed_mesh} != {mesh_elements}"
            )
        if not math.isclose(
            executed_temperature,
            gas_temperature_K,
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            raise BuiltinBenchmarkError(
                "executed gas-temperature marker "
                f"{executed_temperature} != {gas_temperature_K}"
            )
        if reported_rows != profile_checks["rows"]:
            raise BuiltinBenchmarkError(
                f"profile marker rows {reported_rows} != CSV rows "
                f"{profile_checks['rows']}"
            )
        identity = _extract_comsol_identity(
            _read_step_text(compile_result) + "\n" + runtime_text
        )
        if not identity.get("version") or not identity.get("build"):
            raise BuiltinBenchmarkError(
                "COMSOL version/build banner was not found in raw logs"
            )
    except Exception as exc:
        provenance.update(
            {
                "status": "failed",
                "recorded_at_utc": _utc_now(),
                "failure": f"{type(exc).__name__}: {exc}",
                "timing": {
                    "compile_s": (
                        compile_result.get("duration_s")
                        if compile_result is not None
                        else None
                    ),
                    "batch_s": (
                        batch_result.get("duration_s")
                        if batch_result is not None
                        else None
                    ),
                    "end_to_end_s": perf_counter() - end_to_end_started,
                },
                "steps": [
                    result
                    for result in (compile_result, batch_result)
                    if result is not None
                ],
                "artifacts": {
                    **provenance["artifacts"],
                    "class": _file_artifact(paths.class_file),
                    "output_mph": _file_artifact(paths.output_mph),
                    "profile_csv": _file_artifact(paths.profile_csv),
                },
            }
        )
        _write_json(paths.provenance_json, provenance)
        if isinstance(exc, BuiltinBenchmarkError):
            raise
        raise BuiltinBenchmarkError(str(exc)) from exc

    end_to_end_s = perf_counter() - end_to_end_started
    provenance.update(
        {
            "status": "ok",
            "recorded_at_utc": _utc_now(),
            "comsol": {
                **provenance["comsol"],
                "version": identity.get("version"),
                "build": identity.get("build"),
                "progress_version": identity.get("progress_version"),
            },
            "executed_conditions": {
                "voltage_V": executed_voltage,
                "gas_temperature_K": executed_temperature,
                "mesh_elements": executed_mesh,
                "profile_voltage_constant": profile_checks[
                    "applied_voltage_constant"
                ],
                "profile_pressure_constant": profile_checks[
                    "gas_pressure_constant"
                ],
            },
            "profile_checks": profile_checks,
            "timing": {
                "compile_s": compile_result["duration_s"],
                "batch_s": batch_result["duration_s"],
                "solve_s": solve_time_s,
                "end_to_end_s": end_to_end_s,
            },
            "steps": [compile_result, batch_result],
            "artifacts": {
                **provenance["artifacts"],
                "class": _file_artifact(paths.class_file),
                "output_mph": _file_artifact(paths.output_mph),
                "profile_csv": _file_artifact(paths.profile_csv),
                "compile_stdout": _file_artifact(
                    Path(compile_result["stdout"])
                ),
                "compile_stderr": _file_artifact(
                    Path(compile_result["stderr"])
                ),
                "batch_stdout": _file_artifact(Path(batch_result["stdout"])),
                "batch_stderr": _file_artifact(Path(batch_result["stderr"])),
            },
        }
    )
    _write_json(paths.provenance_json, provenance)
    return BuiltinRunResult(
        index=index,
        paths=paths,
        compile_time_s=float(compile_result["duration_s"]),
        batch_time_s=float(batch_result["duration_s"]),
        solve_time_s=solve_time_s,
        end_to_end_time_s=end_to_end_s,
        profile_rows=int(profile_checks["rows"]),
        comsol_version=identity.get("version"),
        comsol_build=identity.get("build"),
    )


def _prepare_run_paths(output_root: Path, index: int) -> RunPaths:
    run_dir = output_root / f"run_{index:02d}"
    if run_dir.exists():
        raise BuiltinBenchmarkError(
            "run directory already exists; refusing reuse so a stale class "
            f"cannot be executed: {run_dir}"
        )
    return RunPaths(
        run_dir=run_dir,
        java=run_dir / f"{CLASS_NAME}.java",
        class_file=run_dir / f"{CLASS_NAME}.class",
        output_mph=run_dir / "positive_column_builtin_solved.mph",
        profile_csv=run_dir / "profile.csv",
        command_json=run_dir / "command.json",
        provenance_json=run_dir / "provenance.json",
    )


def _run_logged(*, name: str, args: Sequence[str], cwd: Path) -> dict[str, Any]:
    stdout_path = cwd / f"{name}_stdout.txt"
    stderr_path = cwd / f"{name}_stderr.txt"
    started = perf_counter()
    try:
        completed = subprocess.run(
            list(args),
            cwd=str(cwd),
            shell=False,
            check=False,
            capture_output=True,
            text=True,
        )
        return_code = int(completed.returncode)
        stdout = _completed_text(completed.stdout)
        stderr = _completed_text(completed.stderr)
    except OSError as exc:
        return_code = -1
        stdout = ""
        stderr = str(exc)
    stdout_path.write_text(stdout, encoding="utf-8")
    stderr_path.write_text(stderr, encoding="utf-8")
    detected = _error_signature(stdout, stderr)
    return {
        "name": name,
        "status": "ok" if return_code == 0 and detected is None else "failed",
        "args": list(args),
        "cwd": str(cwd),
        "return_code": return_code,
        "detected_comsol_error": detected,
        "duration_s": perf_counter() - started,
        "stdout": str(stdout_path),
        "stderr": str(stderr_path),
    }


def _require_success(result: dict[str, Any]) -> None:
    if result["status"] == "ok":
        return
    raise BuiltinBenchmarkError(
        f"COMSOL {result['name']} failed: return_code={result['return_code']}, "
        f"detected_error={result['detected_comsol_error']}; "
        f"stdout={result['stdout']}; stderr={result['stderr']}"
    )


def _require_fresh_class(java: Path, class_file: Path) -> None:
    if not class_file.is_file() or class_file.stat().st_size == 0:
        raise BuiltinBenchmarkError(
            f"COMSOL compile did not create a nonempty class: {class_file}"
        )
    if class_file.stat().st_mtime_ns < java.stat().st_mtime_ns:
        raise BuiltinBenchmarkError(
            f"compiled class is older than its Java source: {class_file}"
        )


def _validate_profile(
    path: Path,
    *,
    voltage_V: float,
    pressure_Pa: float,
) -> dict[str, Any]:
    if not path.is_file() or path.stat().st_size == 0:
        raise BuiltinBenchmarkError(f"profile CSV was not created: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != PROFILE_HEADER:
            raise BuiltinBenchmarkError(
                f"profile columns {reader.fieldnames} != {list(PROFILE_HEADER)}"
            )
        rows = list(reader)
    if not rows:
        raise BuiltinBenchmarkError("profile CSV contains no data rows")
    previous_x = -math.inf
    voltages: list[float] = []
    pressures: list[float] = []
    negative_counts = {
        "electron_density": 0,
        "excitation_source": 0,
        "ionization_source": 0,
    }
    for row_index, row in enumerate(rows, start=2):
        values: dict[str, float] = {}
        for column in PROFILE_HEADER:
            try:
                value = float(row[column])
            except (TypeError, ValueError) as exc:
                raise BuiltinBenchmarkError(
                    f"profile row {row_index} column {column} is not numeric"
                ) from exc
            if not math.isfinite(value):
                raise BuiltinBenchmarkError(
                    f"profile row {row_index} column {column} is non-finite"
                )
            values[column] = value
        if values["x"] < previous_x:
            raise BuiltinBenchmarkError("profile x coordinate is not monotonic")
        previous_x = values["x"]
        voltages.append(values["applied_voltage"])
        pressures.append(values["gas_pressure"])
        for column in negative_counts:
            if values[column] < 0.0:
                negative_counts[column] += 1
    voltage_constant = max(voltages) - min(voltages) <= 1e-10
    pressure_constant = max(pressures) - min(pressures) <= 1e-8
    if not voltage_constant or not math.isclose(
        voltages[0], voltage_V, rel_tol=0.0, abs_tol=1e-9
    ):
        raise BuiltinBenchmarkError("profile applied voltage is not the request")
    if not pressure_constant or not math.isclose(
        pressures[0], pressure_Pa, rel_tol=1e-7, abs_tol=1e-8
    ):
        raise BuiltinBenchmarkError("profile gas pressure is not the request")
    if any(negative_counts.values()):
        raise BuiltinBenchmarkError(
            f"profile contains negative density/source values: {negative_counts}"
        )
    return {
        "rows": len(rows),
        "columns": list(PROFILE_HEADER),
        "all_finite": True,
        "x_monotonic_non_decreasing": True,
        "applied_voltage_constant": voltage_constant,
        "applied_voltage_V": voltages[0],
        "gas_pressure_constant": pressure_constant,
        "gas_pressure_Pa": pressures[0],
        "negative_counts": negative_counts,
    }


def _write_aggregate(output_root: Path, results: Sequence[BuiltinRunResult]) -> None:
    rows = [
        {
            "run": result.index,
            "compile_s": result.compile_time_s,
            "batch_s": result.batch_time_s,
            "solve_s": result.solve_time_s,
            "end_to_end_s": result.end_to_end_time_s,
            "profile_rows": result.profile_rows,
            "comsol_version": result.comsol_version or "",
            "comsol_build": result.comsol_build or "",
            "provenance": str(result.paths.provenance_json),
        }
        for result in results
    ]
    output_root.mkdir(parents=True, exist_ok=True)
    csv_path = output_root / "timing_runs.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    timing_keys = ("compile_s", "batch_s", "solve_s", "end_to_end_s")
    summary = {
        "format_version": 1,
        "status": "ok",
        "run_count": len(results),
        "runs_csv": str(csv_path),
        "statistics": {
            key: {
                "median": statistics.median(float(row[key]) for row in rows),
                "min": min(float(row[key]) for row in rows),
                "max": max(float(row[key]) for row in rows),
            }
            for key in timing_keys
        },
        "provenance": [str(result.paths.provenance_json) for result in results],
    }
    _write_json(output_root / "timing_summary.json", summary)


def _parse_float_marker(text: str, marker: str) -> float:
    match = re.search(re.escape(marker) + r"([0-9eE.+-]+)", text)
    if match is None:
        raise BuiltinBenchmarkError(f"required runtime marker missing: {marker}")
    value = float(match.group(1))
    if not math.isfinite(value):
        raise BuiltinBenchmarkError(f"runtime marker is non-finite: {marker}")
    return value


def _parse_int_marker(text: str, marker: str) -> int:
    match = re.search(re.escape(marker) + r"([0-9]+)", text)
    if match is None:
        raise BuiltinBenchmarkError(f"required runtime marker missing: {marker}")
    return int(match.group(1))


def _extract_comsol_identity(text: str) -> dict[str, str]:
    identity: dict[str, str] = {}
    banner = re.search(
        r"COMSOL Multiphysics\s+([0-9.]+)\s+\(Build:\s*([^)]+)\)",
        text,
        flags=re.IGNORECASE,
    )
    if banner is not None:
        identity["version"] = banner.group(1)
        identity["build"] = banner.group(2).strip()
    progress = re.search(
        r"\*+COMSOL\s+([0-9]+(?:\.[0-9]+)+)\s+progress output file",
        text,
        flags=re.IGNORECASE,
    )
    if progress is not None:
        identity["progress_version"] = progress.group(1)
    return identity


def _error_signature(stdout: str, stderr: str) -> str | None:
    text = f"{stdout}\n{stderr}".lower()
    for name, marker in ERROR_SIGNATURES:
        if marker in text:
            return name
    return None


def _read_step_text(result: dict[str, Any]) -> str:
    return "\n".join(
        Path(result[key]).read_text(encoding="utf-8", errors="replace")
        for key in ("stdout", "stderr")
    )


def _file_artifact(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    payload: dict[str, Any] = {
        "path": str(resolved),
        "exists": resolved.is_file(),
    }
    if not resolved.is_file():
        return payload
    stat = resolved.stat()
    payload.update(
        {
            "size_bytes": stat.st_size,
            "modified_at_utc": datetime.fromtimestamp(
                stat.st_mtime, timezone.utc
            ).isoformat(),
            "sha256": _sha256_file(resolved),
        }
    )
    return payload


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _java_path(path: Path) -> str:
    return str(path.resolve()).replace("\\", "/").replace('"', '\\"')


def _java_number(value: float) -> str:
    if not math.isfinite(value):
        raise BuiltinBenchmarkError("condition values must be finite")
    return format(float(value), ".17g")


def _completed_text(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Run the COMSOL built-in Boltzmann 1D reference in clean, "
            "independent batch processes."
        )
    )
    parser.add_argument(
        "--input-mph",
        type=Path,
        default=Path("Model/positive_column_1d_boltzmann.mph"),
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path(
            "outputs/comsol_swarm_benchmark_2026/builtin_200elem_runs"
        ),
    )
    parser.add_argument("--comsol", type=str, default=None)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--start-index", type=int, default=1)
    parser.add_argument("--voltage-V", type=float, default=200.0)
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
        input_mph=args.input_mph,
        output_root=args.output_root,
        comsol_executable=args.comsol,
        repetitions=args.repetitions,
        start_index=args.start_index,
        voltage_V=args.voltage_V,
        pressure_Pa=args.pressure_Pa,
        gas_temperature_K=args.gas_temperature_K,
        mesh_elements=args.mesh_elements,
    )
    for result in results:
        print(
            f"run_{result.index:02d}: solve={result.solve_time_s:.6g} s, "
            f"end_to_end={result.end_to_end_time_s:.6g} s, "
            f"rows={result.profile_rows}, "
            f"provenance={result.paths.provenance_json}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
