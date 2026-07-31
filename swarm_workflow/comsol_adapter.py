"""Dry-run planning and local COMSOL batch execution helpers."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import os
from pathlib import Path
import re
import shutil
import subprocess
from time import perf_counter
from typing import Any

from ._io import write_json as _write_json
from .comsol_java import (
    DEFAULT_CLASS_NAME,
    VOLTAGE_LOG_PREFIX,
    generate_apply_java_source,
)
from .comsol_mapping import (
    ComsolModelMapping,
    closure_quantity_activity,
    dry_run_mapping_summary,
    load_comsol_mapping,
    validate_comsol_mapping_files,
)


COMSOL_EXECUTABLE_PLACEHOLDER = "<COMSOL_EXECUTABLE>"


class ComsolAdapterError(RuntimeError):
    """Raised when local COMSOL batch execution cannot be planned or run."""

    def __init__(
        self,
        message: str,
        *,
        step_result: dict[str, Any] | None = None,
    ) -> None:
        super().__init__(message)
        self.step_result = step_result


@dataclass(frozen=True, slots=True)
class ApplyComsolPlan:
    mapping: ComsolModelMapping
    java_path: Path
    summary: dict[str, Any]


@dataclass(frozen=True, slots=True)
class JavaWriteSummary:
    mapping: ComsolModelMapping
    java_path: Path
    bytes_written: int


@dataclass(frozen=True, slots=True)
class ResolvedComsolExecutable:
    path: str
    source: str


@dataclass(frozen=True, slots=True)
class ComsolExecutionSummary:
    operation: str
    mapping: ComsolModelMapping
    log_dir: Path
    command_json: Path
    result_json: Path
    stdout_paths: tuple[Path, ...]
    stderr_paths: tuple[Path, ...]
    return_codes: tuple[int, ...]
    java_path: Path | None
    output_mph: Path
    total_time_s: float
    provenance_json: Path | None = None
    comsol_version: str | None = None
    comsol_build: str | None = None
    executed_voltage_sequence_V: tuple[float, ...] | None = None


@dataclass(frozen=True, slots=True)
class _CommandStep:
    name: str
    args: list[str]
    stdout_path: Path
    stderr_path: Path
    result_path: Path


def plan_apply_comsol(
    mapping_path: str | Path,
    java_path: str | Path | None = None,
) -> ApplyComsolPlan:
    mapping = load_comsol_mapping(mapping_path)
    validate_comsol_mapping_files(mapping, require_unit_metadata=True)
    resolved_java_path = (
        Path(java_path).resolve()
        if java_path is not None
        else _default_java_path(mapping)
    )
    summary = dry_run_mapping_summary(mapping)
    summary["java"] = {"path": str(resolved_java_path)}
    return ApplyComsolPlan(
        mapping=mapping,
        java_path=resolved_java_path,
        summary=summary,
    )


def write_apply_comsol_java(
    mapping_path: str | Path,
    java_path: str | Path,
) -> JavaWriteSummary:
    plan = plan_apply_comsol(mapping_path, java_path=java_path)
    source = generate_apply_java_source(plan.mapping)
    plan.java_path.parent.mkdir(parents=True, exist_ok=True)
    plan.java_path.write_text(source, encoding="utf-8")
    return JavaWriteSummary(
        mapping=plan.mapping,
        java_path=plan.java_path,
        bytes_written=len(source.encode("utf-8")),
    )


def execute_apply_comsol(
    mapping_path: str | Path,
    comsol_executable: str | Path | None = None,
    *,
    pressure_Pa: float | None = None,
    gas_temperature_K: float | None = None,
    mesh_elements: int | None = None,
) -> ComsolExecutionSummary:
    """Generate Java, compile it with COMSOL, then run it through batch mode."""

    mapping = load_comsol_mapping(mapping_path)
    validate_comsol_mapping_files(mapping, require_unit_metadata=True)
    _ensure_output_is_not_input(mapping)
    resolved = resolve_comsol_executable(comsol_executable)
    log_dir = _new_operation_log_dir(mapping, "apply")
    java_path = log_dir / f"{DEFAULT_CLASS_NAME}.java"
    source = generate_apply_java_source(
        mapping,
        pressure_Pa=pressure_Pa,
        gas_temperature_K=gas_temperature_K,
        mesh_elements=mesh_elements,
    )
    java_path.write_text(source, encoding="utf-8")
    mapping.model.output_mph.parent.mkdir(parents=True, exist_ok=True)
    class_file = java_path.with_suffix(".class")
    # COMSOL command-line tools support compiling Java model code before batch
    # execution; the command variant depends on whether the user points us at
    # the launcher, comsolbatch, or comsolcompile.
    steps = [
        _step(
            log_dir,
            "compile",
            build_java_compile_command(resolved.path, java_path),
        ),
        _step(
            log_dir,
            "apply",
            build_apply_java_batch_command(resolved.path, class_file),
        ),
    ]
    return _execute_steps(
        operation="apply",
        mapping=mapping,
        resolved=resolved,
        log_dir=log_dir,
        steps=steps,
        java_path=java_path,
        study=mapping.model.study,
    )


def execute_generated_comsol_java(
    mapping: ComsolModelMapping,
    java_path: str | Path,
    *,
    operation: str,
    comsol_executable: str | Path | None = None,
) -> ComsolExecutionSummary:
    resolved = resolve_comsol_executable(comsol_executable)
    resolved_java_path = Path(java_path).resolve()
    if not resolved_java_path.exists():
        raise ComsolAdapterError(f"generated Java source does not exist: {java_path}")
    log_dir = _new_operation_log_dir(mapping, operation)
    staged_java_path = log_dir / resolved_java_path.name
    staged_java_path.write_bytes(resolved_java_path.read_bytes())
    class_file = staged_java_path.with_suffix(".class")
    steps = [
        _step(
            log_dir,
            "compile",
            build_java_compile_command(resolved.path, staged_java_path),
        ),
        _step(
            log_dir,
            operation,
            build_apply_java_batch_command(resolved.path, class_file),
        ),
    ]
    return _execute_steps(
        operation=operation,
        mapping=mapping,
        resolved=resolved,
        log_dir=log_dir,
        steps=steps,
        java_path=staged_java_path,
        source_java_path=resolved_java_path,
        study=mapping.model.study,
    )


def resolve_comsol_executable(
    comsol_executable: str | Path | None = None,
) -> ResolvedComsolExecutable:
    if comsol_executable is not None:
        return _resolve_candidate(str(comsol_executable), source="cli")
    for env_name in ("COMSOL_BATCH", "COMSOL_EXECUTABLE"):
        value = os.environ.get(env_name)
        if value:
            return _resolve_candidate(value, source=f"env:{env_name}")
    for name in ("comsol", "comsolbatch"):
        found = shutil.which(name)
        if found:
            return ResolvedComsolExecutable(
                path=str(Path(found).resolve()),
                source=f"path:{name}",
            )
    raise ComsolAdapterError(
        "COMSOL executable not found; pass --comsol, set COMSOL_BATCH or "
        "COMSOL_EXECUTABLE, or put comsol/comsolbatch on PATH"
    )


def find_comsol_executable(
    comsol_executable: str | Path | None = None,
) -> ResolvedComsolExecutable | None:
    try:
        return resolve_comsol_executable(comsol_executable)
    except ComsolAdapterError:
        return None


def build_java_compile_command(
    comsol_executable: str | Path,
    java_path: str | Path,
) -> list[str]:
    executable = str(comsol_executable)
    stem = Path(executable).stem.lower()
    if stem == "comsolcompile":
        return [executable, str(java_path)]
    if stem == "comsolbatch" or _is_windows_comsol_launcher(executable):
        return [_sibling_command(executable, "comsolcompile"), str(java_path)]
    return [executable, "compile", str(java_path)]


def build_apply_java_batch_command(
    comsol_executable: str | Path,
    class_file: str | Path,
) -> list[str]:
    return [*_batch_prefix(str(comsol_executable)), "-inputfile", str(class_file)]


def format_apply_plan(plan: ApplyComsolPlan) -> str:
    lines = [
        "COMSOL apply dry-run",
        f"mapping: {plan.mapping.path}",
        f"input_mph: {plan.mapping.model.input_mph}",
        f"output_mph: {plan.mapping.model.output_mph}",
        f"bundle: {plan.mapping.bundle.path}",
        f"logs: {plan.mapping.logs.path}",
        f"java_path: {plan.java_path}",
        "interpolation_functions:",
    ]
    for function in plan.mapping.functions:
        lines.append(
            "  - "
            f"{function.name}: tag={function.tag}, file={function.path}, "
            f"column={function.column or ''}, nargs={function.nargs}, "
            f"argunit={function.argunit}, fununit={function.fununit}, "
            f"interp={function.interp}, extrap={function.extrap}"
        )
    closure = plan.mapping.closure
    table = closure.mean_energy
    formulation = closure.mean_energy_formulation
    activity = closure_quantity_activity(formulation)
    lines.append("closure:")
    lines.append(
        "  - written_transport: transport_vs_mean_energy.csv; mean_energy_table: "
        f"{table.table} ({table.x_column}, {table.y_column})"
    )
    lines.append(
        "  - mean_energy_formulation: "
        f"{formulation.mode}; {formulation.property_group}/"
        f"{formulation.property}={formulation.comsol_value}"
    )
    lines.append(
        "  - active: "
        + ", ".join(
            quantity
            for quantity, item in activity.items()
            if item["status"] == "active"
        )
    )
    lines.append(
        "  - inactive_written: "
        + ", ".join(
            quantity
            for quantity, item in activity.items()
            if item["status"] == "inactive"
        )
    )
    if plan.mapping.reaction_lookups:
        lines.append("reaction_lookups:")
        for lookup in plan.mapping.reaction_lookups:
            lines.append(
                "  - "
                f"{lookup.name}: {lookup.physics}.{lookup.feature}, "
                f"form={lookup.form}, source_model={lookup.source_model}, "
                f"file={lookup.path}, "
                f"x={lookup.x_column}, y={lookup.y_column}, "
                f"process_type={lookup.process_type}"
            )
    lines.append("COMSOL command: not executed")
    return "\n".join(lines)


def format_execution_summary(summary: ComsolExecutionSummary) -> str:
    return "\n".join(
        [
            f"COMSOL {summary.operation} completed",
            f"log_dir: {summary.log_dir}",
            f"command_json: {summary.command_json}",
            f"result_json: {summary.result_json}",
            f"provenance_json: {summary.provenance_json or ''}",
            f"generated_java: {summary.java_path or ''}",
            f"output_mph: {summary.output_mph}",
            f"return_codes: {list(summary.return_codes)}",
        ]
    )


def _resolve_candidate(value: str, *, source: str) -> ResolvedComsolExecutable:
    path = Path(value)
    if path.exists():
        return ResolvedComsolExecutable(path=str(path.resolve()), source=source)
    found = shutil.which(value)
    if found:
        return ResolvedComsolExecutable(path=str(Path(found).resolve()), source=source)
    raise ComsolAdapterError(f"COMSOL executable does not exist: {value}")


def _dry_run_executable(comsol_executable: str | Path | None) -> str:
    if comsol_executable is not None:
        return str(comsol_executable)
    for env_name in ("COMSOL_BATCH", "COMSOL_EXECUTABLE"):
        value = os.environ.get(env_name)
        if value:
            return value
    for name in ("comsol", "comsolbatch"):
        found = shutil.which(name)
        if found:
            return str(Path(found).resolve())
    return COMSOL_EXECUTABLE_PLACEHOLDER


def _batch_prefix(executable: str) -> list[str]:
    stem = Path(executable).stem.lower()
    if stem == "comsolbatch":
        return [executable]
    if stem == "comsolcompile" or _is_windows_comsol_launcher(executable):
        return [_sibling_command(executable, "comsolbatch")]
    return [executable, "batch"]


def _is_windows_comsol_launcher(executable: str) -> bool:
    path = Path(executable)
    return path.stem.lower() == "comsol" and (
        os.name == "nt" or path.suffix.lower() == ".exe"
    )


def _sibling_command(executable: str, command_name: str) -> str:
    path = Path(executable)
    return str(path.with_name(command_name + path.suffix))


def _default_java_path(mapping: ComsolModelMapping) -> Path:
    return mapping.model.output_mph.with_suffix(".java")


def _new_operation_log_dir(mapping: ComsolModelMapping, operation: str) -> Path:
    base = mapping.logs.path
    base.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    candidate = base / f"{operation}_{stamp}"
    suffix = 1
    while candidate.exists():
        suffix += 1
        candidate = base / f"{operation}_{stamp}_{suffix}"
    candidate.mkdir(parents=True)
    return candidate


def _step(log_dir: Path, name: str, args: list[str]) -> _CommandStep:
    return _CommandStep(
        name=name,
        args=list(args),
        stdout_path=log_dir / f"{name}_stdout.txt",
        stderr_path=log_dir / f"{name}_stderr.txt",
        result_path=log_dir / f"{name}_result.json",
    )


def _execute_steps(
    *,
    operation: str,
    mapping: ComsolModelMapping,
    resolved: ResolvedComsolExecutable,
    log_dir: Path,
    steps: list[_CommandStep],
    java_path: Path | None,
    study: str,
    source_java_path: Path | None = None,
) -> ComsolExecutionSummary:
    started = perf_counter()
    started_at_utc = _utc_now()
    cwd = mapping.root
    command_json = log_dir / "command.json"
    result_json = log_dir / "result.json"
    provenance_json = log_dir / "provenance.json"
    preexisting_output_mph = _file_artifact(mapping.model.output_mph)
    if java_path is not None:
        _ensure_clean_class_target(java_path, log_dir)
    command_payload = {
        "operation": operation,
        "mapping_path": str(mapping.path),
        "cwd": str(cwd),
        "executable": resolved.path,
        "executable_source": resolved.source,
        "generated_java_path": str(java_path) if java_path is not None else None,
        "source_java_path": (
            str(source_java_path)
            if source_java_path is not None
            else (str(java_path) if java_path is not None else None)
        ),
        "input_mph": str(mapping.model.input_mph),
        "output_mph": str(mapping.model.output_mph),
        "study": study,
        "steps": [
            {
                "name": step.name,
                "args": step.args,
                "stdout": str(step.stdout_path),
                "stderr": str(step.stderr_path),
            }
            for step in steps
        ],
    }
    _write_json(command_json, command_payload)
    _write_execution_provenance(
        provenance_json,
        status="running",
        operation=operation,
        mapping=mapping,
        resolved=resolved,
        java_path=java_path,
        source_java_path=source_java_path,
        started_at_utc=started_at_utc,
        results=[],
        total_time_s=None,
        preexisting_output_mph=preexisting_output_mph,
    )

    results: list[dict[str, Any]] = []
    try:
        for step in steps:
            result = _run_logged_command(step, cwd=cwd, log_dir=log_dir)
            results.append(result)
            if step.name == "compile" and java_path is not None:
                try:
                    _ensure_class_file_exists(java_path, log_dir, result)
                    result["compiled_class"] = _file_artifact(
                        java_path.with_suffix(".class")
                    )
                finally:
                    _write_json(step.result_path, result)
        if operation == "apply":
            try:
                _ensure_apply_output_is_fresh(
                    mapping,
                    log_dir,
                    results[-1],
                    preexisting_output_mph=preexisting_output_mph,
                )
            finally:
                _write_json(steps[-1].result_path, results[-1])
    except ComsolAdapterError as exc:
        if exc.step_result is not None and (
            not results or results[-1] is not exc.step_result
        ):
            results.append(exc.step_result)
        total_time_s = perf_counter() - started
        _write_execution_provenance(
            provenance_json,
            status="failed",
            operation=operation,
            mapping=mapping,
            resolved=resolved,
            java_path=java_path,
            source_java_path=source_java_path,
            started_at_utc=started_at_utc,
            results=results,
            total_time_s=total_time_s,
            preexisting_output_mph=preexisting_output_mph,
        )
        _write_json(
            result_json,
            {
                "status": "failed",
                "operation": operation,
                "log_dir": str(log_dir),
                "provenance": str(provenance_json),
                "steps": results,
                "total_time_s": total_time_s,
            },
        )
        raise

    total_time_s = perf_counter() - started
    provenance = _write_execution_provenance(
        provenance_json,
        status="ok",
        operation=operation,
        mapping=mapping,
        resolved=resolved,
        java_path=java_path,
        source_java_path=source_java_path,
        started_at_utc=started_at_utc,
        results=results,
        total_time_s=total_time_s,
        preexisting_output_mph=preexisting_output_mph,
    )
    _write_json(
        result_json,
        {
            "status": "ok",
            "operation": operation,
            "log_dir": str(log_dir),
            "provenance": str(provenance_json),
            "steps": results,
            "total_time_s": total_time_s,
        },
    )
    return ComsolExecutionSummary(
        operation=operation,
        mapping=mapping,
        log_dir=log_dir,
        command_json=command_json,
        result_json=result_json,
        stdout_paths=tuple(step.stdout_path for step in steps),
        stderr_paths=tuple(step.stderr_path for step in steps),
        return_codes=tuple(int(result["return_code"]) for result in results),
        java_path=java_path,
        output_mph=mapping.model.output_mph,
        total_time_s=total_time_s,
        provenance_json=provenance_json,
        comsol_version=provenance["comsol"]["version"],
        comsol_build=provenance["comsol"]["build"],
        executed_voltage_sequence_V=(
            tuple(provenance["executed_voltage_sequence_V"])
            if provenance["executed_voltage_sequence_V"] is not None
            else None
        ),
    )


def _run_logged_command(
    step: _CommandStep,
    *,
    cwd: Path,
    log_dir: Path,
) -> dict[str, Any]:
    started = perf_counter()
    try:
        completed = subprocess.run(
            step.args,
            cwd=str(cwd),
            shell=False,
            check=True,
            capture_output=True,
            text=True,
        )
        return_code = int(completed.returncode)
        stdout = _completed_text(completed.stdout)
        stderr = _completed_text(completed.stderr)
        detected_error = _comsol_error_signature(stdout, stderr)
        status = "failed" if detected_error is not None else "ok"
    except subprocess.CalledProcessError as exc:
        return_code = int(exc.returncode)
        stdout = _completed_text(exc.stdout)
        stderr = _completed_text(exc.stderr)
        status = "failed"
        detected_error = "nonzero_return_code"
    except OSError as exc:
        return_code = -1
        stdout = ""
        stderr = str(exc)
        status = "failed"
        detected_error = "os_error"

    step.stdout_path.write_text(stdout, encoding="utf-8")
    step.stderr_path.write_text(stderr, encoding="utf-8")
    result = {
        "name": step.name,
        "status": status,
        "args": step.args,
        "cwd": str(cwd),
        "return_code": return_code,
        "stdout": str(step.stdout_path),
        "stderr": str(step.stderr_path),
        "duration_s": perf_counter() - started,
    }
    if detected_error is not None:
        result["detected_comsol_error"] = detected_error
    _write_json(step.result_path, result)
    if status == "failed":
        raise ComsolAdapterError(
            "COMSOL batch step failed: "
            f"{step.name}; log_dir: {log_dir}; "
            f"stdout: {step.stdout_path}; stderr: {step.stderr_path}",
            step_result=result,
        )
    return result


def _completed_text(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def _looks_like_comsol_error(stdout: str, stderr: str) -> bool:
    return _comsol_error_signature(stdout, stderr) is not None


def _comsol_error_signature(stdout: str, stderr: str) -> str | None:
    text = f"{stdout}\n{stderr}".lower()
    signatures = (
        ("compiler_failed", "failed to compile java file"),
        ("compilation_failed", "compilation failed"),
        ("compilation_error", "compilation error on line"),
        ("java_class_error", "error running java class"),
        ("comsol_error_block", "/*****error"),
        ("undefined_variable", "undefined variable"),
        ("unknown_function", "unknown function or operator"),
    )
    for name, marker in signatures:
        if marker in text:
            return name
    return None


def _ensure_clean_class_target(java_path: Path, log_dir: Path) -> None:
    class_file = java_path.with_suffix(".class")
    if not class_file.exists():
        return
    raise ComsolAdapterError(
        "refusing to compile with a pre-existing class file in the per-run "
        f"directory: {class_file}; log_dir: {log_dir}"
    )


def _ensure_class_file_exists(
    java_path: Path,
    log_dir: Path,
    step_result: dict[str, Any],
) -> None:
    class_file = java_path.with_suffix(".class")
    if class_file.is_file() and class_file.stat().st_size > 0:
        return
    step_result["status"] = "failed"
    step_result["missing_artifact"] = str(class_file)
    raise ComsolAdapterError(
        "COMSOL Java compile did not produce a class file: "
        f"{class_file}; log_dir: {log_dir}; "
        f"stdout: {step_result['stdout']}; stderr: {step_result['stderr']}",
        step_result=step_result,
    )


def _ensure_apply_output_is_fresh(
    mapping: ComsolModelMapping,
    log_dir: Path,
    step_result: dict[str, Any],
    *,
    preexisting_output_mph: dict[str, Any] | None,
) -> None:
    output = _file_artifact(mapping.model.output_mph)
    if output is None or not output.get("exists"):
        step_result["status"] = "failed"
        step_result["missing_artifact"] = str(mapping.model.output_mph)
        raise ComsolAdapterError(
            "COMSOL apply did not produce the mapped output .mph: "
            f"{mapping.model.output_mph}; log_dir: {log_dir}; "
            f"stdout: {step_result['stdout']}; stderr: {step_result['stderr']}",
            step_result=step_result,
        )
    if (
        preexisting_output_mph is not None
        and preexisting_output_mph.get("exists")
        and output.get("sha256") == preexisting_output_mph.get("sha256")
        and output.get("modified_at_utc")
        == preexisting_output_mph.get("modified_at_utc")
    ):
        step_result["status"] = "failed"
        step_result["stale_artifact"] = str(mapping.model.output_mph)
        raise ComsolAdapterError(
            "COMSOL apply left the pre-existing output .mph unchanged; "
            f"refusing a stale result: {mapping.model.output_mph}; "
            f"log_dir: {log_dir}; stdout: {step_result['stdout']}; "
            f"stderr: {step_result['stderr']}",
            step_result=step_result,
        )


def _ensure_output_is_not_input(mapping: ComsolModelMapping) -> None:
    if mapping.model.input_mph == mapping.model.output_mph:
        raise ComsolAdapterError(
            "mapping model.output_mph must differ from model.input_mph; "
            f"refusing to overwrite input .mph: {mapping.model.input_mph}"
        )


def _write_execution_provenance(
    path: Path,
    *,
    status: str,
    operation: str,
    mapping: ComsolModelMapping,
    resolved: ResolvedComsolExecutable,
    java_path: Path | None,
    source_java_path: Path | None,
    started_at_utc: str,
    results: list[dict[str, Any]],
    total_time_s: float | None,
    preexisting_output_mph: dict[str, Any] | None,
) -> dict[str, Any]:
    identity = _extract_comsol_identity(results)
    executed_voltages, execution_mode = _executed_voltage_sequence(
        operation,
        results,
    )
    source_origin = source_java_path if source_java_path is not None else java_path
    payload = {
        "format_version": 1,
        "status": status,
        "operation": operation,
        "recorded_at_utc": _utc_now(),
        "started_at_utc": started_at_utc,
        "total_time_s": total_time_s,
        "comsol": {
            "executable": resolved.path,
            "executable_source": resolved.source,
            "version": identity.get("version"),
            "build": identity.get("build"),
            "progress_version": identity.get("progress_version"),
        },
        "inputs": {
            "mapping": _file_artifact(mapping.path),
            "input_mph": _file_artifact(mapping.model.input_mph),
            "bundle_manifest": _file_artifact(mapping.bundle.path / "manifest.json"),
            "preexisting_output_mph": preexisting_output_mph,
        },
        "generated_artifacts": {
            "source_origin": _file_artifact(source_origin),
            "staged_java": _file_artifact(java_path),
            "compiled_class": _file_artifact(
                java_path.with_suffix(".class") if java_path is not None else None
            ),
            "output_mph": _file_artifact(mapping.model.output_mph),
        },
        "study": mapping.model.study,
        "configured_voltage_sequence_V": (
            list(mapping.run.voltages_V)
            if getattr(mapping.run, "voltages_V", None) is not None
            else None
        ),
        "executed_voltage_sequence_V": executed_voltages,
        "voltage_execution_mode": execution_mode,
        "steps": [
            {
                "name": result.get("name"),
                "status": result.get("status"),
                "return_code": result.get("return_code"),
                "detected_comsol_error": result.get("detected_comsol_error"),
            }
            for result in results
        ],
    }
    _write_json(path, payload)
    return payload


def _file_artifact(path: Path | None) -> dict[str, Any] | None:
    if path is None:
        return None
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
            "modified_at_utc": _timestamp_utc(stat.st_mtime),
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


def _extract_comsol_identity(results: list[dict[str, Any]]) -> dict[str, str]:
    texts: list[str] = []
    for result in results:
        for key in ("stdout", "stderr"):
            value = result.get(key)
            if not isinstance(value, str):
                continue
            path = Path(value)
            if path.is_file():
                texts.append(path.read_text(encoding="utf-8", errors="replace"))
    joined = "\n".join(texts)
    identity: dict[str, str] = {}
    banner = re.search(
        r"COMSOL Multiphysics\s+([0-9.]+)\s+\(Build:\s*([^)]+)\)",
        joined,
        flags=re.IGNORECASE,
    )
    if banner is not None:
        identity["version"] = banner.group(1)
        identity["build"] = banner.group(2).strip()
    progress = re.search(
        r"\*+COMSOL\s+([0-9]+(?:\.[0-9]+)+)\s+progress output file",
        joined,
        flags=re.IGNORECASE,
    )
    if progress is not None:
        identity["progress_version"] = progress.group(1)
    return identity


def _executed_voltage_sequence(
    operation: str,
    results: list[dict[str, Any]],
) -> tuple[list[float] | None, str | None]:
    if operation != "positive_column_run" or not results:
        return None, None
    runtime = next(
        (result for result in results if result.get("name") == operation),
        None,
    )
    if runtime is None:
        return None, None
    stdout_path = runtime.get("stdout")
    stdout = ""
    if isinstance(stdout_path, str) and Path(stdout_path).is_file():
        stdout = Path(stdout_path).read_text(
            encoding="utf-8",
            errors="replace",
        )
    executed = [
        float(value)
        for value in re.findall(
            re.escape(VOLTAGE_LOG_PREFIX)
            + r"\s*([+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)"
            + r"(?:[eE][+-]?[0-9]+)?)\s*\[V\]",
            stdout,
        )
    ]
    if not executed:
        return None, None
    if "Direct target-voltage solve failed; using continuation." in stdout:
        return executed, "direct_then_continuation"
    return executed, "direct"


def _timestamp_utc(value: float) -> str:
    return datetime.fromtimestamp(value, tz=timezone.utc).isoformat().replace(
        "+00:00",
        "Z",
    )


def _utc_now() -> str:
    return datetime.now(tz=timezone.utc).isoformat().replace("+00:00", "Z")
