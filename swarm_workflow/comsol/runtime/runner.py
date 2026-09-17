"""Compile and execute model-generated COMSOL Java sources."""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from time import perf_counter
from typing import Any

from ..._io import write_json
from .artifacts import (
    file_artifact,
    new_operation_log_dir,
    utc_now,
    write_execution_provenance,
)
from .commands import build_java_batch_command, build_java_compile_command
from .executable import resolve_comsol_executable
from .process import (
    CommandStep,
    command_step,
    ensure_class_file_exists,
    ensure_clean_class_target,
    ensure_fresh_output,
    run_logged_command,
)
from .types import (
    ComsolAdapterError,
    ComsolExecutionSummary,
    ComsolRunContext,
    ResolvedComsolExecutable,
)


def execute_comsol_java_source(
    context: ComsolRunContext,
    source: str,
    *,
    class_name: str,
    operation: str,
    comsol_executable: str | Path | None = None,
    study: str,
    require_fresh_output: bool = False,
) -> ComsolExecutionSummary:
    """Materialize, compile, and execute a Java source in one run directory."""

    resolved = resolve_comsol_executable(comsol_executable)
    log_dir = new_operation_log_dir(context.log_path, operation)
    java_path = log_dir / f"{class_name}.java"
    java_path.write_text(source, encoding="utf-8")
    context.output_mph.parent.mkdir(parents=True, exist_ok=True)
    steps, compile_targets = _compile_and_run_steps(
        log_dir,
        operation,
        resolved,
        java_path,
    )
    return _execute_steps(
        operation=operation,
        context=context,
        resolved=resolved,
        log_dir=log_dir,
        steps=steps,
        compile_targets=compile_targets,
        java_path=java_path,
        support_java_paths=(),
        source_support_java_paths=(),
        study=study,
        require_fresh_output=require_fresh_output,
    )


def execute_generated_comsol_java(
    context: ComsolRunContext,
    java_path: str | Path,
    *,
    operation: str,
    comsol_executable: str | Path | None = None,
    study: str,
    require_fresh_output: bool = False,
    support_java_paths: Sequence[str | Path] = (),
) -> ComsolExecutionSummary:
    """Stage, compile, and execute an existing generated Java program.

    Supporting public classes are compiled beside the entry point and loaded
    by the same COMSOL batch process.  This is the supported way to compose
    independently generated stages without repeating licensed JVM startup.
    """

    resolved = resolve_comsol_executable(comsol_executable)
    resolved_java_path = Path(java_path).resolve()
    if not resolved_java_path.exists():
        raise ComsolAdapterError(f"generated Java source does not exist: {java_path}")
    resolved_support_paths = tuple(Path(path).resolve() for path in support_java_paths)
    missing_support = next(
        (path for path in resolved_support_paths if not path.is_file()),
        None,
    )
    if missing_support is not None:
        raise ComsolAdapterError(
            f"generated Java support source does not exist: {missing_support}"
        )
    source_names = (resolved_java_path.name,) + tuple(
        path.name for path in resolved_support_paths
    )
    if len(source_names) != len(set(source_names)):
        raise ComsolAdapterError(
            "generated Java entry and support source filenames must be unique"
        )
    log_dir = new_operation_log_dir(context.log_path, operation)
    staged_java_path = log_dir / resolved_java_path.name
    staged_java_path.write_bytes(resolved_java_path.read_bytes())
    staged_support_paths: list[Path] = []
    for support_path in resolved_support_paths:
        staged = log_dir / support_path.name
        staged.write_bytes(support_path.read_bytes())
        staged_support_paths.append(staged)
    steps, compile_targets = _compile_and_run_steps(
        log_dir,
        operation,
        resolved,
        staged_java_path,
        support_java_paths=tuple(staged_support_paths),
    )
    return _execute_steps(
        operation=operation,
        context=context,
        resolved=resolved,
        log_dir=log_dir,
        steps=steps,
        compile_targets=compile_targets,
        java_path=staged_java_path,
        support_java_paths=tuple(staged_support_paths),
        source_support_java_paths=resolved_support_paths,
        source_java_path=resolved_java_path,
        study=study,
        require_fresh_output=require_fresh_output,
    )


def _compile_and_run_steps(
    log_dir: Path,
    operation: str,
    resolved: ResolvedComsolExecutable,
    java_path: Path,
    *,
    support_java_paths: tuple[Path, ...] = (),
) -> tuple[list[CommandStep], dict[str, Path]]:
    class_file = java_path.with_suffix(".class")
    support_steps = [
        command_step(
            log_dir,
            f"compile_support_{index:03d}",
            build_java_compile_command(resolved.path, path),
        )
        for index, path in enumerate(support_java_paths)
    ]
    steps = [
        command_step(
            log_dir,
            "compile",
            build_java_compile_command(resolved.path, java_path),
        ),
        *support_steps,
        command_step(
            log_dir,
            operation,
            build_java_batch_command(
                resolved.path,
                class_file,
                license_path=resolved.license_path,
            ),
        ),
    ]
    compile_targets = {
        step.name: path for step, path in zip(support_steps, support_java_paths)
    }
    compile_targets["compile"] = java_path
    return steps, compile_targets


def _execute_steps(
    *,
    operation: str,
    context: ComsolRunContext,
    resolved: ResolvedComsolExecutable,
    log_dir: Path,
    steps: list[CommandStep],
    compile_targets: dict[str, Path],
    java_path: Path | None,
    support_java_paths: tuple[Path, ...],
    source_support_java_paths: tuple[Path, ...],
    study: str,
    source_java_path: Path | None = None,
    require_fresh_output: bool = False,
) -> ComsolExecutionSummary:
    started = perf_counter()
    started_at_utc = utc_now()
    command_json = log_dir / "command.json"
    result_json = log_dir / "result.json"
    provenance_json = log_dir / "provenance.json"
    preexisting_output_mph = file_artifact(context.output_mph)
    if java_path is not None:
        ensure_clean_class_target(java_path, log_dir)
    for support_path in support_java_paths:
        ensure_clean_class_target(support_path, log_dir)
    command_payload = {
        "operation": operation,
        "mapping_path": str(context.mapping_path),
        "cwd": str(context.root),
        "executable": resolved.path,
        "executable_source": resolved.source,
        "license_file": file_artifact(
            Path(resolved.license_path) if resolved.license_path is not None else None
        ),
        "license_source": resolved.license_source,
        "generated_java_path": str(java_path) if java_path is not None else None,
        "generated_support_java_paths": [str(path) for path in support_java_paths],
        "source_java_path": (
            str(source_java_path)
            if source_java_path is not None
            else (str(java_path) if java_path is not None else None)
        ),
        "source_support_java_paths": [str(path) for path in source_support_java_paths],
        "input_mph": str(context.input_mph),
        "output_mph": str(context.output_mph),
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
    write_json(command_json, command_payload)
    write_execution_provenance(
        provenance_json,
        status="running",
        operation=operation,
        context=context,
        resolved=resolved,
        java_path=java_path,
        source_java_path=source_java_path,
        support_java_paths=support_java_paths,
        source_support_java_paths=source_support_java_paths,
        started_at_utc=started_at_utc,
        results=[],
        total_time_s=None,
        preexisting_output_mph=preexisting_output_mph,
        study=study,
    )

    results: list[dict[str, Any]] = []
    try:
        for step in steps:
            result = run_logged_command(step, cwd=context.root, log_dir=log_dir)
            results.append(result)
            compile_target = compile_targets.get(step.name)
            if compile_target is not None:
                try:
                    ensure_class_file_exists(compile_target, log_dir, result)
                    result["compiled_class"] = file_artifact(
                        compile_target.with_suffix(".class")
                    )
                finally:
                    write_json(step.result_path, result)
        if require_fresh_output:
            try:
                ensure_fresh_output(
                    context.output_mph,
                    log_dir,
                    results[-1],
                    preexisting_output_mph=preexisting_output_mph,
                )
            finally:
                write_json(steps[-1].result_path, results[-1])
    except ComsolAdapterError as exc:
        if exc.step_result is not None and (
            not results or results[-1] is not exc.step_result
        ):
            results.append(exc.step_result)
        total_time_s = perf_counter() - started
        write_execution_provenance(
            provenance_json,
            status="failed",
            operation=operation,
            context=context,
            resolved=resolved,
            java_path=java_path,
            source_java_path=source_java_path,
            support_java_paths=support_java_paths,
            source_support_java_paths=source_support_java_paths,
            started_at_utc=started_at_utc,
            results=results,
            total_time_s=total_time_s,
            preexisting_output_mph=preexisting_output_mph,
            study=study,
        )
        write_json(
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
    provenance = write_execution_provenance(
        provenance_json,
        status="ok",
        operation=operation,
        context=context,
        resolved=resolved,
        java_path=java_path,
        source_java_path=source_java_path,
        support_java_paths=support_java_paths,
        source_support_java_paths=source_support_java_paths,
        started_at_utc=started_at_utc,
        results=results,
        total_time_s=total_time_s,
        preexisting_output_mph=preexisting_output_mph,
        study=study,
    )
    write_json(
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
        log_dir=log_dir,
        command_json=command_json,
        result_json=result_json,
        stdout_paths=tuple(step.stdout_path for step in steps),
        stderr_paths=tuple(step.stderr_path for step in steps),
        return_codes=tuple(int(result["return_code"]) for result in results),
        java_path=java_path,
        output_mph=context.output_mph,
        total_time_s=total_time_s,
        provenance_json=provenance_json,
        comsol_version=provenance["comsol"]["version"],
        comsol_build=provenance["comsol"]["build"],
    )
