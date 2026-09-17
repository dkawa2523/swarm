"""Execute and validate individual COMSOL batch process steps."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import subprocess
from time import perf_counter
from typing import Any

from ..._io import write_json
from .artifacts import file_artifact
from .types import ComsolAdapterError


@dataclass(frozen=True, slots=True)
class CommandStep:
    name: str
    args: list[str]
    stdout_path: Path
    stderr_path: Path
    result_path: Path


def command_step(log_dir: Path, name: str, args: list[str]) -> CommandStep:
    return CommandStep(
        name=name,
        args=list(args),
        stdout_path=log_dir / f"{name}_stdout.txt",
        stderr_path=log_dir / f"{name}_stderr.txt",
        result_path=log_dir / f"{name}_result.json",
    )


def run_logged_command(
    step: CommandStep,
    *,
    cwd: Path,
    log_dir: Path,
) -> dict[str, Any]:
    started = perf_counter()
    return_code = -1
    process_failed = False
    os_error: OSError | None = None
    try:
        with (
            step.stdout_path.open("wb") as stdout_stream,
            step.stderr_path.open("wb") as stderr_stream,
        ):
            completed = subprocess.run(
                step.args,
                cwd=str(cwd),
                shell=False,
                check=True,
                stdout=stdout_stream,
                stderr=stderr_stream,
            )
        return_code = int(completed.returncode)
        _write_returned_process_text(step.stdout_path, completed.stdout)
        _write_returned_process_text(step.stderr_path, completed.stderr)
    except subprocess.CalledProcessError as exc:
        return_code = int(exc.returncode)
        process_failed = True
        _write_returned_process_text(step.stdout_path, exc.stdout)
        _write_returned_process_text(step.stderr_path, exc.stderr)
    except OSError as exc:
        os_error = exc

    stdout = _read_process_log(step.stdout_path)
    stderr = _read_process_log(step.stderr_path)
    if os_error is not None:
        stderr = "\n".join(value for value in (stderr, str(os_error)) if value)
        step.stderr_path.write_text(stderr, encoding="utf-8")
        detected_error = "os_error"
    elif process_failed:
        detected_error = "nonzero_return_code"
    else:
        detected_error = _comsol_error_signature(stdout, stderr)
    status = "failed" if detected_error is not None else "ok"
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
    write_json(step.result_path, result)
    if status == "failed":
        raise ComsolAdapterError(
            "COMSOL batch step failed: "
            f"{step.name}; log_dir: {log_dir}; "
            f"stdout: {step.stdout_path}; stderr: {step.stderr_path}",
            step_result=result,
        )
    return result


def ensure_clean_class_target(java_path: Path, log_dir: Path) -> None:
    class_file = java_path.with_suffix(".class")
    if class_file.exists():
        raise ComsolAdapterError(
            "refusing to compile with a pre-existing class file in the per-run "
            f"directory: {class_file}; log_dir: {log_dir}"
        )


def ensure_class_file_exists(
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


def ensure_fresh_output(
    output_mph: Path,
    log_dir: Path,
    step_result: dict[str, Any],
    *,
    preexisting_output_mph: dict[str, Any] | None,
) -> None:
    output = file_artifact(output_mph)
    if output is None or not output.get("exists"):
        step_result["status"] = "failed"
        step_result["missing_artifact"] = str(output_mph)
        raise ComsolAdapterError(
            "COMSOL operation did not produce the expected output .mph: "
            f"{output_mph}; log_dir: {log_dir}; "
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
        step_result["stale_artifact"] = str(output_mph)
        raise ComsolAdapterError(
            "COMSOL operation left the pre-existing output .mph unchanged; "
            f"refusing a stale result: {output_mph}; log_dir: {log_dir}; "
            f"stdout: {step_result['stdout']}; stderr: {step_result['stderr']}",
            step_result=step_result,
        )


def _read_process_log(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_bytes().decode("utf-8", errors="replace")


def _write_returned_process_text(path: Path, value: object) -> None:
    """Support subprocess doubles while real commands write directly to disk."""

    if value is not None:
        path.write_text(_completed_text(value), encoding="utf-8")


def _completed_text(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


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
