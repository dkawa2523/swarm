"""Build COMSOL Java compilation and batch commands."""

from __future__ import annotations

import os
from pathlib import Path


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


def build_java_batch_command(
    comsol_executable: str | Path,
    class_file: str | Path,
    *,
    license_path: str | Path | None = None,
) -> list[str]:
    license_args = ["-c", str(license_path)] if license_path is not None else []
    return [
        *_batch_prefix(str(comsol_executable)),
        *license_args,
        "-inputfile",
        str(class_file),
    ]


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
