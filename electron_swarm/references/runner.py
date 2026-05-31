"""External reference command execution helpers.

External reference binaries are optional benchmark dependencies.  This module
only runs an explicit user-provided command and verifies that it produced the
declared output file; it does not synthesize BOLSIG+ or MCIG data.
"""

from __future__ import annotations

from pathlib import Path
import subprocess


def _render_command(
    command: str,
    *,
    reference_id: str,
    output: Path,
    input_path: Path | None,
    config_path: Path | None,
) -> str:
    mapping = {
        "reference_id": reference_id,
        "output": str(output),
        "input": "" if input_path is None else str(input_path),
        "config": "" if config_path is None else str(config_path),
    }
    try:
        return command.format(**mapping)
    except KeyError as exc:
        raise ValueError(
            f"unknown external reference command placeholder: {exc.args[0]}"
        ) from exc


def run_external_reference_command(
    *,
    reference_id: str,
    command: str,
    output: Path,
    input_path: Path | None = None,
    config_path: Path | None = None,
    working_directory: Path | None = None,
    timeout_s: float | None = None,
) -> Path:
    """Run an explicit external-reference command and return its output path."""

    if not command.strip():
        raise ValueError("external reference command must not be empty")
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    rendered = _render_command(
        command,
        reference_id=reference_id,
        output=output,
        input_path=input_path.resolve() if input_path is not None else None,
        config_path=config_path.resolve() if config_path is not None else None,
    )
    completed = subprocess.run(
        rendered,
        cwd=None if working_directory is None else working_directory,
        shell=True,
        text=True,
        capture_output=True,
        timeout=timeout_s,
        check=False,
    )
    log_path = output.with_suffix(output.suffix + ".log")
    log_path.write_text(
        "\n".join(
            [
                f"reference_id={reference_id}",
                f"command={rendered}",
                f"returncode={completed.returncode}",
                "",
                "[stdout]",
                completed.stdout,
                "",
                "[stderr]",
                completed.stderr,
            ]
        ),
        encoding="utf-8",
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"{reference_id} external command failed with exit code "
            f"{completed.returncode}; see {log_path}"
        )
    if not output.exists():
        raise FileNotFoundError(
            f"{reference_id} external command did not produce output file: {output}"
        )
    return output
