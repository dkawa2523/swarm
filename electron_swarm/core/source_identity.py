"""Portable numerical-source identities for all canonical solvers.

Only installed runtime sources belong here. Research tests, qualifier tools,
documentation, and local repository paths are not runtime dependencies.
"""

from __future__ import annotations

from hashlib import sha256
import json
from importlib.metadata import version
from pathlib import Path
import platform
from typing import Iterable

from .solver_ids import CANONICAL_SOLVER_IDS


SOLVER_SOURCE_FINGERPRINT_SCHEMA = "swarm.numerical_model_sha256.v2"
EXECUTION_SOURCE_SCHEMA = "swarm.execution_pipeline_sha256.v1"
SOLVER_IMPLEMENTATIONS_KEY = "solver_implementations_json"
PROPAGATOR_SOLVER_SOURCE_METADATA_KEY = "propagator_solver_source_sha256"


def solver_source_fingerprint(
    solver: str,
    repository_root: str | Path | None = None,
) -> dict[str, object]:
    if solver not in CANONICAL_SOLVER_IDS:
        raise ValueError(f"unknown solver: {solver!r}")
    root = (
        Path(repository_root).resolve()
        if repository_root is not None
        else Path(__file__).resolve().parents[2]
    )
    package = root / "electron_swarm"
    directories = [package / name for name in ("physics", "grids")]
    directories.append(package / "solvers" / solver)
    if solver in {"two_term", "multi_term"}:
        directories.extend(
            (package / "collisions", package / "solvers" / "boltzmann_common")
        )
    files = {
        package / "core" / name
        for name in (
            "constants.py",
            "gas.py",
            "numerics.py",
            "cross_sections.py",
            "scattering.py",
            "solver_configs.py",
            "transport.py",
        )
    }
    files.add(package / "solvers" / "base.py")
    for directory in directories:
        if not directory.is_dir():
            raise FileNotFoundError(f"solver source directory is missing: {directory}")
        files.update(directory.rglob("*.py"))
    relative_files = sorted(path.relative_to(root).as_posix() for path in files)
    digest = sha256()
    for relative in relative_files:
        # Git newline conversion must not invalidate otherwise identical code.
        content = (root / relative).read_bytes().replace(b"\r\n", b"\n")
        digest.update(relative.encode("utf-8") + b"\0" + content + b"\0")
    return {
        "schema": SOLVER_SOURCE_FINGERPRINT_SCHEMA,
        "sha256": digest.hexdigest(),
        "file_count": len(files),
        "files": relative_files,
    }


def execution_source_fingerprint(
    repository_root: str | Path | None = None,
) -> dict[str, str]:
    package = (
        Path(repository_root).resolve() / "electron_swarm"
        if repository_root is not None
        else Path(__file__).resolve().parents[1]
    )
    execution_files = [package / "runner.py"]
    execution_files.extend(
        path
        for path in (package / "orchestration").glob("*.py")
        if path.name != "comparison.py"
    )
    execution_files.extend((package / "diagnostics").rglob("*.py"))
    execution_files.extend(
        package / "core" / name
        for name in (
            "config.py",
            "config_parser.py",
            "config_validation.py",
            "capabilities.py",
            "feature_resolution.py",
            "solver_registry.py",
            "physical_context.py",
            "identity.py",
            "source_identity.py",
            "yaml_input.py",
            "observables.py",
            "results.py",
            "result_metadata.py",
        )
    )
    execution_files.extend(
        path
        for path in (package / "core" / "config_sections").glob("*.py")
        if path.name != "output.py"
    )
    digest = sha256()
    for path in sorted(set(execution_files)):
        digest.update(
            path.relative_to(package).as_posix().encode()
            + b"\0"
            + path.read_bytes().replace(b"\r\n", b"\n")
            + b"\0"
        )
    return {"schema": EXECUTION_SOURCE_SCHEMA, "sha256": digest.hexdigest()}


def solver_implementations(solvers: Iterable[str]) -> dict[str, object]:
    execution_sha256 = execution_source_fingerprint()["sha256"]
    runtime = {
        "python": platform.python_version(),
        **{name: version(name) for name in ("numpy", "scipy", "lxcat_data_parser")},
    }
    return {
        solver: {
            "model_schema": SOLVER_SOURCE_FINGERPRINT_SCHEMA,
            "source_sha256": solver_source_fingerprint(solver)["sha256"],
            "execution_schema": EXECUTION_SOURCE_SCHEMA,
            "execution_sha256": execution_sha256,
            "runtime": {
                **runtime,
                **({"numba": version("numba")} if solver == "monte_carlo" else {}),
            },
        }
        for solver in sorted(set(solvers))
    }


def input_file_hashes(config: object) -> dict[str, str]:
    """Bind calculations to supplied collision and angular data, not paths alone."""

    paths = {Path(item.path).resolve() for item in config.cross_sections.files}
    moment_table = config.physics.angular_scattering.moment_table
    if moment_table is not None:
        paths.add(Path(moment_table.path).resolve())
    paths.update(
        path.with_suffix(".provenance.json")
        for path in list(paths)
        if path.with_suffix(".provenance.json").is_file()
    )
    return {str(path): sha256(path.read_bytes()).hexdigest() for path in sorted(paths)}


def input_evidence(config: object, hashes: dict[str, str]) -> dict:
    """Report supplied attribution separately from numerical execution evidence."""
    evidence = {}
    for item in config.cross_sections.files:
        path = Path(item.path).resolve()
        sidecar = path.with_suffix(".provenance.json")
        record = {"attribution_status": "not_provided"}
        if str(sidecar) in hashes:
            payload = json.loads(sidecar.read_text(encoding="utf-8"))
            source = payload.get("source", {})
            matches = payload.get("output", {}).get("sha256") == hashes[str(path)]
            record = {
                "attribution_status": "declared"
                if matches and source.get("url")
                else "unattributed"
                if matches
                else "hash_mismatch",
                "source": source,
                "provenance_file": str(sidecar),
            }
        evidence[str(path)] = record
    return evidence
