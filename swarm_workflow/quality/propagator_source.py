"""Canonical source identity for deterministic Propagator qualification."""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
from typing import Any


PROPAGATOR_SOURCE_FINGERPRINT_SCHEMA = (
    "swarm.propagator_source_tree_sha256.v1"
)
PROPAGATOR_SOLVER_SOURCE_METADATA_KEY = "propagator_solver_source_sha256"


def propagator_qualification_source_fingerprint(
    repository_root: str | Path,
) -> dict[str, Any]:
    """Hash the complete numerical and qualification source boundary."""

    root = Path(repository_root).resolve()
    files = sorted(
        {
            *(root / "electron_swarm" / "solvers" / "propagator").glob("*.py"),
            *(root / "electron_swarm" / "core" / "config_sections").glob("*.py"),
            root / "electron_swarm" / "core" / "config.py",
            root / "electron_swarm" / "core" / "config_parser.py",
            root / "electron_swarm" / "core" / "config_validation.py",
            root / "electron_swarm" / "core" / "constants.py",
            root / "electron_swarm" / "core" / "cross_sections.py",
            root / "electron_swarm" / "core" / "legacy_schema.py",
            root / "electron_swarm" / "core" / "results.py",
            root / "electron_swarm" / "core" / "scattering.py",
            root / "electron_swarm" / "core" / "solver_configs.py",
            root / "electron_swarm" / "core" / "solver_ids.py",
            root / "electron_swarm" / "core" / "transport.py",
            root / "electron_swarm" / "physics" / "angular_scattering.py",
            root / "electron_swarm" / "physics" / "electron_neutral.py",
            root / "electron_swarm" / "physics" / "kinetics.py",
            root / "electron_swarm" / "solvers" / "base.py",
            root / "tests" / "support" / "propagator_reference.py",
            root / "tests" / "test_propagator_operator.py",
            root / "tests" / "test_propagator_solver.py",
            root / "tools" / "qualify_propagator_p1_deterministic.py",
            Path(__file__).resolve(),
        },
        key=lambda item: item.relative_to(root).as_posix(),
    )
    missing = [path for path in files if not path.is_file()]
    if missing:
        names = ", ".join(path.relative_to(root).as_posix() for path in missing)
        raise FileNotFoundError(f"Propagator qualification source is missing: {names}")

    digest = sha256()
    for path in files:
        relative = path.relative_to(root).as_posix()
        digest.update(relative.encode("utf-8") + b"\0" + path.read_bytes() + b"\0")
    return {
        "schema": PROPAGATOR_SOURCE_FINGERPRINT_SCHEMA,
        "sha256": digest.hexdigest(),
        "file_count": len(files),
        "files": [path.relative_to(root).as_posix() for path in files],
    }


__all__ = [
    "PROPAGATOR_SOURCE_FINGERPRINT_SCHEMA",
    "PROPAGATOR_SOLVER_SOURCE_METADATA_KEY",
    "propagator_qualification_source_fingerprint",
]
