"""Deterministic identity of the Monte Carlo numerical implementation."""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path


def monte_carlo_source_sha256() -> str:
    """Hash MC sources and the shared numerical modules they execute."""

    repository = Path(__file__).resolve().parents[3]
    package = Path(__file__).resolve().parent
    files = {
        *package.rglob("*.py"),
        repository / "electron_swarm" / "core" / "constants.py",
        repository / "electron_swarm" / "core" / "cross_sections.py",
        repository / "electron_swarm" / "core" / "scattering.py",
        repository / "electron_swarm" / "core" / "transport.py",
        repository / "electron_swarm" / "diagnostics" / "tail.py",
        repository / "electron_swarm" / "physics" / "angular_scattering.py",
        repository / "electron_swarm" / "physics" / "electron_neutral.py",
        repository / "electron_swarm" / "physics" / "kinetics.py",
    }
    digest = sha256()
    for path in sorted(files, key=lambda item: item.relative_to(repository).as_posix()):
        relative = path.relative_to(repository).as_posix().encode("utf-8")
        digest.update(relative + b"\0" + path.read_bytes() + b"\0")
    return digest.hexdigest()


__all__ = ["monte_carlo_source_sha256"]
