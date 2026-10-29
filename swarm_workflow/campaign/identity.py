"""Portable identity of workflow projection, seeding, persistence, and aggregation."""

from hashlib import sha256
from importlib.metadata import version
from pathlib import Path

WORKFLOW_IMPLEMENTATION_KEY = "workflow_implementation_json"


def workflow_execution_identity(repository_root: Path | None = None) -> dict:
    directory = (
        Path(repository_root) / "swarm_workflow" / "campaign"
        if repository_root is not None
        else Path(__file__).resolve().parent
    )
    files = sorted(directory.glob("*.py"))
    if not files:
        raise FileNotFoundError(f"workflow sources are missing: {directory}")
    digest = sha256()
    for path in files:
        digest.update(
            path.name.encode()
            + b"\0"
            + path.read_bytes().replace(b"\r\n", b"\n")
            + b"\0"
        )
    shared_statistics = (
        directory.parents[1] / "electron_swarm" / "core" / "statistics.py"
    )
    digest.update(
        b"electron_swarm/core/statistics.py\0"
        + shared_statistics.read_bytes().replace(b"\r\n", b"\n")
        + b"\0"
    )
    return {
        "schema": "swarm.workflow_execution_sha256.v2",
        "sha256": digest.hexdigest(),
        "runtime": {"pandas": version("pandas")},
    }
