from pathlib import Path
import shutil

from electron_swarm.core.solver_ids import CANONICAL_SOLVER_IDS
from electron_swarm.core.source_identity import solver_source_fingerprint
from swarm_workflow.quality.solver import (
    propagator_target_qualifier_fingerprint,
    validate_propagator_core_qualification,
)


ROOT = Path(__file__).resolve().parents[1]


def test_runtime_identity_and_qualification_need_no_repository_tools(tmp_path: Path) -> None:
    """Reproduce the installed-wheel layout without examples, tests, or tools."""
    for package in ("electron_swarm", "swarm_workflow"):
        shutil.copytree(ROOT / package, tmp_path / package,
                        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    for solver in CANONICAL_SOLVER_IDS:
        assert solver_source_fingerprint(solver, tmp_path) == solver_source_fingerprint(solver)
    assert propagator_target_qualifier_fingerprint(tmp_path) == (
        propagator_target_qualifier_fingerprint(ROOT)
    )
    validate_propagator_core_qualification(
        ROOT / "docs/dev/results/propagator_p1_core_qualification.json",
        repository_root=tmp_path,
    )
