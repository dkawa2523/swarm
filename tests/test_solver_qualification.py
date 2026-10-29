from __future__ import annotations

import json
from pathlib import Path

import pytest

from swarm_workflow.quality.solver import (
    PROPAGATOR_CORE_QUALIFICATION_SCHEMA,
    PROPAGATOR_TARGET_QUALIFICATION_SCHEMA,
    PropagatorTargetRequirement,
    SolverQualification,
    SolverQualificationError,
    validate_propagator_core_qualification,
    validate_propagator_target_qualification,
)


ROOT = Path(__file__).resolve().parents[1]
QUALIFICATION = (
    ROOT
    / "docs"
    / "dev"
    / "results"
    / "propagator_p1_core_qualification.json"
)
TARGET_QUALIFICATION = (
    ROOT
    / "docs"
    / "dev"
    / "results"
    / "propagator_gec_icp_target_qualification.json"
)
TARGET_REQUIREMENT = PropagatorTargetRequirement(
    target="argon_gec_icp_restricted_local_mean_energy",
    fields_Td=(1500.0, 2500.0, 2744.2936035),
    operating_range_bracket_Td=(1500.0, 2500.0),
    table_support_cap_Td=2744.2936035,
    mixture_id=0,
    medium_grid=(300, 48),
    fine_grid=(600, 72),
    max_scalar_refinement_limit=0.01,
    max_eedf_weighted_L1_limit=0.02,
)


def _target_bound_core(payload: dict[str, object]) -> SolverQualification:
    inputs = payload["inputs"]
    assert isinstance(inputs, dict)
    entry = inputs["core_qualification"]
    assert isinstance(entry, dict)
    fingerprint = entry["implementation_fingerprint"]
    assert isinstance(fingerprint, dict)
    return SolverQualification(
        source_path=QUALIFICATION,
        sha256=str(entry["sha256"]),
        schema=str(entry["schema"]),
        decision=str(entry["decision"]),
        implementation_fingerprint_schema=str(fingerprint["schema"]),
        implementation_fingerprint_sha256=str(fingerprint["sha256"]),
    )


def test_current_propagator_qualification_matches_source_and_inputs() -> None:
    evidence = validate_propagator_core_qualification(QUALIFICATION)

    assert evidence.schema == PROPAGATOR_CORE_QUALIFICATION_SCHEMA
    assert evidence.decision == "p1_deterministic_core_qualified"
    assert len(evidence.sha256) == 64
    assert len(evidence.implementation_fingerprint_sha256) == 64


def test_propagator_qualification_rejects_failed_decision(
    tmp_path: Path,
) -> None:
    payload = json.loads(QUALIFICATION.read_text(encoding="utf-8"))
    payload["decision"]["p1_deterministic_core_qualified"] = False
    path = tmp_path / "failed.json"
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(
        SolverQualificationError,
        match="deterministic P1 core is not qualified",
    ):
        validate_propagator_core_qualification(path)


def test_propagator_qualification_rejects_quick_run(tmp_path: Path) -> None:
    payload = json.loads(QUALIFICATION.read_text(encoding="utf-8"))
    payload["scope"]["quick"] = True
    path = tmp_path / "quick.json"
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(
        SolverQualificationError,
        match="not a full deterministic P0/P1 run",
    ):
        validate_propagator_core_qualification(path)


def test_propagator_qualification_rejects_stale_implementation(
    tmp_path: Path,
) -> None:
    payload = json.loads(QUALIFICATION.read_text(encoding="utf-8"))
    payload["environment"]["implementation_fingerprint"]["sha256"] = "0" * 64
    path = tmp_path / "stale.json"
    path.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(
        SolverQualificationError,
        match="implementation changed after qualification",
    ):
        validate_propagator_core_qualification(path)


def test_propagator_qualification_rejects_altered_input_snapshot(tmp_path: Path) -> None:
    payload = json.loads(QUALIFICATION.read_text(encoding="utf-8"))
    inputs = payload["environment"]["input_artifacts"]
    inputs[next(iter(inputs))] += "\n# changed\n"
    path = tmp_path / "altered-input.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(SolverQualificationError, match="input snapshot changed"):
        validate_propagator_core_qualification(path)


def test_current_propagator_target_qualification_is_internally_valid() -> None:
    payload = json.loads(TARGET_QUALIFICATION.read_text(encoding="utf-8"))
    core = _target_bound_core(payload)
    evidence = validate_propagator_target_qualification(
        TARGET_QUALIFICATION,
        core_qualification=core,
        requirement=TARGET_REQUIREMENT,
    )

    assert evidence.schema == PROPAGATOR_TARGET_QUALIFICATION_SCHEMA
    assert evidence.decision == "target_refinement_qualified"
    assert evidence.fields_Td == (1500.0, 2500.0, 2744.2936035)
    assert evidence.medium_grid == (300, 48)
    assert evidence.fine_grid == (600, 72)


def test_propagator_target_qualification_rejects_failed_decision(
    tmp_path: Path,
) -> None:
    payload = json.loads(TARGET_QUALIFICATION.read_text(encoding="utf-8"))
    payload["decision"]["target_refinement_qualified"] = False
    path = tmp_path / "failed-target.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    core = _target_bound_core(payload)

    with pytest.raises(
        SolverQualificationError,
        match="target refinement is not qualified",
    ):
        validate_propagator_target_qualification(
            path,
            core_qualification=core,
            requirement=TARGET_REQUIREMENT,
        )


def test_propagator_target_qualification_rejects_bundle_provenance() -> None:
    payload = json.loads(TARGET_QUALIFICATION.read_text(encoding="utf-8"))
    core = _target_bound_core(payload)

    with pytest.raises(
        SolverQualificationError,
        match="provenance differs from the downstream manifest",
    ):
        validate_propagator_target_qualification(
            TARGET_QUALIFICATION,
            core_qualification=core,
            requirement=TARGET_REQUIREMENT,
            provenance_hashes={
                "calculation_sha256": "0" * 64,
                "workflow_config_sha256": "0" * 64,
                "cross_sections_sha256": "0" * 64,
            },
        )


@pytest.mark.parametrize("damage", ["nan", "missing_gate", "negative"])
def test_target_refinement_rejects_invalid_numerical_evidence(tmp_path: Path, damage: str) -> None:
    payload = json.loads(TARGET_QUALIFICATION.read_text(encoding="utf-8"))
    row = payload["medium_fine_refinement"][0]
    if damage == "missing_gate":
        row["gated_scalars"] = ["mean_energy"]
    else:
        row["eedf_weighted_L1"] = float("nan") if damage == "nan" else -1.0
    artifact = tmp_path / "invalid-target.json"
    artifact.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(SolverQualificationError, match="refinement (is malformed|failed)"):
        validate_propagator_target_qualification(
            artifact, core_qualification=_target_bound_core(payload),
            requirement=TARGET_REQUIREMENT,
        )
