from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
from typing import Any

import pytest

from swarm_workflow.quality.solver import (
    PROPAGATOR_TARGET_QUALIFICATION_SCHEMA,
    PropagatorTargetRequirement,
    SolverQualification,
    SolverQualificationError,
    validate_propagator_target_qualification,
)
from tools import qualify_propagator_target as target_qualifier


ROOT = Path(__file__).resolve().parents[1]
CORE_QUALIFICATION = (
    ROOT
    / "docs"
    / "dev"
    / "results"
    / "propagator_p1_deterministic_qualification_20260908.json"
)
ICP_QUALIFICATION = (
    ROOT
    / "docs"
    / "dev"
    / "results"
    / "propagator_gec_icp_mean_energy_guarded_qualification_20260913.json"
)


def _icp_spec() -> target_qualifier.TargetRefinementSpec:
    return target_qualifier.TargetRefinementSpec(
        target="argon_gec_icp_restricted_local_mean_energy",
        fields_Td=(1500.0, 2500.0),
        operating_range_bracket_Td=(1500.0, 2500.0),
        medium_grid=(300, 48),
        fine_grid=(600, 72),
    )


def _core_bound_to_target(payload: dict[str, Any]) -> SolverQualification:
    entry = payload["inputs"]["core_qualification"]
    fingerprint = entry["implementation_fingerprint"]
    return SolverQualification(
        source_path=CORE_QUALIFICATION,
        sha256=entry["sha256"],
        schema=entry["schema"],
        decision=entry["decision"],
        implementation_fingerprint_schema=fingerprint["schema"],
        implementation_fingerprint_sha256=fingerprint["sha256"],
    )


def test_target_spec_represents_icp_refinement_without_gec_defaults() -> None:
    spec = _icp_spec()

    assert spec.fields_Td == (1500.0, 2500.0)
    assert spec.operating_range_bracket_Td == (1500.0, 2500.0)
    assert spec.table_support_cap_Td == 2500.0
    assert spec.medium_grid == (300, 48)
    assert spec.fine_grid == (600, 72)
    assert spec.family == "argon_gec_icp_restricted_local_mean_energy"


@pytest.mark.parametrize(
    "changes, match",
    [
        ({"fields_Td": (1500.0, 1500.0)}, "strictly increasing"),
        ({"operating_range_bracket_Td": (1500.0, 2000.0)}, "anchors"),
        ({"fine_grid": (300, 48)}, "refine"),
        ({"fine_grid": (299, 72)}, "coarser"),
        ({"fine_grid": (600, 71)}, "even polar"),
        ({"table_support_cap_Td": 1500.0}, "at or above"),
    ],
)
def test_target_spec_rejects_invalid_numerical_contract(
    changes: dict[str, Any],
    match: str,
) -> None:
    with pytest.raises(ValueError, match=match):
        replace(_icp_spec(), **changes)


def test_generic_cli_forwards_explicit_icp_contract_and_writes_artifact(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output = tmp_path / "nested" / "icp-target.json"
    captured: dict[str, Any] = {}

    def fake_run_target_qualification(**kwargs: Any) -> dict[str, Any]:
        captured.update(kwargs)
        return {
            "schema": target_qualifier.SCHEMA,
            "decision": {
                "target_refinement_qualified": True,
                "blocking_gates": [],
            },
        }

    monkeypatch.setattr(
        target_qualifier,
        "run_target_qualification",
        fake_run_target_qualification,
    )
    result = target_qualifier.main(
        [
            "--target",
            "argon_gec_icp_restricted_local_mean_energy",
            "--fields-Td",
            "1500",
            "2500",
            "--operating-bracket-Td",
            "1500",
            "2500",
            "--medium-grid",
            "300",
            "48",
            "--fine-grid",
            "600",
            "72",
            "--database",
            str(tmp_path / "icp.sqlite"),
            "--config",
            str(tmp_path / "propagator.yaml"),
            "--fine-max-memory-mb",
            "1024",
            "--memory-budget-mb",
            "2048",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    spec = captured["spec"]
    assert spec == _icp_spec()
    assert captured["database"] == (tmp_path / "icp.sqlite").resolve()
    assert captured["config_path"] == (tmp_path / "propagator.yaml").resolve()
    assert json.loads(output.read_text(encoding="utf-8"))["schema"] == (
        target_qualifier.SCHEMA
    )
    assert output.read_bytes().endswith(b"\n")


def test_generic_cli_forwards_explicit_gec_ccp_contract(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output = tmp_path / "gec-ccp-target.json"
    captured: dict[str, Any] = {}

    def fake_run_target_qualification(**kwargs: Any) -> dict[str, Any]:
        captured.update(kwargs)
        return {
            "schema": target_qualifier.SCHEMA,
            "decision": {"target_refinement_qualified": True},
        }

    monkeypatch.setattr(
        target_qualifier,
        "run_target_qualification",
        fake_run_target_qualification,
    )
    result = target_qualifier.main(
        [
            "--target",
            "argon_gec_ccp_restricted_local_mean_energy",
            "--fields-Td",
            "3000",
            "3500",
            "4000",
            "--operating-bracket-Td",
            "3000",
            "3500",
            "--table-support-cap-Td",
            "4000",
            "--medium-grid",
            "300",
            "48",
            "--fine-grid",
            "600",
            "72",
            "--mixture-id",
            "0",
            "--fine-max-memory-mb",
            "1024",
            "--family",
            "argon_gec_ccp_high_field",
            "--database",
            str(tmp_path / "gec.sqlite"),
            "--config",
            str(tmp_path / "gec.yaml"),
            "--workers",
            "3",
            "--memory-budget-mb",
            "3072",
            "--output",
            str(output),
        ]
    )

    assert result == 0
    spec = captured["spec"]
    assert spec.target == "argon_gec_ccp_restricted_local_mean_energy"
    assert spec.fields_Td == (3000.0, 3500.0, 4000.0)
    assert spec.operating_range_bracket_Td == (3000.0, 3500.0)
    assert spec.table_support_cap_Td == 4000.0
    assert spec.medium_grid == (300, 48)
    assert spec.fine_grid == (600, 72)
    assert spec.family == "argon_gec_ccp_high_field"
    assert spec.progress_label == "argon_gec_ccp_restricted_local_mean_energy"
    assert spec.release_claim == (
        "target_specific_P1_numerical_scope_only; downstream closure remains "
        "a separate decision"
    )
    assert json.loads(output.read_text(encoding="utf-8"))["schema"] == (
        target_qualifier.SCHEMA
    )


def test_generic_runner_uses_icp_fields_grids_and_family(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    spec = replace(_icp_spec(), family="argon_gec_icp_high_field")
    config = tmp_path / "propagator.yaml"
    submitted: list[tuple[Any, ...]] = []
    executor_workers: list[int] = []

    class FakeCore:
        source_path = CORE_QUALIFICATION

        def identity_entry(self) -> dict[str, Any]:
            return {
                "schema": "core",
                "sha256": "c" * 64,
                "decision": "p1_deterministic_core_qualified",
                "implementation_fingerprint": {
                    "schema": "core-source",
                    "sha256": "f" * 64,
                },
            }

    class FakeFuture:
        def __init__(self, field: float) -> None:
            self.field = field

        def result(self) -> tuple[dict[str, Any], object]:
            return (
                {
                    "E_over_N_Td": self.field,
                    "status": "ok",
                    "elapsed_s": self.field / 1000.0,
                    "quality_gates_passed": True,
                    "energy_cells": 600,
                    "polar_cells": 72,
                },
                object(),
            )

    class FakeExecutor:
        def __init__(self, *, max_workers: int, mp_context: Any) -> None:
            del mp_context
            executor_workers.append(max_workers)

        def __enter__(self) -> FakeExecutor:
            return self

        def __exit__(self, *_args: Any) -> None:
            return None

        def submit(self, _function: Any, *args: Any) -> FakeFuture:
            submitted.append(args)
            return FakeFuture(float(args[1]))

    monkeypatch.setattr(
        target_qualifier,
        "validate_propagator_core_qualification",
        lambda _path: FakeCore(),
    )
    monkeypatch.setattr(
        target_qualifier,
        "_read_medium_cases",
        lambda _database, selected: (
            {field: object() for field in selected.fields_Td},
                {
                    "base_config_sha256": "a" * 64,
                    "workflow_config_sha256": "b" * 64,
                    "cross_sections_sha256": "d" * 64,
                    "propagator_solver_source_sha256": "f" * 64,
                },
            "e" * 64,
        ),
    )
    monkeypatch.setattr(target_qualifier, "file_sha256", lambda _path: "a" * 64)
    monkeypatch.setattr(target_qualifier, "load_config", lambda _path: object())
    monkeypatch.setattr(
        target_qualifier,
        "_cross_sections_sha256",
        lambda _config: "d" * 64,
    )
    target_fingerprint = {
        "schema": "target-source",
        "sha256": "9" * 64,
        "file_count": 2,
        "files": ["first", "second"],
    }
    monkeypatch.setattr(
        target_qualifier,
        "propagator_target_qualifier_fingerprint",
        lambda _root: target_fingerprint,
    )
    monkeypatch.setattr(
        target_qualifier,
        "_require_unchanged_target_inputs",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(target_qualifier, "ProcessPoolExecutor", FakeExecutor)
    monkeypatch.setattr(
        target_qualifier.multiprocessing,
        "get_context",
        lambda _method: object(),
    )
    monkeypatch.setattr(
        target_qualifier,
        "as_completed",
        lambda futures: reversed(tuple(futures)),
    )
    monkeypatch.setattr(
        target_qualifier,
        "_comparison",
        lambda *_args, **_kwargs: {
            "status": "available",
            "passed": True,
        },
    )

    payload = target_qualifier.run_target_qualification(
        spec=spec,
        database=tmp_path / "icp.sqlite",
        config_path=config,
        core_qualification_path=CORE_QUALIFICATION,
        workers=4,
        memory_budget_mb=2048,
    )

    assert executor_workers == [2]
    assert [row["E_over_N_Td"] for row in payload["fine_runs"]] == [
        1500.0,
        2500.0,
    ]
    assert [args[1:] for args in submitted] == [
        (
            1500.0,
            (600, 72),
            1024,
            "argon_gec_icp_high_field",
            "a" * 64,
            "d" * 64,
            "f" * 64,
            "9" * 64,
        ),
        (
            2500.0,
            (600, 72),
            1024,
            "argon_gec_icp_high_field",
            "a" * 64,
            "d" * 64,
            "f" * 64,
            "9" * 64,
        ),
    ]
    assert payload["scope"]["fields_Td"] == [1500.0, 2500.0]
    assert payload["scope"]["medium_grid"] == [300, 48]
    assert payload["scope"]["fine_grid"] == [600, 72]
    assert payload["scope"]["mixture_id"] == 0
    assert payload["decision"] == {
        "fine_case_quality_passed": True,
        "medium_fine_refinement_passed": True,
        "target_refinement_qualified": True,
        "blocking_gates": [],
        "release_claim": (
            "target_specific_P1_numerical_scope_only; downstream closure "
            "remains a separate decision"
        ),
    }


def test_target_source_guard_rejects_changed_medium_evidence(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    core = SolverQualification(
        source_path=CORE_QUALIFICATION,
        sha256="c" * 64,
        schema="core",
        decision="p1_deterministic_core_qualified",
        implementation_fingerprint_schema="core-source",
        implementation_fingerprint_sha256="f" * 64,
    )
    target_fingerprint = {
        "schema": "target-source",
        "sha256": "9" * 64,
        "file_count": 2,
        "files": ["first", "second"],
    }
    metadata = {
        "base_config_sha256": "a" * 64,
        "cross_sections_sha256": "d" * 64,
    }
    monkeypatch.setattr(
        target_qualifier,
        "validate_propagator_core_qualification",
        lambda _path: core,
    )
    monkeypatch.setattr(target_qualifier, "file_sha256", lambda _path: "a" * 64)
    monkeypatch.setattr(target_qualifier, "load_config", lambda _path: object())
    monkeypatch.setattr(
        target_qualifier,
        "_cross_sections_sha256",
        lambda _config: "d" * 64,
    )
    monkeypatch.setattr(
        target_qualifier,
        "propagator_target_medium_evidence",
        lambda *_args, **_kwargs: ("e" * 64, metadata),
    )
    monkeypatch.setattr(
        target_qualifier,
        "propagator_target_qualifier_fingerprint",
        lambda _root: target_fingerprint,
    )
    starting_state = {
        "core_qualification": core.identity_entry(),
        "base_config_sha256": "a" * 64,
        "cross_sections_sha256": "d" * 64,
        "medium_evidence_sha256": "0" * 64,
        "medium_metadata": metadata,
        "target_qualifier_fingerprint": target_fingerprint,
    }

    with pytest.raises(RuntimeError, match="changed during execution"):
        target_qualifier._require_unchanged_target_inputs(
            starting_state,
            spec=_icp_spec(),
            database=tmp_path / "icp.sqlite",
            config_path=tmp_path / "propagator.yaml",
            core_qualification_path=CORE_QUALIFICATION,
        )


def test_generic_target_artifact_is_validated_without_gec_field_assumptions() -> None:
    payload = json.loads(ICP_QUALIFICATION.read_text(encoding="utf-8"))
    core = _core_bound_to_target(payload)
    requirement = PropagatorTargetRequirement(
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
    with pytest.raises(SolverQualificationError, match="explicit target requirement"):
        validate_propagator_target_qualification(
            ICP_QUALIFICATION,
            core_qualification=core,
        )
    evidence = validate_propagator_target_qualification(
        ICP_QUALIFICATION,
        core_qualification=core,
        requirement=requirement,
    )

    assert evidence.schema == PROPAGATOR_TARGET_QUALIFICATION_SCHEMA
    assert evidence.target == "argon_gec_icp_restricted_local_mean_energy"
    assert evidence.fields_Td == (1500.0, 2500.0, 2744.2936035)
    assert evidence.operating_range_bracket_Td == (1500.0, 2500.0)
    assert evidence.table_support_cap_Td == 2744.2936035
    assert evidence.medium_grid == (300, 48)
    assert evidence.fine_grid == (600, 72)


def test_legacy_gec_target_schema_is_rejected(tmp_path: Path) -> None:
    payload = json.loads(ICP_QUALIFICATION.read_text(encoding="utf-8"))
    core = _core_bound_to_target(payload)
    payload["schema"] = "swarm.propagator_gec_target_qualification.v1"
    artifact = tmp_path / "legacy-target.json"
    artifact.write_text(json.dumps(payload), encoding="utf-8")

    with pytest.raises(SolverQualificationError, match="unsupported.*schema"):
        validate_propagator_target_qualification(
            artifact,
            core_qualification=core,
            allow_self_described_generic=True,
        )
