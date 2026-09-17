from __future__ import annotations

from hashlib import sha256
import json

import pytest

import swarm_workflow.campaign.aggregate as aggregate_module
import swarm_workflow.campaign.quality as database_policy
import swarm_workflow.quality.policy as policy_module
from swarm_workflow.quality.policy import (
    QualityThresholds,
    RequiredRateRse,
    parse_quality_thresholds,
    quality_policy_provenance,
    quality_thresholds_json,
    quality_thresholds_payload,
    quality_thresholds_sha256,
)


def test_quality_policy_has_one_canonical_codec() -> None:
    policy = QualityThresholds(
        mobility_rse=0.1,
        required_rate_rse=RequiredRateRse(
            process_type={"ionization": 0.2},
            process={"Ar excitation": None},
        ),
    )
    payload = quality_thresholds_payload(policy)
    encoded = quality_thresholds_json(policy)

    assert parse_quality_thresholds(payload) == policy
    assert encoded == json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    assert quality_thresholds_sha256(policy) == sha256(
        encoded.encode("utf-8")
    ).hexdigest()
    assert quality_policy_provenance(
        policy,
        policy,
        source_json=encoded,
    ) == {
        "source": payload,
        "source_sha256": sha256(encoded.encode("utf-8")).hexdigest(),
        "evaluation": payload,
        "evaluation_sha256": quality_thresholds_sha256(policy),
        "quality_only_reevaluation": False,
    }


@pytest.mark.parametrize(
    "name",
    (
        "QUALITY_THRESHOLDS_METADATA_KEY",
        "RequiredRateRse",
        "QualityThresholds",
        "parse_quality_thresholds",
        "quality_policy_provenance",
        "quality_thresholds_json",
        "quality_thresholds_payload",
        "quality_thresholds_sha256",
        "resolve_evaluation_quality_thresholds",
        "resolve_quality_thresholds",
        "source_quality_thresholds_json",
        "validate_aggregate_quality_thresholds",
    ),
)
def test_aggregate_does_not_reexport_quality_policy(name: str) -> None:
    assert not hasattr(aggregate_module, name)


@pytest.mark.parametrize(
    "name",
    (
        "resolve_evaluation_quality_thresholds",
        "resolve_quality_thresholds",
        "source_quality_thresholds_json",
        "validate_aggregate_quality_thresholds",
    ),
)
def test_database_quality_policy_has_one_campaign_owner(name: str) -> None:
    assert hasattr(database_policy, name)
    assert not hasattr(policy_module, name)


def test_pure_quality_policy_has_no_persistence_dependencies() -> None:
    assert "sqlite3" not in policy_module.__dict__
    assert "_repository" not in policy_module.__dict__
    assert "_WorkflowSchemaError" not in policy_module.__dict__


def test_required_rate_names_remain_case_insensitively_unique() -> None:
    with pytest.raises(ValueError, match="duplicate case-insensitive name"):
        parse_quality_thresholds(
            {
                "required_rate_rse": {
                    "process": {"Ar excitation": 0.1, "ar EXCITATION": 0.2}
                }
            }
        )
