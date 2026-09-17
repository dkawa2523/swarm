from __future__ import annotations

import pytest

from swarm_workflow.quality.table import (
    COMMON_QUALITY_COLUMNS,
    DETERMINISTIC_QUALITY_EVIDENCE_COLUMNS,
    MONTE_CARLO_QUALITY_EVIDENCE_COLUMNS,
    quality_table_schema,
)


def test_deterministic_quality_schemas_exclude_monte_carlo_evidence() -> None:
    two_term = quality_table_schema("two_term")
    propagator = quality_table_schema("propagator")

    assert two_term.columns == propagator.columns
    assert two_term.evidence_kind == "deterministic_convergence"
    assert set(COMMON_QUALITY_COLUMNS) < set(two_term.columns)
    assert set(DETERMINISTIC_QUALITY_EVIDENCE_COLUMNS) < set(two_term.columns)
    assert set(MONTE_CARLO_QUALITY_EVIDENCE_COLUMNS).isdisjoint(
        two_term.columns
    )


def test_monte_carlo_quality_schema_excludes_deterministic_evidence() -> None:
    schema = quality_table_schema("monte_carlo")

    assert schema.evidence_kind == "independent_replica_statistics"
    assert set(COMMON_QUALITY_COLUMNS) < set(schema.columns)
    assert set(MONTE_CARLO_QUALITY_EVIDENCE_COLUMNS) < set(schema.columns)
    assert set(DETERMINISTIC_QUALITY_EVIDENCE_COLUMNS).isdisjoint(
        schema.columns
    )


def test_quality_schema_rejects_unknown_source() -> None:
    with pytest.raises(ValueError, match="unsupported quality-table source"):
        quality_table_schema("legacy")
