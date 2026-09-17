"""Persistence-bound quality-policy provenance for workflow databases."""

from __future__ import annotations

import json
import sqlite3

from ..quality.policy import (
    QUALITY_THRESHOLDS_METADATA_KEY,
    QualityThresholds,
    parse_quality_thresholds,
    quality_thresholds_json,
)
from . import repository as _repository
from .store import WorkflowSchemaError as _WorkflowSchemaError


def resolve_quality_thresholds(
    connection: sqlite3.Connection,
    requested: QualityThresholds | None = None,
) -> QualityThresholds:
    """Restore and verify the immutable quality policy for a workflow database."""

    metadata = _repository.read_metadata(connection)
    raw = metadata.get(QUALITY_THRESHOLDS_METADATA_KEY)
    if raw is None:
        if _repository.database_has_workflow_data(connection):
            raise _WorkflowSchemaError(
                "Workflow provenance is missing quality_thresholds_json; "
                "regenerate or explicitly migrate the database"
            )
        return requested or QualityThresholds()

    try:
        payload = json.loads(raw)
        stored = parse_quality_thresholds(payload)
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        raise _WorkflowSchemaError(
            "Workflow provenance has invalid quality_thresholds_json"
        ) from exc
    canonical = quality_thresholds_json(stored)
    if raw != canonical:
        raise _WorkflowSchemaError(
            "Workflow provenance quality_thresholds_json is not canonical; "
            "regenerate or explicitly migrate the database"
        )
    if requested is not None and quality_thresholds_json(requested) != canonical:
        raise _WorkflowSchemaError(
            "Requested quality thresholds differ from workflow database provenance"
        )
    return stored


def source_quality_thresholds_json(connection: sqlite3.Connection) -> str:
    """Return the exact validated immutable source-policy encoding."""

    source = resolve_quality_thresholds(connection)
    return _repository.read_metadata(connection).get(
        QUALITY_THRESHOLDS_METADATA_KEY,
        quality_thresholds_json(source),
    )


def resolve_evaluation_quality_thresholds(
    connection: sqlite3.Connection,
) -> QualityThresholds:
    """Return the last explicit aggregate policy, or immutable source policy."""

    source = resolve_quality_thresholds(connection)
    if not _repository.table_exists(connection, "aggregate_quality"):
        return source
    columns = {
        str(row[1])
        for row in connection.execute(
            "PRAGMA table_info(aggregate_quality)"
        ).fetchall()
    }
    if "thresholds_json" not in columns:
        return source
    rows = connection.execute(
        "SELECT DISTINCT thresholds_json FROM aggregate_quality"
    ).fetchall()
    if not rows:
        return source
    if len(rows) != 1 or not isinstance(rows[0][0], str):
        raise _WorkflowSchemaError(
            "aggregate_quality contains inconsistent evaluation policies"
        )
    raw = str(rows[0][0])
    try:
        evaluation = parse_quality_thresholds(json.loads(raw))
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        raise _WorkflowSchemaError(
            "aggregate_quality contains an invalid evaluation policy"
        ) from exc
    if quality_thresholds_json(evaluation) != raw:
        raise _WorkflowSchemaError(
            "aggregate_quality evaluation policy is not canonical"
        )
    return evaluation


def validate_aggregate_quality_thresholds(
    connection: sqlite3.Connection,
    thresholds: QualityThresholds,
) -> None:
    """Validate source and evaluation policies recorded on aggregate rows."""

    if not _repository.table_exists(connection, "aggregate_quality"):
        return
    expected = quality_thresholds_json(thresholds)
    source_expected = source_quality_thresholds_json(connection)
    reevaluated = int(source_expected != expected)
    columns = {
        str(row[1])
        for row in connection.execute(
            "PRAGMA table_info(aggregate_quality)"
        ).fetchall()
    }
    required_columns = {
        "source_thresholds_json",
        "thresholds_json",
        "quality_policy_reevaluated",
    }
    if not required_columns.issubset(columns):
        raise _WorkflowSchemaError(
            "aggregate_quality lacks source/evaluation quality provenance; "
            "rerun aggregation"
        )
    rows = connection.execute(
        """
        SELECT DISTINCT source_thresholds_json, thresholds_json,
                        quality_policy_reevaluated
        FROM aggregate_quality
        """
    ).fetchall()
    for row in rows:
        if (
            row[0] != source_expected
            or row[1] != expected
            or int(row[2]) != reevaluated
        ):
            raise _WorkflowSchemaError(
                "aggregate_quality source/evaluation policy provenance is "
                "inconsistent; rerun aggregation"
            )
