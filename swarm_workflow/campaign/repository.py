"""Shared read-only queries for workflow SQLite repositories."""

from __future__ import annotations

import sqlite3


WORKFLOW_DATA_TABLES = (
    "metadata",
    "mixtures",
    "cases",
    "rates",
    "eedf_bins",
    "aggregate_scalars",
    "aggregate_eedf_bins",
    "aggregate_quality",
)


def table_exists(connection: sqlite3.Connection, name: str) -> bool:
    row = connection.execute(
        "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ?",
        (name,),
    ).fetchone()
    return row is not None


def read_metadata(connection: sqlite3.Connection) -> dict[str, str]:
    if not table_exists(connection, "metadata"):
        return {}
    return {
        str(row[0]): str(row[1])
        for row in connection.execute("SELECT key, value FROM metadata")
    }


def database_has_workflow_data(connection: sqlite3.Connection) -> bool:
    for table_name in WORKFLOW_DATA_TABLES:
        if not table_exists(connection, table_name):
            continue
        if connection.execute(f"SELECT 1 FROM {table_name} LIMIT 1").fetchone():
            return True
    return False
