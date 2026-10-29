"""Typed CSV without inference, index alignment, or silent scalar coercion."""

from dataclasses import dataclass
import csv
from io import StringIO
from math import isfinite
from numbers import Integral, Real
from pathlib import Path

import numpy as np
import pandas as pd


PRODUCT_ARTIFACT_SCHEMA = "swarm.product_run.v3"
FIGURE_ARTIFACT_SCHEMA = "swarm.figures.v1"
_KINDS = {"string", "boolean", "integer", "number"}


@dataclass(frozen=True, slots=True)
class ProductRunArtifact:
    metadata: dict
    tables: dict[str, pd.DataFrame]


def json_scalar(value):
    """Keep numeric/boolean metadata typed instead of converting it to text."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"unsupported artifact metadata: {type(value).__name__}")


def _validate_columns(columns: dict[str, str]) -> None:
    if not columns or any(not isinstance(name, str) or not name for name in columns):
        raise ValueError("product table requires nonempty, unique column names")
    if any(
        not isinstance(kind, str) or kind not in _KINDS for kind in columns.values()
    ):
        raise ValueError("unsupported product column kind")


def _encode(value: object, kind: str) -> str:
    if not pd.api.types.is_scalar(value):
        raise ValueError("product table cells must be scalars")
    if (
        value is None
        or value is pd.NA
        or (not isinstance(value, str) and pd.isna(value))
    ):
        return ""
    if kind == "string":
        if not isinstance(value, str):
            raise ValueError("text column requires text values")
        return value
    if kind == "boolean":
        if not isinstance(value, (bool, np.bool_)):
            raise ValueError("boolean column requires boolean values")
        return "True" if value else "False"
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError("numeric column requires numeric values, not booleans or text")
    if kind == "integer":
        if isinstance(value, Integral):
            return str(int(value))
        number = float(value)
        if not isfinite(number) or not number.is_integer() or abs(number) >= 2**53:
            raise ValueError(
                "integer column requires exact integers; large values need integer storage"
            )
        return str(int(number))
    number = float(value)
    if not isfinite(number):
        raise ValueError("numeric table values must be finite or absent")
    if isinstance(value, Integral) and int(number) != value:
        raise ValueError(
            "number column would lose integer precision; use an integer column"
        )
    return format(number, ".17g")


def write_scalar_table(frame: pd.DataFrame, path: Path, kinds: dict[str, str]) -> dict:
    names = list(frame.columns)
    if len(set(names)) != len(names) or set(names) != set(kinds):
        raise ValueError("product table columns must match the declared kinds exactly")
    ordered = {name: kinds[name] for name in names}
    _validate_columns(ordered)
    stream = StringIO(newline="")
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(names)
    for row_index, row in enumerate(frame.itertuples(index=False, name=None), start=1):
        encoded = []
        for name, value in zip(names, row, strict=True):
            try:
                encoded.append(_encode(value, ordered[name]))
            except ValueError as error:
                raise ValueError(
                    f"column {name!r}, row {row_index}: {error}"
                ) from error
        writer.writerow(encoded)
    # Conversion errors leave any existing file intact. Product publication owns
    # staging, locking, and atomic replacement of the complete generation.
    path.write_bytes(stream.getvalue().encode("utf-8"))
    return {
        "columns": [{"name": name, "kind": kind} for name, kind in ordered.items()],
        "rows": len(frame),
    }


def read_scalar_table(content: bytes, schema: dict) -> pd.DataFrame:
    descriptors = schema.get("columns")
    count = schema.get("rows")
    if (
        not isinstance(descriptors, list)
        or isinstance(count, bool)
        or not isinstance(count, int)
        or count < 0
    ):
        raise ValueError("invalid product table schema")
    try:
        names = [column["name"] for column in descriptors]
        kinds = {column["name"]: column["kind"] for column in descriptors}
    except (KeyError, TypeError) as error:
        raise ValueError("invalid product table schema") from error
    if len(names) != len(kinds):
        raise ValueError("duplicate product table columns")
    _validate_columns(kinds)
    try:
        records = list(
            csv.reader(StringIO(content.decode("utf-8"), newline=""), strict=True)
        )
    except csv.Error as error:
        raise ValueError("invalid product CSV") from error
    if not records or records[0] != names or len(records) - 1 != count:
        raise ValueError(
            "product table does not match its declared columns or row count"
        )
    if any(len(row) != len(names) for row in records[1:]):
        raise ValueError("product CSV row width does not match its declared columns")
    values = {}
    for index, (name, kind) in enumerate(kinds.items()):
        column = [row[index] for row in records[1:]]
        if kind == "string":
            values[name] = pd.Series(column, dtype=object)
        elif kind == "boolean":
            if any(value not in {"", "True", "False"} for value in column):
                raise ValueError("invalid boolean in product table")
            values[name] = pd.Series(
                [pd.NA if value == "" else value == "True" for value in column],
                dtype="boolean",
            )
        elif kind == "integer":
            values[name] = pd.Series(
                [None if value == "" else int(value) for value in column], dtype=object
            )
        else:
            parsed = [float("nan") if value == "" else float(value) for value in column]
            if any(
                value != "" and not isfinite(number)
                for value, number in zip(column, parsed, strict=True)
            ):
                raise ValueError("nonfinite number in product table")
            values[name] = pd.Series(parsed, dtype=float)
    return pd.DataFrame(values)
