"""Validated scalar lookup data with an explicit collision-channel identity."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import math
from pathlib import Path


class LookupValidationError(ValueError):
    """A table cannot define one unambiguous interpolation function."""


@dataclass(frozen=True, slots=True)
class ReactionIdentity:
    species: str
    process: str
    process_type: str
    threshold_eV: float


@dataclass(frozen=True, slots=True)
class ValidatedLookup:
    path: Path
    x_column: str
    y_column: str
    x: tuple[float, ...]
    y: tuple[float, ...]
    reaction: ReactionIdentity | None = None


def read_lookup(
    path: Path,
    x_column: str,
    y_column: str,
    *,
    reaction: ReactionIdentity | None = None,
) -> ValidatedLookup:
    """Select one channel and reject duplicate coordinates instead of overwriting."""

    try:
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            required = {x_column, y_column}
            if reaction is not None:
                required.update({"species", "process", "process_type", "threshold_eV"})
            missing = required - set(reader.fieldnames or ())
            if missing:
                raise LookupValidationError(
                    f"lookup {path.name} is missing columns: {', '.join(sorted(missing))}"
                )
            values: list[tuple[float, float]] = []
            for index, row in enumerate(reader):
                if reaction is not None:
                    if any(
                        row[name] != getattr(reaction, name)
                        for name in ("species", "process", "process_type")
                    ):
                        continue
                    threshold = float(row["threshold_eV"])
                    if not math.isfinite(threshold):
                        raise LookupValidationError(
                            f"nonfinite reaction threshold in {path}"
                        )
                    if not math.isclose(
                        threshold, reaction.threshold_eV, rel_tol=1e-12, abs_tol=1e-12
                    ):
                        continue
                x, y = float(row[x_column]), float(row[y_column])
                if not math.isfinite(x) or not math.isfinite(y) or x < 0 or y < 0:
                    raise LookupValidationError(
                        f"invalid lookup value in {path.name} row {index}"
                    )
                values.append((x, y))
    except (OSError, KeyError, TypeError, ValueError) as exc:
        if isinstance(exc, LookupValidationError):
            raise
        raise LookupValidationError(f"invalid lookup CSV: {path}") from exc
    values.sort()
    if len(values) < 2:
        raise LookupValidationError(
            f"lookup {path.name} needs at least two rows for reaction {reaction!r}"
        )
    if any(b[0] <= a[0] for a, b in zip(values, values[1:])):
        raise LookupValidationError(
            f"lookup {path.name} has duplicate coordinates for reaction {reaction!r}; "
            "ambiguous channels cannot be combined"
        )
    return ValidatedLookup(
        path,
        x_column,
        y_column,
        tuple(x for x, _ in values),
        tuple(y for _, y in values),
        reaction,
    )
