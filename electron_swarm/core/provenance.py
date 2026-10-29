"""Separate executed numerical results from caller-provided or derived data."""

from dataclasses import dataclass, replace
import json
from pathlib import Path

import numpy as np

from .identity import content_digest


def _metadata_json(metadata) -> str:
    def encode(value):
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(f"unsupported execution metadata: {type(value).__name__}")

    return json.dumps(
        {k: v for k, v in metadata.items() if k != "output_paths"},
        sort_keys=True,
        default=encode,
        allow_nan=False,
    )


@dataclass(frozen=True, slots=True)
class ExecutionRecord:
    metadata_json: str
    results_sha256: str

    @property
    def identity(self) -> str:
        return content_digest((self.metadata_json, self.results_sha256))

    def matches(self, result) -> bool:
        return self.metadata_json == _metadata_json(
            result.metadata
        ) and self.results_sha256 == content_digest(result.cases)


def seal_execution(result) -> None:
    result.execution_record = ExecutionRecord(
        _metadata_json(result.metadata), content_digest(result.cases)
    )


def publication_view(result):
    """Return truthful publication context without mutating the caller's result."""
    record = result.execution_record
    if record is not None and record.matches(result):
        return (
            result,
            "solver_execution",
            {"id": record.identity, "results_sha256": record.results_sha256},
        )
    scope = "derived_results" if record is not None else "provided_results"
    metadata = {
        "comparison": {
            "status": "not_evaluated",
            "reason": "results_modified" if record else "execution_not_provided",
        }
    }
    if record is not None:
        metadata["source_execution"] = {
            "id": record.identity,
            "results_sha256": record.results_sha256,
            "metadata": json.loads(record.metadata_json),
        }
    metadata["outcome"] = {
        "status": "derived" if record else "provided",
        "case_count": len(result.cases),
        "comparison_status": metadata.get("comparison", {}).get(
            "status", "not_evaluated"
        ),
    }
    cases = [
        replace(
            case,
            metadata={
                **case.metadata,
                "numerical_status": metadata["outcome"]["status"],
            },
        )
        for case in result.cases
    ]
    return replace(result, cases=cases, metadata=metadata), scope, None
