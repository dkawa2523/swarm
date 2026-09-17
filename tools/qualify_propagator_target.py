"""Qualify a bounded Propagator target range by medium/fine refinement.

The medium cases come from a provenance-bound workflow database.  The fine
cases are independent direct solves at explicitly selected target anchors.
This tool qualifies only the numerical P1 target range; it does not execute or
approve any downstream plasma model.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from hashlib import sha256
import json
import math
import multiprocessing
import os
from pathlib import Path
import re
import sqlite3
import sys
from types import SimpleNamespace
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from electron_swarm import load_config  # noqa: E402
from electron_swarm.core.config import SwarmConfig  # noqa: E402
from swarm_workflow.quality.solver import (  # noqa: E402
    PROPAGATOR_TARGET_QUALIFICATION_SCHEMA,
    file_sha256,
    propagator_target_medium_evidence,
    propagator_target_qualifier_fingerprint,
    validate_propagator_core_qualification,
)
from swarm_workflow.quality.propagator_source import (  # noqa: E402
    PROPAGATOR_SOLVER_SOURCE_METADATA_KEY,
    propagator_qualification_source_fingerprint,
)
from tools.qualify_propagator_p1_deterministic import (  # noqa: E402
    EEDF_L1_REFINEMENT_LIMIT,
    SCALAR_REFINEMENT_LIMIT,
    _comparison,
    _run_case,
    _variant,
)


SCHEMA = PROPAGATOR_TARGET_QUALIFICATION_SCHEMA
DEFAULT_CORE_QUALIFICATION = (
    ROOT
    / "docs"
    / "dev"
    / "results"
    / "propagator_p1_deterministic_qualification_20260908.json"
)
_TARGET_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")


def _cross_sections_sha256(config: SwarmConfig) -> str:
    file_hashes = {
        str(file_config.path): file_sha256(file_config.path)
        for file_config in config.cross_sections.files
    }
    digest = sha256()
    for path, value in sorted(file_hashes.items()):
        digest.update(path.encode("utf-8") + b"\0" + value.encode("ascii") + b"\0")
    return digest.hexdigest()


def _require_unchanged_target_inputs(
    starting_state: dict[str, Any],
    *,
    spec: TargetRefinementSpec,
    database: Path,
    config_path: Path,
    core_qualification_path: Path,
) -> None:
    try:
        core = validate_propagator_core_qualification(core_qualification_path)
        config_sha256 = file_sha256(config_path)
        cross_sections_sha256 = _cross_sections_sha256(load_config(config_path))
        medium_digest, metadata = propagator_target_medium_evidence(
            database,
            spec.fields_Td,
            expected_grid=spec.medium_grid,
            mixture_id=spec.mixture_id,
        )
        current_state = {
            "core_qualification": core.identity_entry(),
            "base_config_sha256": config_sha256,
            "cross_sections_sha256": cross_sections_sha256,
            "medium_evidence_sha256": medium_digest,
            "medium_metadata": metadata,
            "target_qualifier_fingerprint": (
                propagator_target_qualifier_fingerprint(ROOT)
            ),
        }
    except (OSError, TypeError, ValueError, sqlite3.Error) as exc:
        raise RuntimeError(
            "Propagator target qualification inputs changed during execution; "
            "discard the results and rerun"
        ) from exc
    if current_state != starting_state:
        raise RuntimeError(
            "Propagator target qualification inputs changed during execution; "
            "discard the results and rerun"
        )


@dataclass(frozen=True, slots=True)
class TargetRefinementSpec:
    """Explicit numerical contract for one Propagator target range."""

    target: str
    fields_Td: tuple[float, ...]
    operating_range_bracket_Td: tuple[float, float]
    medium_grid: tuple[int, int]
    fine_grid: tuple[int, int]
    table_support_cap_Td: float | None = None
    mixture_id: int = 0
    fine_max_memory_mb: int = 1024
    scalar_refinement_limit: float = SCALAR_REFINEMENT_LIMIT
    eedf_weighted_L1_limit: float = EEDF_L1_REFINEMENT_LIMIT
    family: str | None = None
    progress_label: str | None = None
    release_claim: str = (
        "target_specific_P1_numerical_scope_only; downstream closure remains "
        "a separate decision"
    )

    def __post_init__(self) -> None:
        fields = tuple(float(value) for value in self.fields_Td)
        bracket = tuple(float(value) for value in self.operating_range_bracket_Td)
        medium_grid = _validated_grid(self.medium_grid, "medium grid")
        fine_grid = _validated_grid(self.fine_grid, "fine grid")
        if self.table_support_cap_Td is None:
            cap = max(fields) if fields else None
        else:
            cap = float(self.table_support_cap_Td)
        object.__setattr__(self, "fields_Td", fields)
        object.__setattr__(self, "operating_range_bracket_Td", bracket)
        object.__setattr__(self, "medium_grid", medium_grid)
        object.__setattr__(self, "fine_grid", fine_grid)
        object.__setattr__(self, "table_support_cap_Td", cap)

        if not isinstance(self.target, str) or not _TARGET_PATTERN.fullmatch(
            self.target
        ):
            raise ValueError("target must be a non-empty machine-readable name")
        if not fields or any(
            not math.isfinite(value) or value <= 0 for value in fields
        ):
            raise ValueError("target fields must be finite and positive")
        if any(right <= left for left, right in zip(fields, fields[1:])):
            raise ValueError("target fields must be unique and strictly increasing")
        if len(bracket) != 2 or any(
            not math.isfinite(value) or value <= 0 for value in bracket
        ):
            raise ValueError("operating bracket must contain two positive fields")
        if bracket[0] > bracket[1]:
            raise ValueError("operating bracket must be ordered")
        if any(value not in fields for value in bracket):
            raise ValueError("operating bracket anchors must be target fields")
        if cap is None or not math.isfinite(cap) or cap <= 0:
            raise ValueError("table support cap must be finite and positive")
        if cap not in fields or cap < bracket[1]:
            raise ValueError(
                "table support cap must be a target field at or above the bracket"
            )
        if fine_grid[0] < medium_grid[0] or fine_grid[1] < medium_grid[1]:
            raise ValueError("fine grid cannot be coarser than medium grid")
        if fine_grid == medium_grid:
            raise ValueError("fine grid must refine at least one dimension")
        if (
            isinstance(self.mixture_id, bool)
            or not isinstance(self.mixture_id, int)
            or self.mixture_id < 0
        ):
            raise ValueError("mixture id must be a non-negative integer")
        if (
            isinstance(self.fine_max_memory_mb, bool)
            or not isinstance(self.fine_max_memory_mb, int)
            or self.fine_max_memory_mb < 128
        ):
            raise ValueError("fine-grid memory bound must be at least 128 MB")
        for value, label in (
            (self.scalar_refinement_limit, "scalar refinement limit"),
            (self.eedf_weighted_L1_limit, "EEDF refinement limit"),
        ):
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{label} must be finite and positive")
        if not self.family:
            object.__setattr__(self, "family", self.target)
        if not self.progress_label:
            object.__setattr__(self, "progress_label", self.target)


def _validated_grid(value: tuple[int, int], label: str) -> tuple[int, int]:
    if (
        len(value) != 2
        or any(isinstance(item, bool) or not isinstance(item, int) for item in value)
        or value[0] < 2
        or value[1] < 2
        or value[1] % 2
    ):
        raise ValueError(f"{label} must contain positive cells and even polar cells")
    return int(value[0]), int(value[1])


def _read_medium_cases(
    database: Path,
    spec: TargetRefinementSpec,
) -> tuple[dict[float, Any], dict[str, Any], str]:
    starting_digest, starting_metadata = propagator_target_medium_evidence(
        database,
        spec.fields_Td,
        expected_grid=spec.medium_grid,
        mixture_id=spec.mixture_id,
    )
    cases: dict[float, Any] = {}
    with sqlite3.connect(database) as connection:
        connection.row_factory = sqlite3.Row
        metadata = {
            str(row[0]): str(row[1])
            for row in connection.execute(
                "SELECT key, value FROM metadata ORDER BY key"
            ).fetchall()
        }
        for field in spec.fields_Td:
            row = connection.execute(
                "SELECT mean_energy_eV, drift_velocity_m_s, diagnostics_json "
                "FROM cases WHERE mixture_id=? AND solver='propagator' "
                "AND e_over_n_Td=? AND replicate=0",
                (spec.mixture_id, field),
            ).fetchone()
            if row is None:
                raise ValueError(f"workflow database lacks Propagator {field:g} Td")
            quality = connection.execute(
                "SELECT passed FROM aggregate_quality WHERE mixture_id=? "
                "AND solver='propagator' AND e_over_n_Td=?",
                (spec.mixture_id, field),
            ).fetchone()
            if quality is None or int(quality[0]) != 1:
                raise ValueError(f"workflow database has unqualified {field:g} Td")
            bins = connection.execute(
                "SELECT energy_eV, energy_width_eV, eedf FROM eedf_bins "
                "WHERE mixture_id=? AND solver='propagator' "
                "AND e_over_n_Td=? AND replicate=0 ORDER BY bin_index",
                (spec.mixture_id, field),
            ).fetchall()
            if not bins:
                raise ValueError(f"workflow database lacks EEDF bins at {field:g} Td")
            diagnostics = json.loads(str(row["diagnostics_json"]))
            cases[field] = SimpleNamespace(
                e_over_n_Td=field,
                mean_energy_eV=float(row["mean_energy_eV"]),
                drift_velocity_m_s=float(row["drift_velocity_m_s"]),
                diagnostics=diagnostics,
                energy_eV=np.asarray([item[0] for item in bins], dtype=float),
                energy_widths_eV=np.asarray([item[1] for item in bins], dtype=float),
                eedf=np.asarray([item[2] for item in bins], dtype=float),
            )
    medium_digest, validated_metadata = propagator_target_medium_evidence(
        database,
        spec.fields_Td,
        expected_grid=spec.medium_grid,
        mixture_id=spec.mixture_id,
    )
    if (
        metadata != starting_metadata
        or metadata != validated_metadata
        or medium_digest != starting_digest
    ):
        raise ValueError("workflow evidence changed while target cases were read")
    return cases, metadata, medium_digest


def _fine_worker(
    config_path: str,
    field: float,
    fine_grid: tuple[int, int],
    fine_max_memory_mb: int,
    family: str,
    expected_config_sha256: str,
    expected_cross_sections_sha256: str,
    expected_core_source_sha256: str,
    expected_target_qualifier_sha256: str,
) -> tuple[dict[str, Any], Any | None]:
    if file_sha256(config_path) != expected_config_sha256:
        raise RuntimeError("target qualification config changed before worker start")
    base = load_config(config_path)
    if _cross_sections_sha256(base) != expected_cross_sections_sha256:
        raise RuntimeError(
            "target qualification cross sections changed before worker start"
        )
    if (
        propagator_qualification_source_fingerprint(ROOT)["sha256"]
        != expected_core_source_sha256
    ):
        raise RuntimeError(
            "Propagator implementation changed before target worker start"
        )
    if (
        propagator_target_qualifier_fingerprint(ROOT)["sha256"]
        != expected_target_qualifier_sha256
    ):
        raise RuntimeError(
            "Propagator target qualification algorithm changed before worker start"
        )
    base.solvers.propagator.max_memory_mb = fine_max_memory_mb
    return _run_case(
        _variant(base, field, fine_grid),
        family=family,
        profile="fine",
    )


def run_target_qualification(
    *,
    spec: TargetRefinementSpec,
    database: Path,
    config_path: Path,
    core_qualification_path: Path,
    workers: int,
    memory_budget_mb: int,
) -> dict[str, Any]:
    """Run the independent fine cases and compare them with database cases."""

    if workers < 1:
        raise ValueError("workers must be positive")
    if memory_budget_mb < spec.fine_max_memory_mb:
        raise ValueError("memory budget must cover at least one fine-grid worker bound")
    target_fingerprint = propagator_target_qualifier_fingerprint(ROOT)
    core = validate_propagator_core_qualification(core_qualification_path)
    core_entry = core.identity_entry()
    config_sha256 = file_sha256(config_path)
    config = load_config(config_path)
    cross_sections_sha256 = _cross_sections_sha256(config)
    medium, metadata, medium_digest = _read_medium_cases(database, spec)
    if metadata.get("base_config_sha256") != config_sha256:
        raise ValueError(
            "workflow database base config does not match the target config"
        )
    if metadata.get("cross_sections_sha256") != cross_sections_sha256:
        raise ValueError(
            "workflow database cross sections do not match the target config"
        )
    core_source_sha256 = str(core_entry["implementation_fingerprint"]["sha256"])
    if metadata.get(PROPAGATOR_SOLVER_SOURCE_METADATA_KEY) != core_source_sha256:
        raise ValueError(
            "workflow database Propagator source does not match the qualified core"
        )
    starting_state = {
        "core_qualification": core_entry,
        "base_config_sha256": config_sha256,
        "cross_sections_sha256": cross_sections_sha256,
        "medium_evidence_sha256": medium_digest,
        "medium_metadata": dict(metadata),
        "target_qualifier_fingerprint": target_fingerprint,
    }
    memory_workers = max(1, memory_budget_mb // spec.fine_max_memory_mb)
    effective_workers = max(
        1,
        min(workers, len(spec.fields_Td), os.cpu_count() or 1, memory_workers),
    )

    fine_cases: dict[float, Any | None] = {}
    fine_records: list[dict[str, Any]] = []
    context = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=effective_workers,
        mp_context=context,
    ) as executor:
        futures = {
            executor.submit(
                _fine_worker,
                str(config_path),
                field,
                spec.fine_grid,
                spec.fine_max_memory_mb,
                str(spec.family),
                config_sha256,
                cross_sections_sha256,
                core_source_sha256,
                str(target_fingerprint["sha256"]),
            ): field
            for field in spec.fields_Td
        }
        for future in as_completed(futures):
            field = futures[future]
            record, case = future.result()
            fine_records.append(record)
            fine_cases[field] = case
            print(
                f"{spec.progress_label}/fine/{field:g}Td: {record['status']} "
                f"({record['elapsed_s']:.2f}s)",
                flush=True,
            )
    fine_records.sort(key=lambda item: float(item["E_over_N_Td"]))

    refinement = [
        {
            "E_over_N_Td": field,
            **_comparison(
                medium[field],
                fine_cases.get(field),
                scalar_limit=spec.scalar_refinement_limit,
                eedf_limit=spec.eedf_weighted_L1_limit,
                gate_growth=True,
            ),
        }
        for field in spec.fields_Td
    ]
    run_gate = bool(fine_records) and all(
        row["quality_gates_passed"] for row in fine_records
    )
    refinement_gate = all(row["passed"] for row in refinement)
    blockers = []
    if not run_gate:
        blockers.append("fine_case_quality")
    if not refinement_gate:
        blockers.append("medium_fine_refinement")
    passed = not blockers
    scope: dict[str, Any] = {
        "target": spec.target,
        "included_solver": "propagator",
        "field_model": "homogeneous_dc_B0",
        "fields_Td": list(spec.fields_Td),
        "operating_range_bracket_Td": list(spec.operating_range_bracket_Td),
        "table_support_cap_Td": spec.table_support_cap_Td,
        "medium_grid": list(spec.medium_grid),
        "fine_grid": list(spec.fine_grid),
        "scalar_refinement_limit": spec.scalar_refinement_limit,
        "eedf_weighted_L1_limit": spec.eedf_weighted_L1_limit,
        "mixture_id": spec.mixture_id,
        "excluded": [
            "two_term",
            "multi_term",
            "monte_carlo",
            "COMSOL",
            "P2_bulk_transport_and_diffusion",
        ],
    }
    _require_unchanged_target_inputs(
        starting_state,
        spec=spec,
        database=database,
        config_path=config_path,
        core_qualification_path=core_qualification_path,
    )
    return {
        "schema": SCHEMA,
        "scope": scope,
        "environment": {
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "logical_cpu_count": os.cpu_count(),
            "requested_workers": workers,
            "effective_workers": effective_workers,
            "memory_budget_mb": memory_budget_mb,
            "configured_memory_bound_per_worker_mb": spec.fine_max_memory_mb,
            "target_qualifier_fingerprint": target_fingerprint,
        },
        "inputs": {
            "database": str(database),
            "medium_evidence_sha256": medium_digest,
            "base_config": str(config_path),
            "base_config_sha256": config_sha256,
            "workflow_config_sha256": metadata.get("workflow_config_sha256"),
            "cross_sections_sha256": cross_sections_sha256,
            PROPAGATOR_SOLVER_SOURCE_METADATA_KEY: core_source_sha256,
            "core_qualification": core_entry,
        },
        "fine_runs": fine_records,
        "medium_fine_refinement": refinement,
        "decision": {
            "fine_case_quality_passed": run_gate,
            "medium_fine_refinement_passed": refinement_gate,
            "target_refinement_qualified": passed,
            "blocking_gates": blockers,
            "release_claim": spec.release_claim,
        },
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Qualify a target-specific Propagator medium/fine refinement."
    )
    parser.add_argument("--target", required=True)
    parser.add_argument("--fields-Td", type=float, nargs="+", required=True)
    parser.add_argument("--operating-bracket-Td", type=float, nargs=2, required=True)
    parser.add_argument("--table-support-cap-Td", type=float)
    parser.add_argument("--medium-grid", type=int, nargs=2, default=(300, 48))
    parser.add_argument("--fine-grid", type=int, nargs=2, default=(600, 72))
    parser.add_argument("--mixture-id", type=int, default=0)
    parser.add_argument("--fine-max-memory-mb", type=int, default=1024)
    parser.add_argument("--family")
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument(
        "--core-qualification", type=Path, default=DEFAULT_CORE_QUALIFICATION
    )
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--memory-budget-mb", type=int, default=2048)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        spec = TargetRefinementSpec(
            target=args.target,
            fields_Td=tuple(args.fields_Td),
            operating_range_bracket_Td=tuple(args.operating_bracket_Td),
            table_support_cap_Td=args.table_support_cap_Td,
            medium_grid=tuple(args.medium_grid),
            fine_grid=tuple(args.fine_grid),
            mixture_id=args.mixture_id,
            fine_max_memory_mb=args.fine_max_memory_mb,
            family=args.family,
        )
        if args.workers < 1 or args.memory_budget_mb < spec.fine_max_memory_mb:
            raise ValueError(
                "workers must be positive and memory budget must cover one worker"
            )
    except (TypeError, ValueError) as exc:
        parser.error(str(exc))
    payload = run_target_qualification(
        spec=spec,
        database=args.database.resolve(),
        config_path=args.config.resolve(),
        core_qualification_path=args.core_qualification.resolve(),
        workers=args.workers,
        memory_budget_mb=args.memory_budget_mb,
    )
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(output)
    print(json.dumps(payload["decision"], indent=2, sort_keys=True))
    return 0 if payload["decision"]["target_refinement_qualified"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
