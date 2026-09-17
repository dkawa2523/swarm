"""Run independent propagator validation cases and retain failed evidence.

The campaign runner groups deterministic anchors for warm starts.  A failed
anchor therefore prevents later results in the same group from being written,
and the current propagator exception is not process-pickleable.  This external
validation utility makes every mixture/E/N point an independent bounded job,
catches the original exception in its worker, and writes one auditable JSON
matrix.  It does not alter solver behavior.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import dataclasses
import hashlib
import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from electron_swarm.core.config_parser import load_config  # noqa: E402
from electron_swarm.runner import run  # noqa: E402
from electron_swarm.solvers.propagator.steady import (  # noqa: E402
    PropagatorConvergenceError,
)
from swarm_workflow.campaign.config import MixtureSpec, load_workflow  # noqa: E402
from swarm_workflow.campaign.sweep import _clone_config  # noqa: E402


SCHEMA = "swarm.validation.propagator_case_matrix.v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_value(value: Any) -> Any:
    if dataclasses.is_dataclass(value):
        return _json_value(dataclasses.asdict(value))
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, np.ndarray):
        return [_json_value(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


def _run_one(
    base_config_path: str,
    mixture_id: int,
    fractions: dict[str, float],
    field_td: float,
) -> dict[str, Any]:
    started = time.perf_counter()
    base = load_config(base_config_path)
    mixture = MixtureSpec(mixture_id=mixture_id, fractions=fractions)
    config = _clone_config(
        base,
        mixture=mixture,
        e_over_n_values=(field_td,),
        solver_ids=["propagator"],
        replicate=0,
    )
    try:
        result = run(config, write=False)
        case = result.cases[0]
        return {
            "mixture_id": mixture_id,
            "fractions": fractions,
            "E_over_N_Td": field_td,
            "status": "passed",
            "elapsed_s": time.perf_counter() - started,
            "mean_energy_eV": case.mean_energy_eV,
            "drift_velocity_m_s": case.drift_velocity_m_s,
            "net_ionization_frequency_s_inv": case.net_ionization_frequency_s,
            "energy_eV": case.energy_eV,
            "energy_widths_eV": case.energy_widths_eV,
            "eedf": case.eedf,
            "diagnostics": case.diagnostics,
        }
    except PropagatorConvergenceError as exc:
        return {
            "mixture_id": mixture_id,
            "fractions": fractions,
            "E_over_N_Td": field_td,
            "status": "failed",
            "elapsed_s": time.perf_counter() - started,
            "error_type": type(exc).__name__,
            "error": str(exc),
            "diagnostics": exc.diagnostics,
        }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workflow", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    workflow_path = args.workflow.resolve()
    output = args.output.resolve()
    if output.exists():
        parser.error(f"refusing to overwrite {output}")
    workflow = load_workflow(workflow_path)
    jobs = [
        (
            str(workflow.base_config_path),
            mixture.mixture_id,
            dict(mixture.fractions),
            float(field),
        )
        for mixture in workflow.mixtures
        for field in workflow.e_over_n_Td
    ]
    results: list[dict[str, Any]] = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(_run_one, *job) for job in jobs]
        for future in concurrent.futures.as_completed(futures):
            item = future.result()
            results.append(item)
            print(
                f"mixture {item['mixture_id']} / {item['E_over_N_Td']:g} Td: "
                f"{item['status']} ({item['elapsed_s']:.2f}s)",
                flush=True,
            )
    results.sort(key=lambda item: (item["mixture_id"], item["E_over_N_Td"]))
    payload = {
        "schema": SCHEMA,
        "workflow": str(workflow_path),
        "workflow_sha256": _sha256(workflow_path),
        "base_config": str(workflow.base_config_path),
        "base_config_sha256": _sha256(workflow.base_config_path),
        "workers": args.workers,
        "passed": sum(item["status"] == "passed" for item in results),
        "failed": sum(item["status"] == "failed" for item in results),
        "cases": results,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(_json_value(payload), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(output)


if __name__ == "__main__":
    main()
