"""Thin orchestration for the typed GEC-CCP execution stages."""

from __future__ import annotations

from pathlib import Path

from ..contracts import GecCcpRunSummary
from .postsolve import evaluate_run_quality, run_postsolve_audits
from .preflight import prepare_execution
from .runtime import validate_execution_outputs
from .solver import run_solver_steps
from .status import enforce_quality_acceptance, write_final_status


def execute_gec_ccp_run(
    mapping_path: str | Path,
    *,
    bundle_path: str | Path | None = None,
    comsol_executable: str | Path | None = None,
) -> GecCcpRunSummary:
    """Execute one validated GEC plan and require physics-quality acceptance."""

    prepared = prepare_execution(mapping_path, bundle_path=bundle_path)
    solver = run_solver_steps(
        prepared,
        comsol_executable=comsol_executable,
    )
    comsol_runtime = validate_execution_outputs(prepared, solver)
    audits = run_postsolve_audits(prepared.plan)
    quality = evaluate_run_quality(prepared, solver, audits)
    write_final_status(prepared, solver, audits, quality, comsol_runtime)
    enforce_quality_acceptance(prepared, solver, audits, quality)
    return GecCcpRunSummary(
        prepared.plan,
        solver.executions,
        prepared.status_path,
    )


__all__ = ["execute_gec_ccp_run"]
