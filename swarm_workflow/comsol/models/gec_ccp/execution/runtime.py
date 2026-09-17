"""Validate completed COMSOL outputs and runtime/model identity."""

from __future__ import annotations

import json
import re
from typing import Any, Iterable

from swarm_workflow._io import write_json
from swarm_workflow.comsol.runtime import ComsolExecutionSummary

from ..contracts import GecCcpPlan, GecCcpWorkflowError
from ..mph import inspect_gec_ccp_mph
from .context import PreparedExecution, SolverExecution
from .status import status_contract_metadata


def validate_execution_outputs(
    prepared: PreparedExecution,
    solver: SolverExecution,
) -> dict[str, Any]:
    """Require all result files and one consistent COMSOL runtime identity."""

    plan = prepared.plan
    missing = [
        str(path) for path in plan.expected_result_files if not path.exists()
    ]
    if missing:
        write_json(
            prepared.status_path,
            {
                "stage": "run-gec-ccp",
                "status": "blocked",
                "blocker_code": "missing_result_files",
                "solve_status": "incomplete",
                "quality_status": "not_evaluated",
                "completed_operations": [
                    item.operation for item in solver.executions
                ],
                "missing_results": missing,
                "plan": str(plan.plan_json),
                "plan_sha256": prepared.plan_sha256,
                **status_contract_metadata(
                    plan.mapping,
                    gradient_response=prepared.gradient_response,
                ),
                **prepared.input_provenance,
            },
        )
        raise GecCcpWorkflowError(
            "COMSOL completed without all expected result files: "
            + ", ".join(missing)
        )
    try:
        comsol_runtime = _validated_comsol_runtime(
            solver.executions,
            expected_version=plan.contract.comsol_version,
        )
        saved_versions = _validated_saved_mph_versions(
            plan,
            runtime_version=str(comsol_runtime["version"]),
        )
        comsol_runtime["saved_mph_versions"] = saved_versions
        comsol_runtime["model_upgrade"] = {
            "performed": bool(
                comsol_runtime["input_model_upgrade_required"]
            ),
            "from": plan.contract.comsol_version,
            "to": saved_versions["external"],
        }
    except GecCcpWorkflowError as exc:
        write_json(
            prepared.status_path,
            {
                "stage": "run-gec-ccp",
                "status": "blocked",
                "blocker_code": "comsol_runtime_identity_mismatch",
                "reason": str(exc),
                "solve_status": "completed",
                "quality_status": "not_evaluated",
                "plan_sha256": prepared.plan_sha256,
                **status_contract_metadata(
                    plan.mapping,
                    gradient_response=prepared.gradient_response,
                ),
                **prepared.input_provenance,
            },
        )
        raise
    return comsol_runtime


def _comsol_major_minor(version: str) -> tuple[int, int]:
    match = re.match(r"^\s*(\d+)\.(\d+)", version)
    if match is None:
        raise GecCcpWorkflowError(f"unrecognized COMSOL version: {version}")
    return int(match.group(1)), int(match.group(2))


def _validated_comsol_runtime(
    executions: Iterable[ComsolExecutionSummary],
    *,
    expected_version: str | None,
) -> dict[str, Any]:
    execution_list = list(executions)
    versions = {item.comsol_version for item in execution_list}
    builds = {item.comsol_build for item in execution_list}
    executables: set[str] = set()
    for execution in execution_list:
        if execution.provenance_json is None:
            continue
        provenance = json.loads(
            execution.provenance_json.read_text(encoding="utf-8")
        )
        comsol = provenance.get("comsol")
        if isinstance(comsol, dict) and isinstance(comsol.get("executable"), str):
            executables.add(comsol["executable"])
    if (
        len(versions) != 1
        or None in versions
        or len(builds) != 1
        or None in builds
        or len(executables) != 1
    ):
        raise GecCcpWorkflowError(
            "COMSOL version/build/executable provenance is missing or mixed"
        )
    version = next(iter(versions))
    input_version = str(expected_version) if expected_version is not None else None
    upgraded = False
    if input_version is not None:
        input_major_minor = _comsol_major_minor(input_version)
        runtime_major_minor = _comsol_major_minor(str(version))
        if runtime_major_minor < input_major_minor:
            raise GecCcpWorkflowError(
                f"COMSOL runtime version {version} is older than MPH version "
                f"{input_version}"
            )
        upgraded = runtime_major_minor != input_major_minor
    return {
        "version": str(version),
        "build": str(next(iter(builds))),
        "executable": next(iter(executables)),
        "input_mph_version": input_version,
        "input_model_upgrade_required": upgraded,
    }


def _validated_saved_mph_versions(
    plan: GecCcpPlan,
    *,
    runtime_version: str,
) -> dict[str, str]:
    runtime_major_minor = _comsol_major_minor(runtime_version)
    versions: dict[str, str] = {}
    saved_models = {"external": plan.mapping.model.output_mph}
    if plan.mapping.run.include_builtin_reference:
        if plan.mapping.model.baseline_output_mph is None:
            raise GecCcpWorkflowError(
                "built-in reference output MPH is not configured"
            )
        saved_models["baseline"] = plan.mapping.model.baseline_output_mph
    for name, path in saved_models.items():
        version = inspect_gec_ccp_mph(path).comsol_version
        if version is None:
            raise GecCcpWorkflowError(
                f"saved {name} MPH has no readable COMSOL version"
            )
        if _comsol_major_minor(version) != runtime_major_minor:
            raise GecCcpWorkflowError(
                f"saved {name} MPH version {version} differs from runtime "
                f"version {runtime_version}"
            )
        versions[name] = version
    return versions


__all__ = ["validate_execution_outputs"]
