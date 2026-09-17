"""Independent Function-EEDF adapter for the Argon GEC-ICP model."""

from __future__ import annotations

import json

from .bundle import (
    GecIcpBundleError,
    validate_gec_icp_bundle,
)
from .comparison_analysis import (
    GecIcpComparisonAnalysisError,
    GecIcpComparisonAnalysisSummary,
    generate_gec_icp_saved_solution_comparison,
)
from .comparison_export import (
    GecIcpComparisonExportError,
    GecIcpComparisonExportPlan,
    execute_gec_icp_comparison_export,
    prepare_gec_icp_comparison_export,
)
from .execution import execute_gec_icp_run
from .prepare import prepare_gec_icp_run
from .run_contracts import (
    GecIcpPlan,
    GecIcpQualityError,
    GecIcpRunSummary,
    GecIcpWorkflowError,
)


def format_gec_icp_plan(plan: GecIcpPlan) -> str:
    return "\n".join(
        [
            "GEC ICP COMSOL dry-run ready",
            f"source: {plan.mapping.bundle.expected_source}",
            f"model: {plan.mapping.model_mapping.model.input_mph}",
            f"bundle: {plan.mapping.bundle.path}",
            "closure: Function EEDF + reduced mobility; COMSOL owns "
            "diffusion, energy transport, RF, and chemistry",
            f"final physical time: {plan.mapping.run.final_time_s:.6g} s",
            f"plan: {plan.plan_json}",
            "COMSOL command: not executed",
        ]
    )


def format_gec_icp_summary(summary: GecIcpRunSummary) -> str:
    status = json.loads(summary.status_json.read_text(encoding="utf-8"))
    convergence = status.get("convergence") or {}
    return "\n".join(
        [
            "GEC ICP COMSOL workflow completed",
            f"source: {status.get('source')}",
            f"physical target accepted: {status.get('physical_target_accepted')}",
            f"convergence: {convergence.get('status')}",
            f"status: {summary.status_json}",
            f"results: {summary.plan.output_directory}",
        ]
    )


__all__ = (
    "GecIcpBundleError",
    "GecIcpComparisonAnalysisError",
    "GecIcpComparisonAnalysisSummary",
    "GecIcpComparisonExportError",
    "GecIcpComparisonExportPlan",
    "GecIcpPlan",
    "GecIcpQualityError",
    "GecIcpRunSummary",
    "GecIcpWorkflowError",
    "execute_gec_icp_comparison_export",
    "execute_gec_icp_run",
    "format_gec_icp_plan",
    "format_gec_icp_summary",
    "generate_gec_icp_saved_solution_comparison",
    "prepare_gec_icp_comparison_export",
    "prepare_gec_icp_run",
    "validate_gec_icp_bundle",
)
