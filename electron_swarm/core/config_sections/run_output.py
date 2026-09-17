"""Run selection, comparison policy, and output parsers."""

from __future__ import annotations

from pathlib import Path
from typing import Any, cast

from electron_swarm.core.config import (
    ComparisonConfig,
    DegradedPolicy,
    FeaturePolicyConfig,
    OutputConfig,
    RequestedSolverConfig,
    RunConfig,
    UnsupportedPolicy,
)
from electron_swarm.core.config_sections.common import (
    as_mapping_section,
    as_path,
    bool_field,
    canonical_solver,
    list_float,
    reject_unknown_fields,
    string_value,
    validate_literal,
)
from electron_swarm.core.legacy_schema import MIGRATION_ERROR


def parse_run(raw: dict[str, Any]) -> RunConfig:
    run_raw = as_mapping_section(raw, "run", "run")
    reject_unknown_fields(run_raw, {"solvers", "e_over_n_Td", "case_prefix"}, "run")
    solvers_raw = run_raw.get("solvers", [])
    if not isinstance(solvers_raw, list):
        raise ValueError("run.solvers must be a list")

    solvers: list[RequestedSolverConfig] = []
    for index, item in enumerate(solvers_raw):
        if not isinstance(item, dict) or "id" not in item:
            raise ValueError("run.solvers entries must be mappings with id")
        reject_unknown_fields(item, {"id", "enabled"}, f"run.solvers[{index}]")
        solvers.append(
            RequestedSolverConfig(
                id=canonical_solver(item["id"], f"run.solvers[{index}].id"),
                enabled=bool_field(
                    item,
                    "enabled",
                    True,
                    f"run.solvers[{index}].enabled",
                ),
            )
        )

    return RunConfig(
        solvers=solvers,
        e_over_n_Td=list_float(run_raw.get("e_over_n_Td", [100.0])),
        case_prefix=string_value(
            run_raw.get("case_prefix", "case"),
            "run.case_prefix",
        ),
    )


def parse_comparison(raw: dict[str, Any]) -> ComparisonConfig:
    cmp_raw = as_mapping_section(raw, "comparison", "comparison")
    reject_unknown_fields(
        cmp_raw,
        {
            "enabled",
            "reference_solver",
            "candidate_solvers",
            "compare_eedf",
            "required",
        },
        "comparison",
    )
    reference = cmp_raw.get("reference_solver")
    candidates = cmp_raw.get("candidate_solvers", [])
    if candidates is None:
        candidates = []
    if not isinstance(candidates, list):
        raise ValueError("comparison.candidate_solvers must be a list")
    return ComparisonConfig(
        enabled=bool_field(cmp_raw, "enabled", False, "comparison.enabled"),
        reference_solver=(
            canonical_solver(reference, "comparison.reference_solver")
            if reference is not None
            else None
        ),
        candidate_solvers=[
            canonical_solver(item, "comparison.candidate_solvers")
            for item in candidates
        ],
        compare_eedf=bool_field(
            cmp_raw, "compare_eedf", True, "comparison.compare_eedf"
        ),
        required=bool_field(cmp_raw, "required", False, "comparison.required"),
    )


def parse_feature_policy(raw: dict[str, Any]) -> FeaturePolicyConfig:
    policy_raw = as_mapping_section(raw, "feature_policy", "feature_policy")
    if "allow_unsupported_fallback" in policy_raw:
        raise ValueError(
            f"{MIGRATION_ERROR}; remove feature_policy.allow_unsupported_fallback"
        )
    reject_unknown_fields(policy_raw, {"unsupported", "degraded"}, "feature_policy")
    return FeaturePolicyConfig(
        unsupported=cast(
            UnsupportedPolicy,
            validate_literal(
                string_value(
                    policy_raw.get("unsupported", "fail"),
                    "feature_policy.unsupported",
                ),
                {"fail", "skip_solver"},
                "feature_policy.unsupported",
            ),
        ),
        degraded=cast(
            DegradedPolicy,
            validate_literal(
                string_value(
                    policy_raw.get("degraded", "record"),
                    "feature_policy.degraded",
                ),
                {"fail", "record"},
                "feature_policy.degraded",
            ),
        ),
    )


def parse_output(raw: dict[str, Any], base: Path) -> OutputConfig:
    out_raw = as_mapping_section(raw, "output", "output")
    return OutputConfig(
        directory=as_path(out_raw.get("directory", "outputs"), base)
        or (base / "outputs"),
        base_name=string_value(out_raw.get("base_name", "swarm"), "output.base_name"),
        float_format=string_value(
            out_raw.get("float_format", "%.10e"),
            "output.float_format",
        ),
    )
