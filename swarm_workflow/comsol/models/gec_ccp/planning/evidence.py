"""Assemble the closure evidence required before a GEC-CCP run is planned."""

from __future__ import annotations

from dataclasses import dataclass
import json
from typing import Any

from swarm_workflow.comsol.input.bundle_selection import validate_bundle_selection
from swarm_workflow.selection import ClosureSelectionError

from ..closure import (
    _active_rate_support_metadata,
    _restricted_gradient_response_metadata,
    _uses_external_elastic_energy_loss,
    _uses_external_rates,
    _uses_transport,
)
from ..contracts import (
    FUNCTION_EEDF_PREINTEGRATED_INELASTIC,
    GEC_RESTRICTED_TRANSPORT_CLOSURE,
    GEC_TWO_TERM_HYBRID_TRANSPORT_CLOSURE,
    GecCcpMapping,
    GecCcpWorkflowError,
)
from ..data import _read_csv
from ..validation.bundle import validate_gec_ccp_bundle
from ..validation.bundle_guards import _validate_low_energy_guard_bundle
from ..validation.joint_consistency import (
    _active_closure_mean_energy_support,
    _build_dense_function_eedf_rate_closure,
    _build_two_term_joint_consistency_audit,
    _build_upstream_eedf_rate_consistency_audit,
)
from ..validation.mc_bundle import (
    _mc_rate_censoring_audit,
    _mc_rate_interpolation_audit,
)


@dataclass(frozen=True, slots=True)
class ClosureEvidence:
    bundle: dict[str, Any]
    low_energy_guard: dict[str, Any] | None
    restricted_gradient_response: dict[str, Any]
    preintegrated_rate_rows: list[dict[str, Any]] | None
    preintegrated_rate_summary: dict[str, Any] | None
    active_rate_support: dict[str, Any] | None
    mc_rate_censoring: dict[str, Any] | None
    mc_rate_interpolation: dict[str, Any] | None
    active_closure_support: dict[str, Any] | None
    upstream_eedf_rate_consistency: dict[str, Any] | None
    two_term_joint_consistency: dict[str, Any] | None


def build_closure_evidence(mapping: GecCcpMapping) -> ClosureEvidence:
    bundle = validate_gec_ccp_bundle(mapping)
    try:
        bundle["solver_selection"] = validate_bundle_selection(
            mapping.bundle.path,
            required=mapping.bundle.require_solver_selection,
        )
    except ClosureSelectionError as exc:
        raise GecCcpWorkflowError(str(exc)) from exc
    low_energy_guard = _validate_low_energy_guard_bundle(
        mapping,
        primary_cross_sections_sha256=str(
            bundle["hashes"]["cross_sections_sha256"]
        ),
    )
    bundle["low_energy_guard"] = low_energy_guard
    gradient_response = _restricted_gradient_response_metadata(
        mapping,
        scalar_diffusion_audit=bundle.get("restricted_scalar_diffusion_audit"),
    )
    rate_rows, rate_summary = _preintegrated_rates(mapping)
    active_rate_rows = _active_rate_rows(mapping, rate_rows)
    mc_evidence_rows = active_rate_rows
    if mapping.bundle.expected_source == "monte_carlo":
        mc_evidence_rows = _read_csv(
            mapping.bundle.path / "rates_vs_mean_energy.csv"
        )
    rate_support = (
        _active_rate_support_metadata(mapping, active_rate_rows)
        if active_rate_rows is not None
        else None
    )
    censoring = (
        _mc_rate_censoring_audit(
            mapping,
            mc_evidence_rows,
            required_rate_min_process_peak_fraction=float(
                bundle["quality_thresholds"][
                    "required_rate_min_process_peak_fraction"
                ]
            ),
        )
        if mc_evidence_rows is not None
        else None
    )
    _require_mc_censoring(censoring)
    interpolation = (
        _mc_rate_interpolation_audit(mapping, active_rate_rows)
        if active_rate_rows is not None
        else None
    )
    closure_support = _closure_support(mapping, active_rate_rows)
    upstream = _upstream_rate_consistency(mapping, bundle)
    joint = _two_term_joint_consistency(mapping, bundle, upstream)
    return ClosureEvidence(
        bundle=bundle,
        low_energy_guard=low_energy_guard,
        restricted_gradient_response=gradient_response,
        preintegrated_rate_rows=rate_rows,
        preintegrated_rate_summary=rate_summary,
        active_rate_support=rate_support,
        mc_rate_censoring=censoring,
        mc_rate_interpolation=interpolation,
        active_closure_support=closure_support,
        upstream_eedf_rate_consistency=upstream,
        two_term_joint_consistency=joint,
    )


def _preintegrated_rates(
    mapping: GecCcpMapping,
) -> tuple[list[dict[str, Any]] | None, dict[str, Any] | None]:
    if mapping.closure.reaction_model != FUNCTION_EEDF_PREINTEGRATED_INELASTIC:
        return None, None
    return _build_dense_function_eedf_rate_closure(
        mapping, input_mph=mapping.model.input_mph
    )


def _active_rate_rows(
    mapping: GecCcpMapping,
    preintegrated: list[dict[str, Any]] | None,
) -> list[dict[str, Any]] | None:
    if not _uses_external_rates(mapping.closure):
        return None
    return (
        preintegrated
        if preintegrated is not None
        else _read_csv(mapping.bundle.path / "rates_vs_mean_energy.csv")
    )


def _require_mc_censoring(censoring: dict[str, Any] | None) -> None:
    if not isinstance(censoring, dict) or (
        censoring.get("physics_qualification", {}).get("passed") is True
    ):
        return
    failed = censoring["physics_qualification"].get(
        "failed_censored_anchors", []
    )
    raise GecCcpWorkflowError(
        "Monte Carlo censored-rate qualification failed before COMSOL "
        f"execution: {json.dumps(failed, sort_keys=True, separators=(',', ':'))}"
    )


def _closure_support(
    mapping: GecCcpMapping,
    active_rate_rows: list[dict[str, Any]] | None,
) -> dict[str, Any] | None:
    if not (
        _uses_transport(mapping.closure)
        or _uses_external_rates(mapping.closure)
        or _uses_external_elastic_energy_loss(mapping.closure)
    ):
        return None
    support = _active_closure_mean_energy_support(
        mapping, preintegrated_rate_rows=active_rate_rows
    )
    if not support["passed"]:
        raise GecCcpWorkflowError(
            "active closure tables do not share a usable mean-energy "
            f"support: {support.get('reason')}"
        )
    return support


def _upstream_rate_consistency(
    mapping: GecCcpMapping, bundle: dict[str, Any]
) -> dict[str, Any] | None:
    if not (
        mapping.closure.reaction_model == "external_rates"
        and bundle.get("source") in {"two_term", "monte_carlo"}
        and bundle.get("upstream_eedf_evidence") is not None
    ):
        return None
    audit = _build_upstream_eedf_rate_consistency_audit(
        mapping, input_mph=mapping.model.input_mph
    )
    if not audit["passed"]:
        raise GecCcpWorkflowError(
            "external Swarm projected-EEDF/rate consistency audit failed: "
            + ", ".join(audit["failure_reasons"])
        )
    return audit


def _two_term_joint_consistency(
    mapping: GecCcpMapping,
    bundle: dict[str, Any],
    upstream: dict[str, Any] | None,
) -> dict[str, Any] | None:
    if not (
        bundle.get("source") == "two_term"
        and mapping.closure.electron_transport
        in {
            GEC_TWO_TERM_HYBRID_TRANSPORT_CLOSURE,
            GEC_RESTRICTED_TRANSPORT_CLOSURE,
        }
        and bundle.get("upstream_eedf_evidence") is not None
    ):
        return None
    audit = _build_two_term_joint_consistency_audit(
        mapping,
        input_mph=mapping.model.input_mph,
        upstream_rate_consistency=upstream,
    )
    if not audit["passed"]:
        raise GecCcpWorkflowError(
            "two-term Function-EEDF/transport joint consistency audit failed: "
            + ", ".join(audit["failure_reasons"])
        )
    return audit
