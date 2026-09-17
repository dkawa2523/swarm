"""Build and cross-check coupled GEC-CCP EEDF, rate, and transport closures."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any
from zipfile import ZipFile

import numpy as np
from scipy.interpolate import PchipInterpolator

from electron_swarm.core.constants import E_CHARGE_C, ELECTRON_MASS_KG
from electron_swarm.solvers.two_term.transport import (
    temporal_growth_effective_momentum_frequency,
    TWO_TERM_TEMPORAL_GROWTH_EFFECTIVE_MOMENTUM_MODEL,
)

from swarm_workflow.comsol.input.function_eedf import (
    evaluate_comsol_function_eedf_grid,
)
from swarm_workflow.comsol.input.mean_energy import MeanEnergyArgument
from swarm_workflow.comsol.models.gec_ccp.validation.bundle_guards import (
    GEC_ELASTIC_ENERGY_LOSS_COLUMN,
    GEC_ELASTIC_ENERGY_LOSS_TABLE,
    _validate_two_term_temporal_growth_transport_contract,
)
from swarm_workflow.comsol.models.gec_ccp.closure import (
    _active_rate_table,
    _active_transport_functions,
    _reaction_uses_preintegrated_rate,
    _thermal_diffusion_enabled,
    _uses_external_elastic_energy_loss,
    _uses_external_rates,
    _uses_function_eedf,
    _uses_hybrid_einstein_transport,
    _uses_transport,
)
from swarm_workflow.comsol.models.gec_ccp.contracts import (
    FUNCTION_EEDF_PREINTEGRATED_INELASTIC,
    GEC_LOOKUP_INTERPOLATION,
    GEC_RESTRICTED_TRANSPORT_CLOSURE,
    GEC_TWO_TERM_HYBRID_TRANSPORT_CLOSURE,
    FUNCTION_EEDF_TABLE,
    GecCcpMapping,
    GecCcpWorkflowError,
)
from swarm_workflow.comsol.models.gec_ccp.data import (
    _lookup,
    _read_csv,
    _required_float,
)
from swarm_workflow.comsol.models.gec_ccp.validation.eedf import (
    FUNCTION_EEDF_RATE_SIGNIFICANCE_FRACTION,
    UPSTREAM_EEDF_RATE_AUDIT_TABLE,
    _audit_function_eedf_source_rates,
    _audit_gec_argon_cross_section_identity,
    _integrate_native_function_eedf_rate,
    _reaction_cross_sections_from_model_xml,
    _read_function_eedf,
)
from swarm_workflow.tables.contracts import (
    TWO_TERM_MOMENTUM_CROSS_SECTION_FLOOR_M2,
)


GEC_TWO_TERM_JOINT_TRANSPORT_RELATIVE_ERROR_LIMIT = 0.04
GEC_TWO_TERM_TRANSPORT_PROCESS_TYPES = (
    "elastic",
    "excitation",
    "ionization",
)
# Closure arguments live on a logarithmic mean-energy axis.  Keep accepted
# operating points a few smoothing widths away from every active table edge;
# equality with an endpoint is not independent evidence against constant
# extrapolation.
GEC_CLOSURE_SUPPORT_MARGIN_LOG_FRACTION = 1.0e-6
GEC_CLOSURE_SUPPORT_MARGIN_LOG_ABSOLUTE = 5.0e-6


def _build_dense_function_eedf_rate_closure(
    mapping: GecCcpMapping,
    *,
    input_mph: Path,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Integrate the active Function-EEDF on its native mean-energy grid.

    The bundle rate table remains independent Swarm/MC evidence.  Hybrid
    inelastic rates used by COMSOL are instead evaluated from the exact active
    physical EEDF projection and the cross sections embedded in the input MPH.
    """

    if mapping.closure.reaction_model != FUNCTION_EEDF_PREINTEGRATED_INELASTIC:
        raise GecCcpWorkflowError(
            "dense Function-EEDF rate closure requires the hybrid reaction model"
        )
    function_spec = mapping.closure.function_eedf
    if function_spec is None:
        raise GecCcpWorkflowError(
            "dense Function-EEDF rate closure lacks its active EEDF table"
        )
    function_path = mapping.bundle.path / function_spec.table
    function_grid = _read_function_eedf(function_path)
    inelastic_reactions = tuple(
        reaction
        for reaction in mapping.reactions
        if _reaction_uses_preintegrated_rate(mapping.closure, reaction)
    )
    if not inelastic_reactions:
        raise GecCcpWorkflowError(
            "hybrid Function-EEDF closure has no preintegrated reactions"
        )
    try:
        with ZipFile(input_mph) as archive:
            model_xml = archive.read("dmodel.xml").decode("utf-8", errors="replace")
    except (OSError, KeyError) as exc:
        raise GecCcpWorkflowError(
            f"cannot read embedded reaction cross sections from {input_mph}"
        ) from exc
    cross_sections = _reaction_cross_sections_from_model_xml(
        model_xml,
        inelastic_reactions,
    )

    rows: list[dict[str, Any]] = []
    process_summary: dict[str, Any] = {}
    cross_section_hashes: dict[str, str] = {}
    for reaction in inelastic_reactions:
        cross_energy, cross_sigma = cross_sections[reaction.process_type]
        if (
            np.any(~np.isfinite(cross_energy))
            or np.any(~np.isfinite(cross_sigma))
            or np.any(np.diff(cross_energy) <= 0.0)
            or np.any(cross_sigma < 0.0)
            or not np.any(cross_sigma > 0.0)
        ):
            raise GecCcpWorkflowError(
                "embedded cross section is invalid for dense Function-EEDF "
                f"integration: {reaction.process_type}"
            )
        rates = np.asarray(
            [
                _integrate_native_function_eedf_rate(
                    function_grid,
                    float(mean),
                    cross_energy,
                    cross_sigma,
                )
                for mean in function_grid.mean_energies_eV
            ],
            dtype=float,
        )
        if np.any(~np.isfinite(rates)) or np.any(rates < 0.0):
            raise GecCcpWorkflowError(
                "active Function-EEDF integration produced an invalid rate "
                f"for {reaction.process_type}"
            )
        positive = np.flatnonzero(rates > 0.0)
        if positive.size == 0:
            raise GecCcpWorkflowError(
                "active Function-EEDF integration did not resolve a positive "
                f"rate for {reaction.process_type}; no rate floor is permitted"
            )
        zero = np.flatnonzero(rates == 0.0)
        support_start = int(zero[-1]) + 1 if zero.size else 0
        selected_rates = rates[support_start:]
        selected_means = function_grid.mean_energies_eV[support_start:]
        if selected_rates.size < 2:
            raise GecCcpWorkflowError(
                "active Function-EEDF integration needs at least two positive "
                f"rate points for {reaction.process_type}"
            )
        rows.extend(
            {
                "mean_energy_eV": float(mean),
                "process_type": reaction.process_type,
                "rate_coefficient_m3_s": float(rate),
            }
            for mean, rate in zip(
                selected_means,
                selected_rates,
                strict=True,
            )
        )
        cross_section_payload = {
            "electron_energy_eV": cross_energy.tolist(),
            "cross_section_m2": cross_sigma.tolist(),
        }
        cross_section_hashes[reaction.process_type] = hashlib.sha256(
            json.dumps(
                cross_section_payload,
                sort_keys=True,
                separators=(",", ":"),
                allow_nan=False,
            ).encode("utf-8")
        ).hexdigest()
        process_summary[reaction.process_type] = {
            "feature": reaction.feature,
            "points": int(selected_rates.size),
            "rows_omitted_before_contiguous_positive_support": support_start,
            "zero_rows_observed": int(zero.size),
            "positive_rows_discarded_before_support": int(
                np.count_nonzero(rates[:support_start] > 0.0)
            ),
            "positive_support_mean_energy_eV": [
                float(selected_means[0]),
                float(selected_means[-1]),
            ],
            "rate_range_m3_s": [
                float(np.min(selected_rates)),
                float(np.max(selected_rates)),
            ],
            "cross_section_points": int(cross_energy.size),
            "cross_section_sha256": cross_section_hashes[reaction.process_type],
        }

    canonical_rows = [
        {
            "mean_energy_eV": f"{float(row['mean_energy_eV']):.17e}",
            "process_type": str(row["process_type"]),
            "rate_coefficient_m3_s": (f"{float(row['rate_coefficient_m3_s']):.17e}"),
        }
        for row in rows
    ]
    closure_digest = hashlib.sha256(
        json.dumps(
            {
                "function_eedf_table_sha256": hashlib.sha256(
                    function_path.read_bytes()
                ).hexdigest(),
                "cross_section_sha256": cross_section_hashes,
                "rows": canonical_rows,
            },
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()
    return rows, {
        "source": "active_function_eedf_dense_reintegration",
        "active_function_eedf_table": function_spec.table,
        "active_function_eedf_table_sha256": hashlib.sha256(
            function_path.read_bytes()
        ).hexdigest(),
        "cross_section_source": "input_mph_embedded_xdata_ydata",
        "input_mph_sha256": hashlib.sha256(input_mph.read_bytes()).hexdigest(),
        "integration": ("sigma(E)*v(E)*sqrt(E)*f0(E,meanE) Gauss-Legendre quadrature"),
        "interpolation": GEC_LOOKUP_INTERPOLATION,
        "strictly_positive_without_floor": True,
        "leading_zero_rate_policy": (
            "use_only_the_final_contiguous_strictly_positive_suffix_after_the_"
            "last_unresolved_or_subthreshold_zero; no floor is added; constant_"
            "extrapolation_is_a_Newton_trial_guard_only_and_the_accepted_"
            "solution_must_remain_inside_every_positive_rate_support"
        ),
        "mean_energy_grid": "active_function_eedf_native_rows",
        "mean_energy_grid_points": int(function_grid.mean_energies_eV.size),
        "mean_energy_range_eV": [
            float(function_grid.mean_energies_eV[0]),
            float(function_grid.mean_energies_eV[-1]),
        ],
        "raw_bundle_rates_role": "independent_source_and_uncertainty_evidence",
        "processes": process_summary,
        "closure_sha256": closure_digest,
    }


def _build_upstream_eedf_rate_consistency_audit(
    mapping: GecCcpMapping,
    *,
    input_mph: Path,
) -> dict[str, Any]:
    """Audit inactive projected EEDF evidence against active external rates."""

    if (
        mapping.closure.reaction_model != "external_rates"
        or mapping.bundle.expected_source not in {"two_term", "monte_carlo"}
    ):
        raise GecCcpWorkflowError(
            "upstream EEDF/rate consistency requires two_term or "
            "monte_carlo external_rates"
        )
    function_path = mapping.bundle.path / FUNCTION_EEDF_TABLE
    function_grid = _read_function_eedf(function_path)
    try:
        with ZipFile(input_mph) as archive:
            model_xml = archive.read("dmodel.xml").decode("utf-8", errors="replace")
    except (OSError, KeyError) as exc:
        raise GecCcpWorkflowError(
            f"cannot read embedded cross sections from {input_mph}"
        ) from exc
    cross_sections = _reaction_cross_sections_from_model_xml(
        model_xml,
        mapping.reactions,
    )
    cross_section_identity = _audit_gec_argon_cross_section_identity(
        cross_sections,
        mapping.reactions,
    )
    inelastic_reactions = tuple(
        reaction
        for reaction in mapping.reactions
        if reaction.process_type in {"excitation", "ionization"}
        and _reaction_uses_preintegrated_rate(mapping.closure, reaction)
    )
    active_inelastic_processes = {
        process
        for process in mapping.closure.external_rate_processes
        if process in {"excitation", "ionization"}
    }
    if (
        not active_inelastic_processes
        or {reaction.process_type for reaction in inelastic_reactions}
        != active_inelastic_processes
    ):
        raise GecCcpWorkflowError(
            "external-rate EEDF audit requires at least one active inelastic process"
        )
    rate_consistency = _audit_function_eedf_source_rates(
        mapping,
        function_grid,
        cross_sections,
        inelastic_reactions,
        rate_table_name="rates_vs_mean_energy.csv",
        audit_csv_name=UPSTREAM_EEDF_RATE_AUDIT_TABLE,
    )
    monte_carlo = mapping.bundle.expected_source == "monte_carlo"
    owner_prefix = (
        "monte_carlo_direct_trajectory_rates"
        if monte_carlo
        else "two_term_integrated_rates"
    )
    active_owner = f"{owner_prefix}:rates_vs_mean_energy.csv"
    if active_inelastic_processes != {"excitation", "ionization"}:
        active_owner += ":" + ",".join(sorted(active_inelastic_processes))
    failure_reasons: list[str] = []
    if not cross_section_identity["passed"]:
        failed_processes = [
            name
            for name, result in cross_section_identity.get("processes", {}).items()
            if not result.get("passed", False)
        ]
        failure_reasons.append(
            "cross_section_identity=" + (",".join(failed_processes) or "failed")
        )
    for process_type, result in rate_consistency["processes"].items():
        if result["passed"]:
            continue
        failure_reasons.append(
            f"rate_{process_type}="
            f"max:{result['maximum_relative_error_significant']:.3%},"
            f"p95:{result['p95_relative_error_significant']:.3%},"
            f"nrmse:{result['normalized_rmse']:.3%}"
        )
    passed = bool(cross_section_identity["passed"] and rate_consistency["passed"])
    return {
        "status": "passed" if passed else "failed",
        "passed": passed,
        "source": mapping.bundle.expected_source,
        "projected_eedf": {
            "path": str(function_path),
            "sha256": hashlib.sha256(function_path.read_bytes()).hexdigest(),
            "active_in_comsol": False,
            "role": "inactive_upstream_consistency_evidence",
        },
        "active_inelastic_rate_owner": active_owner,
        "reintegration_role": "fail_closed_consistency_audit_only",
        "reintegration_modifies_active_rates": False,
        "cross_section_identity": cross_section_identity,
        "rate_consistency": rate_consistency,
        "magnitude_policy": {
            "significance_fraction_of_process_peak": (
                FUNCTION_EEDF_RATE_SIGNIFICANCE_FRACTION
            ),
            "below_significance_floor": ("record_only_for_pointwise_relative_error"),
            "artificial_rate_floor": False,
            "censored_zero_handling": (
                "preserve_zero_and_separately_retain_MC_upper_bound_evidence"
                if monte_carlo
                else "not_applicable"
            ),
        },
        "failure_reasons": failure_reasons,
    }


@dataclass(frozen=True)
class _TwoTermAuditContext:
    mapping: GecCcpMapping
    input_mph: Path
    function_path: Path
    transport_path: Path
    function_grid: Any
    transport_rows: list[dict[str, Any]]
    temporal_growth_contract: dict[str, Any] | None
    cross_sections: dict[str, tuple[np.ndarray, np.ndarray]]
    cross_section_identity: dict[str, Any]
    external_rate_consistency: dict[str, Any] | None


@dataclass(frozen=True)
class _TransportQuadrature:
    energy: np.ndarray
    weights: np.ndarray
    left: np.ndarray
    right: np.ndarray
    half_width: np.ndarray
    sample_energy: np.ndarray
    interpolation_fraction: np.ndarray
    sigma_m: np.ndarray
    cross_section_summary: dict[str, Any]


@dataclass(frozen=True)
class _ActiveEedfSamples:
    means: np.ndarray
    f0_rows: np.ndarray
    gas_densities: np.ndarray
    growth_frequencies: np.ndarray


@dataclass(frozen=True)
class _TransportComparison:
    passed: bool
    coefficients: dict[str, Any]
    evidence_sha256: str


def _build_two_term_joint_consistency_audit(
    mapping: GecCcpMapping,
    *,
    input_mph: Path,
    upstream_rate_consistency: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Reintegrate two-term transport from the projected upstream EEDF."""

    context = _load_two_term_audit_context(
        mapping,
        input_mph=input_mph,
        upstream_rate_consistency=upstream_rate_consistency,
    )
    quadrature = _build_transport_quadrature(context)
    samples = _active_eedf_samples(context, quadrature.energy)
    reintegrated = _reintegrate_two_term_transport(
        context,
        quadrature,
        samples,
    )
    comparison = _compare_two_term_transport(
        context,
        samples.means,
        reintegrated,
    )
    kernel_contract, definitions = _transport_kernel_contract(context)
    failure_reasons = _two_term_failure_reasons(
        comparison.coefficients,
        context.cross_section_identity,
        context.external_rate_consistency,
    )
    return {
        "status": "passed" if comparison.passed else "failed",
        "passed": comparison.passed,
        "applicability": ("two_term_external_transport_with_projected_EEDF_evidence"),
        "upstream_projected_eedf_evidence": {
            "path": str(context.function_path),
            "sha256": hashlib.sha256(context.function_path.read_bytes()).hexdigest(),
            "representation": "physical_2d_piecewise_linear_C0",
            "energy_grid_points": int(quadrature.energy.size),
            "mean_energy_grid_points": int(samples.means.size),
        },
        "transport_table": {
            "path": str(context.transport_path),
            "sha256": hashlib.sha256(context.transport_path.read_bytes()).hexdigest(),
            "interpolation": "log_PCHIP_piecewise_cubic_Hermite_C1",
        },
        "cross_section_source": {
            "kind": "input_mph_embedded_electron_impact_xdata_ydata",
            "input_mph_sha256": hashlib.sha256(input_mph.read_bytes()).hexdigest(),
            "processes": quadrature.cross_section_summary,
            "combination": (
                "sum_raw_elastic_excitation_ionization_cross_sections_then_"
                "single_aggregate_floor"
            ),
            "aggregate_floor_m2": TWO_TERM_MOMENTUM_CROSS_SECTION_FLOOR_M2,
        },
        "cross_section_identity": context.cross_section_identity,
        "external_rate_consistency": context.external_rate_consistency,
        "closure_ownership": {
            "projected_eedf": (
                "inactive_upstream_consistency_evidence"
                if context.external_rate_consistency is not None
                else "active_or_non_external_rate_closure"
            ),
            "inelastic_rate_owner": (
                context.external_rate_consistency["active_inelastic_rate_owner"]
                if context.external_rate_consistency is not None
                else "closure_specific"
            ),
            "reintegration_modifies_active_rates": False,
        },
        "integration": {
            "energy_interpolation": "exact_active_piecewise_linear_f0",
            "quadrature": ("8_point_Gauss_Legendre_per_active_energy_interval"),
            "definitions": definitions,
            "gamma": "sqrt(2 electron_charge/electron_mass)",
        },
        "transport_kernel_contract": kernel_contract,
        "temporal_growth_transport": context.temporal_growth_contract,
        "relative_error_limit": (GEC_TWO_TERM_JOINT_TRANSPORT_RELATIVE_ERROR_LIMIT),
        "tolerance_basis": (
            "bounds C0 Function-EEDF projection and derivative sensitivity; "
            "it is not a nonlinear-solver or coefficient-smoothing tolerance"
        ),
        "coefficients": comparison.coefficients,
        "failure_reasons": failure_reasons,
        "evidence_rows_sha256": comparison.evidence_sha256,
        "thermal_diffusion_enabled": _thermal_diffusion_enabled(mapping.closure),
        "scope_limit": (
            "active diagonal f0-moment coefficients only; this does not "
            "identify the full number-density/mean-energy two-gradient "
            "response matrix"
        ),
    }


def _load_two_term_audit_context(
    mapping: GecCcpMapping,
    *,
    input_mph: Path,
    upstream_rate_consistency: dict[str, Any] | None,
) -> _TwoTermAuditContext:
    if (
        mapping.bundle.expected_source != "two_term"
        or mapping.closure.electron_transport
        not in {
            GEC_TWO_TERM_HYBRID_TRANSPORT_CLOSURE,
            GEC_RESTRICTED_TRANSPORT_CLOSURE,
        }
    ):
        raise GecCcpWorkflowError(
            "two-term joint consistency audit is only valid for a two_term "
            "hybrid or full-transport closure"
        )
    function_spec = mapping.closure.function_eedf
    function_path = mapping.bundle.path / (
        function_spec.table if function_spec is not None else FUNCTION_EEDF_TABLE
    )
    transport_path = mapping.bundle.path / "transport_vs_mean_energy.csv"
    function_grid = _read_function_eedf(function_path)
    transport_rows = _read_csv(transport_path)
    manifest = json.loads(
        (mapping.bundle.path / "manifest.json").read_text(encoding="utf-8")
    )
    source_policy = manifest.get("source_policy")
    if not isinstance(source_policy, dict):
        raise GecCcpWorkflowError(
            "two-term joint transport audit requires bundle source_policy"
        )
    temporal_growth_contract = _validate_two_term_temporal_growth_transport_contract(
        mapping,
        source_policy,
        transport_rows,
    )
    try:
        with ZipFile(input_mph) as archive:
            model_xml = archive.read("dmodel.xml").decode("utf-8", errors="replace")
    except (OSError, KeyError) as exc:
        raise GecCcpWorkflowError(
            f"cannot read embedded cross sections from {input_mph}"
        ) from exc
    cross_sections = _reaction_cross_sections_from_model_xml(
        model_xml,
        mapping.reactions,
    )
    if (
        mapping.closure.reaction_model == "external_rates"
        and upstream_rate_consistency is None
    ):
        upstream_rate_consistency = _build_upstream_eedf_rate_consistency_audit(
            mapping,
            input_mph=input_mph,
        )
    cross_section_identity = (
        upstream_rate_consistency["cross_section_identity"]
        if upstream_rate_consistency is not None
        else _audit_gec_argon_cross_section_identity(
            cross_sections,
            mapping.reactions,
        )
    )
    if set(cross_sections) != set(GEC_TWO_TERM_TRANSPORT_PROCESS_TYPES):
        raise GecCcpWorkflowError(
            "two-term GEC transport audit requires exactly the embedded "
            "elastic, excitation, and ionization cross sections"
        )
    external_rate_consistency = None
    if upstream_rate_consistency is not None:
        external_rate_consistency = dict(upstream_rate_consistency["rate_consistency"])
        external_rate_consistency.update(
            {
                "eedf_active_in_comsol": False,
                "eedf_role": "inactive_upstream_consistency_evidence",
                "active_inelastic_rate_owner": upstream_rate_consistency[
                    "active_inelastic_rate_owner"
                ],
                "reintegration_role": ("fail_closed_consistency_audit_only"),
                "active_rate_values_modified": False,
            }
        )
    return _TwoTermAuditContext(
        mapping=mapping,
        input_mph=input_mph,
        function_path=function_path,
        transport_path=transport_path,
        function_grid=function_grid,
        transport_rows=transport_rows,
        temporal_growth_contract=temporal_growth_contract,
        cross_sections=cross_sections,
        cross_section_identity=cross_section_identity,
        external_rate_consistency=external_rate_consistency,
    )


def _build_transport_quadrature(
    context: _TwoTermAuditContext,
) -> _TransportQuadrature:
    energy = np.asarray(
        context.function_grid.electron_energies_eV,
        dtype=float,
    )
    if energy.size < 3 or np.any(np.diff(energy) <= 0.0):
        raise GecCcpWorkflowError(
            "two-term joint consistency audit requires an increasing "
            "Function-EEDF energy grid"
        )
    quadrature_nodes, quadrature_weights = np.polynomial.legendre.leggauss(8)
    left = energy[:-1, None]
    right = energy[1:, None]
    half_width = 0.5 * (right - left)
    sample_energy = 0.5 * (right + left) + half_width * quadrature_nodes
    interpolation_fraction = (sample_energy - left) / (right - left)
    sigma_m = np.zeros_like(sample_energy)
    cross_section_summary: dict[str, Any] = {}
    for process_type in GEC_TWO_TERM_TRANSPORT_PROCESS_TYPES:
        cross_energy, cross_sigma = context.cross_sections[process_type]
        if (
            cross_energy.size < 2
            or cross_energy.shape != cross_sigma.shape
            or np.any(~np.isfinite(cross_energy))
            or np.any(~np.isfinite(cross_sigma))
            or np.any(np.diff(cross_energy) <= 0.0)
            or np.any(cross_sigma < 0.0)
            or cross_energy[0] > energy[0]
            or cross_energy[-1] < energy[-1]
        ):
            raise GecCcpWorkflowError(
                "embedded cross section does not cover the active "
                f"Function-EEDF energy support: {process_type}"
            )
        interpolated = np.interp(
            sample_energy,
            cross_energy,
            cross_sigma,
            left=0.0,
            right=0.0,
        )
        sigma_m += interpolated
        cross_section_payload = {
            "electron_energy_eV": cross_energy.tolist(),
            "cross_section_m2": cross_sigma.tolist(),
        }
        cross_section_summary[process_type] = {
            "points": int(cross_energy.size),
            "energy_range_eV": [
                float(cross_energy[0]),
                float(cross_energy[-1]),
            ],
            "sha256": hashlib.sha256(
                json.dumps(
                    cross_section_payload,
                    sort_keys=True,
                    separators=(",", ":"),
                    allow_nan=False,
                ).encode("utf-8")
            ).hexdigest(),
        }
    sigma_m = np.maximum(
        sigma_m,
        TWO_TERM_MOMENTUM_CROSS_SECTION_FLOOR_M2,
    )
    return _TransportQuadrature(
        energy=energy,
        weights=quadrature_weights,
        left=left,
        right=right,
        half_width=half_width,
        sample_energy=sample_energy,
        interpolation_fraction=interpolation_fraction,
        sigma_m=sigma_m,
        cross_section_summary=cross_section_summary,
    )


def _active_eedf_samples(
    context: _TwoTermAuditContext,
    energy: np.ndarray,
) -> _ActiveEedfSamples:
    if context.temporal_growth_contract is not None:
        means = np.asarray(
            [_required_float(row, "mean_energy_eV") for row in context.transport_rows],
            dtype=float,
        )
        f0_rows = np.asarray(
            [
                evaluate_comsol_function_eedf_grid(
                    context.function_grid,
                    energy,
                    float(mean),
                )
                for mean in means
            ],
            dtype=float,
        )
        gas_densities = np.asarray(
            [
                _required_float(row, "gas_number_density_m3")
                for row in context.transport_rows
            ],
            dtype=float,
        )
        growth_frequencies = np.asarray(
            [
                _required_float(row, "temporal_growth_frequency_s_inv")
                for row in context.transport_rows
            ],
            dtype=float,
        )
    else:
        means = np.asarray(
            context.function_grid.mean_energies_eV,
            dtype=float,
        )
        f0_rows = np.asarray(
            context.function_grid.values_eV_m32,
            dtype=float,
        )
        gas_densities = np.full(means.shape, math.nan)
        growth_frequencies = np.zeros(means.shape, dtype=float)
    if means.size < 2 or np.any(np.diff(means) <= 0.0):
        raise GecCcpWorkflowError(
            "two-term joint consistency audit requires an increasing "
            "Function-EEDF mean-energy grid"
        )
    return _ActiveEedfSamples(
        means=means,
        f0_rows=f0_rows,
        gas_densities=gas_densities,
        growth_frequencies=growth_frequencies,
    )


def _reintegrate_two_term_transport(
    context: _TwoTermAuditContext,
    quadrature: _TransportQuadrature,
    samples: _ActiveEedfSamples,
) -> np.ndarray:
    gamma = math.sqrt(2.0 * E_CHARGE_C / ELECTRON_MASS_KG)
    reintegrated = np.empty((samples.means.size, 4), dtype=float)
    for mean_index, mean_energy in enumerate(samples.means):
        f0 = np.asarray(samples.f0_rows[mean_index], dtype=float)
        f0_left = f0[:-1, None]
        f0_right = f0[1:, None]
        f0_sample = f0_left + (f0_right - f0_left) * (quadrature.interpolation_fraction)
        f0_derivative = (f0_right - f0_left) / (quadrature.right - quadrature.left)
        inverse_sigma = 1.0 / quadrature.sigma_m
        if context.temporal_growth_contract is not None:
            density = float(samples.gas_densities[mean_index])
            speed = gamma * np.sqrt(quadrature.sample_energy)
            base_frequency = density * speed * quadrature.sigma_m
            try:
                effective_frequency = temporal_growth_effective_momentum_frequency(
                    base_frequency.reshape(-1),
                    float(samples.growth_frequencies[mean_index]),
                ).reshape(base_frequency.shape)
            except FloatingPointError as exc:
                raise GecCcpWorkflowError(
                    "two-term PT effective momentum frequency is invalid "
                    f"at mean energy {mean_energy:.17g} eV"
                ) from exc
            inverse_sigma = density * speed / effective_frequency

        def integrate(values: np.ndarray) -> float:
            return float(np.sum(quadrature.half_width * quadrature.weights * values))

        reintegrated[mean_index] = (
            -gamma
            / 3.0
            * integrate(quadrature.sample_energy * inverse_sigma * f0_derivative),
            gamma
            / 3.0
            * integrate(quadrature.sample_energy * inverse_sigma * f0_sample),
            -gamma
            / (3.0 * mean_energy)
            * integrate(
                quadrature.sample_energy
                * quadrature.sample_energy
                * inverse_sigma
                * f0_derivative
            ),
            gamma
            / (3.0 * mean_energy)
            * integrate(
                quadrature.sample_energy
                * quadrature.sample_energy
                * inverse_sigma
                * f0_sample
            ),
        )
    if np.any(~np.isfinite(reintegrated)) or np.any(reintegrated <= 0.0):
        raise GecCcpWorkflowError(
            "active Function-EEDF reintegration produced nonpositive or "
            "nonfinite two-term transport"
        )
    return reintegrated


def _compare_two_term_transport(
    context: _TwoTermAuditContext,
    means: np.ndarray,
    reintegrated: np.ndarray,
) -> _TransportComparison:
    if _uses_hybrid_einstein_transport(context.mapping.closure):
        comparisons = (
            ("muN", "reduced_mobility_m2_V_s_m3", 0),
            (
                "muenN",
                "reduced_electron_energy_mobility_m2_V_s_m3",
                2,
            ),
            (
                "DenN",
                "reduced_electron_energy_diffusion_m2_s_m3",
                3,
            ),
        )
    else:
        comparisons = (
            ("muN", "reduced_mobility_m2_V_s_m3", 0),
            ("DeN_L", "reduced_diffusion_L_m2_s_m3", 1),
            ("DeN_T", "reduced_diffusion_T_m2_s_m3", 1),
            (
                "muenN",
                "reduced_electron_energy_mobility_m2_V_s_m3",
                2,
            ),
            (
                "DenN_L",
                "reduced_electron_energy_diffusion_L_m2_s_m3",
                3,
            ),
            (
                "DenN_T",
                "reduced_electron_energy_diffusion_T_m2_s_m3",
                3,
            ),
        )
    coefficient_summary: dict[str, Any] = {}
    evidence_rows: list[dict[str, str]] = []
    passed = bool(
        context.cross_section_identity["passed"]
        and (
            context.external_rate_consistency is None
            or context.external_rate_consistency["passed"]
        )
    )
    for name, column, reintegrated_index in comparisons:
        transport_means, _ = _lookup(
            context.transport_rows,
            "mean_energy_eV",
            column,
        )
        reference = _independent_log_piecewise_cubic_transport(
            context.transport_rows,
            column,
            means,
            support_minimum_eV=float(transport_means[0]),
            support_maximum_eV=float(transport_means[-1]),
        )
        if np.any(~np.isfinite(reference)) or np.any(reference <= 0.0):
            raise GecCcpWorkflowError(
                f"transport spline produced invalid values for {column}"
            )
        calculated = reintegrated[:, reintegrated_index]
        relative_error = np.abs(calculated / reference - 1.0)
        worst_index = int(np.argmax(relative_error))
        maximum_error = float(relative_error[worst_index])
        coefficient_passed = bool(
            maximum_error <= GEC_TWO_TERM_JOINT_TRANSPORT_RELATIVE_ERROR_LIMIT
        )
        passed = passed and coefficient_passed
        coefficient_summary[name] = {
            "transport_column": column,
            "samples": int(means.size),
            "maximum_relative_error": maximum_error,
            "p95_relative_error": float(np.quantile(relative_error, 0.95)),
            "worst_mean_energy_eV": float(means[worst_index]),
            "passed": coefficient_passed,
        }
        evidence_rows.extend(
            {
                "mean_energy_eV": f"{float(mean):.17e}",
                "coefficient": name,
                "reintegrated": f"{float(actual):.17e}",
                "transport_spline": f"{float(expected):.17e}",
            }
            for mean, actual, expected in zip(
                means,
                calculated,
                reference,
                strict=True,
            )
        )
    evidence_sha256 = hashlib.sha256(
        json.dumps(
            evidence_rows,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()
    return _TransportComparison(
        passed=passed,
        coefficients=coefficient_summary,
        evidence_sha256=evidence_sha256,
    )


def _transport_kernel_contract(
    context: _TwoTermAuditContext,
) -> tuple[dict[str, Any], dict[str, str]]:
    temporal_growth = context.temporal_growth_contract
    denominator = "sigma_m_eff" if temporal_growth is not None else "sigma_m"
    definitions = {
        "muN": f"-gamma/3 integral((E/{denominator}) df0/dE dE)",
        "DeN": f"gamma/3 integral((E/{denominator}) f0 dE)",
        "muenN": (f"-gamma/(3 meanE) integral((E^2/{denominator}) df0/dE dE)"),
        "DenN": (f"gamma/(3 meanE) integral((E^2/{denominator}) f0 dE)"),
    }
    if temporal_growth is not None:
        definitions["sigma_m_eff"] = "(nu_m+growth_frequency)/(N*gamma*sqrt(E))"
    kernel_contract: dict[str, Any] = {
        "schema": "two_term_active_function_transport.v1",
        "process_selection": list(GEC_TWO_TERM_TRANSPORT_PROCESS_TYPES),
        "effective_momentum_model": (
            TWO_TERM_TEMPORAL_GROWTH_EFFECTIVE_MOMENTUM_MODEL
            if temporal_growth is not None
            else "sum_raw_process_cross_sections_then_single_floor"
        ),
        "aggregate_floor_m2": TWO_TERM_MOMENTUM_CROSS_SECTION_FLOOR_M2,
        "temporal_growth_schema": (
            temporal_growth["schema"] if temporal_growth is not None else None
        ),
        "cross_section_interpolation": (
            "piecewise_linear_with_full_active_energy_support_required"
        ),
        "active_f0_interpolation": "piecewise_linear_C0",
        "quadrature": ("8_point_Gauss_Legendre_per_active_energy_interval"),
        "transport_interpolation": ("log_PCHIP_piecewise_cubic_Hermite_C1"),
        "definitions": definitions,
        "thermal_diffusion_enabled": _thermal_diffusion_enabled(
            context.mapping.closure
        ),
    }
    kernel_contract["sha256"] = hashlib.sha256(
        json.dumps(
            kernel_contract,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
    ).hexdigest()
    return kernel_contract, definitions


def _two_term_failure_reasons(
    coefficient_summary: dict[str, Any],
    cross_section_identity: dict[str, Any],
    external_rate_consistency: dict[str, Any] | None,
) -> list[str]:
    failure_reasons = [
        f"{name}={result['maximum_relative_error']:.3%}"
        for name, result in coefficient_summary.items()
        if not result["passed"]
    ]
    if not cross_section_identity["passed"]:
        failed_processes = [
            name
            for name, result in cross_section_identity.get("processes", {}).items()
            if not result.get("passed", False)
        ]
        failure_reasons.append(
            "cross_section_identity=" + (",".join(failed_processes) or "failed")
        )
    if (
        external_rate_consistency is not None
        and not external_rate_consistency["passed"]
    ):
        for process_type, result in external_rate_consistency["processes"].items():
            if result["passed"]:
                continue
            failure_reasons.append(
                f"rate_{process_type}="
                f"max:{result['maximum_relative_error_significant']:.3%},"
                f"p95:{result['p95_relative_error_significant']:.3%},"
                f"nrmse:{result['normalized_rmse']:.3%}"
            )
    return failure_reasons


def _independent_log_piecewise_cubic_transport(
    rows: list[dict[str, str]],
    column: str,
    mean_energy_eV: np.ndarray,
    *,
    support_minimum_eV: float,
    support_maximum_eV: float,
) -> np.ndarray:
    """Evaluate a bundle column without using any saved COMSOL function."""

    x_values, y_values = _lookup(rows, "mean_energy_eV", column)
    if any(value <= 0.0 for value in y_values):
        raise GecCcpWorkflowError(
            f"independent log-piecewise-cubic audit requires positive {column}"
        )
    smooth = MeanEnergyArgument(support_minimum_eV, support_maximum_eV).log_values(
        mean_energy_eV
    )
    interpolator = PchipInterpolator(
        np.log(np.asarray(x_values)),
        np.log(np.asarray(y_values)),
        extrapolate=False,
    )
    return np.exp(interpolator(smooth))


def _active_closure_mean_energy_support(
    mapping: GecCcpMapping,
    *,
    preintegrated_rate_rows: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Return the fail-closed intersection of every active closure table."""

    supports: dict[str, list[float]] = {}
    if _uses_transport(mapping.closure):
        rows = _read_csv(mapping.bundle.path / "transport_vs_mean_energy.csv")
        for _, name, column, _ in _active_transport_functions(
            mapping.closure, source=mapping.bundle.expected_source
        ):
            energies, _ = _lookup(rows, "mean_energy_eV", column)
            supports[f"transport:{name}"] = [min(energies), max(energies)]

    function_support: list[float] | None = None
    if _uses_function_eedf(mapping.closure):
        spec = mapping.closure.function_eedf
        if spec is None:
            raise GecCcpWorkflowError("missing Function-EEDF specification")
        grid = _read_function_eedf(mapping.bundle.path / spec.table)
        function_support = [
            float(grid.mean_energies_eV[0]),
            float(grid.mean_energies_eV[-1]),
        ]
        supports["function_eedf"] = function_support

    if _uses_external_rates(mapping.closure):
        if mapping.closure.reaction_model == FUNCTION_EEDF_PREINTEGRATED_INELASTIC:
            if function_support is None:
                raise GecCcpWorkflowError(
                    "hybrid Function-EEDF rates lack a Function-EEDF support"
                )
            dense_rows = preintegrated_rate_rows
            if dense_rows is None:
                dense_rows, _ = _build_dense_function_eedf_rate_closure(
                    mapping,
                    input_mph=mapping.model.input_mph,
                )
            for reaction in mapping.reactions:
                if _reaction_uses_preintegrated_rate(mapping.closure, reaction):
                    selected, _ = _active_rate_table(
                        mapping,
                        dense_rows,
                        process_type=reaction.process_type,
                    )
                    energies, _ = _lookup(
                        selected, "mean_energy_eV", "rate_coefficient_m3_s"
                    )
                    supports[f"rate:{reaction.process_type}"] = [
                        min(energies),
                        max(energies),
                    ]
        else:
            rows = (
                preintegrated_rate_rows
                if preintegrated_rate_rows is not None
                else _read_csv(mapping.bundle.path / "rates_vs_mean_energy.csv")
            )
            for reaction in mapping.reactions:
                if not _reaction_uses_preintegrated_rate(mapping.closure, reaction):
                    continue
                selected, _ = _active_rate_table(
                    mapping,
                    rows,
                    process_type=reaction.process_type,
                )
                energies, _ = _lookup(
                    selected, "mean_energy_eV", "rate_coefficient_m3_s"
                )
                supports[f"rate:{reaction.process_type}"] = [
                    min(energies),
                    max(energies),
                ]

    elastic_energy_loss_support: list[float] | None = None
    if _uses_external_elastic_energy_loss(mapping.closure):
        elastic_rows = _read_csv(mapping.bundle.path / GEC_ELASTIC_ENERGY_LOSS_TABLE)
        elastic_energies, _ = _lookup(
            elastic_rows,
            "mean_energy_eV",
            GEC_ELASTIC_ENERGY_LOSS_COLUMN,
        )
        elastic_energy_loss_support = [min(elastic_energies), max(elastic_energies)]
        supports["elastic_energy_loss"] = elastic_energy_loss_support

    if not supports:
        raise GecCcpWorkflowError("external closure has no active table support")
    transport_supports = {
        name: value for name, value in supports.items() if name.startswith("transport:")
    }
    rate_supports = {
        name.removeprefix("rate:"): value
        for name, value in supports.items()
        if name.startswith("rate:")
    }
    transport_intersection = (
        [
            max(value[0] for value in transport_supports.values()),
            min(value[1] for value in transport_supports.values()),
        ]
        if transport_supports
        else None
    )
    lower = max(value[0] for value in supports.values())
    upper = min(value[1] for value in supports.values())
    if not 0.0 < lower < upper:
        return {
            "passed": False,
            "supports_mean_energy_eV": supports,
            "reason": "active closure tables have no positive common support",
        }
    log_lower = math.log(lower)
    log_upper = math.log(upper)
    log_margin = max(
        GEC_CLOSURE_SUPPORT_MARGIN_LOG_ABSOLUTE,
        GEC_CLOSURE_SUPPORT_MARGIN_LOG_FRACTION * (log_upper - log_lower),
    )
    if 2.0 * log_margin >= log_upper - log_lower:
        return {
            "passed": False,
            "supports_mean_energy_eV": supports,
            "common_intersection_mean_energy_eV": [lower, upper],
            "reason": "common support is too narrow for an interior margin",
        }
    return {
        "passed": True,
        "supports_mean_energy_eV": supports,
        "transport_intersection_mean_energy_eV": transport_intersection,
        "rate_supports_mean_energy_eV": rate_supports,
        "elastic_energy_loss_support_mean_energy_eV": (elastic_energy_loss_support),
        "guard_policy": (
            "separate_transport_rate_and_elastic_energy_loss_guards"
            if mapping.bundle.expected_source == "monte_carlo"
            else "shared_positive_log_closure_guard"
        ),
        "common_intersection_mean_energy_eV": [lower, upper],
        "accepted_interior_mean_energy_eV": [
            math.exp(log_lower + log_margin),
            math.exp(log_upper - log_margin),
        ],
        "endpoint_margin_log_fraction": (GEC_CLOSURE_SUPPORT_MARGIN_LOG_FRACTION),
        "endpoint_margin_log_absolute": (GEC_CLOSURE_SUPPORT_MARGIN_LOG_ABSOLUTE),
    }
