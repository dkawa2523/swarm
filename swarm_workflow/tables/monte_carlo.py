"""Monte Carlo table qualification and statistical evidence materialization."""

from __future__ import annotations

import json
import math
import sqlite3
from typing import Any

import electron_swarm.solvers.monte_carlo.evidence as _mc_evidence

from . import contracts as _contracts
from . import repository as _repository
from ..campaign.aggregate import RATE_SCALAR
from ..quality.monte_carlo.contracts import (
    MC_EEDF_RATE_ENSEMBLE_RELATIVE_TOLERANCE,
    MC_TRANSPORT_DIFFUSION_LAG_RELATIVE_TOLERANCE,
    MC_TRANSPORT_MIN_ENSEMBLE_REPLICATES,
    MC_TRANSPORT_STATE_RELATIVE_TOLERANCE,
    mc_function_eedf_restricted_lmea_failure_reasons,
    paired_log_ratio_ci95_bound,
    quality_failure_reasons,
)
from ..quality.monte_carlo.policy import failure_axes_for_reasons
from ..quality.policy import QualityThresholds
from ..campaign.repository import read_metadata


def _assess_monte_carlo(
    connection: sqlite3.Connection,
    mixture_id: int,
    *,
    quality_thresholds: QualityThresholds,
    qualification_profile: str,
) -> _contracts.MonteCarloEvidence:
    all_quality = _repository._load_quality(
        connection, _contracts.MC_SOLVER, mixture_id
    )
    estimator_schema = read_metadata(connection).get(
        "mc_transport_estimator_schema_version"
    )
    aggregate_eligible = {
        float(row["E_over_N_Td"])
        for row in all_quality
        if int(row.get("aggregate_quality_passed", 0)) == 1
    }
    full_transport_eligible = {
        float(row["E_over_N_Td"])
        for row in all_quality
        if int(row.get("aggregate_quality_passed", 0)) == 1
        and int(row.get("solver_transport_qualified", 0)) == 1
    }
    active_failure_reasons: dict[float, list[str]] = {}
    if (
        qualification_profile
        == _contracts.MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
    ):
        for row in all_quality:
            e_over_n = float(row["E_over_N_Td"])
            try:
                active_failure_reasons[e_over_n] = (
                    mc_function_eedf_restricted_lmea_failure_reasons(
                        row,
                        thresholds=quality_thresholds,
                        estimator_schema=estimator_schema,
                    )
                )
            except (KeyError, TypeError, ValueError, OverflowError) as exc:
                raise _contracts.TableBuildError(
                    "Monte Carlo restricted-LMEA quality evidence is invalid "
                    f"at {e_over_n:.17g} Td"
                ) from exc
        active_eligible = {
            e_over_n
            for e_over_n, reasons in active_failure_reasons.items()
            if not reasons
        }
    else:
        active_eligible = set(full_transport_eligible)
        for row in all_quality:
            e_over_n = float(row["E_over_N_Td"])
            active_failure_reasons[e_over_n] = quality_failure_reasons(
                row, "failure_reasons_json"
            )
    rate_failures = _mc_eedf_rate_consistency_failures(
        connection,
        mixture_id,
        minimum_process_peak_fraction=(
            quality_thresholds.required_rate_min_process_peak_fraction
        ),
        # Rate closure is an independent observable.  Evaluating it only on
        # anchors that survived mobility/growth/transport gates changes the
        # process peaks, so an unrelated failure can manufacture (or hide) a
        # rate failure at another field.  All structurally usable aggregate
        # anchors must therefore be assessed before profiles compose gates.
        allowed_e=aggregate_eligible,
    )
    censored_failures = _mc_censored_rate_relevance_failures(
        connection,
        mixture_id,
        allowed_e=aggregate_eligible,
        minimum_fraction=(
            quality_thresholds.required_rate_min_process_peak_fraction
        ),
    )
    for anchor, failures in censored_failures.items():
        rate_failures[anchor] = tuple(
            sorted(set(rate_failures.get(anchor, ()) + failures))
        )
    if (
        qualification_profile
        == _contracts.MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
    ):
        # Every active input shares one support. No mobility or Function-EEDF
        # interpolation is allowed to cross an anchor rejected by the direct
        # trajectory-rate consistency evidence.
        transport_eligible = active_eligible.difference(rate_failures)
        eedf_rate_eligible = set(transport_eligible)
    else:
        transport_eligible = set(full_transport_eligible)
        eedf_rate_eligible = aggregate_eligible.difference(rate_failures)
    for row in all_quality:
        anchor = float(row["E_over_N_Td"])
        full_reasons = quality_failure_reasons(row, "failure_reasons_json")
        reasons = sorted(
            set(active_failure_reasons[anchor] + list(rate_failures.get(anchor, ())))
        )
        row.update(
            {
                "failure_axes_json": json.dumps(
                    [axis.value for axis in failure_axes_for_reasons(full_reasons)],
                    separators=(",", ":"),
                ),
                "qualification_profile": qualification_profile,
                "active_closure_quality_passed": int(not reasons),
                "active_closure_failure_reasons_json": json.dumps(
                    reasons, separators=(",", ":")
                ),
                "active_closure_failure_axes_json": json.dumps(
                    [axis.value for axis in failure_axes_for_reasons(reasons)],
                    separators=(",", ":"),
                ),
            }
        )
    return _contracts.MonteCarloEvidence(
        all_quality,
        estimator_schema,
        aggregate_eligible,
        full_transport_eligible,
        active_failure_reasons,
        rate_failures,
        transport_eligible,
        eedf_rate_eligible,
    )


def _build_qualified_monte_carlo(
    connection: sqlite3.Connection,
    mixture_id: int,
    *,
    mc_sampling_plan: dict[str, object] | None,
    quality_thresholds: QualityThresholds,
    qualification_profile: str,
    evidence: _contracts.MonteCarloEvidence | None = None,
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[str, Any],
]:
    """Build a pure-MC table from independently qualified replicas.

    Monte Carlo transport, rates, and the EEDF all remain MC observations.
    No two-term coefficient, EEDF bin, or reaction-rate prior is substituted.
    Qualification is component-specific: a failed transport lag excludes that
    anchor from the six transport curves without discarding independently
    qualified EEDF/rate evidence at the same E/N.
    """

    evidence = evidence or _assess_monte_carlo(
        connection,
        mixture_id,
        quality_thresholds=quality_thresholds,
        qualification_profile=qualification_profile,
    )
    exclusions = _mc_exclusion_maps(evidence, qualification_profile)
    _require_mc_anchor_support(
        evidence,
        qualification_profile=qualification_profile,
        mc_sampling_plan=mc_sampling_plan,
        exclusions=exclusions,
    )
    cases, rates, eedf, quality, eedf_rate_cases = _load_qualified_mc_data(
        connection,
        mixture_id,
        evidence=evidence,
        qualification_profile=qualification_profile,
    )
    _derive_net_townsend_from_rates(cases, rates)
    weighted_growth, transport_policy, transport_convergence = (
        _mc_transport_descriptions(
            evidence.estimator_schema,
            qualification_profile,
        )
    )
    metadata = _mc_table_metadata(
        evidence=evidence,
        qualification_profile=qualification_profile,
        exclusions=exclusions,
        cases=cases,
        eedf_rate_cases=eedf_rate_cases,
        weighted_growth_transport=weighted_growth,
        transport_policy=transport_policy,
        transport_convergence=transport_convergence,
    )
    return cases, rates, eedf, quality, metadata


def _mc_exclusion_maps(
    evidence: _contracts.MonteCarloEvidence,
    qualification_profile: str,
) -> tuple[dict[str, list[str]], dict[str, list[str]], dict[str, list[str]]]:
    transport_excluded: dict[str, list[str]] = {}
    eedf_rate_excluded: dict[str, list[str]] = {}
    full_transport_excluded: dict[str, list[str]] = {}
    for row in evidence.quality:
        e_over_n = float(row["E_over_N_Td"])
        reasons = json.loads(str(row.get("failure_reasons_json", "[]")))
        aggregate_reasons = json.loads(
            str(row.get("aggregate_failure_reasons_json", "[]"))
        )
        if e_over_n not in evidence.full_transport_eligible:
            full_transport_excluded[f"{e_over_n:.17g}"] = list(dict.fromkeys(reasons))
        if e_over_n not in evidence.transport_eligible:
            selected = list(evidence.active_failure_reasons[e_over_n])
            selected.extend(evidence.rate_failures.get(e_over_n, ()))
            transport_excluded[f"{e_over_n:.17g}"] = list(dict.fromkeys(selected))
        if e_over_n not in evidence.eedf_rate_eligible:
            if (
                qualification_profile
                == _contracts.MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
            ):
                aggregate_reasons = list(evidence.active_failure_reasons[e_over_n])
            aggregate_reasons.extend(evidence.rate_failures.get(e_over_n, ()))
            eedf_rate_excluded[f"{e_over_n:.17g}"] = list(
                dict.fromkeys(aggregate_reasons)
            )
    return transport_excluded, eedf_rate_excluded, full_transport_excluded


def _require_mc_anchor_support(
    evidence: _contracts.MonteCarloEvidence,
    *,
    qualification_profile: str,
    mc_sampling_plan: dict[str, object] | None,
    exclusions: tuple[
        dict[str, list[str]],
        dict[str, list[str]],
        dict[str, list[str]],
    ],
) -> None:
    transport_excluded, eedf_rate_excluded, _ = exclusions
    if len(evidence.transport_eligible) < 4:
        qualification_label = (
            "active-profile-"
            if (
                qualification_profile
                == _contracts.MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
            )
            else "aggregate- and direct-transport-"
        )
        raise _contracts.MonteCarloQualificationError(
            "pure monte_carlo tables require at least four "
            f"{qualification_label}qualified E/N anchors; found "
            f"{len(evidence.transport_eligible)}. "
            "Apply the recorded failure-axis-specific sampling remedy instead "
            "of substituting a two_term prior. Excluded anchors: "
            + json.dumps(
                transport_excluded,
                sort_keys=True,
                separators=(",", ":"),
            )
        )
    if len(evidence.eedf_rate_eligible) < 4:
        raise _contracts.MonteCarloQualificationError(
            "pure monte_carlo tables require at least four aggregate- and "
            "EEDF/rate-consistency-qualified E/N anchors; found "
            f"{len(evidence.eedf_rate_eligible)}. Apply the recorded "
            "failure-axis-specific sampling remedy instead of substituting a "
            "two_term prior. "
            "Excluded anchors: "
            + json.dumps(
                eedf_rate_excluded,
                sort_keys=True,
                separators=(",", ":"),
            )
        )
    if not isinstance(mc_sampling_plan, dict) or not isinstance(
        mc_sampling_plan.get("entries"), list
    ):
        raise _contracts.TableBuildError(
            "pure monte_carlo tables require validated sampling-plan provenance"
        )
    entries = mc_sampling_plan["entries"]
    planned_fields = sorted(
        float(entry["e_over_n_Td"]) for entry in entries if isinstance(entry, dict)
    )
    quality_fields = sorted(float(row["E_over_N_Td"]) for row in evidence.quality)
    if (
        len(planned_fields) != len(entries)
        or len(planned_fields) != len(set(planned_fields))
        or planned_fields != quality_fields
    ):
        raise _contracts.TableBuildError(
            "Monte Carlo sampling-plan and aggregate quality anchors disagree"
        )
    qualified_fields = sorted(evidence.transport_eligible)
    first = planned_fields.index(qualified_fields[0])
    last = planned_fields.index(qualified_fields[-1])
    expected_qualified = planned_fields[first : last + 1]
    if qualified_fields != expected_qualified:
        internal_failures = sorted(
            set(expected_qualified).difference(evidence.transport_eligible)
        )
        raise _contracts.MonteCarloQualificationError(
            "Monte Carlo transport support crosses failed internal planned "
            "E/N anchors: " + ", ".join(f"{value:.17g}" for value in internal_failures)
        )
    if any(not reasons for reasons in transport_excluded.values()):
        raise _contracts.TableBuildError(
            "Monte Carlo excluded transport anchors require failure reasons"
        )


def _load_qualified_mc_data(
    connection: sqlite3.Connection,
    mixture_id: int,
    *,
    evidence: _contracts.MonteCarloEvidence,
    qualification_profile: str,
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
]:
    cases = _repository._load_cases(
        connection,
        _contracts.MC_SOLVER,
        mixture_id,
        allowed_e=evidence.transport_eligible,
    )
    _require_case_mean_energy_strictly_increasing(cases)
    eedf_rate_cases = _repository._load_cases(
        connection,
        _contracts.MC_SOLVER,
        mixture_id,
        allowed_e=evidence.eedf_rate_eligible,
    )
    _require_case_mean_energy_strictly_increasing(eedf_rate_cases)
    rates = _repository._load_rates(
        connection,
        _contracts.MC_SOLVER,
        mixture_id,
        eedf_rate_cases,
        allowed_e=evidence.eedf_rate_eligible,
    )
    eedf = _repository._load_eedf(
        connection,
        _contracts.MC_SOLVER,
        mixture_id,
        eedf_rate_cases,
        allowed_e=evidence.eedf_rate_eligible,
    )
    quality: list[dict[str, Any]] = []
    for row in evidence.quality:
        e_over_n = float(row["E_over_N_Td"])
        if e_over_n not in evidence.transport_eligible:
            continue
        active_reasons = list(evidence.active_failure_reasons[e_over_n])
        active_reasons.extend(evidence.rate_failures.get(e_over_n, ()))
        quality.append(
            {
                **dict(row),
                "qualification_profile": qualification_profile,
                "active_closure_quality_passed": int(not active_reasons),
                "active_closure_failure_reasons_json": json.dumps(
                    sorted(set(active_reasons)), separators=(",", ":")
                ),
            }
        )
    return cases, rates, eedf, quality, eedf_rate_cases


def _mc_transport_descriptions(
    estimator_schema: str | None,
    qualification_profile: str,
) -> tuple[bool, str, str]:
    weighted_growth = (
        estimator_schema == _mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
    )
    policy = (
        "direct_mc_synchronized_weighted_growth_flux_mobility_position_"
        "velocity_covariance_diffusion_L_T_residence_energy_mobility_and_"
        "restricted_density_packet_energy_diffusion_L_T"
        if weighted_growth
        else (
            "direct_mc_production_displacement_mobility_block_helfand_"
            "particle_diffusion_L_T_residence_energy_mobility_and_"
            "restricted_density_packet_energy_diffusion_L_T"
        )
    )
    convergence = (
        (
            "restricted_lmea_mean_energy_paired_ci95_stationarity_mobility_"
            "paired_systematic_log_drift_and_independent_replica_student_t_"
            "95_coefficient_precision"
            if (
                qualification_profile
                == _contracts.MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
            )
            else (
                "paired_log_ratio_common_time_flux_moment_stationarity_with_"
                "minimum_per_lag_lineage_ESS_across_all_complete_blocks_and_"
                "case_vs_transport_mean_energy_across_independent_replicas"
            )
        )
        if weighted_growth
        else (
            "paired_log_ratio_mean_energy_origin_stationarity_at_10pct_"
            "fixed_lag_diffusion_plateau_at_25pct_and_case_vs_transport_"
            "mean_energy_at_10pct_across_independent_replicas"
        )
    )
    return weighted_growth, policy, convergence


def _mc_table_metadata(
    *,
    evidence: _contracts.MonteCarloEvidence,
    qualification_profile: str,
    exclusions: tuple[
        dict[str, list[str]],
        dict[str, list[str]],
        dict[str, list[str]],
    ],
    cases: list[dict[str, Any]],
    eedf_rate_cases: list[dict[str, Any]],
    weighted_growth_transport: bool,
    transport_policy: str,
    transport_convergence: str,
) -> dict[str, Any]:
    transport_excluded, eedf_rate_excluded, full_transport_excluded = exclusions
    active_mean_energies = [
        _repository._required_float(case, "mean_energy_eV") for case in cases
    ]
    adjacent_mean_ratios = [
        right / left
        for left, right in zip(active_mean_energies, active_mean_energies[1:])
    ]
    adjacent_mean_gaps = [
        right - left
        for left, right in zip(active_mean_energies, active_mean_energies[1:])
    ]
    return {
        "source": _contracts.MC_SOLVER,
        "postprocess": "none",
        "selection": "component_specific_fail_closed_mc_anchors",
        "qualification_profile": qualification_profile,
        # These are the quantities covered by this statistical qualification.
        # A downstream model adapter independently declares which subset it binds.
        "qualification_outputs": (
            [
                "function_eedf",
                "reduced_mobility",
                "direct_elastic_energy_loss",
            ]
            if qualification_profile
            == _contracts.MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
            else [
                "full_particle_transport",
                "full_restricted_energy_transport",
                "direct_reaction_rates",
                "eedf",
            ]
        ),
        "additional_mc_evidence": (
            [
                "particle_diffusion_L_T",
                "energy_mobility",
                "restricted_energy_diffusion_L_T",
                "direct_reaction_rates",
            ]
            if qualification_profile
            == _contracts.MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
            else []
        ),
        "transport_eligible_anchor_points": len(evidence.transport_eligible),
        "transport_eligible_E_over_N_Td": sorted(evidence.transport_eligible),
        "transport_excluded_E_over_N_Td": transport_excluded,
        "full_transport_eligible_anchor_points": len(evidence.full_transport_eligible),
        "full_transport_eligible_E_over_N_Td": sorted(evidence.full_transport_eligible),
        "full_transport_excluded_E_over_N_Td": full_transport_excluded,
        "active_mean_energy_anchors_eV": active_mean_energies,
        "maximum_adjacent_mean_energy_ratio": max(adjacent_mean_ratios),
        "maximum_adjacent_mean_energy_gap_eV": max(adjacent_mean_gaps),
        "eedf_rate_eligible_anchor_points": len(evidence.eedf_rate_eligible),
        "eedf_rate_eligible_E_over_N_Td": sorted(evidence.eedf_rate_eligible),
        "eedf_rate_excluded_E_over_N_Td": eedf_rate_excluded,
        "component_mean_energy_range_eV": {
            "transport": [
                _repository._required_float(cases[0], "mean_energy_eV"),
                _repository._required_float(cases[-1], "mean_energy_eV"),
            ],
            "eedf_rate": [
                _repository._required_float(eedf_rate_cases[0], "mean_energy_eV"),
                _repository._required_float(eedf_rate_cases[-1], "mean_energy_eV"),
            ],
        },
        "transport_policy": transport_policy,
        "transport_scope": {
            "standard_energy_density_gradient_closure": (
                "restricted_direct_density_packet_energy_flux_correlation"
            ),
            "full_cross_gradient_response_identified": False,
        },
        "eedf_policy": (
            "independent_replica_mc_histogram_then_nonnegative_local_C1_"
            "reconstruction_and_moment_preserving_COMSOL_projection"
        ),
        "reaction_rate_policy": (
            "direct_mc_trajectory_rates_as_independent_function_eedf_tail_"
            "evidence_not_active_COMSOL_input"
            if qualification_profile
            == _contracts.MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
            else (
                "direct_mc_trajectory_rate_coefficients_as_active_input;_"
                "raw_rate_evidence_as_qualification"
            )
        ),
        "smoothing": {
            "transport_and_rates": "none",
            "eedf": (
                "shape_preserving_local_reconstruction_only; no_external_"
                "distribution_prior"
            ),
        },
        "quality_policy": {
            "failed_points": (
                "active_closure_failures_excluded_only_outside_contiguous_"
                "support;inactive_MC_failures_retained_not_repaired"
                if qualification_profile
                == _contracts.MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
                else (
                    "excluded_per_component_only_outside_contiguous_transport_"
                    "support_not_repaired"
                )
            ),
            "transport_convergence": transport_convergence,
            "origin_stationarity_relative_ci95_bound_limit": (
                MC_TRANSPORT_STATE_RELATIVE_TOLERANCE
            ),
            **(
                {
                    "restricted_lmea_mean_energy_stationarity_"
                    "relative_ci95_bound_limit": (
                        MC_TRANSPORT_STATE_RELATIVE_TOLERANCE
                    ),
                    "restricted_lmea_mobility_systematic_absolute_"
                    "log_drift_limit": MC_TRANSPORT_STATE_RELATIVE_TOLERANCE,
                    "restricted_lmea_mobility_independent_replica_"
                    "relative_ci95_half_width_limit": (
                        MC_TRANSPORT_STATE_RELATIVE_TOLERANCE
                    ),
                }
                if qualification_profile
                == _contracts.MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
                else {}
            ),
            **(
                {
                    "diffusion_stationarity_relative_ci95_bound_limit": (
                        MC_TRANSPORT_DIFFUSION_LAG_RELATIVE_TOLERANCE
                    ),
                    "lineage_qualification_lag_rule": "[L/2,L]",
                    "required_effective_lineage_count": (
                        _mc_evidence.WEIGHTED_MC_TRANSPORT_MIN_EFFECTIVE_LINEAGE_COUNT
                    ),
                    "required_effective_lineage_fraction": (
                        _mc_evidence.WEIGHTED_MC_TRANSPORT_MIN_EFFECTIVE_LINEAGE_FRACTION
                    ),
                }
                if weighted_growth_transport
                else {
                    "diffusion_lag_relative_ci95_bound_limit": (
                        MC_TRANSPORT_DIFFUSION_LAG_RELATIVE_TOLERANCE
                    )
                }
            ),
            "transport_mean_energy_relative_ci95_bound_limit": (
                MC_TRANSPORT_STATE_RELATIVE_TOLERANCE
            ),
            "eedf_rate_consistency": (
                "direct_trajectory_rate_vs_same_MC_histogram_convolution"
            ),
            "eedf_rate_relative_ci95_bound_limit": (
                MC_EEDF_RATE_ENSEMBLE_RELATIVE_TOLERANCE
            ),
            "unresolved_rate_tail": (
                "defer_to_rate_RSE_and_pooled_zero_event_upper_bound"
            ),
        },
        "net_ionization_policy": "derive_from_mc_ionization_and_attachment_rates",
    }


def _require_case_mean_energy_strictly_increasing(
    cases: list[dict[str, Any]],
) -> None:
    values = [_repository._required_float(case, "mean_energy_eV") for case in cases]
    if any(right <= left for left, right in zip(values, values[1:])):
        raise _contracts.TableBuildError(
            "qualified MC mean energy must be strictly increasing with E/N; "
            "points are never repaired by monotonicity postprocessing"
        )


def _derive_net_townsend_from_rates(
    cases: list[dict[str, Any]],
    rates: list[dict[str, Any]],
) -> None:
    by_e: dict[float, list[dict[str, Any]]] = {}
    for rate in rates:
        by_e.setdefault(float(rate["E_over_N_Td"]), []).append(rate)
    for case in cases:
        rows = by_e.get(float(case["E_over_N_Td"]), [])
        ionization = 0.0
        attachment = 0.0
        saw_component = False
        for row in rows:
            process_type = str(row["process_type"]).lower()
            value = float(row["mixture_weighted_reduced_townsend_m2"])
            if "ionization" in process_type:
                ionization += value
                saw_component = True
            elif "attachment" in process_type:
                attachment += value
                saw_component = True
        if saw_component:
            case["effective_townsend_m2"] = ionization - attachment


def _load_rate_evidence(
    connection: sqlite3.Connection,
    mixture_id: int,
    *,
    allowed_e: set[float],
) -> list[dict[str, Any]]:
    """Load untouched MC replica statistics and pooled event exposure."""

    raw_cases = _repository._load_cases(
        connection,
        _contracts.MC_SOLVER,
        mixture_id,
        allowed_e=allowed_e,
    )
    mean_by_e = {
        float(case["E_over_N_Td"]): _repository._required_float(case, "mean_energy_eV")
        for case in raw_cases
    }
    event_diagnostics = _load_mc_rate_event_diagnostics(
        connection,
        mixture_id,
        allowed_e=allowed_e,
    )
    rows = connection.execute(
        """
        SELECT e_over_n_Td, species, process, process_type, threshold_eV,
               target_species_fraction, mean, standard_error,
               relative_standard_error, ci95_low, ci95_high,
               ci95_critical_value, estimate_status, valid_replicates,
               uncertainty_available
        FROM aggregate_scalars
        WHERE solver = ? AND mixture_id = ? AND scalar_group = 'rate'
          AND scalar_name = ?
        ORDER BY e_over_n_Td, species, process, process_type, threshold_key
        """,
        (_contracts.MC_SOLVER, mixture_id, RATE_SCALAR),
    ).fetchall()
    evidence: list[dict[str, Any]] = []
    for row in rows:
        e_over_n = float(row["e_over_n_Td"])
        if e_over_n not in allowed_e:
            continue
        if e_over_n not in mean_by_e:
            raise _contracts.TableBuildError(
                f"MC rate evidence lacks mean energy at {e_over_n:.17g} Td"
            )
        key = _repository._rate_key(dict(row))
        target_fraction = _repository._float_or_none(row["target_species_fraction"])
        if target_fraction is None:
            target_fraction = 1.0
        event_count, exposure, upper, event_status = _pooled_zero_event_evidence(
            event_diagnostics.get((e_over_n, key), []),
            expected_replicates=int(row["valid_replicates"]),
            expected_target_fraction=target_fraction,
        )
        evidence.append(
            {
                **_repository._en_fields(e_over_n),
                "mean_energy_eV": mean_by_e[e_over_n],
                "species": str(row["species"]),
                "process": str(row["process"]),
                "process_type": str(row["process_type"]),
                "threshold_eV": _repository._float_or_none(row["threshold_eV"]),
                "target_species_fraction": target_fraction,
                "rate_coefficient_mean_m3_s": _repository._float_or_none(row["mean"]),
                "rate_coefficient_standard_error_m3_s": _repository._float_or_none(
                    row["standard_error"]
                ),
                "rate_coefficient_relative_standard_error": _repository._float_or_none(
                    row["relative_standard_error"]
                ),
                "rate_coefficient_ci95_low_m3_s": _repository._float_or_none(
                    row["ci95_low"]
                ),
                "rate_coefficient_ci95_high_m3_s": _repository._float_or_none(
                    row["ci95_high"]
                ),
                "ci95_critical_value": _repository._float_or_none(
                    row["ci95_critical_value"]
                ),
                "estimate_status": str(row["estimate_status"]),
                "valid_replicates": int(row["valid_replicates"]),
                "uncertainty_available": int(row["uncertainty_available"]),
                "pooled_event_count": event_count,
                "pooled_target_exposure_s_m3": exposure,
                "pooled_all_zero_upper_95_m3_s": upper,
                "pooled_zero_event_status": event_status,
            }
        )
    return evidence


def _validate_mc_censored_rate_relevance(
    rates: list[dict[str, Any]],
    evidence: list[dict[str, Any]],
    *,
    minimum_fraction: float,
    raise_on_failure: bool = True,
) -> dict[str, Any]:
    """Accept a nominal zero only when its exact pooled upper bound is negligible."""

    fraction = float(minimum_fraction)
    peaks: dict[tuple[str, str, str, str], float] = {}
    canonical: dict[tuple[float, tuple[str, str, str, str]], float] = {}
    for row in rates:
        key = _repository._rate_key(row)
        field = _repository._required_float(row, "E_over_N_Td")
        value = _repository._required_float(row, "rate_coefficient_m3_s")
        anchor_key = (field, key)
        if anchor_key in canonical:
            raise _contracts.TableBuildError(
                "Monte Carlo canonical rates contain a duplicate process anchor"
            )
        canonical[anchor_key] = value
        peaks[key] = max(value, peaks.get(key, 0.0))

    censored: list[dict[str, Any]] = []
    failed: list[dict[str, Any]] = []
    for row in evidence:
        key = _repository._rate_key(row)
        field = _repository._required_float(row, "E_over_N_Td")
        anchor_key = (field, key)
        if anchor_key not in canonical:
            raise _contracts.TableBuildError(
                "Monte Carlo rate evidence does not match the canonical rate grid"
            )
        canonical_rate = canonical[anchor_key]
        evidence_rate = _repository._required_float(row, "rate_coefficient_mean_m3_s")
        is_censored = str(row.get("estimate_status", "")) == "censored_all_zero"
        if canonical_rate != 0.0 and not is_censored:
            continue
        count = _repository._required_float(row, "pooled_event_count")
        exposure = _repository._required_float(row, "pooled_target_exposure_s_m3")
        upper = _repository._required_float(row, "pooled_all_zero_upper_95_m3_s")
        recomputed = (
            -math.log(1.0 - _contracts.ZERO_EVENT_CONFIDENCE) / exposure
            if exposure > 0.0
            else math.nan
        )
        if (
            canonical_rate != 0.0
            or evidence_rate != 0.0
            or not is_censored
            or str(row.get("pooled_zero_event_status", "")) != "all_zero_upper_95"
            or str(row.get("uncertainty_available", "")).strip().lower()
            not in {"0", "false"}
            or count != 0.0
            or exposure <= 0.0
            or upper <= 0.0
            or not math.isfinite(recomputed)
            or not math.isclose(upper, recomputed, rel_tol=1.0e-12, abs_tol=0.0)
        ):
            raise _contracts.TableBuildError(
                "Monte Carlo censored rate lacks self-consistent pooled "
                f"zero-event evidence at {field:.17g} Td"
            )
        peak = peaks.get(key, 0.0)
        limit = fraction * peak
        passed = peak > 0.0 and upper < limit
        result = {
            "E_over_N_Td": field,
            "species": key[0],
            "process": key[1],
            "process_type": key[2],
            "process_peak_m3_s": peak,
            "pooled_upper_95_m3_s": upper,
            "required_upper_limit_m3_s": limit,
            "upper_to_process_peak_fraction": (upper / peak if peak > 0.0 else None),
            "passed": passed,
        }
        censored.append(result)
        if not passed:
            failed.append(result)
    if failed and raise_on_failure:
        raise _contracts.MonteCarloQualificationError(
            "Monte Carlo censored-rate pooled 95% upper bound is not strictly "
            "below the configured process-peak fraction: "
            + json.dumps(failed, sort_keys=True, separators=(",", ":"))
        )
    return {
        "status": "qualified" if not failed else "unqualified",
        "criterion": (
            "pooled_upper_95_strictly_below_required_fraction_of_process_peak"
        ),
        "required_rate_min_process_peak_fraction": fraction,
        "zero_event_confidence": _contracts.ZERO_EVENT_CONFIDENCE,
        "censored_anchor_count": len(censored),
        "failed_anchor_count": len(failed),
        "censored_anchors": censored,
        "artificial_rate_floor": False,
    }


def _mc_censored_rate_relevance_failures(
    connection: sqlite3.Connection,
    mixture_id: int,
    *,
    allowed_e: set[float],
    minimum_fraction: float,
) -> dict[float, tuple[str, ...]]:
    """Expose unresolved zero-event tails to the bounded retry policy."""

    if not allowed_e:
        return {}
    cases = _repository._load_cases(
        connection,
        _contracts.MC_SOLVER,
        mixture_id,
        allowed_e=allowed_e,
    )
    rates = _repository._load_rates(
        connection,
        _contracts.MC_SOLVER,
        mixture_id,
        cases,
        allowed_e=allowed_e,
    )
    evidence = _load_rate_evidence(
        connection,
        mixture_id,
        allowed_e=allowed_e,
    )
    assessment = _validate_mc_censored_rate_relevance(
        rates,
        evidence,
        minimum_fraction=minimum_fraction,
        raise_on_failure=False,
    )
    failures: dict[float, tuple[str, ...]] = {}
    for row in assessment["censored_anchors"]:
        if not bool(row["passed"]):
            failures[float(row["E_over_N_Td"])] = (
                "mc_censored_rate_relevance_failed",
            )
    return failures


def _load_mc_rate_event_diagnostics(
    connection: sqlite3.Connection,
    mixture_id: int,
    *,
    allowed_e: set[float],
) -> dict[tuple[float, tuple[str, str, str, str]], list[dict[str, Any]]]:
    grouped: dict[
        tuple[float, tuple[str, str, str, str]],
        list[dict[str, Any]],
    ] = {}
    rows = connection.execute(
        """
        SELECT e_over_n_Td, gas_number_density_m3, diagnostics_json
        FROM cases
        WHERE solver = ? AND mixture_id = ?
        ORDER BY e_over_n_Td, replicate
        """,
        (_contracts.MC_SOLVER, mixture_id),
    ).fetchall()
    for row in rows:
        e_over_n = float(row["e_over_n_Td"])
        if e_over_n not in allowed_e:
            continue
        try:
            diagnostics = json.loads(str(row["diagnostics_json"]))
        except json.JSONDecodeError:
            continue
        section = diagnostics.get("internal_monte_carlo_reaction_rates")
        rate_rows = section.get("rates") if isinstance(section, dict) else None
        if not isinstance(rate_rows, list):
            continue
        density = _repository._float_or_none(row["gas_number_density_m3"])
        for item in rate_rows:
            if not isinstance(item, dict):
                continue
            try:
                key = _repository._rate_key(item)
            except (KeyError, TypeError, ValueError):
                continue
            count = _nonnegative_int_or_none(item.get("event_count"))
            residence_time = _repository._finite_nonnegative(
                item.get("event_observation_residence_time_s")
            )
            fraction = _repository._finite_nonnegative(
                item.get("target_species_fraction")
            )
            confidence = _repository._float_or_none(item.get("zero_event_confidence"))
            valid = (
                density is not None
                and density > 0.0
                and residence_time is not None
                and residence_time > 0.0
                and fraction is not None
                and fraction > 0.0
                and count is not None
                and confidence is not None
                and math.isclose(confidence, _contracts.ZERO_EVENT_CONFIDENCE)
            )
            grouped.setdefault((e_over_n, key), []).append(
                {
                    "valid": valid,
                    "sampled": (
                        item.get("event_sampling_enabled") is True
                        and item.get("event_observation_status")
                        in {
                            "observed",
                            "unobserved_zero_events_upper_bound_only",
                        }
                    ),
                    "event_count": count,
                    "target_species_fraction": fraction,
                    "target_exposure_s_m3": (
                        density * fraction * residence_time if valid else None
                    ),
                }
            )
    return grouped


def _mc_eedf_rate_consistency_failures(
    connection: sqlite3.Connection,
    mixture_id: int,
    *,
    allowed_e: set[float],
    minimum_process_peak_fraction: float = 0.0,
) -> dict[float, tuple[str, ...]]:
    """Check that each MC EEDF reproduces its direct trajectory rates.

    Both estimates come from the same trajectory.  This is therefore a
    representation check on the histogram EEDF, not a second kinetic model
    and not an additional Monte Carlo run.
    """

    rows = connection.execute(
        """
        SELECT e_over_n_Td, diagnostics_json
        FROM cases
        WHERE solver = ? AND mixture_id = ?
        ORDER BY e_over_n_Td, replicate
        """,
        (_contracts.MC_SOLVER, mixture_id),
    ).fetchall()
    expected_replicates: dict[float, int] = {}
    grouped: dict[
        tuple[float, tuple[str, str, str, str]],
        list[tuple[float, float]],
    ] = {}
    invalid: set[float] = set()
    for row in rows:
        e_over_n = float(row["e_over_n_Td"])
        if e_over_n not in allowed_e:
            continue
        expected_replicates[e_over_n] = expected_replicates.get(e_over_n, 0) + 1
        try:
            diagnostics = json.loads(str(row["diagnostics_json"]))
        except json.JSONDecodeError:
            invalid.add(e_over_n)
            continue
        section = diagnostics.get("internal_monte_carlo_reaction_rates")
        rate_rows = section.get("rates") if isinstance(section, dict) else None
        if (
            not isinstance(section, dict)
            or section.get("estimator") != "trajectory_time_average_sigma_v"
            or not isinstance(rate_rows, list)
            or not rate_rows
        ):
            invalid.add(e_over_n)
            continue
        seen: set[tuple[str, str, str, str]] = set()
        for item in rate_rows:
            if not isinstance(item, dict):
                invalid.add(e_over_n)
                continue
            try:
                key = _repository._rate_key(item)
            except (KeyError, TypeError, ValueError):
                invalid.add(e_over_n)
                continue
            direct = _repository._finite_nonnegative(item.get("rate_coefficient_m3_s"))
            histogram = _repository._finite_nonnegative(
                item.get("histogram_convolution_rate_coefficient_m3_s")
            )
            if key in seen or direct is None or histogram is None:
                invalid.add(e_over_n)
                continue
            seen.add(key)
            grouped.setdefault((e_over_n, key), []).append((direct, histogram))

    minimum_fraction = float(minimum_process_peak_fraction)
    if (
        not math.isfinite(minimum_fraction)
        or minimum_fraction < 0.0
        or minimum_fraction > 1.0
    ):
        raise ValueError("MC EEDF-rate peak fraction must be in [0, 1]")
    process_peaks: dict[tuple[str, str, str, str], float] = {}
    for (field, key), values in grouped.items():
        if field not in allowed_e or not values:
            continue
        process_peaks[key] = max(
            process_peaks.get(key, 0.0),
            float(sum(value[0] for value in values) / len(values)),
        )

    failures: dict[float, tuple[str, ...]] = {}
    for e_over_n in allowed_e:
        expected = expected_replicates.get(e_over_n, 0)
        evidence = [
            (key, values)
            for (field, key), values in grouped.items()
            if field == e_over_n
        ]
        if e_over_n in invalid or expected < MC_TRANSPORT_MIN_ENSEMBLE_REPLICATES:
            failures[e_over_n] = ("mc_eedf_rate_evidence_incomplete",)
            continue
        if not evidence or any(len(values) != expected for _key, values in evidence):
            failures[e_over_n] = ("mc_eedf_rate_evidence_incomplete",)
            continue
        consistent = True
        for key, values in evidence:
            direct = [value[0] for value in values]
            histogram = [value[1] for value in values]
            field_mean = float(sum(direct) / len(direct))
            if field_mean < minimum_fraction * process_peaks.get(key, 0.0):
                continue
            if all(left == 0.0 and right == 0.0 for left, right in values):
                continue
            if any(left <= 0.0 or right <= 0.0 for left, right in values):
                # A mixed zero/nonzero ensemble is an unresolved reaction
                # tail. Its qualification belongs to the rate-RSE and pooled
                # zero-event upper-bound policy, not this shape check.
                continue
            bound = paired_log_ratio_ci95_bound(histogram, direct)
            if bound is None or bound > MC_EEDF_RATE_ENSEMBLE_RELATIVE_TOLERANCE:
                consistent = False
                break
        if not consistent:
            failures[e_over_n] = ("mc_eedf_rate_consistency_failed",)
    return failures


def _pooled_zero_event_evidence(
    rows: list[dict[str, Any]],
    *,
    expected_replicates: int,
    expected_target_fraction: float,
) -> tuple[int | None, float | None, float | None, str]:
    if not rows:
        return None, None, None, "diagnostics_unavailable"
    if len(rows) != expected_replicates:
        return None, None, None, "diagnostics_incomplete"
    if not all(bool(row["valid"]) for row in rows):
        return None, None, None, "diagnostics_invalid"
    if not all(
        math.isclose(
            float(row["target_species_fraction"]),
            expected_target_fraction,
            rel_tol=1.0e-12,
            abs_tol=1.0e-15,
        )
        for row in rows
    ):
        return None, None, None, "diagnostics_invalid"
    if not all(bool(row["sampled"]) for row in rows):
        return None, None, None, "not_sampled_by_trajectory_model"
    event_count = sum(int(row["event_count"]) for row in rows)
    exposure = sum(float(row["target_exposure_s_m3"]) for row in rows)
    if event_count > 0:
        return event_count, exposure, None, "observed_events"
    upper = -math.log(1.0 - _contracts.ZERO_EVENT_CONFIDENCE) / exposure
    return event_count, exposure, upper, "all_zero_upper_95"


def _nonnegative_int_or_none(value: object) -> int | None:
    number = _repository._float_or_none(value)
    if number is None or number < 0.0 or not number.is_integer():
        return None
    return int(number)
