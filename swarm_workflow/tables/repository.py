"""SQLite-backed reads and provenance validation for workflow tables."""

from __future__ import annotations

import json
import math
import sqlite3
from typing import Any

from electron_swarm.solvers.two_term.transport import (
    TWO_TERM_TEMPORAL_GROWTH_EFFECTIVE_MOMENTUM_MODEL,
    TWO_TERM_TEMPORAL_GROWTH_TRANSPORT_SCHEMA_VERSION,
)
import electron_swarm.solvers.monte_carlo.evidence as _mc_evidence
from electron_swarm.solvers.monte_carlo.source_identity import (
    monte_carlo_source_sha256,
)

from . import contracts as _contracts
from ..campaign.aggregate import RATE_SCALAR, stable_mc_seed
from ..quality.monte_carlo.direct_transport import summarize_mc_transport_diagnostics
from ..quality.monte_carlo.weighted_transport import (
    summarize_weighted_mc_transport_diagnostics,
)
from ..quality.policy import parse_quality_thresholds
from ..campaign.repository import read_metadata, table_exists
from ..campaign.store import WorkflowSchemaError, mc_sampling_plan_provenance


def _field_source_policy(
    connection: sqlite3.Connection,
    *,
    solver: str,
    mixture_id: int,
    require_mc_solver_physics: bool = False,
) -> dict[str, Any]:
    rows = connection.execute(
        """
        SELECT transport_definition, metadata_json, diagnostics_json
        FROM cases
        WHERE solver = ? AND mixture_id = ?
        ORDER BY e_over_n_Td, replicate
        """,
        (solver, mixture_id),
    ).fetchall()
    payloads = [json.loads(str(row["metadata_json"])) for row in rows]
    treatments = {
        str(payload.get("rf_field_treatment", "none")) for payload in payloads
    }
    frequencies = {
        float(payload["rf_frequency_Hz"])
        for payload in payloads
        if payload.get("rf_frequency_Hz") is not None
    }
    amplitudes = {
        str(payload.get("rf_amplitude_definition", "none")) for payload in payloads
    }
    transport_definitions = {
        str(row["transport_definition"])
        for row in rows
        if row["transport_definition"] not in {None, ""}
    }
    if (
        len(treatments) > 1
        or len(frequencies) > 1
        or len(amplitudes) > 1
        or len(transport_definitions) > 1
    ):
        raise _contracts.TableBuildError(
            "workflow cases have inconsistent field or transport metadata"
        )
    treatment = next(iter(treatments), "none")
    policy: dict[str, Any] = {
        "field_type": "time_dependent" if treatment != "none" else "dc",
        "rf_field_treatment": treatment,
        "rf_frequency_Hz": next(iter(frequencies), None),
        "rf_amplitude_definition": next(iter(amplitudes), "none"),
        "transport_definition": next(iter(transport_definitions), "unspecified"),
    }
    if solver == _contracts.MC_SOLVER and require_mc_solver_physics:
        physics_keys = (
            "angular_model",
            "angular_moment_source",
            "exact_dcs_based",
            "ordinary_integral_xs_closure",
            "ionization_source_model",
            "ionization_source_treatment",
            "electron_electron_treatment",
            "electron_electron_transport_stale",
            "magnetic_field_treatment",
            "tail_refinement_treatment",
        )
        physics_values = {
            key: {payload.get(key) for payload in payloads} for key in physics_keys
        }
        if any(
            len(values) != 1 or None in values for values in physics_values.values()
        ):
            raise _contracts.TableBuildError(
                "Monte Carlo cases have inconsistent solver-physics metadata"
            )
        populations: set[str] = set()
        extrapolations: set[str] = set()
        for row in rows:
            try:
                diagnostics = json.loads(str(row["diagnostics_json"]))
                transport = diagnostics["internal_monte_carlo_transport"]
                provenance = transport["mc_run_provenance"]
                populations.add(str(provenance["population_model"]))
                extrapolations.add(str(provenance["high_energy_extrapolation"]))
            except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
                raise _contracts.TableBuildError(
                    "Monte Carlo cases lack solver-physics run provenance"
                ) from exc
        if len(populations) != 1 or len(extrapolations) != 1:
            raise _contracts.TableBuildError(
                "Monte Carlo cases have inconsistent population or tail models"
            )
        policy["solver_physics"] = {
            **{key: next(iter(values)) for key, values in physics_values.items()},
            "population_model": next(iter(populations)),
            "high_energy_extrapolation": next(iter(extrapolations)),
        }
    return policy


def _mixture_ids(connection: sqlite3.Connection, source: str) -> list[int]:
    solvers = (source,)
    placeholders = ", ".join("?" for _solver in solvers)
    rows = connection.execute(
        f"""
        SELECT DISTINCT mixture_id
        FROM aggregate_scalars
        WHERE solver IN ({placeholders})
        ORDER BY mixture_id
        """,
        solvers,
    ).fetchall()
    return [int(row["mixture_id"]) for row in rows]


def _mixture_rows(
    connection: sqlite3.Connection,
    mixture_id: int,
) -> list[dict[str, str | float]]:
    rows = connection.execute(
        """
        SELECT species, fraction, mass_amu
        FROM mixture_species
        WHERE mixture_id = ?
        ORDER BY species
        """,
        (mixture_id,),
    ).fetchall()
    return [
        {
            "species": str(row["species"]),
            "fraction": float(row["fraction"]),
            "mass_amu": float(row["mass_amu"]),
        }
        for row in rows
    ]


def _load_cases(
    connection: sqlite3.Connection,
    solver: str,
    mixture_id: int,
    *,
    allowed_e: set[float] | None = None,
) -> list[dict[str, Any]]:
    rows = connection.execute(
        """
        SELECT e_over_n_Td, scalar_name, mean, relative_standard_error
        FROM aggregate_scalars
        WHERE solver = ? AND mixture_id = ? AND scalar_group = 'case'
        ORDER BY e_over_n_Td, scalar_name
        """,
        (solver, mixture_id),
    ).fetchall()
    by_e: dict[float, dict[str, Any]] = {}
    for row in rows:
        e_over_n = float(row["e_over_n_Td"])
        if allowed_e is not None and e_over_n not in allowed_e:
            continue
        case = by_e.setdefault(e_over_n, _en_fields(e_over_n))
        case[str(row["scalar_name"])] = _float_or_none(row["mean"])
        case[f"{row['scalar_name']}_rse"] = _float_or_none(
            row["relative_standard_error"]
        )
    return [by_e[e_over_n] for e_over_n in sorted(by_e)]


def _load_two_term_temporal_growth_transport_contract(
    connection: sqlite3.Connection,
    mixture_id: int,
    cases: list[dict[str, Any]],
) -> dict[str, Any] | None:
    """Bind the PT growth eigenvalue to every exported two-term anchor.

    Older ``ignore``/time-dependent two-term runs do not carry this contract
    and remain exportable for their own profiles.  If any anchor claims the PT
    correction, however, every exported anchor must provide the same complete
    converged-eigenvalue provenance; partial evidence is rejected.
    """

    expected_fields = {float(case["E_over_N_Td"]) for case in cases}
    rows = connection.execute(
        """
        SELECT e_over_n_Td, replicate, gas_number_density_m3,
               diagnostics_json
        FROM cases
        WHERE solver = ? AND mixture_id = ?
        ORDER BY e_over_n_Td, replicate
        """,
        (_contracts.TWO_TERM_SOLVER, mixture_id),
    ).fetchall()
    grouped: dict[float, list[tuple[float, float]]] = {}
    applied_flags: list[bool] = []
    invalid_claims: list[str] = []
    for row in rows:
        field = float(row["e_over_n_Td"])
        if field not in expected_fields:
            continue
        try:
            diagnostics = json.loads(str(row["diagnostics_json"]))
        except (TypeError, json.JSONDecodeError):
            diagnostics = None
        two_term = (
            diagnostics.get("two_term") if isinstance(diagnostics, dict) else None
        )
        correction = (
            two_term.get("temporal_growth_momentum_correction")
            if isinstance(two_term, dict)
            else None
        )
        applied = bool(
            isinstance(correction, dict) and correction.get("applied") is True
        )
        applied_flags.append(applied)
        if not applied:
            continue
        density = _float_or_none(row["gas_number_density_m3"])
        growth = _float_or_none(correction.get("growth_frequency_s-1"))
        reported_growth = _float_or_none(
            two_term.get("growth_frequency_s-1") if isinstance(two_term, dict) else None
        )
        valid = bool(
            isinstance(two_term, dict)
            and two_term.get("converged") is True
            and correction.get("model")
            == TWO_TERM_TEMPORAL_GROWTH_EFFECTIVE_MOMENTUM_MODEL
            and correction.get("frequency_source")
            == "converged_temporal_growth_eigenvalue"
            and density is not None
            and density > 0.0
            and growth is not None
            and reported_growth is not None
            and math.isclose(
                growth,
                reported_growth,
                rel_tol=1.0e-12,
                abs_tol=1.0e-12,
            )
        )
        if not valid:
            invalid_claims.append(f"{field:.17g}Td/replicate={int(row['replicate'])}")
            continue
        assert density is not None and growth is not None
        grouped.setdefault(field, []).append((density, growth))

    if not applied_flags or not any(applied_flags):
        return None
    if not all(applied_flags) or invalid_claims or set(grouped) != expected_fields:
        details = invalid_claims or [
            f"{field:.17g}Td" for field in sorted(expected_fields.difference(grouped))
        ]
        raise _contracts.TableBuildError(
            "two-term temporal-growth transport provenance is partial or "
            "invalid: " + ", ".join(details)
        )

    by_field = {float(case["E_over_N_Td"]): case for case in cases}
    growth_values: list[float] = []
    density_values: list[float] = []
    for field in sorted(expected_fields):
        evidence = grouped[field]
        density0, growth0 = evidence[0]
        if any(
            not math.isclose(density, density0, rel_tol=1.0e-12, abs_tol=0.0)
            or not math.isclose(
                growth,
                growth0,
                rel_tol=1.0e-10,
                abs_tol=1.0e-12,
            )
            for density, growth in evidence[1:]
        ):
            raise _contracts.TableBuildError(
                "two-term temporal-growth transport replicas disagree at "
                f"{field:.17g} Td"
            )
        case = by_field[field]
        case["gas_number_density_m3"] = density0
        case["temporal_growth_frequency_s_inv"] = growth0
        case["reduced_temporal_growth_frequency_m3_s"] = growth0 / density0
        density_values.append(density0)
        growth_values.append(growth0)

    return {
        "schema": TWO_TERM_TEMPORAL_GROWTH_TRANSPORT_SCHEMA_VERSION,
        "correction_applied_every_anchor": True,
        "effective_momentum_model": (TWO_TERM_TEMPORAL_GROWTH_EFFECTIVE_MOMENTUM_MODEL),
        "growth_frequency_source": "converged_temporal_growth_eigenvalue",
        "momentum_cross_section_floor_policy": (
            "sum_raw_process_cross_sections_then_single_floor"
        ),
        "minimum_momentum_cross_section_m2": (
            _contracts.TWO_TERM_MOMENTUM_CROSS_SECTION_FLOOR_M2
        ),
        "evidence_columns": list(_contracts.TWO_TERM_TEMPORAL_GROWTH_TRANSPORT_COLUMNS),
        "anchor_points": len(expected_fields),
        "gas_number_density_m3": [
            min(density_values),
            max(density_values),
        ],
        "growth_frequency_s_inv": [
            min(growth_values),
            max(growth_values),
        ],
    }


def _load_rates(
    connection: sqlite3.Connection,
    solver: str,
    mixture_id: int,
    cases: list[dict[str, Any]],
    *,
    allowed_e: set[float] | None = None,
) -> list[dict[str, Any]]:
    by_e = {float(case["E_over_N_Td"]): case for case in cases}
    rows = connection.execute(
        """
        SELECT e_over_n_Td, species, process, process_type, threshold_eV,
               target_species_fraction, energy_loss_eV, mean,
               relative_standard_error, estimate_status
        FROM aggregate_scalars
        WHERE solver = ? AND mixture_id = ? AND scalar_group = 'rate'
          AND scalar_name = ?
        ORDER BY e_over_n_Td, species, process, process_type, threshold_key
        """,
        (solver, mixture_id, RATE_SCALAR),
    ).fetchall()
    rates: list[dict[str, Any]] = []
    for row in rows:
        e_over_n = float(row["e_over_n_Td"])
        if allowed_e is not None and e_over_n not in allowed_e:
            continue
        case = by_e.get(e_over_n)
        if case is None:
            continue
        rate = {
            **_en_fields(e_over_n),
            "mean_energy_eV": case.get("mean_energy_eV"),
            "species": str(row["species"]),
            "process": str(row["process"]),
            "process_type": str(row["process_type"]),
            "threshold_eV": _float_or_none(row["threshold_eV"]),
            "target_species_fraction": _float_or_none(row["target_species_fraction"])
            or 1.0,
            "energy_loss_eV": _float_or_none(row["energy_loss_eV"]),
            "rate_coefficient_m3_s": _float_or_none(row["mean"]),
            "relative_standard_error": _float_or_none(row["relative_standard_error"]),
            "estimate_status": str(row["estimate_status"]),
            "case_drift_velocity_m_s": case.get("drift_velocity_m_s"),
        }
        _fill_rate_derived_columns(rate, drift=case.get("drift_velocity_m_s"))
        rates.append(rate)
    return rates


def _load_eedf(
    connection: sqlite3.Connection,
    solver: str,
    mixture_id: int,
    cases: list[dict[str, Any]],
    *,
    allowed_e: set[float] | None = None,
) -> list[dict[str, Any]]:
    mean_by_e = {
        float(case["E_over_N_Td"]): _float_or_none(case.get("mean_energy_eV"))
        for case in cases
    }
    rows = connection.execute(
        """
        SELECT aggregate.e_over_n_Td, aggregate.energy_eV,
               aggregate.energy_width_eV, aggregate.eedf,
               aggregate.pooled_effective_sample_count
        FROM aggregate_eedf_bins AS aggregate
        WHERE aggregate.solver = ? AND aggregate.mixture_id = ?
        ORDER BY aggregate.e_over_n_Td, aggregate.bin_index
        """,
        (solver, mixture_id),
    ).fetchall()
    eedf: list[dict[str, Any]] = []
    for row in rows:
        e_over_n = float(row["e_over_n_Td"])
        if allowed_e is not None and e_over_n not in allowed_e:
            continue
        if e_over_n not in mean_by_e:
            continue
        eedf.append(
            {
                **_en_fields(e_over_n),
                "mean_energy_eV": mean_by_e[e_over_n],
                "electron_energy_eV": float(row["energy_eV"]),
                "energy_width_eV": float(row["energy_width_eV"]),
                "eedf": float(row["eedf"]),
                "pooled_effective_sample_count": (
                    _float_or_none(row["pooled_effective_sample_count"])
                    if solver == _contracts.MC_SOLVER
                    else None
                ),
            }
        )
    return eedf


def _load_solver_quality_diagnostics(
    connection: sqlite3.Connection,
    solver: str,
    mixture_id: int,
) -> dict[float, dict[str, Any]]:
    """Summarize raw per-case solver diagnostics for table qualification."""

    mc_transport_schema = read_metadata(connection).get(
        "mc_transport_estimator_schema_version"
    )
    mc_eedf_schema = read_metadata(connection).get(
        "mc_eedf_estimator_schema_version"
    )
    if solver == _contracts.MC_SOLVER and mc_transport_schema not in {
        _mc_evidence.DIRECT_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
        _mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
    }:
        raise _contracts.TableBuildError(
            "monte_carlo workflow database lacks the current direct-transport "
            "estimator provenance; regenerate it instead of mixing old cases"
        )
    if (
        solver == _contracts.MC_SOLVER
        and mc_eedf_schema != _mc_evidence.MC_EEDF_ESTIMATOR_SCHEMA_VERSION
    ):
        raise _contracts.TableBuildError(
            "monte_carlo workflow database lacks the current finite-volume "
            "EEDF estimator provenance; regenerate it instead of mixing old cases"
        )

    rows = connection.execute(
        """
        SELECT e_over_n_Td, mean_energy_eV, gas_number_density_m3,
               drift_velocity_m_s, mobility_m2_V_s,
               diffusion_L_m2_s, diffusion_T_m2_s,
               reduced_electron_energy_mobility_m2_V_s_m3,
               reduced_electron_energy_diffusion_L_m2_s_m3,
               reduced_electron_energy_diffusion_T_m2_s_m3,
               diagnostics_json
        FROM cases
        WHERE solver = ? AND mixture_id = ?
        ORDER BY e_over_n_Td, replicate
        """,
        (solver, mixture_id),
    ).fetchall()
    grouped: dict[float, list[dict[str, Any]]] = {}
    for row in rows:
        try:
            payload = json.loads(str(row["diagnostics_json"]))
        except (TypeError, json.JSONDecodeError):
            payload = None
        diagnostic_raw = (
            payload.get("internal_monte_carlo_transport")
            if isinstance(payload, dict) and solver == _contracts.MC_SOLVER
            else payload.get(solver)
            if isinstance(payload, dict)
            else None
        )
        diagnostic = (
            dict(diagnostic_raw)
            if isinstance(diagnostic_raw, dict)
            else {"_diagnostic_invalid": True}
        )
        if solver == _contracts.MC_SOLVER and not diagnostic.get("_diagnostic_invalid"):
            density = _float_or_none(row["gas_number_density_m3"])
            try:
                if density is None or density <= 0.0:
                    raise ValueError("invalid gas density")
                diagnostic["_case_mean_energy_eV"] = float(row["mean_energy_eV"])
                diagnostic["_case_reported_transport"] = {
                    "drift_velocity_m_s": float(row["drift_velocity_m_s"]),
                    "mobility_m2_V_s": float(row["mobility_m2_V_s"]),
                    "diffusion_L_m2_s": float(row["diffusion_L_m2_s"]),
                    "diffusion_T_m2_s": float(row["diffusion_T_m2_s"]),
                    "energy_mobility_m2_V_s": float(
                        row["reduced_electron_energy_mobility_m2_V_s_m3"]
                    )
                    / density,
                    "energy_diffusion_L_m2_s": float(
                        row["reduced_electron_energy_diffusion_L_m2_s_m3"]
                    )
                    / density,
                    "energy_diffusion_T_m2_s": float(
                        row["reduced_electron_energy_diffusion_T_m2_s_m3"]
                    )
                    / density,
                }
            except (TypeError, ValueError, OverflowError):
                diagnostic = {"_diagnostic_invalid": True}
        grouped.setdefault(float(row["e_over_n_Td"]), []).append(diagnostic)

    summaries: dict[float, dict[str, Any]] = {}
    for e_over_n, diagnostics in grouped.items():
        if solver == _contracts.MC_SOLVER:
            summaries[e_over_n] = (
                summarize_weighted_mc_transport_diagnostics(diagnostics)
                if mc_transport_schema
                == _mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
                else summarize_mc_transport_diagnostics(diagnostics)
            )
            continue
        converged = all(bool(item.get("converged", False)) for item in diagnostics)
        iterations = max(int(item.get("iterations", 0)) for item in diagnostics)
        residual = max(float(item.get("residual_L1", math.inf)) for item in diagnostics)
        residual_target = min(
            float(item.get("residual_tolerance", math.nan)) for item in diagnostics
        )
        tail_probability = max(
            float(item.get("tail_probability", math.inf)) for item in diagnostics
        )
        tail_target = min(
            float(item.get("tail_probability_target", math.nan)) for item in diagnostics
        )
        edge_to_peak = max(
            float(item.get("edge_to_peak", math.inf)) for item in diagnostics
        )
        edge_target = min(
            float(item.get("edge_to_peak_target", math.nan)) for item in diagnostics
        )
        grid_max = max(float(item.get("grid_max_eV", math.inf)) for item in diagnostics)
        grid_limit = min(
            float(item.get("grid_max_limit_eV", math.nan)) for item in diagnostics
        )
        thresholds_available = all(
            math.isfinite(value)
            for value in (
                residual_target,
                tail_target,
                edge_target,
                grid_limit,
            )
        )
        configured_ceiling_domain = all(
            item.get("energy_domain_strategy")
            == "one_shot_core_plus_stretched_tail_to_configured_ceiling"
            for item in diagnostics
        )
        grid_limit_hit = bool(
            thresholds_available
            and grid_max >= grid_limit * (1.0 - 1.0e-9)
            and not configured_ceiling_domain
        )
        passed = bool(
            thresholds_available
            and converged
            and residual <= residual_target
            and tail_probability <= tail_target
            and edge_to_peak <= edge_target
            and not grid_limit_hit
        )
        reasons: list[str] = []
        if not thresholds_available:
            reasons.append("solver_diagnostic_thresholds_missing")
        if not converged:
            reasons.append("solver_not_converged")
        if math.isfinite(residual_target) and residual > residual_target:
            reasons.append("solver_residual_above_tolerance")
        if math.isfinite(tail_target) and tail_probability > tail_target:
            reasons.append("eedf_tail_probability_above_target")
        if math.isfinite(edge_target) and edge_to_peak > edge_target:
            reasons.append("eedf_edge_to_peak_above_target")
        if grid_limit_hit:
            reasons.append("energy_grid_limit_hit")
        summaries[e_over_n] = {
            "solver_diagnostics_available": 1,
            "solver_converged": int(converged),
            "solver_iterations_max": iterations,
            "solver_residual_L1_max": residual,
            "solver_residual_tolerance": _float_or_none(residual_target),
            "solver_tail_probability_max": tail_probability,
            "solver_tail_probability_target": _float_or_none(tail_target),
            "solver_edge_to_peak_max": edge_to_peak,
            "solver_edge_to_peak_target": _float_or_none(edge_target),
            "solver_grid_max_eV_max": grid_max,
            "solver_grid_max_limit_eV": _float_or_none(grid_limit),
            "solver_grid_limit_hit": int(grid_limit_hit),
            "solver_diagnostics_passed": int(passed),
            "solver_failure_reasons": reasons,
        }
    return summaries


def _load_quality(
    connection: sqlite3.Connection,
    solver: str,
    mixture_id: int,
    *,
    allowed_e: set[float] | None = None,
) -> list[dict[str, Any]]:
    diagnostics_by_e = _load_solver_quality_diagnostics(connection, solver, mixture_id)
    rows = connection.execute(
        """
        SELECT e_over_n_Td, passed, failure_reasons_json, thresholds_json,
               quality_policy_reevaluated,
               mobility_rse, diffusion_L_rse, diffusion_T_rse,
               energy_mobility_rse, energy_diffusion_L_rse,
               energy_diffusion_T_rse,
               max_major_rate_rse, eedf_normalization_error,
               valid_replicates, uncertainty_available
        FROM aggregate_quality
        WHERE solver = ? AND mixture_id = ?
        ORDER BY e_over_n_Td
        """,
        (solver, mixture_id),
    ).fetchall()
    quality: list[dict[str, Any]] = []
    for row in rows:
        e_over_n = float(row["e_over_n_Td"])
        if allowed_e is not None and e_over_n not in allowed_e:
            continue
        diagnostic = diagnostics_by_e.get(e_over_n)
        try:
            failure_reasons = json.loads(str(row["failure_reasons_json"]))
        except json.JSONDecodeError:
            failure_reasons = [str(row["failure_reasons_json"])]
        if not isinstance(failure_reasons, list):
            failure_reasons = [str(failure_reasons)]
        aggregate_failure_reasons = list(failure_reasons)
        try:
            evaluation_thresholds = parse_quality_thresholds(
                json.loads(str(row["thresholds_json"]))
            )
        except (TypeError, ValueError, json.JSONDecodeError) as exc:
            raise _contracts.TableBuildError(
                "aggregate quality row has invalid evaluation policy"
            ) from exc
        aggregate_passed = bool(row["passed"])
        if diagnostic is not None:
            failure_reasons.extend(diagnostic["solver_failure_reasons"])
            passed = aggregate_passed and bool(diagnostic["solver_diagnostics_passed"])
        else:
            passed = aggregate_passed and solver != _contracts.MC_SOLVER
            missing_reasons = (
                ["mc_transport_replica_diagnostic_invalid"]
                if solver == _contracts.MC_SOLVER
                else []
            )
            failure_reasons.extend(missing_reasons)
            diagnostic = {
                "solver_diagnostics_available": 0,
                "solver_converged": None,
                "solver_transport_qualified": (
                    0 if solver == _contracts.MC_SOLVER else None
                ),
                "solver_iterations_max": None,
                "solver_residual_L1_max": None,
                "solver_residual_tolerance": None,
                "solver_tail_probability_max": None,
                "solver_tail_probability_target": None,
                "solver_edge_to_peak_max": None,
                "solver_edge_to_peak_target": None,
                "solver_grid_max_eV_max": None,
                "solver_grid_max_limit_eV": None,
                "solver_grid_limit_hit": None,
                "solver_diagnostics_passed": None,
                "solver_transport_mean_energy_max_relative_ci95_bound": None,
                "solver_mean_energy_stationarity_relative_ci95_bound": None,
                "solver_mobility_stationarity_absolute_log_drift": None,
                "solver_population_growth_gate_mode": None,
                "solver_population_growth_max_relative_ci95_bound": None,
                "solver_population_growth_poisson_interval_max_ratio": None,
                "solver_population_growth_sparse_max_metric": None,
                "solver_population_growth_pooled_event_count": None,
                "solver_population_growth_pooled_exposure_s": None,
                "solver_population_growth_pooled_frequency_ci95_low_s_inv": None,
                "solver_population_growth_pooled_frequency_ci95_high_s_inv": None,
                "solver_population_growth_pooled_direct_within_ci95": None,
                "solver_population_growth_pooled_population_within_ci95": None,
                "solver_population_growth_relative_tolerance": None,
                "solver_failure_reasons": [],
            }
        quality.append(
            {
                **_en_fields(e_over_n),
                "passed": int(passed),
                "aggregate_quality_passed": int(aggregate_passed),
                "failure_reasons_json": json.dumps(
                    sorted(set(str(value) for value in failure_reasons)),
                    separators=(",", ":"),
                ),
                "aggregate_failure_reasons_json": json.dumps(
                    sorted(set(str(value) for value in aggregate_failure_reasons)),
                    separators=(",", ":"),
                ),
                "thresholds_json": str(row["thresholds_json"]),
                "mobility_rse": _float_or_none(row["mobility_rse"]),
                "diffusion_L_rse": _float_or_none(row["diffusion_L_rse"]),
                "diffusion_T_rse": _float_or_none(row["diffusion_T_rse"]),
                "energy_mobility_rse": _float_or_none(row["energy_mobility_rse"]),
                "energy_diffusion_L_rse": _float_or_none(row["energy_diffusion_L_rse"]),
                "energy_diffusion_T_rse": _float_or_none(row["energy_diffusion_T_rse"]),
                "max_major_rate_rse": _float_or_none(row["max_major_rate_rse"]),
                "eedf_normalization_error": _float_or_none(
                    row["eedf_normalization_error"]
                ),
                "valid_replicates": int(row["valid_replicates"]),
                "uncertainty_available": int(row["uncertainty_available"]),
                "required_rate_min_process_peak_fraction": (
                    evaluation_thresholds.required_rate_min_process_peak_fraction
                ),
                "quality_policy_reevaluated": int(row["quality_policy_reevaluated"]),
                **{
                    key: value
                    for key, value in diagnostic.items()
                    if key != "solver_failure_reasons"
                },
                "quality_source": f"{solver}_aggregate_quality",
            }
        )
    return quality


def _fill_rate_derived_columns(
    row: dict[str, Any], drift: object | None = None
) -> None:
    drift_value = _float_or_none(
        drift if drift is not None else row.get("drift_velocity_m_s")
    )
    if drift_value is None:
        drift_value = _float_or_none(row.get("case_drift_velocity_m_s"))
    if drift_value is None or abs(drift_value) <= 0.0:
        raise _contracts.TableBuildError(
            "cannot compute reduced Townsend with zero drift"
        )
    k = _required_float(row, "rate_coefficient_m3_s")
    fraction = _required_float(row, "target_species_fraction")
    energy_loss = _float_or_none(row.get("energy_loss_eV")) or 0.0
    row["mixture_weighted_rate_m3_s"] = k * fraction
    row["reduced_townsend_m2"] = k / abs(drift_value)
    row["mixture_weighted_reduced_townsend_m2"] = k * fraction / abs(drift_value)
    row["energy_loss_rate_coefficient_eV_m3_s"] = k * energy_loss


def _en_fields(e_over_n_Td: float) -> dict[str, float]:
    return {
        "E_over_N_Td": float(e_over_n_Td),
        "E_over_N_V_m2": float(e_over_n_Td) * _contracts.TD_TO_V_M2,
    }


def _validate_mc_sampling_plan_against_cases(
    connection: sqlite3.Connection,
    metadata: dict[str, str],
) -> dict[str, object]:
    """Bind canonical workflow sampling controls to every stored MC replica."""

    sampling = mc_sampling_plan_provenance(metadata, required=True)
    assert sampling is not None
    entries = sampling["entries"]
    assert isinstance(entries, list)
    by_anchor = {float(entry["e_over_n_Td"]): entry for entry in entries}
    rows = connection.execute(
        """
        SELECT mixture_id, e_over_n_Td, replicate, diagnostics_json
        FROM cases
        WHERE solver = ?
        ORDER BY mixture_id, e_over_n_Td, replicate
        """,
        (_contracts.MC_SOLVER,),
    ).fetchall()
    base_seed = int(metadata["mc_base_seed"]) if "mc_base_seed" in metadata else None
    solver_source = metadata.get("mc_solver_source_sha256")
    if solver_source != monte_carlo_source_sha256():
        raise WorkflowSchemaError(
            "Monte Carlo solver source differs from the workflow database"
        )
    if base_seed is not None:
        for row in rows:
            provenance = (
                json.loads(row["diagnostics_json"])
                .get("internal_monte_carlo_transport", {})
                .get("mc_run_provenance", {})
            )
            expected_seed = stable_mc_seed(
                base_seed=base_seed,
                mixture_id=int(row["mixture_id"]),
                e_over_n_Td=float(row["e_over_n_Td"]),
                replicate=int(row["replicate"]),
            )
            if provenance.get("seed") != expected_seed:
                raise WorkflowSchemaError(
                    "MC base seed does not match per-case run provenance"
                )
            if provenance.get("solver_source_sha256") != solver_source:
                raise WorkflowSchemaError(
                    "MC solver source does not match per-case run provenance"
                )
    if not rows:
        raise WorkflowSchemaError(
            "MC sampling-plan provenance has no stored Monte Carlo cases"
        )

    grouped: dict[tuple[int, float], list[sqlite3.Row]] = {}
    mixture_ids: set[int] = set()
    for row in rows:
        mixture_id = int(row["mixture_id"])
        anchor = float(row["e_over_n_Td"])
        mixture_ids.add(mixture_id)
        grouped.setdefault((mixture_id, anchor), []).append(row)

    expected_anchors = set(by_anchor)
    weighted_growth = (
        metadata.get("mc_transport_estimator_schema_version")
        == _mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
    )
    for mixture_id in sorted(mixture_ids):
        actual_anchors = {
            anchor
            for candidate_mixture, anchor in grouped
            if candidate_mixture == mixture_id
        }
        if actual_anchors != expected_anchors:
            raise WorkflowSchemaError(
                "MC sampling-plan anchors do not match stored cases for "
                f"mixture {mixture_id}"
            )
        for anchor in sorted(expected_anchors):
            entry = by_anchor[anchor]
            anchor_rows = grouped[(mixture_id, anchor)]
            expected_replicates = int(entry["replicas"])
            actual_replicates = [int(row["replicate"]) for row in anchor_rows]
            if actual_replicates != list(range(expected_replicates)):
                raise WorkflowSchemaError(
                    "MC sampling-plan replica count/index does not match stored "
                    f"cases at {anchor:g} Td for mixture {mixture_id}"
                )
            expected_controls = (
                int(entry["particles"]),
                int(entry["warmup_collisions"]),
                int(entry["max_collisions"]),
                int(entry["tail_max_collisions"]),
            )
            expected_lag = int(entry["transport_correlation_lag_barriers"])
            expected_transport_estimator = str(entry["transport_estimator"])
            for row in anchor_rows:
                try:
                    diagnostics = json.loads(str(row["diagnostics_json"]))
                    transport = diagnostics["internal_monte_carlo_transport"]
                    run_provenance = transport["mc_run_provenance"]
                    actual_controls = (
                        int(run_provenance["particles"]),
                        int(run_provenance["warmup_collisions"]),
                        int(run_provenance["production_collisions"]),
                        int(
                            run_provenance.get(
                                "tail_max_collisions",
                                run_provenance["production_collisions"],
                            )
                        ),
                    )
                    actual_lag = (
                        int(
                            transport["block_lag_sampling"][
                                "configured_correlation_lag_barriers"
                            ]
                        )
                        if weighted_growth
                        else expected_lag
                    )
                    actual_transport_estimator = str(
                        run_provenance["transport_estimator"]
                    )
                    if (
                        expected_transport_estimator == "paired_field_parity"
                        and actual_transport_estimator == "paired_field_parity"
                    ):
                        parity_schema = transport["field_parity_response"][
                            "estimator_schema_version"
                        ]
                        if (
                            parity_schema
                            != _mc_evidence.FIELD_PARITY_RESPONSE_ESTIMATOR_SCHEMA_VERSION
                        ):
                            raise ValueError(
                                "field-parity response schema is incompatible"
                            )
                except (
                    KeyError,
                    TypeError,
                    ValueError,
                    json.JSONDecodeError,
                ) as exc:
                    raise WorkflowSchemaError(
                        "MC sampling-plan validation requires per-case run provenance "
                        f"at {anchor:g} Td for mixture {mixture_id}"
                    ) from exc
                if (
                    actual_controls != expected_controls
                    or actual_lag != expected_lag
                    or actual_transport_estimator != expected_transport_estimator
                ):
                    raise WorkflowSchemaError(
                        "MC sampling-plan controls do not match per-case run provenance "
                        f"at {anchor:g} Td for mixture {mixture_id}"
                    )

    return {
        **sampling,
        "status": "validated_against_mc_cases",
        "mixtures": len(mixture_ids),
        "base_seed": base_seed,
    }


def _aggregate_rows_are_current(connection: sqlite3.Connection) -> bool:
    """Return whether aggregates cover every currently committed raw case.

    A resumable MC sweep may append cases after someone produced a diagnostic
    partial aggregate.  Presence of aggregate tables alone is therefore not a
    freshness signal.  The workflow store treats committed case keys as
    immutable, so exact group/replica-count agreement is sufficient here.
    """

    required = (
        "aggregate_quality",
        "aggregate_scalars",
        "aggregate_eedf_bins",
    )
    if any(not table_exists(connection, name) for name in required):
        return False
    aggregate_eedf_columns = {
        str(row[1])
        for row in connection.execute(
            "PRAGMA table_info(aggregate_eedf_bins)"
        ).fetchall()
    }
    if "pooled_effective_sample_count" not in aggregate_eedf_columns:
        return False
    raw = {
        (int(row[0]), str(row[1]), float(row[2])): int(row[3])
        for row in connection.execute(
            """
            SELECT mixture_id, solver, e_over_n_Td, COUNT(*)
            FROM cases
            GROUP BY mixture_id, solver, e_over_n_Td
            """
        )
    }
    aggregated = {
        (int(row[0]), str(row[1]), float(row[2])): int(row[3])
        for row in connection.execute(
            """
            SELECT mixture_id, solver, e_over_n_Td, n
            FROM aggregate_quality
            """
        )
    }
    if not raw or raw != aggregated:
        return False
    scalar_count = connection.execute(
        "SELECT COUNT(*) FROM aggregate_scalars"
    ).fetchone()
    eedf_count = connection.execute(
        "SELECT COUNT(*) FROM aggregate_eedf_bins"
    ).fetchone()
    return bool(
        scalar_count
        and int(scalar_count[0]) > 0
        and eedf_count
        and int(eedf_count[0]) > 0
    )


def _rate_key(row: dict[str, Any]) -> tuple[str, str, str, str]:
    threshold = row.get("threshold_eV")
    threshold_key = "" if threshold is None else f"{float(threshold):.17g}"
    return (
        str(row["species"]),
        str(row["process"]),
        str(row["process_type"]),
        threshold_key,
    )


def _finite_nonnegative(value: object) -> float | None:
    number = _float_or_none(value)
    if number is None or number < 0.0:
        return None
    return number


def _required_float(row: dict[str, Any], name: str) -> float:
    value = _float_or_none(row.get(name))
    if value is None:
        raise _contracts.TableBuildError(f"missing required value {name}")
    return value


def _float_or_none(value: object) -> float | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None
