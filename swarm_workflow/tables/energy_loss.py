"""Solver-owned elastic energy-loss table materialization."""

from __future__ import annotations

import json
import sqlite3
from typing import Any

from . import contracts as _contracts
from . import repository as _repository
from ..campaign.statistics import summarize_replicates
from .math import strictly_monotonic

def _load_elastic_energy_loss(
    connection: sqlite3.Connection,
    solver: str,
    mixture_id: int,
    cases: list[dict[str, Any]],
    *,
    allowed_e: set[float],
) -> tuple[list[dict[str, Any]] | None, dict[str, Any] | None]:
    """Load an all-or-none solver-native elastic energy-loss closure."""

    mean_by_e = {
        float(case["E_over_N_Td"]): _repository._required_float(case, "mean_energy_eV")
        for case in cases
    }
    rows = connection.execute(
        """
        SELECT e_over_n_Td, replicate, diagnostics_json
        FROM cases
        WHERE solver = ? AND mixture_id = ?
        ORDER BY e_over_n_Td, replicate
        """,
        (solver, mixture_id),
    ).fetchall()
    relevant = [row for row in rows if float(row["e_over_n_Td"]) in allowed_e]
    if not relevant:
        return None, None

    extracted: list[
        tuple[float, int, float, dict[str, Any], int | None] | None
    ] = []
    for row in relevant:
        try:
            payload = json.loads(str(row["diagnostics_json"]))
        except (TypeError, json.JSONDecodeError) as exc:
            raise _contracts.TableBuildError(
                "elastic energy-loss closure has invalid solver diagnostics"
            ) from exc
        if not isinstance(payload, dict):
            extracted.append(None)
            continue
        if solver == _contracts.TWO_TERM_SOLVER:
            record = _two_term_elastic_energy_loss_record(payload)
        elif solver == _contracts.PROPAGATOR_SOLVER:
            record = _propagator_elastic_energy_loss_record(payload)
        else:
            record = _mc_elastic_energy_loss_record(payload)
        if record is None:
            extracted.append(None)
            continue
        coefficient, contract = record
        seed = _mc_replica_seed(payload) if solver == _contracts.MC_SOLVER else None
        extracted.append(
            (
                float(row["e_over_n_Td"]),
                int(row["replicate"]),
                coefficient,
                contract,
                seed,
            )
        )

    if all(item is None for item in extracted):
        return None, None
    if any(item is None for item in extracted):
        raise _contracts.TableBuildError(
            "elastic energy-loss closure is present for only part of the "
            f"selected {solver} replicas; regenerate one consistent database"
        )

    records = [item for item in extracted if item is not None]
    contracts = [item[3] for item in records]
    canonical_contract = contracts[0]
    if any(contract != canonical_contract for contract in contracts[1:]):
        raise _contracts.TableBuildError(
            "elastic energy-loss closure model differs across selected replicas"
        )
    values_by_e: dict[float, list[float]] = {}
    seeds_by_e: dict[float, list[int]] = {}
    for e_over_n, _replicate, value, _contract, seed in records:
        values_by_e.setdefault(e_over_n, []).append(value)
        if seed is not None:
            seeds_by_e.setdefault(e_over_n, []).append(seed)
    if solver == _contracts.MC_SOLVER:
        for e_over_n, values in values_by_e.items():
            seeds = seeds_by_e.get(e_over_n, [])
            if len(seeds) != len(values) or len(set(seeds)) != len(seeds):
                raise _contracts.TableBuildError(
                    "monte_carlo elastic energy-loss uncertainty requires "
                    f"distinct replica seeds at {e_over_n:.17g} Td"
                )
    if set(values_by_e) != allowed_e or set(mean_by_e) != allowed_e:
        raise _contracts.TableBuildError(
            "elastic energy-loss closure does not cover every selected E/N anchor"
        )

    output: list[dict[str, Any]] = []
    for e_over_n in sorted(allowed_e):
        stats = summarize_replicates(values_by_e[e_over_n])
        if stats.mean is None:
            raise _contracts.TableBuildError("elastic energy-loss aggregate is unavailable")
        if solver == _contracts.MC_SOLVER and not stats.uncertainty_available:
            raise _contracts.TableBuildError(
                "monte_carlo elastic energy loss requires independent-replica "
                f"uncertainty at {e_over_n:.17g} Td"
            )
        output.append(
            {
                "mean_energy_eV": mean_by_e[e_over_n],
                **_repository._en_fields(e_over_n),
                "elastic_energy_loss_rate_coefficient_eV_m3_s": stats.mean,
                "elastic_energy_loss_standard_error_eV_m3_s": (
                    stats.standard_error
                ),
                "elastic_energy_loss_relative_standard_error": (
                    stats.relative_standard_error
                ),
                "elastic_energy_loss_ci95_low_eV_m3_s": stats.ci95_low,
                "elastic_energy_loss_ci95_high_eV_m3_s": stats.ci95_high,
                "ci95_critical_value": stats.ci95_critical_value,
                "estimate_status": (
                    "deterministic_operator_moment"
                    if solver in {_contracts.TWO_TERM_SOLVER, _contracts.PROPAGATOR_SOLVER}
                    else stats.estimate_status
                ),
                "valid_replicates": stats.valid_replicates,
                "uncertainty_available": int(stats.uncertainty_available),
            }
        )
    if not strictly_monotonic(row["mean_energy_eV"] for row in output):
        raise _contracts.TableBuildError(
            "elastic energy-loss closure cannot be represented as a "
            "single-valued function of mean energy"
        )
    return output, canonical_contract


def _mc_replica_seed(diagnostics: dict[str, Any]) -> int:
    transport = diagnostics.get("internal_monte_carlo_transport")
    provenance = (
        transport.get("mc_run_provenance")
        if isinstance(transport, dict)
        else None
    )
    seed = provenance.get("seed") if isinstance(provenance, dict) else None
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise _contracts.TableBuildError(
            "monte_carlo elastic energy-loss artifact lacks an integer "
            "replica seed"
        )
    return seed


def _propagator_elastic_energy_loss_record(
    diagnostics: dict[str, Any],
) -> tuple[float, dict[str, Any]] | None:
    section = diagnostics.get(_contracts.PROPAGATOR_SOLVER)
    if not isinstance(section, dict):
        return None
    artifact = section.get("elastic_energy_loss")
    if artifact is None:
        return None
    if not isinstance(artifact, dict):
        raise _contracts.TableBuildError(
            "propagator elastic energy-loss artifact is invalid"
        )
    required = {
        "schema": "swarm.elastic_energy_loss.v1",
        "status": "available",
        "symbol": "K_epsilon_el",
        "estimator": "same_discrete_finite_temperature_elastic_sg_operator",
        "operator": (
            "finite_volume_scharfetter_gummel_elastic_energy_generator"
        ),
        "source_eedf": "same_solved_energy_angle_distribution",
        "neutral_thermal_motion_model": "finite_temperature_fokker_planck",
        "gas_temperature_terms_included": True,
        "sign_convention": "positive_is_net_electron_energy_loss",
        "uncertainty_status": "deterministic_kinetic_solver",
    }
    if any(artifact.get(key) != value for key, value in required.items()):
        raise _contracts.TableBuildError(
            "propagator elastic energy-loss artifact does not match its "
            "finite-temperature discrete collision operator"
        )
    value = _repository._float_or_none(artifact.get("rate_coefficient_eV_m3_s"))
    gas_temperature = _repository._float_or_none(artifact.get("gas_temperature_K"))
    if value is None or gas_temperature is None or gas_temperature <= 0.0:
        raise _contracts.TableBuildError(
            "propagator elastic energy-loss artifact lacks a finite coefficient "
            "or gas temperature"
        )
    return value, {
        "schema": required["schema"],
        "symbol": required["symbol"],
        "estimator": required["estimator"],
        "operator": required["operator"],
        "source_eedf": required["source_eedf"],
        "neutral_thermal_motion_model": required[
            "neutral_thermal_motion_model"
        ],
        "gas_temperature_terms_included": True,
        "gas_temperature_K": gas_temperature,
        "sign_convention": required["sign_convention"],
        "uncertainty": "deterministic_kinetic_solver",
    }


def _two_term_elastic_energy_loss_record(
    diagnostics: dict[str, Any],
) -> tuple[float, dict[str, Any]] | None:
    section = diagnostics.get(_contracts.TWO_TERM_SOLVER)
    if not isinstance(section, dict):
        return None
    artifact = section.get("elastic_energy_loss")
    if artifact is None:
        return None
    if not isinstance(artifact, dict):
        raise _contracts.TableBuildError("two_term elastic energy-loss artifact is invalid")
    required = {
        "schema": "swarm.elastic_energy_loss.v1",
        "status": "available",
        "symbol": "K_epsilon_el",
        "estimator": "same_discrete_elastic_collision_operator_energy_moment",
        "operator": (
            "native_finite_volume_scharfetter_gummel_elastic_A_D_zero_field"
        ),
        "operator_isolation": (
            "zero_field_reassembly_from_same_elastic_A_D_coefficients"
        ),
        "source_eedf": "same_solved_eedf",
        "neutral_thermal_motion_model": "finite_temperature_fokker_planck",
        "gas_temperature_terms_included": True,
        "sign_convention": "positive_is_net_electron_energy_loss",
    }
    if any(artifact.get(key) != value for key, value in required.items()):
        raise _contracts.TableBuildError(
            "two_term elastic energy-loss artifact is not the same discrete "
            "finite-temperature collision-operator moment"
        )
    value = _repository._float_or_none(artifact.get("rate_coefficient_eV_m3_s"))
    gas_temperature = _repository._float_or_none(artifact.get("gas_temperature_K"))
    if value is None or gas_temperature is None or gas_temperature <= 0.0:
        raise _contracts.TableBuildError(
            "two_term elastic energy-loss artifact lacks a finite coefficient "
            "or gas temperature"
        )
    return value, {
        "schema": "swarm.elastic_energy_loss.v1",
        "symbol": "K_epsilon_el",
        "estimator": required["estimator"],
        "operator": required["operator"],
        "operator_isolation": required["operator_isolation"],
        "source_eedf": required["source_eedf"],
        "neutral_thermal_motion_model": required[
            "neutral_thermal_motion_model"
        ],
        "gas_temperature_terms_included": True,
        "gas_temperature_K": gas_temperature,
        "sign_convention": required["sign_convention"],
        "uncertainty": "deterministic_kinetic_solver",
    }


def _mc_elastic_energy_loss_record(
    diagnostics: dict[str, Any],
) -> tuple[float, dict[str, Any]] | None:
    section = diagnostics.get("internal_monte_carlo_reaction_rates")
    if not isinstance(section, dict):
        return None
    rates = section.get("rates")
    if not isinstance(rates, list):
        return None
    transport = diagnostics.get("internal_monte_carlo_transport")
    provenance = (
        transport.get("mc_run_provenance")
        if isinstance(transport, dict)
        else None
    )
    gas_temperature = _repository._float_or_none(
        provenance.get("gas_temperature_K")
        if isinstance(provenance, dict)
        else None
    )
    scattering_types = {"elastic", "effective", "momentum"}
    supported: list[tuple[float, float, str]] = []
    saw_sampled_scattering = False
    missing_sampled_loss = False
    for item in rates:
        if not isinstance(item, dict):
            raise _contracts.TableBuildError("monte_carlo reaction-rate diagnostic is invalid")
        if str(item.get("process_type", "")).lower() not in scattering_types:
            continue
        if item.get("event_sampling_enabled") is not True:
            continue
        saw_sampled_scattering = True
        loss = item.get("energy_loss")
        if not isinstance(loss, dict):
            missing_sampled_loss = True
            continue
        model = str(loss.get("model", ""))
        if (
            loss.get("status") != "direct_event_estimator"
            or loss.get("estimator")
            != "trajectory_event_energy_change_per_target_density_residence_time"
            or model
            not in {
                "maxwellian_target_exact_binary_collision_isotropic",
                "maxwellian_target_exact_binary_collision_maxent_p1",
            }
        ):
            raise _contracts.TableBuildError(
                "monte_carlo elastic energy loss is not the direct "
                "finite-temperature binary-collision event estimator"
            )
        neutral_model = loss.get("neutral_thermal_motion_model")
        gas_temperature_terms = loss.get("gas_temperature_terms_included")
        if (
            neutral_model
            != "maxwellian_relative_speed_exact_binary_collision"
            or gas_temperature_terms is not True
        ):
            raise _contracts.TableBuildError(
                "monte_carlo elastic energy-loss thermal model is inconsistent"
            )
        value = _repository._float_or_none(loss.get("rate_coefficient_eV_m3_s"))
        fraction = _repository._float_or_none(item.get("target_species_fraction"))
        if value is None or fraction is None or fraction <= 0.0:
            raise _contracts.TableBuildError(
                "monte_carlo elastic energy-loss coefficient or target fraction "
                "is invalid"
            )
        supported.append((fraction, value, model))
    if not saw_sampled_scattering:
        return None
    if gas_temperature is None or gas_temperature <= 0.0:
        raise _contracts.TableBuildError(
            "monte_carlo exact thermal elastic energy loss lacks gas temperature"
        )
    if not supported:
        return None
    if missing_sampled_loss:
        raise _contracts.TableBuildError(
            "monte_carlo direct elastic energy loss is present for only part "
            "of the sampled scattering channels"
        )
    return sum(fraction * value for fraction, value, _model in supported), {
        "schema": "swarm.elastic_energy_loss.v1",
        "symbol": "K_epsilon_el",
        "estimator": (
            "trajectory_event_energy_change_per_target_density_residence_time"
        ),
        "aggregation": "sum_target_fraction_times_process_coefficient",
        "source_eedf": "same_sampled_trajectory_event_and_residence_measure",
        "event_models": sorted({model for _fraction, _value, model in supported}),
        "neutral_thermal_motion_model": (
            "maxwellian_relative_speed_exact_binary_collision"
        ),
        "gas_temperature_terms_included": True,
        "gas_temperature_K": gas_temperature,
        "sign_convention": "positive_is_net_electron_energy_loss",
        "uncertainty": "independent_replica_student_t_95",
    }
