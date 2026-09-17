"""Validate the external elastic-energy-loss table and uncertainty contract."""

from __future__ import annotations

import math
from typing import Any

from ..contracts import (
    GEC_EXTERNAL_ELASTIC_ENERGY_LOSS,
    GecCcpMapping,
    GecCcpWorkflowError,
)
from ..data import _float_or_none, _read_csv, _required_float
from . import bundle_guards as _bundle_guards
from .bundle_context import BundleUsage, ManifestEvidence
from .table_inputs import positive_numeric_values as _positive_numeric_values


def validate_elastic_energy_loss(
    mapping: GecCcpMapping,
    usage: BundleUsage,
    manifest: ManifestEvidence,
) -> dict[str, Any] | None:
    if not usage.external_elastic_loss:
        return None
    entry = manifest.tables[_bundle_guards.GEC_ELASTIC_ENERGY_LOSS_TABLE]
    physics_contract = entry.get("physics_contract")
    if not isinstance(physics_contract, dict):
        raise GecCcpWorkflowError(
            "external elastic energy-loss table lacks physics_contract"
        )
    _validate_elastic_physics_contract(
        source=manifest.source,
        physics_contract=physics_contract,
    )
    elastic_rows = _read_csv(
        mapping.bundle.path / _bundle_guards.GEC_ELASTIC_ENERGY_LOSS_TABLE
    )
    elastic_energy = _positive_numeric_values(
        elastic_rows,
        "mean_energy_eV",
    )
    coefficients = _positive_numeric_values(
        elastic_rows,
        _bundle_guards.GEC_ELASTIC_ENERGY_LOSS_COLUMN,
    )
    if len(elastic_rows) < 2 or any(
        right <= left for left, right in zip(elastic_energy, elastic_energy[1:])
    ):
        raise GecCcpWorkflowError(
            "external elastic energy-loss table requires at least two "
            "strictly increasing mean-energy anchors"
        )
    statistical_evidence = _elastic_statistical_evidence(
        source=manifest.source,
        rows=elastic_rows,
        coefficients=coefficients,
        quality_thresholds=manifest.quality_thresholds,
    )
    return {
        "model": GEC_EXTERNAL_ELASTIC_ENERGY_LOSS,
        "table": _bundle_guards.GEC_ELASTIC_ENERGY_LOSS_TABLE,
        "table_sha256": manifest.artifact_digests[
            _bundle_guards.GEC_ELASTIC_ENERGY_LOSS_TABLE
        ],
        "source": manifest.source,
        "physics_contract": physics_contract,
        "mean_energy_range_eV": [
            float(elastic_energy[0]),
            float(elastic_energy[-1]),
        ],
        "coefficient_range_eV_m3_s": [
            float(min(coefficients)),
            float(max(coefficients)),
        ],
        "points": len(elastic_rows),
        "statistical_evidence": statistical_evidence,
        "energy_equation_owner": "GeneralPowerDeposition:swElLoss",
        "COMSOL_elastic_reaction_owner": "disabled:eir1",
    }


def _validate_elastic_physics_contract(
    *,
    source: str,
    physics_contract: dict[str, Any],
) -> None:
    common_contract = {
        "schema": "swarm.elastic_energy_loss.v1",
        "symbol": "K_epsilon_el",
        "sign_convention": "positive_is_net_electron_energy_loss",
    }
    if source == "two_term":
        required_contract: dict[str, Any] = {
            **common_contract,
            "estimator": ("same_discrete_elastic_collision_operator_energy_moment"),
            "operator": (
                "native_finite_volume_scharfetter_gummel_elastic_A_D_zero_field"
            ),
            "operator_isolation": (
                "zero_field_reassembly_from_same_elastic_A_D_coefficients"
            ),
            "source_eedf": "same_solved_eedf",
            "neutral_thermal_motion_model": ("finite_temperature_fokker_planck"),
            "gas_temperature_terms_included": True,
            "uncertainty": "deterministic_kinetic_solver",
        }
        candidate = dict(physics_contract)
        gas_temperature = _float_or_none(candidate.pop("gas_temperature_K", None))
        contract_valid = bool(
            candidate == required_contract
            and gas_temperature is not None
            and gas_temperature > 0.0
        )
    elif source == "monte_carlo":
        required_contract = {
            **common_contract,
            "estimator": (
                "trajectory_event_energy_change_per_target_density_residence_time"
            ),
            "aggregation": "sum_target_fraction_times_process_coefficient",
            "source_eedf": ("same_sampled_trajectory_event_and_residence_measure"),
            "neutral_thermal_motion_model": (
                "maxwellian_relative_speed_exact_binary_collision"
            ),
            "gas_temperature_terms_included": True,
            "uncertainty": "independent_replica_student_t_95",
        }
        event_models = physics_contract.get("event_models")
        allowed_event_models = {
            "maxwellian_target_exact_binary_collision_isotropic",
            "maxwellian_target_exact_binary_collision_maxent_p1",
        }
        candidate = dict(physics_contract)
        candidate.pop("event_models", None)
        gas_temperature = _float_or_none(candidate.pop("gas_temperature_K", None))
        contract_valid = bool(
            candidate == required_contract
            and gas_temperature is not None
            and gas_temperature > 0.0
            and isinstance(event_models, list)
            and event_models
            and len(event_models) == len(set(event_models))
            and set(event_models).issubset(allowed_event_models)
        )
    elif source == "propagator":
        required_contract = {
            **common_contract,
            "estimator": ("same_discrete_finite_temperature_elastic_sg_operator"),
            "operator": ("finite_volume_scharfetter_gummel_elastic_energy_generator"),
            "source_eedf": "same_solved_energy_angle_distribution",
            "neutral_thermal_motion_model": ("finite_temperature_fokker_planck"),
            "gas_temperature_terms_included": True,
            "uncertainty": "deterministic_kinetic_solver",
        }
        candidate = dict(physics_contract)
        gas_temperature = _float_or_none(candidate.pop("gas_temperature_K", None))
        contract_valid = bool(
            candidate == required_contract
            and gas_temperature is not None
            and gas_temperature > 0.0
        )
    else:  # The parser rejects this before bundle validation.
        contract_valid = False
    if not contract_valid:
        raise GecCcpWorkflowError(
            "external elastic energy-loss physics_contract is "
            f"incompatible with source={source}"
        )


def _elastic_statistical_evidence(
    *,
    source: str,
    rows: list[dict[str, Any]],
    coefficients: list[float],
    quality_thresholds: dict[str, Any],
) -> dict[str, Any]:
    if source != "monte_carlo":
        return {
            "kind": "deterministic_kinetic_solver",
            "uncertainty_applicability": "not_applicable",
        }
    elastic_rse_limit = float(quality_thresholds["major_rate_rse"])
    relative_standard_errors: list[float] = []
    valid_replicates: list[int] = []
    estimate_statuses: set[str] = set()
    for row, coefficient in zip(rows, coefficients, strict=True):
        standard_error = _required_float(
            row,
            "elastic_energy_loss_standard_error_eV_m3_s",
        )
        relative_standard_error = _required_float(
            row,
            "elastic_energy_loss_relative_standard_error",
        )
        ci_low = _required_float(
            row,
            "elastic_energy_loss_ci95_low_eV_m3_s",
        )
        ci_high = _required_float(
            row,
            "elastic_energy_loss_ci95_high_eV_m3_s",
        )
        critical = _required_float(row, "ci95_critical_value")
        try:
            replicate_count = int(str(row["valid_replicates"]))
            uncertainty_available = int(str(row["uncertainty_available"]))
        except (KeyError, TypeError, ValueError) as exc:
            raise GecCcpWorkflowError(
                "Monte Carlo elastic energy-loss uncertainty columns "
                "must contain integer replicate/availability values"
            ) from exc
        estimate_status = str(row.get("estimate_status", "")).strip()
        uncertainty_valid = bool(
            standard_error >= 0.0
            and relative_standard_error >= 0.0
            and ci_low <= coefficient <= ci_high
            and critical > 0.0
            and replicate_count >= 2
            and uncertainty_available == 1
            and relative_standard_error <= elastic_rse_limit
            and estimate_status
            and estimate_status != "unavailable"
            and math.isclose(
                relative_standard_error,
                standard_error / coefficient,
                rel_tol=1.0e-8,
                abs_tol=1.0e-15,
            )
            and math.isclose(
                ci_low,
                coefficient - critical * standard_error,
                rel_tol=1.0e-8,
                abs_tol=1.0e-12 * coefficient,
            )
            and math.isclose(
                ci_high,
                coefficient + critical * standard_error,
                rel_tol=1.0e-8,
                abs_tol=1.0e-12 * coefficient,
            )
        )
        if not uncertainty_valid:
            raise GecCcpWorkflowError(
                "Monte Carlo elastic energy-loss table lacks a "
                "self-consistent independent-replica CI/RSE contract"
            )
        relative_standard_errors.append(relative_standard_error)
        valid_replicates.append(replicate_count)
        estimate_statuses.add(estimate_status)
    return {
        "kind": "independent_replica_student_t_95",
        "uncertainty_available_at_every_anchor": True,
        "minimum_valid_replicates": min(valid_replicates),
        "maximum_relative_standard_error": max(relative_standard_errors),
        "relative_standard_error_limit": elastic_rse_limit,
        "relative_standard_error_qualified": True,
        "estimate_statuses": sorted(estimate_statuses),
    }
