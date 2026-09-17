"""Validate active GEC-CCP coefficient-table contracts and values."""

from __future__ import annotations

import math
from typing import Any, Iterable

from swarm_workflow.tables.contracts import ZERO_EVENT_CONFIDENCE

from ..closure import (
    _active_transport_functions,
    _uses_hybrid_einstein_transport,
    _validate_restricted_closure_policy,
)
from ..contracts import (
    FUNCTION_EEDF_TABLE,
    GEC_RESTRICTED_TRANSPORT_CLOSURE,
    GecCcpMapping,
    GecCcpWorkflowError,
)
from ..data import _read_csv
from . import bundle_guards as _bundle_guards
from .bundle_context import BundleUsage, ManifestEvidence, TableEvidence
from .mc_bundle import (
    REQUIRED_MC_QUALITY_COLUMNS as _REQUIRED_MC_QUALITY_COLUMNS,
)


REQUIRED_TABLE_COLUMNS = {
    "transport_vs_mean_energy.csv": (
        "mean_energy_eV",
        "E_over_N_Td",
        "reduced_mobility_m2_V_s_m3",
        "reduced_diffusion_L_m2_s_m3",
        "reduced_diffusion_T_m2_s_m3",
        "reduced_electron_energy_mobility_m2_V_s_m3",
        "reduced_electron_energy_diffusion_m2_s_m3",
        "reduced_electron_energy_diffusion_L_m2_s_m3",
        "reduced_electron_energy_diffusion_T_m2_s_m3",
    ),
    "rates_vs_mean_energy.csv": (
        "mean_energy_eV",
        "process_type",
        "rate_coefficient_m3_s",
    ),
    "quality.csv": (
        "E_over_N_Td",
        "passed",
        "eedf_normalization_error",
    ),
    "elastic_energy_loss_vs_mean_energy.csv": (
        "mean_energy_eV",
        "E_over_N_Td",
        "E_over_N_V_m2",
        "elastic_energy_loss_rate_coefficient_eV_m3_s",
        "elastic_energy_loss_standard_error_eV_m3_s",
        "elastic_energy_loss_relative_standard_error",
        "elastic_energy_loss_ci95_low_eV_m3_s",
        "elastic_energy_loss_ci95_high_eV_m3_s",
        "ci95_critical_value",
        "estimate_status",
        "valid_replicates",
        "uncertainty_available",
    ),
}
REQUIRED_RATE_EVIDENCE_COLUMNS = (
    "E_over_N_Td",
    "mean_energy_eV",
    "species",
    "process",
    "process_type",
    "rate_coefficient_mean_m3_s",
    "rate_coefficient_ci95_low_m3_s",
    "rate_coefficient_ci95_high_m3_s",
    "estimate_status",
    "uncertainty_available",
    "pooled_event_count",
    "pooled_target_exposure_s_m3",
    "pooled_all_zero_upper_95_m3_s",
    "pooled_zero_event_status",
)


def validate_active_tables(
    mapping: GecCcpMapping,
    usage: BundleUsage,
    manifest: ManifestEvidence,
) -> TableEvidence:
    required = {"quality.csv": REQUIRED_TABLE_COLUMNS["quality.csv"]}
    if manifest.source == "monte_carlo":
        required["quality.csv"] += _REQUIRED_MC_QUALITY_COLUMNS
    if usage.transport:
        required["transport_vs_mean_energy.csv"] = REQUIRED_TABLE_COLUMNS[
            "transport_vs_mean_energy.csv"
        ]
    if usage.rates or usage.function_eedf:
        required["rates_vs_mean_energy.csv"] = REQUIRED_TABLE_COLUMNS[
            "rates_vs_mean_energy.csv"
        ]
    if usage.external_elastic_loss:
        required[_bundle_guards.GEC_ELASTIC_ENERGY_LOSS_TABLE] = REQUIRED_TABLE_COLUMNS[
            _bundle_guards.GEC_ELASTIC_ENERGY_LOSS_TABLE
        ]
    source_eedf_table = FUNCTION_EEDF_TABLE
    if usage.external_elastic_loss or usage.function_eedf:
        required[source_eedf_table] = (
            "electron_energy_eV",
            "mean_energy_eV",
            "eepf_eV_m32",
        )
    for name, columns in required.items():
        entry = manifest.tables.get(name)
        listed = entry.get("columns", []) if isinstance(entry, dict) else []
        missing = [column for column in columns if column not in listed]
        if missing:
            raise GecCcpWorkflowError(
                f"bundle table {name} missing columns: {', '.join(missing)}"
            )
    roles = _active_artifact_roles(mapping, usage, source_eedf_table)
    for name in required:
        if manifest.tables[name].get("artifact_role") != roles[name]:
            raise GecCcpWorkflowError(f"bundle table has invalid artifact role: {name}")
    if manifest.source == "monte_carlo":
        _validate_mc_rate_evidence(manifest.tables)
    rate_rows = (
        _read_csv(mapping.bundle.path / "rates_vs_mean_energy.csv")
        if usage.rates or usage.function_eedf
        else []
    )
    _validate_reaction_coverage(mapping, usage, rate_rows)
    monotonic = manifest.manifest.get("monotonicity", {}).get(
        "mean_energy_strictly_monotonic"
    )
    if (
        usage.transport or usage.rates or usage.function_eedf
    ) and monotonic is not True:
        raise GecCcpWorkflowError(
            "GEC CCP lookup requires strictly monotonic mean energy"
        )
    transport_rows = (
        _read_csv(mapping.bundle.path / "transport_vs_mean_energy.csv")
        if usage.transport
        else []
    )
    kernel = (
        _bundle_guards._validate_two_term_temporal_growth_transport_contract(
            mapping, manifest.source_policy, transport_rows
        )
        if usage.transport
        else None
    )
    if usage.transport:
        positive_numeric_values(transport_rows, "mean_energy_eV")
    _validate_transport_values(mapping, manifest.source, transport_rows)
    scalar_diffusion = _validate_restricted_closure_policy(
        mapping,
        source=manifest.source,
        transport_rows=transport_rows if usage.transport else None,
    )
    return TableEvidence(
        source_eedf_table=source_eedf_table,
        expected_active_roles=roles,
        rate_rows=rate_rows,
        transport_rows=transport_rows,
        two_term_transport_kernel=kernel,
        scalar_diffusion_audit=scalar_diffusion,
    )


def _active_artifact_roles(
    mapping: GecCcpMapping,
    usage: BundleUsage,
    source_eedf_table: str,
) -> dict[str, str]:
    roles = {
        "quality.csv": "swarm_quality_audit",
        "transport_vs_mean_energy.csv": "canonical_comsol_coefficient_input",
        "rates_vs_mean_energy.csv": "canonical_comsol_coefficient_input",
        _bundle_guards.GEC_ELASTIC_ENERGY_LOSS_TABLE: (
            "canonical_comsol_elastic_energy_loss_input"
        ),
    }
    if usage.function_eedf and mapping.closure.function_eedf is not None:
        roles[mapping.closure.function_eedf.table] = (
            "canonical_comsol_function_eedf_input"
        )
    if usage.external_elastic_loss:
        roles[source_eedf_table] = "canonical_comsol_function_eedf_input"
    return roles


def _validate_mc_rate_evidence(tables: dict[str, dict[str, Any]]) -> None:
    evidence = tables.get("rate_evidence.csv")
    columns = evidence.get("columns", []) if isinstance(evidence, dict) else []
    expected_statistics = {
        "replicate_interval": "two_sided_student_t_95",
        "zero_event_confidence": ZERO_EVENT_CONFIDENCE,
        "zero_event_upper_bound": ("poisson_zero_count_over_pooled_target_exposure"),
    }
    missing = [
        column for column in REQUIRED_RATE_EVIDENCE_COLUMNS if column not in columns
    ]
    if (
        missing
        or not isinstance(evidence, dict)
        or evidence.get("artifact_role") != "raw_swarm_rate_evidence"
        or evidence.get("statistics") != expected_statistics
    ):
        raise GecCcpWorkflowError(
            "pure Monte Carlo GEC input requires canonical raw trajectory "
            "rate_evidence.csv with no substituted rates"
        )


def _validate_reaction_coverage(
    mapping: GecCcpMapping,
    usage: BundleUsage,
    rate_rows: list[dict[str, str]],
) -> None:
    if not (usage.rates or usage.function_eedf):
        return
    present = {row.get("process_type") for row in rate_rows}
    missing = [
        reaction.process_type
        for reaction in mapping.reactions
        if not (reaction.process_type == "elastic" and usage.external_elastic_loss)
        if reaction.process_type not in present
    ]
    if missing:
        raise GecCcpWorkflowError(
            "bundle lacks mapped reaction process types: " + ", ".join(missing)
        )


def _validate_transport_values(
    mapping: GecCcpMapping,
    source: str,
    rows: list[dict[str, str]],
) -> None:
    if mapping.closure.electron_transport == "swarm_mobility_einstein":
        required = ("reduced_mobility_m2_V_s_m3",)
    elif _uses_hybrid_einstein_transport(mapping.closure):
        required = tuple(
            item[2]
            for item in _active_transport_functions(mapping.closure, source=source)
        )
    elif mapping.closure.electron_transport == GEC_RESTRICTED_TRANSPORT_CLOSURE:
        required = tuple(
            item[2]
            for item in _active_transport_functions(mapping.closure, source=source)
        )
        if source == "two_term":
            required += (
                "reduced_diffusion_T_m2_s_m3",
                "reduced_electron_energy_diffusion_L_m2_s_m3",
                "reduced_electron_energy_diffusion_T_m2_s_m3",
            )
    else:
        required = ()
    for column in required:
        positive_numeric_values(rows, column)


def positive_numeric_values(rows: Iterable[dict[str, str]], column: str) -> list[float]:
    values = nonnegative_numeric_values(rows, column)
    if any(value <= 0.0 for value in values):
        raise GecCcpWorkflowError(f"bundle column {column} must be positive")
    return values


def nonnegative_numeric_values(
    rows: Iterable[dict[str, str]], column: str
) -> list[float]:
    values: list[float] = []
    for row in rows:
        try:
            value = float(row[column])
        except (KeyError, TypeError, ValueError) as exc:
            raise GecCcpWorkflowError(
                f"bundle column {column} contains a nonnumeric value"
            ) from exc
        if not math.isfinite(value) or value < 0.0:
            raise GecCcpWorkflowError(
                f"bundle column {column} contains an invalid value"
            )
        values.append(value)
    if not values:
        raise GecCcpWorkflowError(f"bundle column {column} has no values")
    return values
