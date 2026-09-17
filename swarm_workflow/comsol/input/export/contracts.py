"""Contracts and artifact names for canonical COMSOL bundle export."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from swarm_workflow.selection import (
    MC_QUALIFICATION_TABLE,
)
from swarm_workflow.tables.contracts import (
    COLLISION_RATE_KERNEL_TABLE,
    ELASTIC_ENERGY_LOSS_TABLE,
    RATE_EVIDENCE_TABLE,
)


FORMAT_VERSION = 1


BASE_TABLES = (
    "mean_energy_vs_en.csv",
    "transport_vs_en.csv",
    "rates_vs_en.csv",
    "townsend_vs_en.csv",
    "quality.csv",
)


OPTIONAL_TABLES = (
    "energy_loss.csv",
    MC_QUALIFICATION_TABLE,
    "rates_vs_mean_energy.csv",
    "transport_vs_mean_energy.csv",
    "eedf.csv",
    RATE_EVIDENCE_TABLE,
    ELASTIC_ENERGY_LOSS_TABLE,
    COLLISION_RATE_KERNEL_TABLE,
)


COMSOL_FUNCTION_TABLES = (
    (
        "comsol_functions/sw_meanE.csv",
        "mean_energy_vs_en.csv",
        "mean_energy_eV",
        "mean_energy_eV",
    ),
    (
        "comsol_functions/sw_muN.csv",
        "transport_vs_en.csv",
        "reduced_mobility_m2_V_s_m3",
        "reduced_mobility_m2_V_s_m3",
    ),
)


FUNCTION_EEDF_TABLE = "eedf_f0_vs_mean_energy.csv"


COMSOL_FUNCTION_EEDF_TABLE = "eedf_f0_comsol_2d.csv"


COPIED_TABLE_ROLES = {
    "energy_loss.csv": "swarm_coefficient_evidence",
    MC_QUALIFICATION_TABLE: "all_planned_mc_anchor_qualification",
    "mean_energy_vs_en.csv": "swarm_coefficient_evidence",
    "transport_vs_en.csv": "swarm_coefficient_evidence",
    "rates_vs_en.csv": "swarm_coefficient_evidence",
    "townsend_vs_en.csv": "swarm_coefficient_evidence",
    "quality.csv": "swarm_quality_audit",
    "rates_vs_mean_energy.csv": "canonical_comsol_coefficient_input",
    "transport_vs_mean_energy.csv": "canonical_comsol_coefficient_input",
    "eedf.csv": "raw_swarm_eedf_audit_and_canonicalization_source",
    RATE_EVIDENCE_TABLE: "raw_swarm_rate_evidence",
    ELASTIC_ENERGY_LOSS_TABLE: ("canonical_comsol_elastic_energy_loss_input"),
    COLLISION_RATE_KERNEL_TABLE: "canonical_eedf_projection_rate_kernels",
}


DERIVED_FUNCTION_ROLE = "derived_comsol_interpolation_input"


class ComsolExportError(RuntimeError):
    """Raised when built tables cannot be exported safely for COMSOL."""


@dataclass(frozen=True, slots=True)
class ComsolExportSummary:
    table_directory: Path
    output_directory: Path
    bundles_written: int
