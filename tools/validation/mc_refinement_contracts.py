"""Shared contracts for Monte Carlo population-refinement evaluation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from swarm_workflow.plots.eedf_contracts import EedfCase


SCHEMA = "swarm.validation.mc_population_refinement.v2"
DEFAULT_MIXTURES = ("ar_o2", "ar_n2", "ar_cl2")
DEFAULT_POPULATIONS = (512, 1024, 2048)
REQUIRED_ANCHORS_TD = (1.0, 10.0, 100.0, 1000.0)
LOW_FIELD_TAIL_ANCHORS_TD = (1.0, 10.0)
HIGH_FIELD_TD = 1000.0
TAIL_THRESHOLDS_EV = (80.0, 100.0, 120.0)
MAX_REPRESENTATIONAL_MERGE_MASS = 1.0e-14
MAX_REPRESENTATIONAL_MOMENT_CHANGE_EV = 1.0e-12
MAX_AGGREGATE_RECONSTRUCTION_TV = 1.0e-12
MAX_AGGREGATE_MEAN_RELATIVE_ERROR = 1.0e-12
REFINEMENT_CONTEXT_METADATA = (
    "base_config_sha256",
    "cross_section_files_json",
    "cross_sections_sha256",
    "mc_base_seed",
    "mc_eedf_estimator_schema_version",
    "mc_tail_estimator_schema_version",
    "mc_transport_estimator_schema_version",
    "physical_context_json",
    "quality_thresholds_json",
)


class RefinementEvaluationError(RuntimeError):
    """Raised when MC refinement evidence is incomplete or inconsistent."""


@dataclass(frozen=True, slots=True)
class DatabaseInput:
    mixture: str
    population: int
    path: Path


def csv_bool(value: bool | None) -> str | int:
    return "" if value is None else int(value)


@dataclass(frozen=True, slots=True)
class QualificationStatus:
    profile: str | None
    overall_passed: bool | None
    aggregate_quality_passed: bool | None
    active_closure_passed: bool | None
    solver_transport_passed: bool | None
    explicit_eedf_passed: bool | None
    eedf_normalization_passed: bool | None

    def csv_values(self) -> dict[str, object]:
        return {
            "qualification_profile": self.profile,
            "qualification_overall_passed": csv_bool(self.overall_passed),
            "qualification_aggregate_quality_passed": csv_bool(
                self.aggregate_quality_passed
            ),
            "qualification_active_closure_passed": csv_bool(
                self.active_closure_passed
            ),
            "qualification_solver_transport_passed": csv_bool(
                self.solver_transport_passed
            ),
            "qualification_explicit_eedf_passed": csv_bool(
                self.explicit_eedf_passed
            ),
            "qualification_eedf_normalization_passed": csv_bool(
                self.eedf_normalization_passed
            ),
        }


@dataclass(slots=True)
class LoadedDatabase:
    source: DatabaseInput
    source_sha256: str
    fractions: dict[str, float]
    tail_rows: list[dict[str, object]]
    aggregate_high_field: EedfCase
    final_replica_cases: dict[float, list[tuple[int, EedfCase]]]
    projection_events: list[dict[str, object]]
    metadata: dict[str, str]
    refinement_context_sha256: str
    sampling_contract_sha256: str
    replicas_by_anchor: dict[float, int]
    campaign_limits: dict[str, int]


__all__ = [
    "DEFAULT_MIXTURES",
    "DEFAULT_POPULATIONS",
    "DatabaseInput",
    "HIGH_FIELD_TD",
    "LOW_FIELD_TAIL_ANCHORS_TD",
    "LoadedDatabase",
    "MAX_AGGREGATE_MEAN_RELATIVE_ERROR",
    "MAX_AGGREGATE_RECONSTRUCTION_TV",
    "MAX_REPRESENTATIONAL_MERGE_MASS",
    "MAX_REPRESENTATIONAL_MOMENT_CHANGE_EV",
    "QualificationStatus",
    "REFINEMENT_CONTEXT_METADATA",
    "REQUIRED_ANCHORS_TD",
    "RefinementEvaluationError",
    "SCHEMA",
    "TAIL_THRESHOLDS_EV",
]
