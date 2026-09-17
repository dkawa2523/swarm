"""Typed state passed through the GEC-CCP bundle-validation pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ..closure import (
    _uses_external_elastic_energy_loss,
    _uses_external_rates,
    _uses_function_eedf,
    _uses_transport,
)
from ..contracts import GecCcpMapping


GEC_PROVENANCE_HASH_KEYS = (
    "workflow_config_sha256",
    "base_config_sha256",
    "cross_sections_sha256",
)
GEC_OPTIONAL_PROVENANCE_HASH_KEYS = (
    "mc_sampling_plan_json",
    "mc_sampling_plan_sha256",
    "mc_transport_estimator_schema_version",
    "mc_eedf_estimator_schema_version",
    "mc_tail_estimator_schema_version",
    "mc_seed_derivation_schema_version",
    "mc_solver_source_sha256",
    "propagator_solver_source_sha256",
    "source_database_sha256",
)


@dataclass(frozen=True, slots=True)
class BundleUsage:
    transport: bool
    rates: bool
    function_eedf: bool
    external_elastic_loss: bool

    @classmethod
    def from_mapping(cls, mapping: GecCcpMapping) -> BundleUsage:
        return cls(
            transport=_uses_transport(mapping.closure),
            rates=_uses_external_rates(mapping.closure),
            function_eedf=_uses_function_eedf(mapping.closure),
            external_elastic_loss=_uses_external_elastic_energy_loss(mapping.closure),
        )


@dataclass(frozen=True, slots=True)
class ManifestEvidence:
    path: Path
    manifest: dict[str, Any]
    source: str
    solver_qualification: dict[str, Any] | None
    target_qualification: dict[str, Any] | None
    source_policy: dict[str, Any]
    expected_mc_profile: str | None
    hashes: dict[str, Any]
    mc_sampling_plan: list[dict[str, Any]] | None
    mc_estimator_schema: str | None
    mc_eedf_estimator_schema: str | None
    mc_tail_estimator_schema: str | None
    quality_thresholds: dict[str, Any]
    tables: dict[str, dict[str, Any]]
    artifact_digests: dict[str, str]


@dataclass(frozen=True, slots=True)
class TableEvidence:
    source_eedf_table: str
    expected_active_roles: dict[str, str]
    rate_rows: list[dict[str, str]]
    transport_rows: list[dict[str, str]]
    two_term_transport_kernel: dict[str, Any] | None
    scalar_diffusion_audit: dict[str, Any] | None
