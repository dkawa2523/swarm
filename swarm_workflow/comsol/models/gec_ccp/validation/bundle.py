"""Compose the fail-closed GEC-CCP bundle-validation pipeline."""

from __future__ import annotations

import hashlib
from typing import Any

from ..contracts import GecCcpMapping
from . import bundle_context as _context
from . import eedf_artifacts as _eedf_artifacts
from . import elastic_energy as _elastic_energy
from . import provenance as _provenance
from . import quality as _quality
from . import table_inputs as _table_inputs


def validate_gec_ccp_bundle(mapping: GecCcpMapping) -> dict[str, Any]:
    usage = _context.BundleUsage.from_mapping(mapping)
    manifest = _provenance.validate_manifest(mapping, uses_transport=usage.transport)
    tables = _table_inputs.validate_active_tables(mapping, usage, manifest)
    elastic_energy_loss = _elastic_energy.validate_elastic_energy_loss(
        mapping, usage, manifest
    )
    upstream_eedf_evidence = _eedf_artifacts.validate_upstream_eedf_evidence(
        mapping, usage, manifest, tables
    )
    function_eedf = _eedf_artifacts.validate_active_function_eedf(
        mapping, usage, manifest, tables
    )
    quality_audit = _quality.validate_bundle_quality(mapping, usage, manifest, tables)
    return {
        "path": str(mapping.bundle.path),
        "manifest_path": str(manifest.path),
        "manifest_sha256": hashlib.sha256(manifest.path.read_bytes()).hexdigest(),
        "source": manifest.source,
        "expected_source": mapping.bundle.expected_source,
        "solver_qualification": manifest.solver_qualification,
        "target_qualification": manifest.target_qualification,
        "mc_transport_estimator_schema": manifest.mc_estimator_schema,
        "mc_eedf_estimator_schema": manifest.mc_eedf_estimator_schema,
        "expected_mc_transport_estimator_schema": (
            mapping.bundle.expected_mc_transport_estimator_schema
        ),
        "mc_tail_estimator_schema": manifest.mc_tail_estimator_schema,
        "expected_mc_tail_estimator_schema": (
            mapping.bundle.expected_mc_tail_estimator_schema
        ),
        "mc_qualification_profile": manifest.source_policy.get("qualification_profile"),
        "expected_mc_qualification_profile": manifest.expected_mc_profile,
        "hashes": {
            name: manifest.hashes[name] for name in _context.GEC_PROVENANCE_HASH_KEYS
        },
        "mc_sampling_plan_sha256": (
            str(manifest.hashes["mc_sampling_plan_sha256"])
            if manifest.source == "monte_carlo"
            else None
        ),
        "mc_solver_source_sha256": (
            str(manifest.hashes["mc_solver_source_sha256"])
            if manifest.source == "monte_carlo"
            else None
        ),
        "quality_thresholds": manifest.quality_thresholds,
        "independent_quality_audit": quality_audit,
        "source_policy": manifest.source_policy,
        "two_term_transport_kernel": tables.two_term_transport_kernel,
        "restricted_scalar_diffusion_audit": tables.scalar_diffusion_audit,
        "valid_ranges": manifest.manifest.get("valid_ranges", {}),
        "artifact_verification": {
            "status": "passed",
            "verified_files": len(manifest.artifact_digests),
            "sha256": dict(sorted(manifest.artifact_digests.items())),
        },
        "function_eedf": function_eedf,
        "elastic_energy_loss": elastic_energy_loss,
        "upstream_eedf_evidence": upstream_eedf_evidence,
        "eedf_artifacts": manifest.manifest.get("eedf_artifacts"),
    }


__all__ = ["validate_gec_ccp_bundle"]
