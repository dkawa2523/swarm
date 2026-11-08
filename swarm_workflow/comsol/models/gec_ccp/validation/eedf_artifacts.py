"""Validate projected and active Function-EEDF bundle artifacts."""

from __future__ import annotations

from typing import Any

from swarm_workflow.tables.contracts import PROPAGATOR_SOLVER

from ..closure import _function_eedf_rate_reactions
from ..contracts import GecCcpMapping, GecCcpWorkflowError
from .bundle_context import BundleUsage, ManifestEvidence, TableEvidence
from .eedf import (
    _audit_upstream_eedf_artifact as _audit_eedf_artifact,
    _audit_function_eedf_source_rates,
    _read_function_eedf as _load_function_eedf_for_preflight,
    read_mph_reaction_cross_sections,
)


def validate_upstream_eedf_evidence(
    mapping: GecCcpMapping,
    usage: BundleUsage,
    manifest: ManifestEvidence,
    table: TableEvidence,
) -> dict[str, Any] | None:
    tables = manifest.tables
    source_eedf_table = table.source_eedf_table
    artifact_digests = manifest.artifact_digests
    uses_transport = usage.transport
    uses_rates = usage.rates
    uses_external_elastic_loss = usage.external_elastic_loss
    uses_function_eedf = usage.function_eedf
    eedf_entry = tables.get(source_eedf_table)
    upstream_eedf_evidence = (
        {
            "table": source_eedf_table,
            "table_sha256": artifact_digests.get(source_eedf_table),
            "bundle_artifact_role": eedf_entry.get("artifact_role"),
            "bundle_comsol_import_capable": eedf_entry.get("canonical_comsol_input"),
            "representation": eedf_entry.get("representation"),
            "mean_energy_range_eV": eedf_entry.get("mean_energy_range_eV"),
            "grid_shape": eedf_entry.get("grid_shape"),
            "active_in_comsol": uses_function_eedf,
            "binding_status": (
                "active_comsol_function_eedf"
                if uses_function_eedf
                else "inactive_evidence_only"
            ),
            "role": (
                "upstream_projected_eedf_evidence"
                if not uses_function_eedf
                else "active_comsol_function_eedf"
            ),
        }
        if isinstance(eedf_entry, dict)
        else None
    )
    if (
        uses_transport or uses_rates or uses_external_elastic_loss
    ) and not uses_function_eedf:
        if not isinstance(eedf_entry, dict):
            raise GecCcpWorkflowError(
                "external Swarm closure requires its projected EEDF evidence"
            )
        upstream_moment_audit = _audit_eedf_artifact(
            mapping.bundle.path / source_eedf_table,
            eedf_entry,
            expected_sha256=artifact_digests[source_eedf_table],
        )
        if upstream_eedf_evidence is None:
            raise GecCcpWorkflowError("projected EEDF evidence metadata is unavailable")
        upstream_eedf_evidence["recomputed_moment_audit"] = upstream_moment_audit
        upstream_eedf_evidence["evidence_verified_from_csv"] = True
    return upstream_eedf_evidence


def validate_active_function_eedf(
    mapping: GecCcpMapping,
    usage: BundleUsage,
    manifest: ManifestEvidence,
    table: TableEvidence,
) -> dict[str, Any] | None:
    uses_function_eedf = usage.function_eedf
    tables = manifest.tables
    artifact_digests = manifest.artifact_digests
    expected_active_roles = table.expected_active_roles
    function_eedf_summary: dict[str, Any] | None = None
    if uses_function_eedf:
        spec = mapping.closure.function_eedf
        if spec is None:  # Parser validation should make this unreachable.
            raise GecCcpWorkflowError("missing Function-EEDF specification")
        function_path = mapping.bundle.path / spec.table
        entry = tables.get(spec.table)
        if not isinstance(entry, dict) or not function_path.exists():
            raise GecCcpWorkflowError(f"bundle lacks Function-EEDF table: {spec.table}")
        if entry.get("artifact_role") != expected_active_roles[spec.table]:
            raise GecCcpWorkflowError(
                "Function-EEDF table is not the canonical COMSOL input"
            )
        grid = _load_function_eedf_for_preflight(function_path, entry)
        expected_shape = list(grid.values_eV_m32.shape)
        projected_audit = {
            "integration": "exact_piecewise_linear_energy_segments",
            "rows": expected_shape[0],
            "normalization_error_max": grid.normalization_error_max,
            "mean_energy_relative_error_max": grid.mean_energy_relative_error_max,
            "nonnegative_minimum": grid.minimum_value,
        }
        if (
            max(grid.normalization_error_max, grid.mean_energy_relative_error_max)
            > 1e-8
        ):
            raise GecCcpWorkflowError(
                "Function-EEDF normalization or mean-energy error exceeds 1e-8"
            )
        energy_range = [
            float(grid.electron_energies_eV[0]),
            float(grid.electron_energies_eV[-1]),
        ]
        mean_range = [float(grid.mean_energies_eV[0]), float(grid.mean_energies_eV[-1])]
        if (
            entry.get("argument_order") != ["electron_energy_eV", "mean_energy_eV"]
            or entry.get("electron_energy_range_eV") != energy_range
            or entry.get("mean_energy_range_eV") != mean_range
        ):
            raise GecCcpWorkflowError(
                "Function-EEDF manifest arguments or ranges disagree with the grid"
            )

        # Validate source-derived rates before running a model that integrates
        # its Function EEDF. This also rejects older inconsistent bundles.
        rate_consistency = None
        if mapping.closure.reaction_model == "function_eedf" and (
            manifest.source == PROPAGATOR_SOLVER
            or entry.get("source_rate_projection") is not None
        ):
            reactions = tuple(
                reaction for _, reaction in _function_eedf_rate_reactions(mapping)
            )
            cross_sections = read_mph_reaction_cross_sections(
                mapping.model.input_mph, reactions
            )
            rate_consistency = _audit_function_eedf_source_rates(
                mapping, grid, cross_sections, reactions, audit_csv_name=None
            )
            if not rate_consistency["passed"]:
                failed = ", ".join(
                    name
                    for name, result in rate_consistency["processes"].items()
                    if not result["passed"]
                )
                raise GecCcpWorkflowError(
                    f"Function-EEDF source rate mismatch before COMSOL ({failed}). Rebuild the bundle from the qualified tables."
                )
        actual_sha256 = artifact_digests[spec.table]
        function_eedf_summary = {
            "table": spec.table,
            "function_tag": spec.function_tag,
            "representation": spec.interpolation,
            "comsol_feature_operation": "Interpolation",
            "continuity": "C0",
            "table_sha256": actual_sha256,
            "electron_energy_range_eV": energy_range,
            "mean_energy_range_eV": mean_range,
            "grid_shape": expected_shape,
            "interpolation": "linear",
            "data_structure": "spreadsheet",
            "extrapolation": "constant",
            "projected_moment_audit": projected_audit,
        }
        if rate_consistency is not None:
            function_eedf_summary["source_rate_consistency"] = rate_consistency
    return function_eedf_summary
