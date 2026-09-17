"""Validate projected and active Function-EEDF bundle artifacts."""

from __future__ import annotations

import math
from typing import Any

from swarm_workflow.comsol.input.function_eedf import ComsolEedfImportContract

from ..contracts import GecCcpMapping, GecCcpWorkflowError
from ..data import _float_or_none
from .bundle_context import BundleUsage, ManifestEvidence, TableEvidence
from .eedf import (
    _audit_upstream_eedf_artifact as _audit_eedf_artifact,
    _native_function_eedf_moment_audit as _audit_native_eedf_moments,
    _read_function_eedf as _load_function_eedf_for_preflight,
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
        expected_columns = [
            "electron_energy_eV",
            "mean_energy_eV",
            "eepf_eV_m32",
        ]
        expected_import = ComsolEedfImportContract().import_settings()
        if (
            entry.get("canonical_comsol_input") is not True
            or entry.get("format") != "csv"
            or entry.get("representation")
            != "physical_2d_adaptive_moment_rate_projected"
            or entry.get("argument_order") != ["electron_energy_eV", "mean_energy_eV"]
            or entry.get("columns") != expected_columns
            or entry.get("grid_axis_order") != ["electron_energy_eV", "mean_energy_eV"]
            or entry.get("row_order") != "mean_energy_major_then_electron_energy"
        ):
            raise GecCcpWorkflowError(
                "Function-EEDF active artifact must be the canonical physical "
                "adaptive linear moment/rate projection"
            )
        comsol_import = entry.get("comsol_import")
        if not isinstance(comsol_import, dict) or any(
            comsol_import.get(name) != value for name, value in expected_import.items()
        ):
            raise GecCcpWorkflowError(
                "Function-EEDF manifest does not match the canonical "
                "two-argument COMSOL Interpolation contract"
            )
        grid = _load_function_eedf_for_preflight(function_path)
        if ComsolEedfImportContract.from_path(function_path).structure != "spreadsheet":
            raise GecCcpWorkflowError(
                "Function-EEDF file serialization differs from its import contract"
            )
        expected_shape = [
            len(grid.mean_energies_eV),
            len(grid.electron_energies_eV),
        ]
        if (
            entry.get("grid_shape") != expected_shape
            or entry.get("mean_energy_grid_points") != expected_shape[0]
            or entry.get("energy_grid_points") != expected_shape[1]
        ):
            raise GecCcpWorkflowError(
                "Function-EEDF manifest grid dimensions disagree with the "
                "active structured grid"
            )
        projected_audit = _audit_native_eedf_moments(grid)
        norm_error = _float_or_none(entry.get("projected_normalization_error_max"))
        mean_error = _float_or_none(
            entry.get("projected_mean_energy_relative_error_max")
        )
        minimum = _float_or_none(entry.get("projected_nonnegative_minimum"))
        if (
            norm_error is None
            or norm_error > 1.0e-8
            or projected_audit["normalization_error_max"] > 1.0e-8
        ):
            raise GecCcpWorkflowError(
                "projected Function-EEDF normalization audit is missing or exceeds 1e-8"
            )
        if (
            mean_error is None
            or mean_error > 1.0e-8
            or projected_audit["mean_energy_relative_error_max"] > 1.0e-8
        ):
            raise GecCcpWorkflowError(
                "projected Function-EEDF mean-energy audit is missing or exceeds 1e-8"
            )
        if (
            minimum is None
            or minimum < 0.0
            or projected_audit["nonnegative_minimum"] < 0.0
        ):
            raise GecCcpWorkflowError(
                "projected Function-EEDF nonnegativity audit is missing or failed"
            )
        for reported, actual, name in (
            (
                norm_error,
                projected_audit["normalization_error_max"],
                "normalization",
            ),
            (
                mean_error,
                projected_audit["mean_energy_relative_error_max"],
                "mean energy",
            ),
            (
                minimum,
                projected_audit["nonnegative_minimum"],
                "nonnegativity",
            ),
        ):
            if not math.isclose(reported, actual, rel_tol=1.0e-9, abs_tol=1.0e-13):
                raise GecCcpWorkflowError(
                    f"Function-EEDF manifest {name} audit disagrees with "
                    "the active structured grid"
                )
        energy_range = [
            float(grid.electron_energies_eV[0]),
            float(grid.electron_energies_eV[-1]),
        ]
        mean_range = [
            float(grid.mean_energies_eV[0]),
            float(grid.mean_energies_eV[-1]),
        ]
        if (
            entry.get("electron_energy_range_eV") != energy_range
            or entry.get("mean_energy_range_eV") != mean_range
        ):
            raise GecCcpWorkflowError(
                "Function-EEDF manifest ranges disagree with the active structured grid"
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
    return function_eedf_summary
