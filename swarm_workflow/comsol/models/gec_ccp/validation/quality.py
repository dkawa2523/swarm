"""Apply final statistical and reaction-rate gates to a GEC-CCP bundle."""

from __future__ import annotations

from typing import Any

from ..contracts import GecCcpMapping, GecCcpWorkflowError
from ..data import _read_csv
from .bundle_context import BundleUsage, ManifestEvidence, TableEvidence
from .mc_quality import independent_bundle_quality_audit
from .table_inputs import nonnegative_numeric_values


def validate_bundle_quality(
    mapping: GecCcpMapping,
    usage: BundleUsage,
    manifest: ManifestEvidence,
    table: TableEvidence,
) -> dict[str, Any]:
    quality_rows = _read_csv(mapping.bundle.path / "quality.csv")
    audit = independent_bundle_quality_audit(
        quality_rows,
        source=manifest.source,
        thresholds=manifest.quality_thresholds,
        mc_sampling_plan=manifest.mc_sampling_plan,
        mc_estimator_schema=manifest.mc_estimator_schema,
        source_policy=manifest.source_policy,
        transport_rows=(
            table.transport_rows
            if manifest.source == "monte_carlo" and usage.transport
            else None
        ),
        valid_ranges=(
            manifest.manifest.get("valid_ranges")
            if manifest.source == "monte_carlo"
            and isinstance(manifest.manifest.get("valid_ranges"), dict)
            else None
        ),
    )
    failed = [str(row["E_over_N_Td"]) for row in audit["rows"] if not row["passed"]]
    if not audit["passed"]:
        raise GecCcpWorkflowError(
            "bundle has failed quality points: " + ", ".join(failed)
        )
    for reaction in mapping.reactions if usage.rates or usage.function_eedf else ():
        if reaction.process_type == "elastic" and usage.external_elastic_loss:
            continue
        selected = [
            row
            for row in table.rate_rows
            if row.get("process_type") == reaction.process_type
        ]
        nonnegative_numeric_values(selected, "rate_coefficient_m3_s")
    return audit
