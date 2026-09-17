"""Validate guarded GEC-CCP bundle contracts."""

from __future__ import annotations

import hashlib
import json
import math
import re
from typing import Any

from electron_swarm.solvers.two_term.transport import (
    TWO_TERM_TEMPORAL_GROWTH_EFFECTIVE_MOMENTUM_MODEL,
    TWO_TERM_TEMPORAL_GROWTH_TRANSPORT_SCHEMA_VERSION,
)

from swarm_workflow.comsol.models.gec_ccp.contracts import (
    GecCcpMapping,
    GecCcpWorkflowError,
)
from swarm_workflow.comsol.models.gec_ccp.data import (
    _lookup,
    _read_csv,
    _required_float,
)
from swarm_workflow.comsol.models.gec_ccp.validation.mc_quality import (
    independent_bundle_quality_audit,
)
from swarm_workflow.quality.policy import (
    parse_quality_thresholds,
    quality_thresholds_payload,
)
from swarm_workflow.tables.contracts import (
    TWO_TERM_MOMENTUM_CROSS_SECTION_FLOOR_M2,
    TWO_TERM_TEMPORAL_GROWTH_TRANSPORT_COLUMNS,
)


GEC_LOW_ENERGY_GUARD_RELATIVE_INFLUENCE_LIMIT = 1.0e-2
GEC_ELASTIC_ENERGY_LOSS_TABLE = "elastic_energy_loss_vs_mean_energy.csv"
GEC_ELASTIC_ENERGY_LOSS_COLUMN = "elastic_energy_loss_rate_coefficient_eV_m3_s"


def _validate_two_term_temporal_growth_transport_contract(
    mapping: GecCcpMapping,
    source_policy: dict[str, Any],
    transport_rows: list[dict[str, Any]],
) -> dict[str, Any] | None:
    """Validate the per-anchor PT eigenvalue used by two-term transport."""

    expected = mapping.bundle.expected_two_term_transport_kernel_schema
    raw = source_policy.get("two_term_transport_kernel")
    if raw is None:
        if expected is not None:
            raise GecCcpWorkflowError(
                "two-term bundle lacks the required temporal-growth "
                "transport-kernel provenance"
            )
        return None
    if mapping.bundle.expected_source != "two_term" or not isinstance(raw, dict):
        raise GecCcpWorkflowError(
            "two-term transport-kernel provenance is invalid for this source"
        )
    schema = raw.get("schema")
    if schema != TWO_TERM_TEMPORAL_GROWTH_TRANSPORT_SCHEMA_VERSION:
        raise GecCcpWorkflowError(
            f"unsupported bundle two-term transport-kernel schema: {schema}"
        )
    if expected is not None and schema != expected:
        raise GecCcpWorkflowError(
            "COMSOL bundle two-term transport-kernel schema mismatch: "
            f"expected {expected}, found {schema}"
        )
    required_contract = {
        "correction_applied_every_anchor": True,
        "effective_momentum_model": (TWO_TERM_TEMPORAL_GROWTH_EFFECTIVE_MOMENTUM_MODEL),
        "growth_frequency_source": "converged_temporal_growth_eigenvalue",
        "momentum_cross_section_floor_policy": (
            "sum_raw_process_cross_sections_then_single_floor"
        ),
        "minimum_momentum_cross_section_m2": (TWO_TERM_MOMENTUM_CROSS_SECTION_FLOOR_M2),
        "evidence_columns": list(TWO_TERM_TEMPORAL_GROWTH_TRANSPORT_COLUMNS),
    }
    if any(raw.get(name) != value for name, value in required_contract.items()):
        raise GecCcpWorkflowError(
            "two-term temporal-growth transport-kernel contract is incomplete"
        )
    if raw.get("anchor_points") != len(transport_rows) or len(transport_rows) < 2:
        raise GecCcpWorkflowError(
            "two-term temporal-growth evidence does not cover every transport anchor"
        )
    densities: list[float] = []
    growths: list[float] = []
    for row in transport_rows:
        density = _required_float(row, "gas_number_density_m3")
        growth = _required_float(row, "temporal_growth_frequency_s_inv")
        reduced = _required_float(row, "reduced_temporal_growth_frequency_m3_s")
        if density <= 0.0 or not math.isclose(
            reduced,
            growth / density,
            rel_tol=1.0e-10,
            abs_tol=1.0e-30,
        ):
            raise GecCcpWorkflowError(
                "two-term temporal-growth frequency/density evidence is "
                "nonphysical or inconsistent"
            )
        densities.append(density)
        growths.append(growth)
    reported_density = raw.get("gas_number_density_m3")
    reported_growth = raw.get("growth_frequency_s_inv")
    if reported_density != [min(densities), max(densities)] or reported_growth != [
        min(growths),
        max(growths),
    ]:
        raise GecCcpWorkflowError(
            "two-term temporal-growth source-policy ranges disagree with "
            "the transport table"
        )
    return dict(raw)


def _validate_low_energy_guard_bundle(
    mapping: GecCcpMapping,
    *,
    primary_cross_sections_sha256: str,
) -> dict[str, Any] | None:
    """Validate the independent two-term bundle used only for sensitivity.

    The guard never replaces a primary MC coefficient. It supplies a
    qualified alternative below the primary MC support so the converged
    application can bound the influence of constant endpoint continuation.
    """

    guard_root = mapping.run.low_energy_guard_bundle
    if mapping.run.support_policy == "strict":
        if guard_root is not None:
            raise GecCcpWorkflowError(
                "strict support policy must not declare a low-energy guard"
            )
        return None
    if guard_root is None:
        raise GecCcpWorkflowError("low-energy guard bundle is missing")
    manifest_path = guard_root / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise GecCcpWorkflowError(
            f"cannot read low-energy guard manifest: {manifest_path}"
        ) from exc
    source_policy = manifest.get("source_policy")
    hashes = manifest.get("hashes")
    if (
        manifest.get("status") != "ok"
        or manifest.get("source") != "two_term"
        or not isinstance(source_policy, dict)
        or source_policy.get("postprocess") != "none"
        or source_policy.get("field_type") != "dc"
        or not isinstance(hashes, dict)
        or hashes.get("cross_sections_sha256") != primary_cross_sections_sha256
        or manifest.get("quality_summary", {}).get("passed") is not True
    ):
        raise GecCcpWorkflowError(
            "low-energy guard must be a passed, unmodified, steady-DC "
            "two_term bundle with the same cross sections as the primary MC bundle"
        )
    tables = manifest.get("tables")
    required = (
        "quality.csv",
        "transport_vs_mean_energy.csv",
        "rates_vs_mean_energy.csv",
        GEC_ELASTIC_ENERGY_LOSS_TABLE,
    )
    if not isinstance(tables, dict) or any(name not in tables for name in required):
        raise GecCcpWorkflowError(
            "low-energy guard lacks a required quality or coefficient table"
        )
    artifact_sha256: dict[str, str] = {}
    for name in required:
        entry = tables[name]
        path = guard_root / name
        expected = entry.get("sha256") if isinstance(entry, dict) else None
        if (
            not path.is_file()
            or re.fullmatch(r"[0-9a-f]{64}", str(expected or "")) is None
        ):
            raise GecCcpWorkflowError(
                f"low-energy guard artifact metadata is invalid: {name}"
            )
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise GecCcpWorkflowError(
                f"low-energy guard artifact SHA-256 mismatch: {name}"
            )
        artifact_sha256[name] = actual

    raw_thresholds = manifest.get("quality_thresholds")
    if not isinstance(raw_thresholds, dict):
        raise GecCcpWorkflowError("low-energy guard lacks quality thresholds")
    thresholds = quality_thresholds_payload(parse_quality_thresholds(raw_thresholds))
    if raw_thresholds != thresholds:
        raise GecCcpWorkflowError(
            "low-energy guard quality thresholds are noncanonical"
        )
    quality = independent_bundle_quality_audit(
        _read_csv(guard_root / "quality.csv"),
        source="two_term",
        thresholds=thresholds,
        mc_sampling_plan=None,
        mc_estimator_schema=None,
    )
    if not quality["passed"]:
        raise GecCcpWorkflowError(
            "low-energy guard contains independently failed quality points"
        )
    failed_quality_points = sum(not bool(row.get("passed")) for row in quality["rows"])

    transport_rows = _read_csv(guard_root / "transport_vs_mean_energy.csv")
    transport_energy, transport_mobility = _lookup(
        transport_rows,
        "mean_energy_eV",
        "reduced_mobility_m2_V_s_m3",
    )
    if (
        len(transport_energy) < 4
        or any(
            right <= left for left, right in zip(transport_energy, transport_energy[1:])
        )
        or any(value <= 0.0 for value in transport_mobility)
    ):
        raise GecCcpWorkflowError(
            "low-energy guard mobility table is not positive and monotonic in mean energy"
        )
    rate_rows = _read_csv(guard_root / "rates_vs_mean_energy.csv")
    rate_supports: dict[str, list[float]] = {}
    for process_type in ("excitation", "ionization"):
        selected = [row for row in rate_rows if row.get("process_type") == process_type]
        energy, rate = _lookup(selected, "mean_energy_eV", "rate_coefficient_m3_s")
        if (
            len(energy) < 4
            or any(right <= left for left, right in zip(energy, energy[1:]))
            or any(value <= 0.0 for value in rate)
        ):
            raise GecCcpWorkflowError(
                f"low-energy guard {process_type} rate table is invalid"
            )
        rate_supports[process_type] = [float(energy[0]), float(energy[-1])]
    elastic_rows = _read_csv(guard_root / GEC_ELASTIC_ENERGY_LOSS_TABLE)
    elastic_energy, elastic_loss = _lookup(
        elastic_rows, "mean_energy_eV", GEC_ELASTIC_ENERGY_LOSS_COLUMN
    )
    if (
        len(elastic_energy) < 4
        or any(right <= left for left, right in zip(elastic_energy, elastic_energy[1:]))
        or any(value <= 0.0 for value in elastic_loss)
    ):
        raise GecCcpWorkflowError("low-energy guard elastic-loss table is invalid")
    return {
        "status": "qualified",
        "role": "frozen_field_low_energy_sensitivity_only",
        "active_comsol_input": False,
        "source": "two_term",
        "path": str(guard_root),
        "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "artifact_sha256": artifact_sha256,
        "cross_sections_sha256": primary_cross_sections_sha256,
        "quality_audit": {
            "passed": True,
            "decision": quality["decision"],
            "points": len(quality["rows"]),
            "failed_points": failed_quality_points,
        },
        "supports_mean_energy_eV": {
            "mobility": [float(transport_energy[0]), float(transport_energy[-1])],
            "rates": rate_supports,
            "elastic_energy_loss": [
                float(elastic_energy[0]),
                float(elastic_energy[-1]),
            ],
        },
        "relative_influence_limit": (GEC_LOW_ENERGY_GUARD_RELATIVE_INFLUENCE_LIMIT),
    }
