"""Fail-closed bundle and physical-context validation for positive column."""

from __future__ import annotations

from electron_swarm.core.physical_context import PHYSICAL_CONTEXT_SCHEMA

import csv
from dataclasses import asdict
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from typing import TYPE_CHECKING

from swarm_workflow.comsol.input.bundle_selection import validate_bundle_selection
from swarm_workflow.comsol.input.export.contracts import TABLE_ROLES, FORMAT_VERSION
from swarm_workflow.comsol.input.validation import (
    ComsolInputArtifactError,
    validate_manifest_tables,
)
from swarm_workflow.selection import ClosureSelectionError

from .contracts import BundleValidationSummary, PositiveColumnWorkflowError
from .lookup import LookupValidationError, read_lookup

if TYPE_CHECKING:
    from .config import ComsolModelMapping

REQUIRED_TABLES = (
    "mean_energy_vs_en.csv",
    "transport_vs_en.csv",
    "transport_vs_mean_energy.csv",
    "rates_vs_mean_energy.csv",
    "quality.csv",
)
TRANSPORT_COLUMNS = (
    "reduced_mobility_m2_V_s_m3",
    "reduced_diffusion_L_m2_s_m3",
    "reduced_electron_energy_mobility_m2_V_s_m3",
    "reduced_electron_energy_diffusion_m2_s_m3",
)


def validate_positive_column_bundle(
    bundle_path: str | Path,
    *,
    mapping: ComsolModelMapping | None = None,
    gas_temperature_K: float | None = None,
) -> BundleValidationSummary:
    """Validate portable evidence before preparation or any licensed execution."""

    bundle = Path(bundle_path).resolve()
    try:
        manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise PositiveColumnWorkflowError(
            f"cannot read COMSOL bundle manifest: {bundle}"
        ) from exc
    if (
        not isinstance(manifest, dict)
        or manifest.get("format_version") != FORMAT_VERSION
        or isinstance(manifest.get("format_version"), bool)
    ):
        raise PositiveColumnWorkflowError(
            "COMSOL bundle format_version must be 2; regenerate obsolete bundles"
        )
    if manifest.get("stage") != "export-comsol" or manifest.get("status") != "ok":
        raise PositiveColumnWorkflowError(
            "COMSOL bundle must have stage=export-comsol and status=ok"
        )
    source = manifest.get("source")
    if source != "two_term" or (
        mapping is not None and source != mapping.bundle.source
    ):
        raise PositiveColumnWorkflowError(
            "COMSOL bundle source does not match the full-transport mapping"
        )
    policy = manifest.get("source_policy")
    if (
        not isinstance(policy, dict)
        or policy.get("source") != source
        or policy.get("postprocess") != "none"
    ):
        raise PositiveColumnWorkflowError(
            "COMSOL bundle requires a pure-solver source_policy with postprocess=none"
        )
    hashes = manifest.get("hashes")
    for name in (
        "workflow_config_sha256",
        "calculation_sha256",
        "cross_sections_sha256",
        "source_database_sha256",
    ):
        if (
            not isinstance(hashes, dict)
            or not isinstance(hashes.get(name), str)
            or not re.fullmatch(r"[0-9a-f]{64}", hashes[name])
        ):
            raise PositiveColumnWorkflowError(
                f"COMSOL bundle lacks canonical source provenance hash: {name}"
            )
    context = manifest.get("physical_context")
    if (
        not isinstance(context, dict)
        or context.get("schema") != PHYSICAL_CONTEXT_SCHEMA
    ):
        raise PositiveColumnWorkflowError(
            "COMSOL bundle lacks physical_context provenance"
        )
    temperature = _positive(
        context.get("gas_temperature_K"), "physical_context.gas_temperature_K"
    )
    _positive(
        context.get("gas_number_density_m3"), "physical_context.gas_number_density_m3"
    )
    field = context.get("field")
    if (
        not isinstance(field, dict)
        or field.get("type") != "dc"
        or not isinstance(field.get("magnetic_field"), dict)
    ):
        raise PositiveColumnWorkflowError(
            "positive-column bundle requires explicit homogeneous DC/B=0 context"
        )
    magnetic = field["magnetic_field"]
    if (
        not isinstance(magnetic.get("enabled"), bool)
        or isinstance(magnetic.get("B_T"), bool)
        or not isinstance(magnetic.get("B_T"), (float, int))
        or not math.isfinite(magnetic["B_T"])
        or magnetic["B_T"] < 0
        or (magnetic["enabled"] and magnetic["B_T"] != 0.0)
    ):
        raise PositiveColumnWorkflowError(
            "positive-column bundle does not support magnetic swarm inputs"
        )
    for key in ("angular_scattering", "electron_electron", "ionization", "finite_k"):
        if not isinstance(context.get(key), dict):
            raise PositiveColumnWorkflowError(f"physical_context.{key} is missing")
    if (
        context["electron_electron"].get("enabled") is not False
        or context["finite_k"].get("enabled") is not False
    ):
        raise PositiveColumnWorkflowError(
            "positive-column bundle requires unprocessed, local swarm inputs"
        )
    if context.get("high_energy_extrapolation") not in {"zero", "hold", "error"}:
        raise PositiveColumnWorkflowError(
            "physical_context.high_energy_extrapolation is missing"
        )
    mixture = manifest.get("mixture")
    species = mixture.get("species") if isinstance(mixture, dict) else None
    expected_species = (
        [asdict(item) for item in mapping.physical_context.gas_mixture]
        if mapping
        else [{"species": "Ar", "fraction": 1.0, "mass_amu": 39.948}]
    )
    if species != expected_species:
        raise PositiveColumnWorkflowError(
            "COMSOL bundle gas mixture does not match the pure-Ar model"
        )
    requested_temperature = (
        gas_temperature_K
        if gas_temperature_K is not None
        else (mapping.physical_context.gas_temperature_K if mapping else temperature)
    )
    if not math.isclose(
        temperature, requested_temperature, rel_tol=1e-12, abs_tol=1e-12
    ):
        raise PositiveColumnWorkflowError(
            "COMSOL model gas temperature does not match bundle physical_context"
        )
    summary = manifest.get("quality_summary")
    if (
        not isinstance(summary, dict)
        or summary.get("passed") is not True
        or summary.get("failed_points") != 0
        or not isinstance(summary.get("total_points"), int)
        or isinstance(summary.get("total_points"), bool)
        or summary["total_points"] <= 0
    ):
        raise PositiveColumnWorkflowError(
            "COMSOL bundle quality_summary must certify all table points"
        )
    try:
        artifacts = validate_manifest_tables(
            bundle, manifest, required_names=REQUIRED_TABLES, require_roles=True
        )
        validate_bundle_selection(bundle)
        for name in REQUIRED_TABLES:
            if manifest["tables"][name]["artifact_role"] != TABLE_ROLES[name]:
                raise PositiveColumnWorkflowError(
                    f"wrong artifact role for positive-column table: {name}"
                )
        mean = read_lookup(
            artifacts["mean_energy_vs_en.csv"].path, "E_over_N_Td", "mean_energy_eV"
        )
        for table in ("transport_vs_en.csv", "transport_vs_mean_energy.csv"):
            for column in TRANSPORT_COLUMNS:
                lookup = read_lookup(
                    artifacts[table].path,
                    "E_over_N_Td" if table.endswith("_en.csv") else "mean_energy_eV",
                    column,
                )
                if any(value <= 0 for value in lookup.y):
                    raise PositiveColumnWorkflowError(
                        f"COMSOL transport coefficient must be positive: {table}:{column}"
                    )
    except (
        ComsolInputArtifactError,
        ClosureSelectionError,
        LookupValidationError,
    ) as exc:
        raise PositiveColumnWorkflowError(str(exc)) from exc
    ranges = manifest.get("valid_ranges")
    bounds = ranges.get("E_over_N_Td") if isinstance(ranges, dict) else None
    if (
        not isinstance(bounds, list)
        or len(bounds) != 2
        or any(
            isinstance(v, bool)
            or not isinstance(v, (float, int))
            or not math.isfinite(v)
            or v <= 0
            for v in bounds
        )
        or bounds[0] > bounds[1]
    ):
        raise PositiveColumnWorkflowError(
            "COMSOL bundle valid_ranges.E_over_N_Td must be finite and ordered"
        )
    if not math.isclose(bounds[0], mean.x[0], rel_tol=1e-12) or not math.isclose(
        bounds[1], mean.x[-1], rel_tol=1e-12
    ):
        raise PositiveColumnWorkflowError(
            "COMSOL bundle valid range differs from coefficient tables"
        )
    with artifacts["quality.csv"].path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    if (
        not rows
        or len(rows) != summary["total_points"]
        or any(row.get("passed") != "1" for row in rows)
    ):
        raise PositiveColumnWorkflowError(
            "COMSOL bundle quality.csv does not certify every table point"
        )
    try:
        quality_fields = sorted(float(row["E_over_N_Td"]) for row in rows)
    except (KeyError, ValueError) as exc:
        raise PositiveColumnWorkflowError(
            "invalid quality.csv E/N coordinates"
        ) from exc
    if quality_fields != list(mean.x):
        raise PositiveColumnWorkflowError(
            "COMSOL bundle quality.csv coverage differs from coefficient tables"
        )
    return BundleValidationSummary(
        bundle,
        (float(bounds[0]), float(bounds[1])),
        True,
        REQUIRED_TABLES,
        sha256((bundle / "manifest.json").read_bytes()).hexdigest(),
        tuple(sorted((name, artifact.sha256) for name, artifact in artifacts.items())),
    )


def assert_bundle_unchanged(summary: BundleValidationSummary) -> None:
    """The execution stages must consume the exact generation prepared earlier."""
    for name, expected in (
        ("manifest.json", summary.manifest_sha256),
        *summary.table_sha256,
    ):
        path = summary.bundle_path / name
        if not path.is_file() or sha256(path.read_bytes()).hexdigest() != expected:
            raise PositiveColumnWorkflowError(
                f"positive-column prepared bundle changed: {name}"
            )


def _positive(value: object, name: str) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value <= 0
    ):
        raise PositiveColumnWorkflowError(f"{name} must be finite and positive")
    return float(value)
