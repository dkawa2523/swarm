"""Validate canonical Swarm bundles for the Argon GEC-ICP closure."""

from __future__ import annotations

import csv
from hashlib import sha256
import json
import math
from pathlib import Path, PurePosixPath
import re
from typing import Any

import electron_swarm.solvers.monte_carlo.evidence as _mc_evidence
from electron_swarm.solvers.monte_carlo.source_identity import (
    monte_carlo_source_sha256,
)

from swarm_workflow.selection import (
    ANCHOR_FALLBACK_SCHEMA,
    ANCHOR_FALLBACK_SOURCE,
    ANCHOR_FALLBACK_STAGE,
    COMPOSITE_QUALITY_COLUMNS,
    COMPOSITE_QUALITY_SCHEMA,
    MC_QUALIFICATION_TABLE,
)
from swarm_workflow.comsol.input.export.contracts import (
    COMSOL_FUNCTION_EEDF_TABLE,
    FORMAT_VERSION,
)
from swarm_workflow.comsol.input.function_eedf import (
    COMSOL_EEDF_COLUMNS,
    ComsolEedfImportContract,
    FunctionEedfError,
    read_comsol_function_eedf_grid,
)
from swarm_workflow.quality.table import quality_table_schema
from swarm_workflow.quality.solver import (
    PROPAGATOR_TARGET_QUALIFICATION_SCHEMA,
    PropagatorTargetRequirement,
    SolverQualificationError,
    validate_copied_qualification,
    validate_propagator_core_qualification,
    validate_propagator_target_qualification,
)
from swarm_workflow.tables.contracts import (
    MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA,
    MC_SOLVER,
    PROPAGATOR_SOLVER,
    SOURCE_CHOICES,
)

from .contracts import GecIcpContractError
from .cross_sections import canonical_cross_sections_combined_sha256


_BOLTZMANN_J_K = 1.380649e-23
_ARGON_MASS_AMU = 39.948
_TRANSPORT_TABLE = "transport_vs_mean_energy.csv"
_QUALITY_TABLE = "quality.csv"
_ICP_SOURCES = (*SOURCE_CHOICES, ANCHOR_FALLBACK_SOURCE)
_COMPOSITE_EVIDENCE_ROLES = {
    "anchor_fallback_plan.json": "anchor_fallback_plan",
    "mc_solver_decision.json": "terminal_mc_policy_decision",
    "monte_carlo_qualification.csv": "monte_carlo_anchor_qualification",
    "monte_carlo_source_manifest.json": "monte_carlo_source_manifest",
    "two_term_quality.csv": "two_term_anchor_qualification",
    "two_term_source_manifest.json": "two_term_source_manifest",
}
_REQUIRED_PROVENANCE_HASHES = (
    "workflow_config_sha256",
    "base_config_sha256",
    "cross_sections_sha256",
)
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
_ICP_PROPAGATOR_TARGET_REQUIREMENT = PropagatorTargetRequirement(
    target="argon_gec_icp_restricted_local_mean_energy",
    fields_Td=(1500.0, 2500.0, 2744.2936035),
    operating_range_bracket_Td=(1500.0, 2500.0),
    table_support_cap_Td=2744.2936035,
    mixture_id=0,
    medium_grid=(300, 48),
    fine_grid=(600, 72),
    max_scalar_refinement_limit=0.01,
    max_eedf_weighted_L1_limit=0.02,
)


class GecIcpBundleError(ValueError):
    """Raised when a Swarm bundle cannot safely drive the GEC-ICP closure."""


def validate_gec_icp_bundle(
    bundle_path: str | Path,
    expected_source: str,
    expected_pressure_Pa: float = 2.66644,
    expected_temperature_K: float = 300.0,
    expected_transport_definition: str | None = None,
    expected_mc_qualification_profile: str | None = None,
) -> dict[str, Any]:
    """Validate one model-independent bundle for the restricted ICP closure."""

    root = Path(bundle_path).resolve()
    if not root.is_dir():
        raise GecIcpBundleError(f"GEC-ICP bundle directory does not exist: {root}")
    if expected_source not in _ICP_SOURCES:
        raise GecIcpBundleError(f"unsupported GEC-ICP bundle source: {expected_source}")
    expected_pressure = _positive_number(expected_pressure_Pa, "expected pressure")
    expected_temperature = _positive_number(
        expected_temperature_K, "expected temperature"
    )

    manifest_path = root / "manifest.json"
    manifest = _read_manifest(manifest_path)
    source, transport_definition, source_policy, hashes = _validate_manifest_identity(
        manifest,
        expected_source=expected_source,
        expected_transport_definition=expected_transport_definition,
    )
    pressure, temperature, mixture = _validate_physical_context(
        manifest,
        expected_pressure_Pa=expected_pressure,
        expected_temperature_K=expected_temperature,
    )
    tables, table_paths, digests = _verify_table_artifacts(root, manifest)
    evidence_paths, evidence_digests = _verify_composite_evidence(
        root, source, manifest
    )
    digests.update(evidence_digests)
    function_path, function_support = _validate_active_function_eedf(
        source, manifest, tables, table_paths
    )
    transport_path, transport_support = _validate_transport_table(
        manifest, tables, table_paths
    )
    if not all(
        math.isclose(left, right, rel_tol=1.0e-10, abs_tol=1.0e-12)
        for left, right in zip(function_support, transport_support, strict=True)
    ):
        raise GecIcpBundleError(
            "active Function-EEDF and mean-energy transport supports differ"
        )
    quality_path, quality_summary = _validate_quality(
        source, manifest, tables, table_paths
    )
    solver_qualification, target_qualification = _validate_solver_qualifications(
        source, root, manifest, hashes=hashes
    )
    mc_qualification = _validate_mc_qualification(
        source,
        manifest,
        source_policy,
        tables,
        table_paths,
        evidence_paths,
        expected_profile=expected_mc_qualification_profile,
    )

    selected_paths = {
        "function_eedf": function_path,
        "transport": transport_path,
        "quality": quality_path,
    }
    if source == MC_SOLVER:
        selected_paths["mc_qualification"] = table_paths[MC_QUALIFICATION_TABLE]
    if source == ANCHOR_FALLBACK_SOURCE:
        selected_paths.update(
            {f"evidence:{name}": path for name, path in evidence_paths.items()}
        )
    return {
        "path": str(root),
        "manifest_path": str(manifest_path),
        "manifest_sha256": _file_sha256(manifest_path),
        "source": source,
        "expected_source": expected_source,
        "mc_solver_source_sha256": (
            mc_qualification.get("solver_source_sha256")
            if mc_qualification is not None
            else None
        ),
        "transport_definition": transport_definition,
        "pressure_Pa": pressure,
        "temperature_K": temperature,
        "mean_energy_support_eV": list(transport_support),
        "support_min_mean_energy_eV": transport_support[0],
        "support_max_mean_energy_eV": transport_support[1],
        "table_paths": {name: str(path) for name, path in selected_paths.items()},
        "artifact_verification": {
            "status": "passed",
            "verified_files": len(digests),
            "sha256": dict(sorted(digests.items())),
        },
        "provenance": {
            "format_version": manifest["format_version"],
            "stage": manifest["stage"],
            "hashes": hashes,
            "physical_context": manifest["physical_context"],
            "mixture": mixture,
            "source_policy": source_policy,
            "quality_summary": quality_summary,
            "solver_qualification": solver_qualification,
            "target_qualification": target_qualification,
            "mc_qualification": mc_qualification,
            "source_composition": manifest.get("source_composition"),
            "evidence": manifest.get("evidence"),
        },
    }


def _read_manifest(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, json.JSONDecodeError) as exc:
        raise GecIcpBundleError(f"cannot read GEC-ICP bundle manifest: {path}") from exc
    if not isinstance(value, dict):
        raise GecIcpBundleError("GEC-ICP bundle manifest must be a JSON object")
    return value


def _validate_manifest_identity(
    manifest: dict[str, Any],
    *,
    expected_source: str,
    expected_transport_definition: str | None,
) -> tuple[str, str, dict[str, Any], dict[str, Any]]:
    if (
        manifest.get("format_version") != FORMAT_VERSION
        or manifest.get("stage") != "export-comsol"
        or manifest.get("status") != "ok"
        or manifest.get("missing_required_coefficients") != []
    ):
        raise GecIcpBundleError("bundle is not a successful canonical COMSOL export")
    source = manifest.get("source")
    if source != expected_source:
        raise GecIcpBundleError(
            f"GEC-ICP bundle source mismatch: expected {expected_source}, found {source}"
        )
    source_policy = manifest.get("source_policy")
    if not isinstance(source_policy, dict):
        raise GecIcpBundleError("bundle lacks source_policy provenance")
    policy_source = source_policy.get("source")
    rf_frequency = source_policy.get("rf_frequency_Hz")
    if (
        (policy_source is not None and policy_source != source)
        or source_policy.get("postprocess") != "none"
        or source_policy.get("field_type") != "dc"
        or (rf_frequency is not None and rf_frequency != "")
    ):
        raise GecIcpBundleError("bundle source policy is not an unmodified DC source")
    transport_definition = source_policy.get("transport_definition")
    if not isinstance(transport_definition, str) or not transport_definition.strip():
        raise GecIcpBundleError("bundle lacks a transport_definition")
    if expected_transport_definition is not None and (
        not isinstance(expected_transport_definition, str)
        or not expected_transport_definition.strip()
        or transport_definition != expected_transport_definition
    ):
        raise GecIcpBundleError(
            "GEC-ICP transport_definition mismatch: expected "
            f"{expected_transport_definition}, found {transport_definition}"
        )
    hashes = manifest.get("hashes")
    if not isinstance(hashes, dict) or any(
        _SHA256_PATTERN.fullmatch(str(hashes.get(name, ""))) is None
        for name in _REQUIRED_PROVENANCE_HASHES
    ):
        raise GecIcpBundleError(
            "bundle lacks canonical workflow, config, or cross-section SHA-256 provenance"
        )
    try:
        expected_cross_sections = canonical_cross_sections_combined_sha256()
    except GecIcpContractError as exc:
        raise GecIcpBundleError(str(exc)) from exc
    if hashes["cross_sections_sha256"] != expected_cross_sections:
        raise GecIcpBundleError(
            "bundle cross_sections_sha256 does not match the repository canonical "
            "examples/cross_sections/argon_application_library.csv combined hash"
        )
    return str(source), transport_definition, source_policy, hashes


def _validate_physical_context(
    manifest: dict[str, Any],
    *,
    expected_pressure_Pa: float,
    expected_temperature_K: float,
) -> tuple[float, float, dict[str, Any]]:
    context = manifest.get("physical_context")
    if (
        not isinstance(context, dict)
        or context.get("schema") != "swarm_physical_context.v1"
    ):
        raise GecIcpBundleError("bundle lacks canonical physical-context provenance")
    field = context.get("field")
    if not isinstance(field, dict) or field.get("type") != "dc":
        raise GecIcpBundleError("GEC-ICP requires a homogeneous DC Swarm source")
    temperature = _positive_number(
        context.get("gas_temperature_K"), "bundle gas temperature"
    )
    density = _positive_number(
        context.get("gas_number_density_m3"), "bundle gas number density"
    )
    pressure = density * _BOLTZMANN_J_K * temperature
    if not math.isclose(
        temperature, expected_temperature_K, rel_tol=1.0e-10, abs_tol=1.0e-10
    ):
        raise GecIcpBundleError("bundle gas temperature does not match GEC-ICP")
    if not math.isclose(
        pressure, expected_pressure_Pa, rel_tol=1.0e-10, abs_tol=1.0e-10
    ):
        raise GecIcpBundleError("bundle pressure does not match GEC-ICP")

    mixture = manifest.get("mixture")
    species = mixture.get("species") if isinstance(mixture, dict) else None
    if not isinstance(species, list) or len(species) != 1:
        raise GecIcpBundleError("GEC-ICP requires one pure-argon mixture entry")
    argon = species[0]
    if not isinstance(argon, dict) or argon.get("species") != "Ar":
        raise GecIcpBundleError("GEC-ICP bundle mixture must be Ar")
    fraction = _positive_number(argon.get("fraction"), "argon fraction")
    mass = _positive_number(argon.get("mass_amu"), "argon mass")
    if not math.isclose(
        fraction, 1.0, rel_tol=0.0, abs_tol=1.0e-12
    ) or not math.isclose(mass, _ARGON_MASS_AMU, rel_tol=0.0, abs_tol=1.0e-6):
        raise GecIcpBundleError(
            "GEC-ICP bundle must contain pure Ar with canonical mass"
        )
    return pressure, temperature, mixture


def _verify_table_artifacts(
    root: Path, manifest: dict[str, Any]
) -> tuple[dict[str, dict[str, Any]], dict[str, Path], dict[str, str]]:
    raw_tables = manifest.get("tables")
    if not isinstance(raw_tables, dict) or not raw_tables:
        raise GecIcpBundleError("bundle manifest.tables must be a nonempty mapping")
    tables: dict[str, dict[str, Any]] = {}
    paths: dict[str, Path] = {}
    digests: dict[str, str] = {}
    for name, metadata in raw_tables.items():
        path = _contained_artifact(root, name)
        role = metadata.get("artifact_role") if isinstance(metadata, dict) else None
        if (
            not isinstance(metadata, dict)
            or not isinstance(role, str)
            or not role.strip()
        ):
            raise GecIcpBundleError(f"bundle artifact lacks a role: {name}")
        expected = metadata.get("sha256")
        actual = _file_sha256(path)
        if _SHA256_PATTERN.fullmatch(str(expected or "")) is None or expected != actual:
            raise GecIcpBundleError(f"bundle artifact SHA-256 mismatch: {name}")
        tables[name] = metadata
        paths[name] = path
        digests[name] = actual
    return tables, paths, digests


def _verify_composite_evidence(
    root: Path,
    source: str,
    manifest: dict[str, Any],
) -> tuple[dict[str, Path], dict[str, str]]:
    composition = manifest.get("source_composition")
    evidence = manifest.get("evidence")
    if source != ANCHOR_FALLBACK_SOURCE:
        if composition is not None or evidence is not None:
            raise GecIcpBundleError(
                "non-composite bundle contains composite source evidence"
            )
        return {}, {}
    composition_evidence = (
        composition.get("evidence") if isinstance(composition, dict) else None
    )
    if (
        not isinstance(evidence, dict)
        or set(evidence) != set(_COMPOSITE_EVIDENCE_ROLES)
        or composition_evidence != evidence
    ):
        raise GecIcpBundleError(
            "composite bundle evidence inventory is incomplete or inconsistent"
        )
    paths: dict[str, Path] = {}
    digests: dict[str, str] = {}
    for name, expected_role in _COMPOSITE_EVIDENCE_ROLES.items():
        metadata = evidence.get(name)
        if (
            not isinstance(metadata, dict)
            or set(metadata) != {"role", "sha256"}
            or metadata.get("role") != expected_role
        ):
            raise GecIcpBundleError(f"composite evidence metadata is invalid: {name}")
        path = _contained_artifact(root, name)
        actual = _file_sha256(path)
        if metadata.get("sha256") != actual:
            raise GecIcpBundleError(f"composite evidence SHA-256 mismatch: {name}")
        paths[name] = path
        digests[f"evidence/{name}"] = actual
    return paths, digests


def _validate_active_function_eedf(
    source: str,
    manifest: dict[str, Any],
    tables: dict[str, dict[str, Any]],
    paths: dict[str, Path],
) -> tuple[Path, tuple[float, float]]:
    expected_name = COMSOL_FUNCTION_EEDF_TABLE
    reference = manifest.get("eedf_artifacts")
    reference = reference.get("comsol_input") if isinstance(reference, dict) else None
    canonical = [
        name
        for name, metadata in tables.items()
        if metadata.get("artifact_role") == "canonical_comsol_function_eedf_input"
        and metadata.get("canonical_comsol_input") is True
    ]
    if not isinstance(reference, dict) or canonical != [expected_name]:
        raise GecIcpBundleError(
            "bundle does not identify one canonical active Function-EEDF"
        )
    metadata = tables[expected_name]
    if (
        reference.get("file") != expected_name
        or reference.get("sha256") != metadata.get("sha256")
        or reference.get("role") != metadata.get("artifact_role")
    ):
        raise GecIcpBundleError("active Function-EEDF reference is inconsistent")
    path = paths[expected_name]
    try:
        contract = ComsolEedfImportContract.from_path(path)
        grid = read_comsol_function_eedf_grid(path)
    except FunctionEedfError as exc:
        raise GecIcpBundleError(f"invalid active Function-EEDF: {exc}") from exc
    if (
        contract.structure != "spreadsheet"
        or metadata.get("representation") != contract.representation
        or metadata.get("format") != contract.format
        or metadata.get("columns") != list(COMSOL_EEDF_COLUMNS)
        or metadata.get("argument_order") != list(COMSOL_EEDF_COLUMNS[:2])
        or metadata.get("grid_axis_order") != list(COMSOL_EEDF_COLUMNS[:2])
        or metadata.get("row_order") != contract.row_order
        or not isinstance(metadata.get("comsol_import"), dict)
        or any(
            metadata["comsol_import"].get(name) != value
            for name, value in contract.import_settings().items()
        )
    ):
        raise GecIcpBundleError(
            "active Function-EEDF metadata violates its import contract"
        )
    shape = [len(grid.mean_energies_eV), len(grid.electron_energies_eV)]
    support = (
        float(grid.mean_energies_eV[0]),
        float(grid.mean_energies_eV[-1]),
    )
    if (
        metadata.get("grid_shape") != shape
        or metadata.get("mean_energy_grid_points") != shape[0]
        or metadata.get("energy_grid_points") != shape[1]
        or not _numeric_pair_matches(metadata.get("mean_energy_range_eV"), support)
        or grid.normalization_error_max > 1.0e-8
        or grid.mean_energy_relative_error_max > 1.0e-8
        or grid.minimum_value < 0.0
        or not _reported_metric_matches(
            metadata,
            "projected_normalization_error_max",
            grid.normalization_error_max,
        )
        or not _reported_metric_matches(
            metadata,
            "projected_mean_energy_relative_error_max",
            grid.mean_energy_relative_error_max,
        )
        or not _reported_metric_matches(
            metadata,
            "projected_nonnegative_minimum",
            grid.minimum_value,
        )
    ):
        raise GecIcpBundleError("active Function-EEDF grid or moment audit failed")
    return path, support


def _validate_transport_table(
    manifest: dict[str, Any],
    tables: dict[str, dict[str, Any]],
    paths: dict[str, Path],
) -> tuple[Path, tuple[float, float]]:
    metadata = tables.get(_TRANSPORT_TABLE)
    if not isinstance(metadata, dict):
        raise GecIcpBundleError("bundle lacks mean-energy transport table")
    required_columns = {
        "mean_energy_eV",
        "E_over_N_Td",
        "reduced_mobility_m2_V_s_m3",
    }
    listed_columns = metadata.get("columns")
    if (
        metadata.get("artifact_role") != "canonical_comsol_coefficient_input"
        or metadata.get("argument") != "mean_energy_eV"
        or not isinstance(listed_columns, list)
        or not required_columns.issubset(listed_columns)
    ):
        raise GecIcpBundleError("mean-energy transport metadata is not canonical")
    rows = _read_csv(paths[_TRANSPORT_TABLE])
    means = [
        _positive_number(row.get("mean_energy_eV"), "transport mean energy")
        for row in rows
    ]
    fields = [_positive_number(row.get("E_over_N_Td"), "transport E/N") for row in rows]
    mobilities = [
        _positive_number(row.get("reduced_mobility_m2_V_s_m3"), "reduced mobility")
        for row in rows
    ]
    if len(means) < 2 or any(right <= left for left, right in zip(means, means[1:])):
        raise GecIcpBundleError(
            "mean-energy transport support is not strictly increasing"
        )
    if len(fields) != len(means) or len(mobilities) != len(means):
        raise GecIcpBundleError(
            "mean-energy transport columns have inconsistent lengths"
        )
    monotonicity = manifest.get("monotonicity")
    if (
        not isinstance(monotonicity, dict)
        or monotonicity.get("mean_energy_strictly_monotonic") is not True
    ):
        raise GecIcpBundleError("bundle does not certify monotonic mean-energy lookup")
    support = (means[0], means[-1])
    valid_ranges = manifest.get("valid_ranges")
    if not isinstance(valid_ranges, dict) or not _numeric_pair_matches(
        valid_ranges.get("mean_energy_eV"), support
    ):
        raise GecIcpBundleError(
            "manifest mean-energy range disagrees with transport table"
        )
    return paths[_TRANSPORT_TABLE], support


def _validate_quality(
    source: str,
    manifest: dict[str, Any],
    tables: dict[str, dict[str, Any]],
    paths: dict[str, Path],
) -> tuple[Path, dict[str, Any]]:
    metadata = tables.get(_QUALITY_TABLE)
    if (
        not isinstance(metadata, dict)
        or metadata.get("artifact_role") != "swarm_quality_audit"
    ):
        raise GecIcpBundleError("bundle lacks canonical quality evidence")
    rows = _read_csv(paths[_QUALITY_TABLE])
    gate = "active_closure_quality_passed" if source == MC_SOLVER else "passed"
    expected_columns = (
        COMPOSITE_QUALITY_COLUMNS
        if source == ANCHOR_FALLBACK_SOURCE
        else quality_table_schema(source).columns
    )
    if tuple(metadata.get("columns") or ()) != expected_columns:
        raise GecIcpBundleError(
            "bundle quality table does not use its canonical schema"
        )
    if source == ANCHOR_FALLBACK_SOURCE and (
        metadata.get("schema") != COMPOSITE_QUALITY_SCHEMA
        or metadata.get("evidence_kind") != "per_anchor_selected_source_qualification"
    ):
        raise GecIcpBundleError(
            "composite quality table lacks its anchor-selection schema"
        )
    if not rows or gate not in (metadata.get("columns") or []):
        raise GecIcpBundleError("bundle quality evidence lacks its acceptance gate")
    if source == MC_SOLVER and any(
        row.get("qualification_profile")
        != MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
        for row in rows
    ):
        raise GecIcpBundleError(
            "MC quality rows do not use restricted-LMEA qualification"
        )
    if any(not _accepted_flag(row.get(gate)) for row in rows):
        raise GecIcpBundleError("bundle contains unaccepted quality points")
    summary = manifest.get("quality_summary")
    if (
        not isinstance(summary, dict)
        or summary.get("passed") is not True
        or summary.get("failed_points") != 0
        or summary.get("total_points") != len(rows)
    ):
        raise GecIcpBundleError("bundle quality summary disagrees with quality.csv")
    return paths[_QUALITY_TABLE], summary


def _validate_solver_qualifications(
    source: str,
    root: Path,
    manifest: dict[str, Any],
    *,
    hashes: dict[str, Any],
) -> tuple[dict[str, Any] | None, dict[str, Any] | None]:
    entry = manifest.get("solver_qualification")
    target_entry = manifest.get("target_qualification")
    if source != PROPAGATOR_SOLVER:
        if entry is not None or target_entry is not None:
            raise GecIcpBundleError(
                "solver and target qualification are only valid for propagator"
            )
        return None, None
    if not isinstance(entry, dict):
        raise GecIcpBundleError("propagator bundle lacks core solver qualification")
    try:
        path = validate_copied_qualification(root, entry)
        validated = validate_propagator_core_qualification(path)
    except SolverQualificationError as exc:
        raise GecIcpBundleError(str(exc)) from exc
    expected = validated.manifest_entry(file=path.name)
    if entry != expected:
        raise GecIcpBundleError(
            "propagator qualification manifest entry is inconsistent"
        )
    if target_entry is None:
        raise GecIcpBundleError(
            "propagator bundle lacks required generic ICP target qualification"
        )
    if not isinstance(target_entry, dict):
        raise GecIcpBundleError("propagator target qualification entry is malformed")
    try:
        target_path = validate_copied_qualification(root, target_entry)
        if target_entry.get("schema") != PROPAGATOR_TARGET_QUALIFICATION_SCHEMA:
            raise GecIcpBundleError(
                "GEC-ICP requires the generic propagator target qualification schema"
            )
        target = validate_propagator_target_qualification(
            target_path,
            core_qualification=validated,
            provenance_hashes=hashes,
            requirement=_ICP_PROPAGATOR_TARGET_REQUIREMENT,
        )
    except SolverQualificationError as exc:
        raise GecIcpBundleError(str(exc)) from exc
    if target.schema != PROPAGATOR_TARGET_QUALIFICATION_SCHEMA:
        raise GecIcpBundleError(
            "GEC-ICP requires the generic propagator target qualification schema"
        )
    expected_target = target.manifest_entry(file=target_path.name)
    if target_entry != expected_target:
        raise GecIcpBundleError(
            "propagator target qualification manifest entry is inconsistent"
        )
    return expected, expected_target


def _validate_mc_qualification(
    source: str,
    manifest: dict[str, Any],
    source_policy: dict[str, Any],
    tables: dict[str, dict[str, Any]],
    paths: dict[str, Path],
    evidence_paths: dict[str, Path],
    *,
    expected_profile: str | None,
) -> dict[str, Any] | None:
    if source == ANCHOR_FALLBACK_SOURCE:
        return _validate_composite_qualification(
            manifest,
            source_policy,
            paths,
            evidence_paths,
            expected_profile=expected_profile,
        )
    if source != MC_SOLVER:
        if manifest.get("mc_qualification") is not None or expected_profile is not None:
            raise GecIcpBundleError("MC qualification is only valid for monte_carlo")
        return None
    profile = MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
    if expected_profile not in {None, profile}:
        raise GecIcpBundleError(
            f"unsupported GEC-ICP MC qualification profile: {expected_profile}"
        )
    if source_policy.get("qualification_profile") != profile or source_policy.get(
        "qualification_outputs"
    ) != ["function_eedf", "reduced_mobility", "direct_elastic_energy_loss"]:
        raise GecIcpBundleError(
            "MC bundle does not use the restricted-LMEA closure profile"
        )
    entry = tables.get(MC_QUALIFICATION_TABLE)
    qualification = manifest.get("mc_qualification")
    if (
        not isinstance(entry, dict)
        or entry.get("artifact_role") != "all_planned_mc_anchor_qualification"
        or not isinstance(qualification, dict)
        or qualification.get("file") != MC_QUALIFICATION_TABLE
    ):
        raise GecIcpBundleError("MC bundle lacks canonical qualification evidence")
    rows = _read_csv(paths[MC_QUALIFICATION_TABLE])
    if not rows or any(
        row.get("qualification_profile") != profile
        or not _accepted_flag(row.get("active_closure_quality_passed"))
        for row in rows
    ):
        raise GecIcpBundleError("MC restricted-LMEA qualification did not pass")
    if (
        qualification.get("all_planned_anchors") != len(rows)
        or qualification.get("qualified_anchors") != len(rows)
        or qualification.get("coefficient_tables_available") is not True
        or qualification.get("table_build_failure") is not None
    ):
        raise GecIcpBundleError("MC qualification summary is incomplete or failed")
    hashes = manifest["hashes"]
    raw_plan = hashes.get("mc_sampling_plan_json")
    plan_digest = hashes.get("mc_sampling_plan_sha256")
    eedf_estimator_schema = hashes.get("mc_eedf_estimator_schema_version")
    solver_source_sha256 = hashes.get("mc_solver_source_sha256")
    recorded_plan = manifest.get("mc_sampling_plan")
    try:
        decoded = json.loads(raw_plan) if isinstance(raw_plan, str) else None
    except json.JSONDecodeError as exc:
        raise GecIcpBundleError("MC sampling-plan provenance is invalid") from exc
    if (
        not isinstance(decoded, list)
        or _SHA256_PATTERN.fullmatch(str(plan_digest or "")) is None
        or sha256(raw_plan.encode("utf-8")).hexdigest() != plan_digest
        or not isinstance(recorded_plan, dict)
        or recorded_plan.get("status") != "validated_against_mc_cases"
        or recorded_plan.get("sha256") != plan_digest
        or recorded_plan.get("entries") != decoded
        or eedf_estimator_schema != _mc_evidence.MC_EEDF_ESTIMATOR_SCHEMA_VERSION
        or solver_source_sha256 != monte_carlo_source_sha256()
    ):
        raise GecIcpBundleError(
            "MC sampling-plan, EEDF estimator, or solver-source provenance is "
            "inconsistent"
        )
    return {
        **qualification,
        "qualification_profile": profile,
        "sampling_plan_sha256": plan_digest,
        "eedf_estimator_schema_version": str(eedf_estimator_schema),
        "solver_source_sha256": str(solver_source_sha256),
    }


def _validate_composite_qualification(
    manifest: dict[str, Any],
    source_policy: dict[str, Any],
    table_paths: dict[str, Path],
    evidence_paths: dict[str, Path],
    *,
    expected_profile: str | None,
) -> dict[str, Any]:
    profile = MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
    if expected_profile not in {None, profile}:
        raise GecIcpBundleError(
            f"unsupported GEC-ICP composite qualification profile: {expected_profile}"
        )
    if (
        source_policy.get("source") != ANCHOR_FALLBACK_SOURCE
        or source_policy.get("primary_solver") != MC_SOLVER
        or source_policy.get("fallback_solver") != "two_term"
        or source_policy.get("selection")
        != "qualified_mc_else_exact_qualified_two_term_anchor"
        or source_policy.get("qualification_profile") != profile
        or source_policy.get("qualification_outputs")
        != ["function_eedf", "reduced_mobility"]
        or source_policy.get("component_mixing_within_anchor") is not False
        or source_policy.get("coefficient_repair") is not False
    ):
        raise GecIcpBundleError(
            "composite bundle does not declare the restricted anchor-fallback closure"
        )

    composition = manifest.get("source_composition")
    if (
        not isinstance(composition, dict)
        or composition.get("schema") != ANCHOR_FALLBACK_SCHEMA
        or composition.get("primary_solver") != MC_SOLVER
        or composition.get("fallback_solver") != "two_term"
        or composition.get("scope") != "complete_anchor_replacement"
        or composition.get("selection_rule")
        != "qualified_monte_carlo_else_exact_qualified_two_term"
        or composition.get("component_mixing_within_anchor") is not False
        or composition.get("postprocess_repair") is not False
        or composition.get("decision_relationship", {}).get("status")
        != "superseded_for_composite"
    ):
        raise GecIcpBundleError("composite source_composition contract is invalid")

    hashes = manifest.get("hashes")
    if not isinstance(hashes, dict):
        raise GecIcpBundleError("composite bundle lacks provenance hashes")
    mc_manifest_hash = _file_sha256(evidence_paths["monte_carlo_source_manifest.json"])
    two_term_manifest_hash = _file_sha256(
        evidence_paths["two_term_source_manifest.json"]
    )
    plan_hash = _file_sha256(evidence_paths["anchor_fallback_plan.json"])
    if (
        hashes.get("monte_carlo_manifest_sha256") != mc_manifest_hash
        or hashes.get("two_term_manifest_sha256") != two_term_manifest_hash
        or hashes.get("anchor_fallback_plan_sha256") != plan_hash
        or composition.get("mc_database_sha256")
        != hashes.get("monte_carlo_database_sha256")
    ):
        raise GecIcpBundleError(
            "composite source manifests, plan, or MC database binding changed"
        )

    anchors, mc_fields, fallback_fields = _composite_anchor_selection(
        composition,
        mc_manifest_hash=mc_manifest_hash,
        two_term_manifest_hash=two_term_manifest_hash,
    )
    quality_rows = _read_csv(table_paths[_QUALITY_TABLE])
    _validate_composite_quality_rows(quality_rows, anchors)

    mc_rows = _read_csv(evidence_paths["monte_carlo_qualification.csv"])
    mc_by_field = _unique_rows_by_field(mc_rows, "Monte Carlo qualification")
    if set(mc_by_field) != set(anchors):
        raise GecIcpBundleError(
            "Monte Carlo qualification anchors differ from composite anchors"
        )
    for field, anchor in anchors.items():
        row = mc_by_field[field]
        try:
            reasons = _string_list_json(
                row.get("active_closure_failure_reasons_json"),
                "Monte Carlo active-closure failure reasons",
            )
        except GecIcpBundleError as exc:
            raise GecIcpBundleError(f"{exc} at {field:.17g} Td") from exc
        passed = _accepted_flag(row.get("active_closure_quality_passed"))
        expected_passed = anchor["effective_source"] == MC_SOLVER
        if (
            row.get("qualification_profile") != profile
            or passed != expected_passed
            or reasons != anchor["mc_failure_reasons"]
            or passed == bool(reasons)
        ):
            raise GecIcpBundleError(
                f"Monte Carlo qualification disagrees at {field:.17g} Td"
            )

    fallback_rows = _unique_rows_by_field(
        _read_csv(evidence_paths["two_term_quality.csv"]),
        "two-term qualification",
    )
    if any(
        field not in fallback_rows
        or not _accepted_flag(fallback_rows[field].get("passed"))
        for field in fallback_fields
    ):
        raise GecIcpBundleError(
            "two-term fallback is not qualified at every replacement anchor"
        )

    mc_manifest = _read_json_object(
        evidence_paths["monte_carlo_source_manifest.json"],
        "Monte Carlo source manifest",
    )
    two_term_manifest = _read_json_object(
        evidence_paths["two_term_source_manifest.json"],
        "two-term source manifest",
    )
    mc_solver_source_sha256 = mc_manifest.get("hashes", {}).get(
        "mc_solver_source_sha256"
    )
    if (
        mc_manifest.get("source") != MC_SOLVER
        or two_term_manifest.get("source") != "two_term"
        or mc_manifest.get("physical_context") != manifest.get("physical_context")
        or two_term_manifest.get("physical_context") != manifest.get("physical_context")
        or mc_manifest.get("mixture") != manifest.get("mixture")
        or two_term_manifest.get("mixture") != manifest.get("mixture")
        or mc_manifest.get("hashes", {}).get("cross_sections_sha256")
        != hashes.get("cross_sections_sha256")
        or two_term_manifest.get("hashes", {}).get("cross_sections_sha256")
        != hashes.get("cross_sections_sha256")
        or mc_manifest.get("hashes", {}).get("workflow_config_sha256")
        != hashes.get("workflow_config_sha256")
        or mc_manifest.get("hashes", {}).get("base_config_sha256")
        != hashes.get("base_config_sha256")
        or mc_manifest.get("mc_sampling_plan") != manifest.get("mc_sampling_plan")
        or mc_solver_source_sha256 != monte_carlo_source_sha256()
    ):
        raise GecIcpBundleError(
            "composite source snapshots disagree with bundle provenance"
        )

    decision_hash = _file_sha256(evidence_paths["mc_solver_decision.json"])
    decision = _read_json_object(
        evidence_paths["mc_solver_decision.json"], "MC policy decision"
    )
    plan = _read_json_object(
        evidence_paths["anchor_fallback_plan.json"], "anchor-fallback plan"
    )
    plan_composition = plan.get("source_composition")
    if (
        plan.get("format_version") != 1
        or plan.get("stage") != ANCHOR_FALLBACK_STAGE
        or plan.get("status") != "selected"
        or plan.get("source") != ANCHOR_FALLBACK_SOURCE
        or plan.get("quality_scope") != profile
        or plan.get("physical_context") != manifest.get("physical_context")
        or plan.get("mixture") != manifest.get("mixture")
        or plan.get("cross_sections_sha256") != hashes.get("cross_sections_sha256")
        or not isinstance(plan_composition, dict)
        or plan_composition.get("anchors") != composition.get("anchors")
        or plan_composition.get("decision_relationship")
        != composition.get("decision_relationship")
        or plan.get("decision", {}).get("sha256") != decision_hash
        or plan.get("inputs", {}).get(MC_SOLVER, {}).get("manifest_sha256")
        != mc_manifest_hash
        or plan.get("inputs", {}).get("two_term", {}).get("manifest_sha256")
        != two_term_manifest_hash
        or plan.get("inputs", {}).get(MC_SOLVER, {}).get("quality_sha256")
        != _file_sha256(evidence_paths["monte_carlo_qualification.csv"])
        or plan.get("inputs", {}).get("two_term", {}).get("quality_sha256")
        != _file_sha256(evidence_paths["two_term_quality.csv"])
    ):
        raise GecIcpBundleError(
            "anchor-fallback plan does not bind the exported composite"
        )
    failed_fields = sorted(fallback_fields)
    if (
        decision.get("status") != "selected"
        or decision.get("action") != "select_two_term"
        or decision.get("selected_solver") != "two_term"
        or decision.get("selection_scope") != "whole_comsol_closure"
        or decision.get("quality_scope") != profile
        or composition.get("decision_relationship")
        != {
            "terminal_scope": "whole_comsol_closure",
            "terminal_selection": "two_term",
            "status": "superseded_for_composite",
            "replacement_scope": "complete_anchor_replacement",
            "trigger": "explicit_per_anchor_fallback_request",
        }
        or not _field_lists_match(decision.get("failed_anchors_Td"), failed_fields)
        or decision.get("fallback", {}).get("qualified_for_required_anchors")
        is not True
        or decision.get("mc_run_identity", {}).get("mc_solver_source_sha256")
        != mc_solver_source_sha256
    ):
        raise GecIcpBundleError(
            "terminal MC policy decision does not authorize the selected fallback"
        )

    qualification = manifest.get("mc_qualification")
    if (
        not isinstance(qualification, dict)
        or qualification.get("profile") != profile
        or qualification.get("all_planned_anchors") != len(anchors)
        or qualification.get("qualified_monte_carlo_anchors") != len(mc_fields)
        or qualification.get("two_term_fallback_anchors") != len(fallback_fields)
        or not _field_lists_match(
            qualification.get("failed_monte_carlo_anchors_Td"), failed_fields
        )
        or qualification.get("coefficient_tables_available") is not True
        or qualification.get("fallback_applied") is not True
    ):
        raise GecIcpBundleError(
            "composite MC qualification summary is incomplete or inconsistent"
        )
    return {
        **qualification,
        "qualification_profile": profile,
        "primary_solver": MC_SOLVER,
        "fallback_solver": "two_term",
        "evidence_verified": True,
        "solver_source_sha256": str(mc_solver_source_sha256),
        "anchor_sources": {
            f"{field:.17g}": anchors[field]["effective_source"]
            for field in sorted(anchors)
        },
    }


def _composite_anchor_selection(
    composition: dict[str, Any],
    *,
    mc_manifest_hash: str,
    two_term_manifest_hash: str,
) -> tuple[dict[float, dict[str, Any]], set[float], set[float]]:
    raw_anchors = composition.get("anchors")
    if not isinstance(raw_anchors, list) or not raw_anchors:
        raise GecIcpBundleError("composite source has no anchor selections")
    anchors: dict[float, dict[str, Any]] = {}
    order: list[float] = []
    mc_fields: set[float] = set()
    fallback_fields: set[float] = set()
    for raw in raw_anchors:
        if not isinstance(raw, dict):
            raise GecIcpBundleError("composite anchor selection is malformed")
        field = _positive_number(raw.get("E_over_N_Td"), "composite E/N anchor")
        effective_source = raw.get("effective_source")
        try:
            reasons = _string_list(raw.get("mc_failure_reasons"))
        except GecIcpBundleError as exc:
            raise GecIcpBundleError(f"{exc} at {field:.17g} Td") from exc
        expected_hash = (
            mc_manifest_hash
            if effective_source == MC_SOLVER
            else two_term_manifest_hash
        )
        expected_reason = (
            "mc_anchor_qualified"
            if effective_source == MC_SOLVER
            else "mc_anchor_unqualified"
        )
        if (
            effective_source not in {MC_SOLVER, "two_term"}
            or field in anchors
            or raw.get("selection_reason") != expected_reason
            or raw.get("source_manifest_sha256") != expected_hash
            or (effective_source == MC_SOLVER and reasons)
            or (effective_source == "two_term" and not reasons)
        ):
            raise GecIcpBundleError(
                f"composite source selection is inconsistent at {field:.17g} Td"
            )
        entry = dict(raw)
        entry["mc_failure_reasons"] = reasons
        anchors[field] = entry
        order.append(field)
        (mc_fields if effective_source == MC_SOLVER else fallback_fields).add(field)
    if order != sorted(order) or not mc_fields or not fallback_fields:
        raise GecIcpBundleError(
            "composite anchors must be increasing and include MC plus fallback"
        )
    return anchors, mc_fields, fallback_fields


def _validate_composite_quality_rows(
    rows: list[dict[str, str]],
    anchors: dict[float, dict[str, Any]],
) -> None:
    by_field = _unique_rows_by_field(rows, "composite quality")
    if set(by_field) != set(anchors):
        raise GecIcpBundleError(
            "composite quality anchors differ from source_composition"
        )
    for field, anchor in anchors.items():
        row = by_field[field]
        reasons = _string_list_json(
            row.get("mc_failure_reasons_json"),
            "composite MC failure reasons",
        )
        expected_source = anchor["effective_source"]
        expected_gate = (
            "active_closure_quality_passed"
            if expected_source == MC_SOLVER
            else "passed"
        )
        if (
            not _accepted_flag(row.get("passed"))
            or row.get("effective_source") != expected_source
            or row.get("source_quality_gate") != expected_gate
            or row.get("selection_reason") != anchor["selection_reason"]
            or row.get("source_manifest_sha256") != anchor["source_manifest_sha256"]
            or reasons != anchor["mc_failure_reasons"]
        ):
            raise GecIcpBundleError(
                f"composite quality evidence disagrees at {field:.17g} Td"
            )


def _unique_rows_by_field(
    rows: list[dict[str, str]], label: str
) -> dict[float, dict[str, str]]:
    result: dict[float, dict[str, str]] = {}
    for row in rows:
        field = _positive_number(row.get("E_over_N_Td"), f"{label} E/N")
        if field in result:
            raise GecIcpBundleError(f"{label} repeats {field:.17g} Td")
        result[field] = row
    if not result:
        raise GecIcpBundleError(f"{label} is empty")
    return result


def _string_list_json(value: object, label: str) -> list[str]:
    try:
        decoded = json.loads(value) if isinstance(value, str) else None
    except json.JSONDecodeError as exc:
        raise GecIcpBundleError(f"{label} is invalid") from exc
    try:
        return _string_list(decoded)
    except GecIcpBundleError as exc:
        raise GecIcpBundleError(f"{label} is invalid") from exc


def _string_list(value: object) -> list[str]:
    if (
        not isinstance(value, list)
        or any(not isinstance(item, str) or not item for item in value)
        or len(value) != len(set(value))
    ):
        raise GecIcpBundleError("failure-reason list is invalid")
    return list(value)


def _field_lists_match(value: object, expected: list[float]) -> bool:
    if not isinstance(value, list) or len(value) != len(expected):
        return False
    try:
        actual = [float(item) for item in value]
    except (TypeError, ValueError):
        return False
    return all(
        math.isfinite(item)
        and math.isclose(item, target, rel_tol=1.0e-12, abs_tol=1.0e-12)
        for item, target in zip(actual, expected, strict=True)
    )


def _read_json_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8-sig"))
    except (OSError, json.JSONDecodeError) as exc:
        raise GecIcpBundleError(f"cannot read {label}") from exc
    if not isinstance(value, dict):
        raise GecIcpBundleError(f"{label} must be a JSON object")
    return value


def _contained_artifact(root: Path, name: object) -> Path:
    relative = PurePosixPath(name) if isinstance(name, str) else None
    if (
        relative is None
        or not name
        or "\\" in name
        or relative.is_absolute()
        or relative.as_posix() != name
        or any(part in {"", ".", ".."} for part in relative.parts)
    ):
        raise GecIcpBundleError("bundle manifest contains an unsafe artifact path")
    path = (root / Path(*relative.parts)).resolve()
    if not path.is_relative_to(root) or not path.is_file():
        raise GecIcpBundleError(f"bundle artifact does not exist: {name}")
    return path


def _read_csv(path: Path) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            if not reader.fieldnames:
                raise GecIcpBundleError(f"bundle table has no header: {path.name}")
            return list(reader)
    except (OSError, csv.Error) as exc:
        raise GecIcpBundleError(f"cannot read bundle table: {path.name}") from exc


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _positive_number(value: object, label: str) -> float:
    if isinstance(value, bool):
        raise GecIcpBundleError(f"{label} must be positive and finite")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise GecIcpBundleError(f"{label} must be positive and finite") from exc
    if not math.isfinite(result) or result <= 0.0:
        raise GecIcpBundleError(f"{label} must be positive and finite")
    return result


def _numeric_pair_matches(value: object, expected: tuple[float, float]) -> bool:
    if not isinstance(value, list) or len(value) != 2:
        return False
    try:
        actual = (float(value[0]), float(value[1]))
    except (TypeError, ValueError):
        return False
    return all(
        math.isfinite(item)
        and math.isclose(item, reference, rel_tol=1.0e-10, abs_tol=1.0e-12)
        for item, reference in zip(actual, expected, strict=True)
    )


def _accepted_flag(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true"}


def _reported_metric_matches(
    metadata: dict[str, Any], name: str, actual: float
) -> bool:
    try:
        reported = float(metadata[name])
    except (KeyError, TypeError, ValueError):
        return False
    return math.isfinite(reported) and math.isclose(
        reported, actual, rel_tol=1.0e-9, abs_tol=1.0e-13
    )
