"""Orchestrate provenance-bound COMSOL bundle export."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from swarm_workflow._io import write_json as _write_json
from swarm_workflow.comsol.input.bundle_selection import (
    SELECTION_FILE,
    SOURCE_MANIFEST_FILE,
    validate_bundle_selection,
)
from swarm_workflow.selection import (
    ANCHOR_FALLBACK_SOURCE,
    ClosureSelectionError,
    validate_selected_tables,
)
from swarm_workflow.quality.solver import (
    PROPAGATOR_CORE_QUALIFICATION_FILE,
    PROPAGATOR_TARGET_QUALIFICATION_FILE,
    SolverQualificationError,
    validate_copied_qualification,
    validate_propagator_core_qualification,
    validate_propagator_target_qualification,
)
from swarm_workflow.quality.propagator_source import (
    PROPAGATOR_SOLVER_SOURCE_METADATA_KEY,
)
from swarm_workflow.tables.contracts import (
    COLLISION_RATE_KERNEL_TABLE,
    ELASTIC_ENERGY_LOSS_TABLE,
    RATE_EVIDENCE_TABLE,
)

from .contracts import (
    BASE_TABLES,
    FORMAT_VERSION,
    OPTIONAL_TABLES,
    ComsolExportError,
    ComsolExportSummary,
)
from .function_tables import _write_eedf_mean_energy_function_table
from .manifest import (
    _bundle_manifest,
    _canonical_eedf_source,
    _is_root_manifest,
    _read_manifest,
    _sha256_file,
    _table_metadata,
    _validate_source_artifact,
    _validated_quality_policy,
    _validated_quality_thresholds,
)
from .writers import (
    _clear_known_outputs,
    _csv_header,
    _missing_required_coefficients,
    _write_comsol_function_tables,
)


def export_comsol_bundle(
    table_directory: str | Path,
    output_directory: str | Path,
    *,
    selection_path: str | Path | None = None,
) -> ComsolExportSummary:
    table_dir = Path(table_directory).resolve()
    output_dir = Path(output_directory).resolve()
    if output_dir == table_dir or table_dir.is_relative_to(output_dir):
        raise ComsolExportError("export output must be separate from the source tables")
    if selection_path is not None and Path(selection_path).resolve().is_relative_to(
        output_dir
    ):
        raise ComsolExportError(
            "selection evidence must be kept outside the export output"
        )
    manifest = _read_manifest(table_dir)
    if _is_root_manifest(manifest):
        if selection_path is not None:
            raise ComsolExportError(
                "a whole-closure selection exports one mixture directory at a time"
            )
        root_quality = _validated_quality_thresholds(manifest)
        root_quality_policy = _validated_quality_policy(manifest)
        entries = manifest["mixtures"]
        output_dir.mkdir(parents=True, exist_ok=True)
        bundle_entries = []
        for entry in entries:
            source_manifest = table_dir / str(entry["path"])
            source_dir = source_manifest.parent
            child_quality = _validated_quality_thresholds(_read_manifest(source_dir))
            if child_quality != root_quality:
                raise ComsolExportError(
                    "mixture table quality thresholds differ from root manifest"
                )
            if (
                _validated_quality_policy(_read_manifest(source_dir))
                != root_quality_policy
            ):
                raise ComsolExportError(
                    "mixture table quality policy differs from root manifest"
                )
            target_dir = output_dir / f"mixture_{int(entry['mixture_id']):04d}"
            _export_single(source_dir, target_dir)
            bundle_entries.append(
                {
                    "mixture_id": int(entry["mixture_id"]),
                    "path": f"mixture_{int(entry['mixture_id']):04d}/manifest.json",
                }
            )
        _write_json(
            output_dir / "manifest.json",
            {
                "format_version": FORMAT_VERSION,
                "stage": "export-comsol",
                "source": manifest.get("source"),
                "hashes": manifest.get("hashes", {}),
                "mc_sampling_plan": manifest.get("mc_sampling_plan"),
                "quality_thresholds": root_quality,
                "quality_policy": root_quality_policy,
                "bundles": bundle_entries,
            },
        )
        return ComsolExportSummary(table_dir, output_dir, len(bundle_entries))

    output_dir.mkdir(parents=True, exist_ok=True)
    _export_single(table_dir, output_dir, selection_path=selection_path)
    return ComsolExportSummary(table_dir, output_dir, 1)


def _export_single(
    table_dir: Path,
    output_dir: Path,
    *,
    selection_path: str | Path | None = None,
) -> None:
    source_manifest = _read_manifest(table_dir)
    if source_manifest.get("status") == "unqualified":
        raise ComsolExportError(
            "MC qualification evidence has no exportable coefficient closure"
        )
    if selection_path is not None:
        try:
            validate_selected_tables(selection_path, table_dir)
        except ClosureSelectionError as exc:
            raise ComsolExportError(str(exc)) from exc
    _validated_quality_thresholds(source_manifest)
    _validated_quality_policy(source_manifest)
    output_dir.mkdir(parents=True, exist_ok=True)
    _clear_known_outputs(output_dir)
    listed_tables = source_manifest.get("tables")
    if not isinstance(listed_tables, dict):
        raise ComsolExportError("source table manifest.tables must be a mapping")
    required_tables = list(BASE_TABLES)
    if source_manifest.get("source") == "monte_carlo":
        required_tables.append(RATE_EVIDENCE_TABLE)
    missing_manifest_tables = [
        name for name in required_tables if name not in listed_tables
    ]
    if missing_manifest_tables:
        raise ComsolExportError(
            "source table manifest omits required tables: "
            + ", ".join(missing_manifest_tables)
        )
    qualification_entry = source_manifest.get("solver_qualification")
    propagator_core_qualification = None
    if source_manifest.get("source") == "propagator":
        if not isinstance(qualification_entry, dict):
            raise ComsolExportError(
                "propagator table manifest lacks bound solver qualification"
            )
        try:
            qualification_source = validate_copied_qualification(
                table_dir, qualification_entry
            )
            propagator_core_qualification = validate_propagator_core_qualification(
                qualification_source
            )
        except SolverQualificationError as exc:
            raise ComsolExportError(str(exc)) from exc
        expected_qualification_entry = propagator_core_qualification.manifest_entry(
            file=PROPAGATOR_CORE_QUALIFICATION_FILE
        )
        if qualification_entry != expected_qualification_entry:
            raise ComsolExportError(
                "propagator solver qualification manifest entry is noncanonical"
            )
        provenance_hashes = source_manifest.get("hashes")
        if (
            not isinstance(provenance_hashes, dict)
            or provenance_hashes.get(PROPAGATOR_SOLVER_SOURCE_METADATA_KEY)
            != propagator_core_qualification.implementation_fingerprint_sha256
        ):
            raise ComsolExportError(
                "propagator source provenance differs from the qualified implementation"
            )
        qualification_target = output_dir / PROPAGATOR_CORE_QUALIFICATION_FILE
        qualification_target.write_bytes(qualification_source.read_bytes())
        try:
            validate_copied_qualification(output_dir, expected_qualification_entry)
        except SolverQualificationError as exc:
            raise ComsolExportError(str(exc)) from exc
    elif qualification_entry is not None:
        raise ComsolExportError(
            "solver qualification is only valid for a propagator source"
        )
    target_qualification_entry = source_manifest.get("target_qualification")
    if target_qualification_entry is not None:
        if source_manifest.get("source") != "propagator" or not isinstance(
            target_qualification_entry, dict
        ):
            raise ComsolExportError(
                "target qualification is only valid for a propagator source"
            )
        try:
            target_qualification_source = validate_copied_qualification(
                table_dir, target_qualification_entry
            )
            if propagator_core_qualification is None:
                raise SolverQualificationError(
                    "Propagator target qualification lacks validated core evidence"
                )
            target_qualification = validate_propagator_target_qualification(
                target_qualification_source,
                core_qualification=propagator_core_qualification,
                provenance_hashes=source_manifest.get("hashes"),
                allow_self_described_generic=True,
            )
        except SolverQualificationError as exc:
            raise ComsolExportError(str(exc)) from exc
        expected_target_entry = target_qualification.manifest_entry(
            file=PROPAGATOR_TARGET_QUALIFICATION_FILE
        )
        if target_qualification_entry != expected_target_entry:
            raise ComsolExportError(
                "propagator target qualification manifest entry is noncanonical"
            )
        target_qualification_target = output_dir / PROPAGATOR_TARGET_QUALIFICATION_FILE
        target_qualification_target.write_bytes(
            target_qualification_source.read_bytes()
        )
        try:
            validate_copied_qualification(output_dir, expected_target_entry)
        except SolverQualificationError as exc:
            raise ComsolExportError(str(exc)) from exc
    _copy_composite_evidence(table_dir, output_dir, source_manifest)
    missing = _missing_required_coefficients(table_dir)
    if missing:
        failed = _bundle_manifest(
            source_manifest,
            status="failed",
            tables={},
            missing_required_coefficients=missing,
        )
        _write_json(output_dir / "manifest.json", failed)
        raise ComsolExportError(
            "COMSOL export missing required coefficients: " + ", ".join(missing)
        )

    tables: dict[str, dict[str, Any]] = {}
    for name in (*BASE_TABLES, *OPTIONAL_TABLES):
        if name not in listed_tables:
            continue
        source = table_dir / name
        if not source.is_file():
            raise ComsolExportError(
                f"manifest-listed source table does not exist: {name}"
            )
        if name == "eedf.csv":
            _validate_source_artifact(
                source,
                listed_tables[name],
                expected_roles={
                    "raw_swarm_eedf_audit_and_canonicalization_source",
                },
            )
        elif name == RATE_EVIDENCE_TABLE:
            _validate_source_artifact(
                source,
                listed_tables[name],
                expected_roles={"raw_swarm_rate_evidence"},
            )
        elif name == ELASTIC_ENERGY_LOSS_TABLE:
            _validate_source_artifact(
                source,
                listed_tables[name],
                expected_roles={"canonical_comsol_elastic_energy_loss_input"},
            )
        elif name == COLLISION_RATE_KERNEL_TABLE:
            _validate_source_artifact(
                source,
                listed_tables[name],
                expected_roles={"canonical_eedf_projection_rate_kernels"},
            )
            expected_cross_sections = source_manifest.get("hashes", {}).get(
                "cross_sections_sha256"
            )
            if (
                listed_tables[name].get("cross_sections_sha256")
                != expected_cross_sections
            ):
                raise ComsolExportError(
                    "collision-rate kernels are not bound to source cross sections"
                )
        columns = _csv_header(source)
        target = output_dir / name
        target.write_bytes(source.read_bytes())
        metadata = _table_metadata(columns, source_manifest, name)
        metadata["sha256"] = _sha256_file(target)
        if name == ELASTIC_ENERGY_LOSS_TABLE:
            physics_contract = listed_tables[name].get("physics_contract")
            if not isinstance(physics_contract, dict):
                raise ComsolExportError(
                    "elastic energy-loss table lacks its physics contract"
                )
            metadata["physics_contract"] = physics_contract
        if name == COLLISION_RATE_KERNEL_TABLE:
            metadata["cross_sections_sha256"] = listed_tables[name][
                "cross_sections_sha256"
            ]
        if name == "eedf.csv":
            metadata.update(
                {
                    "artifact_role": (
                        "raw_swarm_eedf_audit_and_canonicalization_source"
                    ),
                    "sha256": _sha256_file(target),
                }
            )
        tables[name] = metadata
    tables.update(_write_comsol_function_tables(output_dir))
    function_eedf_source = _canonical_eedf_source(
        table_dir,
        output_dir,
        source_manifest,
        listed_tables,
    )
    tables.update(
        _write_eedf_mean_energy_function_table(
            output_dir,
            source=function_eedf_source,
            source_kind=str(source_manifest.get("source", "")),
            source_composition=source_manifest.get("source_composition"),
            rate_kernel_path=(
                output_dir / COLLISION_RATE_KERNEL_TABLE
                if COLLISION_RATE_KERNEL_TABLE in tables
                else None
            ),
        )
    )
    bundle_manifest = _bundle_manifest(
        source_manifest,
        status="ok",
        tables=tables,
        missing_required_coefficients=[],
    )
    if selection_path is not None:
        try:
            selection = validate_selected_tables(selection_path, table_dir)
        except ClosureSelectionError as exc:
            raise ComsolExportError(str(exc)) from exc
        (output_dir / SELECTION_FILE).write_bytes(Path(selection_path).read_bytes())
        (output_dir / SOURCE_MANIFEST_FILE).write_bytes(
            (table_dir / "manifest.json").read_bytes()
        )
        bundle_manifest["solver_selection"] = {
            "file": SELECTION_FILE,
            "sha256": _sha256_file(output_dir / SELECTION_FILE),
            "source_manifest_sha256": _sha256_file(output_dir / SOURCE_MANIFEST_FILE),
            "selected_solver": selection["selected_solver"],
        }
    _write_json(output_dir / "manifest.json", bundle_manifest)
    if selection_path is not None:
        try:
            validate_bundle_selection(output_dir, required=True)
        except ClosureSelectionError as exc:
            raise ComsolExportError(str(exc)) from exc


def _copy_composite_evidence(
    table_dir: Path,
    output_dir: Path,
    source_manifest: dict[str, Any],
) -> dict[str, dict[str, str]]:
    """Copy immutable anchor-selection evidence without treating it as a table."""

    if source_manifest.get("source") != ANCHOR_FALLBACK_SOURCE:
        return {}
    evidence = source_manifest.get("evidence")
    composition = source_manifest.get("source_composition")
    composition_evidence = (
        composition.get("evidence") if isinstance(composition, dict) else None
    )
    if (
        not isinstance(evidence, dict)
        or not evidence
        or composition_evidence != evidence
    ):
        raise ComsolExportError(
            "composite source lacks one consistent evidence inventory"
        )
    copied: dict[str, dict[str, str]] = {}
    for name, metadata in sorted(evidence.items()):
        if (
            not isinstance(name, str)
            or not name
            or Path(name).name != name
            or not isinstance(metadata, dict)
        ):
            raise ComsolExportError("composite evidence path is unsafe")
        role = metadata.get("role")
        expected_hash = metadata.get("sha256")
        if (
            not isinstance(role, str)
            or not role.strip()
            or not isinstance(expected_hash, str)
            or len(expected_hash) != 64
            or any(character not in "0123456789abcdef" for character in expected_hash)
        ):
            raise ComsolExportError(f"composite evidence metadata is invalid: {name}")
        source = table_dir / name
        if not source.is_file() or _sha256_file(source) != expected_hash:
            raise ComsolExportError(f"composite evidence hash mismatch: {name}")
        target = output_dir / name
        if target.is_symlink():
            raise ComsolExportError(
                f"composite evidence output may not be a symlink: {name}"
            )
        target.write_bytes(source.read_bytes())
        if _sha256_file(target) != expected_hash:
            raise ComsolExportError(f"copied composite evidence hash mismatch: {name}")
        copied[name] = {"role": role, "sha256": expected_hash}
    return copied
