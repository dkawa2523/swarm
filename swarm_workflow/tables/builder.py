"""Build COMSOL-ready workflow table directories from aggregate data."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import shutil
import sqlite3
import tempfile
from typing import Any

from .._io import write_csv as _write_csv_rows
from .._io import write_json as _write_json
from . import contracts as _contracts
from . import collision_kernels as _collision_kernels
from . import energy_loss as _energy_loss
from . import monte_carlo as _mc_tables
from . import repository as _repository
from ..campaign.aggregate import aggregate_workflow_results
from ..campaign.quality import (
    resolve_evaluation_quality_thresholds,
    resolve_quality_thresholds,
    source_quality_thresholds_json,
    validate_aggregate_quality_thresholds,
)
from ..campaign.store import (
    WorkflowSchemaError,
    provenance_hash_manifest,
)
from ..selection import MC_QUALIFICATION_TABLE, context_from_metadata
from ..quality.solver import (
    PROPAGATOR_CORE_QUALIFICATION_FILE,
    PROPAGATOR_TARGET_QUALIFICATION_FILE,
    PropagatorTargetQualification,
    SolverQualification,
    SolverQualificationError,
    validate_propagator_core_qualification,
    validate_propagator_target_qualification,
)
from ..quality.propagator_source import PROPAGATOR_SOLVER_SOURCE_METADATA_KEY
from ..quality.table import COMMON_QUALITY_COLUMNS, quality_table_schema
from ..quality.policy import (
    QualityThresholds,
    quality_policy_provenance,
    quality_thresholds_json,
    quality_thresholds_payload,
)
from .math import finite_range, strictly_monotonic
from ..campaign.repository import read_metadata


def build_tables(
    database_path: str | Path,
    output_directory: str | Path,
    *,
    source: str,
    thresholds: QualityThresholds | None = None,
    mc_qualification_profile: str = _contracts.MC_QUALIFICATION_FULL_TRANSPORT,
    allow_unqualified_mc: bool = False,
    solver_qualification_path: str | Path | None = None,
    target_qualification_path: str | Path | None = None,
) -> _contracts.TableBuildSummary:
    if source not in _contracts.SOURCE_CHOICES:
        raise _contracts.TableBuildError(f"unsupported table source: {source}")
    if allow_unqualified_mc and source != _contracts.MC_SOLVER:
        raise _contracts.TableBuildError(
            "allow_unqualified_mc requires source=monte_carlo"
        )
    if mc_qualification_profile not in _contracts.MC_QUALIFICATION_PROFILES:
        raise _contracts.TableBuildError(
            f"unsupported Monte Carlo qualification profile: {mc_qualification_profile}"
        )
    if (
        source != _contracts.MC_SOLVER
        and mc_qualification_profile != _contracts.MC_QUALIFICATION_FULL_TRANSPORT
    ):
        raise _contracts.TableBuildError(
            "Monte Carlo qualification profiles are only valid for a "
            "monte_carlo table source"
        )
    db_path = Path(database_path)
    solver_qualification: SolverQualification | None = None
    target_qualification: PropagatorTargetQualification | None = None
    if source == _contracts.PROPAGATOR_SOLVER:
        if solver_qualification_path is None:
            raise _contracts.TableBuildError(
                "propagator tables require --solver-qualification so the "
                "numerical release evidence is bound into the bundle"
            )
        try:
            solver_qualification = validate_propagator_core_qualification(
                solver_qualification_path
            )
        except SolverQualificationError as exc:
            raise _contracts.TableBuildError(str(exc)) from exc
        if target_qualification_path is not None:
            try:
                target_qualification = validate_propagator_target_qualification(
                    target_qualification_path,
                    core_qualification=solver_qualification,
                    database_path=db_path,
                    allow_self_described_generic=True,
                )
            except SolverQualificationError as exc:
                raise _contracts.TableBuildError(str(exc)) from exc
    elif solver_qualification_path is not None or target_qualification_path is not None:
        raise _contracts.TableBuildError(
            "solver and target qualification paths are valid only for source=propagator"
        )
    output = Path(output_directory)
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(
            prefix=f".{output.name}.staging-",
            dir=output.parent,
        )
    )

    try:
        connection = sqlite3.connect(db_path)
        try:
            connection.row_factory = sqlite3.Row
            metadata = read_metadata(connection)
            if solver_qualification is not None:
                expected_source = (
                    solver_qualification.implementation_fingerprint_sha256
                )
                actual_source = metadata.get(
                    PROPAGATOR_SOLVER_SOURCE_METADATA_KEY
                )
                if actual_source != expected_source:
                    raise _contracts.TableBuildError(
                        "propagator database solver source does not match "
                        "the qualified implementation"
                    )
            mc_sampling_plan = (
                _repository._validate_mc_sampling_plan_against_cases(
                    connection, metadata
                )
                if source == _contracts.MC_SOLVER
                else None
            )
            source_thresholds = resolve_quality_thresholds(connection)
            source_policy_json = source_quality_thresholds_json(connection)
            if not _repository._aggregate_rows_are_current(connection):
                if thresholds is not None:
                    resolve_quality_thresholds(connection, thresholds)
                aggregate_workflow_results(connection, source_thresholds)
            evaluation_thresholds = resolve_evaluation_quality_thresholds(connection)
            if thresholds is not None and quality_thresholds_json(
                thresholds
            ) != quality_thresholds_json(evaluation_thresholds):
                raise WorkflowSchemaError(
                    "requested table quality thresholds differ from the "
                    "saved aggregate evaluation policy"
                )
            validate_aggregate_quality_thresholds(connection, evaluation_thresholds)
            # Aggregation is the last write performed by this command.  Bind every
            # emitted table manifest to these exact database bytes so a later
            # comparison cannot pair the qualification CSV with another campaign.
            connection.commit()
            source_database_sha256 = _sha256_file(db_path.resolve())
            mixture_ids = _repository._mixture_ids(connection, source)
            if not mixture_ids:
                raise _contracts.TableBuildError(
                    f"database has no usable {source} data"
                )
            entries: list[dict[str, Any]] = []
            for mixture_id in mixture_ids:
                mixture_dir = staging / f"mixture_{mixture_id:04d}"
                evidence = (
                    _mc_tables._assess_monte_carlo(
                        connection,
                        mixture_id,
                        quality_thresholds=evaluation_thresholds,
                        qualification_profile=mc_qualification_profile,
                    )
                    if source == _contracts.MC_SOLVER
                    else None
                )
                qualification_error: str | None = None
                try:
                    dataset = _build_dataset(
                        connection,
                        mixture_id,
                        source,
                        mc_sampling_plan=mc_sampling_plan,
                        quality_thresholds=evaluation_thresholds,
                        mc_qualification_profile=mc_qualification_profile,
                        mc_evidence=evidence,
                    )
                    _write_dataset_tables(
                        dataset,
                        mixture_dir,
                        metadata=metadata,
                        source_quality_thresholds=source_thresholds,
                        source_quality_json=source_policy_json,
                        quality_thresholds=evaluation_thresholds,
                        mc_sampling_plan=mc_sampling_plan,
                        solver_qualification=solver_qualification,
                        target_qualification=target_qualification,
                        source_database_sha256=source_database_sha256,
                    )
                except _contracts.MonteCarloQualificationError as exc:
                    if not allow_unqualified_mc:
                        raise
                    qualification_error = str(exc)
                if evidence is not None:
                    _write_mc_qualification_evidence(
                        mixture_dir,
                        evidence,
                        connection=connection,
                        mixture_id=mixture_id,
                        metadata=metadata,
                        sampling_plan=mc_sampling_plan,
                        qualification_profile=mc_qualification_profile,
                        qualification_error=qualification_error,
                        source_database_sha256=source_database_sha256,
                    )
                entries.append(
                    {
                        "mixture_id": mixture_id,
                        "path": f"mixture_{mixture_id:04d}/manifest.json",
                        "status": "unqualified" if qualification_error else "ok",
                    }
                )
        finally:
            connection.close()

        _write_json(
            staging / "manifest.json",
            {
                "format_version": _contracts.FORMAT_VERSION,
                "stage": "build-tables",
                "source": source,
                "mc_qualification_profile": (
                    mc_qualification_profile if source == _contracts.MC_SOLVER else None
                ),
                "database_path": str(db_path),
                "hashes": {
                    **provenance_hash_manifest(metadata, source=source),
                    "source_database_sha256": source_database_sha256,
                },
                "solver_qualification": (
                    solver_qualification.manifest_entry(
                        file=PROPAGATOR_CORE_QUALIFICATION_FILE
                    )
                    if solver_qualification is not None
                    else None
                ),
                "target_qualification": (
                    target_qualification.manifest_entry(
                        file=PROPAGATOR_TARGET_QUALIFICATION_FILE
                    )
                    if target_qualification is not None
                    else None
                ),
                "mc_sampling_plan": mc_sampling_plan,
                "quality_thresholds": quality_thresholds_payload(evaluation_thresholds),
                "quality_policy": quality_policy_provenance(
                    source_thresholds,
                    evaluation_thresholds,
                    source_json=source_policy_json,
                ),
                "mixtures": entries,
            },
        )
        _replace_directory(staging, output)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return _contracts.TableBuildSummary(
        db_path,
        output,
        source,
        len(entries),
        sum(entry["status"] == "unqualified" for entry in entries),
    )


def _write_mc_qualification_evidence(
    directory: Path,
    evidence: _contracts.MonteCarloEvidence,
    *,
    connection: sqlite3.Connection,
    mixture_id: int,
    metadata: dict[str, str],
    sampling_plan: dict[str, object] | None,
    qualification_profile: str,
    qualification_error: str | None,
    source_database_sha256: str,
) -> None:
    """Keep failed anchors available to policy even without coefficient tables."""
    path = directory / "manifest.json"
    if path.is_file():
        manifest = json.loads(path.read_text(encoding="utf-8"))
    else:
        manifest = {
            "format_version": _contracts.FORMAT_VERSION,
            "stage": "build-tables",
            "status": "unqualified",
            "source": _contracts.MC_SOLVER,
            "physical_context": context_from_metadata(metadata),
            "hashes": {
                **provenance_hash_manifest(
                    metadata, source=_contracts.MC_SOLVER
                ),
                "source_database_sha256": source_database_sha256,
            },
            "mc_sampling_plan": sampling_plan,
            "mc_campaign": json.loads(metadata.get("mc_campaign_json", "null")),
            "mixture": {
                "mixture_id": mixture_id,
                "species": _repository._mixture_rows(connection, mixture_id),
            },
            "source_policy": {
                **_repository._field_source_policy(
                    connection,
                    solver=_contracts.MC_SOLVER,
                    mixture_id=mixture_id,
                    require_mc_solver_physics=(
                        qualification_profile
                        == _contracts.MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
                    ),
                ),
                "qualification_profile": qualification_profile,
            },
            "tables": {},
        }
    columns = tuple(dict.fromkeys(key for row in evidence.quality for key in row))
    entry = _write_csv(
        directory / MC_QUALIFICATION_TABLE,
        columns,
        evidence.quality,
        argument="E_over_N_Td",
    )
    entry["artifact_role"] = "all_planned_mc_anchor_qualification"
    manifest["tables"][MC_QUALIFICATION_TABLE] = entry
    manifest["mc_qualification"] = {
        "file": MC_QUALIFICATION_TABLE,
        "all_planned_anchors": len(evidence.quality),
        "qualified_anchors": len(evidence.transport_eligible),
        "coefficient_tables_available": qualification_error is None,
        "table_build_failure": qualification_error,
    }
    _write_json(path, manifest)


def _replace_directory(staging: Path, target: Path) -> None:
    """Replace a complete table generation while retaining rollback state."""

    if not target.exists():
        staging.replace(target)
        return
    backup = Path(
        tempfile.mkdtemp(
            prefix=f".{target.name}.backup-",
            dir=target.parent,
        )
    )
    backup.rmdir()
    target.replace(backup)
    try:
        staging.replace(target)
    except Exception:
        backup.replace(target)
        raise
    shutil.rmtree(backup)


def _build_dataset(
    connection: sqlite3.Connection,
    mixture_id: int,
    source: str,
    *,
    mc_sampling_plan: dict[str, object] | None,
    quality_thresholds: QualityThresholds,
    mc_qualification_profile: str,
    mc_evidence: _contracts.MonteCarloEvidence | None = None,
) -> _contracts.TableDataset:
    mixture = _repository._mixture_rows(connection, mixture_id)
    if source != _contracts.MC_SOLVER:
        cases = _repository._load_cases(connection, source, mixture_id)
        rates = _repository._load_rates(connection, source, mixture_id, cases)
        eedf = _repository._load_eedf(connection, source, mixture_id, cases)
        quality = _repository._load_quality(connection, source, mixture_id)
        source_policy = {
            "source": source,
            "postprocess": "none",
            "selection": "all_aggregate_points",
        }
        if source == _contracts.TWO_TERM_SOLVER:
            two_term_transport_kernel = (
                _repository._load_two_term_temporal_growth_transport_contract(
                    connection,
                    mixture_id,
                    cases,
                )
            )
            if two_term_transport_kernel is not None:
                source_policy["two_term_transport_kernel"] = two_term_transport_kernel
    else:
        cases, rates, eedf, quality, source_policy = (
            _mc_tables._build_qualified_monte_carlo(
                connection,
                mixture_id,
                mc_sampling_plan=mc_sampling_plan,
                quality_thresholds=quality_thresholds,
                qualification_profile=mc_qualification_profile,
                evidence=mc_evidence,
            )
        )
    rate_evidence = (
        _mc_tables._load_rate_evidence(
            connection,
            mixture_id,
            allowed_e={float(row["E_over_N_Td"]) for row in eedf},
        )
        if source == _contracts.MC_SOLVER
        else None
    )
    if source == _contracts.MC_SOLVER:
        assert rate_evidence is not None
        source_policy["rate_censoring_qualification"] = (
            _mc_tables._validate_mc_censored_rate_relevance(
                rates,
                rate_evidence,
                minimum_fraction=(
                    quality_thresholds.required_rate_min_process_peak_fraction
                ),
            )
        )
    if not cases:
        raise _contracts.TableBuildError(
            f"mixture {mixture_id} has no cases for {source}"
        )
    elastic_allowed = (
        {
            float(value)
            for value in source_policy.get("eedf_rate_eligible_E_over_N_Td", ())
        }
        if source == _contracts.MC_SOLVER
        else {float(case["E_over_N_Td"]) for case in cases}
    )
    elastic_cases = _repository._load_cases(
        connection,
        source,
        mixture_id,
        allowed_e=elastic_allowed,
    )
    elastic_energy_loss, elastic_energy_loss_metadata = (
        _energy_loss._load_elastic_energy_loss(
            connection,
            source,
            mixture_id,
            elastic_cases,
            allowed_e=elastic_allowed,
        )
    )
    source_policy.update(
        _repository._field_source_policy(
            connection,
            solver=source,
            mixture_id=mixture_id,
            require_mc_solver_physics=(
                source == _contracts.MC_SOLVER
                and mc_qualification_profile
                == _contracts.MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
            ),
        )
    )
    return _contracts.TableDataset(
        source=source,
        mixture_id=mixture_id,
        mixture=mixture,
        cases=cases,
        rates=rates,
        eedf=eedf,
        rate_evidence=rate_evidence,
        elastic_energy_loss=elastic_energy_loss,
        elastic_energy_loss_metadata=elastic_energy_loss_metadata,
        quality=quality,
        source_policy=source_policy,
    )


def _write_dataset_tables(
    dataset: _contracts.TableDataset,
    directory: Path,
    *,
    metadata: dict[str, str],
    source_quality_thresholds: QualityThresholds,
    source_quality_json: str,
    quality_thresholds: QualityThresholds,
    mc_sampling_plan: dict[str, object] | None,
    solver_qualification: SolverQualification | None,
    target_qualification: PropagatorTargetQualification | None,
    source_database_sha256: str,
) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    qualification_entry: dict[str, Any] | None = None
    if solver_qualification is not None:
        qualification_path = directory / PROPAGATOR_CORE_QUALIFICATION_FILE
        shutil.copy2(solver_qualification.source_path, qualification_path)
        qualification_entry = solver_qualification.manifest_entry(
            file=PROPAGATOR_CORE_QUALIFICATION_FILE
        )
    target_qualification_entry: dict[str, Any] | None = None
    if target_qualification is not None:
        target_path = directory / PROPAGATOR_TARGET_QUALIFICATION_FILE
        shutil.copy2(target_qualification.source_path, target_path)
        target_qualification_entry = target_qualification.manifest_entry(
            file=PROPAGATOR_TARGET_QUALIFICATION_FILE
        )
    _validate_positive_values(dataset)
    monotonic = strictly_monotonic(
        [float(case["mean_energy_eV"]) for case in dataset.cases]
    )
    monotonic_reason = None if monotonic else "mean_energy_not_strictly_monotonic"

    tables: dict[str, dict[str, Any]] = {}
    temporal_growth_transport = isinstance(
        dataset.source_policy.get("two_term_transport_kernel"), dict
    )
    transport_columns = (
        *_contracts.CASE_COLUMNS,
        *(
            _contracts.TWO_TERM_TEMPORAL_GROWTH_TRANSPORT_COLUMNS
            if temporal_growth_transport
            else ()
        ),
    )
    tables["mean_energy_vs_en.csv"] = _write_csv(
        directory / "mean_energy_vs_en.csv",
        ("E_over_N_Td", "E_over_N_V_m2", "mean_energy_eV"),
        dataset.cases,
        argument="E_over_N_Td",
    )
    tables["transport_vs_en.csv"] = _write_csv(
        directory / "transport_vs_en.csv",
        transport_columns,
        dataset.cases,
        argument="E_over_N_Td",
    )
    tables["rates_vs_en.csv"] = _write_csv(
        directory / "rates_vs_en.csv",
        _contracts.RATE_COLUMNS,
        dataset.rates,
        argument="E_over_N_Td",
    )
    tables["townsend_vs_en.csv"] = _write_csv(
        directory / "townsend_vs_en.csv",
        ("E_over_N_Td", "E_over_N_V_m2", "effective_townsend_m2"),
        dataset.cases,
        argument="E_over_N_Td",
    )
    tables["energy_loss.csv"] = _write_csv(
        directory / "energy_loss.csv",
        _contracts.ENERGY_LOSS_COLUMNS,
        dataset.rates,
        argument="E_over_N_Td",
    )
    quality_schema = quality_table_schema(dataset.source)
    tables["quality.csv"] = _write_csv(
        directory / "quality.csv",
        quality_schema.columns,
        dataset.quality,
        argument="E_over_N_Td",
    )
    tables["quality.csv"].update(
        {
            "schema": "solver_quality.v2",
            "common_columns": list(COMMON_QUALITY_COLUMNS),
            "solver_evidence": {
                "kind": quality_schema.evidence_kind,
                "columns": list(quality_schema.evidence_columns),
            },
        }
    )
    tables["eedf.csv"] = _write_csv(
        directory / "eedf.csv",
        _contracts.EEDF_COLUMNS,
        dataset.eedf,
        argument="electron_energy_eV,E_over_N_Td",
    )
    tables["eedf.csv"].update(
        {
            "artifact_role": "raw_swarm_eedf_audit_and_canonicalization_source",
            "sha256": _sha256_file(directory / "eedf.csv"),
        }
    )
    collision_kernels = _collision_kernels.active_collision_kernel_rows(
        metadata,
        dataset.rates,
        mixture=dataset.mixture,
    )
    if collision_kernels is not None:
        columns, rows = collision_kernels
        name = _contracts.COLLISION_RATE_KERNEL_TABLE
        tables[name] = _write_csv(
            directory / name,
            columns,
            rows,
            argument="electron_energy_eV,process",
        )
        tables[name].update(
            {
                "artifact_role": "canonical_eedf_projection_rate_kernels",
                "cross_sections_sha256": metadata["cross_sections_sha256"],
            }
        )
    if dataset.rate_evidence is not None:
        tables[_contracts.RATE_EVIDENCE_TABLE] = _write_csv(
            directory / _contracts.RATE_EVIDENCE_TABLE,
            _contracts.RATE_EVIDENCE_COLUMNS,
            dataset.rate_evidence,
            argument="E_over_N_Td,process",
        )
        tables[_contracts.RATE_EVIDENCE_TABLE].update(
            {
                "artifact_role": "raw_swarm_rate_evidence",
                "sha256": _sha256_file(directory / _contracts.RATE_EVIDENCE_TABLE),
                "statistics": {
                    "replicate_interval": "two_sided_student_t_95",
                    "zero_event_confidence": _contracts.ZERO_EVENT_CONFIDENCE,
                    "zero_event_upper_bound": (
                        "poisson_zero_count_over_pooled_target_exposure"
                    ),
                },
            }
        )
    if dataset.elastic_energy_loss is not None:
        if dataset.elastic_energy_loss_metadata is None:
            raise _contracts.TableBuildError(
                "elastic energy-loss rows lack their physics contract"
            )
        tables[_contracts.ELASTIC_ENERGY_LOSS_TABLE] = _write_csv(
            directory / _contracts.ELASTIC_ENERGY_LOSS_TABLE,
            _contracts.ELASTIC_ENERGY_LOSS_COLUMNS,
            dataset.elastic_energy_loss,
            argument="mean_energy_eV",
        )
        tables[_contracts.ELASTIC_ENERGY_LOSS_TABLE].update(
            {
                "artifact_role": ("canonical_comsol_elastic_energy_loss_input"),
                "sha256": _sha256_file(
                    directory / _contracts.ELASTIC_ENERGY_LOSS_TABLE
                ),
                "physics_contract": dataset.elastic_energy_loss_metadata,
            }
        )
    if monotonic:
        tables["transport_vs_mean_energy.csv"] = _write_csv(
            directory / "transport_vs_mean_energy.csv",
            (
                "mean_energy_eV",
                "E_over_N_Td",
                "E_over_N_V_m2",
                *_contracts.CASE_COLUMNS[3:],
                *(
                    _contracts.TWO_TERM_TEMPORAL_GROWTH_TRANSPORT_COLUMNS
                    if temporal_growth_transport
                    else ()
                ),
            ),
            dataset.cases,
            argument="mean_energy_eV",
        )
        tables["rates_vs_mean_energy.csv"] = _write_csv(
            directory / "rates_vs_mean_energy.csv",
            ("mean_energy_eV", *_contracts.RATE_COLUMNS),
            dataset.rates,
            argument="mean_energy_eV",
        )
    _write_json(
        directory / "manifest.json",
        {
            "format_version": _contracts.FORMAT_VERSION,
            "stage": "build-tables",
            "source": dataset.source,
            "physical_context": context_from_metadata(metadata),
            "mc_campaign": json.loads(metadata.get("mc_campaign_json", "null")),
            "hashes": {
                **provenance_hash_manifest(metadata, source=dataset.source),
                "source_database_sha256": source_database_sha256,
            },
            "solver_qualification": qualification_entry,
            "target_qualification": target_qualification_entry,
            "mc_sampling_plan": mc_sampling_plan,
            "quality_thresholds": quality_thresholds_payload(quality_thresholds),
            "quality_policy": quality_policy_provenance(
                source_quality_thresholds,
                quality_thresholds,
                source_json=source_quality_json,
            ),
            "mixture": {
                "mixture_id": dataset.mixture_id,
                "species": dataset.mixture,
            },
            "source_policy": dataset.source_policy,
            "valid_ranges": _valid_ranges(dataset.cases),
            "units": _manifest_units(tables),
            "table_argument": {
                "primary": "E_over_N_Td",
                "secondary": "mean_energy_eV" if monotonic else None,
                "eedf": "electron_energy_eV,E_over_N_Td",
            },
            "monotonicity": {
                "mean_energy_strictly_monotonic": monotonic,
                "reason": monotonic_reason,
            },
            "quality_summary": _quality_summary(dataset.quality),
            "tables": tables,
        },
    )


def _validate_positive_values(dataset: _contracts.TableDataset) -> None:
    positive_case = [
        "mean_energy_eV",
        "reduced_mobility_m2_V_s_m3",
    ]
    if dataset.source != _contracts.PROPAGATOR_SOLVER:
        positive_case.extend(
            [
                "reduced_diffusion_L_m2_s_m3",
                "reduced_diffusion_T_m2_s_m3",
            ]
        )
    optional_positive_case = [
        "reduced_electron_energy_mobility_m2_V_s_m3",
        "reduced_electron_energy_diffusion_m2_s_m3",
        "reduced_electron_energy_diffusion_L_m2_s_m3",
        "reduced_electron_energy_diffusion_T_m2_s_m3",
    ]
    positive_rate = [
        "rate_coefficient_m3_s",
        "mixture_weighted_rate_m3_s",
        "reduced_townsend_m2",
        "mixture_weighted_reduced_townsend_m2",
    ]
    for case in dataset.cases:
        for name in positive_case:
            _require_nonnegative_finite(case.get(name), name)
        if dataset.source == _contracts.MC_SOLVER:
            for name in (
                "reduced_electron_energy_mobility_m2_V_s_m3",
                "reduced_electron_energy_diffusion_L_m2_s_m3",
                "reduced_electron_energy_diffusion_T_m2_s_m3",
            ):
                _require_positive_finite(case.get(name), name)
        for name in optional_positive_case:
            if case.get(name) is not None:
                _require_nonnegative_finite(case.get(name), name)
    for rate in dataset.rates:
        for name in positive_rate:
            _require_nonnegative_finite(rate.get(name), name)
    for row in dataset.eedf:
        _require_nonnegative_finite(row.get("eedf"), "raw EEDF eedf")
        _require_nonnegative_finite(
            row.get("energy_width_eV"), "raw EEDF energy_width_eV"
        )


def _require_nonnegative_finite(value: object, name: str) -> None:
    if value is None:
        raise _contracts.TableBuildError(f"missing positive value {name}")
    number = float(value)
    if not math.isfinite(number) or number < 0.0:
        raise _contracts.TableBuildError(f"invalid positive value {name}: {value}")


def _require_positive_finite(value: object, name: str) -> None:
    if value is None:
        raise _contracts.TableBuildError(f"missing positive value {name}")
    number = float(value)
    if not math.isfinite(number) or number <= 0.0:
        raise _contracts.TableBuildError(f"invalid positive value {name}: {value}")


def _write_csv(
    path: Path,
    columns: tuple[str, ...],
    rows: list[dict[str, Any]],
    *,
    argument: str,
) -> dict[str, Any]:
    _write_csv_rows(path, columns, rows)
    return {
        "columns": list(columns),
        "sha256": _sha256_file(path),
        "units": {
            column: _contracts.UNITS[column]
            for column in columns
            if column in _contracts.UNITS
        },
        "argument": argument,
    }


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _valid_ranges(cases: list[dict[str, Any]]) -> dict[str, list[float] | None]:
    return {
        "E_over_N_Td": finite_range(case.get("E_over_N_Td") for case in cases),
        "E_over_N_V_m2": finite_range(case.get("E_over_N_V_m2") for case in cases),
        "mean_energy_eV": finite_range(case.get("mean_energy_eV") for case in cases),
    }


def _quality_summary(quality: list[dict[str, Any]]) -> dict[str, Any]:
    failed = [
        row
        for row in quality
        if int(
            row.get(
                "active_closure_quality_passed",
                row.get("passed", 0),
            )
        )
        == 0
    ]
    return {
        "passed": not failed,
        "failed_points": len(failed),
        "total_points": len(quality),
    }


def _manifest_units(tables: dict[str, dict[str, Any]]) -> dict[str, str]:
    columns = {
        column
        for table in tables.values()
        for column in table.get("columns", [])
        if column in _contracts.UNITS
    }
    return {column: _contracts.UNITS[column] for column in sorted(columns)}
