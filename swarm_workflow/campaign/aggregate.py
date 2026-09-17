"""Aggregate workflow SQLite results for table-building stages."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from hashlib import sha256
import json
import math
from pathlib import Path
import sqlite3
from typing import Iterable

import numpy as np

import electron_swarm.solvers.monte_carlo.evidence as _mc_evidence

from .._io import write_json as _write_json
from ..quality import policy as _quality_policy
from . import quality as _quality_database
from .statistics import (
    ScalarStats as _ScalarStats,
    conservative_rebin_probability_mass as _conservative_rebin_probability_mass,
    energy_edges_from_cells,
    energy_edges_from_nodes,
    summarize_replicates as _summarize_replicates,
)
from .store import (
    mc_sampling_plan_provenance,
    provenance_hash_manifest,
    validate_workflow_schema,
)
from .repository import read_metadata


CASE_SCALARS = (
    "mean_energy_eV",
    "drift_velocity_m_s",
    "reduced_mobility_m2_V_s_m3",
    "reduced_diffusion_L_m2_s_m3",
    "reduced_diffusion_T_m2_s_m3",
    "reduced_electron_energy_mobility_m2_V_s_m3",
    "reduced_electron_energy_diffusion_m2_s_m3",
    "reduced_electron_energy_diffusion_L_m2_s_m3",
    "reduced_electron_energy_diffusion_T_m2_s_m3",
    "effective_townsend_m2",
)
OPTIONAL_CASE_SCALARS = {
    "reduced_diffusion_L_m2_s_m3",
    "reduced_diffusion_T_m2_s_m3",
    "reduced_electron_energy_mobility_m2_V_s_m3",
    "reduced_electron_energy_diffusion_m2_s_m3",
    "reduced_electron_energy_diffusion_L_m2_s_m3",
    "reduced_electron_energy_diffusion_T_m2_s_m3",
}
RATE_SCALAR = "rate_coefficient_m3_s"
MC_SOLVER = "monte_carlo"
FINITE_VOLUME_EEDF_SOLVERS = frozenset({MC_SOLVER, "propagator"})
RateKey = tuple[str, str, str, str, float | None]

AGGREGATE_CSVS = (
    "aggregate_scalars.csv",
    "aggregate_eedf_bins.csv",
    "aggregate_quality.csv",
)

__all__ = [
    "AGGREGATE_CSVS",
    "CASE_SCALARS",
    "AggregateSummary",
    "aggregate_database",
    "aggregate_workflow_results",
    "stable_mc_seed",
    "write_aggregate_outputs",
]


@dataclass(frozen=True, slots=True)
class AggregateSummary:
    database_path: Path
    output_directory: Path | None
    solvers: tuple[str, ...]
    mixtures: int
    e_over_n_points: int


def stable_mc_seed(
    *,
    base_seed: int | None,
    mixture_id: int,
    e_over_n_Td: float,
    replicate: int | None = None,
    solver: str = MC_SOLVER,
    replica_index: int | None = None,
) -> int | None:
    """Return an execution-order-independent Monte Carlo seed."""

    if base_seed is None:
        return None
    if replicate is None:
        if replica_index is None:
            raise TypeError("stable_mc_seed requires replicate")
        replicate = replica_index
    material = {
        "seed_derivation": _mc_evidence.MC_WORKFLOW_REPLICA_SEED_DERIVATION,
        "base_seed": int(base_seed),
        "mixture_id": int(mixture_id),
        "E_over_N_Td": float(e_over_n_Td),
        "replicate": int(replicate),
        "solver": str(solver),
    }
    encoded = json.dumps(material, sort_keys=True, separators=(",", ":")).encode(
        "utf-8"
    )
    return int.from_bytes(sha256(encoded).digest()[:4], "big", signed=False)


def aggregate_database(
    database_path: str | Path,
    output_directory: str | Path | None = None,
    thresholds: _quality_policy.QualityThresholds | None = None,
) -> AggregateSummary:
    db_path = Path(database_path)
    connection = sqlite3.connect(db_path)
    try:
        connection.row_factory = sqlite3.Row
        aggregate_workflow_results(connection, thresholds)
        output_path = Path(output_directory) if output_directory is not None else None
        if output_path is not None:
            write_aggregate_outputs(connection, output_path, database_path=db_path)
        solvers = tuple(
            str(row["solver"])
            for row in connection.execute(
                "SELECT DISTINCT solver FROM aggregate_quality ORDER BY solver"
            )
        )
        mixtures = int(
            connection.execute(
                "SELECT COUNT(DISTINCT mixture_id) FROM aggregate_quality"
            ).fetchone()[0]
            or 0
        )
        e_points = int(
            connection.execute(
                "SELECT COUNT(*) FROM aggregate_quality"
            ).fetchone()[0]
            or 0
        )
    finally:
        connection.close()
    return AggregateSummary(db_path, output_path, solvers, mixtures, e_points)


def aggregate_workflow_results(
    connection: sqlite3.Connection,
    thresholds: _quality_policy.QualityThresholds | None = None,
    *,
    solvers: tuple[str, ...] | None = None,
    _allow_evaluation_policy_override: bool = False,
) -> None:
    connection.row_factory = sqlite3.Row
    validate_workflow_schema(connection)
    source_thresholds = _quality_database.resolve_quality_thresholds(connection)
    source_policy_json = _quality_database.source_quality_thresholds_json(connection)
    if thresholds is None:
        thresholds = _quality_database.resolve_evaluation_quality_thresholds(connection)
    elif not _allow_evaluation_policy_override:
        _quality_database.resolve_quality_thresholds(connection, thresholds)
    policy_reevaluated = (
        _quality_policy.quality_thresholds_json(source_thresholds)
        != _quality_policy.quality_thresholds_json(thresholds)
    )
    rate_mean_peaks = _rate_mean_peaks(connection, solvers=solvers)
    _reset_tables(connection)
    solver_filter = ""
    params: tuple[object, ...] = ()
    if solvers is not None:
        placeholders = ", ".join("?" for _item in solvers)
        solver_filter = f"WHERE solver IN ({placeholders})"
        params = tuple(solvers)
    groups = connection.execute(
        f"""
        SELECT solver, mixture_id, e_over_n_Td,
               COUNT(DISTINCT replicate) AS replicate_count
        FROM cases
        {solver_filter}
        GROUP BY solver, mixture_id, e_over_n_Td
        ORDER BY solver, mixture_id, e_over_n_Td
        """,
        params,
    ).fetchall()

    with connection:
        for group in groups:
            solver = str(group["solver"])
            mixture_id = int(group["mixture_id"])
            e_over_n_Td = float(group["e_over_n_Td"])
            replicate_count = int(group["replicate_count"])
            require_uncertainty = solver == MC_SOLVER
            case_stats = _case_scalar_stats(connection, solver, mixture_id, e_over_n_Td)
            rate_stats, rate_metadata = _rate_scalar_stats(
                connection, solver, mixture_id, e_over_n_Td
            )

            for name, stats in case_stats.items():
                _insert_scalar(
                    connection,
                    solver=solver,
                    mixture_id=mixture_id,
                    e_over_n_Td=e_over_n_Td,
                    scalar_group="case",
                    scalar_name=name,
                    stats=stats,
                )

            for key, stats in rate_stats.items():
                species, process, process_type, threshold_key, threshold = key
                metadata = rate_metadata[key]
                _insert_scalar(
                    connection,
                    solver=solver,
                    mixture_id=mixture_id,
                    e_over_n_Td=e_over_n_Td,
                    scalar_group="rate",
                    scalar_name=RATE_SCALAR,
                    species=species,
                    process=process,
                    process_type=process_type,
                    threshold_key=threshold_key,
                    threshold_eV=threshold,
                    target_species_fraction=metadata["target_species_fraction"],
                    energy_loss_eV=metadata["energy_loss_eV"],
                    stats=stats,
                )

            eedf_normalization_error, eedf_finite_nonnegative = _aggregate_eedf(
                connection, solver, mixture_id, e_over_n_Td
            )
            (
                passed,
                failure_reasons,
                details,
                uncertainty_available,
            ) = _quality_for_group(
                case_stats=case_stats,
                rate_stats=rate_stats,
                rate_mean_peaks=rate_mean_peaks.get(
                    (solver, mixture_id), {}
                ),
                eedf_normalization_error=eedf_normalization_error,
                eedf_finite_nonnegative=eedf_finite_nonnegative,
                thresholds=thresholds,
                require_uncertainty=require_uncertainty,
            )
            connection.execute(
                """
                INSERT INTO aggregate_quality(
                    mixture_id, solver, e_over_n_Td, passed,
                    failure_reasons_json, source_thresholds_json,
                    thresholds_json, quality_policy_reevaluated, mobility_rse,
                    diffusion_L_rse, diffusion_T_rse, max_major_rate_rse,
                    energy_mobility_rse, energy_diffusion_L_rse,
                    energy_diffusion_T_rse,
                    eedf_normalization_error, valid_replicates, n,
                    uncertainty_available
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    mixture_id,
                    solver,
                    e_over_n_Td,
                    int(passed),
                    _json_dumps(list(failure_reasons)),
                    source_policy_json,
                    _quality_policy.quality_thresholds_json(thresholds),
                    int(policy_reevaluated),
                    _finite_or_none(details["mobility_rse"]),
                    _finite_or_none(details["diffusion_L_rse"]),
                    _finite_or_none(details["diffusion_T_rse"]),
                    _finite_or_none(details["max_major_rate_rse"]),
                    _finite_or_none(details["energy_mobility_rse"]),
                    _finite_or_none(details["energy_diffusion_L_rse"]),
                    _finite_or_none(details["energy_diffusion_T_rse"]),
                    _finite_or_none(details["eedf_normalization_error"]),
                    replicate_count,
                    replicate_count,
                    int(uncertainty_available),
                ),
            )


def write_aggregate_outputs(
    connection: sqlite3.Connection,
    output_directory: str | Path,
    *,
    database_path: str | Path | None = None,
) -> None:
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    for table_name, filename in [
        ("aggregate_scalars", "aggregate_scalars.csv"),
        ("aggregate_eedf_bins", "aggregate_eedf_bins.csv"),
        ("aggregate_quality", "aggregate_quality.csv"),
    ]:
        select = (
            "SELECT *, 'raw_replica_statistical_screen' AS quality_scope "
            "FROM aggregate_quality"
            if table_name == "aggregate_quality"
            else f"SELECT * FROM {table_name}"
        )
        _write_query_csv(
            connection,
            output / filename,
            f"{select} ORDER BY { _order_columns(table_name) }",
        )
    metadata = read_metadata(connection)
    source_thresholds = _quality_database.resolve_quality_thresholds(connection)
    source_policy_json = _quality_database.source_quality_thresholds_json(connection)
    evaluation_thresholds = _quality_database.resolve_evaluation_quality_thresholds(connection)
    _quality_database.validate_aggregate_quality_thresholds(
        connection, evaluation_thresholds
    )
    solvers = [
        str(row["solver"])
        for row in connection.execute(
            "SELECT DISTINCT solver FROM aggregate_quality ORDER BY solver"
        )
    ]
    _write_json(
        output / "manifest.json",
        {
            "format_version": 1,
            "stage": "aggregate",
            "database_path": str(database_path) if database_path is not None else None,
            "hashes": provenance_hash_manifest(metadata),
            "mc_sampling_plan": mc_sampling_plan_provenance(metadata),
            "quality_thresholds": _quality_policy.quality_thresholds_payload(
                evaluation_thresholds
            ),
            "quality_policy": _quality_policy.quality_policy_provenance(
                source_thresholds,
                evaluation_thresholds,
                source_json=source_policy_json,
            ),
            "quality_scope": {
                "aggregate_quality.csv": "raw_replica_statistical_screen",
                "comsol_input_qualification": (
                    "performed_later_by_build_tables_with_solver_specific_"
                    "diagnostics"
                ),
            },
            "outputs": list(AGGREGATE_CSVS),
            "solvers": solvers,
            "mixtures": [
                int(row["mixture_id"])
                for row in connection.execute(
                    "SELECT DISTINCT mixture_id FROM aggregate_quality ORDER BY mixture_id"
                )
            ],
            "e_over_n_points": int(
                connection.execute("SELECT COUNT(*) FROM aggregate_quality").fetchone()[
                    0
                ]
                or 0
            ),
        },
    )


def _reset_tables(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        DROP TABLE IF EXISTS aggregate_scalars;
        DROP TABLE IF EXISTS aggregate_eedf_bins;
        DROP TABLE IF EXISTS aggregate_quality;

        CREATE TABLE aggregate_scalars (
            mixture_id INTEGER NOT NULL,
            solver TEXT NOT NULL,
            e_over_n_Td REAL NOT NULL,
            scalar_group TEXT NOT NULL,
            scalar_name TEXT NOT NULL,
            species TEXT NOT NULL,
            process TEXT NOT NULL,
            process_type TEXT NOT NULL,
            threshold_key TEXT NOT NULL,
            threshold_eV REAL,
            target_species_fraction REAL,
            energy_loss_eV REAL,
            mean REAL,
            sample_stddev REAL,
            standard_error REAL,
            relative_standard_error REAL,
            ci95_low REAL,
            ci95_high REAL,
            ci95_critical_value REAL,
            estimate_status TEXT NOT NULL,
            valid_replicates INTEGER NOT NULL,
            n INTEGER NOT NULL,
            uncertainty_available INTEGER NOT NULL,
            PRIMARY KEY (
                mixture_id, solver, e_over_n_Td, scalar_group, scalar_name,
                species, process, process_type, threshold_key
            )
        );
        CREATE TABLE aggregate_eedf_bins (
            mixture_id INTEGER NOT NULL,
            solver TEXT NOT NULL,
            e_over_n_Td REAL NOT NULL,
            bin_index INTEGER NOT NULL,
            energy_left_eV REAL NOT NULL,
            energy_right_eV REAL NOT NULL,
            energy_eV REAL NOT NULL,
            energy_width_eV REAL NOT NULL,
            mean_probability_mass REAL NOT NULL,
            probability_mass_sample_stddev REAL,
            probability_mass_standard_error REAL,
            eedf REAL NOT NULL,
            pooled_effective_sample_count REAL,
            valid_replicates INTEGER NOT NULL,
            n INTEGER NOT NULL,
            uncertainty_available INTEGER NOT NULL,
            PRIMARY KEY (mixture_id, solver, e_over_n_Td, bin_index)
        );
        CREATE TABLE aggregate_quality (
            mixture_id INTEGER NOT NULL,
            solver TEXT NOT NULL,
            e_over_n_Td REAL NOT NULL,
            passed INTEGER NOT NULL,
            failure_reasons_json TEXT NOT NULL,
            source_thresholds_json TEXT NOT NULL,
            thresholds_json TEXT NOT NULL,
            quality_policy_reevaluated INTEGER NOT NULL,
            mobility_rse REAL,
            diffusion_L_rse REAL,
            diffusion_T_rse REAL,
            energy_mobility_rse REAL,
            energy_diffusion_L_rse REAL,
            energy_diffusion_T_rse REAL,
            max_major_rate_rse REAL,
            eedf_normalization_error REAL,
            valid_replicates INTEGER NOT NULL,
            n INTEGER NOT NULL,
            uncertainty_available INTEGER NOT NULL,
            PRIMARY KEY (mixture_id, solver, e_over_n_Td)
        );
        """
    )


def _rate_mean_peaks(
    connection: sqlite3.Connection,
    *,
    solvers: tuple[str, ...] | None,
) -> dict[tuple[str, int], dict[RateKey, float]]:
    """Calculate each reaction's peak raw replica mean across E/N."""

    solver_filter = ""
    params: tuple[object, ...] = ()
    if solvers is not None:
        placeholders = ", ".join("?" for _solver in solvers)
        solver_filter = f"WHERE solver IN ({placeholders})"
        params = tuple(solvers)
    rows = connection.execute(
        f"""
        SELECT solver, mixture_id, e_over_n_Td, species, process,
               process_type, threshold_eV,
               AVG(rate_coefficient_m3_s) AS raw_mean
        FROM rates
        {solver_filter}
        GROUP BY solver, mixture_id, e_over_n_Td, species, process,
                 process_type, threshold_eV
        """,
        params,
    ).fetchall()
    peaks: dict[tuple[str, int], dict[RateKey, float]] = {}
    for row in rows:
        mean = _finite_or_none(row["raw_mean"])
        if mean is None or mean < 0.0:
            continue
        group = (str(row["solver"]), int(row["mixture_id"]))
        key = _rate_key(row)
        by_reaction = peaks.setdefault(group, {})
        by_reaction[key] = max(mean, by_reaction.get(key, 0.0))
    return peaks


def _case_scalar_stats(
    connection: sqlite3.Connection,
    solver: str,
    mixture_id: int,
    e_over_n_Td: float,
) -> dict[str, _ScalarStats]:
    rows = connection.execute(
        f"""
        SELECT {", ".join(CASE_SCALARS)}
        FROM cases
        WHERE solver = ? AND mixture_id = ? AND e_over_n_Td = ?
        ORDER BY replicate
        """,
        (solver, mixture_id, e_over_n_Td),
    ).fetchall()
    return {
        scalar_name: _summarize_replicates(row[scalar_name] for row in rows)
        for scalar_name in CASE_SCALARS
    }


def _rate_scalar_stats(
    connection: sqlite3.Connection,
    solver: str,
    mixture_id: int,
    e_over_n_Td: float,
) -> tuple[
    dict[tuple[str, str, str, str, float | None], _ScalarStats],
    dict[tuple[str, str, str, str, float | None], dict[str, float | None]],
]:
    rows = connection.execute(
        """
        SELECT replicate, species, process, process_type, threshold_eV,
               rate_coefficient_m3_s, target_species_fraction, energy_loss_eV
        FROM rates
        WHERE solver = ? AND mixture_id = ? AND e_over_n_Td = ?
        ORDER BY species, process, process_type, threshold_eV, replicate, rate_index
        """,
        (solver, mixture_id, e_over_n_Td),
    ).fetchall()
    grouped: dict[tuple[str, str, str, str, float | None], list[float]] = {}
    target_fractions: dict[tuple[str, str, str, str, float | None], list[float | None]] = {}
    energy_losses: dict[tuple[str, str, str, str, float | None], list[float | None]] = {}
    for row in rows:
        key = _rate_key(row)
        grouped.setdefault(key, []).append(float(row[RATE_SCALAR]))
        target_fractions.setdefault(key, []).append(_finite_or_none(row["target_species_fraction"]))
        energy_losses.setdefault(key, []).append(_finite_or_none(row["energy_loss_eV"]))

    metadata: dict[tuple[str, str, str, str, float | None], dict[str, float | None]] = {}
    for key in grouped:
        target_fraction = _constant_optional_number(target_fractions[key])
        energy_loss = _constant_optional_number(energy_losses[key])
        metadata[key] = {
            "target_species_fraction": target_fraction,
            "energy_loss_eV": energy_loss,
        }
    return {
        key: _summarize_replicates(
            values, zero_is_censored=(solver == MC_SOLVER)
        )
        for key, values in grouped.items()
    }, metadata


def _aggregate_eedf(
    connection: sqlite3.Connection,
    solver: str,
    mixture_id: int,
    e_over_n_Td: float,
) -> tuple[float | None, bool]:
    rows = connection.execute(
        """
        SELECT replicate, bin_index, energy_eV, energy_width_eV, eedf,
               effective_sample_count
        FROM eedf_bins
        WHERE solver = ? AND mixture_id = ? AND e_over_n_Td = ?
        ORDER BY replicate, bin_index
        """,
        (solver, mixture_id, e_over_n_Td),
    ).fetchall()
    by_replicate: dict[int, list[sqlite3.Row]] = {}
    for row in rows:
        by_replicate.setdefault(int(row["replicate"]), []).append(row)
    if not by_replicate:
        return None, False

    replicate_edges: dict[int, np.ndarray] = {}
    replicate_masses: dict[int, np.ndarray] = {}
    replicate_effective_counts: dict[int, np.ndarray] = {}
    for replicate, replicate_rows in by_replicate.items():
        energy = np.asarray([row["energy_eV"] for row in replicate_rows], dtype=float)
        widths = np.asarray([row["energy_width_eV"] for row in replicate_rows], dtype=float)
        eedf = np.asarray([row["eedf"] for row in replicate_rows], dtype=float)
        effective_counts = np.asarray(
            [
                float(row["effective_sample_count"])
                if row["effective_sample_count"] is not None
                else math.nan
                for row in replicate_rows
            ],
            dtype=float,
        )
        if (
            len(energy) == 0
            or np.any(~np.isfinite(energy))
            or np.any(~np.isfinite(widths))
            or np.any(widths <= 0.0)
            or np.any(~np.isfinite(eedf))
        ):
            raise ValueError(
                "EEDF cells must contain finite centers, widths, and densities "
                f"for {solver} mixture {mixture_id} at {e_over_n_Td:g} Td "
                f"replicate {replicate}"
            )
        try:
            edges = (
                energy_edges_from_cells(energy, widths)
                if solver in FINITE_VOLUME_EEDF_SOLVERS
                else energy_edges_from_nodes(energy)
            )
        except ValueError as exc:
            raise ValueError(
                "Invalid finite-volume EEDF grid for "
                f"{solver} mixture {mixture_id} at {e_over_n_Td:g} Td "
                f"replicate {replicate}: {exc}"
            ) from exc
        reconstructed_widths = np.diff(edges)
        # ``center +/- width/2`` necessarily loses a few ulps when a narrow
        # physical cell sits at a large absolute energy.  Use the same
        # scale-aware IEEE-754 bound as ``energy_edges_from_cells`` instead of
        # a second, width-only tolerance that rejects an otherwise valid
        # contiguous finite-volume grid.
        roundtrip_scale = np.maximum.reduce(
            (np.ones_like(energy), np.abs(energy), widths)
        )
        roundtrip_tolerance = (
            64.0 * np.finfo(float).eps * roundtrip_scale
        )
        if np.any(np.abs(reconstructed_widths - widths) > roundtrip_tolerance):
            raise ValueError(
                "Stored EEDF widths do not match their grid representation for "
                f"{solver} mixture {mixture_id} at {e_over_n_Td:g} Td "
                f"replicate {replicate}"
            )
        mass = eedf * widths
        replicate_edges[replicate] = edges
        replicate_masses[replicate] = mass
        replicate_effective_counts[replicate] = effective_counts
    if not replicate_edges:
        return None, False

    ordered_replicates = sorted(replicate_edges)
    reference_edges = replicate_edges[ordered_replicates[0]]
    common_grid = all(
        np.array_equal(reference_edges, replicate_edges[replicate])
        for replicate in ordered_replicates[1:]
    )
    if common_grid:
        target_edges = reference_edges
        rebinned = np.vstack(
            [replicate_masses[replicate] for replicate in ordered_replicates]
        )
    else:
        target_edges = _union_energy_edges(
            replicate_edges[replicate] for replicate in ordered_replicates
        )
        rebinned = np.vstack(
            [
                _conservative_rebin_probability_mass(
                    replicate_edges[replicate],
                    replicate_masses[replicate],
                    target_edges,
                )
                for replicate in ordered_replicates
            ]
        )

    pooled_effective_counts: np.ndarray | None = None
    if all(
        np.all(np.isfinite(replicate_effective_counts[replicate]))
        and np.all(replicate_effective_counts[replicate] >= 0.0)
        for replicate in ordered_replicates
    ):
        if common_grid:
            support = np.vstack(
                [
                    replicate_effective_counts[replicate]
                    for replicate in ordered_replicates
                ]
            )
        else:
            support = np.vstack(
                [
                    _conservative_rebin_probability_mass(
                        replicate_edges[replicate],
                        replicate_effective_counts[replicate],
                        target_edges,
                    )
                    for replicate in ordered_replicates
                ]
            )
        pooled_effective_counts = np.sum(support, axis=0)

    stats = [
        _summarize_replicates(rebinned[:, index])
        for index in range(rebinned.shape[1])
    ]
    widths = np.diff(target_edges)
    mean_mass = np.asarray([stat.mean or 0.0 for stat in stats], dtype=float)
    normalization_error = float(abs(np.sum(mean_mass) - 1.0))
    finite_nonnegative = bool(
        np.all(np.isfinite(rebinned))
        and np.all(rebinned >= 0.0)
        and np.all(np.isfinite(mean_mass))
        and np.all(mean_mass >= 0.0)
        and np.all(np.isfinite(widths))
        and np.all(widths > 0.0)
    )
    connection.executemany(
        """
        INSERT INTO aggregate_eedf_bins(
            mixture_id, solver, e_over_n_Td, bin_index, energy_left_eV,
            energy_right_eV, energy_eV, energy_width_eV, mean_probability_mass,
            probability_mass_sample_stddev, probability_mass_standard_error,
            eedf, pooled_effective_sample_count, valid_replicates, n,
            uncertainty_available
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        [
            (
                mixture_id,
                solver,
                e_over_n_Td,
                index,
                float(target_edges[index]),
                float(target_edges[index + 1]),
                float(0.5 * (target_edges[index] + target_edges[index + 1])),
                float(widths[index]),
                float(mean_mass[index]),
                _finite_or_none(stats[index].sample_stddev),
                _finite_or_none(stats[index].standard_error),
                float(mean_mass[index] / widths[index]),
                (
                    float(pooled_effective_counts[index])
                    if pooled_effective_counts is not None
                    else None
                ),
                stats[index].valid_replicates,
                stats[index].valid_replicates,
                int(stats[index].uncertainty_available),
            )
            for index in range(len(widths))
            if widths[index] > 0.0
        ],
    )
    return normalization_error, finite_nonnegative


def _union_energy_edges(edge_sets: Iterable[np.ndarray]) -> np.ndarray:
    """Build a stable union of exact cell endpoints for conservative rebinning."""

    values = np.sort(
        np.concatenate([np.asarray(edges, dtype=float) for edges in edge_sets])
    )
    if len(values) < 2 or np.any(~np.isfinite(values)):
        raise ValueError("EEDF edge union requires finite cell endpoints")
    merged = [float(values[0])]
    for raw_value in values[1:]:
        value = float(raw_value)
        tolerance = 64.0 * np.finfo(float).eps * max(
            1.0, abs(merged[-1]), abs(value)
        )
        if value - merged[-1] <= tolerance:
            merged[-1] = 0.5 * (merged[-1] + value)
        else:
            merged.append(value)
    result = np.asarray(merged, dtype=float)
    if len(result) < 2 or np.any(np.diff(result) <= 0.0):
        raise ValueError("EEDF edge union must be strictly increasing")
    return result


def _quality_for_group(
    *,
    case_stats: dict[str, _ScalarStats],
    rate_stats: dict[RateKey, _ScalarStats],
    rate_mean_peaks: dict[RateKey, float],
    eedf_normalization_error: float | None,
    eedf_finite_nonnegative: bool,
    thresholds: _quality_policy.QualityThresholds,
    require_uncertainty: bool,
) -> tuple[bool, tuple[str, ...], dict[str, float | None], bool]:
    failures: list[str] = []
    finite_nonnegative = True

    for name, stats in case_stats.items():
        if stats.mean is None and name in OPTIONAL_CASE_SCALARS:
            continue
        if stats.mean is None or not math.isfinite(stats.mean) or stats.mean < 0.0:
            finite_nonnegative = False
            failures.append(f"{name}_not_finite_nonnegative")

    mobility = case_stats["reduced_mobility_m2_V_s_m3"]
    diffusion_l = case_stats["reduced_diffusion_L_m2_s_m3"]
    diffusion_t = case_stats["reduced_diffusion_T_m2_s_m3"]
    energy_mobility = case_stats[
        "reduced_electron_energy_mobility_m2_V_s_m3"
    ]
    energy_diffusion_l = case_stats[
        "reduced_electron_energy_diffusion_L_m2_s_m3"
    ]
    energy_diffusion_t = case_stats[
        "reduced_electron_energy_diffusion_T_m2_s_m3"
    ]
    if not _scalar_meets_quality(
        mobility,
        rse_limit=thresholds.mobility_rse,
        require_uncertainty=require_uncertainty,
    ):
        failures.append("reduced_mobility_rse_unavailable_or_exceeds_threshold")
    if (
        diffusion_l.mean is not None or require_uncertainty
    ) and not _scalar_meets_quality(
        diffusion_l,
        rse_limit=thresholds.diffusion_rse,
        require_uncertainty=require_uncertainty,
    ):
        failures.append("reduced_diffusion_L_rse_unavailable_or_exceeds_threshold")
    if (
        diffusion_t.mean is not None or require_uncertainty
    ) and not _scalar_meets_quality(
        diffusion_t,
        rse_limit=thresholds.diffusion_rse,
        require_uncertainty=require_uncertainty,
    ):
        failures.append("reduced_diffusion_T_rse_unavailable_or_exceeds_threshold")
    energy_stats = (
        (energy_mobility, thresholds.mobility_rse, "energy_mobility"),
        (energy_diffusion_l, thresholds.diffusion_rse, "energy_diffusion_L"),
        (energy_diffusion_t, thresholds.diffusion_rse, "energy_diffusion_T"),
    )
    for stats, limit, label in energy_stats:
        if stats.mean is None:
            continue
        if not _scalar_meets_quality(
            stats,
            rse_limit=limit,
            require_uncertainty=require_uncertainty,
        ):
            failures.append(f"reduced_{label}_rse_unavailable_or_exceeds_threshold")

    rate_means = [
        stats.mean
        for stats in rate_stats.values()
        if stats.mean is not None and math.isfinite(stats.mean) and stats.mean >= 0.0
    ]
    for key, stats in rate_stats.items():
        species, process, _process_type, _threshold_key_value, _threshold = key
        if stats.mean is None or not math.isfinite(stats.mean) or stats.mean < 0.0:
            finite_nonnegative = False
            failures.append(f"rate_not_finite_nonnegative:{species}:{process}")
    max_rate = max(rate_means) if rate_means else 0.0
    cutoff = max_rate * thresholds.major_rate_fraction
    major_rses: list[float] = []
    required_uncertainty_available = True
    for key, stats in rate_stats.items():
        species, process, process_type, _threshold_key_value, _threshold = key
        required_limit = _required_rate_rse_limit(
            thresholds, process=process, process_type=process_type
        )
        if required_limit is not None and not _required_rate_is_relevant(
            stats,
            process_peak=rate_mean_peaks.get(key),
            minimum_fraction=(
                thresholds.required_rate_min_process_peak_fraction
            ),
        ):
            required_limit = None
        is_major = (
            max_rate > 0.0
            and stats.mean is not None
            and stats.mean >= cutoff
        )
        if is_major and stats.relative_standard_error is not None:
            major_rses.append(stats.relative_standard_error)
        if required_limit is None and not is_major:
            continue

        gate_name = "required_rate" if required_limit is not None else "major_rate"
        rse_limit = (
            required_limit
            if required_limit is not None
            else thresholds.major_rate_rse
        )
        label = f"{species}:{process}:{process_type}"
        if not require_uncertainty:
            continue
        if stats.estimate_status == "censored_all_zero":
            failures.append(f"{gate_name}_censored_all_zero:{label}")
            if required_limit is not None:
                required_uncertainty_available = False
            continue
        if (
            not stats.uncertainty_available
            or stats.relative_standard_error is None
        ):
            failures.append(f"{gate_name}_rse_unavailable:{label}")
            if required_limit is not None:
                required_uncertainty_available = False
        elif stats.relative_standard_error > rse_limit:
            failures.append(f"{gate_name}_rse_exceeds_threshold:{label}")

    if eedf_normalization_error is None:
        failures.append("eedf_normalization_unavailable")
    elif eedf_normalization_error > thresholds.eedf_normalization_error:
        failures.append("eedf_normalization_error_exceeds_threshold")
    if not eedf_finite_nonnegative:
        finite_nonnegative = False
        failures.append("eedf_not_finite_nonnegative")

    uncertainty_available = all(
        stats.uncertainty_available
        for stats in [mobility, diffusion_l, diffusion_t]
    )
    reported_energy_stats = [
        stats for stats, _limit, _label in energy_stats if stats.mean is not None
    ]
    if reported_energy_stats:
        uncertainty_available = uncertainty_available and all(
            stats.uncertainty_available for stats in reported_energy_stats
        )
    if max_rate > 0.0:
        uncertainty_available = uncertainty_available and not any(
            reason.startswith("major_rate_rse_unavailable")
            or reason.startswith("major_rate_censored")
            for reason in failures
        )
    uncertainty_available = (
        uncertainty_available and required_uncertainty_available
    )

    details = {
        "mobility_rse": mobility.relative_standard_error,
        "diffusion_L_rse": diffusion_l.relative_standard_error,
        "diffusion_T_rse": diffusion_t.relative_standard_error,
        "energy_mobility_rse": energy_mobility.relative_standard_error,
        "energy_diffusion_L_rse": energy_diffusion_l.relative_standard_error,
        "energy_diffusion_T_rse": energy_diffusion_t.relative_standard_error,
        "max_major_rate_rse": max(major_rses) if major_rses else None,
        "eedf_normalization_error": eedf_normalization_error,
    }
    return (
        len(failures) == 0 and finite_nonnegative,
        tuple(failures),
        details,
        uncertainty_available,
    )


def _required_rate_is_relevant(
    stats: _ScalarStats,
    *,
    process_peak: float | None,
    minimum_fraction: float,
) -> bool:
    if minimum_fraction == 0.0:
        return True
    if (
        stats.mean is None
        or not math.isfinite(stats.mean)
        or stats.mean < 0.0
        or process_peak is None
        or not math.isfinite(process_peak)
        or process_peak <= 0.0
    ):
        return True
    return stats.mean >= process_peak * minimum_fraction


def _required_rate_rse_limit(
    thresholds: _quality_policy.QualityThresholds,
    *,
    process: str,
    process_type: str,
) -> float | None:
    """Return the mandatory RSE limit, with process overriding process type."""

    selectors = thresholds.required_rate_rse
    for configured_name, configured_limit in selectors.process.items():
        if process.casefold() == configured_name.casefold():
            return (
                thresholds.major_rate_rse
                if configured_limit is None
                else float(configured_limit)
            )
    for configured_name, configured_limit in selectors.process_type.items():
        if process_type.casefold() == configured_name.casefold():
            return (
                thresholds.major_rate_rse
                if configured_limit is None
                else float(configured_limit)
            )
    return None


def _scalar_meets_quality(
    stats: _ScalarStats,
    *,
    rse_limit: float | None,
    require_uncertainty: bool,
) -> bool:
    if stats.mean is None or not math.isfinite(stats.mean) or stats.mean < 0.0:
        return False
    if not require_uncertainty:
        return True
    return (
        stats.uncertainty_available
        and stats.relative_standard_error is not None
        and (rse_limit is None or stats.relative_standard_error <= rse_limit)
    )


def _insert_scalar(
    connection: sqlite3.Connection,
    *,
    solver: str,
    mixture_id: int,
    e_over_n_Td: float,
    scalar_group: str,
    scalar_name: str,
    stats: _ScalarStats,
    species: str = "",
    process: str = "",
    process_type: str = "",
    threshold_key: str = "",
    threshold_eV: float | None = None,
    target_species_fraction: float | None = None,
    energy_loss_eV: float | None = None,
) -> None:
    connection.execute(
        """
        INSERT INTO aggregate_scalars(
            mixture_id, solver, e_over_n_Td, scalar_group, scalar_name,
            species, process, process_type, threshold_key, threshold_eV,
            target_species_fraction, energy_loss_eV, mean, sample_stddev,
            standard_error, relative_standard_error, ci95_low, ci95_high,
            ci95_critical_value, estimate_status, valid_replicates, n,
            uncertainty_available
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            mixture_id,
            solver,
            e_over_n_Td,
            scalar_group,
            scalar_name,
            species,
            process,
            process_type,
            threshold_key,
            threshold_eV,
            _finite_or_none(target_species_fraction),
            _finite_or_none(energy_loss_eV),
            _finite_or_none(stats.mean),
            _finite_or_none(stats.sample_stddev),
            _finite_or_none(stats.standard_error),
            _finite_or_none(stats.relative_standard_error),
            _finite_or_none(stats.ci95_low),
            _finite_or_none(stats.ci95_high),
            _finite_or_none(stats.ci95_critical_value),
            stats.estimate_status,
            stats.valid_replicates,
            stats.valid_replicates,
            int(stats.uncertainty_available),
        ),
    )


def _threshold_key(value: object) -> tuple[str, float | None]:
    if value is None:
        return "", None
    number = float(value)
    return f"{number:.17g}", number


def _rate_key(row: sqlite3.Row) -> tuple[str, str, str, str, float | None]:
    key, threshold = _threshold_key(row["threshold_eV"])
    return (
        str(row["species"]),
        str(row["process"]),
        str(row["process_type"]),
        key,
        threshold,
    )


def _constant_optional_number(values: Iterable[float | None]) -> float | None:
    finite = [_finite_or_none(value) for value in values]
    finite = [value for value in finite if value is not None]
    if not finite:
        return None
    first = finite[0]
    if any(not math.isclose(first, value, rel_tol=1.0e-12, abs_tol=1.0e-30) for value in finite[1:]):
        raise ValueError("rate static metadata differs across replicates")
    return first


def _order_columns(table_name: str) -> str:
    if table_name == "aggregate_scalars":
        return (
            "solver, mixture_id, e_over_n_Td, scalar_group, scalar_name, "
            "species, process, process_type, threshold_key"
        )
    if table_name == "aggregate_eedf_bins":
        return "solver, mixture_id, e_over_n_Td, bin_index"
    return "solver, mixture_id, e_over_n_Td"


def _write_query_csv(
    connection: sqlite3.Connection,
    path: Path,
    query: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = connection.execute(query).fetchall()
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    columns = rows[0].keys()
    with path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=list(columns))
        writer.writeheader()
        for row in rows:
            writer.writerow({column: row[column] for column in columns})


def _finite_or_none(value: object) -> float | None:
    if value is None:
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def _json_dumps(value: object) -> str:
    return json.dumps(value, sort_keys=True)
