"""Aggregate workflow SQLite results for table-building stages."""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
import math
from pathlib import Path
import sqlite3
from typing import Iterable

import numpy as np

from ._io import write_json as _write_json
from .aggregate_stats import (
    ScalarStats,
    conservative_rebin_probability_mass,
    energy_edges_from_centers,
    summarize_replicates,
)
from .store import validate_workflow_schema


CASE_SCALARS = (
    "mean_energy_eV",
    "drift_velocity_m_s",
    "reduced_mobility_m2_V_s_m3",
    "reduced_diffusion_L_m2_s_m3",
    "reduced_diffusion_T_m2_s_m3",
    "reduced_electron_energy_mobility_m2_V_s_m3",
    "reduced_electron_energy_diffusion_m2_s_m3",
    "effective_townsend_m2",
)
RATE_SCALAR = "rate_coefficient_m3_s"
MC_SOLVER = "monte_carlo"

AGGREGATE_CSVS = (
    "aggregate_scalars.csv",
    "aggregate_eedf_bins.csv",
    "aggregate_quality.csv",
)

__all__ = [
    "AGGREGATE_CSVS",
    "CASE_SCALARS",
    "QualityThresholds",
    "AggregateSummary",
    "ScalarStats",
    "aggregate_database",
    "aggregate_monte_carlo_results",
    "aggregate_workflow_results",
    "conservative_rebin_probability_mass",
    "parse_quality_thresholds",
    "stable_mc_seed",
    "summarize_replicates",
    "write_aggregate_outputs",
]


@dataclass(frozen=True, slots=True)
class QualityThresholds:
    mobility_rse: float = 0.02
    diffusion_rse: float = 0.05
    major_rate_rse: float = 0.05
    major_rate_fraction: float = 0.01
    eedf_normalization_error: float = 1.0e-6


@dataclass(frozen=True, slots=True)
class AggregateSummary:
    database_path: Path
    output_directory: Path | None
    solvers: tuple[str, ...]
    mixtures: int
    e_over_n_points: int


def parse_quality_thresholds(raw: object) -> QualityThresholds:
    if raw is None:
        return QualityThresholds()
    if not isinstance(raw, dict):
        raise ValueError("workflow quality must be a mapping")
    allowed = {
        "mobility_rse",
        "diffusion_rse",
        "major_rate_rse",
        "major_rate_fraction",
        "eedf_normalization_error",
    }
    unknown = set(raw) - allowed
    if unknown:
        raise ValueError(f"Unsupported workflow quality fields: {sorted(unknown)}")
    values = asdict(QualityThresholds())
    for key, value in raw.items():
        number = float(value)
        if not math.isfinite(number) or number < 0.0:
            raise ValueError(f"workflow quality.{key} must be finite and nonnegative")
        if key == "major_rate_fraction" and number > 1.0:
            raise ValueError("workflow quality.major_rate_fraction must be in [0, 1]")
        values[key] = number
    return QualityThresholds(**values)


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
    thresholds: QualityThresholds | None = None,
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


def aggregate_monte_carlo_results(
    connection: sqlite3.Connection,
    thresholds: QualityThresholds | None = None,
) -> None:
    aggregate_workflow_results(connection, thresholds, solvers=(MC_SOLVER,))


def aggregate_workflow_results(
    connection: sqlite3.Connection,
    thresholds: QualityThresholds | None = None,
    *,
    solvers: tuple[str, ...] | None = None,
) -> None:
    thresholds = thresholds or QualityThresholds()
    connection.row_factory = sqlite3.Row
    validate_workflow_schema(connection)
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
                eedf_normalization_error=eedf_normalization_error,
                eedf_finite_nonnegative=eedf_finite_nonnegative,
                thresholds=thresholds,
                require_uncertainty=require_uncertainty,
            )
            connection.execute(
                """
                INSERT INTO aggregate_quality(
                    mixture_id, solver, e_over_n_Td, passed,
                    failure_reasons_json, thresholds_json, mobility_rse,
                    diffusion_L_rse, diffusion_T_rse, max_major_rate_rse,
                    eedf_normalization_error, valid_replicates, n,
                    uncertainty_available
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    mixture_id,
                    solver,
                    e_over_n_Td,
                    int(passed),
                    _json_dumps(list(failure_reasons)),
                    _json_dumps(asdict(thresholds)),
                    _finite_or_none(details["mobility_rse"]),
                    _finite_or_none(details["diffusion_L_rse"]),
                    _finite_or_none(details["diffusion_T_rse"]),
                    _finite_or_none(details["max_major_rate_rse"]),
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
        _write_query_csv(
            connection,
            output / filename,
            f"SELECT * FROM {table_name} ORDER BY { _order_columns(table_name) }",
        )
    metadata = _metadata(connection)
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
            "hashes": _hash_manifest(metadata),
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
            thresholds_json TEXT NOT NULL,
            mobility_rse REAL,
            diffusion_L_rse REAL,
            diffusion_T_rse REAL,
            max_major_rate_rse REAL,
            eedf_normalization_error REAL,
            valid_replicates INTEGER NOT NULL,
            n INTEGER NOT NULL,
            uncertainty_available INTEGER NOT NULL,
            PRIMARY KEY (mixture_id, solver, e_over_n_Td)
        );
        """
    )


def _case_scalar_stats(
    connection: sqlite3.Connection,
    solver: str,
    mixture_id: int,
    e_over_n_Td: float,
) -> dict[str, ScalarStats]:
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
        scalar_name: summarize_replicates(row[scalar_name] for row in rows)
        for scalar_name in CASE_SCALARS
    }


def _rate_scalar_stats(
    connection: sqlite3.Connection,
    solver: str,
    mixture_id: int,
    e_over_n_Td: float,
) -> tuple[
    dict[tuple[str, str, str, str, float | None], ScalarStats],
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
    return {key: summarize_replicates(values) for key, values in grouped.items()}, metadata


def _aggregate_eedf(
    connection: sqlite3.Connection,
    solver: str,
    mixture_id: int,
    e_over_n_Td: float,
) -> tuple[float | None, bool]:
    rows = connection.execute(
        """
        SELECT replicate, bin_index, energy_eV, energy_width_eV, eedf
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
    endpoints: list[float] = []
    for replicate, replicate_rows in by_replicate.items():
        energy = np.asarray([row["energy_eV"] for row in replicate_rows], dtype=float)
        widths = np.asarray([row["energy_width_eV"] for row in replicate_rows], dtype=float)
        eedf = np.asarray([row["eedf"] for row in replicate_rows], dtype=float)
        if (
            len(energy) == 0
            or np.any(~np.isfinite(energy))
            or np.any(~np.isfinite(widths))
            or np.any(widths <= 0.0)
            or np.any(~np.isfinite(eedf))
        ):
            continue
        edges = energy_edges_from_centers(energy)
        mass = eedf * widths
        replicate_edges[replicate] = edges
        replicate_masses[replicate] = mass
        endpoints.extend(float(value) for value in edges)
    if not replicate_edges:
        return None, False

    target_edges = np.unique(np.asarray(endpoints, dtype=float))
    target_edges = target_edges[np.isfinite(target_edges)]
    target_edges.sort()
    target_edges = target_edges[np.concatenate(([True], np.diff(target_edges) > 0.0))]
    if len(target_edges) < 2:
        return None, False

    rebinned = np.vstack(
        [
            conservative_rebin_probability_mass(
                replicate_edges[replicate],
                replicate_masses[replicate],
                target_edges,
            )
            for replicate in sorted(replicate_edges)
        ]
    )
    stats = [summarize_replicates(rebinned[:, index]) for index in range(rebinned.shape[1])]
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
            eedf, valid_replicates, n, uncertainty_available
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
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
                stats[index].valid_replicates,
                stats[index].valid_replicates,
                int(stats[index].uncertainty_available),
            )
            for index in range(len(widths))
            if widths[index] > 0.0
        ],
    )
    return normalization_error, finite_nonnegative


def _quality_for_group(
    *,
    case_stats: dict[str, ScalarStats],
    rate_stats: dict[tuple[str, str, str, str, float | None], ScalarStats],
    eedf_normalization_error: float | None,
    eedf_finite_nonnegative: bool,
    thresholds: QualityThresholds,
    require_uncertainty: bool,
) -> tuple[bool, tuple[str, ...], dict[str, float | None], bool]:
    failures: list[str] = []
    finite_nonnegative = True

    for name, stats in case_stats.items():
        if stats.mean is None or not math.isfinite(stats.mean) or stats.mean < 0.0:
            finite_nonnegative = False
            failures.append(f"{name}_not_finite_nonnegative")

    mobility = case_stats["reduced_mobility_m2_V_s_m3"]
    diffusion_l = case_stats["reduced_diffusion_L_m2_s_m3"]
    diffusion_t = case_stats["reduced_diffusion_T_m2_s_m3"]
    if not _scalar_meets_quality(
        mobility,
        rse_limit=thresholds.mobility_rse,
        require_uncertainty=require_uncertainty,
    ):
        failures.append("reduced_mobility_rse_unavailable_or_exceeds_threshold")
    if not _scalar_meets_quality(
        diffusion_l,
        rse_limit=thresholds.diffusion_rse,
        require_uncertainty=require_uncertainty,
    ):
        failures.append("reduced_diffusion_L_rse_unavailable_or_exceeds_threshold")
    if not _scalar_meets_quality(
        diffusion_t,
        rse_limit=thresholds.diffusion_rse,
        require_uncertainty=require_uncertainty,
    ):
        failures.append("reduced_diffusion_T_rse_unavailable_or_exceeds_threshold")

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
    major_rses: list[float] = []
    if max_rate > 0.0:
        cutoff = max_rate * thresholds.major_rate_fraction
        for key, stats in rate_stats.items():
            species, process, process_type, _threshold_key_value, _threshold = key
            if stats.mean is None or stats.mean < cutoff:
                continue
            if stats.relative_standard_error is not None:
                major_rses.append(stats.relative_standard_error)
            if require_uncertainty and (
                not stats.uncertainty_available
                or stats.relative_standard_error is None
            ):
                failures.append(
                    f"major_rate_rse_unavailable:{species}:{process}:{process_type}"
                )
            elif (
                stats.relative_standard_error is not None
                and stats.relative_standard_error > thresholds.major_rate_rse
            ):
                failures.append(
                    f"major_rate_rse_exceeds_threshold:{species}:{process}:{process_type}"
                )

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
    if max_rate > 0.0:
        uncertainty_available = uncertainty_available and not any(
            reason.startswith("major_rate_rse_unavailable") for reason in failures
        )

    details = {
        "mobility_rse": mobility.relative_standard_error,
        "diffusion_L_rse": diffusion_l.relative_standard_error,
        "diffusion_T_rse": diffusion_t.relative_standard_error,
        "max_major_rate_rse": max(major_rses) if major_rses else None,
        "eedf_normalization_error": eedf_normalization_error,
    }
    return (
        len(failures) == 0 and finite_nonnegative,
        tuple(failures),
        details,
        uncertainty_available,
    )


def _scalar_meets_quality(
    stats: ScalarStats,
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
    stats: ScalarStats,
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
            valid_replicates, n, uncertainty_available
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
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


def _metadata(connection: sqlite3.Connection) -> dict[str, str]:
    if not _table_exists(connection, "metadata"):
        return {}
    return {
        str(row["key"]): str(row["value"])
        for row in connection.execute("SELECT key, value FROM metadata")
    }


def _hash_manifest(metadata: dict[str, str]) -> dict[str, str | None]:
    return {
        "base_config_sha256": metadata.get("base_config_sha256"),
        "cross_sections_sha256": metadata.get("cross_sections_sha256"),
    }


def _table_exists(connection: sqlite3.Connection, name: str) -> bool:
    row = connection.execute(
        "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ?",
        (name,),
    ).fetchone()
    return row is not None


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
