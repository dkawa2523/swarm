"""Build COMSOL-ready workflow table directories from aggregate data."""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path
import sqlite3
from typing import Any

from ._io import write_csv as _write_csv_rows
from ._io import write_json as _write_json
from .aggregate import (
    MC_SOLVER,
    RATE_SCALAR,
    QualityThresholds,
    aggregate_workflow_results,
)
from .table_math import (
    finite_range,
    interpolate_no_extrapolate,
    positive_log_interpolate,
    strictly_monotonic,
)


FORMAT_VERSION = 1
TD_TO_V_M2 = 1.0e-21
TWO_TERM_SOLVER = "two_term"

SOURCE_CHOICES = ("two_term", "monte_carlo", "hybrid")

CASE_COLUMNS = (
    "E_over_N_Td",
    "E_over_N_V_m2",
    "mean_energy_eV",
    "drift_velocity_m_s",
    "reduced_mobility_m2_V_s_m3",
    "reduced_diffusion_L_m2_s_m3",
    "reduced_diffusion_T_m2_s_m3",
    "reduced_electron_energy_mobility_m2_V_s_m3",
    "reduced_electron_energy_diffusion_m2_s_m3",
)
RATE_COLUMNS = (
    "E_over_N_Td",
    "E_over_N_V_m2",
    "species",
    "process",
    "process_type",
    "threshold_eV",
    "target_species_fraction",
    "rate_coefficient_m3_s",
    "mixture_weighted_rate_m3_s",
    "reduced_townsend_m2",
    "mixture_weighted_reduced_townsend_m2",
)
ENERGY_LOSS_COLUMNS = (
    "E_over_N_Td",
    "E_over_N_V_m2",
    "species",
    "process",
    "process_type",
    "threshold_eV",
    "target_species_fraction",
    "energy_loss_eV",
    "energy_loss_rate_coefficient_eV_m3_s",
)
QUALITY_COLUMNS = (
    "E_over_N_Td",
    "E_over_N_V_m2",
    "passed",
    "failure_reasons_json",
    "mobility_rse",
    "diffusion_L_rse",
    "diffusion_T_rse",
    "max_major_rate_rse",
    "eedf_normalization_error",
    "valid_replicates",
    "uncertainty_available",
    "solver_diagnostics_available",
    "solver_converged",
    "solver_iterations_max",
    "solver_residual_L1_max",
    "solver_residual_tolerance",
    "solver_tail_probability_max",
    "solver_tail_probability_target",
    "solver_edge_to_peak_max",
    "solver_edge_to_peak_target",
    "solver_grid_max_eV_max",
    "solver_grid_max_limit_eV",
    "solver_grid_limit_hit",
    "solver_diagnostics_passed",
    "quality_source",
)
EEDF_COLUMNS = (
    "electron_energy_eV",
    "energy_width_eV",
    "E_over_N_Td",
    "E_over_N_V_m2",
    "mean_energy_eV",
    "eedf",
)

BUILT_TABLES = (
    "mean_energy_vs_en.csv",
    "transport_vs_en.csv",
    "rates_vs_en.csv",
    "townsend_vs_en.csv",
    "energy_loss.csv",
    "quality.csv",
    "eedf.csv",
)

UNITS = {
    "E_over_N_Td": "Td",
    "E_over_N_V_m2": "V m^2",
    "mean_energy_eV": "eV",
    "electron_energy_eV": "eV",
    "energy_width_eV": "eV",
    "drift_velocity_m_s": "m/s",
    "reduced_mobility_m2_V_s_m3": "1/(V m s)",
    "reduced_diffusion_L_m2_s_m3": "1/(m s)",
    "reduced_diffusion_T_m2_s_m3": "1/(m s)",
    "reduced_electron_energy_mobility_m2_V_s_m3": "1/(V m s)",
    "reduced_electron_energy_diffusion_m2_s_m3": "1/(m s)",
    "effective_townsend_m2": "m^2",
    "rate_coefficient_m3_s": "m^3/s",
    "mixture_weighted_rate_m3_s": "m^3/s",
    "reduced_townsend_m2": "m^2",
    "mixture_weighted_reduced_townsend_m2": "m^2",
    "threshold_eV": "eV",
    "energy_loss_eV": "eV",
    "energy_loss_rate_coefficient_eV_m3_s": "eV m^3/s",
    "eedf": "eV^-1",
}


class TableBuildError(RuntimeError):
    """Raised when workflow tables cannot be built safely."""


@dataclass(frozen=True, slots=True)
class TableBuildSummary:
    database_path: Path
    output_directory: Path
    source: str
    mixtures: int


@dataclass(slots=True)
class _Dataset:
    source: str
    mixture_id: int
    mixture: list[dict[str, str | float]]
    cases: list[dict[str, Any]]
    rates: list[dict[str, Any]]
    eedf: list[dict[str, Any]]
    quality: list[dict[str, Any]]
    source_policy: dict[str, Any]


def build_tables(
    database_path: str | Path,
    output_directory: str | Path,
    *,
    source: str,
    thresholds: QualityThresholds | None = None,
) -> TableBuildSummary:
    if source not in SOURCE_CHOICES:
        raise TableBuildError(f"unsupported table source: {source}")
    db_path = Path(database_path)
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)

    connection = sqlite3.connect(db_path)
    try:
        connection.row_factory = sqlite3.Row
        aggregate_workflow_results(connection, thresholds)
        metadata = _metadata(connection)
        mixture_ids = _mixture_ids(connection, source)
        if not mixture_ids:
            raise TableBuildError(f"database has no usable {source} data")
        entries: list[dict[str, Any]] = []
        for mixture_id in mixture_ids:
            dataset = _build_dataset(connection, mixture_id, source)
            mixture_dir = output / f"mixture_{mixture_id:04d}"
            _write_dataset_tables(dataset, mixture_dir, metadata=metadata)
            entries.append(
                {
                    "mixture_id": mixture_id,
                    "path": f"mixture_{mixture_id:04d}/manifest.json",
                    "status": "ok",
                }
            )
    finally:
        connection.close()

    _write_json(
        output / "manifest.json",
        {
            "format_version": FORMAT_VERSION,
            "stage": "build-tables",
            "source": source,
            "database_path": str(db_path),
            "hashes": _hash_manifest(metadata),
            "mixtures": entries,
        },
    )
    return TableBuildSummary(db_path, output, source, len(entries))


def _build_dataset(
    connection: sqlite3.Connection,
    mixture_id: int,
    source: str,
) -> _Dataset:
    mixture = _mixture_rows(connection, mixture_id)
    if source == "two_term":
        cases = _load_cases(connection, TWO_TERM_SOLVER, mixture_id)
        rates = _load_rates(connection, TWO_TERM_SOLVER, mixture_id, cases)
        eedf = _load_eedf(connection, TWO_TERM_SOLVER, mixture_id, cases)
        quality = _load_quality(connection, TWO_TERM_SOLVER, mixture_id)
        source_policy = {"source": "two_term", "selection": "all_aggregate_points"}
    elif source == "monte_carlo":
        allowed = _passing_e_over_n(connection, MC_SOLVER, mixture_id)
        if not allowed:
            raise TableBuildError(
                f"mixture {mixture_id} has no quality-passing monte_carlo points"
            )
        cases = _load_cases(connection, MC_SOLVER, mixture_id, allowed_e=allowed)
        rates = _load_rates(connection, MC_SOLVER, mixture_id, cases, allowed_e=allowed)
        eedf = _load_eedf(connection, MC_SOLVER, mixture_id, cases, allowed_e=allowed)
        quality = _load_quality(connection, MC_SOLVER, mixture_id, allowed_e=allowed)
        source_policy = {
            "source": "monte_carlo",
            "selection": "quality_passing_aggregate_points_only",
        }
    else:
        cases, rates, eedf, quality, source_policy = _build_hybrid(
            connection, mixture_id
        )
    if not cases:
        raise TableBuildError(f"mixture {mixture_id} has no cases for {source}")
    return _Dataset(source, mixture_id, mixture, cases, rates, eedf, quality, source_policy)


def _build_hybrid(
    connection: sqlite3.Connection,
    mixture_id: int,
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[str, Any],
]:
    base_cases = _load_cases(connection, TWO_TERM_SOLVER, mixture_id)
    base_rates = _load_rates(connection, TWO_TERM_SOLVER, mixture_id, base_cases)
    base_eedf = _load_eedf(connection, TWO_TERM_SOLVER, mixture_id, base_cases)
    if not base_cases:
        raise TableBuildError("hybrid source requires two_term baseline")
    mc_allowed = _passing_e_over_n(connection, MC_SOLVER, mixture_id)
    mc_cases = (
        _load_cases(connection, MC_SOLVER, mixture_id, allowed_e=mc_allowed)
        if mc_allowed
        else []
    )
    mc_rates = (
        _load_rates(connection, MC_SOLVER, mixture_id, mc_cases, allowed_e=mc_allowed)
        if mc_allowed
        else []
    )
    mc_quality = _load_quality(connection, MC_SOLVER, mixture_id, allowed_e=mc_allowed)
    thresholds = _thresholds_from_quality(mc_quality)

    cases = [dict(row) for row in base_cases]
    for scalar in [
        "mean_energy_eV",
        "drift_velocity_m_s",
        "reduced_mobility_m2_V_s_m3",
        "reduced_diffusion_L_m2_s_m3",
        "reduced_diffusion_T_m2_s_m3",
        "reduced_electron_energy_mobility_m2_V_s_m3",
        "reduced_electron_energy_diffusion_m2_s_m3",
    ]:
        _apply_case_log_correction(
            cases,
            base_cases,
            mc_cases,
            scalar,
            thresholds=thresholds,
        )

    rates = _hybrid_rates(base_rates, mc_rates, thresholds)
    _derive_net_townsend_from_rates(cases, rates)
    quality = _hybrid_quality_rows(base_cases, mc_quality)
    source_policy = {
        "source": "hybrid",
        "baseline": "two_term",
        "correction_source": "monte_carlo_quality_passing_points",
        "positive_quantity_policy": "log_delta_with_rse_shrinkage",
        "net_ionization_policy": "derive_from_ionization_and_attachment_rates",
        "eedf_policy": "two_term_baseline_not_mc_corrected_phase_9_1",
        "mc_points_used": len(mc_allowed),
    }
    return cases, rates, base_eedf, quality, source_policy


def _apply_case_log_correction(
    cases: list[dict[str, Any]],
    base_cases: list[dict[str, Any]],
    mc_cases: list[dict[str, Any]],
    scalar: str,
    *,
    thresholds: QualityThresholds,
) -> None:
    known_x: list[float] = []
    known_delta: list[float] = []
    for mc_case in mc_cases:
        x = float(mc_case["E_over_N_Td"])
        mc_value = _finite_positive(mc_case.get(scalar))
        base_value = positive_log_interpolate(base_cases, scalar, x)
        if mc_value is None or base_value is None:
            continue
        weight = _case_weight(mc_case, scalar, thresholds)
        known_x.append(x)
        known_delta.append(math.log(mc_value / base_value) * weight)
    for case in cases:
        base_value = _finite_positive(case.get(scalar))
        if base_value is None:
            continue
        delta = interpolate_no_extrapolate(
            float(case["E_over_N_Td"]), known_x, known_delta
        )
        case[scalar] = base_value * math.exp(delta)


def _hybrid_rates(
    base_rates: list[dict[str, Any]],
    mc_rates: list[dict[str, Any]],
    thresholds: QualityThresholds,
) -> list[dict[str, Any]]:
    by_key: dict[tuple[str, str, str, str], list[dict[str, Any]]] = {}
    for rate in mc_rates:
        by_key.setdefault(_rate_key(rate), []).append(rate)
    base_by_key: dict[tuple[str, str, str, str], list[dict[str, Any]]] = {}
    for rate in base_rates:
        base_by_key.setdefault(_rate_key(rate), []).append(rate)

    corrected: list[dict[str, Any]] = []
    for rate in base_rates:
        row = dict(rate)
        key = _rate_key(rate)
        known_x: list[float] = []
        known_delta: list[float] = []
        for mc_rate in by_key.get(key, []):
            x = float(mc_rate["E_over_N_Td"])
            mc_value = _finite_positive(mc_rate.get("rate_coefficient_m3_s"))
            base_value = positive_log_interpolate(
                base_by_key.get(key, []), "rate_coefficient_m3_s", x
            )
            if mc_value is None or base_value is None:
                continue
            rse = _finite_nonnegative(mc_rate.get("relative_standard_error"))
            weight = _rse_weight(rse, thresholds.major_rate_rse)
            known_x.append(x)
            known_delta.append(math.log(mc_value / base_value) * weight)
        base_value = _finite_positive(row["rate_coefficient_m3_s"])
        if base_value is not None:
            delta = interpolate_no_extrapolate(
                float(row["E_over_N_Td"]), known_x, known_delta
            )
            row["rate_coefficient_m3_s"] = base_value * math.exp(delta)
        _fill_rate_derived_columns(row)
        corrected.append(row)
    return corrected


def _derive_net_townsend_from_rates(
    cases: list[dict[str, Any]],
    rates: list[dict[str, Any]],
) -> None:
    by_e: dict[float, list[dict[str, Any]]] = {}
    for rate in rates:
        by_e.setdefault(float(rate["E_over_N_Td"]), []).append(rate)
    for case in cases:
        rows = by_e.get(float(case["E_over_N_Td"]), [])
        ionization = 0.0
        attachment = 0.0
        saw_component = False
        for row in rows:
            process_type = str(row["process_type"]).lower()
            value = float(row["mixture_weighted_reduced_townsend_m2"])
            if "ionization" in process_type:
                ionization += value
                saw_component = True
            elif "attachment" in process_type:
                attachment += value
                saw_component = True
        if saw_component:
            case["effective_townsend_m2"] = ionization - attachment


def _write_dataset_tables(
    dataset: _Dataset,
    directory: Path,
    *,
    metadata: dict[str, str],
) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    _validate_positive_values(dataset)
    monotonic = strictly_monotonic(
        [float(case["mean_energy_eV"]) for case in dataset.cases]
    )
    monotonic_reason = None if monotonic else "mean_energy_not_strictly_monotonic"

    tables: dict[str, dict[str, Any]] = {}
    tables["mean_energy_vs_en.csv"] = _write_csv(
        directory / "mean_energy_vs_en.csv",
        ("E_over_N_Td", "E_over_N_V_m2", "mean_energy_eV"),
        dataset.cases,
        argument="E_over_N_Td",
    )
    tables["transport_vs_en.csv"] = _write_csv(
        directory / "transport_vs_en.csv",
        CASE_COLUMNS,
        dataset.cases,
        argument="E_over_N_Td",
    )
    tables["rates_vs_en.csv"] = _write_csv(
        directory / "rates_vs_en.csv",
        RATE_COLUMNS,
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
        ENERGY_LOSS_COLUMNS,
        dataset.rates,
        argument="E_over_N_Td",
    )
    tables["quality.csv"] = _write_csv(
        directory / "quality.csv",
        QUALITY_COLUMNS,
        dataset.quality,
        argument="E_over_N_Td",
    )
    tables["eedf.csv"] = _write_csv(
        directory / "eedf.csv",
        EEDF_COLUMNS,
        dataset.eedf,
        argument="electron_energy_eV,E_over_N_Td",
    )
    if monotonic:
        tables["transport_vs_mean_energy.csv"] = _write_csv(
            directory / "transport_vs_mean_energy.csv",
            ("mean_energy_eV", "E_over_N_Td", "E_over_N_V_m2", *CASE_COLUMNS[3:]),
            dataset.cases,
            argument="mean_energy_eV",
        )
        tables["rates_vs_mean_energy.csv"] = _write_csv(
            directory / "rates_vs_mean_energy.csv",
            ("mean_energy_eV", *RATE_COLUMNS),
            dataset.rates,
            argument="mean_energy_eV",
        )
    _write_json(
        directory / "manifest.json",
        {
            "format_version": FORMAT_VERSION,
            "stage": "build-tables",
            "source": dataset.source,
            "hashes": _hash_manifest(metadata),
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


def _mixture_ids(connection: sqlite3.Connection, source: str) -> list[int]:
    if source == "two_term":
        solvers = (TWO_TERM_SOLVER,)
    elif source == "monte_carlo":
        solvers = (MC_SOLVER,)
    else:
        solvers = (TWO_TERM_SOLVER,)
    placeholders = ", ".join("?" for _solver in solvers)
    rows = connection.execute(
        f"""
        SELECT DISTINCT mixture_id
        FROM aggregate_scalars
        WHERE solver IN ({placeholders})
        ORDER BY mixture_id
        """,
        solvers,
    ).fetchall()
    return [int(row["mixture_id"]) for row in rows]


def _mixture_rows(
    connection: sqlite3.Connection,
    mixture_id: int,
) -> list[dict[str, str | float]]:
    rows = connection.execute(
        """
        SELECT species, fraction, mass_amu
        FROM mixture_species
        WHERE mixture_id = ?
        ORDER BY species
        """,
        (mixture_id,),
    ).fetchall()
    return [
        {
            "species": str(row["species"]),
            "fraction": float(row["fraction"]),
            "mass_amu": float(row["mass_amu"]),
        }
        for row in rows
    ]


def _load_cases(
    connection: sqlite3.Connection,
    solver: str,
    mixture_id: int,
    *,
    allowed_e: set[float] | None = None,
) -> list[dict[str, Any]]:
    rows = connection.execute(
        """
        SELECT e_over_n_Td, scalar_name, mean, relative_standard_error
        FROM aggregate_scalars
        WHERE solver = ? AND mixture_id = ? AND scalar_group = 'case'
        ORDER BY e_over_n_Td, scalar_name
        """,
        (solver, mixture_id),
    ).fetchall()
    by_e: dict[float, dict[str, Any]] = {}
    for row in rows:
        e_over_n = float(row["e_over_n_Td"])
        if allowed_e is not None and e_over_n not in allowed_e:
            continue
        case = by_e.setdefault(e_over_n, _en_fields(e_over_n))
        case[str(row["scalar_name"])] = _float_or_none(row["mean"])
        case[f"{row['scalar_name']}_rse"] = _float_or_none(
            row["relative_standard_error"]
        )
    return [by_e[e_over_n] for e_over_n in sorted(by_e)]


def _load_rates(
    connection: sqlite3.Connection,
    solver: str,
    mixture_id: int,
    cases: list[dict[str, Any]],
    *,
    allowed_e: set[float] | None = None,
) -> list[dict[str, Any]]:
    by_e = {float(case["E_over_N_Td"]): case for case in cases}
    rows = connection.execute(
        """
        SELECT e_over_n_Td, species, process, process_type, threshold_eV,
               target_species_fraction, energy_loss_eV, mean,
               relative_standard_error
        FROM aggregate_scalars
        WHERE solver = ? AND mixture_id = ? AND scalar_group = 'rate'
          AND scalar_name = ?
        ORDER BY e_over_n_Td, species, process, process_type, threshold_key
        """,
        (solver, mixture_id, RATE_SCALAR),
    ).fetchall()
    rates: list[dict[str, Any]] = []
    for row in rows:
        e_over_n = float(row["e_over_n_Td"])
        if allowed_e is not None and e_over_n not in allowed_e:
            continue
        case = by_e.get(e_over_n)
        if case is None:
            continue
        rate = {
            **_en_fields(e_over_n),
            "mean_energy_eV": case.get("mean_energy_eV"),
            "species": str(row["species"]),
            "process": str(row["process"]),
            "process_type": str(row["process_type"]),
            "threshold_eV": _float_or_none(row["threshold_eV"]),
            "target_species_fraction": _float_or_none(
                row["target_species_fraction"]
            )
            or 1.0,
            "energy_loss_eV": _float_or_none(row["energy_loss_eV"]),
            "rate_coefficient_m3_s": _float_or_none(row["mean"]),
            "relative_standard_error": _float_or_none(row["relative_standard_error"]),
            "case_drift_velocity_m_s": case.get("drift_velocity_m_s"),
        }
        _fill_rate_derived_columns(rate, drift=case.get("drift_velocity_m_s"))
        rates.append(rate)
    return rates


def _load_eedf(
    connection: sqlite3.Connection,
    solver: str,
    mixture_id: int,
    cases: list[dict[str, Any]],
    *,
    allowed_e: set[float] | None = None,
) -> list[dict[str, Any]]:
    mean_by_e = {
        float(case["E_over_N_Td"]): _float_or_none(case.get("mean_energy_eV"))
        for case in cases
    }
    rows = connection.execute(
        """
        SELECT e_over_n_Td, energy_eV, energy_width_eV, eedf
        FROM aggregate_eedf_bins
        WHERE solver = ? AND mixture_id = ?
        ORDER BY e_over_n_Td, bin_index
        """,
        (solver, mixture_id),
    ).fetchall()
    eedf: list[dict[str, Any]] = []
    for row in rows:
        e_over_n = float(row["e_over_n_Td"])
        if allowed_e is not None and e_over_n not in allowed_e:
            continue
        if e_over_n not in mean_by_e:
            continue
        eedf.append(
            {
                **_en_fields(e_over_n),
                "mean_energy_eV": mean_by_e[e_over_n],
                "electron_energy_eV": float(row["energy_eV"]),
                "energy_width_eV": float(row["energy_width_eV"]),
                "eedf": float(row["eedf"]),
            }
        )
    return eedf


def _load_solver_quality_diagnostics(
    connection: sqlite3.Connection,
    solver: str,
    mixture_id: int,
) -> dict[float, dict[str, Any]]:
    """Summarize raw per-case solver diagnostics for table qualification."""

    rows = connection.execute(
        """
        SELECT e_over_n_Td, diagnostics_json
        FROM cases
        WHERE solver = ? AND mixture_id = ?
        ORDER BY e_over_n_Td, replicate
        """,
        (solver, mixture_id),
    ).fetchall()
    grouped: dict[float, list[dict[str, Any]]] = {}
    for row in rows:
        payload = json.loads(str(row["diagnostics_json"]))
        diagnostic = payload.get(solver)
        if isinstance(diagnostic, dict):
            grouped.setdefault(float(row["e_over_n_Td"]), []).append(diagnostic)

    summaries: dict[float, dict[str, Any]] = {}
    for e_over_n, diagnostics in grouped.items():
        converged = all(bool(item.get("converged", False)) for item in diagnostics)
        iterations = max(int(item.get("iterations", 0)) for item in diagnostics)
        residual = max(float(item.get("residual_L1", math.inf)) for item in diagnostics)
        residual_target = min(
            float(item.get("residual_tolerance", math.nan))
            for item in diagnostics
        )
        tail_probability = max(
            float(item.get("tail_probability", math.inf))
            for item in diagnostics
        )
        tail_target = min(
            float(item.get("tail_probability_target", math.nan))
            for item in diagnostics
        )
        edge_to_peak = max(
            float(item.get("edge_to_peak", math.inf)) for item in diagnostics
        )
        edge_target = min(
            float(item.get("edge_to_peak_target", math.nan))
            for item in diagnostics
        )
        grid_max = max(
            float(item.get("grid_max_eV", math.inf)) for item in diagnostics
        )
        grid_limit = min(
            float(item.get("grid_max_limit_eV", math.nan))
            for item in diagnostics
        )
        thresholds_available = all(
            math.isfinite(value)
            for value in (
                residual_target,
                tail_target,
                edge_target,
                grid_limit,
            )
        )
        grid_limit_hit = bool(
            thresholds_available and grid_max >= grid_limit * (1.0 - 1.0e-9)
        )
        passed = bool(
            thresholds_available
            and converged
            and residual <= residual_target
            and tail_probability <= tail_target
            and edge_to_peak <= edge_target
            and not grid_limit_hit
        )
        reasons: list[str] = []
        if not thresholds_available:
            reasons.append("solver_diagnostic_thresholds_missing")
        if not converged:
            reasons.append("solver_not_converged")
        if math.isfinite(residual_target) and residual > residual_target:
            reasons.append("solver_residual_above_tolerance")
        if math.isfinite(tail_target) and tail_probability > tail_target:
            reasons.append("eedf_tail_probability_above_target")
        if math.isfinite(edge_target) and edge_to_peak > edge_target:
            reasons.append("eedf_edge_to_peak_above_target")
        if grid_limit_hit:
            reasons.append("energy_grid_limit_hit")
        summaries[e_over_n] = {
            "solver_diagnostics_available": 1,
            "solver_converged": int(converged),
            "solver_iterations_max": iterations,
            "solver_residual_L1_max": residual,
            "solver_residual_tolerance": _float_or_none(residual_target),
            "solver_tail_probability_max": tail_probability,
            "solver_tail_probability_target": _float_or_none(tail_target),
            "solver_edge_to_peak_max": edge_to_peak,
            "solver_edge_to_peak_target": _float_or_none(edge_target),
            "solver_grid_max_eV_max": grid_max,
            "solver_grid_max_limit_eV": _float_or_none(grid_limit),
            "solver_grid_limit_hit": int(grid_limit_hit),
            "solver_diagnostics_passed": int(passed),
            "solver_failure_reasons": reasons,
        }
    return summaries


def _load_quality(
    connection: sqlite3.Connection,
    solver: str,
    mixture_id: int,
    *,
    allowed_e: set[float] | None = None,
) -> list[dict[str, Any]]:
    diagnostics_by_e = _load_solver_quality_diagnostics(
        connection, solver, mixture_id
    )
    rows = connection.execute(
        """
        SELECT e_over_n_Td, passed, failure_reasons_json, thresholds_json,
               mobility_rse, diffusion_L_rse, diffusion_T_rse,
               max_major_rate_rse, eedf_normalization_error,
               valid_replicates, uncertainty_available
        FROM aggregate_quality
        WHERE solver = ? AND mixture_id = ?
        ORDER BY e_over_n_Td
        """,
        (solver, mixture_id),
    ).fetchall()
    quality: list[dict[str, Any]] = []
    for row in rows:
        e_over_n = float(row["e_over_n_Td"])
        if allowed_e is not None and e_over_n not in allowed_e:
            continue
        diagnostic = diagnostics_by_e.get(e_over_n)
        try:
            failure_reasons = json.loads(str(row["failure_reasons_json"]))
        except json.JSONDecodeError:
            failure_reasons = [str(row["failure_reasons_json"])]
        if not isinstance(failure_reasons, list):
            failure_reasons = [str(failure_reasons)]
        aggregate_passed = bool(row["passed"])
        if diagnostic is not None:
            failure_reasons.extend(diagnostic["solver_failure_reasons"])
            passed = aggregate_passed and bool(
                diagnostic["solver_diagnostics_passed"]
            )
        else:
            passed = aggregate_passed
            diagnostic = {
                "solver_diagnostics_available": 0,
                "solver_converged": None,
                "solver_iterations_max": None,
                "solver_residual_L1_max": None,
                "solver_residual_tolerance": None,
                "solver_tail_probability_max": None,
                "solver_tail_probability_target": None,
                "solver_edge_to_peak_max": None,
                "solver_edge_to_peak_target": None,
                "solver_grid_max_eV_max": None,
                "solver_grid_max_limit_eV": None,
                "solver_grid_limit_hit": None,
                "solver_diagnostics_passed": None,
                "solver_failure_reasons": [],
            }
        quality.append(
            {
                **_en_fields(e_over_n),
                "passed": int(passed),
                "failure_reasons_json": json.dumps(
                    sorted(set(str(value) for value in failure_reasons)),
                    separators=(",", ":"),
                ),
                "thresholds_json": str(row["thresholds_json"]),
                "mobility_rse": _float_or_none(row["mobility_rse"]),
                "diffusion_L_rse": _float_or_none(row["diffusion_L_rse"]),
                "diffusion_T_rse": _float_or_none(row["diffusion_T_rse"]),
                "max_major_rate_rse": _float_or_none(row["max_major_rate_rse"]),
                "eedf_normalization_error": _float_or_none(
                    row["eedf_normalization_error"]
                ),
                "valid_replicates": int(row["valid_replicates"]),
                "uncertainty_available": int(row["uncertainty_available"]),
                **{
                    key: value
                    for key, value in diagnostic.items()
                    if key != "solver_failure_reasons"
                },
                "quality_source": f"{solver}_aggregate_quality",
            }
        )
    return quality


def _passing_e_over_n(
    connection: sqlite3.Connection,
    solver: str,
    mixture_id: int,
) -> set[float]:
    rows = connection.execute(
        """
        SELECT e_over_n_Td
        FROM aggregate_quality
        WHERE solver = ? AND mixture_id = ? AND passed = 1
        """,
        (solver, mixture_id),
    ).fetchall()
    return {float(row["e_over_n_Td"]) for row in rows}


def _hybrid_quality_rows(
    base_cases: list[dict[str, Any]],
    mc_quality: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    mc_by_e = {float(row["E_over_N_Td"]): row for row in mc_quality}
    quality = []
    for case in base_cases:
        e_over_n = float(case["E_over_N_Td"])
        mc_row = mc_by_e.get(e_over_n)
        quality.append(
            {
                **_en_fields(e_over_n),
                "passed": 1,
                "failure_reasons_json": "[]",
                "mobility_rse": None if mc_row is None else mc_row.get("mobility_rse"),
                "diffusion_L_rse": None
                if mc_row is None
                else mc_row.get("diffusion_L_rse"),
                "diffusion_T_rse": None
                if mc_row is None
                else mc_row.get("diffusion_T_rse"),
                "max_major_rate_rse": None
                if mc_row is None
                else mc_row.get("max_major_rate_rse"),
                "eedf_normalization_error": None
                if mc_row is None
                else mc_row.get("eedf_normalization_error"),
                "valid_replicates": 1
                if mc_row is None
                else mc_row.get("valid_replicates"),
                "uncertainty_available": 0
                if mc_row is None
                else mc_row.get("uncertainty_available"),
                "quality_source": "hybrid_two_term_baseline_with_mc_correction",
            }
        )
    return quality


def _fill_rate_derived_columns(row: dict[str, Any], drift: object | None = None) -> None:
    drift_value = _float_or_none(drift if drift is not None else row.get("drift_velocity_m_s"))
    if drift_value is None:
        drift_value = _float_or_none(row.get("case_drift_velocity_m_s"))
    if drift_value is None or abs(drift_value) <= 0.0:
        raise TableBuildError("cannot compute reduced Townsend with zero drift")
    k = _required_float(row, "rate_coefficient_m3_s")
    fraction = _required_float(row, "target_species_fraction")
    energy_loss = _float_or_none(row.get("energy_loss_eV")) or 0.0
    row["mixture_weighted_rate_m3_s"] = k * fraction
    row["reduced_townsend_m2"] = k / abs(drift_value)
    row["mixture_weighted_reduced_townsend_m2"] = k * fraction / abs(drift_value)
    row["energy_loss_rate_coefficient_eV_m3_s"] = k * energy_loss


def _validate_positive_values(dataset: _Dataset) -> None:
    positive_case = [
        "mean_energy_eV",
        "reduced_mobility_m2_V_s_m3",
        "reduced_diffusion_L_m2_s_m3",
        "reduced_diffusion_T_m2_s_m3",
    ]
    optional_positive_case = [
        "reduced_electron_energy_mobility_m2_V_s_m3",
        "reduced_electron_energy_diffusion_m2_s_m3",
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
        for name in optional_positive_case:
            if case.get(name) is not None:
                _require_nonnegative_finite(case.get(name), name)
    for rate in dataset.rates:
        for name in positive_rate:
            _require_nonnegative_finite(rate.get(name), name)
    for row in dataset.eedf:
        _require_nonnegative_finite(row.get("eedf"), "eedf")
        _require_nonnegative_finite(row.get("energy_width_eV"), "energy_width_eV")


def _require_nonnegative_finite(value: object, name: str) -> None:
    if value is None:
        raise TableBuildError(f"missing positive value {name}")
    number = float(value)
    if not math.isfinite(number) or number < 0.0:
        raise TableBuildError(f"invalid positive value {name}: {value}")


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
        "units": {column: UNITS[column] for column in columns if column in UNITS},
        "argument": argument,
    }


def _en_fields(e_over_n_Td: float) -> dict[str, float]:
    return {
        "E_over_N_Td": float(e_over_n_Td),
        "E_over_N_V_m2": float(e_over_n_Td) * TD_TO_V_M2,
    }


def _valid_ranges(cases: list[dict[str, Any]]) -> dict[str, list[float] | None]:
    return {
        "E_over_N_Td": finite_range(case.get("E_over_N_Td") for case in cases),
        "E_over_N_V_m2": finite_range(case.get("E_over_N_V_m2") for case in cases),
        "mean_energy_eV": finite_range(case.get("mean_energy_eV") for case in cases),
    }


def _quality_summary(quality: list[dict[str, Any]]) -> dict[str, Any]:
    failed = [row for row in quality if int(row.get("passed", 0)) == 0]
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
        if column in UNITS
    }
    return {column: UNITS[column] for column in sorted(columns)}


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


def _thresholds_from_quality(quality: list[dict[str, Any]]) -> QualityThresholds:
    for row in quality:
        raw = row.get("thresholds_json")
        if raw:
            try:
                data = json.loads(str(raw))
            except json.JSONDecodeError:
                continue
            return QualityThresholds(
                mobility_rse=float(data.get("mobility_rse", QualityThresholds().mobility_rse)),
                diffusion_rse=float(data.get("diffusion_rse", QualityThresholds().diffusion_rse)),
                major_rate_rse=float(data.get("major_rate_rse", QualityThresholds().major_rate_rse)),
                major_rate_fraction=float(
                    data.get("major_rate_fraction", QualityThresholds().major_rate_fraction)
                ),
                eedf_normalization_error=float(
                    data.get(
                        "eedf_normalization_error",
                        QualityThresholds().eedf_normalization_error,
                    )
                ),
            )
    return QualityThresholds()


def _case_weight(
    case: dict[str, Any],
    scalar: str,
    thresholds: QualityThresholds,
) -> float:
    if scalar == "reduced_mobility_m2_V_s_m3":
        return _rse_weight(case.get(f"{scalar}_rse"), thresholds.mobility_rse)
    if scalar == "reduced_diffusion_L_m2_s_m3":
        return _rse_weight(case.get(f"{scalar}_rse"), thresholds.diffusion_rse)
    if scalar == "reduced_diffusion_T_m2_s_m3":
        return _rse_weight(case.get(f"{scalar}_rse"), thresholds.diffusion_rse)
    return 1.0


def _rse_weight(value: object, limit: float) -> float:
    rse = _finite_nonnegative(value)
    if rse is None:
        return 0.0
    if limit <= 0.0:
        return 1.0 if rse <= 0.0 else 0.0
    return max(0.0, min(1.0, 1.0 - rse / limit))


def _rate_key(row: dict[str, Any]) -> tuple[str, str, str, str]:
    threshold = row.get("threshold_eV")
    threshold_key = "" if threshold is None else f"{float(threshold):.17g}"
    return (
        str(row["species"]),
        str(row["process"]),
        str(row["process_type"]),
        threshold_key,
    )


def _finite_positive(value: object) -> float | None:
    number = _float_or_none(value)
    if number is None or number <= 0.0:
        return None
    return number


def _finite_nonnegative(value: object) -> float | None:
    number = _float_or_none(value)
    if number is None or number < 0.0:
        return None
    return number


def _required_float(row: dict[str, Any], name: str) -> float:
    value = _float_or_none(row.get(name))
    if value is None:
        raise TableBuildError(f"missing required value {name}")
    return value


def _float_or_none(value: object) -> float | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None
