"""Validated canonical-table and Monte Carlo EEDF readers."""

from __future__ import annotations

import csv
from contextlib import closing
import math
from pathlib import Path
import sqlite3
from typing import Any, Sequence

import numpy as np
from scipy.stats import t as student_t

from swarm_workflow.campaign.statistics import energy_edges_from_cells
from swarm_workflow.plots.eedf_contracts import (
    EedfCase,
    EedfComparisonError,
    EedfDataset,
    MAX_SOURCE_NORMALIZATION_ERROR,
)
from swarm_workflow.plots.eedf_provenance import (
    manifest_provenance as _manifest_provenance,
    read_source_manifest as _read_source_manifest,
    sha256,
    validate_manifest_artifact as validate_manifest_artifact,
    validate_mc_database_provenance,
    verify_manifest_file as _verify_manifest_file,
)


def validate_mc_qualification_manifest(
    manifest_path: str | Path,
    qualification_path: str | Path,
) -> None:
    """Bind a qualification CSV to its declared MC profile and manifest."""

    resolved_manifest = Path(manifest_path).resolve()
    resolved_qualification = Path(qualification_path).resolve()
    manifest = _read_source_manifest(resolved_manifest, solver="monte_carlo")
    qualification = manifest.get("mc_qualification")
    expected_file = (
        qualification.get("file") if isinstance(qualification, dict) else None
    )
    if expected_file != resolved_qualification.name:
        raise EedfComparisonError(
            "Monte Carlo manifest does not identify the supplied qualification CSV"
        )
    policy = manifest.get("source_policy")
    expected_profile = (
        policy.get("qualification_profile") if isinstance(policy, dict) else None
    )
    if not isinstance(expected_profile, str) or not expected_profile:
        raise EedfComparisonError(
            "Monte Carlo manifest lacks a qualification profile"
        )
    profiles = {
        str(row.get("qualification_profile", ""))
        for row in csv_rows(resolved_qualification)
    }
    if profiles != {expected_profile}:
        raise EedfComparisonError(
            "Monte Carlo qualification CSV profile differs from its manifest"
        )


def finite_float(value: object, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise EedfComparisonError(f"{label} must be numeric") from exc
    if not math.isfinite(result):
        raise EedfComparisonError(f"{label} must be finite")
    return result


def validated_case(
    *,
    solver: str,
    e_over_n_td: float,
    energy_eV: Sequence[float],
    widths_eV: Sequence[float],
    density_eV_inv: Sequence[float],
    reported_mean_energy_eV: float,
    ci95_half_density_eV_inv: Sequence[float] | None = None,
    uncertainty_label: str | None = None,
) -> EedfCase:
    energy = np.asarray(energy_eV, dtype=float)
    widths = np.asarray(widths_eV, dtype=float)
    density = np.asarray(density_eV_inv, dtype=float)
    if (
        energy.ndim != 1
        or len(energy) == 0
        or widths.shape != energy.shape
        or density.shape != energy.shape
        or np.any(~np.isfinite(density))
        or np.any(density < 0.0)
    ):
        raise EedfComparisonError(
            f"{solver} {e_over_n_td:g} Td EEDF must be finite and nonnegative"
        )
    try:
        edges = energy_edges_from_cells(energy, widths)
    except ValueError as exc:
        raise EedfComparisonError(
            f"{solver} {e_over_n_td:g} Td has an invalid cell grid: {exc}"
        ) from exc
    normalization = float(np.sum(density * widths))
    mean = finite_float(reported_mean_energy_eV, "reported mean energy")
    if normalization <= 0.0 or mean <= 0.0:
        raise EedfComparisonError(
            f"{solver} {e_over_n_td:g} Td EEDF has no positive probability mass"
        )
    if abs(normalization - 1.0) > MAX_SOURCE_NORMALIZATION_ERROR:
        raise EedfComparisonError(
            f"{solver} {e_over_n_td:g} Td EEDF normalization error exceeds "
            f"{MAX_SOURCE_NORMALIZATION_ERROR:g}: {normalization:.17g}"
        )
    half: np.ndarray | None = None
    if ci95_half_density_eV_inv is not None:
        half = np.asarray(ci95_half_density_eV_inv, dtype=float)
        if (
            half.shape != density.shape
            or np.any(~np.isfinite(half))
            or np.any(half < 0.0)
        ):
            raise EedfComparisonError(
                f"{solver} {e_over_n_td:g} Td has an invalid uncertainty band"
            )
    return EedfCase(
        solver=solver,
        e_over_n_td=float(e_over_n_td),
        edges_eV=edges,
        density_eV_inv=density,
        reported_mean_energy_eV=mean,
        ci95_half_density_eV_inv=half,
        uncertainty_label=uncertainty_label,
    )


def csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise EedfComparisonError(f"missing EEDF evidence: {path}")
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle))
    except OSError as exc:
        raise EedfComparisonError(f"could not read EEDF evidence: {path}") from exc


def load_table_eedf(
    directory: str | Path,
    solver: str,
    *,
    require_manifest: bool = False,
) -> EedfDataset:
    """Load canonical ``eedf.csv`` and ``mean_energy_vs_en.csv`` tables."""

    root = Path(directory).resolve()
    eedf_path = root / "eedf.csv"
    mean_path = root / "mean_energy_vs_en.csv"
    source_provenance: dict[str, object] = {}
    manifest_path = root / "manifest.json"
    if require_manifest or manifest_path.is_file():
        manifest = _read_source_manifest(manifest_path, solver=solver)
        _verify_manifest_file(manifest, root=root, name="eedf.csv", solver=solver)
        _verify_manifest_file(
            manifest,
            root=root,
            name="mean_energy_vs_en.csv",
            solver=solver,
        )
        source_provenance = _manifest_provenance(manifest_path, manifest)
    means: dict[float, float] = {}
    for row in csv_rows(mean_path):
        field = finite_float(row.get("E_over_N_Td"), "E/N")
        if field in means:
            raise EedfComparisonError(f"duplicate mean-energy anchor at {field:g} Td")
        means[field] = finite_float(row.get("mean_energy_eV"), "mean energy")

    grouped: dict[float, list[dict[str, str]]] = {}
    for row in csv_rows(eedf_path):
        field = finite_float(row.get("E_over_N_Td"), "E/N")
        grouped.setdefault(field, []).append(row)
    if set(grouped) != set(means):
        raise EedfComparisonError(
            f"{solver} EEDF and mean-energy anchor sets do not match"
        )

    cases: dict[float, EedfCase] = {}
    for field, rows in grouped.items():
        rows.sort(key=lambda row: finite_float(row.get("electron_energy_eV"), "energy"))
        cases[field] = validated_case(
            solver=solver,
            e_over_n_td=field,
            energy_eV=[row["electron_energy_eV"] for row in rows],
            widths_eV=[row["energy_width_eV"] for row in rows],
            density_eV_inv=[row["eedf"] for row in rows],
            reported_mean_energy_eV=means[field],
        )
    return EedfDataset(
        solver=solver,
        cases=cases,
        provenance={
            **source_provenance,
            "eedf_path": str(eedf_path),
            "eedf_sha256": sha256(eedf_path),
            "mean_energy_path": str(mean_path),
            "mean_energy_sha256": sha256(mean_path),
        },
    )


def load_mc_eedf_database(
    database_path: str | Path,
    *,
    mixture_id: int = 0,
    source_manifest_path: str | Path | None = None,
) -> EedfDataset:
    """Load aggregate MC EEDFs and pointwise replica uncertainty from SQLite."""

    path = Path(database_path).resolve()
    if not path.is_file():
        raise EedfComparisonError(f"missing Monte Carlo database: {path}")
    source_provenance: dict[str, object] = {}
    manifest: dict[str, Any] | None = None
    manifest_path: Path | None = None
    if source_manifest_path is not None:
        manifest_path = Path(source_manifest_path).resolve()
        manifest = _read_source_manifest(manifest_path, solver="monte_carlo")
        expected_database_sha256 = manifest["hashes"].get(
            "source_database_sha256"
        )
        if (
            not isinstance(expected_database_sha256, str)
            or sha256(path) != expected_database_sha256
        ):
            raise EedfComparisonError(
                "Monte Carlo database does not match its source manifest"
            )
    metadata: dict[str, str] = {}
    mixture_row: sqlite3.Row | None = None
    mixture_species_rows: list[sqlite3.Row] = []
    try:
        with closing(
            sqlite3.connect(f"file:{path.as_posix()}?mode=ro", uri=True)
        ) as connection:
            connection.row_factory = sqlite3.Row
            if manifest is not None:
                metadata = {
                    str(row[0]): str(row[1])
                    for row in connection.execute("SELECT key, value FROM metadata")
                }
                mixture_row = connection.execute(
                    "SELECT fractions_json FROM mixtures WHERE mixture_id = ?",
                    (int(mixture_id),),
                ).fetchone()
                mixture_species_rows = connection.execute(
                    """
                    SELECT species, fraction, mass_amu
                    FROM mixture_species
                    WHERE mixture_id = ?
                    ORDER BY species
                    """,
                    (int(mixture_id),),
                ).fetchall()
            mean_rows = connection.execute(
                """
                SELECT e_over_n_Td, mean
                FROM aggregate_scalars
                WHERE solver = 'monte_carlo' AND mixture_id = ?
                  AND scalar_group = 'case' AND scalar_name = 'mean_energy_eV'
                ORDER BY e_over_n_Td
                """,
                (int(mixture_id),),
            ).fetchall()
            bin_rows = connection.execute(
                """
                SELECT e_over_n_Td, bin_index, energy_eV, energy_width_eV, eedf,
                       probability_mass_standard_error, valid_replicates,
                       uncertainty_available
                FROM aggregate_eedf_bins
                WHERE solver = 'monte_carlo' AND mixture_id = ?
                ORDER BY e_over_n_Td, bin_index
                """,
                (int(mixture_id),),
            ).fetchall()
    except sqlite3.Error as exc:
        raise EedfComparisonError(
            "Monte Carlo aggregate EEDF tables are unavailable"
        ) from exc

    if manifest is not None and manifest_path is not None:
        source_provenance = validate_mc_database_provenance(
            manifest_path=manifest_path,
            manifest=manifest,
            metadata=metadata,
            mixture_id=mixture_id,
            mixture_row=mixture_row,
            mixture_species_rows=mixture_species_rows,
        )

    means = {float(row["e_over_n_Td"]): float(row["mean"]) for row in mean_rows}
    grouped: dict[float, list[sqlite3.Row]] = {}
    for row in bin_rows:
        grouped.setdefault(float(row["e_over_n_Td"]), []).append(row)
    if not grouped or set(grouped) != set(means):
        raise EedfComparisonError(
            "Monte Carlo EEDF and mean-energy anchors do not match"
        )

    cases: dict[float, EedfCase] = {}
    for field, rows in grouped.items():
        half_widths: list[float] = []
        uncertainty_available = True
        replicate_counts: set[int] = set()
        for row in rows:
            count = int(row["valid_replicates"])
            standard_error_mass = row["probability_mass_standard_error"]
            width = float(row["energy_width_eV"])
            if (
                int(row["uncertainty_available"]) != 1
                or count < 2
                or standard_error_mass is None
            ):
                uncertainty_available = False
                half_widths.append(0.0)
            else:
                replicate_counts.add(count)
                critical = float(student_t.ppf(0.975, df=count - 1))
                half_widths.append(critical * float(standard_error_mass) / width)
        if uncertainty_available and len(replicate_counts) != 1:
            raise EedfComparisonError(
                f"Monte Carlo {field:g} Td uses inconsistent replica counts by bin"
            )
        cases[field] = validated_case(
            solver="monte_carlo",
            e_over_n_td=field,
            energy_eV=[row["energy_eV"] for row in rows],
            widths_eV=[row["energy_width_eV"] for row in rows],
            density_eV_inv=[row["eedf"] for row in rows],
            reported_mean_energy_eV=means[field],
            ci95_half_density_eV_inv=(half_widths if uncertainty_available else None),
            uncertainty_label=(
                "pointwise 95% Student-t interval across independent replicas"
                if uncertainty_available
                else None
            ),
        )
    return EedfDataset(
        solver="monte_carlo",
        cases=cases,
        provenance={
            **source_provenance,
            "database_path": str(path),
            "database_sha256": sha256(path),
            "mixture_id": int(mixture_id),
            "uncertainty": (
                "pointwise 95% Student-t interval from aggregate probability-mass "
                "standard errors; it is not a simultaneous confidence region"
            ),
        },
    )


def validate_compatible_physical_sources(*datasets: EedfDataset) -> None:
    """Require one cross-section, mixture, and physical context across solvers."""

    if not datasets:
        raise EedfComparisonError("EEDF comparison has no source datasets")
    required = ("physical_context", "mixture", "cross_sections_sha256")
    reference = datasets[0]
    for key in required:
        if key not in reference.provenance:
            raise EedfComparisonError(
                f"{reference.solver} comparison source lacks {key} provenance"
            )
    def physical_value(dataset: EedfDataset, key: str) -> object:
        value = dataset.provenance[key]
        if key != "mixture":
            return value
        species = value.get("species") if isinstance(value, dict) else None
        if not isinstance(species, list):
            raise EedfComparisonError(
                f"{dataset.solver} comparison source has invalid mixture provenance"
            )
        return sorted(
            (
                str(item["species"]),
                float(item["fraction"]),
                float(item["mass_amu"]),
            )
            for item in species
        )

    for candidate in datasets[1:]:
        for key in required:
            if physical_value(candidate, key) != physical_value(reference, key):
                raise EedfComparisonError(
                    f"{candidate.solver} comparison source has mismatched {key}"
                )


def qualification_by_field(path: Path, column: str) -> dict[float, bool]:
    statuses: dict[float, bool] = {}
    for row in csv_rows(path):
        field = finite_float(row.get("E_over_N_Td"), "qualification E/N")
        raw = str(row.get(column, "")).strip().lower()
        if raw not in {"0", "1", "false", "true"}:
            raise EedfComparisonError(
                f"qualification column {column!r} is missing or non-boolean"
            )
        if field in statuses:
            raise EedfComparisonError(f"duplicate qualification anchor at {field:g} Td")
        statuses[field] = raw in {"1", "true"}
    if not statuses:
        raise EedfComparisonError("Monte Carlo qualification evidence is empty")
    return statuses


def case_at(dataset: EedfDataset, field: float) -> EedfCase:
    matches = [
        case
        for value, case in dataset.cases.items()
        if math.isclose(value, field, rel_tol=1e-12, abs_tol=1e-12)
    ]
    if len(matches) != 1:
        raise EedfComparisonError(
            f"{dataset.solver} does not contain exactly one {field:g} Td EEDF"
        )
    return matches[0]
