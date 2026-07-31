"""SQLite persistence for external swarm workflow runs."""

from __future__ import annotations

import json
import math
import sqlite3
from pathlib import Path
import numpy as np

from electron_swarm import SwarmCaseResult, eepf_from_eedf, widths_from_centers


LEGACY_ENERGY_TRANSPORT_COLUMNS = frozenset(
    {
        "reduced_electron_energy_mobility_eV_m2_V_s_m3",
        "reduced_electron_energy_diffusion_eV_m2_s_m3",
    }
)
CANONICAL_ENERGY_TRANSPORT_COLUMNS = frozenset(
    {
        "reduced_electron_energy_mobility_m2_V_s_m3",
        "reduced_electron_energy_diffusion_m2_s_m3",
    }
)


class WorkflowSchemaError(RuntimeError):
    """Raised when an existing workflow database uses an incompatible schema."""


def validate_workflow_schema(connection: sqlite3.Connection) -> None:
    """Reject legacy energy-transport columns instead of silently reinterpreting them."""

    exists = connection.execute(
        "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'cases'"
    ).fetchone()
    if exists is None:
        return
    columns = {
        str(row[1])
        for row in connection.execute("PRAGMA table_info(cases)").fetchall()
    }
    legacy = sorted(columns & LEGACY_ENERGY_TRANSPORT_COLUMNS)
    if legacy:
        raise WorkflowSchemaError(
            "Existing workflow database uses obsolete energy-transport columns "
            f"{legacy}. Their names incorrectly encode an eV factor. "
            "Schema v2 does not silently reinterpret or alias these fields; "
            "regenerate the database with the current workflow."
        )
    missing = sorted(CANONICAL_ENERGY_TRANSPORT_COLUMNS - columns)
    if missing:
        raise WorkflowSchemaError(
            "Existing workflow database is missing canonical schema-v2 "
            f"energy-transport columns {missing}; regenerate the database."
        )


def _json_default(value: object) -> object:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _json_dumps(value: object) -> str:
    return json.dumps(value, default=_json_default, sort_keys=True)


def _finite_or_none(value: object) -> float | None:
    if value is None:
        return None
    number = float(value)
    if not math.isfinite(number):
        return None
    return number


class WorkflowStore:
    """Small SQLite wrapper for workflow results."""

    def __init__(self, path: Path) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.connection = sqlite3.connect(self.path)
        self.connection.execute("PRAGMA foreign_keys = ON")
        try:
            validate_workflow_schema(self.connection)
            self._create_schema()
        except Exception:
            self.connection.close()
            raise

    def close(self) -> None:
        self.connection.close()

    def __enter__(self) -> "WorkflowStore":
        return self

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        self.close()

    def _create_schema(self) -> None:
        self.connection.executescript(
            """
            CREATE TABLE IF NOT EXISTS metadata (
                key TEXT PRIMARY KEY,
                value TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS mixtures (
                mixture_id INTEGER PRIMARY KEY,
                fractions_json TEXT NOT NULL
            );
            CREATE TABLE IF NOT EXISTS mixture_species (
                mixture_id INTEGER NOT NULL,
                species TEXT NOT NULL,
                fraction REAL NOT NULL,
                mass_amu REAL NOT NULL,
                PRIMARY KEY (mixture_id, species),
                FOREIGN KEY (mixture_id) REFERENCES mixtures(mixture_id)
            );
            CREATE TABLE IF NOT EXISTS cases (
                mixture_id INTEGER NOT NULL,
                solver TEXT NOT NULL,
                e_over_n_Td REAL NOT NULL,
                replicate INTEGER NOT NULL,
                case_id TEXT NOT NULL,
                mean_energy_eV REAL NOT NULL,
                gas_number_density_m3 REAL NOT NULL,
                transport_definition TEXT NOT NULL,
                drift_velocity_m_s REAL NOT NULL,
                mobility_m2_V_s REAL NOT NULL,
                reduced_mobility_m2_V_s_m3 REAL NOT NULL,
                diffusion_L_m2_s REAL NOT NULL,
                diffusion_T_m2_s REAL NOT NULL,
                reduced_diffusion_L_m2_s_m3 REAL NOT NULL,
                reduced_diffusion_T_m2_s_m3 REAL NOT NULL,
                reduced_electron_energy_mobility_m2_V_s_m3 REAL,
                reduced_electron_energy_diffusion_m2_s_m3 REAL,
                net_ionization_frequency_s REAL NOT NULL,
                effective_townsend_m2 REAL NOT NULL,
                schema_version TEXT NOT NULL,
                metadata_json TEXT NOT NULL,
                diagnostics_json TEXT NOT NULL,
                PRIMARY KEY (mixture_id, solver, e_over_n_Td, replicate),
                FOREIGN KEY (mixture_id) REFERENCES mixtures(mixture_id)
            );
            CREATE TABLE IF NOT EXISTS rates (
                mixture_id INTEGER NOT NULL,
                solver TEXT NOT NULL,
                e_over_n_Td REAL NOT NULL,
                replicate INTEGER NOT NULL,
                rate_index INTEGER NOT NULL,
                species TEXT NOT NULL,
                process TEXT NOT NULL,
                process_type TEXT NOT NULL,
                threshold_eV REAL,
                rate_coefficient_m3_s REAL NOT NULL,
                target_species_fraction REAL NOT NULL,
                energy_loss_eV REAL,
                energy_loss_rate_coefficient_eV_m3_s REAL NOT NULL,
                mixture_weighted_rate_m3_s REAL NOT NULL,
                frequency_s_inv REAL,
                power_loss_eV_s REAL,
                tail_fraction REAL,
                PRIMARY KEY (
                    mixture_id, solver, e_over_n_Td, replicate, rate_index
                ),
                FOREIGN KEY (
                    mixture_id, solver, e_over_n_Td, replicate
                ) REFERENCES cases(mixture_id, solver, e_over_n_Td, replicate)
            );
            CREATE TABLE IF NOT EXISTS eedf_bins (
                mixture_id INTEGER NOT NULL,
                solver TEXT NOT NULL,
                e_over_n_Td REAL NOT NULL,
                replicate INTEGER NOT NULL,
                bin_index INTEGER NOT NULL,
                energy_eV REAL NOT NULL,
                energy_width_eV REAL NOT NULL,
                eedf REAL NOT NULL,
                eepf REAL NOT NULL,
                sample_count INTEGER,
                effective_sample_count REAL,
                relative_standard_error REAL,
                PRIMARY KEY (
                    mixture_id, solver, e_over_n_Td, replicate, bin_index
                ),
                FOREIGN KEY (
                    mixture_id, solver, e_over_n_Td, replicate
                ) REFERENCES cases(mixture_id, solver, e_over_n_Td, replicate)
            );
            """
        )
        self.connection.commit()

    def set_provenance(self, values: dict[str, str]) -> None:
        with self.connection:
            for key, value in values.items():
                row = self.connection.execute(
                    "SELECT value FROM metadata WHERE key = ?",
                    (key,),
                ).fetchone()
                if row is not None and row[0] != value:
                    raise ValueError(
                        f"Existing workflow database has different {key}; "
                        "use a separate database"
                    )
                self.connection.execute(
                    "INSERT OR IGNORE INTO metadata(key, value) VALUES (?, ?)",
                    (key, value),
                )

    def write_mixture(
        self,
        mixture_id: int,
        species_rows: list[tuple[str, float, float]],
    ) -> None:
        fractions_json = _json_dumps(
            {species: fraction for species, fraction, _mass in species_rows}
        )
        existing = self.connection.execute(
            "SELECT fractions_json FROM mixtures WHERE mixture_id = ?",
            (mixture_id,),
        ).fetchone()
        if existing is not None and existing[0] != fractions_json:
            raise ValueError(
                f"Existing workflow database has different mixture_id {mixture_id}"
            )
        with self.connection:
            self.connection.execute(
                "INSERT OR IGNORE INTO mixtures(mixture_id, fractions_json) VALUES (?, ?)",
                (mixture_id, fractions_json),
            )
            self.connection.executemany(
                """
                INSERT OR REPLACE INTO mixture_species(
                    mixture_id, species, fraction, mass_amu
                ) VALUES (?, ?, ?, ?)
                """,
                [
                    (mixture_id, species, float(fraction), float(mass))
                    for species, fraction, mass in species_rows
                ],
            )

    def write_case(
        self,
        *,
        mixture_id: int,
        replicate: int,
        case: SwarmCaseResult,
    ) -> None:
        key = (mixture_id, case.solver, float(case.e_over_n_Td), int(replicate))
        widths = case.energy_widths_eV
        if widths is None or len(widths) != len(case.energy_eV):
            widths = widths_from_centers(case.energy_eV)
        eepf = eepf_from_eedf(case.energy_eV, case.eedf)
        counts = case.eedf_counts
        effective_counts = case.eedf_effective_counts

        with self.connection:
            self.connection.execute(
                """
                INSERT OR REPLACE INTO cases(
                    mixture_id, solver, e_over_n_Td, replicate, case_id,
                    mean_energy_eV, gas_number_density_m3, transport_definition,
                    drift_velocity_m_s, mobility_m2_V_s,
                    reduced_mobility_m2_V_s_m3, diffusion_L_m2_s,
                    diffusion_T_m2_s, reduced_diffusion_L_m2_s_m3,
                    reduced_diffusion_T_m2_s_m3,
                    reduced_electron_energy_mobility_m2_V_s_m3,
                    reduced_electron_energy_diffusion_m2_s_m3,
                    net_ionization_frequency_s, effective_townsend_m2,
                    schema_version, metadata_json, diagnostics_json
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    *key,
                    case.case_id,
                    float(case.mean_energy_eV),
                    float(case.transport.gas_number_density_m3),
                    case.transport.definition,
                    float(case.drift_velocity_m_s),
                    float(case.mobility_m2_V_s),
                    float(case.reduced_mobility_m2_V_s_m3),
                    float(case.diffusion_L_m2_s),
                    float(case.diffusion_T_m2_s),
                    float(case.reduced_diffusion_L_m2_s_m3),
                    float(case.reduced_diffusion_T_m2_s_m3),
                    _finite_or_none(
                        case.reduced_electron_energy_mobility_m2_V_s_m3
                    ),
                    _finite_or_none(
                        case.reduced_electron_energy_diffusion_m2_s_m3
                    ),
                    float(case.net_ionization_frequency_s),
                    float(case.effective_townsend_m2),
                    case.schema_version,
                    _json_dumps(case.metadata),
                    _json_dumps(case.diagnostics),
                ),
            )
            self.connection.execute(
                """
                DELETE FROM rates
                WHERE mixture_id = ? AND solver = ? AND e_over_n_Td = ?
                  AND replicate = ?
                """,
                key,
            )
            self.connection.executemany(
                """
                INSERT INTO rates(
                    mixture_id, solver, e_over_n_Td, replicate, rate_index,
                    species, process, process_type, threshold_eV,
                    rate_coefficient_m3_s, target_species_fraction, energy_loss_eV,
                    energy_loss_rate_coefficient_eV_m3_s,
                    mixture_weighted_rate_m3_s, frequency_s_inv, power_loss_eV_s,
                    tail_fraction
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [
                    (
                        *key,
                        index,
                        rate.species,
                        rate.process,
                        rate.process_type,
                        _finite_or_none(rate.threshold_eV),
                        float(rate.rate_coefficient_m3_s),
                        float(rate.target_species_fraction),
                        _finite_or_none(rate.energy_loss_eV),
                        float(rate.energy_loss_rate_coefficient_eV_m3_s),
                        float(rate.mixture_weighted_rate_m3_s),
                        _finite_or_none(rate.frequency_s_inv),
                        _finite_or_none(rate.power_loss_eV_s),
                        _finite_or_none(rate.tail_fraction),
                    )
                    for index, rate in enumerate(case.rates)
                ],
            )
            self.connection.execute(
                """
                DELETE FROM eedf_bins
                WHERE mixture_id = ? AND solver = ? AND e_over_n_Td = ?
                  AND replicate = ?
                """,
                key,
            )
            self.connection.executemany(
                """
                INSERT INTO eedf_bins(
                    mixture_id, solver, e_over_n_Td, replicate, bin_index,
                    energy_eV, energy_width_eV, eedf, eepf, sample_count,
                    effective_sample_count, relative_standard_error
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [
                    (
                        *key,
                        index,
                        float(energy),
                        float(widths[index]),
                        float(case.eedf[index]),
                        float(eepf[index]),
                        (
                            int(counts[index])
                            if counts is not None and index < len(counts)
                            else None
                        ),
                        (
                            float(effective_counts[index])
                            if effective_counts is not None
                            and index < len(effective_counts)
                            else None
                        ),
                        (
                            float(1.0 / math.sqrt(effective_counts[index]))
                            if effective_counts is not None
                            and index < len(effective_counts)
                            and effective_counts[index] > 0.0
                            else None
                        ),
                    )
                    for index, energy in enumerate(case.energy_eV)
                ],
            )
