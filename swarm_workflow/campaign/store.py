"""SQLite persistence for external swarm workflow runs."""

from __future__ import annotations

import json
from hashlib import sha256
import math
import sqlite3
from pathlib import Path
import numpy as np

from electron_swarm import SwarmCaseResult, eepf_from_eedf, widths_from_centers

from ..quality.propagator_source import PROPAGATOR_SOLVER_SOURCE_METADATA_KEY
from .repository import read_metadata


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
        "reduced_electron_energy_diffusion_L_m2_s_m3",
        "reduced_electron_energy_diffusion_T_m2_s_m3",
    }
)
NULLABLE_PARTICLE_TRANSPORT_COLUMNS = frozenset(
    {
        "diffusion_L_m2_s",
        "diffusion_T_m2_s",
        "reduced_diffusion_L_m2_s_m3",
        "reduced_diffusion_T_m2_s_m3",
    }
)
MC_SAMPLING_PLAN_METADATA_KEY = "mc_sampling_plan_json"
MC_SAMPLING_PLAN_FIELDS = frozenset(
    {
        "e_over_n_Td",
        "particles",
        "warmup_collisions",
        "max_collisions",
        "tail_max_collisions",
        "replicas",
        "transport_correlation_lag_barriers",
        "transport_estimator",
    }
)


class WorkflowSchemaError(RuntimeError):
    """Raised when an existing workflow database uses an incompatible schema."""


def canonical_mc_sampling_plan_json(value: object) -> str:
    """Validate and serialize the resolved per-anchor MC sampling plan."""

    if not isinstance(value, list) or not value:
        raise WorkflowSchemaError(
            "MC sampling-plan provenance must be a non-empty list"
        )
    normalized: list[dict[str, int | float | str]] = []
    anchors: set[float] = set()
    for index, item in enumerate(value):
        field_name = f"MC sampling-plan provenance row {index}"
        if not isinstance(item, dict) or set(item) != MC_SAMPLING_PLAN_FIELDS:
            raise WorkflowSchemaError(
                f"{field_name} must contain exactly {sorted(MC_SAMPLING_PLAN_FIELDS)}"
            )
        raw_anchor = item["e_over_n_Td"]
        if isinstance(raw_anchor, bool) or not isinstance(raw_anchor, (int, float)):
            raise WorkflowSchemaError(f"{field_name}.e_over_n_Td must be numeric")
        anchor = float(raw_anchor)
        if not math.isfinite(anchor) or anchor <= 0.0:
            raise WorkflowSchemaError(
                f"{field_name}.e_over_n_Td must be finite and positive"
            )
        if anchor in anchors:
            raise WorkflowSchemaError(
                f"MC sampling-plan provenance has duplicate anchor {anchor:g} Td"
            )
        anchors.add(anchor)

        controls: dict[str, int] = {}
        for name in (
            "particles",
            "warmup_collisions",
            "max_collisions",
            "tail_max_collisions",
            "replicas",
        ):
            raw_control = item[name]
            if isinstance(raw_control, bool) or not isinstance(raw_control, int):
                raise WorkflowSchemaError(f"{field_name}.{name} must be an integer")
            control = int(raw_control)
            minimum = (
                0 if name in {"warmup_collisions", "tail_max_collisions"} else 1
            )
            if control < minimum:
                qualifier = "nonnegative" if minimum == 0 else "positive"
                raise WorkflowSchemaError(f"{field_name}.{name} must be {qualifier}")
            controls[name] = control
        raw_lag = item["transport_correlation_lag_barriers"]
        if isinstance(raw_lag, bool) or not isinstance(raw_lag, int):
            raise WorkflowSchemaError(
                f"{field_name}.transport_correlation_lag_barriers must be an integer"
            )
        lag = int(raw_lag)
        if lag < 4 or lag & (lag - 1):
            raise WorkflowSchemaError(
                f"{field_name}.transport_correlation_lag_barriers must be a "
                "power of two greater than or equal to 4"
            )
        controls["transport_correlation_lag_barriers"] = lag
        transport_estimator = item["transport_estimator"]
        if transport_estimator not in {"single_field", "paired_field_parity"}:
            raise WorkflowSchemaError(
                f"{field_name}.transport_estimator must be single_field or "
                "paired_field_parity"
            )
        normalized.append(
            {
                "e_over_n_Td": anchor,
                **controls,
                "transport_estimator": transport_estimator,
            }
        )
    return json.dumps(
        normalized,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def mc_sampling_plan_provenance(
    metadata: dict[str, str],
    *,
    required: bool = False,
) -> dict[str, object] | None:
    """Return canonical sampling-plan provenance and reject altered metadata."""

    raw = metadata.get(MC_SAMPLING_PLAN_METADATA_KEY)
    if raw is None:
        if required:
            raise WorkflowSchemaError(
                "Workflow provenance is missing mc_sampling_plan_json; "
                "regenerate the Monte Carlo database"
            )
        return None
    try:
        payload = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise WorkflowSchemaError(
            "Workflow provenance has invalid mc_sampling_plan_json"
        ) from exc
    canonical = canonical_mc_sampling_plan_json(payload)
    if raw != canonical:
        raise WorkflowSchemaError(
            "Workflow mc_sampling_plan_json provenance is not canonical"
        )
    return {
        "status": "canonical",
        "sha256": sha256(canonical.encode("utf-8")).hexdigest(),
        "entries": payload,
    }


def provenance_hash_manifest(
    metadata: dict[str, str],
    *,
    source: str | None = None,
) -> dict[str, str | None]:
    """Return the canonical hash chain for a database or selected solver."""

    workflow_key = "workflow_config_sha256"
    if metadata and workflow_key not in metadata:
        raise WorkflowSchemaError(
            "Workflow provenance is missing workflow_config_sha256; "
            "regenerate or explicitly migrate the database"
        )
    manifest = {
        workflow_key: metadata.get(workflow_key),
        "base_config_sha256": metadata.get("base_config_sha256"),
        "cross_sections_sha256": metadata.get("cross_sections_sha256"),
    }
    include_mc = source in (None, "monte_carlo")
    include_propagator = source in (None, "propagator")
    estimator_key = "mc_transport_estimator_schema_version"
    if include_mc and estimator_key in metadata:
        manifest[estimator_key] = metadata[estimator_key]
        sampling = mc_sampling_plan_provenance(metadata, required=True)
        assert sampling is not None
        manifest[MC_SAMPLING_PLAN_METADATA_KEY] = metadata[
            MC_SAMPLING_PLAN_METADATA_KEY
        ]
        manifest["mc_sampling_plan_sha256"] = str(sampling["sha256"])
    eedf_estimator_key = "mc_eedf_estimator_schema_version"
    if include_mc and eedf_estimator_key in metadata:
        manifest[eedf_estimator_key] = metadata[eedf_estimator_key]
    tail_estimator_key = "mc_tail_estimator_schema_version"
    if include_mc and tail_estimator_key in metadata:
        manifest[tail_estimator_key] = metadata[tail_estimator_key]
    seed_derivation_key = "mc_seed_derivation_schema_version"
    if include_mc and seed_derivation_key in metadata:
        manifest[seed_derivation_key] = metadata[seed_derivation_key]
    solver_source_key = "mc_solver_source_sha256"
    if include_mc and solver_source_key in metadata:
        manifest[solver_source_key] = metadata[solver_source_key]
    propagator_source_key = PROPAGATOR_SOLVER_SOURCE_METADATA_KEY
    if include_propagator and propagator_source_key in metadata:
        manifest[propagator_source_key] = metadata[propagator_source_key]
    return manifest


def validate_workflow_schema(connection: sqlite3.Connection) -> None:
    """Reject legacy energy-transport columns instead of silently reinterpreting them."""

    exists = connection.execute(
        "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = 'cases'"
    ).fetchone()
    if exists is None:
        return
    table_info = connection.execute("PRAGMA table_info(cases)").fetchall()
    columns = {str(row[1]) for row in table_info}
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
    nonnullable = sorted(
        str(row[1])
        for row in table_info
        if str(row[1]) in NULLABLE_PARTICLE_TRANSPORT_COLUMNS and bool(row[3])
    )
    if nonnullable:
        raise WorkflowSchemaError(
            "Existing workflow database requires particle diffusion values in "
            f"{nonnullable}. Propagator P1 reports these quantities as "
            "unavailable; regenerate the database with the nullable schema."
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

    def metadata(self) -> dict[str, str]:
        return read_metadata(self.connection)

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
                diffusion_L_m2_s REAL,
                diffusion_T_m2_s REAL,
                reduced_diffusion_L_m2_s_m3 REAL,
                reduced_diffusion_T_m2_s_m3 REAL,
                reduced_electron_energy_mobility_m2_V_s_m3 REAL,
                reduced_electron_energy_diffusion_m2_s_m3 REAL,
                reduced_electron_energy_diffusion_L_m2_s_m3 REAL,
                reduced_electron_energy_diffusion_T_m2_s_m3 REAL,
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
            CREATE TABLE IF NOT EXISTS energy_angle_bins (
                mixture_id INTEGER NOT NULL,
                solver TEXT NOT NULL,
                e_over_n_Td REAL NOT NULL,
                replicate INTEGER NOT NULL,
                energy_index INTEGER NOT NULL,
                polar_index INTEGER NOT NULL,
                energy_eV REAL NOT NULL,
                energy_width_eV REAL NOT NULL,
                mu REAL NOT NULL,
                mu_width REAL NOT NULL,
                energy_angle_density_eV_inv REAL NOT NULL,
                PRIMARY KEY (
                    mixture_id, solver, e_over_n_Td, replicate,
                    energy_index, polar_index
                ),
                FOREIGN KEY (
                    mixture_id, solver, e_over_n_Td, replicate
                ) REFERENCES cases(mixture_id, solver, e_over_n_Td, replicate)
            );
            """
        )
        self.connection.commit()

    def set_provenance(self, values: dict[str, str]) -> None:
        if MC_SAMPLING_PLAN_METADATA_KEY in values:
            mc_sampling_plan_provenance(values, required=True)
        existing = {
            str(row[0]): str(row[1])
            for row in self.connection.execute(
                "SELECT key, value FROM metadata"
            ).fetchall()
        }
        has_workflow_rows = any(
            self.connection.execute(f"SELECT 1 FROM {table_name} LIMIT 1").fetchone()
            is not None
            for table_name in ("mixtures", "cases")
        )
        required_keys = [
            "workflow_config_sha256",
            "quality_thresholds_json",
        ]
        if "mc_transport_estimator_schema_version" in values:
            required_keys.append("mc_transport_estimator_schema_version")
        if "mc_eedf_estimator_schema_version" in values:
            required_keys.append("mc_eedf_estimator_schema_version")
        if "mc_tail_estimator_schema_version" in values:
            required_keys.append("mc_tail_estimator_schema_version")
        if "mc_seed_derivation_schema_version" in values:
            required_keys.append("mc_seed_derivation_schema_version")
        if "mc_solver_source_sha256" in values:
            required_keys.append("mc_solver_source_sha256")
        if PROPAGATOR_SOLVER_SOURCE_METADATA_KEY in values:
            required_keys.append(PROPAGATOR_SOLVER_SOURCE_METADATA_KEY)
        if "mc_sampling_plan_json" in values:
            required_keys.append("mc_sampling_plan_json")
        for required_key in required_keys:
            if (
                required_key in values
                and required_key not in existing
                and (existing or has_workflow_rows)
            ):
                raise WorkflowSchemaError(
                    "Existing workflow database is missing "
                    f"{required_key} provenance and cannot be resumed safely; "
                    "regenerate or explicitly migrate the database"
                )
        for key, value in values.items():
            if key in existing and existing[key] != value:
                raise ValueError(
                    f"Existing workflow database has different {key}; "
                    "use a separate database"
                )

        with self.connection:
            for key, value in values.items():
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

    def existing_case_keys(self) -> set[tuple[int, str, float, int]]:
        """Return committed case keys that are safe resume checkpoints."""

        return {
            (int(row[0]), str(row[1]), float(row[2]), int(row[3]))
            for row in self.connection.execute(
                """
                SELECT mixture_id, solver, e_over_n_Td, replicate
                FROM cases
                """
            )
        }

    def import_cases(
        self,
        source: sqlite3.Connection,
        keys: list[tuple[int, str, float, int]],
    ) -> int:
        """Import complete, already-qualified raw cases from another v2 store.

        Eligibility belongs to the sweep planner because it depends on the
        requested physics and MC sampling plan.  This method only performs the
        storage operation, preserving each parent case and all rate/EEDF rows
        in one transaction.
        """

        if not keys:
            return 0
        tables = ("cases", "rates", "eedf_bins", "energy_angle_bins")
        columns: dict[str, list[str]] = {}
        for table in tables:
            target_columns = [
                str(row[1])
                for row in self.connection.execute(
                    f"PRAGMA table_info({table})"
                ).fetchall()
            ]
            source_columns = [
                str(row[1])
                for row in source.execute(f"PRAGMA table_info({table})").fetchall()
            ]
            if source_columns != target_columns:
                raise WorkflowSchemaError(
                    f"Reusable workflow database has incompatible {table} schema"
                )
            columns[table] = target_columns

        imported = 0
        with self.connection:
            for key in keys:
                parent = source.execute(
                    """
                    SELECT * FROM cases
                    WHERE mixture_id = ? AND solver = ? AND e_over_n_Td = ?
                      AND replicate = ?
                    """,
                    key,
                ).fetchone()
                if parent is None:
                    continue
                for table in tables:
                    rows = (
                        [parent]
                        if table == "cases"
                        else source.execute(
                            f"""
                            SELECT * FROM {table}
                            WHERE mixture_id = ? AND solver = ? AND e_over_n_Td = ?
                              AND replicate = ?
                            ORDER BY {
                                "rate_index"
                                if table == "rates"
                                else "bin_index"
                                if table == "eedf_bins"
                                else "energy_index, polar_index"
                            }
                            """,
                            key,
                        ).fetchall()
                    )
                    placeholders = ", ".join("?" for _ in columns[table])
                    self.connection.executemany(
                        f"INSERT INTO {table} VALUES ({placeholders})",
                        rows,
                    )
                imported += 1
        return imported

    def has_case(
        self,
        *,
        mixture_id: int,
        solver: str,
        e_over_n_Td: float,
        replicate: int,
    ) -> bool:
        key = (int(mixture_id), solver, float(e_over_n_Td), int(replicate))
        return (
            self.connection.execute(
                """
            SELECT 1 FROM cases
            WHERE mixture_id = ? AND solver = ? AND e_over_n_Td = ?
              AND replicate = ?
            """,
                key,
            ).fetchone()
            is not None
        )

    def write_case(
        self,
        *,
        mixture_id: int,
        replicate: int,
        case: SwarmCaseResult,
    ) -> None:
        """Atomically commit a case and all of its rate/EEDF children.

        A later workflow connection sees the parent row only after this
        transaction commits, so a visible ``cases`` key is the complete-case
        checkpoint.
        """

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
                    reduced_electron_energy_diffusion_L_m2_s_m3,
                    reduced_electron_energy_diffusion_T_m2_s_m3,
                    net_ionization_frequency_s, effective_townsend_m2,
                    schema_version, metadata_json, diagnostics_json
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
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
                    _finite_or_none(case.diffusion_L_m2_s),
                    _finite_or_none(case.diffusion_T_m2_s),
                    _finite_or_none(case.reduced_diffusion_L_m2_s_m3),
                    _finite_or_none(case.reduced_diffusion_T_m2_s_m3),
                    _finite_or_none(case.reduced_electron_energy_mobility_m2_V_s_m3),
                    _finite_or_none(case.reduced_electron_energy_diffusion_m2_s_m3),
                    _finite_or_none(case.reduced_electron_energy_diffusion_L_m2_s_m3),
                    _finite_or_none(case.reduced_electron_energy_diffusion_T_m2_s_m3),
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
            self.connection.execute(
                """
                DELETE FROM energy_angle_bins
                WHERE mixture_id = ? AND solver = ? AND e_over_n_Td = ?
                  AND replicate = ?
                """,
                key,
            )
            distribution = case.energy_angle_distribution
            if distribution is not None:
                expected_shape = (
                    len(distribution.energy_eV),
                    len(distribution.mu),
                )
                if distribution.density_eV_inv.shape != expected_shape:
                    raise ValueError(
                        "energy-angle density shape does not match its coordinates"
                    )
                self.connection.executemany(
                    """
                    INSERT INTO energy_angle_bins(
                        mixture_id, solver, e_over_n_Td, replicate,
                        energy_index, polar_index, energy_eV,
                        energy_width_eV, mu, mu_width,
                        energy_angle_density_eV_inv
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    [
                        (
                            *key,
                            energy_index,
                            polar_index,
                            float(energy),
                            float(distribution.energy_widths_eV[energy_index]),
                            float(mu),
                            float(distribution.mu_widths[polar_index]),
                            float(
                                distribution.density_eV_inv[energy_index, polar_index]
                            ),
                        )
                        for energy_index, energy in enumerate(distribution.energy_eV)
                        for polar_index, mu in enumerate(distribution.mu)
                    ],
                )
