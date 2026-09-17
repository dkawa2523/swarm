from __future__ import annotations

from pathlib import Path
import json
import math
import sqlite3

import yaml

from electron_swarm.solvers.monte_carlo.evidence import (
    DIRECT_MC_TRANSPORT_CADENCE_TRIAL_PERIODS,
    DIRECT_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
    DIRECT_MC_TRANSPORT_LAG_PLANES,
    DIRECT_MC_TRANSPORT_OBSERVATION_PLANES,
    MC_EEDF_ESTIMATOR_SCHEMA_VERSION,
)
from electron_swarm.solvers.monte_carlo.source_identity import (
    monte_carlo_source_sha256,
)
from swarm_workflow.quality.policy import QualityThresholds, quality_thresholds_json
from swarm_workflow.campaign.store import (
    MC_SAMPLING_PLAN_METADATA_KEY,
    canonical_mc_sampling_plan_json,
)

ROOT = Path(__file__).resolve().parents[1]
XS_PATH = ROOT / "examples" / "cross_sections" / "argon_minimal.csv"


def ensure_workflow_provenance(
    conn: sqlite3.Connection,
    quality: QualityThresholds | None = None,
) -> None:
    """Attach complete synthetic provenance before inserting workflow rows."""

    values = {
        "workflow_config_sha256": "0" * 64,
        "quality_thresholds_json": quality_thresholds_json(
            quality or QualityThresholds()
        ),
    }
    for key, value in values.items():
        conn.execute(
            "INSERT OR IGNORE INTO metadata(key, value) VALUES (?, ?)",
            (key, value),
        )


def _refresh_mc_sampling_plan_provenance(conn: sqlite3.Connection) -> None:
    """Keep synthetic MC fixtures aligned with production provenance rules."""

    entries: list[dict[str, int | float]] = []
    for e_over_n, replicas in conn.execute(
        """
        SELECT e_over_n_Td, COUNT(*)
        FROM cases
        WHERE solver = 'monte_carlo'
        GROUP BY e_over_n_Td
        ORDER BY e_over_n_Td
        """
    ):
        diagnostic_row = conn.execute(
            """
            SELECT diagnostics_json
            FROM cases
            WHERE solver = 'monte_carlo' AND e_over_n_Td = ?
            ORDER BY replicate
            LIMIT 1
            """,
            (e_over_n,),
        ).fetchone()
        provenance: dict[str, object] = {}
        if diagnostic_row is not None:
            diagnostics = json.loads(str(diagnostic_row[0]))
            transport = diagnostics.get("internal_monte_carlo_transport", {})
            provenance = transport.get("mc_run_provenance", {})
        entries.append(
            {
                "e_over_n_Td": float(e_over_n),
                "particles": int(provenance.get("particles", 128)),
                "warmup_collisions": int(provenance.get("warmup_collisions", 100)),
                "max_collisions": int(provenance.get("production_collisions", 2000)),
                "tail_max_collisions": int(
                    provenance.get(
                        "tail_max_collisions",
                        provenance.get("production_collisions", 2000),
                    )
                ),
                "replicas": int(replicas),
                "transport_correlation_lag_barriers": int(
                    transport.get("block_lag_sampling", {}).get(
                        "configured_correlation_lag_barriers",
                        64,
                    )
                ),
                "transport_estimator": str(
                    provenance.get("transport_estimator", "single_field")
                ),
            }
        )
    if entries:
        conn.execute(
            "INSERT OR REPLACE INTO metadata(key, value) VALUES (?, ?)",
            (
                MC_SAMPLING_PLAN_METADATA_KEY,
                canonical_mc_sampling_plan_json(entries),
            ),
        )


def base_product_config(tmp_path: Path, solvers: list[str] | None = None) -> dict:
    return {
        "schema_version": 2,
        "run": {
            "solvers": [{"id": solver} for solver in (solvers or ["two_term"])],
            "e_over_n_Td": [50.0],
            "case_prefix": "Prod",
        },
        "conditions": {
            "gas_temperature_K": 300.0,
            "pressure_Pa": 100.0,
            "gas_mixture": [{"species": "Ar", "fraction": 1.0, "mass_amu": 39.948}],
        },
        "cross_sections": {
            "format": "csv",
            "high_energy_extrapolation": "zero",
            "files": [{"path": XS_PATH.as_posix(), "species": "Ar", "format": "csv"}],
        },
        "physics": {
            "field": {
                "type": "dc",
                "magnetic_field": {"enabled": False, "B_T": 0.0, "angle_EB_deg": 0.0},
            },
            "angular_scattering": {
                "model": "isotropic",
                "higher_moment_closure": "zero",
            },
            "electron_electron": {"enabled": False, "model": "none"},
            "ionization": {"energy_sharing": "equal"},
            "energy_grid_policy": {
                "adaptive": True,
                "threshold_refinement": False,
                "tail_probability_target": 1.0e-8,
                "tail_metrics": True,
                "tail_threshold_eV": None,
                "tail_rate_warning_fraction": 0.05,
                "max_eV_limit": 1000.0,
            },
        },
        "solvers": {
            "two_term": {"backend": "native_sg"},
            "multi_term": {
                "method": "pn_closure_direct",
                "lmax": 3,
            },
            "monte_carlo": {},
        },
        "comparison": {"enabled": False},
        "feature_policy": {
            "unsupported": "fail",
            "degraded": "record",
        },
        "output": {
            "directory": tmp_path.as_posix(),
            "base_name": "prod",
        },
    }


def write_config(tmp_path: Path, data: dict, name: str = "case.yaml") -> Path:
    path = tmp_path / name
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")
    return path


def write_moment_table(
    tmp_path: Path,
    name: str = "moments.csv",
    m1: float | None = None,
) -> Path:
    path = tmp_path / name
    if m1 is not None:
        m2 = m1 * m1
        m3 = m2 * m1
        rows = [
            "energy_eV,m0,m1,m2,m3",
            f"0,1,{m1},{m2},{m3}",
            f"10,1,{m1},{m2},{m3}",
            f"100,1,{m1},{m2},{m3}",
            f"1000,1,{m1},{m2},{m3}",
        ]
    else:
        rows = [
            "energy_eV,m0,m1,m2,m3",
            "0,1,0.0,0.0,0.0",
            "10,1,0.2,0.04,0.008",
            "100,1,0.4,0.16,0.064",
            "1000,1,0.5,0.25,0.125",
        ]
    path.write_text(
        "\n".join(rows) + "\n",
        encoding="utf-8",
    )
    return path


def ensure_workflow_mixture(
    conn: sqlite3.Connection,
    *,
    mixture_id: int = 0,
    species: str = "Ar",
    fraction: float = 1.0,
    mass_amu: float = 39.948,
) -> None:
    conn.execute(
        "INSERT OR IGNORE INTO mixtures(mixture_id, fractions_json) VALUES (?, ?)",
        (mixture_id, f'{{"{species}": {fraction}}}'),
    )
    conn.execute(
        """
        INSERT OR IGNORE INTO mixture_species(
            mixture_id, species, fraction, mass_amu
        ) VALUES (?, ?, ?, ?)
        """,
        (mixture_id, species, fraction, mass_amu),
    )


def insert_workflow_case(
    conn: sqlite3.Connection,
    *,
    mixture_id: int = 0,
    solver: str = "monte_carlo",
    e_over_n: float = 10.0,
    replicate: int = 0,
    mean_energy: float = 2.0,
    drift: float = 5.0,
    reduced_mobility: float = 12.0,
    reduced_diffusion_l: float = 22.0,
    reduced_diffusion_t: float = 33.0,
    reduced_energy_mobility: float | None = 7.0,
    reduced_energy_diffusion: float | None = 8.0,
    reduced_energy_diffusion_l: float | None = 8.0,
    reduced_energy_diffusion_t: float | None = 8.0,
    effective_townsend: float = 9.0,
    metadata: dict[str, object] | None = None,
    diagnostics: dict[str, object] | None = None,
) -> None:
    ensure_workflow_provenance(conn)
    solver_source_sha256 = (
        monte_carlo_source_sha256() if solver == "monte_carlo" else None
    )
    if solver == "monte_carlo":
        conn.execute(
            """
            INSERT OR IGNORE INTO metadata(key, value) VALUES (?, ?)
            """,
            (
                "mc_transport_estimator_schema_version",
                DIRECT_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
            ),
        )
        conn.execute(
            """
            INSERT OR IGNORE INTO metadata(key, value) VALUES (?, ?)
            """,
            ("mc_solver_source_sha256", solver_source_sha256),
        )
        conn.execute(
            """
            INSERT OR IGNORE INTO metadata(key, value) VALUES (?, ?)
            """,
            (
                "mc_eedf_estimator_schema_version",
                MC_EEDF_ESTIMATOR_SCHEMA_VERSION,
            ),
        )
    ensure_workflow_mixture(conn, mixture_id=mixture_id)
    diagnostic_payload = dict(diagnostics or {})
    direct_values = {
        "drift_velocity_m_s": drift,
        "mobility_m2_V_s": reduced_mobility / 1.0e22,
        "diffusion_L_m2_s": reduced_diffusion_l / 1.0e22,
        "diffusion_T_m2_s": reduced_diffusion_t / 1.0e22,
        "energy_mobility_m2_V_s": (
            None
            if reduced_energy_mobility is None
            else reduced_energy_mobility / 1.0e22
        ),
        "energy_diffusion_L_m2_s": (
            None
            if reduced_energy_diffusion_l is None
            else reduced_energy_diffusion_l / 1.0e22
        ),
        "energy_diffusion_T_m2_s": (
            None
            if reduced_energy_diffusion_t is None
            else reduced_energy_diffusion_t / 1.0e22
        ),
        "mean_energy_eV": mean_energy,
    }
    if (
        solver == "monte_carlo"
        and "internal_monte_carlo_transport" not in diagnostic_payload
        and all(
            value is not None and float(value) > 0.0 for value in direct_values.values()
        )
    ):
        cadence_s = 1.0e-9
        lag_s = [float(value) * cadence_s for value in DIRECT_MC_TRANSPORT_LAG_PLANES]
        diagnostic_payload["internal_monte_carlo_transport"] = {
            "estimator": "direct_mc_fixed_lag_block_helfand_energy_flux",
            "estimator_schema_version": DIRECT_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
            "energy_transport_status": "direct_mc_candidate_for_ensemble_qualification",
            "finite": True,
            "strictly_positive": True,
            "multiple_time_origins": True,
            "complete_observation": True,
            "origin_stationarity_window_mean_early": {
                "mean_energy_eV": float(mean_energy) * 0.99,
            },
            "origin_stationarity_window_mean_late": {
                "mean_energy_eV": float(mean_energy),
            },
            "component_estimators": {
                "drift_mobility": (
                    "production_trajectory_displacement_over_residence_time"
                ),
                "particle_diffusion": "block_helfand_mean_square_displacement",
                "energy_mobility": "production_residence_energy_current",
                "energy_diffusion": ("restricted_density_packet_energy_flux_fixed_lag"),
            },
            "production_estimates": {
                name: float(direct_values[name])
                for name in (
                    "drift_velocity_m_s",
                    "mobility_m2_V_s",
                    "energy_mobility_m2_V_s",
                    "mean_energy_eV",
                )
            },
            "lag_scan": {
                "lag_planes": list(DIRECT_MC_TRANSPORT_LAG_PLANES),
                "lag_s": lag_s,
                "estimates": {
                    name: [
                        float(value) * 0.97,
                        float(value) * 0.98,
                        float(value) * 0.99,
                        float(value),
                    ]
                    for name, value in direct_values.items()
                    if name
                    in {
                        "diffusion_L_m2_s",
                        "diffusion_T_m2_s",
                        "energy_diffusion_L_m2_s",
                        "energy_diffusion_T_m2_s",
                    }
                },
                "completed_origins": [
                    DIRECT_MC_TRANSPORT_OBSERVATION_PLANES - value
                    for value in DIRECT_MC_TRANSPORT_LAG_PLANES
                ],
                "selected_lag_index": len(DIRECT_MC_TRANSPORT_LAG_PLANES) - 1,
            },
            "transport_sampling": {
                "lag_grid_source": (
                    "fixed_physical_time_independent_of_production_length"
                ),
                "origin_cadence_s": cadence_s,
                "origin_cadence_trial_periods": (
                    DIRECT_MC_TRANSPORT_CADENCE_TRIAL_PERIODS
                ),
                "rolling_origins": True,
                "production_length_controls_lag": False,
                "warmup_accumulators_reset": True,
                "observation_planes_requested": (
                    DIRECT_MC_TRANSPORT_OBSERVATION_PLANES
                ),
                "observation_planes_complete": (DIRECT_MC_TRANSPORT_OBSERVATION_PLANES),
                "minimum_observation_planes": (DIRECT_MC_TRANSPORT_OBSERVATION_PLANES),
            },
            "mc_run_provenance": {
                "seed": 1000 + int(replicate),
                "particles": 128,
                "warmup_collisions": 100,
                "production_collisions": 2000,
                "tail_max_collisions": 2000,
                "tail_rate_rse_trigger": 0.25,
                "tail_collisions_executed": 0,
                "population_model": "fixed_particle_single_daughter",
                "transport_estimator": "single_field",
                "transport_field_legs": 1,
                "transport_budget_interpretation": ("configured_budget_per_field_leg"),
                "gas_number_density_m3": 1.0e22,
                "electric_field_V_m": float(e_over_n) * 10.0,
                "trial_collision_frequency_s_inv": 1.0e9,
                "collision_clock": "global",
                "magnetic_enabled": False,
                "magnetic_B_T": 0.0,
                "magnetic_angle_EB_deg": 0.0,
                "angular_scattering_model": "isotropic",
                "ionization_source_model": "equal",
                "high_energy_extrapolation": "zero",
                "max_energy_limit_eV": 1000.0,
            },
        }
    if solver == "monte_carlo":
        assert solver_source_sha256 is not None
        transport = diagnostic_payload.get("internal_monte_carlo_transport")
        if isinstance(transport, dict):
            transport = dict(transport)
            run_provenance = dict(transport.get("mc_run_provenance", {}))
            run_provenance.setdefault(
                "solver_source_sha256",
                solver_source_sha256,
            )
            transport["mc_run_provenance"] = run_provenance
            diagnostic_payload["internal_monte_carlo_transport"] = transport
    conn.execute(
        """
        INSERT INTO cases(
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
            mixture_id,
            solver,
            e_over_n,
            replicate,
            f"{solver}_{e_over_n}_{replicate}",
            mean_energy,
            1.0e22,
            "fixed_particle_flux",
            drift,
            reduced_mobility / 1.0e22,
            reduced_mobility,
            reduced_diffusion_l / 1.0e22,
            reduced_diffusion_t / 1.0e22,
            reduced_diffusion_l,
            reduced_diffusion_t,
            reduced_energy_mobility,
            reduced_energy_diffusion,
            reduced_energy_diffusion_l,
            reduced_energy_diffusion_t,
            0.0,
            effective_townsend,
            "2",
            json.dumps(metadata or {}, sort_keys=True),
            json.dumps(diagnostic_payload, sort_keys=True),
        ),
    )
    if solver == "monte_carlo":
        _refresh_mc_sampling_plan_provenance(conn)


def insert_workflow_rate(
    conn: sqlite3.Connection,
    *,
    mixture_id: int = 0,
    solver: str = "monte_carlo",
    e_over_n: float = 10.0,
    replicate: int = 0,
    value: float = 1.0e-15,
    target_fraction: float = 1.0,
    process: str = "ionization",
    process_type: str = "ionization",
    threshold_eV: float = 15.0,
    rate_index: int = 0,
) -> None:
    energy_loss = threshold_eV
    conn.execute(
        """
        INSERT INTO rates(
            mixture_id, solver, e_over_n_Td, replicate, rate_index,
            species, process, process_type, threshold_eV, rate_coefficient_m3_s,
            target_species_fraction, energy_loss_eV,
            energy_loss_rate_coefficient_eV_m3_s, mixture_weighted_rate_m3_s,
            frequency_s_inv, power_loss_eV_s, tail_fraction
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            mixture_id,
            solver,
            e_over_n,
            replicate,
            rate_index,
            "Ar",
            process,
            process_type,
            threshold_eV,
            value,
            target_fraction,
            energy_loss,
            value * energy_loss,
            value * target_fraction,
            None,
            None,
            None,
        ),
    )
    if solver == "monte_carlo":
        case = conn.execute(
            """
            SELECT diagnostics_json FROM cases
            WHERE mixture_id = ? AND solver = ? AND e_over_n_Td = ?
              AND replicate = ?
            """,
            (mixture_id, solver, e_over_n, replicate),
        ).fetchone()
        if case is not None:
            diagnostics = json.loads(str(case[0]))
            section = diagnostics.setdefault(
                "internal_monte_carlo_reaction_rates",
                {
                    "estimator": "trajectory_time_average_sigma_v",
                    "rates": [],
                },
            )
            rows = section.setdefault("rates", [])
            item = next(
                (
                    row
                    for row in rows
                    if row.get("species") == "Ar"
                    and row.get("process") == process
                    and row.get("process_type") == process_type
                    and float(row.get("threshold_eV", math.nan)) == threshold_eV
                ),
                None,
            )
            if item is None:
                item = {
                    "species": "Ar",
                    "process": process,
                    "process_type": process_type,
                    "threshold_eV": threshold_eV,
                }
                rows.append(item)
            item.setdefault("target_species_fraction", target_fraction)
            item.setdefault("event_sampling_enabled", True)
            item.setdefault("event_count", 1)
            item.setdefault("event_observation_residence_time_s", 1.0e-6)
            item.setdefault("zero_event_confidence", 0.95)
            item["rate_coefficient_m3_s"] = value
            item["histogram_convolution_rate_coefficient_m3_s"] = value
            conn.execute(
                """
                UPDATE cases SET diagnostics_json = ?
                WHERE mixture_id = ? AND solver = ? AND e_over_n_Td = ?
                  AND replicate = ?
                """,
                (
                    json.dumps(diagnostics, sort_keys=True),
                    mixture_id,
                    solver,
                    e_over_n,
                    replicate,
                ),
            )


def insert_workflow_eedf(
    conn: sqlite3.Connection,
    *,
    mixture_id: int = 0,
    solver: str = "monte_carlo",
    e_over_n: float = 10.0,
    replicate: int = 0,
    probability_masses: tuple[float, ...] = (0.4, 0.6),
) -> None:
    for index, mass in enumerate(probability_masses):
        energy = 0.5 + index
        conn.execute(
            """
            INSERT INTO eedf_bins(
                mixture_id, solver, e_over_n_Td, replicate, bin_index,
                energy_eV, energy_width_eV, eedf, eepf, sample_count,
                effective_sample_count, relative_standard_error
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                mixture_id,
                solver,
                e_over_n,
                replicate,
                index,
                energy,
                1.0,
                mass,
                mass,
                None,
                None,
                None,
            ),
        )
