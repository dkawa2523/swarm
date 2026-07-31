from __future__ import annotations

from pathlib import Path
import sqlite3

import yaml

ROOT = Path(__file__).resolve().parents[1]
XS_PATH = ROOT / "examples" / "cross_sections" / "argon_minimal.csv"


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
            "gas_mixture": [
                {"species": "Ar", "fraction": 1.0, "mass_amu": 39.948}
            ],
        },
        "cross_sections": {
            "format": "csv",
            "high_energy_extrapolation": "zero",
            "files": [
                {"path": XS_PATH.as_posix(), "species": "Ar", "format": "csv"}
            ],
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
    effective_townsend: float = 9.0,
) -> None:
    ensure_workflow_mixture(conn, mixture_id=mixture_id)
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
            net_ionization_frequency_s, effective_townsend_m2,
            schema_version, metadata_json, diagnostics_json
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
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
            0.0,
            effective_townsend,
            "2",
            "{}",
            "{}",
        ),
    )


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
