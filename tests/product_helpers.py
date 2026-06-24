from __future__ import annotations

from pathlib import Path

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
                "tail_rate_fraction_target": 1.0e-3,
                "tail_metrics": True,
                "tail_threshold_eV": None,
                "tail_rate_warning_fraction": 0.05,
                "max_eV_limit": 1000.0,
            },
        },
        "solvers": {
            "two_term": {"backend": "native_sg"},
            "multi_term": {
                "formulation": "PN",
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
