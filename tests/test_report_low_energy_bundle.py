from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "reports"
    / "comsol_swarm_benchmark_2026"
    / "repro"
    / "prepare_low_energy_bundle.py"
)


def _load_module():
    spec = importlib.util.spec_from_file_location(
        "prepare_low_energy_bundle_for_test", SCRIPT
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_low_energy_floor_uses_zero_drift_at_zero_reduced_field(
    tmp_path: Path,
) -> None:
    module = _load_module()
    source = tmp_path / "source"
    output = tmp_path / "derived"
    source.mkdir()
    (source / "manifest.json").write_text(
        json.dumps(
            {
                "valid_ranges": {"mean_energy_eV": [0.5, 10.0]},
            }
        ),
        encoding="utf-8",
    )
    (source / "transport_vs_mean_energy.csv").write_text(
        "mean_energy_eV,E_over_N_Td,E_over_N_V_m2,drift_velocity_m_s,"
        "reduced_mobility_m2_V_s_m3\n"
        "0.5,0.05,5e-23,1000,2e25\n"
        "1.0,0.1,1e-22,1500,1.5e25\n",
        encoding="utf-8",
    )
    (source / "rates_vs_mean_energy.csv").write_text(
        "mean_energy_eV,E_over_N_Td,E_over_N_V_m2,process_type,"
        "rate_coefficient_m3_s,mixture_weighted_rate_m3_s,"
        "reduced_townsend_m2,mixture_weighted_reduced_townsend_m2\n"
        "0.5,0.05,5e-23,elastic,1,1,2,2\n"
        "0.5,0.05,5e-23,excitation,3,3,4,4\n"
        "1.0,0.1,1e-22,elastic,5,5,6,6\n"
        "1.0,0.1,1e-22,excitation,7,7,8,8\n",
        encoding="utf-8",
    )

    module.create_bundle(source, output)

    rows = (
        (output / "transport_vs_mean_energy.csv")
        .read_text(encoding="utf-8")
        .splitlines()
    )
    header = rows[0].split(",")
    velocity_index = header.index("drift_velocity_m_s")
    field_index = header.index("E_over_N_Td")
    for row in rows[1:6]:
        values = row.split(",")
        assert float(values[field_index]) == 0.0
        assert float(values[velocity_index]) == 0.0

    manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    policy = manifest["postprocess"]["low_mean_energy_floor"]
    assert policy["drift_velocity_policy"].startswith("zero on derived E/N=0")
    assert manifest["derived_table_hashes_sha256"][
        "transport_vs_mean_energy.csv"
    ] == _sha256(output / "transport_vs_mean_energy.csv")


def test_low_energy_bundle_refuses_existing_output(tmp_path: Path) -> None:
    module = _load_module()
    source = tmp_path / "source"
    output = tmp_path / "derived"
    source.mkdir()
    output.mkdir()
    (source / "manifest.json").write_text("{}", encoding="utf-8")

    with pytest.raises(FileExistsError, match="output already exists"):
        module.create_bundle(source, output)
