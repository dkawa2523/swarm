from __future__ import annotations

import csv
import json
import sqlite3
from pathlib import Path

import yaml

from electron_swarm import load_config
from swarm_workflow.cli import main as workflow_cli_main
from swarm_workflow.comsol_adapter import plan_apply_comsol
from swarm_workflow.comsol_mapping import load_comsol_mapping
from swarm_workflow.comsol_verify import plan_verify_comsol_functions
from swarm_workflow.sweep import load_workflow


ROOT = Path(__file__).resolve().parents[1]


def test_examples_workflow_ar_o2_loads_as_schema_v2_workflow() -> None:
    workflow = load_workflow(ROOT / "examples" / "workflow_ar_o2.yaml")

    assert workflow.base_config_path.name == "ar_o2_base.yaml"
    assert workflow.e_over_n_Td == (20.0, 50.0)
    assert len(workflow.mixtures) == 2
    assert workflow.mixtures[1].fractions["O2"] == 0.01


def test_comsol_benchmark_workflow_matches_positive_column_conditions() -> None:
    config_path = ROOT / "examples" / "argon_comsol_benchmark_2026.yaml"
    workflow_path = (
        ROOT / "examples" / "workflow_argon_comsol_benchmark_2026.yaml"
    )

    config = load_config(config_path)
    workflow = load_workflow(workflow_path)

    assert [str(solver.id) for solver in config.run.solvers] == ["two_term"]
    assert config.conditions.gas_temperature_K == 293.15
    assert config.conditions.pressure_Pa == 13.3322
    assert min(config.run.e_over_n_Td) <= 0.05
    assert max(config.run.e_over_n_Td) >= 5000.0
    assert workflow.base_config_path == config_path.resolve()
    assert workflow.e_over_n_Td == tuple(config.run.e_over_n_Td)


def test_swarm_to_comsol_dry_run_end_to_end(tmp_path: Path) -> None:
    workflow_path = _write_workflow_repo(tmp_path)
    db_path = tmp_path / "outputs" / "swarm.sqlite"
    aggregate_dir = tmp_path / "outputs" / "aggregate"
    tables_dir = tmp_path / "outputs" / "tables"
    bundle_root = tmp_path / "outputs" / "comsol_bundle"

    workflow_cli_main(["sweep", str(workflow_path)])
    assert db_path.exists()
    with sqlite3.connect(db_path) as conn:
        assert conn.execute("SELECT COUNT(*) FROM cases").fetchone()[0] == 4

    workflow_cli_main(["aggregate", str(db_path), "--output", str(aggregate_dir)])
    assert (aggregate_dir / "manifest.json").exists()
    assert (aggregate_dir / "aggregate_scalars.csv").exists()

    workflow_cli_main(
        [
            "build-tables",
            str(db_path),
            "--output",
            str(tables_dir),
            "--source",
            "two_term",
        ]
    )
    assert (tables_dir / "manifest.json").exists()

    workflow_cli_main(["export-comsol", str(tables_dir), "--output", str(bundle_root)])
    bundle = bundle_root / "mixture_0000"
    assert (bundle / "manifest.json").exists()
    manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["status"] == "ok"

    mapping_path = _write_mapping(tmp_path, bundle)
    mapping = load_comsol_mapping(mapping_path, validate_files=True)
    assert mapping.bundle.path == bundle.resolve()

    apply_plan = plan_apply_comsol(mapping_path)
    verify_plan = plan_verify_comsol_functions(mapping_path)

    assert apply_plan.mapping.model.input_mph.exists()
    assert [function.tag for function in apply_plan.mapping.functions]
    assert verify_plan.points

    workflow_cli_main(
        ["run-comsol", str(mapping_path), "--bundle", str(bundle), "--dry-run"]
    )


def test_example_positive_column_mapping_parses_without_local_files() -> None:
    mapping = load_comsol_mapping(
        ROOT / "Model" / "maps" / "positive_column_external.yaml"
    )

    assert mapping.model.input_mph.name == "positive_column_1d.mph"
    assert mapping.bundle.path.name == "mixture_0000"
    assert any(function.tag == "sw_muN" for function in mapping.functions)


def _write_workflow_repo(root: Path) -> Path:
    (root / "pyproject.toml").write_text("[project]\nname='tmp'\n", encoding="utf-8")
    xs_path = root / "cross_sections" / "ar_o2_minimal.csv"
    xs_path.parent.mkdir()
    _write_cross_sections(xs_path)
    base_path = root / "base.yaml"
    base_path.write_text(
        yaml.safe_dump(
            {
                "schema_version": 2,
                "run": {
                    "solvers": [{"id": "two_term"}],
                    "e_over_n_Td": [50.0],
                    "case_prefix": "e2e",
                },
                "conditions": {
                    "gas_temperature_K": 300.0,
                    "pressure_Pa": 13.3,
                    "gas_mixture": [
                        {"species": "Ar", "fraction": 1.0, "mass_amu": 39.948},
                        {"species": "O2", "fraction": 0.0, "mass_amu": 31.998},
                    ],
                },
                "cross_sections": {
                    "format": "csv",
                    "high_energy_extrapolation": "zero",
                    "files": [{"path": xs_path.as_posix(), "format": "csv"}],
                },
                "physics": {
                    "field": {
                        "type": "dc",
                        "magnetic_field": {
                            "enabled": False,
                            "B_T": 0.0,
                            "angle_EB_deg": 0.0,
                        },
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
                        "tail_rate_warning_fraction": 0.05,
                        "max_eV_limit": 1000.0,
                    },
                },
                "solvers": {"two_term": {"backend": "native_sg"}},
                "comparison": {"enabled": False},
                "feature_policy": {"unsupported": "fail", "degraded": "record"},
                "output": {"directory": (root / "outputs").as_posix(), "base_name": "e2e"},
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    workflow_path = root / "workflow.yaml"
    workflow_path.write_text(
        yaml.safe_dump(
            {
                "base_config": base_path.name,
                "database": "outputs/swarm.sqlite",
                "e_over_n_Td": [15.0, 20.0],
                "mixtures": [
                    {"Ar": 1.0, "O2": 0.0},
                    {"Ar": 0.99, "O2": 0.01},
                ],
                "mc": {"replicas": 1, "base_seed": 12345},
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    return workflow_path


def _write_mapping(root: Path, bundle: Path) -> Path:
    model_dir = root / "model" / "work"
    model_dir.mkdir(parents=True)
    input_mph = model_dir / "positive_column_1d.mph"
    input_mph.write_bytes(b"fake mph")
    mapping_path = root / "mapping.yaml"
    mapping_path.write_text(
        f"""
model:
  input_mph: {input_mph.as_posix()}
  output_mph: {(model_dir / "positive_column_1d_swarm_tables.mph").as_posix()}
  study: std1
  component: comp1
  physics: plas
bundle:
  path: {bundle.as_posix()}
logs:
  path: {(root / "model" / "logs").as_posix()}
verify:
  output_path: {(root / "outputs" / "comsol_verify").as_posix()}
functions:
  mean_energy_vs_en:
    tag: sw_meanE
    file: mean_energy_vs_en.csv
    column: mean_energy_eV
    nargs: 1
    argunit: V*m^2
    fununit: eV
    interp: linear
    extrap: const
  muN:
    tag: sw_muN
    file: transport_vs_en.csv
    column: reduced_mobility_m2_V_s_m3
    nargs: 1
    argunit: V*m^2
    fununit: 1/(V*m*s)
    interp: linear
    extrap: const
closure:
  feature: pes1
  mean_energy_formulation:
    mode: local_energy
    property_group: ElectronProperties
    property: MeanElectronEnergyModel
    comsol_value: LocalEnergyApproximationE
  mean_energy:
    table: mean_energy_vs_en.csv
reaction_lookups:
  excitation:
    feature: eir2
    form: townsend
    process_type: excitation
  ionization:
    feature: eir4
    form: townsend
    process_type: ionization
run:
  voltage_feature: mct1
  voltages_V: [20, 50, 100, 200]
results:
  probes:
    - name: x
      expression: x
      unit: m
    - name: electron_density
      expression: plas.ne
      unit: 1/m^3
    - name: mean_electron_energy
      expression: plas.meanE
      unit: eV
    - name: electric_potential
      expression: V
      unit: V
    - name: electron_current_density
      expression: plas.Jex
      unit: A/m^2
    - name: ion_current_density
      expression: plas.Jix
      unit: A/m^2
    - name: total_current_density
      expression: plas.Jex+plas.Jix
      unit: A/m^2
    - name: excitation_source
      expression: plas.Rexc
      unit: 1/(m^3*s)
    - name: ionization_source
      expression: plas.Ri
      unit: 1/(m^3*s)
    - name: E_over_N
      expression: plas.EN
      unit: V*m^2
    - name: reduced_mobility
      expression: plas.muexx*plas.Nn
      unit: 1/(V*m*s)
    - name: reduced_diffusion_L
      expression: plas.Dexx*plas.Nn
      unit: 1/(m*s)
    - name: excitation_townsend
      expression: plas.eir2.alpha(plas.ebar)
      unit: m^2
    - name: ionization_townsend
      expression: plas.eir4.alpha(plas.ebar)
      unit: m^2
    - name: applied_voltage
      expression: V0
      unit: V
    - name: gas_pressure
      expression: p0
      unit: Pa
""".lstrip(),
        encoding="utf-8",
    )
    return mapping_path


def _write_cross_sections(path: Path) -> None:
    rows = [
        ["species", "process", "type", "threshold_eV", "mass_amu", "energy_eV", "cross_section_m2"],
    ]
    for species, mass, scale in [("Ar", 39.948, 1.0), ("O2", 31.998, 1.4)]:
        rows.extend(
            [
                [species, f"{species} momentum", "momentum", "0.0", str(mass), "0.0", f"{scale * 1.0e-20:.6e}"],
                [species, f"{species} momentum", "momentum", "0.0", str(mass), "1.0", f"{scale * 1.2e-20:.6e}"],
                [species, f"{species} momentum", "momentum", "0.0", str(mass), "10.0", f"{scale * 1.4e-20:.6e}"],
                [species, f"{species} momentum", "momentum", "0.0", str(mass), "100.0", f"{scale * 1.1e-20:.6e}"],
                [species, f"{species} excitation", "excitation", "8.0", str(mass), "0.0", "0.0"],
                [species, f"{species} excitation", "excitation", "8.0", str(mass), "20.0", f"{scale * 2.0e-21:.6e}"],
                [species, f"{species} excitation", "excitation", "8.0", str(mass), "100.0", f"{scale * 1.0e-21:.6e}"],
                [species, f"{species} ionization", "ionization", "15.0", str(mass), "0.0", "0.0"],
                [species, f"{species} ionization", "ionization", "15.0", str(mass), "30.0", f"{scale * 1.0e-21:.6e}"],
                [species, f"{species} ionization", "ionization", "15.0", str(mass), "100.0", f"{scale * 0.8e-21:.6e}"],
            ]
        )
    with path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.writer(fp)
        writer.writerows(rows)
