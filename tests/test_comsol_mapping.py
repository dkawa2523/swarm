from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from swarm_workflow.comsol_mapping import (
    ComsolMappingError,
    load_comsol_mapping,
    validate_comsol_mapping_files,
)


def test_mapping_parses_fixed_closure_and_voltage_continuation(tmp_path: Path) -> None:
    mapping_path = _write_mapping_repo(tmp_path)

    mapping = load_comsol_mapping(mapping_path, validate_files=True)

    assert mapping.closure.mean_energy.table == "mean_energy_vs_en.csv"
    assert mapping.run.voltage_feature == "mct1"
    assert mapping.run.voltages_V == (20.0, 50.0, 100.0, 200.0)
    assert {item.process_type for item in mapping.reaction_lookups} == {
        "excitation",
        "ionization",
    }


@pytest.mark.parametrize("legacy_key", ["plasma_closure", "parameters"])
def test_mapping_rejects_removed_generic_sections(
    tmp_path: Path,
    legacy_key: str,
) -> None:
    mapping_path = _write_mapping_repo(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8") + f"\n{legacy_key}: {{}}\n",
        encoding="utf-8",
    )

    with pytest.raises(ComsolMappingError, match="unsupported COMSOL mapping"):
        load_comsol_mapping(mapping_path)


def test_mapping_rejects_nonincreasing_voltage_continuation(tmp_path: Path) -> None:
    mapping_path = _write_mapping_repo(tmp_path)
    text = mapping_path.read_text(encoding="utf-8").replace(
        "voltages_V: [20, 50, 100, 200]",
        "voltages_V: [20, 100, 50, 200]",
    )
    mapping_path.write_text(text, encoding="utf-8")

    with pytest.raises(ComsolMappingError, match="strictly increasing"):
        load_comsol_mapping(mapping_path)


def test_mapping_validation_requires_mean_energy_columns(tmp_path: Path) -> None:
    mapping_path = _write_mapping_repo(tmp_path)
    table = tmp_path / "outputs" / "bundle" / "mean_energy_vs_en.csv"
    table.write_text("E_over_N_V_m2\n1e-20\n", encoding="utf-8")
    mapping = load_comsol_mapping(mapping_path)

    with pytest.raises(ComsolMappingError, match="mean_energy_eV"):
        validate_comsol_mapping_files(mapping)


def _write_mapping_repo(root: Path) -> Path:
    (root / "pyproject.toml").write_text("[project]\nname='tmp'\n", encoding="utf-8")
    (root / "model" / "maps").mkdir(parents=True)
    (root / "model" / "base.mph").write_bytes(b"mph")
    bundle = root / "outputs" / "bundle"
    bundle.mkdir(parents=True)
    tables = {
        "mean_energy_vs_en.csv": (
            ["E_over_N_Td", "E_over_N_V_m2", "mean_energy_eV"],
            [["10", "1e-20", "2"], ["100", "1e-19", "8"]],
        ),
        "transport_vs_en.csv": (
            [
                "E_over_N_Td",
                "E_over_N_V_m2",
                "mean_energy_eV",
                "reduced_mobility_m2_V_s_m3",
                "reduced_diffusion_L_m2_s_m3",
                "reduced_electron_energy_mobility_m2_V_s_m3",
                "reduced_electron_energy_diffusion_m2_s_m3",
            ],
            [
                ["10", "1e-20", "2", "1e24", "2e24", "3e24", "4e24"],
                ["100", "1e-19", "8", "2e24", "3e24", "4e24", "5e24"],
            ],
        ),
        "rates_vs_mean_energy.csv": (
            ["mean_energy_eV", "process_type", "reduced_townsend_m2"],
            [
                ["2", "excitation", "1e-22"],
                ["8", "excitation", "2e-22"],
                ["2", "ionization", "1e-23"],
                ["8", "ionization", "2e-23"],
            ],
        ),
    }
    tables["transport_vs_mean_energy.csv"] = tables["transport_vs_en.csv"]
    metadata = {}
    for name, (columns, rows) in tables.items():
        with (bundle / name).open("w", encoding="utf-8", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(columns)
            writer.writerows(rows)
        metadata[name] = {"columns": columns}
    (bundle / "manifest.json").write_text(
        json.dumps({"tables": metadata}),
        encoding="utf-8",
    )
    mapping = root / "model" / "maps" / "positive_column.yaml"
    mapping.write_text(
        """
model:
  input_mph: model/base.mph
  output_mph: model/work/external.mph
  study: std1
  component: comp1
  physics: plas
bundle:
  path: outputs/bundle
functions:
  mean_energy:
    tag: sw_meanE
    file: mean_energy_vs_en.csv
    column: mean_energy_eV
    nargs: 1
    argunit: V*m^2
    fununit: eV
    interp: linear
    extrap: const
  mobility:
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
""".lstrip(),
        encoding="utf-8",
    )
    return mapping
