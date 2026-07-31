from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path
from zipfile import ZipFile

from swarm_workflow.comsol_gec_ccp import (
    generate_apply_java,
    generate_run_java,
    inspect_gec_ccp_mph,
    prepare_gec_ccp_run,
)


def test_gec_ccp_contract_and_java_preserve_solved_mean_energy(
    tmp_path: Path,
) -> None:
    mapping, bundle = _write_fixture(tmp_path)

    plan = prepare_gec_ccp_run(mapping, bundle_path=bundle)
    source = generate_apply_java(plan.mapping)

    assert plan.contract.original_eedf == "Druyvesteyn"
    assert plan.contract.physics_operation == "ColdPlasmaTimePeriodic"
    assert set(plan.contract.reaction_features) == {"eir1", "eir2", "eir3"}
    assert '"SpecifyElectronDensityAndEnergy", "SpecifyAll"' in source
    assert "ptp.en/max(ptp.ne,1[1/m^3])" in source
    assert "sw_muN_e" in source
    assert "sw_DeN_e" in source
    assert "sw_muenN_e" in source
    assert "sw_DenN_e" in source
    assert "/max(ptp.Nn,1[1/m^3])" in source
    assert "nojac" not in source
    assert '"SourceStabilization", true' in source
    assert '"ReactionSourceStabilization", true' in source
    assert "SpecifyMeanElectronEnergy" not in source
    assert source.count('"RateConstantForm", "UseRate"') == 3
    assert source.count('"SpecifyReactionUsing", "RateConstant"') == 3
    assert source.count(
        "6.02214076e23[1/mol]*exp(sw_logk_"
    ) == 3
    assert source.count('*1[m^3/s]"') == 3
    assert source.count('.set("interp", "piecewisecubic")') == 7
    scale_match = re.search(
        r"shiftPeriodicLogEnergySolution\(model,\s*([0-9.eE+-]+)\);",
        source,
    )
    assert scale_match is not None
    assert math.isclose(
        float(scale_match.group(1)),
        1.216872309669,
        rel_tol=1.0e-14,
        abs_tol=0.0,
    )
    assert '"comp1.En_per".equals(names[nameIndex])' in source
    assert "periodic.createSolution()" in source
    assert '"SpecifyReactionUsing", "UseLookupTable"' not in source
    assert '"SpecifyReactionUsing", "UseCrossSectionData"' not in source
    assert '"xratedata"' not in source
    assert '"yratedata"' not in source
    assert "xtownratedata" not in source
    assert "sw_mobility_blend" not in source
    run_source = generate_run_java(
        plan.mapping,
        class_name="Run",
        input_mph=plan.mapping.model.output_mph,
        output_mph=plan.mapping.model.output_mph,
        robust_nonlinear=True,
    )
    assert 'model.param().set("P0", "1[W]")' in run_source
    assert run_source.count('model.study("std1").run()') == 1
    assert run_source.count('model.study("std2").run()') == 1
    assert "for (" not in run_source
    assert "blend" not in run_source
    assert '.set("dtech", "hnlin")' in run_source
    assert '.set("minsteph", 1.0e-12)' in run_source
    assert '.set("useminsteprecovery", "off")' in run_source
    assert '.set("maxiter", 200)' in run_source
    assert plan.plan_json.exists()
    plan_data = json.loads(plan.plan_json.read_text(encoding="utf-8"))
    assert plan_data["closure"]["coefficient_sweep"] is False
    assert plan_data["closure"]["power_W"] == 1.0
    assert (
        plan_data["closure"]["initial_guess"]["electron_mean_energy_scale"]
        == 1.216872309669
    )


def test_repository_gec_ccp_model_reports_druyvesteyn() -> None:
    model = Path(__file__).parents[1] / "comsol_modes" / "argon_gec_ccp.mph"
    contract = inspect_gec_ccp_mph(model)

    assert contract.original_eedf == "Druyvesteyn"
    assert contract.physics_tag == "ptp"
    assert contract.plasma_feature == "pes1"


def _write_fixture(root: Path) -> tuple[Path, Path]:
    maps = root / "maps"
    maps.mkdir()
    model = root / "argon_gec_ccp.mph"
    xml = """
<Model>
  <Physics op="ColdPlasmaTimePeriodic" tag="ptp"/>
  <PhysicsFeature op="ElectronImpactReaction" tag="eir1">
    <param param="eedf" value="1|1,'FromPhysicsInterfaceProperty'"/>
  </PhysicsFeature>
  <PhysicsFeature op="ElectronImpactReaction" tag="eir2"/>
  <PhysicsFeature op="ElectronImpactReaction" tag="eir3"/>
  <PhysicsFeature op="PlasmaEsModel" tag="pes1"/>
  <PhysicsProp>
    <param param="eedf" value="1|1,'Druyvesteyn'"/>
  </PhysicsProp>
  <Study tag="std1"/>
  <Study tag="std2"/>
  <DatasetFeature tag="dset1"/>
  <DatasetFeature tag="dset2"/>
  <DatasetFeature tag="dset3"/>
  <DatasetFeature tag="cln1"/>
  <DatasetFeature tag="cln2"/>
</Model>
""".strip()
    with ZipFile(model, "w") as archive:
        archive.writestr("dmodel.xml", xml)

    bundle = root / "bundle"
    bundle.mkdir()
    transport_columns = [
        "mean_energy_eV",
        "reduced_mobility_m2_V_s_m3",
        "reduced_diffusion_L_m2_s_m3",
        "reduced_electron_energy_mobility_m2_V_s_m3",
        "reduced_electron_energy_diffusion_m2_s_m3",
    ]
    _write_csv(
        bundle / "transport_vs_mean_energy.csv",
        transport_columns,
        [
            [1, 10, 20, 30, 40],
            [2, 11, 21, 31, 41],
        ],
    )
    rate_columns = [
        "mean_energy_eV",
        "process_type",
        "rate_coefficient_m3_s",
    ]
    _write_csv(
        bundle / "rates_vs_mean_energy.csv",
        rate_columns,
        [
            [1, "elastic", 1e-14],
            [2, "elastic", 2e-14],
            [1, "excitation", 1e-20],
            [2, "excitation", 2e-20],
            [1, "ionization", 1e-22],
            [2, "ionization", 2e-22],
        ],
    )
    eedf_columns = [
        "electron_energy_eV",
        "mean_energy_eV",
        "E_over_N_Td",
        "eepf_eV_m32",
    ]
    _write_csv(
        bundle / "eedf_f0.csv",
        eedf_columns,
        [[0.1, 1, 10, 1], [0.2, 2, 20, 2]],
    )
    quality_columns = [
        "E_over_N_Td",
        "passed",
        "eedf_normalization_error",
    ]
    _write_csv(
        bundle / "quality.csv",
        quality_columns,
        [[10, 1, 0], [20, 1, 1e-16]],
    )
    (bundle / "manifest.json").write_text(
        json.dumps(
            {
                "status": "ok",
                "source": "two_term",
                "valid_ranges": {"mean_energy_eV": [1, 2]},
                "monotonicity": {"mean_energy_strictly_monotonic": True},
                "tables": {
                    "transport_vs_mean_energy.csv": {
                        "columns": transport_columns
                    },
                    "rates_vs_mean_energy.csv": {"columns": rate_columns},
                    "eedf_f0.csv": {"columns": eedf_columns},
                    "quality.csv": {"columns": quality_columns},
                },
            }
        ),
        encoding="utf-8",
    )
    mapping = maps / "gec.yaml"
    mapping.write_text(
        """
schema_version: 2
model:
  input_mph: ../argon_gec_ccp.mph
  baseline_output_mph: ../work/baseline.mph
  external_output_mph: ../work/external.mph
  component: comp1
  physics: ptp
  plasma_feature: pes1
  time_periodic_study: std1
  conversion_study: std2
  expected_original_eedf: Druyvesteyn
bundle:
  path: ../bundle
reactions:
  - {name: elastic, feature: eir1, process_type: elastic}
  - {name: excitation, feature: eir2, process_type: excitation}
  - {name: ionization, feature: eir3, process_type: ionization}
run:
  power_parameter: P0
  power_W: 1.0
  initial_mean_energy_scale: 1.216872309669
  axis_dataset: cln1
  radial_dataset: cln2
  period_dataset: dset1
  phase_dataset: dset3
  waveform_dataset: dset2
results:
  output_directory: ../results
logs:
  path: ../logs
""".lstrip(),
        encoding="utf-8",
    )
    return mapping, bundle


def _write_csv(path: Path, columns: list[str], rows: list[list[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        writer.writerows(rows)
