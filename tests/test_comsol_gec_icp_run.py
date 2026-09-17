from __future__ import annotations

import json
from pathlib import Path
import re
from types import SimpleNamespace
from typing import Any

import pytest
import yaml

from swarm_workflow.comsol.models.gec_icp.convergence import (
    assess_gec_icp_convergence,
)
from swarm_workflow.comsol.models.gec_icp.java import (
    TRANSPORT_TABLE,
    generate_apply_java,
    generate_solve_java,
    output_times,
)
from swarm_workflow.comsol.models.gec_icp.run_contracts import (
    GecIcpWorkflowError,
)
from swarm_workflow.comsol.models.gec_icp.run_mapping import (
    load_gec_icp_run_mapping,
)


_SOURCE_EEDF = {
    "two_term": (
        "eedf_f0_comsol_2d.csv",
        "structured_spreadsheet_linear_projection",
        None,
    ),
    "monte_carlo": (
        "eedf_f0_comsol_2d.csv",
        "structured_spreadsheet_linear_projection",
        "function_eedf_restricted_lmea",
    ),
    "propagator": (
        "eedf_f0_comsol_2d.csv",
        "structured_spreadsheet_linear_projection",
        None,
    ),
}


def _model_mapping_payload() -> dict[str, Any]:
    return {
        "schema_version": 2,
        "model": {
            "input_mph": "input.mph",
            "output_mph": "unused-contract-output.mph",
            "component": "comp1",
            "geometry": "geom1",
            "mesh": "mesh1",
            "geometry_dimension": 2,
            "axisymmetric": True,
            "plasma_physics": "plas",
            "plasma_feature": "pes1",
            "magnetic_physics": "mf",
            "coil_feature": "coil1",
            "plasma_conductivity_coupling": "pcc1",
            "electron_heat_source_coupling": "ehs1",
            "study": "std1",
            "study_feature": "ftrans",
            "solution": "sol1",
            "dataset": "dset1",
        },
        "parameters": {
            "power": "Psp",
            "gas_temperature": "T0",
            "pressure": "p0",
        },
        "species": {"ground": "Ar", "excited": "Ars", "ion": "Ar_1p"},
        "reactions": {
            "elastic": "eir1",
            "excitation": "eir2",
            "superelastic": "eir3",
            "ionization": "eir4",
            "stepwise_ionization": "eir5",
        },
    }


def _run_mapping_payload(source: str) -> dict[str, Any]:
    table, interpolation, mc_profile = _SOURCE_EEDF[source]
    bundle: dict[str, Any] = {
        "path": "bundle",
        "expected_source": source,
        "expected_field_type": "dc",
        "expected_transport_definition": "flux",
    }
    if mc_profile is not None:
        bundle["expected_mc_qualification_profile"] = mc_profile
    return {
        "schema_version": 2,
        "model_mapping": "model.yaml",
        "bundle": bundle,
        "closure": {
            "source_field": "steady_dc",
            "electron_transport": "swarm_mobility_comsol_einstein",
            "reaction_model": "function_eedf",
            "elastic_energy_loss_model": "comsol_cross_section_integral",
            "chemistry_owner": "comsol_embedded_cross_sections",
            "rf_owner": "comsol_frequency_transient",
            "function_eedf": {
                "table": table,
                "function_tag": f"sw_icp_eedf_{source}",
                "interpolation": interpolation,
                "extrapolation": "constant",
            },
        },
        "run": {
            "output_mph": "work/result.mph",
            "power_W": 1500.0,
            "frequency_Hz": 1.356e7,
            "gas_temperature_K": 300.0,
            "pressure_Pa": 2.66644,
            "final_time_s": 0.01,
            "output_points_per_decade": 4,
            "clear_saved_solution": True,
            "convergence": {
                "electron_inventory_relative_change": 0.02,
                "ion_inventory_relative_change": 0.02,
                "metastable_inventory_relative_change": 0.03,
                "mean_energy_relative_change": 0.01,
                "absorbed_power_relative_change": 0.01,
                "coil_power_relative_error": 0.005,
                "mean_energy_support_relative_guard": 0.10,
            },
        },
        "results": {"path": "results"},
        "logs": {"path": "logs"},
    }


def _write_run_mapping(
    tmp_path: Path,
    payload: dict[str, Any],
) -> Path:
    root = tmp_path / "maps"
    root.mkdir(parents=True, exist_ok=True)
    (root / "input.mph").touch()
    (root / "model.yaml").write_text(
        yaml.safe_dump(_model_mapping_payload(), sort_keys=False),
        encoding="utf-8",
    )
    path = root / "run.yaml"
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    return path


@pytest.mark.parametrize(
    ("source", "expected_table", "expected_interpolation", "expected_profile"),
    [
        (
            "two_term",
            "eedf_f0_comsol_2d.csv",
            "structured_spreadsheet_linear_projection",
            None,
        ),
        (
            "monte_carlo",
            "eedf_f0_comsol_2d.csv",
            "structured_spreadsheet_linear_projection",
            "function_eedf_restricted_lmea",
        ),
        (
            "propagator",
            "eedf_f0_comsol_2d.csv",
            "structured_spreadsheet_linear_projection",
            None,
        ),
    ],
)
def test_run_mapping_loads_the_common_function_eedf_contract(
    tmp_path: Path,
    source: str,
    expected_table: str,
    expected_interpolation: str,
    expected_profile: str | None,
) -> None:
    mapping = load_gec_icp_run_mapping(
        _write_run_mapping(tmp_path, _run_mapping_payload(source))
    )

    assert mapping.bundle.expected_source == source
    assert mapping.bundle.expected_mc_qualification_profile == expected_profile
    assert mapping.closure.function_eedf.table == expected_table
    assert mapping.closure.function_eedf.interpolation == expected_interpolation
    assert mapping.run.output_mph == (tmp_path / "maps/work/result.mph").resolve()
    assert mapping.output_directory == (tmp_path / "maps/results").resolve()


@pytest.mark.parametrize(
    ("source", "field", "wrong_value"),
    [
        ("two_term", "table", "eedf_f0_comsol_grid.txt"),
        (
            "two_term",
            "interpolation",
            "structured_grid_linear_projection",
        ),
        ("monte_carlo", "table", "eedf_f0_comsol_grid.txt"),
        (
            "monte_carlo",
            "interpolation",
            "structured_grid_linear_projection",
        ),
        ("propagator", "table", "eedf_f0_comsol_grid.txt"),
        (
            "propagator",
            "interpolation",
            "structured_grid_linear_projection",
        ),
    ],
)
def test_run_mapping_rejects_noncanonical_function_eedf_contract(
    tmp_path: Path,
    source: str,
    field: str,
    wrong_value: str,
) -> None:
    payload = _run_mapping_payload(source)
    function = payload["closure"]["function_eedf"]
    assert isinstance(function, dict)
    function[field] = wrong_value

    with pytest.raises(
        GecIcpWorkflowError,
        match="Function-EEDF must use",
    ):
        load_gec_icp_run_mapping(_write_run_mapping(tmp_path, payload))


@pytest.mark.parametrize(
    ("source", "profile", "message"),
    [
        (
            "monte_carlo",
            None,
            "Monte Carlo ICP mapping requires the "
            "function_eedf_restricted_lmea qualification profile",
        ),
        (
            "two_term",
            "function_eedf_restricted_lmea",
            "expected_mc_qualification_profile is only valid for "
            "monte_carlo or composite",
        ),
    ],
)
def test_run_mapping_enforces_source_specific_mc_qualification(
    tmp_path: Path,
    source: str,
    profile: str | None,
    message: str,
) -> None:
    payload = _run_mapping_payload(source)
    bundle = payload["bundle"]
    assert isinstance(bundle, dict)
    if profile is None:
        bundle.pop("expected_mc_qualification_profile", None)
    else:
        bundle["expected_mc_qualification_profile"] = profile

    with pytest.raises(GecIcpWorkflowError, match=message):
        load_gec_icp_run_mapping(_write_run_mapping(tmp_path, payload))


@pytest.mark.parametrize(
    ("owner", "field", "value"),
    [
        ("run", "power_W", float("nan")),
        ("run", "frequency_Hz", float("inf")),
        ("convergence", "mean_energy_relative_change", float("nan")),
        ("convergence", "coil_power_relative_error", float("inf")),
    ],
)
def test_run_mapping_rejects_nonfinite_numeric_values(
    tmp_path: Path,
    owner: str,
    field: str,
    value: float,
) -> None:
    payload = _run_mapping_payload("two_term")
    run = payload["run"]
    assert isinstance(run, dict)
    target = run if owner == "run" else run["convergence"]
    assert isinstance(target, dict)
    target[field] = value

    with pytest.raises(GecIcpWorkflowError, match="must be a positive number"):
        load_gec_icp_run_mapping(_write_run_mapping(tmp_path, payload))


def _loaded_mapping(tmp_path: Path, source: str = "two_term") -> Any:
    mapping = load_gec_icp_run_mapping(
        _write_run_mapping(tmp_path, _run_mapping_payload(source))
    )
    mapping.bundle.path.mkdir(parents=True, exist_ok=True)
    (mapping.bundle.path / TRANSPORT_TABLE).write_text(
        "mean_energy_eV,reduced_mobility_m2_V_s_m3\n1.0,1.0e24\n10.0,2.0e24\n",
        encoding="utf-8",
    )
    return mapping


def test_generated_java_applies_restricted_closure_and_one_cold_solve(
    tmp_path: Path,
) -> None:
    mapping = _loaded_mapping(tmp_path)

    apply_java = generate_apply_java(mapping)
    solve_java = generate_solve_java(mapping)

    assert apply_java.count('"SpecifyMueOnly"') == 1
    for reaction in ("eir1", "eir2", "eir3", "eir4", "eir5"):
        target = f'model.component("comp1").physics("plas").feature("{reaction}")'
        assert f'{target}.set("eedf", "FromPhysicsInterfaceProperty");' in apply_java
    assert solve_java.count('model.study("std1").run();') == 1
    assert 'model.sol("sol1").clearSolutionData();' in solve_java
    assert solve_java.index('model.sol("sol1").clearSolutionData();') < (
        solve_java.index('model.study("std1").run();')
    )
    assert "Parametric" not in solve_java
    assert "plistarr" not in solve_java
    assert 'model.param().set("Psp", "1500[W]");' in solve_java

    match = re.search(r'\.set\("tlist", "([^"]+)"\);', solve_java)
    assert match is not None
    times = [float(value) for value in match.group(1).split()]
    assert len(times) == 26
    assert times[:2] == pytest.approx([0.0, 1.0e-8])
    assert times[-1] == pytest.approx(0.01)
    assert all(right > left for left, right in zip(times, times[1:]))
    assert solve_java.count("model.result().export().create(") == 4
    for name in (
        "convergence_volume.csv",
        "convergence_mean_energy_min.csv",
        "convergence_mean_energy_max.csv",
        "convergence_coil_power.csv",
    ):
        assert name in solve_java


def test_generated_java_exports_self_describing_full_precision_csv(
    tmp_path: Path,
) -> None:
    solve_java = generate_solve_java(_loaded_mapping(tmp_path))

    assert (
        '.set("expr", new String[]{"1", "plas.ne", '
        '"plas.ne*e_const*plas.ebar", "plas.n_wArs", '
        '"plas.n_wAr_1p", "mf.Qrh"});'
    ) in solve_java
    assert '.set("unit", new String[]{"m^3", "1", "eV", "1", "1", "W"});' in solve_java
    assert solve_java.count('.set("expr", new String[]{"e_const*plas.ebar"});') == 2
    assert (
        '.set("descr", new String[]{"axisymmetric plasma volume", '
        '"electron inventory", "electron mean-energy inventory", '
        '"argon metastable inventory", "argon ion inventory", '
        '"plasma absorbed RF power"});'
    ) in solve_java
    assert '.set("descr", new String[]{"minimum electron mean energy"});' in solve_java
    assert '.set("descr", new String[]{"maximum electron mean energy"});' in solve_java
    assert (
        '.set("descr", new String[]{"configured coil-power readback"});' in solve_java
    )
    for tag in (
        "swIcpVolumeExport",
        "swIcpMinEnergyExport",
        "swIcpMaxEnergyExport",
        "swIcpCoilPowerExport",
    ):
        owner = f'model.result().export("{tag}")'
        assert f'{owner}.set("header", "on");' in solve_java
        assert f'{owner}.set("prec", "full");' in solve_java
        assert f'{owner}.set("ifexists", "overwrite");' in solve_java
        assert f'{owner}.set("separator"' not in solve_java


def test_generated_monte_carlo_java_uses_the_common_spreadsheet_eedf(
    tmp_path: Path,
) -> None:
    apply_java = generate_apply_java(_loaded_mapping(tmp_path, "monte_carlo"))

    assert "eedf_f0_comsol_2d.csv" in apply_java
    assert "eedf_f0_binding_seed.csv" not in apply_java
    assert "eedf_f0_comsol_grid.txt" not in apply_java
    assert '.set("struct", "spreadsheet")' in apply_java
    assert '.set("struct", "grid")' not in apply_java
    assert "discardData" not in apply_java
    assert (
        apply_java.count('.func().create("sw_icp_eedf_monte_carlo", "Interpolation")')
        == 1
    )
    assert apply_java.count('.set("nargs", 2)') == 1
    assert apply_java.count('.func("sw_icp_eedf_monte_carlo").importData();') == 1
    assert '.set("argunit", "eV,eV")' in apply_java


def _write_numeric_export(path: Path, rows: list[tuple[float, ...]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "% synthetic COMSOL export\n"
        + "\n".join(",".join(f"{value:.17g}" for value in row) for row in rows)
        + "\n",
        encoding="utf-8",
    )


def _convergence_plan(tmp_path: Path, *, passing: bool) -> Any:
    mapping = _loaded_mapping(tmp_path)
    output = mapping.output_directory
    times = output_times(mapping)
    if passing:
        late_electrons = (100.0, 100.5, 100.8, 101.0)
        late_metastables = (50.0, 50.1, 50.3, 50.5)
        late_ions = (100.0, 100.4, 100.7, 101.0)
        late_absorbed_power = (1475.0, 1480.0, 1485.0, 1490.0)
        late_coil_power = (1490.0, 1492.0, 1494.0, 1495.0)
        maximum_energy = 4.0
    else:
        late_electrons = (100.0, 105.0, 115.0, 130.0)
        late_metastables = (50.0, 50.1, 50.3, 50.5)
        late_ions = (100.0, 100.4, 100.7, 101.0)
        late_absorbed_power = (1490.0, 1495.0, 1498.0, 1500.0)
        late_coil_power = (1490.0, 1450.0, 1350.0, 1200.0)
        maximum_energy = 12.0
    volume_rows = [
        (time, 1.0, 100.0, 300.0, 50.0, 100.0, 1490.0) for time in times[:-4]
    ]
    volume_rows.extend(
        (
            time,
            1.0,
            electrons,
            3.0 * electrons,
            metastables,
            ions,
            absorbed_power,
        )
        for time, electrons, metastables, ions, absorbed_power in zip(
            times[-4:],
            late_electrons,
            late_metastables,
            late_ions,
            late_absorbed_power,
            strict=True,
        )
    )
    _write_numeric_export(
        output / "convergence_volume.csv",
        volume_rows,
    )
    _write_numeric_export(
        output / "convergence_mean_energy_min.csv",
        [(time, 2.5) for time in times],
    )
    _write_numeric_export(
        output / "convergence_mean_energy_max.csv",
        [(time, maximum_energy) for time in times],
    )
    _write_numeric_export(
        output / "convergence_coil_power.csv",
        [
            *((time, 1490.0) for time in times[:-4]),
            *zip(times[-4:], late_coil_power, strict=True),
        ],
    )
    return SimpleNamespace(mapping=mapping, output_directory=output)


def test_convergence_accepts_stable_supported_synthetic_exports(
    tmp_path: Path,
) -> None:
    plan = _convergence_plan(tmp_path, passing=True)

    result = assess_gec_icp_convergence(plan)

    assert result["schema"] == "swarm.gec_icp_convergence.v1"
    assert result["status"] == "passed"
    assert result["passed"] is True
    assert all(gate["passed"] for gate in result["gates"].values())
    written = json.loads(
        (plan.output_directory / "convergence.json").read_text(encoding="utf-8")
    )
    assert written == result


def test_convergence_rejects_unstable_or_unsupported_synthetic_exports(
    tmp_path: Path,
) -> None:
    result = assess_gec_icp_convergence(_convergence_plan(tmp_path, passing=False))

    assert result["status"] == "failed"
    assert result["passed"] is False
    assert result["gates"]["electron_inventory_relative_change"]["passed"] is False
    assert result["gates"]["coil_power_relative_error"]["passed"] is False
    assert result["gates"]["mean_energy_support"]["passed"] is False


def test_convergence_requires_guarded_upper_mean_energy_support(
    tmp_path: Path,
) -> None:
    plan = _convergence_plan(tmp_path, passing=True)
    times = output_times(plan.mapping)
    _write_numeric_export(
        plan.output_directory / "convergence_mean_energy_max.csv",
        [(time, 9.5) for time in times],
    )

    result = assess_gec_icp_convergence(plan)

    gate = result["gates"]["mean_energy_support"]
    assert result["passed"] is False
    assert gate["passed"] is False
    assert gate["relative_guard"] == pytest.approx(0.10)
    assert gate["required_upper_eV"] == pytest.approx(10.45)


def test_convergence_rejects_an_unplanned_but_shared_time_axis(
    tmp_path: Path,
) -> None:
    plan = _convergence_plan(tmp_path, passing=True)
    filenames = (
        "convergence_volume.csv",
        "convergence_mean_energy_min.csv",
        "convergence_mean_energy_max.csv",
        "convergence_coil_power.csv",
    )
    for filename in filenames:
        path = plan.output_directory / filename
        rows = [
            tuple(float(value) for value in line.split(","))
            for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("%")
        ]
        rows[-2] = ((rows[-3][0] + rows[-1][0]) / 2.0, *rows[-2][1:])
        _write_numeric_export(path, rows)

    with pytest.raises(
        GecIcpWorkflowError,
        match="does not match the planned time axis",
    ):
        assess_gec_icp_convergence(plan)


def test_convergence_rejects_a_csv_with_an_extra_column(
    tmp_path: Path,
) -> None:
    plan = _convergence_plan(tmp_path, passing=True)
    path = plan.output_directory / "convergence_volume.csv"
    lines = path.read_text(encoding="utf-8").splitlines()
    lines[-1] += ",123"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    with pytest.raises(
        GecIcpWorkflowError,
        match="incomplete COMSOL convergence export",
    ):
        assess_gec_icp_convergence(plan)
