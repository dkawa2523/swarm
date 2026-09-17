from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from swarm_workflow.comsol.models.gec_icp.comparison_export import (
    EXPECTED_EXPORT_NAMES,
    FIELD_EXPRESSIONS,
    GecIcpComparisonExportError,
    VOLUME_EXPRESSIONS,
    execute_gec_icp_comparison_export,
    generate_comparison_export_java,
    prepare_gec_icp_comparison_export,
)
from swarm_workflow.comsol.models.gec_icp.contracts import (
    GecIcpMapping,
    GecIcpModelSpec,
    GecIcpParameterTags,
    GecIcpReactionTags,
    GecIcpSpeciesTags,
)


def _model(input_mph: Path, output_mph: Path) -> GecIcpModelSpec:
    return GecIcpModelSpec(
        input_mph=input_mph,
        output_mph=output_mph,
        component="comp1",
        geometry="geom1",
        mesh="mesh1",
        geometry_dimension=2,
        axisymmetric=True,
        plasma_physics="plas",
        plasma_feature="pes1",
        magnetic_physics="mf",
        coil_feature="coil1",
        plasma_conductivity_coupling="pcc1",
        electron_heat_source_coupling="ehs1",
        study="std1",
        study_feature="ftrans",
        solution="sol1",
        dataset="dset1",
    )


def _mapping(tmp_path: Path) -> GecIcpMapping:
    mapping_path = tmp_path / "argon_gec_icp_model.yaml"
    mapping_path.write_text("schema_version: 2\n", encoding="utf-8")
    source = tmp_path / "source.mph"
    return GecIcpMapping(
        path=mapping_path,
        root=tmp_path,
        model=_model(source, tmp_path / "work/result.mph"),
        parameters=GecIcpParameterTags(
            power="Psp", gas_temperature="T0", pressure="p0"
        ),
        species=GecIcpSpeciesTags(ground="Ar", excited="Ars", ion="Ar_1p"),
        reactions=GecIcpReactionTags(
            elastic="eir1",
            excitation="eir2",
            superelastic="eir3",
            ionization="eir4",
            stepwise_ionization="eir5",
        ),
    )


def test_comparison_java_is_read_only_and_selects_exact_common_time(
    tmp_path: Path,
) -> None:
    source = tmp_path / "saved.mph"
    source.write_bytes(b"saved solution")
    java = generate_comparison_export_java(
        _model(source, tmp_path / "unused.mph"),
        input_mph=source,
        output_directory=tmp_path / "comparison",
        common_time_s=1.0e-3,
    )

    assert "ModelUtil.loadCopy" in java
    assert 'model.sol("sol1").getPVals()' in java
    assert "SWARM_GEC_ICP_SOLUTION_TIME" in java
    assert "PrintWriter" not in java
    assert java.count('.set("innerinput", "interp")') == 3
    assert java.count('.set("innerinput", "last")') == 3
    assert java.count('.set("innerinput", "all")') == 2
    assert java.count('.set("t", new double[]{1.00000000000000002e-03})') == 3
    assert "study.run(" not in java
    assert "clearSolutionData(" not in java
    assert "model.save(" not in java
    for expression in (*FIELD_EXPRESSIONS, *VOLUME_EXPRESSIONS, "mf.PCoil_1"):
        assert expression in java


def test_prepare_comparison_export_records_immutable_inputs(tmp_path: Path) -> None:
    mapping = _mapping(tmp_path)
    source = tmp_path / "saved.mph"
    source.write_bytes(b"saved solution")
    output = tmp_path / "comparison" / "original"

    plan = prepare_gec_icp_comparison_export(
        mapping,
        case_id="original",
        input_mph=source,
        output_directory=output,
    )

    assert plan.input_mph_sha256 == hashlib.sha256(b"saved solution").hexdigest()
    assert plan.java_path.is_file()
    assert tuple(path.name for path in plan.expected_outputs) == EXPECTED_EXPORT_NAMES
    manifest = json.loads(plan.manifest_path.read_text(encoding="utf-8"))
    assert manifest["read_only"] is True
    assert manifest["input_mph"]["sha256"] == plan.input_mph_sha256
    assert manifest["time_selection"]["common"] == {
        "mode": "exact_transient_interpolation",
        "time_s": 1.0e-3,
    }
    assert manifest["time_selection"]["terminal"] == {"mode": "last_saved_solution"}
    assert "model.save" in manifest["mutation_contract"]
    assert manifest["expected_outputs"] == [
        str(output.resolve() / name) for name in EXPECTED_EXPORT_NAMES
    ]


@pytest.mark.parametrize(
    ("case_id", "common_time_s", "message"),
    (
        ("Original", 1.0e-3, "case_id"),
        ("original", float("nan"), "finite and nonnegative"),
        ("original", -1.0, "finite and nonnegative"),
    ),
)
def test_prepare_comparison_export_rejects_invalid_requests(
    tmp_path: Path,
    case_id: str,
    common_time_s: float,
    message: str,
) -> None:
    mapping = _mapping(tmp_path)
    source = tmp_path / "saved.mph"
    source.write_bytes(b"saved solution")

    with pytest.raises(GecIcpComparisonExportError, match=message):
        prepare_gec_icp_comparison_export(
            mapping,
            case_id=case_id,
            input_mph=source,
            output_directory=tmp_path / "comparison",
            common_time_s=common_time_s,
        )


def test_prepare_comparison_export_rejects_missing_mph(tmp_path: Path) -> None:
    mapping = _mapping(tmp_path)

    with pytest.raises(GecIcpComparisonExportError, match="existing MPH"):
        prepare_gec_icp_comparison_export(
            mapping,
            case_id="original",
            input_mph=tmp_path / "missing.mph",
            output_directory=tmp_path / "comparison",
        )


def test_execute_comparison_export_verifies_fresh_outputs_and_source_hash(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping = _mapping(tmp_path)
    source = tmp_path / "saved.mph"
    source.write_bytes(b"saved solution")
    plan = prepare_gec_icp_comparison_export(
        mapping,
        case_id="original",
        input_mph=source,
        output_directory=tmp_path / "comparison" / "original",
    )
    plan.expected_outputs[0].write_text("stale", encoding="utf-8")

    def fake_execute(context, java_path, **kwargs):
        assert not plan.expected_outputs[0].exists()
        assert context.input_mph == source.resolve()
        assert context.output_mph == source.resolve()
        assert kwargs["require_fresh_output"] is False
        for path in plan.expected_outputs:
            if path.name != "solution_times.csv":
                path.write_text("fresh\n", encoding="utf-8")
        stdout = tmp_path / "batch_stdout.txt"
        stdout.write_text(
            "COMSOL noise\n"
            "SWARM_GEC_ICP_SOLUTION_TIME\t1\t0.0\n"
            "SWARM_GEC_ICP_SOLUTION_TIME\t2\t0.001\n",
            encoding="utf-8",
        )
        return SimpleNamespace(
            operation=kwargs["operation"],
            log_dir=tmp_path / "logs",
            command_json=tmp_path / "command.json",
            result_json=tmp_path / "result.json",
            provenance_json=tmp_path / "provenance.json",
            stdout_paths=(tmp_path / "compile_stdout.txt", stdout),
            return_codes=(0, 0),
            total_time_s=1.25,
            comsol_version="6.4",
            comsol_build="test",
        )

    monkeypatch.setattr(
        "swarm_workflow.comsol.models.gec_icp.comparison_export."
        "execute_generated_comsol_java",
        fake_execute,
    )
    execute_gec_icp_comparison_export(mapping, plan)

    assert source.read_bytes() == b"saved solution"
    manifest = json.loads(plan.manifest_path.read_text(encoding="utf-8"))
    assert manifest["status"] == "completed"
    assert (
        manifest["input_mph"]["sha256"]
        == manifest["input_mph"]["post_execution_sha256"]
    )
    assert len(manifest["outputs"]) == len(plan.expected_outputs)


def test_execute_comparison_export_rejects_changed_source(tmp_path: Path) -> None:
    mapping = _mapping(tmp_path)
    source = tmp_path / "saved.mph"
    source.write_bytes(b"saved solution")
    plan = prepare_gec_icp_comparison_export(
        mapping,
        case_id="original",
        input_mph=source,
        output_directory=tmp_path / "comparison",
    )
    source.write_bytes(b"changed")

    with pytest.raises(GecIcpComparisonExportError, match="changed after preparation"):
        execute_gec_icp_comparison_export(mapping, plan)
