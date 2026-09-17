from __future__ import annotations

import json
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from swarm_workflow._io import write_json
from swarm_workflow.comsol.input.mean_energy import MeanEnergyArgument
from swarm_workflow.comsol.models.gec_icp import closure_readback as readback
from swarm_workflow.comsol.models.gec_icp import execution
from swarm_workflow.comsol.models.gec_icp.java import generate_apply_java
from swarm_workflow.comsol.models.gec_icp.run_contracts import (
    GecIcpClosureSpec,
    GecIcpFunctionEedfSpec,
    GecIcpQualityError,
)


def _mapping(tmp_path: Path) -> Any:
    bundle = tmp_path / "bundle"
    bundle.mkdir(parents=True, exist_ok=True)
    (bundle / "transport_vs_mean_energy.csv").write_text(
        "mean_energy_eV,reduced_mobility_m2_V_s_m3\n1.0,1.0e24\n10.0,2.0e24\n",
        encoding="utf-8",
    )
    model = SimpleNamespace(
        input_mph=tmp_path / "input.mph",
        component="comp1",
        plasma_physics="plas",
        plasma_feature="pes1",
        dataset="dset1",
        study="std1",
    )
    reactions = SimpleNamespace(
        elastic="eir1",
        excitation="eir2",
        superelastic="eir3",
        ionization="eir4",
        stepwise_ionization="eir5",
    )
    closure = GecIcpClosureSpec(
        source_field="steady_dc",
        electron_transport="swarm_mobility_comsol_einstein",
        reaction_model="function_eedf",
        elastic_energy_loss_model="comsol_cross_section_integral",
        chemistry_owner="comsol_embedded_cross_sections",
        rf_owner="comsol_frequency_transient",
        function_eedf=GecIcpFunctionEedfSpec(
            table="eedf_f0_comsol_2d.csv",
            function_tag="sw_icp_eedf_two_term",
            interpolation="structured_spreadsheet_linear_projection",
            extrapolation="constant",
        ),
    )
    return SimpleNamespace(
        root=tmp_path,
        path=tmp_path / "run.yaml",
        model_mapping=SimpleNamespace(model=model, reactions=reactions),
        bundle=SimpleNamespace(path=bundle, expected_source="two_term"),
        closure=closure,
        run=SimpleNamespace(output_mph=tmp_path / "result.mph"),
        output_directory=tmp_path / "results",
        log_path=tmp_path / "logs",
    )


def _plan(tmp_path: Path) -> Any:
    mapping = _mapping(tmp_path)
    return readback.prepare_gec_icp_closure_readback(
        mapping,
        output_directory=mapping.output_directory,
        write_java=True,
    )


def _write_readback_log(
    plan: Any,
    *,
    overrides: dict[str, str] | None = None,
    omit: set[str] | None = None,
    table_delta: float = 0.0,
    value_delta: float = 0.0,
    duplicate_key: str | None = None,
) -> Path:
    energies = [1.0, 10.0]
    mobilities = [1.0e24, 2.0e24]
    contract = readback._expected_contract(plan, energies)
    contract.update(
        {
            "saved_solution.mean_energy_relative_error_max": "1e-14",
            "saved_solution.support_log_difference_max": "2e-14",
        }
    )
    contract.update(overrides or {})
    for key in omit or set():
        contract.pop(key, None)
    lines = ["COMSOL batch output"]
    lines.extend(
        f"SWARM_ICP_CLOSURE_CONTRACT\t{key}\t{value}"
        for key, value in sorted(contract.items())
    )
    if duplicate_key is not None:
        lines.append(
            f"SWARM_ICP_CLOSURE_CONTRACT\t{duplicate_key}\t{contract[duplicate_key]}"
        )
    for index, (energy, mobility) in enumerate(zip(energies, mobilities, strict=True)):
        log_energy = math.log(energy)
        log_mobility = math.log(mobility)
        lines.append(
            "SWARM_ICP_MOBILITY_TABLE\t"
            f"{index}\t{log_energy:.17g}\t{log_mobility + table_delta:.17g}"
        )
        lines.append(
            "SWARM_ICP_MOBILITY_VALUE\t"
            f"{index}\t{log_energy:.17g}\t{log_mobility + value_delta:.17g}"
        )
    path = plan.report_path.parent / "closure_stdout.txt"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_closure_readback_java_reopens_saved_model_and_reads_every_binding(
    tmp_path: Path,
) -> None:
    plan = _plan(tmp_path)
    source = plan.java_path.read_text(encoding="utf-8")

    assert "ModelUtil.loadCopy" in source
    assert str(plan.mapping.run.output_mph).replace("\\", "/") in source
    for name in (
        "ReducedProps",
        "TensorElectronProps",
        "IncludeThermalDiffusion",
        "MeanElectronEnergyModel",
    ):
        assert f'.getBoolean("{name}")' in source or f'.getString("{name}")' in source
    assert '.feature("pes1").getStringArray("muN")' in source
    assert '.variable("swIcpClosureVars").get("sw_icp_logeps")' in source
    assert 'mobility.getStringMatrix("table")' in source
    assert 'model.param().evaluate("sw_icp_readback_value")' in source
    assert source.count('getString("SpecifyReactionUsing")') == 5
    assert source.count('getString("eedf")') == 5
    assert source.count('getBoolean("UseTownsendTwoTermsBoltzmann")') == 5
    assert 'numerical().create(identityTag, "MaxSurface")' in source
    assert "exp(En-Ne)*1[V]-plas.ebar" in source
    assert "SWARM_ICP_CLOSURE_CONTRACT" in source
    assert "SWARM_ICP_MOBILITY_TABLE" in source
    assert "SWARM_ICP_MOBILITY_VALUE" in source

    apply = generate_apply_java(plan.mapping)
    assert (
        '.prop("ElectronProperties").set("MeanElectronEnergyModel", '
        '"LocalEnergyApproximationE");'
    ) in apply
    assert MeanEnergyArgument(1.0, 10.0).expression("En-Ne") in apply


def test_closure_readback_accepts_complete_saved_contract_and_all_anchors(
    tmp_path: Path,
) -> None:
    plan = _plan(tmp_path)

    result = readback.audit_gec_icp_closure_readback(
        plan,
        _write_readback_log(plan),
    )

    assert result["schema"] == "swarm.gec_icp_closure_readback.v1"
    assert result["status"] == "passed"
    assert result["passed"] is True
    assert all(result["gates"].values())
    assert len(result["mobility_anchor_audit"]["rows"]) == 2
    assert json.loads(plan.report_path.read_text(encoding="utf-8")) == result


@pytest.mark.parametrize(
    ("log_options", "reason"),
    [
        (
            {"overrides": {"ElectronProperties.ReducedProps": "false"}},
            "contract_value_mismatch",
        ),
        (
            {"omit": {"reaction.elastic.eir1.eedf"}},
            "missing_contract_keys",
        ),
        ({"table_delta": 1.0e-4}, "mobility_table_anchor_mismatch"),
        ({"value_delta": 1.0e-4}, "mobility_evaluated_anchor_mismatch"),
        (
            {"overrides": {"saved_solution.mean_energy_relative_error_max": "1e-4"}},
            "saved_solution_mean_energy_identity_mismatch",
        ),
    ],
)
def test_closure_readback_fails_closed_on_contract_or_numeric_mismatch(
    tmp_path: Path,
    log_options: dict[str, Any],
    reason: str,
) -> None:
    plan = _plan(tmp_path)

    result = readback.audit_gec_icp_closure_readback(
        plan,
        _write_readback_log(plan, **log_options),
    )

    assert result["status"] == "failed"
    assert result["passed"] is False
    assert any(item.startswith(reason) for item in result["failure_reasons"])
    assert plan.report_path.is_file()


def test_closure_readback_rejects_duplicate_sentinel_and_changed_table(
    tmp_path: Path,
) -> None:
    plan = _plan(tmp_path)
    duplicate = _write_readback_log(plan, duplicate_key="transport_sha256")

    duplicate_result = readback.audit_gec_icp_closure_readback(plan, duplicate)
    assert duplicate_result["passed"] is False
    assert duplicate_result["failure_reasons"][0].startswith("invalid_readback_log")

    log = _write_readback_log(plan)
    plan.transport_path.write_text(
        "mean_energy_eV,reduced_mobility_m2_V_s_m3\n1.0,1.0e24\n10.0,2.1e24\n",
        encoding="utf-8",
    )
    changed_result = readback.audit_gec_icp_closure_readback(plan, log)
    assert changed_result["passed"] is False
    assert "transport_sha256_changed" in changed_result["failure_reasons"]


def _execution_plan(tmp_path: Path) -> Any:
    mapping = _mapping(tmp_path)
    closure_plan = readback.prepare_gec_icp_closure_readback(
        mapping,
        output_directory=mapping.output_directory,
        write_java=False,
    )
    native = SimpleNamespace(
        java_path=mapping.output_directory / "Native.java",
        values_path=mapping.output_directory / "native_values.csv",
        contract_path=mapping.output_directory / "native_contract.tsv",
    )
    return SimpleNamespace(
        mapping=mapping,
        mph_contract=None,
        bundle_evidence={},
        output_directory=mapping.output_directory,
        plan_json=mapping.output_directory / "plan.json",
        closure_java=mapping.output_directory / "Closure.java",
        closure_support_java=(
            mapping.output_directory / "Apply.java",
            mapping.output_directory / "Readback.java",
            mapping.output_directory / "Native.java",
        ),
        closure_readback=closure_plan,
        native_eedf_audit=native,
        solve_java=mapping.output_directory / "Solve.java",
        expected_result_files=(
            closure_plan.report_path,
            mapping.output_directory / "comsol_eedf_audit.json",
            mapping.output_directory / "convergence.json",
        ),
    )


def _mock_execution(
    monkeypatch: pytest.MonkeyPatch,
    plan: Any,
    *,
    closure_passed: bool,
) -> list[str]:
    operations: list[str] = []

    def fake_execute(*args: Any, operation: str, **kwargs: Any) -> Any:
        del args, kwargs
        operations.append(operation)
        stdout = plan.output_directory / f"{operation}.txt"
        stdout.parent.mkdir(parents=True, exist_ok=True)
        stdout.write_text("synthetic\n", encoding="utf-8")
        return SimpleNamespace(
            operation=operation,
            total_time_s=1.0,
            log_dir=plan.output_directory,
            result_json=plan.output_directory / f"{operation}_result.json",
            provenance_json=plan.output_directory / f"{operation}_provenance.json",
            comsol_version="6.4",
            comsol_build="synthetic",
            stdout_paths=(stdout,),
        )

    def fake_readback(*args: Any, **kwargs: Any) -> dict[str, Any]:
        del args, kwargs
        result = {
            "schema": "swarm.gec_icp_closure_readback.v1",
            "status": "passed" if closure_passed else "failed",
            "passed": closure_passed,
            "failure_reasons": [] if closure_passed else ["synthetic_mismatch"],
        }
        write_json(plan.closure_readback.report_path, result)
        return result

    monkeypatch.setattr(execution, "prepare_gec_icp_run", lambda *a, **k: plan)
    monkeypatch.setattr(execution, "execute_generated_comsol_java", fake_execute)
    monkeypatch.setattr(execution, "audit_gec_icp_closure_readback", fake_readback)
    monkeypatch.setattr(
        execution,
        "extract_comsol_eedf_audit_log",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        execution,
        "analyze_comsol_eedf_audit",
        lambda *args, **kwargs: {"passed": True},
    )
    monkeypatch.setattr(
        execution,
        "assess_gec_icp_convergence",
        lambda *args, **kwargs: {"passed": True, "gates": {}},
    )
    monkeypatch.setattr(execution, "_repo_root", lambda path: plan.mapping.root)
    return operations


def test_execution_orders_readback_before_native_eedf_and_requires_all_gates(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    plan = _execution_plan(tmp_path)
    operations = _mock_execution(monkeypatch, plan, closure_passed=True)

    summary = execution.execute_gec_icp_run(plan.mapping.path)

    assert operations == ["gec_icp_closure", "gec_icp_solve"]
    status = json.loads(summary.status_json.read_text(encoding="utf-8"))
    assert status["closure_readback"]["passed"] is True
    assert status["native_function_eedf"]["passed"] is True
    assert status["convergence"]["passed"] is True
    assert status["physical_target_accepted"] is True
    assert status["quality_accepted_for_declared_closure"] is True
    assert status["quality_acceptance_scope"] == "declared_restricted_lmea_closure"
    assert status["expected_artifacts"]["closure_readback"] == str(
        plan.closure_readback.report_path
    )


def test_execution_rejects_before_native_eedf_when_readback_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    plan = _execution_plan(tmp_path)
    operations = _mock_execution(monkeypatch, plan, closure_passed=False)

    with pytest.raises(GecIcpQualityError, match="post-apply readback"):
        execution.execute_gec_icp_run(plan.mapping.path)

    assert operations == ["gec_icp_closure"]
    status = json.loads(
        (plan.output_directory / "run_status.json").read_text(encoding="utf-8")
    )
    assert status["status"] == "failed"
    assert status["solve_status"] == "not_started"
    assert status["closure_readback"]["passed"] is False
    assert status["physical_target_accepted"] is False
    assert status["quality_accepted_for_declared_closure"] is False
