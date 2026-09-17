from __future__ import annotations

from copy import deepcopy
from dataclasses import replace
from pathlib import Path
import sqlite3
from types import SimpleNamespace

import pytest
import yaml

from electron_swarm import load_config, run
from swarm_workflow.campaign import mean_energy_support
from swarm_workflow.campaign.config import load_workflow
import swarm_workflow.campaign.provenance as provenance_module
import swarm_workflow.campaign.sweep as sweep_module
from product_helpers import base_product_config, write_config


def _write_support_workflow(
    tmp_path: Path,
    *,
    required_mean_eV: float | None,
    solver: str = "two_term",
) -> Path:
    tmp_path.mkdir(parents=True, exist_ok=True)
    data = base_product_config(tmp_path, [solver])
    data["run"]["e_over_n_Td"] = [10.0]
    base = write_config(tmp_path, data, "base.yaml")
    workflow: dict[str, object] = {
        "base_config": base.name,
        "database": "swarm.sqlite",
        "e_over_n_Td": [10.0, 20.0],
        "mixtures": [{"Ar": 1.0}],
    }
    if solver == "monte_carlo":
        workflow["mc"] = {"base_seed": 7}
    if required_mean_eV is not None:
        workflow["mean_energy_support"] = {
            "required_max_mean_energy_eV": required_mean_eV,
            "relative_guard": 0.0,
            "maximum_steps": 1,
            "maximum_e_over_n_Td": 100.0,
        }
    path = tmp_path / "workflow.yaml"
    path.write_text(yaml.safe_dump(workflow, sort_keys=False), encoding="utf-8")
    return path


def test_upper_extension_targets_guard_when_required_value_is_already_supported() -> (
    None
):
    request = mean_energy_support.propose_upper_mean_energy_extension(
        [1500.0, 2500.0],
        [21.05741938486011, 35.45146815872865],
        35.44309,
        relative_guard=0.1,
    )

    assert request is not None
    assert request.guarded_target_mean_energy_eV == pytest.approx(38.987399)
    assert request.suggested_E_over_N_Td == pytest.approx(2744.293603496929)


def test_upper_extension_applies_step_bound_before_flat_response_overflows() -> None:
    request = mean_energy_support.propose_upper_mean_energy_extension(
        [1.0, 2.0],
        [1.0, 1.0 + 1.0e-12],
        2.0,
        maximum_field_step_factor=3.0,
    )

    assert request is not None
    assert request.suggested_E_over_N_Td == pytest.approx(6.0)


def test_upper_extension_uses_record_envelope_across_nonmonotonic_mean_energy() -> None:
    request = mean_energy_support.propose_upper_mean_energy_extension(
        [1.0, 2.0, 3.0, 4.0],
        [1.0, 2.0, 1.8, 2.2],
        2.4,
        relative_guard=0.0,
    )

    assert request is not None
    assert request.current_max_mean_energy_eV == pytest.approx(2.2)
    assert request.suggested_E_over_N_Td > 4.0


def test_upper_extension_rejects_a_nonmonotonic_upper_endpoint() -> None:
    with pytest.raises(ValueError, match="endpoint is not locally invertible"):
        mean_energy_support.propose_upper_mean_energy_extension(
            [10.0, 100.0, 300.0],
            [2.0, 8.0, 7.0],
            9.0,
        )


def test_deterministic_continuation_reuses_solver_and_recalculates_next_step(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = load_config("examples/argon_gec_icp_propagator.yaml")
    requested_fields: list[float] = []

    def fake_run(step_config: object, *, write: bool) -> object:
        assert write is False
        run_config = step_config.run  # type: ignore[attr-defined]
        assert [item.id for item in run_config.solvers] == ["propagator"]
        field = float(run_config.e_over_n_Td[0])
        requested_fields.append(field)
        # Deliberately under-shoot the first inverse prediction so the API must
        # use the newly calculated point for one bounded continuation step.
        mean = 3.0 if len(requested_fields) == 1 else 4.2
        return SimpleNamespace(
            cases=[
                SimpleNamespace(
                    solver="propagator",
                    e_over_n_Td=field,
                    mean_energy_eV=mean,
                )
            ]
        )

    monkeypatch.setattr(mean_energy_support, "_run_product_config", fake_run)
    result = mean_energy_support.continue_deterministic_mean_energy_support(
        config,
        "propagator",
        [1.0, 2.0],
        [1.0, 2.0],
        3.5,
        relative_guard=0.1,
        maximum_steps=2,
    )

    assert len(requested_fields) == 2
    assert requested_fields[0] == pytest.approx(3.85)
    assert requested_fields[1] > requested_fields[0]
    assert result.reached_target is True
    assert result.final_max_mean_energy_eV == pytest.approx(4.2)
    assert len(result.extension_cases) == 2


def test_monte_carlo_support_cannot_bypass_statistical_campaign() -> None:
    config = load_config("examples/argon_gec_icp_propagator.yaml")

    with pytest.raises(ValueError, match="statistical campaign"):
        mean_energy_support.continue_deterministic_mean_energy_support(
            config,
            "monte_carlo",
            [1.0, 2.0],
            [1.0, 2.0],
            3.0,
        )


def test_workflow_support_config_is_typed_hashed_and_rejects_monte_carlo(
    tmp_path: Path,
) -> None:
    workflow_path = _write_support_workflow(
        tmp_path / "deterministic",
        required_mean_eV=None,
    )
    plain = load_workflow(workflow_path)
    raw = yaml.safe_load(workflow_path.read_text(encoding="utf-8"))
    raw["mean_energy_support"] = {
        "required_max_mean_energy_eV": 100.0,
        "relative_guard": 0.0,
        "maximum_steps": 1,
        "maximum_e_over_n_Td": 100.0,
    }
    workflow_path.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
    supported = load_workflow(workflow_path)

    assert supported.mean_energy_support is not None
    assert supported.mean_energy_support.required_max_mean_energy_eV == 100.0
    assert supported.mean_energy_support.maximum_steps == 1
    assert provenance_module.workflow_config_sha256(plain) != (
        provenance_module.workflow_config_sha256(supported)
    )
    inactive_mc_change = replace(
        plain,
        mc_sampling_plan=tuple(
            replace(row, transport_estimator="paired_field_parity")
            for row in plain.mc_sampling_plan
        ),
    )
    assert plain.mc_enabled is False
    assert provenance_module.workflow_config_sha256(plain) == (
        provenance_module.workflow_config_sha256(inactive_mc_change)
    )

    mc_path = _write_support_workflow(
        tmp_path / "mc",
        required_mean_eV=100.0,
        solver="monte_carlo",
    )
    with pytest.raises(ValueError, match="cannot extend monte_carlo"):
        load_workflow(mc_path)


def test_sweep_adds_one_bounded_support_anchor_and_reuses_it_on_resume(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    workflow_path = _write_support_workflow(
        tmp_path,
        required_mean_eV=100.0,
    )
    original_run = run
    continuation_calls = 0

    def fake_continuation(
        config: object,
        solver_id: str,
        fields: list[float],
        means: list[float],
        required: float,
        **kwargs: object,
    ) -> mean_energy_support.MeanEnergyContinuationResult:
        nonlocal continuation_calls
        continuation_calls += 1
        assert solver_id == "two_term"
        assert fields == [10.0, 20.0]
        assert len(means) == 2
        field = 30.0
        extension_config = deepcopy(config)
        extension_config.run.e_over_n_Td = [field]
        case = original_run(extension_config, write=False).cases[0]
        case.mean_energy_eV = required
        return mean_energy_support.MeanEnergyContinuationResult(
            solver_id="two_term",
            guarded_target_mean_energy_eV=required,
            initial_max_mean_energy_eV=max(means),
            final_max_mean_energy_eV=required,
            reached_target=True,
            extension_cases=(case,),
        )

    monkeypatch.setattr(
        sweep_module,
        "continue_deterministic_mean_energy_support",
        fake_continuation,
    )
    first = sweep_module.run_sweep(workflow_path)
    second = sweep_module.run_sweep(workflow_path)

    assert first.cases_written == 3
    assert second.cases_written == 0
    assert continuation_calls == 1
    with sqlite3.connect(first.database_path) as connection:
        assert connection.execute(
            "SELECT e_over_n_Td FROM cases ORDER BY e_over_n_Td"
        ).fetchall() == [(10.0,), (20.0,), (30.0,)]
        assert connection.execute(
            "SELECT COUNT(*) FROM aggregate_quality"
        ).fetchone() == (3,)


def test_sweep_support_noop_and_unreached_continuation_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    noop_path = _write_support_workflow(
        tmp_path / "noop",
        required_mean_eV=1.0e-6,
    )

    def reject_call(*_args: object, **_kwargs: object) -> object:
        raise AssertionError("satisfied support must not launch continuation")

    monkeypatch.setattr(
        sweep_module,
        "continue_deterministic_mean_energy_support",
        reject_call,
    )
    assert sweep_module.run_sweep(noop_path).cases_written == 2

    failed_path = _write_support_workflow(
        tmp_path / "failed",
        required_mean_eV=100.0,
    )

    def failed_continuation(
        _config: object,
        solver_id: str,
        _fields: list[float],
        means: list[float],
        required: float,
        **_kwargs: object,
    ) -> mean_energy_support.MeanEnergyContinuationResult:
        return mean_energy_support.MeanEnergyContinuationResult(
            solver_id=solver_id,
            guarded_target_mean_energy_eV=required,
            initial_max_mean_energy_eV=max(means),
            final_max_mean_energy_eV=max(means),
            reached_target=False,
            extension_cases=(),
        )

    monkeypatch.setattr(
        sweep_module,
        "continue_deterministic_mean_energy_support",
        failed_continuation,
    )
    with pytest.raises(ValueError, match="after bounded continuation"):
        sweep_module.run_sweep(failed_path)
