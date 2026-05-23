from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

import pytest

from electron_swarm import load_config, run
from electron_swarm.core.capabilities import (
    SolverCapabilities,
    SupportLevel,
    get_solver_capabilities,
)
from electron_swarm.core.config import CANONICAL_SOLVER_IDS
from electron_swarm.orchestration.plan import build_solve_plan
import electron_swarm.orchestration.plan as plan_module

from product_helpers import ROOT, base_product_config, write_config


def test_schema_v2_valid_config_and_canonical_solver_ids(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["two_term", "multi_term", "monte_carlo"])
    data["feature_policy"]["allow_unsupported_fallback"] = True
    cfg = load_config(write_config(tmp_path, data))
    assert cfg.schema_version == 2
    assert [item.id for item in cfg.run.solvers] == [
        "two_term",
        "multi_term",
        "monte_carlo",
    ]
    assert not hasattr(cfg, "boltzmann_two_term")
    assert not hasattr(cfg, "multiterm_boltzmann")
    assert not hasattr(cfg, "monte_carlo")
    assert cfg.internal.two_term.backend
    assert cfg.internal.multi_term.product_method == "pn_closure_surrogate"
    assert cfg.feature_policy.allow_unsupported_fallback is True


def test_electron_electron_schema_is_typed(tmp_path: Path) -> None:
    data = base_product_config(tmp_path)
    data["physics"]["electron_electron"] = {
        "enabled": True,
        "model": "relaxation_postprocess",
        "relaxation_fraction": 0.05,
        "conserve_mean_energy": True,
    }
    cfg = load_config(write_config(tmp_path, data))
    assert cfg.physics.electron_electron.model == "relaxation_postprocess"

    data["physics"]["electron_electron"]["model"] = "relaxation"
    with pytest.raises(ValueError, match="relaxation_postprocess"):
        load_config(write_config(tmp_path, data))

    data["physics"]["electron_electron"] = {
        "enabled": True,
        "model": "relaxation_postprocess",
        "relaxation_fraction": 1.5,
    }
    with pytest.raises(ValueError, match="relaxation_fraction"):
        load_config(write_config(tmp_path, data))


def test_tail_metrics_schema_validation(tmp_path: Path) -> None:
    data = base_product_config(tmp_path)
    data["physics"]["energy_grid_policy"]["tail_metrics"] = False
    data["physics"]["energy_grid_policy"]["tail_threshold_eV"] = 12.0
    data["physics"]["energy_grid_policy"]["tail_rate_warning_fraction"] = 0.2
    cfg = load_config(write_config(tmp_path, data))
    policy = cfg.physics.energy_grid_policy
    assert policy.tail_metrics is False
    assert policy.tail_threshold_eV == 12.0
    assert policy.tail_rate_warning_fraction == 0.2

    data["physics"]["energy_grid_policy"]["tail_threshold_eV"] = -1.0
    with pytest.raises(ValueError, match="tail_threshold_eV"):
        load_config(write_config(tmp_path, data))

    data["physics"]["energy_grid_policy"]["tail_threshold_eV"] = None
    data["physics"]["energy_grid_policy"]["tail_rate_warning_fraction"] = 1.5
    with pytest.raises(ValueError, match="tail_rate_warning_fraction"):
        load_config(write_config(tmp_path, data))

    data["physics"]["energy_grid_policy"]["tail_rate_warning_fraction"] = 0.2
    data["physics"]["energy_grid_policy"]["tail_metrics"] = "false"
    with pytest.raises(ValueError, match="tail_metrics"):
        load_config(write_config(tmp_path, data))

    data["physics"]["energy_grid_policy"]["tail_metrics"] = False
    data["physics"]["energy_grid_policy"]["tail_metric"] = False
    with pytest.raises(ValueError, match="Unsupported physics.energy_grid_policy fields"):
        load_config(write_config(tmp_path, data))


def test_schema_v2_rejects_non_boolean_and_nested_unknown_fields(tmp_path: Path) -> None:
    data = base_product_config(tmp_path)
    data["feature_policy"]["allow_unsupported_fallback"] = "false"
    with pytest.raises(ValueError, match="allow_unsupported_fallback"):
        load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["comparison"]["enabled"] = "false"
    with pytest.raises(ValueError, match="comparison.enabled"):
        load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["physics"]["field"]["magnetic_field"]["enabled"] = "false"
    with pytest.raises(ValueError, match="magnetic_field.enabled"):
        load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["physics"]["field"]["magnetic_field"]["solver"] = "full_orbit"
    with pytest.raises(ValueError, match="Unsupported physics.field.magnetic_field"):
        load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["comparison"]["capabilities"] = True
    with pytest.raises(ValueError, match="Unsupported comparison fields"):
        load_config(write_config(tmp_path, data))


def test_schema_v2_rejects_old_public_fields(tmp_path: Path) -> None:
    data = base_product_config(tmp_path)
    data.pop("schema_version")
    with pytest.raises(ValueError, match="schema v2 is required"):
        load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["run"]["mode"] = "both"
    with pytest.raises(ValueError, match="schema v2 is required"):
        load_config(write_config(tmp_path, data))

    for old_id in ["both", "all", "boltzmann_two_term", "multiterm_boltzmann"]:
        data = base_product_config(tmp_path, [old_id])
        with pytest.raises(ValueError, match="obsolete solver id"):
            load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["solvers"]["boltzmann_two_term"] = {}
    with pytest.raises(ValueError, match="obsolete solver id"):
        load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["output"]["write_plots"] = False
    with pytest.raises(ValueError, match="fixed canonical outputs"):
        load_config(write_config(tmp_path, data))


def test_multi_term_product_method_validation(tmp_path: Path) -> None:
    for old_method in ["moment_closure", "operator", "hybrid"]:
        data = base_product_config(tmp_path, ["multi_term"])
        data["solvers"]["multi_term"]["method"] = old_method
        with pytest.raises(ValueError, match="solvers.multi_term.method"):
            load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path, ["multi_term"])
    data["solvers"]["multi_term"]["allow_experimental_operator"] = True
    with pytest.raises(ValueError, match="Unsupported solvers.multi_term fields"):
        load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path, ["multi_term"])
    data["solvers"]["multi_term"]["method"] = "pn_closure_direct"
    cfg = load_config(write_config(tmp_path, data))
    assert cfg.solvers.multi_term.method == "pn_closure_direct"
    with pytest.raises(NotImplementedError, match="direct PN block operator"):
        run(cfg, write=False)

    data = base_product_config(tmp_path, ["multi_term"])
    data["solvers"]["multi_term"]["method"] = "pn_dcs"
    cfg = load_config(write_config(tmp_path, data))
    assert cfg.solvers.multi_term.method == "pn_dcs"
    with pytest.raises(NotImplementedError, match="DCS angular"):
        run(cfg, write=False)


def test_direct_pn_roadmap_records_implementation_gate() -> None:
    roadmap = ROOT / "docs" / "roadmap" / "direct_pn_closure_operator.md"
    text = roadmap.read_text(encoding="utf-8")
    required = [
        "Candidate Operator Equation",
        "Unknown Vector Layout",
        "Field Coupling Block",
        "Collision Damping Block",
        "Source And Sink Treatment",
        "Boundary Conditions",
        "Normalization Constraint",
        "lmax=1 Regression Condition",
        "Required Tests",
        "Physically Uncertain Parts",
        "two-term native Scharfetter-Gummel",
        "direct_pn_operator=false",
    ]
    for phrase in required:
        assert phrase in text


def test_capability_matrix_is_canonical_and_minimal() -> None:
    expected = {
        "solver",
        "electron_neutral",
        "angular_scattering",
        "electron_electron",
        "magnetic_field",
        "tail_refinement",
        "bulk_transport",
    }
    for solver in CANONICAL_SOLVER_IDS:
        caps = get_solver_capabilities(solver)
        assert set(asdict(caps)) == expected
        assert all(isinstance(value, SupportLevel) for value in asdict(caps).values() if value != solver)


def test_solver_plan_and_unsupported_feature_policy(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["two_term", "multi_term"])
    cfg = load_config(write_config(tmp_path, data))
    plan = build_solve_plan(cfg)
    assert [item.solver for item in plan] == ["two_term", "multi_term"]
    assert all(item.runnable for item in plan)

    data["physics"]["field"]["magnetic_field"] = {
        "enabled": True,
        "B_T": 0.05,
        "angle_EB_deg": 90.0,
    }
    with pytest.raises(ValueError, match="magnetic field"):
        run(load_config(write_config(tmp_path, data)), write=False)

    data["feature_policy"]["unsupported"] = "skip_solver"
    cfg = load_config(write_config(tmp_path, data))
    skipped = build_solve_plan(cfg)
    assert all(item.skipped for item in skipped)
    result = run(cfg, write=False)
    assert result.cases == []
    assert all(row["skipped"] for row in result.metadata["solver_plan"])

    data["feature_policy"]["unsupported"] = "approximate"
    with pytest.raises(ValueError, match="cannot be silently approximated"):
        run(load_config(write_config(tmp_path, data)), write=False)

    data["feature_policy"]["allow_unsupported_fallback"] = True
    cfg = load_config(write_config(tmp_path, data))
    fallback_plan = build_solve_plan(cfg)
    assert all(item.runnable for item in fallback_plan)
    assert all(item.degraded for item in fallback_plan)
    assert {
        item.effective_physics["magnetic_field"] for item in fallback_plan
    } == {"ignored_fallback"}
    result = run(cfg, write=False)
    assert result.cases
    assert {
        case.metadata["magnetic_field_treatment"] for case in result.cases
    } == {"ignored_fallback"}


def test_degraded_policy_applies_to_angular_and_tail_features(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["multi_term"])
    data["feature_policy"]["degraded"] = "fail"
    with pytest.raises(ValueError, match="angular_scattering"):
        build_solve_plan(load_config(write_config(tmp_path, data)))


def test_degraded_policy_applies_to_tail_refinement(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data = base_product_config(tmp_path, ["two_term"])
    data["physics"]["energy_grid_policy"]["adaptive"] = True
    data["feature_policy"]["degraded"] = "fail"
    cfg = load_config(write_config(tmp_path, data))

    def fake_capabilities(solver: str) -> SolverCapabilities:
        return SolverCapabilities(
            solver=solver,
            electron_neutral=SupportLevel.EXACT,
            angular_scattering=SupportLevel.EXACT,
            electron_electron=SupportLevel.EXACT,
            magnetic_field=SupportLevel.EXACT,
            tail_refinement=SupportLevel.APPROXIMATE,
            bulk_transport=SupportLevel.EXACT,
        )

    monkeypatch.setattr(plan_module, "get_solver_capabilities", fake_capabilities)
    with pytest.raises(ValueError, match="tail_refinement"):
        build_solve_plan(cfg)


def test_electron_electron_unsupported_solver_policy(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["physics"]["electron_electron"] = {
        "enabled": True,
        "model": "relaxation_postprocess",
        "relaxation_fraction": 0.05,
    }
    with pytest.raises(ValueError, match="electron electron"):
        run(load_config(write_config(tmp_path, data)), write=False)

    data["feature_policy"]["unsupported"] = "skip_solver"
    cfg = load_config(write_config(tmp_path, data))
    [item] = build_solve_plan(cfg)
    assert item.skipped
    assert item.effective_physics["electron_electron"] == "unsupported"
