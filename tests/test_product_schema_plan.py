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
    assert cfg.solvers.monte_carlo.angular_scattering == "external"
    assert cfg.feature_policy.allow_unsupported_fallback is True


def test_monte_carlo_angular_scattering_schema(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {}
    data["solvers"]["monte_carlo"]["angular_scattering"] = "same_as_physics"
    cfg = load_config(write_config(tmp_path, data))
    assert cfg.solvers.monte_carlo.angular_scattering == "same_as_physics"
    assert cfg.internal.monte_carlo.angular_scattering == "same_as_physics"

    data["solvers"]["monte_carlo"]["angular_scattering"] = "isotropic"
    with pytest.raises(ValueError, match="solvers.monte_carlo.angular_scattering"):
        load_config(write_config(tmp_path, data))


def test_monte_carlo_internal_backend_schema(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "backend": "internal",
        "angular_scattering": "same_as_physics",
        "particles": 8,
        "max_collisions": 4,
        "seed": 123,
    }
    cfg = load_config(write_config(tmp_path, data))
    assert cfg.solvers.monte_carlo.backend == "internal"
    assert cfg.internal.monte_carlo.backend == "internal"
    assert cfg.solvers.monte_carlo.seed == 123

    data["solvers"]["monte_carlo"]["command"] = "echo nope"
    with pytest.raises(ValueError, match="backend=internal"):
        load_config(write_config(tmp_path, data))

    data["solvers"]["monte_carlo"] = {"backend": "internal"}
    with pytest.raises(ValueError, match="same_as_physics"):
        load_config(write_config(tmp_path, data))


def test_monte_carlo_same_as_physics_sampler_policy(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {"angular_scattering": "same_as_physics"}
    cfg = load_config(write_config(tmp_path, data))
    [item] = build_solve_plan(cfg)
    assert item.runnable
    assert (
        item.effective_physics["angular_scattering"]
        == "same_as_physics:isotropic:sampler_supported"
    )

    data["physics"]["angular_scattering"] = {
        "model": "momentum_power",
        "higher_moment_closure": "power",
    }
    with pytest.raises(ValueError, match="no product MC sampler"):
        build_solve_plan(load_config(write_config(tmp_path, data)))

    data["feature_policy"]["unsupported"] = "skip_solver"
    cfg = load_config(write_config(tmp_path, data))
    [item] = build_solve_plan(cfg)
    assert item.skipped
    assert item.effective_physics["angular_scattering"].endswith("unsupported_sampler")

    data["feature_policy"]["unsupported"] = "approximate"
    with pytest.raises(ValueError, match="cannot be silently approximated"):
        build_solve_plan(load_config(write_config(tmp_path, data)))

    data["feature_policy"]["allow_unsupported_fallback"] = True
    cfg = load_config(write_config(tmp_path, data))
    [item] = build_solve_plan(cfg)
    assert item.runnable
    assert item.degraded
    assert item.effective_physics["angular_scattering"].endswith(
        "external_metadata_validation_fallback"
    )

    data["solvers"]["monte_carlo"]["backend"] = "internal"
    with pytest.raises(ValueError, match="internal monte_carlo"):
        build_solve_plan(load_config(write_config(tmp_path, data)))


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

    data["physics"]["electron_electron"] = {
        "enabled": True,
        "model": "fp_energy",
        "strength_model": "simple_relaxation",
        "relaxation_fraction": 0.05,
    }
    cfg = load_config(write_config(tmp_path, data))
    assert cfg.physics.electron_electron.model == "fp_energy"
    assert cfg.physics.electron_electron.strength_model == "simple_relaxation"

    data["physics"]["electron_electron"]["strength_model"] = "density_based"
    cfg = load_config(write_config(tmp_path, data))
    with pytest.raises(NotImplementedError, match="density_based"):
        build_solve_plan(cfg)

    data["physics"]["electron_electron"] = {
        "enabled": True,
        "model": "relaxation_postprocess",
        "strength_model": "density_based",
    }
    cfg = load_config(write_config(tmp_path, data))
    with pytest.raises(NotImplementedError, match="density_based"):
        build_solve_plan(cfg)

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
    data["solvers"]["multiterm_boltzmann"] = {}
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
    data["solvers"]["multi_term"]["lmax"] = 1
    cfg = load_config(write_config(tmp_path, data))
    assert cfg.solvers.multi_term.method == "pn_closure_direct"
    with pytest.raises(NotImplementedError, match="not an independent PN block solve"):
        run(cfg, write=False)

    data = base_product_config(tmp_path, ["multi_term"])
    data["solvers"]["multi_term"]["method"] = "pn_dcs"
    cfg = load_config(write_config(tmp_path, data))
    assert cfg.solvers.multi_term.method == "pn_dcs"
    with pytest.raises(NotImplementedError, match="model=moment_table"):
        run(cfg, write=False)


def test_direct_pn_roadmap_records_implementation_gate() -> None:
    roadmap = ROOT / "docs" / "dev" / "direct_pn_closure_operator.md"
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
        "planned independent coupled block PN solver",
        "Required Tests",
        "Physically Uncertain Parts",
        "two-term native Scharfetter-Gummel",
        "direct_pn_operator=true is not emitted",
        "current shared-SG reduction is not an independent PN block solve",
        "coupled f0/f1 sparse block solve",
        "no two-term scalar operator reuse",
        "no post-hoc f1",
        "PN-equation-derived energy-space field coupling",
        "`ell >= 1` collision damping",
        "`l>0` source/sink treatment",
        "`l>0` boundary conditions",
        "lmax=1 two-term Scharfetter-Gummel regression harness",
    ]
    for phrase in required:
        assert phrase in text


def test_direct_pn_gate_has_no_tracked_placeholder_modules() -> None:
    multiterm_dir = ROOT / "electron_swarm" / "solvers" / "multi_term"
    forbidden = {
        "direct.py",
        "direct_pn.py",
        "direct_pn_operator.py",
        "operator.py",
        "operator_core.py",
        "pn_direct.py",
    }
    present = {path.name for path in multiterm_dir.glob("*.py")}
    assert not (present & forbidden)


def test_capability_matrix_is_canonical_and_minimal() -> None:
    expected = {
        "solver",
        "electron_neutral",
        "angular_scattering",
        "ionization_source",
        "electron_electron",
        "magnetic_field",
        "tail_refinement",
        "bulk_transport",
    }
    for solver in CANONICAL_SOLVER_IDS:
        caps = get_solver_capabilities(solver)
        assert set(asdict(caps)) == expected
        assert all(isinstance(value, SupportLevel) for value in asdict(caps).values() if value != solver)
    assert get_solver_capabilities("two_term").electron_electron == SupportLevel.APPROXIMATE
    assert get_solver_capabilities("multi_term").electron_electron == SupportLevel.APPROXIMATE
    assert all(
        get_solver_capabilities(solver).bulk_transport == SupportLevel.UNSUPPORTED
        for solver in CANONICAL_SOLVER_IDS
    )


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


def test_magnetic_policy_skips_unsupported_pn_but_runs_internal_mc(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["two_term", "multi_term", "monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "backend": "internal",
        "angular_scattering": "same_as_physics",
        "particles": 8,
        "max_collisions": 4,
        "seed": 5,
    }
    data["physics"]["field"]["magnetic_field"] = {
        "enabled": True,
        "B_T": 0.01,
        "angle_EB_deg": 90.0,
    }
    data["feature_policy"]["unsupported"] = "skip_solver"
    cfg = load_config(write_config(tmp_path, data))
    plan = build_solve_plan(cfg)
    by_solver = {item.solver: item for item in plan}
    assert by_solver["two_term"].skipped
    assert by_solver["multi_term"].skipped
    assert by_solver["monte_carlo"].runnable
    assert by_solver["monte_carlo"].degraded
    assert by_solver["monte_carlo"].effective_physics["magnetic_field"] == "approximate"


def test_magnetic_field_schema_validation(tmp_path: Path) -> None:
    data = base_product_config(tmp_path)
    data["physics"]["field"]["magnetic_field"]["B_T"] = -0.1
    with pytest.raises(ValueError, match="B_T"):
        load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["physics"]["field"]["magnetic_field"]["B_T"] = ".inf"
    with pytest.raises(ValueError, match="B_T"):
        load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["physics"]["field"]["magnetic_field"]["angle_EB_deg"] = ".nan"
    with pytest.raises(ValueError, match="angle_EB_deg"):
        load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["physics"]["field"]["magnetic_field"]["angle_EB_deg"] = 181.0
    with pytest.raises(ValueError, match="angle_EB_deg"):
        load_config(write_config(tmp_path, data))


def test_finite_k_is_schema_validated_and_policy_handled(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["two_term", "multi_term"])
    data["physics"]["finite_k"] = {"enabled": True, "k_m_inv": 100.0}
    with pytest.raises(ValueError, match="finite k"):
        run(load_config(write_config(tmp_path, data)), write=False)

    data["feature_policy"]["unsupported"] = "skip_solver"
    cfg = load_config(write_config(tmp_path, data))
    plan = build_solve_plan(cfg)
    assert all(item.skipped for item in plan)
    assert {item.effective_physics["finite_k"] for item in plan} == {"unsupported"}
    result = run(cfg, write=False)
    assert result.cases == []
    assert all(row["effective_finite_k"] == "unsupported" for row in result.metadata["solver_plan"])

    data["feature_policy"]["unsupported"] = "approximate"
    with pytest.raises(ValueError, match="cannot be silently approximated"):
        run(load_config(write_config(tmp_path, data)), write=False)


@pytest.mark.parametrize("k_m_inv", [None, 0.0, -1.0, ".inf", ".nan"])
def test_finite_k_requires_positive_finite_wavenumber(
    tmp_path: Path,
    k_m_inv: object,
) -> None:
    data = base_product_config(tmp_path)
    data["physics"]["finite_k"] = {"enabled": True, "k_m_inv": k_m_inv}
    with pytest.raises(ValueError, match="k_m_inv"):
        load_config(write_config(tmp_path, data))


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
            ionization_source=SupportLevel.EXACT,
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
    data["solvers"]["monte_carlo"] = {}
    data["solvers"]["monte_carlo"]["python_api"] = "product_helpers:fake_mc_missing"
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

    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "python_api": "product_helpers:fake_mc_missing",
    }
    data["physics"]["electron_electron"] = {
        "enabled": True,
        "model": "fp_energy",
        "strength_model": "simple_relaxation",
    }
    with pytest.raises(ValueError, match="electron electron"):
        run(load_config(write_config(tmp_path, data)), write=False)

    data["feature_policy"]["unsupported"] = "skip_solver"
    cfg = load_config(write_config(tmp_path, data))
    [item] = build_solve_plan(cfg)
    assert item.skipped
    assert item.effective_physics["electron_electron"] == "unsupported"


def test_ionization_source_model_policy_and_metadata(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["two_term"])
    data["physics"]["ionization"] = {
        "energy_sharing": "primary_secondary",
        "secondary_electron_energy_eV": 0.2,
    }
    cfg = load_config(write_config(tmp_path, data))
    [item] = build_solve_plan(cfg)
    assert item.runnable
    assert item.effective_physics["ionization_source"] == "primary_secondary"
    result = run(cfg, write=False)
    [case] = result.cases
    assert case.metadata["ionization_source_model"] == "primary_secondary"
    assert case.metadata["ionization_source_treatment"] == "primary_secondary"
    assert case.metadata["ionization_secondary_electron_energy_eV"] == 0.2

    data = base_product_config(tmp_path, ["multi_term"])
    data["physics"]["ionization"] = {
        "energy_sharing": "primary_secondary",
        "secondary_electron_energy_eV": 0.2,
    }
    with pytest.raises(ValueError, match="ionization source"):
        run(load_config(write_config(tmp_path, data)), write=False)

    data["feature_policy"]["unsupported"] = "skip_solver"
    cfg = load_config(write_config(tmp_path, data))
    [item] = build_solve_plan(cfg)
    assert item.skipped
    assert item.effective_physics["ionization_source"] == "unsupported"


def test_ionization_secondary_energy_validation(tmp_path: Path) -> None:
    data = base_product_config(tmp_path)
    data["physics"]["ionization"] = {
        "energy_sharing": "primary_secondary",
        "secondary_electron_energy_eV": -0.1,
    }
    with pytest.raises(ValueError, match="secondary_electron_energy_eV"):
        load_config(write_config(tmp_path, data))

    data["physics"]["ionization"] = {
        "energy_sharing": "equal",
        "secondary_electron_energy_eV": 0.1,
    }
    with pytest.raises(ValueError, match="only used with primary_secondary"):
        load_config(write_config(tmp_path, data))
