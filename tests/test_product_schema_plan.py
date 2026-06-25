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
from electron_swarm.core.result_metadata import PRODUCT_CASE_METADATA_KEYS
from electron_swarm.core.solver_configs import build_internal_solver_configs
from electron_swarm.orchestration.plan import build_solve_plan, solver_plan_metadata
import electron_swarm.orchestration.plan as plan_module

from product_helpers import ROOT, base_product_config, write_config, write_moment_table


def _plan_row(item: object) -> dict[str, object]:
    return solver_plan_metadata([item])[0]


def test_schema_v2_valid_config_and_canonical_solver_ids(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["two_term", "multi_term", "monte_carlo"])
    cfg = load_config(write_config(tmp_path, data))
    assert cfg.schema_version == 2
    assert [item.id for item in cfg.run.solvers] == [
        "two_term",
        "multi_term",
        "monte_carlo",
    ]
    internal = build_internal_solver_configs(cfg.solvers, cfg.physics)
    assert internal.two_term.backend
    assert internal.multi_term.product_method == "pn_closure_direct"
    assert cfg.solvers.monte_carlo.population_model == "fixed_particle_single_daughter"
    assert cfg.feature_policy.degraded == "record"


def test_monte_carlo_rejects_unknown_product_fields(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {"command": "run external mc"}
    with pytest.raises(ValueError, match="Unsupported solvers.monte_carlo fields"):
        load_config(write_config(tmp_path, data))


def test_monte_carlo_internal_schema(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "particles": 8,
        "warmup_collisions": 2,
        "max_collisions": 4,
        "seed": 123,
    }
    cfg = load_config(write_config(tmp_path, data))
    internal = build_internal_solver_configs(cfg.solvers, cfg.physics)
    assert cfg.solvers.monte_carlo.seed == 123
    assert cfg.solvers.monte_carlo.warmup_collisions == 2
    assert internal.monte_carlo.warmup_collisions == 2
    assert cfg.solvers.monte_carlo.population_model == "fixed_particle_single_daughter"

    data["solvers"]["monte_carlo"] = {
        "population_model": "weighted_branching",
        "particles": 8,
    }
    cfg = load_config(write_config(tmp_path, data))
    assert cfg.solvers.monte_carlo.population_model == "weighted_branching"
    internal = build_internal_solver_configs(cfg.solvers, cfg.physics)
    assert internal.monte_carlo.population_model == "weighted_branching"

    data["solvers"]["monte_carlo"] = {
        "population_model": "branching_weighted",
        "particles": 8,
    }
    with pytest.raises(ValueError, match="population_model"):
        load_config(write_config(tmp_path, data))

    data["solvers"]["monte_carlo"] = {
        "population_control": "systematic_resampling",
    }
    with pytest.raises(ValueError, match="Unsupported solvers.monte_carlo fields"):
        load_config(write_config(tmp_path, data))

    data["solvers"]["monte_carlo"] = {
        "target_particles": 10,
        "max_particles": 20,
    }
    with pytest.raises(ValueError, match="Unsupported solvers.monte_carlo fields"):
        load_config(write_config(tmp_path, data))

    data["solvers"]["monte_carlo"] = {
        "warmup_collisions": -1,
    }
    with pytest.raises(ValueError, match="warmup_collisions"):
        load_config(write_config(tmp_path, data))


def test_monte_carlo_same_as_physics_sampler_policy(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {}
    cfg = load_config(write_config(tmp_path, data))
    [item] = build_solve_plan(cfg)
    row = _plan_row(item)
    assert item.runnable
    assert row["effective_angular_scattering"] == (
        "same_as_physics:isotropic:sampler_supported"
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
    row = _plan_row(item)
    assert item.skipped
    assert str(row["effective_angular_scattering"]).endswith("unsupported_sampler")
    assert item.skip_reason is not None

    data["feature_policy"]["unsupported"] = "approximate"
    with pytest.raises(ValueError, match="feature_policy.unsupported"):
        load_config(write_config(tmp_path, data))

    data["feature_policy"]["allow_unsupported_fallback"] = True
    with pytest.raises(ValueError, match="allow_unsupported_fallback"):
        load_config(write_config(tmp_path, data))

    data["feature_policy"].pop("allow_unsupported_fallback")
    data["feature_policy"]["unsupported"] = "fail"
    with pytest.raises(ValueError, match="no product MC sampler"):
        build_solve_plan(load_config(write_config(tmp_path, data)))

    table = write_moment_table(tmp_path)
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {}
    data["physics"]["angular_scattering"] = {
        "model": "moment_table",
        "moment_table": {
            "path": table.as_posix(),
            "format": "normalized_legendre_moments",
            "provenance": "model_derived",
            "extrapolation": "error",
        },
    }
    with pytest.raises(ValueError, match="no product MC sampler"):
        build_solve_plan(load_config(write_config(tmp_path, data, "mc_table.yaml")))


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

    data = base_product_config(tmp_path)
    data["run"]["solvers"] = ["two_term"]
    with pytest.raises(ValueError, match="run.solvers entries must be mappings"):
        load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["run"]["solvers"][0]["label"] = "legacy label"
    with pytest.raises(ValueError, match="Unsupported run.solvers\\[0\\] fields"):
        load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["run"]["E_over_N_Td"] = [50.0]
    with pytest.raises(ValueError, match="Unsupported run fields"):
        load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["conditions"]["unused"] = True
    with pytest.raises(ValueError, match="Unsupported conditions fields"):
        load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["cross_sections"]["unused"] = True
    with pytest.raises(ValueError, match="Unsupported cross_sections fields"):
        load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["feature_policy"]["unsupported"] = "approximate"
    with pytest.raises(ValueError, match="feature_policy.unsupported"):
        load_config(write_config(tmp_path, data))

    for degraded in ("warn", "record_only"):
        data = base_product_config(tmp_path)
        data["feature_policy"]["degraded"] = degraded
        with pytest.raises(ValueError, match="feature_policy.degraded"):
            load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["unexpected"] = True
    with pytest.raises(ValueError, match="Unsupported top-level fields"):
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

    data = base_product_config(tmp_path)
    data["run"]["solvers"] = [{"id": "boltzmann_two_term"}]
    with pytest.raises(ValueError, match="obsolete solver id"):
        load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["unexpected"] = True
    with pytest.raises(ValueError, match="Unsupported top-level fields"):
        load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["output"]["write_plots"] = False
    with pytest.raises(ValueError, match="Unsupported output fields"):
        load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["references"] = {"external": []}
    with pytest.raises(ValueError, match="Unsupported top-level fields"):
        load_config(write_config(tmp_path, data))


def test_multi_term_product_method_validation(tmp_path: Path) -> None:
    for old_method in ["moment_closure", "operator", "hybrid", "pn_closure_surrogate"]:
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
    result = run(cfg, write=False)
    [case] = result.cases
    assert case.metadata["solver_method"] == "pn_closure_direct"
    assert case.metadata["direct_pn_operator"] is True

    data["solvers"]["multi_term"]["lmax"] = 2
    cfg = load_config(write_config(tmp_path, data, name="direct_l2.yaml"))
    result = run(cfg, write=False)
    [case] = result.cases
    assert case.metadata["solver_method"] == "pn_closure_direct"
    assert case.metadata["lmax"] == 2
    assert case.metadata["transport_definition"] == "f0_gradient_reconstruction"
    assert set(case.metadata) <= PRODUCT_CASE_METADATA_KEYS

    data = base_product_config(tmp_path, ["multi_term"])
    data["solvers"]["multi_term"]["method"] = "pn_dcs"
    cfg = load_config(write_config(tmp_path, data))
    assert cfg.solvers.multi_term.method == "pn_dcs"
    with pytest.raises(NotImplementedError, match="model=moment_table"):
        run(cfg, write=False)


def test_direct_pn_roadmap_records_implementation_gate() -> None:
    roadmap = ROOT / "docs" / "dev" / "direct_pn_closure_operator.md"
    text = roadmap.read_text(encoding="utf-8")
    higher_l_gate = (
        ROOT / "docs" / "dev" / "direct_pn_lmax_gt1_gate.md"
    ).read_text(encoding="utf-8")
    required = [
        "Unknown Vector Layout",
        "Current lmax=1 Equation",
        "lmax>1 Coefficient-Space Scope",
        "lmax=1 Regression Condition",
        "Required Tests",
        "SG energy-flux auxiliary",
        "physical Legendre coefficient",
        "no post-hoc `G1` construction",
        "`ell>=2` damping",
        "sink-only",
        "Ar/BOLSIG lmax=1 two-term SG regression harness",
    ]
    for phrase in required:
        assert phrase in text
    for phrase in [
        "Implemented Model",
        "energy-flux auxiliary",
        "`ell>=2` damping",
        "Fail-Fast Conditions",
    ]:
        assert phrase in higher_l_gate


def test_capability_matrix_is_canonical_and_minimal() -> None:
    expected = {
        "solver",
        "angular_scattering",
        "ionization_source",
        "electron_electron",
        "magnetic_field",
        "tail_refinement",
    }
    for solver in CANONICAL_SOLVER_IDS:
        caps = get_solver_capabilities(solver)
        assert set(asdict(caps)) == expected
        assert all(isinstance(value, SupportLevel) for value in asdict(caps).values() if value != solver)
    assert get_solver_capabilities("two_term").electron_electron == SupportLevel.APPROXIMATE
    assert get_solver_capabilities("multi_term").electron_electron == SupportLevel.APPROXIMATE


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
    with pytest.raises(ValueError, match="feature_policy.unsupported"):
        load_config(write_config(tmp_path, data))

    data["feature_policy"]["allow_unsupported_fallback"] = True
    with pytest.raises(ValueError, match="allow_unsupported_fallback"):
        load_config(write_config(tmp_path, data))


def test_magnetic_policy_skips_unsupported_pn_but_runs_internal_mc(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["two_term", "multi_term", "monte_carlo"])
    data["solvers"]["monte_carlo"] = {
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
    assert (
        _plan_row(by_solver["monte_carlo"])["effective_magnetic_field"]
        == "boris_lorentz_push"
    )


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
    assert {
        row["effective_finite_k"] for row in solver_plan_metadata(plan)
    } == {"unsupported"}
    result = run(cfg, write=False)
    assert result.cases == []
    assert all(row["effective_finite_k"] == "unsupported" for row in result.metadata["solver_plan"])

    data["feature_policy"]["unsupported"] = "approximate"
    with pytest.raises(ValueError, match="feature_policy.unsupported"):
        load_config(write_config(tmp_path, data))


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
            angular_scattering=SupportLevel.EXACT,
            ionization_source=SupportLevel.EXACT,
            electron_electron=SupportLevel.EXACT,
            magnetic_field=SupportLevel.EXACT,
            tail_refinement=SupportLevel.APPROXIMATE,
        )

    monkeypatch.setattr(plan_module, "get_solver_capabilities", fake_capabilities)
    with pytest.raises(ValueError, match="tail_refinement"):
        build_solve_plan(cfg)


def test_electron_electron_unsupported_solver_policy(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {}
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
    assert _plan_row(item)["effective_electron_electron"] == "unsupported"

    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {}
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
    assert _plan_row(item)["effective_electron_electron"] == "unsupported"


def test_ionization_source_model_policy_and_metadata(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["two_term"])
    data["physics"]["ionization"] = {
        "energy_sharing": "primary_secondary",
        "secondary_electron_energy_eV": 0.2,
    }
    cfg = load_config(write_config(tmp_path, data))
    [item] = build_solve_plan(cfg)
    assert item.runnable
    assert _plan_row(item)["effective_ionization_source"] == "primary_secondary"
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
    assert _plan_row(item)["effective_ionization_source"] == "unsupported"


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
