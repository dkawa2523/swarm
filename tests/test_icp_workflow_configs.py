from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from electron_swarm import load_config
from electron_swarm.orchestration.plan import build_solve_plan
from swarm_workflow.campaign.config import load_workflow
from swarm_workflow.selection import physical_context
from swarm_workflow.quality.monte_carlo.policy import sampling_budget_provenance


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES = ROOT / "examples"
COMMON_ANCHORS = (
    1.0,
    2.0,
    3.0,
    5.0,
    10.0,
    20.0,
    30.0,
    50.0,
    100.0,
    300.0,
    1000.0,
    1500.0,
    2500.0,
)
SOLVER_SUPPORT_EXTENSIONS = {
    "two_term": (3000.0,),
    "monte_carlo": (2726.6505889,),
    "propagator": (2744.2936035,),
}
BASE_CONFIGS = {
    "two_term": "argon_gec_icp_two_term.yaml",
    "monte_carlo": "argon_gec_icp_monte_carlo_weighted_branching.yaml",
    "propagator": "argon_gec_icp_propagator.yaml",
}
WORKFLOWS = {
    "two_term": "workflow_argon_gec_icp_two_term.yaml",
    "monte_carlo": "workflow_argon_gec_icp_monte_carlo.yaml",
    "propagator": "workflow_argon_gec_icp_propagator.yaml",
}


@pytest.mark.parametrize("solver_id", tuple(BASE_CONFIGS))
def test_icp_base_configs_share_one_physical_context(solver_id: str) -> None:
    config = load_config(EXAMPLES / BASE_CONFIGS[solver_id])

    assert config.schema_version == 2
    assert [item.id for item in config.run.solvers] == [solver_id]
    assert config.conditions.gas_temperature_K == pytest.approx(300.0)
    assert config.conditions.pressure_Pa == pytest.approx(2.66644)
    assert [item.species for item in config.conditions.gas_mixture] == ["Ar"]
    assert config.conditions.gas_mixture[0].fraction == pytest.approx(1.0)
    assert config.cross_sections.files[0].path.name == "argon_application_library.csv"
    assert config.cross_sections.high_energy_extrapolation == "zero"
    assert config.physics.field.type == "dc"
    assert config.physics.field.magnetic_field.enabled is False
    assert config.physics.angular_scattering.model == "isotropic"
    assert config.physics.electron_electron.enabled is False
    [planned] = build_solve_plan(config)
    assert planned.solver == solver_id
    assert planned.runnable is True
    assert planned.skipped is False


@pytest.mark.parametrize("solver_id", tuple(WORKFLOWS))
def test_icp_workflows_share_comparison_anchors_and_bound_support_extensions(
    solver_id: str,
) -> None:
    workflow = load_workflow(EXAMPLES / WORKFLOWS[solver_id])

    assert workflow.e_over_n_Td == (
        COMMON_ANCHORS + SOLVER_SUPPORT_EXTENSIONS[solver_id]
    )
    assert workflow.mixtures[0].fractions == {"Ar": 1.0}
    assert workflow.database_path.parent.name == "argon_gec_icp"
    config = load_config(workflow.base_config_path)
    assert [item.id for item in config.run.solvers] == [solver_id]


def test_icp_solver_sources_have_compatible_physical_contexts() -> None:
    contexts = {
        solver_id: physical_context(load_config(EXAMPLES / filename))
        for solver_id, filename in BASE_CONFIGS.items()
    }

    assert contexts["two_term"] == contexts["monte_carlo"]
    assert contexts["two_term"] == contexts["propagator"]


def test_icp_monte_carlo_plan_has_explicit_finite_budget() -> None:
    workflow = load_workflow(EXAMPLES / WORKFLOWS["monte_carlo"])
    gec_workflow = load_workflow(EXAMPLES / "workflow_argon_gec_ccp_monte_carlo.yaml")
    config = load_config(workflow.base_config_path)
    assert workflow.mc_e_over_n_Td == (
        COMMON_ANCHORS + SOLVER_SUPPORT_EXTENSIONS["monte_carlo"]
    )
    assert workflow.mc_replicas == 4
    assert len(workflow.mc_sampling_plan) == len(workflow.mc_e_over_n_Td)
    gec_rows = {row.e_over_n_Td: row for row in gec_workflow.mc_sampling_plan}
    for row in workflow.mc_sampling_plan:
        if row.e_over_n_Td not in gec_rows or row.e_over_n_Td == 5.0:
            continue
        expected = gec_rows[row.e_over_n_Td]
        assert (
            replace(row, transport_estimator=expected.transport_estimator) == expected
        )
    rows = {row.e_over_n_Td: row for row in workflow.mc_sampling_plan}
    assert rows[5.0].warmup_collisions == 524_288
    assert rows[5.0].max_collisions == 1_048_576
    assert rows[2726.6505889].particles == 1024
    assert rows[2726.6505889].max_collisions == 32_768
    assert [
        row.e_over_n_Td
        for row in workflow.mc_sampling_plan
        if row.transport_estimator == "paired_field_parity"
    ] == [1.0]
    assert all(
        row.transport_estimator == "single_field"
        for row in gec_workflow.mc_sampling_plan
    )
    assert workflow.database_path.name == "monte_carlo_function_eedf.sqlite"
    assert workflow.mc_reuse_database is None
    assert workflow.mc_previous_decision is None
    assert config.solvers.monte_carlo.population_model == "weighted_branching"
    assert config.solvers.monte_carlo.numeric_kernel == "auto"

    budget = sampling_budget_provenance(workflow.mc_sampling_plan)
    assert budget["maximum_per_replica_particle_barriers"] == 545_259_520
    assert budget["total_particle_barriers"] == 9_990_832_128
    assert workflow.mc_convergence is not None
    assert workflow.mc_convergence.maximum_particle_barriers_per_replica == 600_000_000
    assert workflow.mc_convergence.maximum_total_particle_barriers == 10_000_000_000


def test_icp_propagator_uses_bounded_high_field_profile() -> None:
    workflow = load_workflow(EXAMPLES / WORKFLOWS["propagator"])
    config = load_config(workflow.base_config_path)
    propagator = config.solvers.propagator

    assert propagator.energy_cells == 300
    assert propagator.polar_cells == 48
    assert propagator.max_iterations == 2000
    assert propagator.convergence_tolerance == pytest.approx(1.0e-8)
    assert propagator.max_memory_mb == 256
    assert workflow.deterministic_execution.workers == 4
    assert workflow.deterministic_execution.global_memory_budget_mb == 2048
