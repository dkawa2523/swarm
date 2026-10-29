"""Snapshot requested inputs, resolve defaults, and bind an execution to them."""

from copy import deepcopy
from dataclasses import asdict, dataclass

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.config_validation import validate_config
from electron_swarm.core.cross_sections import (
    ActiveMixtureInputs,
    load_active_mixture_inputs,
)
from electron_swarm.core.solver_configs import (
    InternalSolverConfigs,
    build_internal_solver_configs,
    resolve_monte_carlo_sampling,
)
from electron_swarm.core.source_identity import (
    input_evidence,
    input_file_hashes,
    solver_implementations,
)
from .plan import SolverPlanItem, build_solve_plan, solver_plan_metadata
from electron_swarm.core.physical_context import (
    PhysicalContext,
    collision_data_identity,
    physical_context,
)


@dataclass(frozen=True, slots=True)
class PreparedRun:
    config: SwarmConfig
    inputs: ActiveMixtureInputs
    internal: InternalSolverConfigs
    plan: list[SolverPlanItem]
    configuration: dict
    implementations: dict
    input_files: dict[str, str]
    resolved_solvers: dict
    input_evidence: dict
    collision_identity: str
    physical_cases: tuple[PhysicalContext, ...]

    def verify_unchanged(self) -> None:
        if asdict(self.config) != self.configuration:
            raise RuntimeError(
                "solver configuration changed during execution; discard the run"
            )
        if input_file_hashes(self.config) != self.input_files:
            raise RuntimeError(
                "solver input files changed during execution; discard the run and retry"
            )
        if collision_data_identity(self.inputs) != self.collision_identity:
            raise RuntimeError(
                "loaded collision data changed during execution; discard the run"
            )
        if solver_implementations(self.implementations) != self.implementations:
            raise RuntimeError(
                "solver implementation or numerical runtime changed during execution; discard the run"
            )

    def description(self) -> dict:
        runnable = sum(item.runnable for item in self.plan)
        return {
            "status": "planned" if runnable else "skipped",
            "solver_plan": solver_plan_metadata(self.plan),
            "resolved_solvers": self.resolved_solvers,
            "expected_case_count": runnable * len(self.config.run.e_over_n_Td),
            "physical_cases": [asdict(case) for case in self.physical_cases],
            "export": self.config.output.export_plan(),
            "figures": self.config.output.figure_plan(
                runnable=[item.solver for item in self.plan if item.runnable],
                field_type=self.config.physics.field.type,
            ),
        }


def prepare_run(config: SwarmConfig, *, collect_audit: bool = False) -> PreparedRun:
    working = deepcopy(config)
    validate_config(working)
    configuration = asdict(working)
    input_files = input_file_hashes(working)
    inputs = load_active_mixture_inputs(working.cross_sections, working.conditions)
    internal = build_internal_solver_configs(working.solvers, working.physics)
    internal.monte_carlo.collect_audit = collect_audit
    plan = build_solve_plan(working, inputs)
    runnable = [item.solver for item in plan if item.runnable]
    working.output.figure_plan(runnable=runnable, field_type=working.physics.field.type)
    if (
        working.output.export is not None
        and working.output.export.solver not in runnable
    ):
        raise ValueError(
            "output.export.solver is not runnable under the requested feature policy"
        )
    if "monte_carlo" in runnable and working.solvers.monte_carlo.seed is None:
        raise ValueError(
            "solvers.monte_carlo.seed is required before planning Monte Carlo execution"
        )
    collision_identity = collision_data_identity(inputs)
    resolved = {solver: asdict(getattr(working.solvers, solver)) for solver in runnable}
    if "propagator" in resolved:
        resolved["propagator"]["energy_cells_semantics"] = (
            "core_cells; thresholds_and_bounded_tail_add_reported_cells"
        )
    if "multi_term" in resolved:
        resolved["multi_term"]["retained_angular_moments"] = (
            working.solvers.multi_term.lmax + 1
        )
    if "monte_carlo" in resolved:
        resolved["monte_carlo"].update(
            asdict(resolve_monte_carlo_sampling(working.solvers.monte_carlo))
        )
        resolved["monte_carlo"].update(
            {
                "sampling_unit": "null_collision_trial_periods",
                "production_legs": 2
                if working.solvers.monte_carlo.transport_estimator
                == "paired_field_parity"
                else 1,
                "tail_legs": 1 if resolved["monte_carlo"]["tail_trial_periods"] else 0,
            }
        )
    prepared = PreparedRun(
        working,
        inputs,
        internal,
        plan,
        configuration,
        solver_implementations(runnable),
        input_files,
        resolved,
        input_evidence(working, input_files),
        collision_identity,
        tuple(
            physical_context(working, collision_identity, value)
            for value in working.run.e_over_n_Td
        ),
    )
    prepared.verify_unchanged()
    return prepared
