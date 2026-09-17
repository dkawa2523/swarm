"""Product adapter for the homogeneous DC energy-angle propagator."""

from __future__ import annotations

from time import perf_counter

import numpy as np

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.constants import E_CHARGE_C, ELECTRON_MASS_KG, TOWNSEND
from electron_swarm.core.cross_sections import ActiveMixtureInputs
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.core.solver_configs import PropagatorInternalConfig
from electron_swarm.physics.kinetics import gas_number_density
from electron_swarm.solvers.base import IndependentCaseSolver
from electron_swarm.solvers.propagator.collisions import (
    build_collision_operator,
)
from electron_swarm.solvers.propagator.case_result import build_case_result
from electron_swarm.solvers.propagator.grid import (
    build_propagator_grid,
)
from electron_swarm.solvers.propagator.memory import preflight_memory_bytes
from electron_swarm.solvers.propagator.operator import PropagatorOperator
from electron_swarm.solvers.propagator.qualification import qualify_p1_inputs
from electron_swarm.solvers.propagator.steady import (
    SteadyControls,
    solve_steady_state,
)


class PropagatorSolver(IndependentCaseSolver):
    name = "propagator"

    def __init__(
        self,
        config: SwarmConfig,
        cross_sections: ActiveMixtureInputs,
        solver_config: PropagatorInternalConfig,
    ) -> None:
        super().__init__(config, cross_sections)
        self.solver_config = solver_config
        self._qualification = qualify_p1_inputs(config, self.cross_sections)
        self._operators: dict[float, PropagatorOperator] = {}
        estimate = preflight_memory_bytes(
            energy_cells=solver_config.energy_cells,
            polar_cells=solver_config.polar_cells,
            angular_model=config.physics.angular_scattering.model,
        )
        limit = int(solver_config.max_memory_mb) * 1024 * 1024
        if estimate > limit:
            raise MemoryError(
                "propagator preflight estimate "
                f"{estimate} bytes exceeds max_memory_mb="
                f"{solver_config.max_memory_mb}"
            )

    def _operator(
        self,
        max_eV: float,
    ) -> tuple[PropagatorOperator, dict[str, float]]:
        key = float(max_eV)
        cached = self._operators.get(key)
        if cached is not None:
            return cached, {
                "grid_build": 0.0,
                "collision_kernel_build": 0.0,
                "operator_cache_hit": 1.0,
            }
        started = perf_counter()
        grid = build_propagator_grid(
            self.solver_config,
            self.cross_sections,
            max_eV=max_eV,
        )
        grid_time = perf_counter() - started
        estimate = preflight_memory_bytes(
            energy_cells=grid.energy_cells,
            polar_cells=grid.polar_cells,
            angular_model=self.config.physics.angular_scattering.model,
        )
        limit = int(self.solver_config.max_memory_mb) * 1024 * 1024
        if estimate > limit:
            raise MemoryError(
                f"propagator grid estimate {estimate} bytes exceeds "
                f"max_memory_mb={self.solver_config.max_memory_mb}"
            )
        collision_started = perf_counter()
        collisions = build_collision_operator(
            self.config,
            self.cross_sections,
            grid,
        )
        collision_time = perf_counter() - collision_started
        operator = PropagatorOperator.build(grid, collisions)
        if operator.memory_bytes > limit:
            raise MemoryError(
                f"propagator operator uses {operator.memory_bytes} bytes, "
                f"exceeding max_memory_mb={self.solver_config.max_memory_mb}"
            )
        self._operators[key] = operator
        return operator, {
            "grid_build": grid_time,
            "collision_kernel_build": collision_time,
            "operator_cache_hit": 0.0,
        }

    def _solve_case_with_state(
        self,
        e_over_n_Td: float,
        case_id: str,
        warm_operator: PropagatorOperator | None,
        warm_population: np.ndarray | None,
    ) -> tuple[SwarmCaseResult, PropagatorOperator, np.ndarray]:
        total_started = perf_counter()
        density = gas_number_density(self.config)
        electric_field = float(e_over_n_Td) * TOWNSEND * density
        acceleration = E_CHARGE_C * electric_field / ELECTRON_MASS_KG
        max_eV = (
            self.solver_config.max_eV_limit
            if self.solver_config.adaptive_grid
            else min(
                self.solver_config.base_max_eV,
                self.solver_config.max_eV_limit,
            )
        )
        timing = {
            "grid_build": 0.0,
            "collision_kernel_build": 0.0,
            "steady_iteration": 0.0,
            "observables": 0.0,
        }

        operator, build_timing = self._operator(max_eV)
        timing["grid_build"] += build_timing["grid_build"]
        timing["collision_kernel_build"] += build_timing[
            "collision_kernel_build"
        ]
        initial = (
            warm_population
            if warm_operator is not None
            and np.array_equal(
                warm_operator.grid.energy_edges_eV,
                operator.grid.energy_edges_eV,
            )
            else None
        )
        steady_started = perf_counter()
        solution = solve_steady_state(
            operator,
            acceleration,
            SteadyControls(
                max_iterations=int(self.solver_config.max_iterations),
                tolerance=self.solver_config.convergence_tolerance,
                tail_probability_target=(
                    self.solver_config.tail_probability_target
                ),
                tail_energy_start_eV=0.8 * float(max_eV),
            ),
            initial_population=initial,
        )
        timing["steady_iteration"] += perf_counter() - steady_started

        case = build_case_result(
            config=self.config,
            cross_sections=self.cross_sections,
            solver_config=self.solver_config,
            qualification=self._qualification,
            operator=operator,
            solution=solution,
            gas_number_density_m3=density,
            e_over_n_Td=e_over_n_Td,
            case_id=case_id,
            solver_name=self.name,
            timing=timing,
            total_started=total_started,
        )
        return case, operator, solution.population

    def solve_case(self, e_over_n_Td: float, case_id: str) -> SwarmCaseResult:
        case, _operator, _population = self._solve_case_with_state(
            e_over_n_Td,
            case_id,
            None,
            None,
        )
        return case

    def solve_all(self) -> list[SwarmCaseResult]:
        cases: list[SwarmCaseResult] = []
        operator: PropagatorOperator | None = None
        population: np.ndarray | None = None
        for index, e_over_n_Td in enumerate(self.config.run.e_over_n_Td):
            case, operator, population = self._solve_case_with_state(
                e_over_n_Td,
                f"{self.config.run.case_prefix}_{index:04d}",
                operator,
                population,
            )
            cases.append(case)
        return cases
