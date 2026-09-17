"""Public adapter and orchestration for the native two-term solver."""

from __future__ import annotations

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.cross_sections import ActiveMixtureInputs
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.core.solver_configs import TwoTermInternalConfig
from electron_swarm.solvers.base import IndependentCaseSolver
from electron_swarm.solvers.two_term.models import (
    NativeDistributionResult,
    TimePeriodicDistributionResult,
)
from electron_swarm.solvers.two_term.observables import build_case_result
from electron_swarm.solvers.two_term.steady import (
    solve_native_distribution as solve_stationary_distribution,
)
from electron_swarm.solvers.two_term.time_periodic import (
    solve_time_periodic_distribution as solve_periodic_distribution,
)


class TwoTermSolver(IndependentCaseSolver):
    """Unified electron Boltzmann two-term solver."""

    name = "two_term"

    def __init__(
        self,
        config: SwarmConfig,
        cross_sections: ActiveMixtureInputs,
        solver_config: TwoTermInternalConfig,
    ) -> None:
        super().__init__(config, cross_sections)
        self.solver_config = solver_config

    def solve_case(self, e_over_n_Td: float, case_id: str) -> SwarmCaseResult:
        field_type = self.config.physics.field.type
        if field_type == "time_dependent":
            return self._solve_case_time_periodic(e_over_n_Td, case_id)
        if field_type != "dc":
            raise NotImplementedError(
                "two_term only implements dc and time_dependent fields; "
                "stationary high-frequency rf is not implemented"
            )
        return self._solve_case_native(e_over_n_Td, case_id)

    def _solve_case_native(
        self,
        e_over_n_Td: float,
        case_id: str,
    ) -> SwarmCaseResult:
        native = self.solve_native_distribution(e_over_n_Td)
        if not native.diagnostics.converged:
            raise RuntimeError(
                "two_term stationary distribution did not converge; refusing "
                "to publish EEDF, rates, or transport from the last iterate"
            )
        return build_case_result(
            self.config,
            self.cross_sections,
            self.solver_config,
            e_over_n_Td,
            case_id,
            native.operator_block,
            native.eedf_eV_inv,
            metadata=native.metadata,
            solver_name=self.name,
            temporal_growth_frequency_s_inv=(
                native.diagnostics.growth_frequency_s
                if self.solver_config.nonconservative_model == "growth"
                else None
            ),
        )

    def _solve_case_time_periodic(
        self,
        e_over_n_rms_Td: float,
        case_id: str,
    ) -> SwarmCaseResult:
        periodic = self.solve_time_periodic_distribution(
            e_over_n_rms_Td,
            case_id=case_id,
        )
        if periodic.diagnostics.get("converged") is not True:
            raise RuntimeError(
                "two_term time-periodic distribution did not converge; refusing "
                "to publish cycle-averaged EEDF, rates, or transport"
            )
        case = build_case_result(
            self.config,
            self.cross_sections,
            self.solver_config,
            e_over_n_rms_Td,
            case_id,
            periodic.operator_block,
            periodic.cycle_averaged_eedf_eV_inv,
            metadata=periodic.diagnostics,
            solver_name=self.name,
        )
        case.rf_phase = periodic.phase
        return case

    def solve_native_distribution(
        self,
        e_over_n_Td: float,
    ) -> NativeDistributionResult:
        return solve_stationary_distribution(
            self.config,
            self.cross_sections,
            self.solver_config,
            e_over_n_Td,
        )

    def solve_time_periodic_distribution(
        self,
        e_over_n_rms_Td: float,
        *,
        case_id: str = "rf_case",
    ) -> TimePeriodicDistributionResult:
        return solve_periodic_distribution(
            self.config,
            self.cross_sections,
            self.solver_config,
            e_over_n_rms_Td,
            case_id=case_id,
            solver_name=self.name,
        )
