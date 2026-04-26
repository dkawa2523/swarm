"""Lightweight multi-term Boltzmann adapter.

This module provides the merge-safe multi-term entry point used by the unified
runner.  The current implementation reuses the validated native two-term
finite-volume block as the l=0/l=1 base and annotates the result with explicit
multi-term approximation metadata.  This keeps the public result/output contract
stable while leaving room for a future block-coupled finite-volume l>=2 solver.
"""

from __future__ import annotations

from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.solvers.boltzmann_two_term import BoltzmannTwoTermSolver


class MultiTermBoltzmannSolver(BoltzmannTwoTermSolver):
    """Minimal multi-term-compatible solver facade.

    The implementation intentionally avoids a second divergent Boltzmann code
    path.  It uses the native two-term solver as the parity baseline and records
    that the current closure is lmax=1-equivalent.  Downstream diagnostics,
    benchmark comparison, transport long-table output, state-resolved
    superelastic generation, and e-e relaxation therefore work identically for
    two-term and multi-term runs.
    """

    name = "multiterm_boltzmann"

    def __init__(self, config: SwarmConfig, cross_sections: CrossSectionSet) -> None:
        super().__init__(config, cross_sections)
        self.backend = "native_bolsig"

    def solve_case(self, e_over_n_Td: float, case_id: str) -> SwarmCaseResult:
        case = self._solve_case_native(e_over_n_Td, case_id)
        case.solver = self.name
        case.metadata.setdefault("multiterm_status", "lmax1_adapter")
        case.metadata.setdefault("multiterm_lmax_effective", 1)
        case.metadata.setdefault(
            "multiterm_notes",
            "Uses validated two-term finite-volume block until block-coupled l>=2 solver is enabled.",
        )
        return case
