"""Fail-closed stationary eigenmode solve for the full PN system."""

from __future__ import annotations

import warnings

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as spla

from electron_swarm.core.solver_configs import ConvergenceConfig
from electron_swarm.solvers.boltzmann_common.observables import (
    mean_energy_from_eedf,
    negative_mass_fraction,
    weighted_integral,
)

from .diagnostics import SolverDiagnostics
from .models import PNOperator, StationaryPNState


_NEGATIVE_F0_MASS_LIMIT = 1.0e-8


def _split(operator: PNOperator, values: np.ndarray) -> tuple[np.ndarray, ...]:
    return tuple(values[block] for block in operator.moment_slices)


def _moment_l1(operator: PNOperator, values: np.ndarray) -> float:
    blocks = _split(operator, values)
    return float(
        sum(
            np.sum(np.abs(block) * weight)
            for block, weight in zip(
                blocks,
                operator.moment_weights_eV,
                strict=True,
            )
        )
    )


def _weight_vector(operator: PNOperator) -> np.ndarray:
    return np.concatenate(operator.moment_weights_eV)


def _weighted_operator_l1_norm(operator: PNOperator) -> float:
    """Induced L1 norm for the staggered moment quadrature."""

    weights = _weight_vector(operator)
    column_sums = np.asarray(abs(operator.matrix).T @ weights).reshape(-1)
    return float(np.max(column_sums / weights))


def _full_relative_residual(
    operator: PNOperator,
    state: np.ndarray,
    growth_s_inv: float,
) -> float:
    """Weighted backward error of ``A f = growth f`` over every PN block."""

    residual_state = operator.matrix @ state - growth_s_inv * state
    state_norm = _moment_l1(operator, state)
    scale = (
        _weighted_operator_l1_norm(operator) + abs(growth_s_inv)
    ) * state_norm
    return _moment_l1(operator, residual_state) / max(scale, 1.0e-300)


def _normalized_linear_solve(
    operator: PNOperator,
    growth_s_inv: float,
) -> np.ndarray:
    """Solve the eigen-equations with one F0 row replaced by normalization."""

    dimension = operator.matrix.shape[0]
    shifted = (
        operator.matrix
        - growth_s_inv * sparse.identity(dimension, format="csr")
    ).tolil()
    n = len(operator.widths_eV)
    diagonal = np.abs(operator.matrix[:n, :n].diagonal())
    row = int(np.argmax(diagonal)) if np.any(diagonal > 0.0) else n - 1
    shifted[row, :] = 0.0
    shifted[row, operator.moment_slices[0]] = operator.widths_eV
    rhs = np.zeros(dimension, dtype=float)
    rhs[row] = 1.0
    with warnings.catch_warnings():
        warnings.simplefilter("error", spla.MatrixRankWarning)
        try:
            state = spla.spsolve(shifted.tocsr(), rhs)
        except (spla.MatrixRankWarning, RuntimeError) as exc:
            raise RuntimeError("stationary PN matrix is singular") from exc
    if state.shape != (dimension,) or not np.all(np.isfinite(state)):
        raise FloatingPointError("stationary PN solve returned non-finite values")
    f0 = state[operator.moment_slices[0]]
    normalization = weighted_integral(f0, operator.widths_eV)
    if not np.isfinite(normalization) or abs(normalization) < 1.0e-300:
        raise FloatingPointError("stationary PN F0 is not normalizable")
    if normalization < 0.0:
        state = -state
        normalization = -normalization
    return state / normalization


def solve_stationary_pn(
    operator: PNOperator,
    convergence: ConvergenceConfig,
) -> StationaryPNState:
    """Find and verify the normalized temporal eigenmode of all PN moments."""

    growth = 0.0
    previous: np.ndarray | None = None
    residual = np.inf
    shape_change = np.inf
    eigenvalue_change = np.inf
    state = np.empty(operator.matrix.shape[0], dtype=float)

    for iteration in range(1, convergence.max_iterations + 1):
        state = _normalized_linear_solve(operator, growth)
        applied = operator.matrix @ state
        applied_blocks = _split(operator, applied)
        new_growth = weighted_integral(
            applied_blocks[0], operator.widths_eV
        )
        residual = _full_relative_residual(operator, state, new_growth)
        if previous is not None:
            shape_change = _moment_l1(operator, state - previous) / max(
                _moment_l1(operator, state),
                1.0e-300,
            )
        eigenvalue_change = abs(new_growth - growth) / max(
            abs(new_growth), abs(growth), 1.0
        )

        if (
            previous is not None
            and shape_change <= convergence.tolerance
            and eigenvalue_change <= convergence.eigenvalue_tolerance
            and residual <= convergence.residual_tolerance
        ):
            growth = new_growth
            break
        previous = state.copy()
        growth = new_growth
    else:
        raise RuntimeError(
            "stationary PN solve did not converge: "
            f"iterations={convergence.max_iterations}, "
            f"full_residual={residual:.3e}, shape_change={shape_change:.3e}, "
            f"eigenvalue_change={eigenvalue_change:.3e}"
        )

    residual = _full_relative_residual(operator, state, growth)
    if not np.isfinite(residual) or residual > convergence.residual_tolerance:
        raise RuntimeError(
            "stationary PN full-system residual exceeds tolerance: "
            f"{residual:.3e} > {convergence.residual_tolerance:.3e}"
        )

    raw_coefficients = _split(operator, state)
    negative_mass = negative_mass_fraction(
        raw_coefficients[0], operator.widths_eV
    )
    if negative_mass > _NEGATIVE_F0_MASS_LIMIT:
        raise RuntimeError(
            "stationary PN converged to a nonphysical negative F0 mass "
            f"({negative_mass:.3e})"
        )
    coefficients = np.vstack(
        [
            projection @ coefficient
            for projection, coefficient in zip(
                operator.center_projections,
                raw_coefficients,
                strict=True,
            )
        ]
    )
    f0 = coefficients[0]
    tail_mask = operator.energy_eV >= 0.9 * operator.edges_eV[-1]
    diagnostics = SolverDiagnostics(
        iterations=iteration,
        growth_frequency_s_inv=float(growth),
        full_pn_relative_residual=float(residual),
        shape_change=float(shape_change),
        eigenvalue_change=float(eigenvalue_change),
        mean_energy_eV=mean_energy_from_eedf(
            operator.energy_eV,
            operator.widths_eV,
            f0,
        ),
        eedf_tail_fraction=float(
            np.sum(f0[tail_mask] * operator.widths_eV[tail_mask])
        ),
    )
    return StationaryPNState(
        coefficients=coefficients,
        raw_coefficients=raw_coefficients,
        diagnostics=diagnostics,
    )


__all__ = ["solve_stationary_pn"]
