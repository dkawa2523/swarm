"""Positive Perron solve for the stationary response propagator."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq
from scipy.sparse.linalg import ArpackNoConvergence, LinearOperator, eigs

from electron_swarm.solvers.propagator.models import (
    PropagatorConvergenceDiagnostics,
    PropagatorSteadySolution,
)
from electron_swarm.solvers.propagator.observables import (
    flux_drift_velocity_m_s,
    mean_energy_eV,
    normalized_population,
)
from electron_swarm.solvers.propagator.operator import (
    PropagatorOperator,
    RESPONSE_BUILD_BLAS_THREADS,
)


class PropagatorConvergenceError(RuntimeError):
    def __init__(
        self,
        message: str,
        diagnostics: PropagatorConvergenceDiagnostics,
    ) -> None:
        super().__init__(message)
        self.diagnostics = diagnostics


@dataclass(frozen=True, slots=True)
class SteadyControls:
    max_iterations: int
    tolerance: float
    tail_probability_target: float
    tail_energy_start_eV: float | None = None


@dataclass(slots=True)
class _WorkBudget:
    maximum: int
    used: int = 0

    def consume(self) -> None:
        if self.used >= self.maximum:
            raise _SweepBudgetExceeded
        self.used += 1


@dataclass(frozen=True, slots=True)
class _PerronEvaluation:
    growth_s_inv: float
    spectral_radius: float
    state: np.ndarray
    response_residual_L1: float
    system_build_time_s: float
    system_memory_bytes: int
    outer_escape_s_inv: float
    boundary_residual_L1: float
    krylov_dimension: int
    ritz_pairs: int
    cone_compatible_ritz_pairs: int
    perron_real_part_gap: float
    shell_rational_error_estimate: float
    shell_coefficient_error_estimate: float
    shell_maximum_coefficient_segments: int


class _SweepBudgetExceeded(RuntimeError):
    pass


def maxwellian_initial_population(
    operator: PropagatorOperator,
    *,
    electron_temperature_eV: float | None = None,
) -> np.ndarray:
    grid = operator.grid
    if electron_temperature_eV is None:
        energy_mass = operator.collisions.elastic_equilibrium_mass.copy()
    else:
        temperature = max(float(electron_temperature_eV), 1.0e-5)
        eedf = np.sqrt(np.maximum(grid.energy_centers_eV, 0.0)) * np.exp(
            -grid.energy_centers_eV / temperature
        )
        energy_mass = eedf * grid.energy_widths_eV
        energy_mass /= float(np.sum(energy_mass))
    return energy_mass[:, None] * operator.collisions.isotropic_weights[None, :]


def krylov_subspace_dimension(cell_count: int) -> int:
    return min(
        256,
        max(48, int(np.sqrt(max(int(cell_count), 1)))),
        max(4, int(cell_count) - 1),
    )


def _base_krylov_subspace_dimension(cell_count: int) -> int:
    return min(
        128,
        max(48, int(np.sqrt(max(int(cell_count), 1)))),
        max(4, int(cell_count) - 1),
    )


def _response_krylov_subspace_dimension(
    operator: PropagatorOperator,
    acceleration_m_s2: float,
) -> int:
    """Size Arnoldi from the characteristic optical depth at the origin."""

    cells = operator.grid.cells
    base = _base_krylov_subspace_dimension(cells)
    expanded = krylov_subspace_dimension(cells)
    if expanded == base:
        return base
    origin_frequency = float(
        sum(
            item.collision_frequency_s_inv[0]
            for item in operator.collisions.elastic
        )
    )
    grazing_cosine = float(np.min(np.abs(operator.grid.mu_centers)))
    optical_depth = (
        origin_frequency
        * float(operator.grid.speed_edges_m_s[1])
        / max(float(acceleration_m_s2) * grazing_cosine, 1.0e-300)
    )
    angular_cluster_threshold = max(
        16.0,
        0.25 * float(operator.grid.polar_cells),
    )
    return expanded if optical_depth > angular_cluster_threshold else base


def _positive_perron_evaluation(
    operator: PropagatorOperator,
    acceleration_m_s2: float,
    growth_s_inv: float,
    shift_s_inv: float,
    initial: np.ndarray,
    tolerance: float,
    budget: _WorkBudget,
) -> _PerronEvaluation:
    system = operator.build_response_system(
        acceleration_m_s2,
        growth_s_inv,
        shift_s_inv,
    )
    shape = initial.shape
    cells = initial.size

    def action(vector: np.ndarray) -> np.ndarray:
        budget.consume()
        return system.apply(np.asarray(vector, dtype=float).reshape(shape)).population.reshape(-1)

    linear = LinearOperator((cells, cells), matvec=action, dtype=float)
    state = normalized_population(initial).reshape(-1)
    # A few genuine positive iterations select the Perron communicating class
    # before Arnoldi accelerates the remaining spectral separation.
    for _ in range(3):
        mapped = action(state)
        mass = float(np.sum(mapped))
        if mass <= 0.0 or np.any(~np.isfinite(mapped)) or np.min(mapped) < -1.0e-13 * max(float(np.max(np.abs(mapped))), 1.0):
            raise FloatingPointError("positive response map produced a nonpositive iterate")
        state = mapped / mass
    positive_seed = state.copy()
    krylov_dimension = _response_krylov_subspace_dimension(
        operator,
        acceleration_m_s2,
    )
    eigensolver_tolerance = min(
        max(tolerance * 0.001, 1.0e-12),
        1.0e-9,
    )
    maximum_ritz_pairs = min(4, max(1, cells - 2))

    def arnoldi(
        requested_pairs: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        try:
            return eigs(
                linear,
                k=requested_pairs,
                which="LR",
                v0=state,
                tol=eigensolver_tolerance,
                maxiter=max(2, budget.maximum - budget.used),
                ncv=krylov_dimension,
            )
        except _SweepBudgetExceeded:
            raise
        except ArpackNoConvergence as exc:
            raise RuntimeError(
                "positive response Arnoldi iteration did not converge"
            ) from exc

    def cone_compatible_candidates(
        values: np.ndarray,
        vectors: np.ndarray,
    ) -> list[tuple[float, int, complex, np.ndarray, float]]:
        candidates: list[
            tuple[float, int, complex, np.ndarray, float]
        ] = []
        for index, raw_eigenvalue in enumerate(values):
            eigenvalue = complex(raw_eigenvalue)
            vector = np.asarray(vectors[:, index])
            imaginary_fraction = float(
                np.linalg.norm(np.imag(vector), ord=1)
                / max(np.linalg.norm(vector, ord=1), 1.0e-300)
            )
            if (
                not np.isfinite(eigenvalue.real)
                or not np.isfinite(eigenvalue.imag)
                or eigenvalue.real <= 0.0
                or abs(eigenvalue.imag)
                > 1.0e-8 * max(abs(eigenvalue.real), 1.0)
                or imaginary_fraction > 1.0e-8
            ):
                continue
            real_vector = np.real(vector)
            if float(np.sum(real_vector)) < 0.0:
                real_vector = -real_vector
            absolute_mass = max(
                float(np.sum(np.abs(real_vector))),
                1.0e-300,
            )
            negative_fraction = float(
                np.sum(np.maximum(-real_vector, 0.0))
            ) / absolute_mass
            if negative_fraction <= 1.0e-6:
                candidates.append(
                    (
                        float(eigenvalue.real),
                        index,
                        eigenvalue,
                        real_vector,
                        negative_fraction,
                    )
                )
        candidates.sort(key=lambda item: item[0], reverse=True)
        return candidates

    # A cone-compatible real Ritz vector is already a Perron certificate for
    # the positive response map; requesting several clustered slow modes in
    # every scalar-root evaluation can consume the complete finite budget.
    # Resolve one principal candidate first and expand the same Arnoldi
    # subspace only when the real/nonnegative certificate is absent.
    ritz_pairs = 1
    eigenvalues, eigenvectors = arnoldi(ritz_pairs)
    cone_candidates = cone_compatible_candidates(
        eigenvalues,
        eigenvectors,
    )
    if not cone_candidates and maximum_ritz_pairs > 1:
        ritz_pairs = maximum_ritz_pairs
        eigenvalues, eigenvectors = arnoldi(ritz_pairs)
        cone_candidates = cone_compatible_candidates(
            eigenvalues,
            eigenvectors,
        )
    if not cone_candidates:
        raise RuntimeError(
            "positive response Arnoldi subspace contained no real "
            "cone-compatible Perron vector"
        )
    _, selected_index, eigenvalue, state, negative_fraction = (
        cone_candidates[0]
    )
    largest_returned_real_part = float(np.max(np.real(eigenvalues)))
    if largest_returned_real_part - eigenvalue.real > max(
        20.0 * eigensolver_tolerance * max(abs(eigenvalue.real), 1.0),
        1.0e-11,
    ):
        raise RuntimeError(
            "cone-compatible response mode is not the largest-real Ritz mode"
        )
    if len(cone_candidates) > 1:
        second_cone_eigenvalue = cone_candidates[1][0]
        if eigenvalue.real - second_cone_eigenvalue <= max(
            20.0
            * eigensolver_tolerance
            * max(abs(eigenvalue.real), 1.0),
            1.0e-11,
        ):
            raise RuntimeError(
                "positive response has an unresolved non-unique Perron branch"
            )
    other_real_parts = [
        float(np.real(value))
        for index, value in enumerate(eigenvalues)
        if index != selected_index and np.isfinite(np.real(value))
    ]
    perron_real_part_gap = (
        float("inf")
        if not other_real_parts
        else (
            float(eigenvalue.real) - max(other_real_parts)
        )
        / max(abs(float(eigenvalue.real)), 1.0e-300)
    )

    # Applying the positive physical map is a Perron iteration, not clipping.
    # It removes roundoff-scale signed Ritz components while retaining the
    # eigenproblem itself.
    state /= float(np.sum(state))
    if float(np.min(state)) < 0.0:
        interior = 0.5 * positive_seed + 0.5 / cells
        interior /= float(np.sum(interior))
        active = state < 0.0
        denominator = interior[active] - state[active]
        if np.any(denominator <= 0.0):
            raise RuntimeError("Perron cone line search has no positive direction")
        fraction = float(np.max(-state[active] / denominator))
        fraction = min(1.0, fraction * (1.0 + 1.0e-10) + 1.0e-15)
        state = (1.0 - fraction) * state + fraction * interior
        if float(np.min(state)) < -2.0e-15:
            raise RuntimeError("Perron cone line search did not enter the cone")
    radius = float(eigenvalue.real)
    response_residual = float("inf")
    application = None
    for _ in range(8):
        budget.consume()
        application = system.apply(state.reshape(shape))
        mapped = application.population.reshape(-1)
        radius = float(np.sum(mapped))
        if radius <= 0.0 or not np.isfinite(radius):
            raise FloatingPointError("positive response has a nonpositive spectral radius")
        candidate = mapped / radius
        scale = max(float(np.max(np.abs(candidate))), 1.0)
        if float(np.min(candidate)) < -1.0e-8 * scale:
            raise RuntimeError(
                "positive response iteration left the nonnegative cone: "
                f"growth={growth_s_inv:.9e}, minimum={float(np.min(candidate)):.3e}, "
                f"scale={scale:.3e}"
            )
        response_residual = float(np.sum(np.abs(candidate - state)))
        state = candidate
        if response_residual <= max(tolerance * 0.1, 2.0e-11):
            break
    assert application is not None
    budget.consume()
    final_application = system.apply(state.reshape(shape))
    mapped = final_application.population.reshape(-1)
    radius = float(np.sum(mapped))
    final_negative = float(np.sum(np.maximum(-mapped, 0.0)))
    if final_negative > 1.0e-10 * max(float(np.sum(np.abs(mapped))), 1.0e-300):
        raise RuntimeError("positive response Perron polish retained negative mass")
    response_residual = float(np.sum(np.abs(mapped - radius * state))) / max(
        abs(radius), 1.0e-300
    )
    return _PerronEvaluation(
        growth_s_inv=float(growth_s_inv),
        spectral_radius=radius,
        state=state.reshape(shape),
        response_residual_L1=response_residual,
        system_build_time_s=system.build_time_s,
        system_memory_bytes=system.memory_bytes,
        outer_escape_s_inv=final_application.outer_escape_s_inv,
        boundary_residual_L1=final_application.boundary_residual_L1,
        krylov_dimension=krylov_dimension,
        ritz_pairs=ritz_pairs,
        cone_compatible_ritz_pairs=len(cone_candidates),
        perron_real_part_gap=perron_real_part_gap,
        shell_rational_error_estimate=(
            system.maximum_shell_rational_error_estimate
        ),
        shell_coefficient_error_estimate=(
            system.maximum_shell_coefficient_error_estimate
        ),
        shell_maximum_coefficient_segments=(
            system.maximum_shell_coefficient_segments
        ),
    )


def _diagnostics(
    operator: PropagatorOperator,
    evaluation: _PerronEvaluation,
    budget: _WorkBudget,
    controls: SteadyControls,
) -> PropagatorConvergenceDiagnostics:
    state = evaluation.state
    growth = evaluation.growth_s_inv
    reaction = operator.reaction_number_rate_s_inv(state)
    number_defect = reaction - evaluation.outer_escape_s_inv - growth
    angular_frequency = np.zeros(operator.grid.energy_cells, dtype=float)
    for transfer in operator.collisions.elastic:
        angular_frequency += transfer.collision_frequency_s_inv
    activity = float(
        np.sum(
            (
                angular_frequency
                + operator.collisions.nonlocal_outflow_s_inv
            )[:, None]
            * state
        )
    )
    activity = max(activity + abs(growth), 1.0)
    tail_start = controls.tail_energy_start_eV
    if tail_start is None:
        tail_start = 0.8 * float(operator.grid.energy_edges_eV[-1])
    tail = float(
        np.sum(
            np.abs(
                state[
                    operator.grid.energy_centers_eV >= float(tail_start), :
                ]
            )
        )
    )
    negative = float(np.sum(np.maximum(-state, 0.0)))
    return PropagatorConvergenceDiagnostics(
        converged=False,
        iterations=budget.used,
        shape_change_L1=evaluation.response_residual_L1,
        growth_relative_change=abs(evaluation.spectral_radius - 1.0),
        mean_energy_relative_change=0.0,
        drift_relative_change=0.0,
        operator_residual_L1=evaluation.response_residual_L1,
        normalization_error=abs(float(np.sum(state)) - 1.0),
        negative_population_mass=negative,
        tail_probability=tail,
        outer_acceleration_flux_fraction=(
            evaluation.outer_escape_s_inv / activity
        ),
        number_balance_residual=abs(number_defect) / activity,
        stop_reason="iterating",
        timings_s={
            "response_system_build_s": evaluation.system_build_time_s,
            "response_boundary_residual_L1": evaluation.boundary_residual_L1,
            "response_spectral_radius_error": abs(
                evaluation.spectral_radius - 1.0
            ),
            "response_operator_memory_bytes": float(
                evaluation.system_memory_bytes
            ),
            "response_build_blas_threads": float(
                RESPONSE_BUILD_BLAS_THREADS
            ),
            "response_krylov_dimension": float(evaluation.krylov_dimension),
            "response_ritz_pairs": float(evaluation.ritz_pairs),
            "response_cone_compatible_ritz_pairs": float(
                evaluation.cone_compatible_ritz_pairs
            ),
            "response_perron_real_part_gap": (
                evaluation.perron_real_part_gap
            ),
            "response_shell_rational_error_estimate": (
                evaluation.shell_rational_error_estimate
            ),
            "response_shell_coefficient_error_estimate": (
                evaluation.shell_coefficient_error_estimate
            ),
            "response_shell_maximum_coefficient_segments": float(
                evaluation.shell_maximum_coefficient_segments
            ),
            "mean_energy_eV": mean_energy_eV(operator.grid, state),
            "drift_velocity_m_s": flux_drift_velocity_m_s(
                operator.grid, state
            ),
        },
    )


def solve_steady_state(
    operator: PropagatorOperator,
    acceleration_m_s2: float,
    controls: SteadyControls,
    *,
    initial_population: np.ndarray | None = None,
) -> PropagatorSteadySolution:
    """Solve the nonlinear response eigencondition on its positive branch."""

    if controls.max_iterations < 1:
        raise ValueError("max_iterations must be positive")
    tolerance = float(controls.tolerance)
    if tolerance <= 0.0 or not np.isfinite(tolerance):
        raise ValueError("convergence tolerance must be finite and positive")
    shape = (operator.grid.energy_cells, operator.grid.polar_cells)
    if initial_population is None:
        initial = maxwellian_initial_population(operator)
    else:
        initial = normalized_population(initial_population)
        if initial.shape != shape:
            raise ValueError("warm-start population shape does not match grid")
    budget = _WorkBudget(int(controls.max_iterations))
    bounds = operator.growth_bounds(acceleration_m_s2)
    # Electron-number balance is a scalar conservation law and can be solved
    # substantially more tightly than the requested shape tolerance.  Leaving
    # the response radius only at O(tolerance) makes the birth/escape balance
    # depend on the very large elastic activity scale at high field.
    root_target = max(tolerance * 0.001, 2.0e-12)

    def growth_balance_defect(evaluation: _PerronEvaluation) -> float:
        return float(
            operator.reaction_number_rate_s_inv(evaluation.state)
            - evaluation.outer_escape_s_inv
            - evaluation.growth_s_inv
        )

    def growth_balance_target(evaluation: _PerronEvaluation) -> float:
        reaction = operator.reaction_number_rate_s_inv(evaluation.state)
        scale = max(
            abs(reaction),
            abs(evaluation.outer_escape_s_inv),
            abs(evaluation.growth_s_inv),
            1.0,
        )
        # Keep the scalar Schur solve one decade tighter than the public
        # normalized number-balance gate after the different reaction and
        # total-collision activity scales are accounted for.
        return max(tolerance * 0.001, 1.0e-12) * scale

    def evaluation_score(evaluation: _PerronEvaluation) -> float:
        return max(
            abs(evaluation.spectral_radius - 1.0) / root_target,
            abs(growth_balance_defect(evaluation))
            / growth_balance_target(evaluation),
        )

    def nonlinear_fixed_point_converged(
        evaluation: _PerronEvaluation,
    ) -> bool:
        """Accept an interior or physical-boundary response fixed point.

        A number-conserving model has its exact growth root at the physical
        upper bound ``g=0``.  Such a root cannot have a sign-changing bracket
        outside the admissible interval.  It is nevertheless a genuine root
        when both the Perron equation and the independent number balance are
        satisfied to the solver tolerance.  Keep the tighter ``root_target``
        for resolving bracketed interior roots, but use the public nonlinear
        tolerance for this endpoint existence test.
        """

        fixed_point_tolerance = max(tolerance, 1.0e-8)
        return bool(
            abs(evaluation.spectral_radius - 1.0)
            <= fixed_point_tolerance
            and evaluation.response_residual_L1
            <= fixed_point_tolerance
            and abs(growth_balance_defect(evaluation))
            <= growth_balance_target(evaluation)
        )

    try:
        seed_growth = float(
            np.clip(
                operator.reaction_number_rate_s_inv(initial),
                bounds.lower_s_inv,
                bounds.upper_s_inv,
            )
        )
        seed = _positive_perron_evaluation(
            operator,
            acceleration_m_s2,
            seed_growth,
            operator.minimum_response_shift(seed_growth),
            initial,
            tolerance,
            budget,
        )
        evaluations = [seed]
        best = seed
        macro_state = seed

        def bracket_from_evaluations() -> tuple[
            _PerronEvaluation | None, _PerronEvaluation | None
        ]:
            positive = [
                item
                for item in evaluations
                if growth_balance_defect(item) >= 0.0
            ]
            negative = [
                item
                for item in evaluations
                if growth_balance_defect(item) <= 0.0
            ]
            lower_item = max(
                positive, key=lambda item: item.growth_s_inv, default=None
            )
            upper_item = min(
                negative, key=lambda item: item.growth_s_inv, default=None
            )
            if (
                lower_item is not None
                and upper_item is not None
                and lower_item.growth_s_inv <= upper_item.growth_s_inv
            ):
                return lower_item, upper_item
            return None, None

        macro_secant_steps = 0
        for _ in range(8):
            if nonlinear_fixed_point_converged(best):
                break
            fallback_growth = float(
                np.clip(
                    macro_state.growth_s_inv
                    + growth_balance_defect(macro_state),
                    bounds.lower_s_inv,
                    bounds.upper_s_inv,
                )
            )
            balance_growth = fallback_growth
            if len(evaluations) >= 2:
                previous = evaluations[-2]
                delta_growth = (
                    macro_state.growth_s_inv - previous.growth_s_inv
                )
                delta_defect = (
                    growth_balance_defect(macro_state)
                    - growth_balance_defect(previous)
                )
                if delta_growth != 0.0:
                    slope = delta_defect / delta_growth
                    if np.isfinite(slope) and slope < 0.0:
                        proposed = float(
                            macro_state.growth_s_inv
                            - growth_balance_defect(macro_state) / slope
                        )
                        lower, upper = bracket_from_evaluations()
                        proposed_lower = (
                            bounds.lower_s_inv
                            if lower is None
                            else lower.growth_s_inv
                        )
                        proposed_upper = (
                            bounds.upper_s_inv
                            if upper is None
                            else upper.growth_s_inv
                        )
                        if (
                            np.isfinite(proposed)
                            and proposed_lower <= proposed <= proposed_upper
                        ):
                            balance_growth = proposed
                            macro_secant_steps += 1
            if abs(balance_growth - macro_state.growth_s_inv) <= (
                growth_balance_target(macro_state)
            ):
                break
            balanced = _positive_perron_evaluation(
                operator,
                acceleration_m_s2,
                balance_growth,
                operator.minimum_response_shift(balance_growth),
                macro_state.state,
                tolerance,
                budget,
            )
            evaluations.append(balanced)
            macro_state = balanced
            if evaluation_score(balanced) < evaluation_score(best):
                best = balanced

        lower, upper = bracket_from_evaluations()
        expansion = max(
            abs(
                evaluations[-1].growth_s_inv
                - evaluations[0].growth_s_inv
            ),
            1.0,
        )
        for _ in range(12):
            if nonlinear_fixed_point_converged(best):
                break
            lower, upper = bracket_from_evaluations()
            if lower is not None and upper is not None:
                break
            if all(
                growth_balance_defect(item) > 0.0
                for item in evaluations
            ):
                anchor = max(evaluations, key=lambda item: item.growth_s_inv)
                candidate_growth = min(
                    bounds.upper_s_inv,
                    anchor.growth_s_inv + expansion,
                )
            elif all(
                growth_balance_defect(item) < 0.0
                for item in evaluations
            ):
                anchor = min(evaluations, key=lambda item: item.growth_s_inv)
                candidate_growth = max(
                    bounds.lower_s_inv,
                    anchor.growth_s_inv - expansion,
                )
            else:
                break
            if any(
                candidate_growth == item.growth_s_inv for item in evaluations
            ):
                break
            candidate = _positive_perron_evaluation(
                operator,
                acceleration_m_s2,
                candidate_growth,
                operator.minimum_response_shift(candidate_growth),
                anchor.state,
                tolerance,
                budget,
            )
            evaluations.append(candidate)
            if evaluation_score(candidate) < evaluation_score(best):
                best = candidate
            expansion *= 4.0
    except _SweepBudgetExceeded as exc:
        diagnostics = PropagatorConvergenceDiagnostics(
            iterations=budget.used,
            stop_reason="maximum_response_applications",
        )
        raise PropagatorConvergenceError(
            "propagator exceeded its finite response-application budget",
            diagnostics,
        ) from exc
    except (FloatingPointError, RuntimeError, ValueError) as exc:
        diagnostics = PropagatorConvergenceDiagnostics(
            iterations=budget.used,
            stop_reason="perron_branch_failure",
        )
        raise PropagatorConvergenceError(str(exc), diagnostics) from exc

    lower, upper = bracket_from_evaluations()
    if (
        not nonlinear_fixed_point_converged(best)
        and (lower is None or upper is None)
    ):
        diagnostics = _diagnostics(
            operator,
            best,
            budget,
            controls,
        )
        diagnostics.stop_reason = "no_normalizable_principal_mode"
        raise PropagatorConvergenceError(
            "propagator physical growth bounds do not bracket a normalizable mode",
            diagnostics,
        )

    root_iterations = 0
    if not nonlinear_fixed_point_converged(best):
        assert lower is not None and upper is not None
        evaluation_by_growth = {
            item.growth_s_inv: item for item in evaluations
        }

        def evaluate_number_balance(candidate_growth: float) -> float:
            nonlocal best, root_iterations
            growth = float(candidate_growth)
            cached = evaluation_by_growth.get(growth)
            if cached is not None:
                return growth_balance_defect(cached)
            warm = min(
                evaluations,
                key=lambda item: abs(item.growth_s_inv - growth),
            ).state
            candidate = _positive_perron_evaluation(
                operator,
                acceleration_m_s2,
                growth,
                operator.minimum_response_shift(growth),
                warm,
                tolerance,
                budget,
            )
            evaluations.append(candidate)
            evaluation_by_growth[growth] = candidate
            root_iterations += 1
            if evaluation_score(candidate) < evaluation_score(best):
                best = candidate
            return growth_balance_defect(candidate)

        lower_defect = growth_balance_defect(lower)
        upper_defect = growth_balance_defect(upper)
        if lower_defect < 0.0 or upper_defect > 0.0:
            diagnostics = _diagnostics(operator, best, budget, controls)
            diagnostics.stop_reason = "nonmonotone_number_balance_bracket"
            raise PropagatorConvergenceError(
                "propagator number-balance root lost its physical bracket",
                diagnostics,
            )
        balance_absolute_tolerance = max(
            0.05
            * min(
                growth_balance_target(lower),
                growth_balance_target(upper),
            ),
            1.0e-8,
        )
        try:
            root_growth = float(
                brentq(
                    evaluate_number_balance,
                    lower.growth_s_inv,
                    upper.growth_s_inv,
                    xtol=balance_absolute_tolerance,
                    rtol=4.0 * np.finfo(float).eps,
                    maxiter=32,
                    disp=True,
                )
            )
            evaluate_number_balance(root_growth)
        except _SweepBudgetExceeded as exc:
            diagnostics = _diagnostics(operator, best, budget, controls)
            diagnostics.stop_reason = "maximum_response_applications"
            raise PropagatorConvergenceError(
                "propagator exceeded its finite response-application budget",
                diagnostics,
            ) from exc
        except (FloatingPointError, RuntimeError, ValueError) as exc:
            diagnostics = _diagnostics(operator, best, budget, controls)
            diagnostics.stop_reason = "perron_branch_failure"
            raise PropagatorConvergenceError(str(exc), diagnostics) from exc

        lower, upper = bracket_from_evaluations()
    diagnostics = _diagnostics(operator, best, budget, controls)
    diagnostics.timings_s["growth_root_iterations"] = float(root_iterations)
    diagnostics.timings_s["growth_response_system_builds"] = float(
        len(evaluations)
    )
    diagnostics.timings_s["growth_secant_steps"] = float(
        macro_secant_steps
    )
    diagnostics.timings_s["growth_bracket_width_s_inv"] = float(
        0.0
        if lower is None or upper is None
        else upper.growth_s_inv - lower.growth_s_inv
    )
    core_ok = (
        abs(best.spectral_radius - 1.0) <= max(tolerance, 1.0e-8)
        and best.response_residual_L1 <= max(tolerance, 1.0e-8)
        and diagnostics.normalization_error <= 1.0e-12
        and diagnostics.negative_population_mass <= 1.0e-12
        and diagnostics.number_balance_residual <= max(tolerance, 1.0e-8)
        and best.boundary_residual_L1 <= 1.0e-11
    )
    if not core_ok:
        diagnostics.stop_reason = "response_eigenproblem_not_converged"
        raise PropagatorConvergenceError(
            "propagator stationary response eigenproblem did not converge",
            diagnostics,
        )
    tail_ok = (
        diagnostics.tail_probability <= controls.tail_probability_target
        and diagnostics.outer_acceleration_flux_fraction
        <= max(controls.tail_probability_target, 1.0e-10)
    )
    if not tail_ok:
        diagnostics.stop_reason = "non_normalizable_or_truncated_tail"
        raise PropagatorConvergenceError(
            "propagator outer escape does not meet the finite-domain bound",
            diagnostics,
        )
    diagnostics.converged = True
    diagnostics.stop_reason = "converged"
    return PropagatorSteadySolution(
        population=best.state,
        growth_frequency_s_inv=best.growth_s_inv,
        diagnostics=diagnostics,
    )
