"""Reciprocal finite-volume angular collision kernels.

The angular transition matrix acts on cell populations.  Its columns are
conditional destination probabilities.  For a rotational phase law the
isotropic source--destination joint measure is symmetric.  A finite-volume
maximum-entropy projection constructs the discrete law inside the conservative
and reciprocal manifold and gives the requested P1 relaxation exactly.  This
is part of collision-operator assembly; no solution or output is corrected.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares
from scipy.special import i0e

from electron_swarm.physics.angular_scattering import inverse_langevin
from electron_swarm.solvers.propagator.models import PropagatorGrid


def isotropic_destination_weights(grid: PropagatorGrid) -> np.ndarray:
    weights = np.asarray(grid.solid_angles_sr, dtype=float) / (4.0 * np.pi)
    if not np.isclose(float(np.sum(weights)), 1.0, rtol=0.0, atol=1.0e-14):
        raise FloatingPointError("isotropic angular weights do not sum to one")
    return weights


def _cell_quadrature(
    grid: PropagatorGrid,
    order: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    local_nodes, local_weights = np.polynomial.legendre.leggauss(order)
    cosine_edges = np.cos(grid.theta_edges_rad)
    nodes = np.empty(grid.polar_cells * order, dtype=float)
    weights = np.empty_like(nodes)
    groups = np.repeat(np.arange(grid.polar_cells, dtype=np.int32), order)
    for index in range(grid.polar_cells):
        lower = float(cosine_edges[index + 1])
        upper = float(cosine_edges[index])
        section = slice(index * order, (index + 1) * order)
        nodes[section] = (
            0.5 * (upper - lower) * local_nodes
            + 0.5 * (upper + lower)
        )
        weights[section] = 0.5 * (upper - lower) * local_weights
    return nodes, weights, groups


def _maxent_joint_measure(
    grid: PropagatorGrid,
    mean_cosine: float,
    *,
    cell_quadrature_order: int,
) -> np.ndarray:
    """Integrate the continuous maxent law over both angular cells.

    The azimuthal integral is analytic.  The remaining two cosine integrals
    use the same Gauss rule in source and destination cells, making the raw
    joint measure symmetric to roundoff.
    """

    value = float(mean_cosine)
    weights = isotropic_destination_weights(grid)
    if abs(value) <= 1.0e-14:
        return np.outer(weights, weights)

    kappa = float(inverse_langevin(value))
    mu, quadrature, groups = _cell_quadrature(
        grid,
        cell_quadrature_order,
    )
    transverse = np.sqrt(np.maximum(1.0 - mu * mu, 0.0))
    magnitude = abs(kappa)
    argument = magnitude * np.outer(transverse, transverse)
    if magnitude < 350.0:
        log_two_sinh = magnitude + np.log1p(-np.exp(-2.0 * magnitude))
    else:
        log_two_sinh = magnitude
    log_density = (
        np.log(magnitude)
        - log_two_sinh
        + kappa * np.outer(mu, mu)
        + np.log(i0e(argument))
        + argument
    )
    point_joint = (
        0.5
        * quadrature[:, None]
        * quadrature[None, :]
        * np.exp(log_density)
    )
    incidence = np.zeros((grid.polar_cells, mu.size), dtype=float)
    incidence[groups, np.arange(mu.size)] = 1.0
    joint = incidence @ point_joint @ incidence.T
    return 0.5 * (joint + joint.T)


def _project_reversible_p1_joint(
    base_joint: np.ndarray,
    grid: PropagatorGrid,
    mean_cosine: float,
    initial_dual: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the closest positive reciprocal joint with exact P0/P1 moments.

    The dual exponential family is the KL projection of the cell-integrated
    continuous maxent law.  A positive numerical floor represents entries that
    are mathematically nonzero but underflow for a strongly forward phase law.
    Continuation in ``mean_cosine`` keeps the solve on the same interior branch.
    """

    cells = grid.polar_cells
    mu = np.asarray(grid.mu_centers, dtype=float)
    weights = isotropic_destination_weights(grid)
    base = np.maximum(np.asarray(base_joint, dtype=float), np.finfo(float).tiny)
    if base.shape != (cells, cells):
        raise ValueError("maxent joint measure does not match the angular grid")
    base = 0.5 * (base + base.T)
    if cells % 2:
        raise ValueError("reversible maxent projection requires an even grid")
    half = cells // 2
    start = (
        np.zeros(cells, dtype=float)
        if initial_dual is None
        else np.asarray(initial_dual, dtype=float).copy()
    )
    if start.shape != (cells,):
        raise ValueError("maxent continuation state has the wrong size")

    def joint_from_dual(dual: np.ndarray) -> np.ndarray:
        # The prior and constraints are centrosymmetric.  Work directly in
        # that invariant subspace: alpha is even and beta is odd under
        # mu -> -mu.  This removes the redundant nonsymmetric dual directions
        # that otherwise drift when solutions are continued across many
        # closely spaced radial bands.
        alpha_half = dual[:half]
        beta_half = dual[half:]
        alpha = np.concatenate((alpha_half, alpha_half[::-1]))
        beta = np.concatenate((beta_half, -beta_half[::-1]))
        exponent = (
            alpha[:, None]
            + alpha[None, :]
            + beta[:, None] * mu[None, :]
            + mu[:, None] * beta[None, :]
        )
        return base * np.exp(np.clip(exponent, -700.0, 700.0))

    def residual(dual: np.ndarray) -> np.ndarray:
        joint = joint_from_dual(dual)
        margin = (np.sum(joint, axis=0) - weights) / weights
        moment = (
            mu @ joint - float(mean_cosine) * mu * weights
        ) / weights
        # Centrosymmetry makes the negative-half equations exact mirrors.
        return np.concatenate((margin[:half], moment[:half]))

    def jacobian(dual: np.ndarray) -> np.ndarray:
        joint = joint_from_dual(dual)
        transposed = joint.T
        margin = np.sum(joint, axis=0)
        first = mu @ joint
        second = (mu * mu) @ joint

        margin_alpha = transposed + np.diag(margin)
        moment_alpha = (
            transposed * mu[None, :] + np.diag(first)
        )
        margin_beta = transposed * mu[:, None] + np.diag(first)
        moment_beta = (
            transposed
            * mu[:, None]
            * mu[None, :]
            + np.diag(second)
        )
        reflected = np.arange(cells - 1, half - 1, -1)
        margin_rows = np.hstack(
            (
                margin_alpha[:, :half] + margin_alpha[:, reflected],
                margin_beta[:, :half] - margin_beta[:, reflected],
            )
        )
        moment_rows = np.hstack(
            (
                moment_alpha[:, :half] + moment_alpha[:, reflected],
                moment_beta[:, :half] - moment_beta[:, reflected],
            )
        )
        return np.vstack(
            (
                margin_rows[:half] / weights[:half, None],
                moment_rows[:half] / weights[:half, None],
            )
        )

    solution = least_squares(
        residual,
        start,
        jac=jacobian,
        xtol=1.0e-12,
        ftol=1.0e-12,
        gtol=1.0e-12,
        max_nfev=600,
        x_scale="jac",
    )
    scaled_error = float(np.max(np.abs(residual(solution.x))))
    if not np.isfinite(scaled_error) or scaled_error > 2.0e-9:
        raise FloatingPointError(
            "reversible maxent finite-volume projection did not converge: "
            f"scaled residual={scaled_error:.3e}"
        )
    return joint_from_dual(solution.x), solution.x


def _continuation_points(previous: float, target: float) -> tuple[float, ...]:
    anchors = (0.5, 0.8, 0.9, 0.95, 0.98, 0.99, 0.995, 0.997, 0.998)
    points = [
        value
        for value in anchors
        if previous + 1.0e-14 < value < target - 1.0e-14
    ]
    points.append(target)
    return tuple(points)


def maxent_p1_transition_matrices(
    grid: PropagatorGrid,
    mean_cosine: np.ndarray,
    *,
    cell_quadrature_order: int = 12,
) -> np.ndarray:
    """Build exact-P1 reversible kernels for an arbitrary vector of moments."""

    if cell_quadrature_order < 4:
        raise ValueError("angular cell quadrature order must be at least four")
    values = np.asarray(mean_cosine, dtype=float)
    if values.ndim != 1 or np.any(~np.isfinite(values)):
        raise ValueError("elastic first angular moments must be a finite vector")
    tolerance = 1.0e-12
    if np.any((values < -1.0 - tolerance) | (values > 1.0 + tolerance)):
        raise ValueError("elastic first angular moment is outside [-1, 1]")
    values = np.clip(values, -1.0, 1.0)
    result = np.empty(
        (values.size, grid.polar_cells, grid.polar_cells), dtype=float
    )
    weights = isotropic_destination_weights(grid)
    zero = np.outer(weights, np.ones(grid.polar_cells, dtype=float))
    identity = np.eye(grid.polar_cells, dtype=float)
    cache: dict[float, np.ndarray] = {0.0: zero, 1.0: identity}

    # Cell quadrature can leave otherwise identical cross-section ratios a few
    # ulps apart.  Sharing those kernels avoids redundant nonlinear solves;
    # the 12-decimal key is far below the enforced P1 moment tolerance.
    magnitude_keys = np.round(np.abs(values), decimals=12)
    requested = sorted(
        {
            float(key)
            for key in magnitude_keys
            if tolerance < float(key) < 1.0 - tolerance
        }
    )
    dual: np.ndarray | None = None
    previous = 0.0
    for target in requested:
        for point in _continuation_points(previous, target):
            base = _maxent_joint_measure(
                grid,
                point,
                cell_quadrature_order=cell_quadrature_order,
            )
            joint, dual = _project_reversible_p1_joint(
                base,
                grid,
                point,
                initial_dual=dual,
            )
            matrix = joint / weights[None, :]
            if point == target:
                cache[target] = matrix
            previous = point

    for index, requested_value in enumerate(values):
        value = float(requested_value)
        magnitude = float(abs(value))
        if magnitude <= tolerance:
            matrix = zero
        elif magnitude >= 1.0 - tolerance:
            matrix = identity
        else:
            matrix = cache[float(magnitude_keys[index])]
        if value < 0.0:
            matrix = matrix[::-1, :]
        result[index] = matrix

    mu = np.asarray(grid.mu_centers, dtype=float)
    for value, matrix in zip(values, result, strict=True):
        joint = matrix * weights[None, :]
        scale = max(float(np.max(matrix)), 1.0)
        if np.any(matrix < -2.0e-14 * scale):
            raise FloatingPointError("maxent-P1 angular kernel is negative")
        column_error = float(np.max(np.abs(np.sum(matrix, axis=0) - 1.0)))
        equilibrium_error = float(np.sum(np.abs(matrix @ weights - weights)))
        reciprocity_error = float(np.max(np.abs(joint - joint.T)))
        moment_error = float(np.max(np.abs(mu @ matrix - value * mu)))
        if column_error > 5.0e-9 or equilibrium_error > 5.0e-9:
            raise FloatingPointError(
                "maxent-P1 projection does not preserve probability"
            )
        if reciprocity_error > 5.0e-9:
            raise FloatingPointError("maxent-P1 angular kernel is not reciprocal")
        if moment_error > 5.0e-9:
            raise FloatingPointError(
                "maxent-P1 angular kernel does not realize momentum relaxation"
            )
    return result


def maxent_p1_transition_stack(
    grid: PropagatorGrid,
    mean_cosine: np.ndarray,
) -> np.ndarray:
    values = np.asarray(mean_cosine, dtype=float)
    if values.shape != (grid.energy_cells,):
        raise ValueError("mean cosine must match the propagator energy grid")
    return maxent_p1_transition_matrices(grid, values)
