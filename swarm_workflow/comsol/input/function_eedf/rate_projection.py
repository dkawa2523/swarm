"""Match solver rate anchors using the physical linear EEPF basis."""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from .contracts import CollisionRateKernel, FunctionEedfError
from .kernels import linear_rate_weights
from .moments import linear_moment_weights


SOURCE_RATE_SIGNIFICANCE_FRACTION = 1e-5
SOURCE_RATE_MATCH_TOLERANCE = 1e-7


def match_source_rates(
    energies: np.ndarray,
    means: np.ndarray,
    rows: np.ndarray,
    kernels: tuple[CollisionRateKernel, ...],
    anchor_indices: np.ndarray,
    anchor_rates: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Correct the shaped rows and return the maximum source-anchor residual.

    The degrees of freedom are probability contributions of nonnegative linear
    basis functions. Correction is a weighted least-distance projection on the
    original support. The scale for selecting negligible rates is never added
    to a rate or to the EEPF.
    """
    if not kernels or anchor_rates.shape != (len(anchor_indices), len(kernels)):
        raise FunctionEedfError("source reaction rates do not cover the EEDF anchors")
    if np.any(~np.isfinite(anchor_rates)) or np.any(anchor_rates < 0):
        raise FunctionEedfError(
            "source reaction coefficients must be finite and nonnegative"
        )
    mass_weights = linear_moment_weights(energies, 0.5)
    energy_per_mass = linear_moment_weights(energies, 1.5) / mass_weights
    rate_per_mass = (
        np.stack([linear_rate_weights(energies, kernel) for kernel in kernels])
        / mass_weights
    )
    original_rates = (rows * mass_weights) @ rate_per_mass.T
    floors = (
        np.maximum(anchor_rates.max(axis=0), original_rates.max(axis=0))
        * SOURCE_RATE_SIGNIFICANCE_FRACTION
    )
    # Keep the nonlinear tail development of the shaped EEDF between source
    # points. Only the multiplicative correction is interpolated in mean E.
    baseline_at_anchors = original_rates[anchor_indices]
    significant = np.maximum(anchor_rates, baseline_at_anchors) >= floors
    factors = np.ones_like(anchor_rates)
    np.divide(
        anchor_rates,
        baseline_at_anchors,
        out=factors,
        where=(baseline_at_anchors > 0) & significant,
    )
    factors[anchor_rates == 0] = 0.0
    correction = np.column_stack(
        [np.interp(means, means[anchor_indices], factor) for factor in factors.T]
    )
    targets = original_rates * correction
    targets[anchor_indices] = anchor_rates
    fitted = []
    for mean, row, target, estimate in zip(
        means, rows, targets, original_rates, strict=True
    ):
        relevant = (np.maximum(target, estimate) >= floors) & (
            np.maximum(target, estimate) > 0
        )
        for _ in range(len(kernels) + 1):
            positive = relevant & (target > 0)
            mass = row * mass_weights
            mass[np.any(rate_per_mass[target == 0] > 0, axis=0)] = 0.0
            operator = np.vstack(
                (
                    np.ones(len(energies)),
                    energy_per_mass / mean,
                    rate_per_mass[positive] / target[positive, None],
                )
            )
            corrected = _match_moments(mass, operator, float(mean))
            estimate = rate_per_mass @ corrected
            newly_relevant = (estimate >= floors) & (estimate > 0) & ~relevant
            if not np.any(newly_relevant):
                break
            relevant |= newly_relevant
        fitted.append(corrected / mass_weights)

    values = np.asarray(fitted)
    projected_rates = (values * mass_weights) @ rate_per_mass.T
    # Re-evaluate relevance after fitting to catch a newly significant rate.
    relevant = (np.maximum(targets, projected_rates) >= floors) & (
        np.maximum(targets, projected_rates) > 0
    )
    residual = np.zeros_like(targets)
    np.divide(
        np.abs(projected_rates - targets),
        np.maximum(targets, floors),
        out=residual,
        where=relevant,
    )
    maximum = float(residual.max())
    if maximum > SOURCE_RATE_MATCH_TOLERANCE:
        anchor, process = np.unravel_index(np.argmax(residual), residual.shape)
        raise FunctionEedfError(
            f"source-rate projection for {kernels[process].identity} at {means[anchor]:g} eV "
            f"has scaled residual {maximum:.3g}"
        )
    return values, float(residual[anchor_indices].max())


def _match_moments(mass: np.ndarray, operator: np.ndarray, mean: float) -> np.ndarray:
    """Solve the few-dimensional dual of the nonnegative quadratic projection."""
    residual = operator @ mass - 1.0
    if np.max(np.abs(residual)) <= 1e-10:
        return mass
    rms = np.sqrt(np.square(operator) @ mass)
    if np.any(~np.isfinite(rms)) or np.any(rms <= 0.0):
        raise FunctionEedfError(
            f"source rates at {mean:g} eV lie outside the EEDF support"
        )
    basis = operator / rms[:, None]

    def factors(dual: np.ndarray) -> np.ndarray:
        return 1.0 + dual @ basis

    def mismatch(dual: np.ndarray) -> np.ndarray:
        return operator @ (mass * np.maximum(factors(dual), 0.0)) - 1.0

    def jacobian(dual: np.ndarray) -> np.ndarray:
        return (operator * (mass * (factors(dual) > 0))) @ basis.T

    fit = least_squares(
        mismatch,
        np.zeros(len(operator)),
        jac=jacobian,
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-11,
        max_nfev=80,
    )
    if np.max(np.abs(fit.fun)) > 1e-9:
        raise FunctionEedfError(
            f"cannot match normalization, mean and source rates at {mean:g} eV"
        )
    return mass * np.maximum(factors(fit.x), 0.0)
