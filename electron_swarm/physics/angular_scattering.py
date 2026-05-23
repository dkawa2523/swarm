"""Angular-scattering moment models for product schema v2."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

AngularModelName = Literal["isotropic", "momentum_power", "maxent_p1"]
HigherMomentClosure = Literal["zero", "power", "maxent"]


def _energy_vector(energy_eV: np.ndarray | float) -> np.ndarray:
    energy = np.asarray(energy_eV, dtype=float).reshape(-1)
    if not np.all(np.isfinite(energy)):
        raise ValueError("energy_eV must be finite")
    return energy


def _base_moments(energy_eV: np.ndarray | float, lmax: int) -> np.ndarray:
    if lmax < 0:
        raise ValueError("lmax must be >= 0")
    energy = _energy_vector(energy_eV)
    out = np.zeros((lmax + 1, energy.size), dtype=float)
    out[0, :] = 1.0
    return out


def _cross_section_vector(
    value: np.ndarray | float | None,
    energy: np.ndarray,
    name: str,
) -> np.ndarray:
    if value is None:
        raise ValueError(f"{name} is required")
    vector = np.asarray(value, dtype=float).reshape(-1)
    if vector.shape != energy.shape:
        raise ValueError(f"{name} must match energy_eV shape")
    if not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must be finite")
    return vector


def _m1_from_cross_sections(
    energy_eV: np.ndarray | float,
    sigma_total: np.ndarray | float | None,
    sigma_momentum: np.ndarray | float | None,
) -> tuple[np.ndarray, np.ndarray]:
    energy = _energy_vector(energy_eV)
    total = _cross_section_vector(sigma_total, energy, "sigma_total")
    momentum = _cross_section_vector(sigma_momentum, energy, "sigma_momentum")
    if np.any(total <= 0.0):
        raise ValueError("sigma_total must be positive")
    m1 = 1.0 - momentum / total
    m1 = np.clip(m1, -1.0, 1.0)
    if not np.all(np.isfinite(m1)):
        raise ValueError("computed m1 must be finite")
    return energy, m1


def _finalize_moments(out: np.ndarray) -> np.ndarray:
    out = np.nan_to_num(out, nan=0.0, posinf=1.0, neginf=-1.0)
    out = np.clip(out, -1.0, 1.0)
    out[0, :] = 1.0
    return out


def _legendre_values(lmax: int, mu: np.ndarray) -> np.ndarray:
    values = np.zeros((lmax + 1, mu.size), dtype=float)
    values[0, :] = 1.0
    if lmax >= 1:
        values[1, :] = mu
    for ell in range(1, lmax):
        values[ell + 1, :] = (
            (2 * ell + 1) * mu * values[ell, :] - ell * values[ell - 1, :]
        ) / (ell + 1)
    return values


def _langevin(kappa: float) -> float:
    abs_kappa = abs(kappa)
    if abs_kappa < 1.0e-6:
        return kappa / 3.0 - kappa**3 / 45.0 + 2.0 * kappa**5 / 945.0
    if abs_kappa > 50.0:
        return float(np.sign(kappa) * (1.0 - 1.0 / abs_kappa))
    return float(1.0 / np.tanh(kappa) - 1.0 / kappa)


def _inverse_langevin(target: float) -> float:
    target = float(np.clip(target, -0.999999, 0.999999))
    if abs(target) < 1.0e-10:
        return 0.0
    sign = 1.0 if target > 0.0 else -1.0
    abs_target = abs(target)
    lo = 0.0
    hi = max(1.0, 1.0 / max(1.0e-12, 1.0 - abs_target))
    while _langevin(hi) < abs_target and hi < 1.0e6:
        hi *= 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if _langevin(mid) < abs_target:
            lo = mid
        else:
            hi = mid
    return sign * 0.5 * (lo + hi)


@dataclass(frozen=True, slots=True)
class IsotropicAngularModel:
    name: str = "isotropic"

    def metadata(self) -> dict[str, str]:
        return {
            "angular_model": self.name,
            "angular_moment_source": "isotropic_closure",
            "angular_higher_moment_closure": "zero",
            "angular_closure_assumption": "isotropic_scattering",
        }

    def moments(
        self,
        energy_eV: np.ndarray | float,
        lmax: int,
        *,
        sigma_total: np.ndarray | float | None = None,
        sigma_momentum: np.ndarray | float | None = None,
    ) -> np.ndarray:
        return _base_moments(energy_eV, lmax)


@dataclass(frozen=True, slots=True)
class MomentumPowerAngularModel:
    name: str = "momentum_power"

    def metadata(self) -> dict[str, str]:
        return {
            "angular_model": self.name,
            "angular_moment_source": "ordinary_integral_xs_closure",
            "angular_higher_moment_closure": "power",
            "angular_closure_assumption": "m1_from_total_and_momentum_xs_power_closure",
        }

    def moments(
        self,
        energy_eV: np.ndarray | float,
        lmax: int,
        *,
        sigma_total: np.ndarray | float | None = None,
        sigma_momentum: np.ndarray | float | None = None,
    ) -> np.ndarray:
        energy, m1 = _m1_from_cross_sections(energy_eV, sigma_total, sigma_momentum)
        out = _base_moments(energy, lmax)
        if lmax >= 1:
            out[1, :] = m1
        for ell in range(2, lmax + 1):
            out[ell, :] = np.power(m1, ell)
        return _finalize_moments(out)


@dataclass(frozen=True, slots=True)
class MaxEntP1AngularModel:
    quadrature_order: int = 64
    name: str = "maxent_p1"

    def metadata(self) -> dict[str, str]:
        return {
            "angular_model": self.name,
            "angular_moment_source": "ordinary_integral_xs_closure",
            "angular_higher_moment_closure": "maxent",
            "angular_closure_assumption": "maximum_entropy_p1_closure",
        }

    def moments(
        self,
        energy_eV: np.ndarray | float,
        lmax: int,
        *,
        sigma_total: np.ndarray | float | None = None,
        sigma_momentum: np.ndarray | float | None = None,
    ) -> np.ndarray:
        energy, m1 = _m1_from_cross_sections(energy_eV, sigma_total, sigma_momentum)
        out = _base_moments(energy, lmax)
        if lmax == 0:
            return out

        mu, weights = np.polynomial.legendre.leggauss(self.quadrature_order)
        legendre = _legendre_values(lmax, mu)
        targets = np.clip(m1, -0.999999, 0.999999)
        for index, target in enumerate(targets):
            if abs(target) < 1.0e-10:
                continue
            kappa = _inverse_langevin(float(target))
            exponent = kappa * mu
            exponent -= float(np.max(exponent))
            weighted = weights * np.exp(exponent)
            norm = float(np.sum(weighted))
            if norm <= 0.0 or not np.isfinite(norm):
                raise ValueError("maxent_p1 quadrature produced invalid norm")
            out[:, index] = (legendre @ weighted) / norm
        return _finalize_moments(out)


def build_angular_model(
    config: object,
) -> IsotropicAngularModel | MomentumPowerAngularModel | MaxEntP1AngularModel:
    angular = getattr(getattr(config, "physics", None), "angular_scattering", config)
    model = str(getattr(angular, "model", "isotropic")).lower()
    closure = str(getattr(angular, "higher_moment_closure", "")).lower()

    if model == "isotropic":
        if closure not in {"", "zero"}:
            raise ValueError("isotropic angular scattering requires zero closure")
        return IsotropicAngularModel()
    if model == "momentum_power":
        if closure not in {"", "power"}:
            raise ValueError("momentum_power angular scattering requires power closure")
        return MomentumPowerAngularModel()
    if model == "maxent_p1":
        if closure not in {"", "maxent"}:
            raise ValueError("maxent_p1 angular scattering requires maxent closure")
        return MaxEntP1AngularModel()
    raise ValueError(
        "physics.angular_scattering.model must be isotropic, momentum_power, or maxent_p1"
    )
