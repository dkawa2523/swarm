"""Angular-scattering moment models for product schema v2."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Literal, Protocol

import numpy as np

AngularModelName = Literal["isotropic", "momentum_power", "maxent_p1", "moment_table"]
HigherMomentClosure = Literal["zero", "power", "maxent", "table"]
ANGULAR_METADATA_KEYS = (
    "angular_model",
    "angular_moment_source",
    "exact_dcs_based",
    "ordinary_integral_xs_closure",
)


class AngularMomentProvider(Protocol):
    name: str

    def metadata(self) -> dict[str, str]:
        ...

    def moments(
        self,
        energy_eV: np.ndarray | float,
        lmax: int,
        *,
        sigma_total: np.ndarray | float | None = None,
        sigma_momentum: np.ndarray | float | None = None,
    ) -> np.ndarray:
        ...

    def sample_mu(
        self,
        energy_eV: np.ndarray | float,
        rng: np.random.Generator,
        size: int | tuple[int, ...] | None = None,
        *,
        sigma_total: np.ndarray | float | None = None,
        sigma_momentum: np.ndarray | float | None = None,
    ) -> np.ndarray | float:
        ...


def _energy_vector(energy_eV: np.ndarray | float) -> np.ndarray:
    energy = np.asarray(energy_eV, dtype=float).reshape(-1)
    if not np.all(np.isfinite(energy)):
        raise ValueError("energy_eV must be finite")
    return energy


def _sample_shape(
    energy_eV: np.ndarray | float,
    energy: np.ndarray,
    size: int | tuple[int, ...] | None,
) -> tuple[int, ...]:
    extra = () if size is None else (size,) if isinstance(size, int) else tuple(size)
    if np.asarray(energy_eV).ndim == 0:
        return extra
    return (energy.size, *extra)


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


def _sample_maxent_mu(
    kappa: float,
    rng: np.random.Generator,
    size: tuple[int, ...],
) -> np.ndarray | float:
    if abs(kappa) < 1.0e-10:
        return rng.uniform(-1.0, 1.0, size=size)
    if kappa < 0.0:
        return -_sample_maxent_mu(-kappa, rng, size)
    u = rng.random(size=size)
    exp_neg2k = np.exp(-2.0 * min(kappa, 350.0))
    values = 1.0 + np.log(u + (1.0 - u) * exp_neg2k) / kappa
    return np.clip(values, -1.0, 1.0)


@dataclass(frozen=True, slots=True)
class IsotropicAngularModel:
    name: str = "isotropic"

    def metadata(self) -> dict[str, str]:
        return {
            "angular_model": self.name,
            "angular_moment_source": "isotropic_closure",
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

    def sample_mu(
        self,
        energy_eV: np.ndarray | float,
        rng: np.random.Generator,
        size: int | tuple[int, ...] | None = None,
        *,
        sigma_total: np.ndarray | float | None = None,
        sigma_momentum: np.ndarray | float | None = None,
    ) -> np.ndarray | float:
        energy = _energy_vector(energy_eV)
        return rng.uniform(-1.0, 1.0, size=_sample_shape(energy_eV, energy, size))


@dataclass(frozen=True, slots=True)
class MomentumPowerAngularModel:
    name: str = "momentum_power"

    def metadata(self) -> dict[str, str]:
        return {
            "angular_model": self.name,
            "angular_moment_source": "ordinary_integral_xs_closure",
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

    def sample_mu(
        self,
        energy_eV: np.ndarray | float,
        rng: np.random.Generator,
        size: int | tuple[int, ...] | None = None,
        *,
        sigma_total: np.ndarray | float | None = None,
        sigma_momentum: np.ndarray | float | None = None,
    ) -> np.ndarray | float:
        raise NotImplementedError(
            "momentum_power angular closure does not define a unique MC sampler"
        )


@dataclass(frozen=True, slots=True)
class MaxEntP1AngularModel:
    quadrature_order: int = 64
    name: str = "maxent_p1"

    def metadata(self) -> dict[str, str]:
        return {
            "angular_model": self.name,
            "angular_moment_source": "ordinary_integral_xs_closure",
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

    def sample_mu(
        self,
        energy_eV: np.ndarray | float,
        rng: np.random.Generator,
        size: int | tuple[int, ...] | None = None,
        *,
        sigma_total: np.ndarray | float | None = None,
        sigma_momentum: np.ndarray | float | None = None,
    ) -> np.ndarray | float:
        energy, m1 = _m1_from_cross_sections(energy_eV, sigma_total, sigma_momentum)
        output_shape = _sample_shape(energy_eV, energy, size)
        if np.asarray(energy_eV).ndim == 0:
            kappa = _inverse_langevin(float(m1[0]))
            return _sample_maxent_mu(kappa, rng, output_shape)
        extra_shape = output_shape[1:]
        out = np.empty(output_shape, dtype=float)
        for index, target in enumerate(m1):
            out[index] = _sample_maxent_mu(
                _inverse_langevin(float(target)),
                rng,
                extra_shape,
            )
        return out


def _moment_column_index(name: str) -> int | None:
    if len(name) < 2 or not name.startswith("m"):
        return None
    suffix = name[1:]
    if not suffix.isdigit():
        return None
    return int(suffix)


def _load_moment_table(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open("r", encoding="utf-8", newline="") as fp:
        reader = csv.DictReader(fp)
        if reader.fieldnames is None:
            raise ValueError("moment_table CSV must have a header row")
        fieldnames = [name.strip() for name in reader.fieldnames]
        if "energy_eV" not in fieldnames:
            raise ValueError("moment_table CSV requires energy_eV column")
        moment_indices = {
            index: name
            for name in fieldnames
            if (index := _moment_column_index(name)) is not None
        }
        if 0 not in moment_indices:
            raise ValueError("moment_table CSV requires m0 column")
        expected = set(range(max(moment_indices) + 1))
        if set(moment_indices) != expected:
            raise ValueError("moment_table CSV moment columns must be contiguous from m0")
        energies: list[float] = []
        rows: list[list[float]] = []
        for row in reader:
            energies.append(float(row["energy_eV"]))
            rows.append([float(row[moment_indices[ell]]) for ell in sorted(moment_indices)])
    if not energies:
        raise ValueError("moment_table CSV must contain at least one row")
    energy = np.asarray(energies, dtype=float)
    moments = np.asarray(rows, dtype=float).T
    if energy.ndim != 1 or moments.ndim != 2 or moments.shape[1] != energy.size:
        raise ValueError("moment_table CSV shape is invalid")
    if not np.all(np.isfinite(energy)) or np.any(energy < 0.0):
        raise ValueError("moment_table energy_eV values must be finite and nonnegative")
    if np.any(np.diff(energy) <= 0.0):
        raise ValueError("moment_table energy_eV values must be strictly increasing")
    if not np.all(np.isfinite(moments)):
        raise ValueError("moment_table moments must be finite")
    if not np.allclose(moments[0], 1.0, rtol=1.0e-8, atol=1.0e-8):
        raise ValueError("moment_table m0 values must be 1")
    if np.any((moments[1:] < -1.0) | (moments[1:] > 1.0)):
        raise ValueError("moment_table moments must be in [-1, 1]")
    moments[0, :] = 1.0
    return energy, moments


@dataclass(frozen=True, slots=True)
class MomentTableAngularModel:
    path: Path
    provenance: Literal["precomputed_moments", "dcs_derived"] = "precomputed_moments"
    extrapolation: Literal["error"] = "error"
    name: str = "moment_table"

    def metadata(self) -> dict[str, str]:
        return {
            "angular_model": self.name,
            "angular_moment_source": "moment_table",
        }

    def moments(
        self,
        energy_eV: np.ndarray | float,
        lmax: int,
        *,
        sigma_total: np.ndarray | float | None = None,
        sigma_momentum: np.ndarray | float | None = None,
    ) -> np.ndarray:
        if self.extrapolation != "error":
            raise ValueError("moment_table only supports extrapolation=error")
        energy = _energy_vector(energy_eV)
        table_energy, table_moments = _load_moment_table(self.path)
        if lmax >= table_moments.shape[0]:
            raise ValueError("moment_table does not contain enough moment columns for lmax")
        if energy.size and (
            energy[0] < table_energy[0] - 1.0e-12
            or energy[-1] > table_energy[-1] + 1.0e-12
        ):
            raise ValueError("moment_table energy range does not cover solver energy grid")
        out = _base_moments(energy, lmax)
        for ell in range(1, lmax + 1):
            out[ell] = np.interp(energy, table_energy, table_moments[ell])
        return _finalize_moments(out)

    def sample_mu(
        self,
        energy_eV: np.ndarray | float,
        rng: np.random.Generator,
        size: int | tuple[int, ...] | None = None,
        *,
        sigma_total: np.ndarray | float | None = None,
        sigma_momentum: np.ndarray | float | None = None,
    ) -> np.ndarray | float:
        raise NotImplementedError(
            "moment_table moments do not define a unique MC sampler"
        )


def build_angular_model(
    config: object,
) -> (
    IsotropicAngularModel
    | MomentumPowerAngularModel
    | MaxEntP1AngularModel
    | MomentTableAngularModel
):
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
    if model == "moment_table":
        if closure not in {"", "table"}:
            raise ValueError("moment_table angular scattering requires table closure")
        table = getattr(angular, "moment_table", None)
        if table is None:
            raise ValueError("moment_table angular scattering requires moment_table config")
        return MomentTableAngularModel(
            path=Path(table.path),
            provenance=table.provenance,
            extrapolation=table.extrapolation,
        )
    raise ValueError(
        "physics.angular_scattering.model must be isotropic, momentum_power, maxent_p1, or moment_table"
    )


def expected_angular_metadata(config: object) -> dict[str, object]:
    angular = getattr(config.physics, "angular_scattering")
    metadata: dict[str, object] = dict(build_angular_model(config).metadata())
    is_moment_table = angular.model == "moment_table"
    table = getattr(angular, "moment_table", None)
    metadata.update(
        {
            "exact_dcs_based": bool(
                is_moment_table
                and table is not None
                and table.provenance == "dcs_derived"
            ),
            "ordinary_integral_xs_closure": not is_moment_table,
        }
    )
    return metadata
