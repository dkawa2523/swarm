"""Moment-closure backend for the multi-term Boltzmann entry point."""

from __future__ import annotations

import numpy as np
from scipy import optimize

from electron_swarm.core.constants import (
    BOLTZMANN_J_K,
    E_CHARGE_C,
    ELECTRON_MASS_KG,
    EV_TO_J,
)
from electron_swarm.core.transport import (
    BulkTransport,
    FluxTransport,
    TransportMetadata,
    TransportSet,
)
from electron_swarm.physics.angular_scattering import build_angular_model

from .angular import LegendreBasis
from .diagnostics import SolverDiagnostics
from .grid import EnergyGrid
from .models import MultiTermCase, MultiTermSolution, RateSet
from .projection import (
    ProjectedCollisionData,
    compute_rate_set,
    elastic_power_loss_eV_s,
    project_collision_data,
    total_momentum_frequency,
    validity_warnings,
)


def maxwellian_energy_pdf(grid: EnergyGrid, mean_energy_eV: float) -> np.ndarray:
    theta = max(2.0 * float(mean_energy_eV) / 3.0, 1.0e-5)
    e = grid.centers_eV
    exponent = np.clip(-e / theta, -700.0, 80.0)
    with np.errstate(under="ignore"):
        F = (
            2.0
            / np.sqrt(np.pi)
            * np.sqrt(np.maximum(e, 0.0))
            / theta**1.5
            * np.exp(exponent)
        )
    return grid.normalize_energy_pdf(F)


def druyvesteyn_energy_pdf(grid: EnergyGrid, mean_energy_eV: float) -> np.ndarray:
    e = grid.centers_eV
    alpha = max(float(mean_energy_eV), 1.0e-5)
    exponent = np.clip(-0.55 * (e / alpha) ** 2, -700.0, 80.0)
    with np.errstate(under="ignore"):
        F = np.sqrt(np.maximum(e, 0.0)) * np.exp(exponent)
    return grid.normalize_energy_pdf(F)


def _angular_cross_sections(
    case: MultiTermCase,
    collisions: ProjectedCollisionData,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sigma_total = np.zeros_like(case.grid.centers_eV)
    sigma_momentum = np.zeros_like(case.grid.centers_eV)
    for item in collisions.processes:
        weighted_sigma = item.fraction * item.sigma_m2
        sigma_total += weighted_sigma
        if item.contributes_to_momentum:
            sigma_momentum += weighted_sigma
    valid = np.isfinite(sigma_total) & np.isfinite(sigma_momentum) & (sigma_total > 0.0)
    return sigma_total, sigma_momentum, valid


def _angular_moment_closure(
    case: MultiTermCase,
    collisions: ProjectedCollisionData,
    energy_pdf_eV_inv: np.ndarray,
) -> tuple[np.ndarray, dict[str, object]]:
    cfg = case.config.internal.multi_term
    angular_model = build_angular_model(case.config)
    energy = case.grid.centers_eV

    if angular_model.name == "isotropic":
        moments = angular_model.moments(energy, cfg.lmax)
        valid = np.ones_like(energy, dtype=bool)
    else:
        sigma_total, sigma_momentum, valid = _angular_cross_sections(case, collisions)
        safe_total = np.where(valid, sigma_total, 1.0)
        safe_momentum = np.where(valid, sigma_momentum, 1.0)
        moments = angular_model.moments(
            energy,
            cfg.lmax,
            sigma_total=safe_total,
            sigma_momentum=safe_momentum,
        )
        if not np.all(valid):
            moments[:, ~valid] = 0.0
            moments[0, ~valid] = 1.0

    weights = case.grid.normalize_energy_pdf(energy_pdf_eV_inv) * case.grid.widths_eV
    averaged = np.sum(moments * weights[np.newaxis, :], axis=1)
    metadata = {
        "angular_moment_lmax": int(cfg.lmax),
        "angular_moment_valid_fraction": float(np.mean(valid)),
        "angular_moment_invalid_cells": int(np.size(valid) - np.count_nonzero(valid)),
        "angular_moment_m1_average": float(averaged[1]) if cfg.lmax >= 1 else 0.0,
        "angular_moment_highest_average": float(averaged[-1]),
    }
    return moments, metadata


def mean_energy_from_power_balance(
    case: MultiTermCase,
    collisions: ProjectedCollisionData,
    shape: str,
    case_id: str,
    solver_name: str,
) -> tuple[float, np.ndarray, RateSet, float, float, float]:
    grid = case.grid
    electric_field = case.electric_field_V_m
    thermal_mean = (
        1.5 * BOLTZMANN_J_K * case.config.conditions.gas_temperature_K / EV_TO_J
    )

    def eedf_for(mean_e: float) -> np.ndarray:
        if shape == "druyvesteyn":
            return druyvesteyn_energy_pdf(grid, mean_e)
        return maxwellian_energy_pdf(grid, mean_e)

    def power_terms(mean_e: float):
        F = eedf_for(mean_e)
        nu_m = total_momentum_frequency(case, collisions, F)
        mobility = E_CHARGE_C / (ELECTRON_MASS_KG * max(nu_m, 1.0))
        drift = mobility * electric_field
        p_in = abs(electric_field * drift)
        rates = compute_rate_set(case, collisions, F, case_id, solver_name)
        p_loss = sum(
            r.power_loss_eV_s or 0.0
            for r in rates.rates
            if r.process_type in {"excitation", "ionization", "superelastic"}
        )
        p_loss += elastic_power_loss_eV_s(case, collisions, F, thermal_mean)
        return p_in, p_loss, rates, drift, mobility

    def balance(mean_e: float) -> float:
        p_in, p_loss, *_ = power_terms(mean_e)
        return p_in - p_loss

    lo = max(thermal_mean * 0.5, 1.0e-3)
    hi = min(max(grid.edges_eV[-1] * 0.8, 5.0), grid.edges_eV[-1] * 0.98)
    points = np.geomspace(lo, hi, 48) if lo > 0.0 else np.linspace(lo, hi, 48)
    values = np.array([balance(float(x)) for x in points], dtype=float)
    finite = np.isfinite(values)
    mean: float | None = None
    if np.any(finite):
        good_points = points[finite]
        good_values = values[finite]
        for left, right, b_left, b_right in zip(
            good_points[:-1],
            good_points[1:],
            good_values[:-1],
            good_values[1:],
            strict=False,
        ):
            if b_left == 0.0:
                mean = float(left)
                break
            if b_left * b_right < 0.0:
                mean = float(optimize.brentq(balance, left, right, maxiter=100))
                break
    if mean is None:
        res = optimize.minimize_scalar(
            lambda x: abs(balance(x)),
            bounds=(lo, hi),
            method="bounded",
            options={"xatol": 1.0e-4},
        )
        if not res.success:
            raise RuntimeError("moment-closure power balance did not converge")
        mean = float(res.x)
        if mean <= lo * 1.01 or mean >= hi * 0.99:
            raise RuntimeError(
                "moment-closure power balance minimized at the energy-grid boundary"
            )
    F = eedf_for(mean)
    p_in, p_loss, rates, drift, mobility = power_terms(mean)
    residual = abs(p_in - p_loss) / max(abs(p_in), abs(p_loss), 1.0e-300)
    if not np.isfinite(residual) or residual > 5.0e-2:
        raise RuntimeError(
            "moment-closure power balance residual is too large "
            f"({residual:.3e})"
        )
    nu_m = total_momentum_frequency(case, collisions, F)
    mobility = E_CHARGE_C / (ELECTRON_MASS_KG * max(nu_m, 1.0))
    drift = mobility * electric_field
    return mean, F, rates, drift, mobility, float(residual)


class MomentClosureEngine:
    def __init__(self, solver_name: str) -> None:
        self.solver_name = solver_name

    def solve(
        self,
        case: MultiTermCase,
        case_id: str,
        previous: MultiTermSolution | None = None,
    ) -> MultiTermSolution:
        cfg = case.config.internal.multi_term
        basis = LegendreBasis(cfg.lmax)
        collisions = project_collision_data(case)
        mean_e, F0, rates, drift, mobility, residual = mean_energy_from_power_balance(
            case, collisions, cfg.eedf_shape, case_id, self.solver_name
        )
        coeff = np.zeros((basis.n_terms, case.grid.n_cells), dtype=float)
        coeff[0] = F0
        dFde = np.gradient(F0, case.grid.centers_eV, edge_order=1)
        anis = -dFde
        norm = np.max(np.abs(anis)) or 1.0
        anis = anis / norm * np.max(F0)
        if basis.lmax >= 1:
            moment = np.sum(case.grid.speeds_m_s * anis * case.grid.widths_eV / 3.0)
            scale = drift / moment if abs(moment) > 1.0e-300 else 0.0
            coeff[1] = anis * scale
        for ell in range(2, basis.n_terms):
            coeff[ell] = coeff[ell - 1] * min(0.35, 0.8 / (ell + 1))
        angular_moments, angular_metadata = _angular_moment_closure(
            case, collisions, F0
        )

        diffusion = max((2.0 / 3.0) * mean_e * abs(mobility), 0.0)
        flux = FluxTransport.from_drift_and_field(
            drift,
            case.electric_field_V_m,
            diffusion_longitudinal_m2_s=diffusion,
            diffusion_transverse_m2_s=diffusion,
        )
        estimated_bulk = None
        if cfg.hydrodynamic:
            nu_eff = rates.effective_growth_frequency_s_inv
            corr_W = diffusion * nu_eff / max(
                abs(drift),
                np.sqrt(abs(diffusion * max(abs(nu_eff), 1.0))),
                1.0,
            )
            estimated_bulk = BulkTransport.from_drift_and_field(
                drift + corr_W,
                case.electric_field_V_m,
                diffusion_longitudinal_m2_s=max(
                    diffusion
                    * (1.0 + 0.1 * np.tanh(nu_eff / max(abs(drift), 1.0))),
                    0.0,
                ),
                diffusion_transverse_m2_s=diffusion,
            )
        notes = (
            ("moment_closure_bulk_estimate_available",)
            if estimated_bulk is not None
            else ()
        )
        transport = TransportSet.from_flux_only(
            flux,
            rates.ionization_frequency_s_inv,
            rates.attachment_frequency_s_inv,
            TransportMetadata(
                solver=self.solver_name,
                coefficient_definition="flux",
                swarm_condition="local_flux",
                notes=notes,
            ),
        )
        tail_mask = case.grid.centers_eV >= 0.9 * case.grid.edges_eV[-1]
        with np.errstate(under="ignore"):
            tail = float(np.sum(F0[tail_mask] * case.grid.widths_eV[tail_mask]))
        diagnostics = SolverDiagnostics(
            warnings=validity_warnings(case, collisions, F0),
            mean_energy_eV=float(mean_e),
            eedf_tail_fraction=tail,
            power_balance_relative_residual=residual,
        )
        return MultiTermSolution(
            case.grid.centers_eV,
            case.grid.widths_eV,
            coeff,
            F0,
            rates,
            transport,
            estimated_bulk,
            diagnostics,
            method_used="moment_closure",
            metadata=angular_metadata,
            angular_moments=angular_moments,
        )
