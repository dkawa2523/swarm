"""Moment-closure backend for the multi-term Boltzmann entry point."""

from __future__ import annotations

import numpy as np
from scipy import optimize

from electron_swarm.collisions.ee_fp_energy import apply_fp_energy_operator
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
    MOMENTUM_TYPES,
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


def _angular_moments_for_grid(
    case: MultiTermCase,
    collisions: ProjectedCollisionData,
) -> np.ndarray:
    cfg = case.config.internal.multi_term
    angular_model = build_angular_model(case.config)
    energy = case.grid.centers_eV

    if angular_model.name in {"isotropic", "moment_table"}:
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

    return moments


def _pn_dcs_momentum_frequency(
    collisions: ProjectedCollisionData,
    angular_moments: np.ndarray,
) -> np.ndarray:
    m1 = np.clip(angular_moments[1], -1.0, 1.0)
    transport_factor = np.clip(1.0 - m1, 1.0e-6, 2.0)
    elastic_like = np.zeros_like(collisions.momentum_frequency_s_inv)
    for item in collisions.processes:
        if item.contributes_to_momentum and item.process.process_type in MOMENTUM_TYPES:
            elastic_like += item.frequency_s_inv
    inelastic_surrogate = np.maximum(
        collisions.momentum_frequency_s_inv - elastic_like,
        0.0,
    )
    return elastic_like * transport_factor + inelastic_surrogate


def _total_momentum_frequency(
    case: MultiTermCase,
    collisions: ProjectedCollisionData,
    F: np.ndarray,
    momentum_frequency_s_inv: np.ndarray | None,
) -> float:
    if momentum_frequency_s_inv is None:
        return total_momentum_frequency(case, collisions, F)
    with np.errstate(under="ignore"):
        return float(
            np.sum(
                np.maximum(momentum_frequency_s_inv, 1.0)
                * F
                * case.grid.widths_eV
            )
        )


def mean_energy_from_power_balance(
    case: MultiTermCase,
    collisions: ProjectedCollisionData,
    shape: str,
    case_id: str,
    solver_name: str,
    momentum_frequency_s_inv: np.ndarray | None = None,
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
        nu_m = _total_momentum_frequency(case, collisions, F, momentum_frequency_s_inv)
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
    nu_m = _total_momentum_frequency(case, collisions, F, momentum_frequency_s_inv)
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
        angular_moments = _angular_moments_for_grid(case, collisions)
        momentum_frequency = None
        if case.config.solvers.multi_term.method == "pn_dcs":
            momentum_frequency = _pn_dcs_momentum_frequency(
                collisions,
                angular_moments,
            )
        mean_e, F0, rates, drift, mobility, residual = mean_energy_from_power_balance(
            case,
            collisions,
            cfg.eedf_shape,
            case_id,
            self.solver_name,
            momentum_frequency,
        )
        transport_mean_e = mean_e
        solution_metadata: dict[str, object] = {}
        ee = case.config.physics.electron_electron
        if ee.enabled and ee.model == "fp_energy":
            if ee.strength_model == "density_based":
                raise NotImplementedError(
                    "electron_electron fp_energy strength_model='density_based' is not implemented"
                )
            F0, ee_metadata = apply_fp_energy_operator(
                case.grid.centers_eV,
                case.grid.widths_eV,
                F0,
                relaxation_fraction=ee.relaxation_fraction,
                conserve_mean_energy=ee.conserve_mean_energy,
                fallback_temperature_eV=ee.fallback_temperature_eV,
            )
            rates = compute_rate_set(
                case,
                collisions,
                F0,
                case_id,
                self.solver_name,
            )
            mean_e = float(np.sum(case.grid.centers_eV * F0 * case.grid.widths_eV))
            solution_metadata.update(ee_metadata)
            solution_metadata.update(
                {
                    "electron_electron_treatment": "fp_energy",
                    "electron_electron_affects_eedf": True,
                    "electron_electron_affects_rates": True,
                    "electron_electron_affects_transport": False,
                    "electron_electron_transport_stale": True,
                }
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

        diffusion = max((2.0 / 3.0) * transport_mean_e * abs(mobility), 0.0)
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
            metadata=solution_metadata,
            angular_moments=angular_moments,
        )
