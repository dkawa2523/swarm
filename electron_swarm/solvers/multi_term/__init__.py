"""Axisymmetric multi-term Boltzmann entry point integrated with electron_swarm."""

from __future__ import annotations

from time import perf_counter

import numpy as np

from electron_swarm.collisions.ee_fp_energy import apply_fp_energy_operator
from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.constants import BOLTZMANN_J_K, TOWNSEND
from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.solvers.base import SwarmSolver
from electron_swarm.solvers.kinetic import (
    compute_rates_from_eedf,
    eepf_from_eedf,
    mean_energy_from_eedf,
    weighted_integral,
)

from .angular import LegendreBasis
from .direct import solve_direct_lmax1, solve_pn_dcs
from .grid import EnergyGrid, electron_speed_m_s, make_energy_grid
from .models import MultiTermCase, MultiTermSolution, RateSet

MULTITERM_SOLVER_NAME = "multi_term"
def _gas_number_density(config: SwarmConfig) -> float:
    cond = config.conditions
    if cond.gas_number_density_m3 is not None:
        return float(cond.gas_number_density_m3)
    assert cond.pressure_Pa is not None
    return float(cond.pressure_Pa / (BOLTZMANN_J_K * cond.gas_temperature_K))


def build_multiterm_case(
    config: SwarmConfig, cross_sections: CrossSectionSet, e_over_n_Td: float
) -> MultiTermCase:
    """Convert unified config/cross sections into the multi-term internal case."""

    number_density = _gas_number_density(config)
    return MultiTermCase(
        config=config,
        cross_sections=cross_sections,
        grid=make_energy_grid(config, cross_sections),
        e_over_n_Td=float(e_over_n_Td),
        gas_number_density_m3=number_density,
        electric_field_V_m=float(e_over_n_Td) * TOWNSEND * number_density,
    )


class MultiTermSolver(SwarmSolver):
    """Product multi-term entry point for ordinary-XS PN closure runs."""

    name = MULTITERM_SOLVER_NAME

    def _require_supported_product_method(self) -> None:
        method = self.config.solvers.multi_term.method
        if method == "pn_closure_direct":
            return
        if method == "pn_dcs":
            if self.config.physics.angular_scattering.model != "moment_table":
                raise NotImplementedError(
                    "multi_term method 'pn_dcs' requires "
                    "physics.angular_scattering.model=moment_table"
                )
            return
        raise NotImplementedError(f"Unsupported multi_term product method {method!r}")

    def solve_case(self, e_over_n_Td: float, case_id: str) -> SwarmCaseResult:
        started = perf_counter()
        self._require_supported_product_method()
        case = build_multiterm_case(self.config, self.cross_sections, e_over_n_Td)
        cfg = self.config.internal.multi_term
        if cfg.method != "moment_closure":
            raise NotImplementedError(
                f"Unsupported multi_term internal method {cfg.method!r}"
            )
        if self.config.solvers.multi_term.method == "pn_closure_direct":
            solution = solve_direct_lmax1(case, case_id, self.name)
        elif self.config.solvers.multi_term.method == "pn_dcs":
            solution = solve_pn_dcs(case, case_id, self.name)
        else:
            raise NotImplementedError(
                f"Unsupported multi_term product method "
                f"{self.config.solvers.multi_term.method!r}"
            )
        run_time_s = perf_counter() - started
        return self._to_case_result(case, solution, case_id, run_time_s)

    def _to_case_result(
        self,
        case: MultiTermCase,
        solution: MultiTermSolution,
        case_id: str,
        run_time_s: float,
    ) -> SwarmCaseResult:
        transport = solution.transport
        flux = transport.flux
        number_density = case.gas_number_density_m3
        diffusion_l = (
            flux.diffusion_longitudinal_m2_s
            if flux.diffusion_longitudinal_m2_s is not None
            else float("nan")
        )
        diffusion_t = (
            flux.diffusion_transverse_m2_s
            if flux.diffusion_transverse_m2_s is not None
            else diffusion_l
        )
        eedf = np.asarray(solution.eedf_eV_inv, dtype=float)
        rates = solution.rates
        mean_energy = float(solution.diagnostics.mean_energy_eV)
        ee_metadata: dict[str, object] = {}
        ee = case.config.physics.electron_electron
        if bool(ee.enabled) and ee.model == "fp_energy":
            if ee.strength_model == "density_based":
                raise NotImplementedError(
                    "electron_electron fp_energy strength_model='density_based' "
                    "is not implemented"
                )
            eedf, ee_operator_metadata = apply_fp_energy_operator(
                solution.energy_eV,
                solution.widths_eV,
                eedf,
                relaxation_fraction=ee.relaxation_fraction,
                conserve_mean_energy=ee.conserve_mean_energy,
                fallback_temperature_eV=ee.fallback_temperature_eV,
            )
            rate_data = compute_rates_from_eedf(
                case.config,
                case.cross_sections,
                solution.energy_eV,
                solution.widths_eV,
                eedf,
                case_id=case_id,
                e_over_n_Td=case.e_over_n_Td,
                solver_name=self.name,
            )
            rates = RateSet(
                tuple(rate_data.rates),
                number_density * rate_data.ionization_rate_m3_s,
                number_density * rate_data.attachment_rate_m3_s,
            )
            mean_energy = mean_energy_from_eedf(
                solution.energy_eV,
                solution.widths_eV,
                eedf,
            )
            ee_metadata = {
                "electron_electron_treatment": "fp_energy",
                "electron_electron_affects_eedf": True,
                "electron_electron_affects_rates": True,
                "electron_electron_affects_transport": False,
                "electron_electron_transport_stale": True,
                **ee_operator_metadata,
            }
        net_freq = rates.effective_growth_frequency_s_inv
        effective_townsend = net_freq / max(
            abs(flux.drift_velocity_m_s) * number_density, 1.0e-300
        )
        eepf = eepf_from_eedf(solution.energy_eV, eedf)
        energy_grid = case.config.internal.multi_term.energy_grid
        metadata = {
            "backend": case.config.solvers.multi_term.method,
            "run_label": case_id,
            "sweep_param": "E_over_N_Td",
            "sweep_value": float(case.e_over_n_Td),
            "run_time_s": float(run_time_s),
            "gas_number_density_m-3": float(number_density),
            "electric_field_V_m": float(case.electric_field_V_m),
            "lmax": int(case.config.internal.multi_term.lmax),
            "direct_pn_operator": False,
            "transport_definition": "flux",
            "grid_n_cells": int(len(solution.energy_eV)),
            "grid_min_eV": float(np.min(solution.energy_eV)),
            "grid_max_eV": float(np.max(solution.energy_eV)),
            "grid_spacing": energy_grid.spacing,
            "threshold_refined": bool(
                energy_grid.refine.enabled and len(solution.energy_eV) > energy_grid.n
            ),
            "cross_section_high_energy_extrapolation": (
                case.config.cross_sections.high_energy_extrapolation
            ),
            "tail_fraction": solution.diagnostics.eedf_tail_fraction,
            "power_balance_relative_residual": (
                solution.diagnostics.power_balance_relative_residual
            ),
            "warnings": "; ".join(solution.diagnostics.warnings),
            "normalization_integral": weighted_integral(eedf, solution.widths_eV),
            "convolution_ionization_rate_coefficient_m3_s": float(
                sum(
                    r.mixture_weighted_rate_m3_s
                    for r in rates.rates
                    if r.process_type == "ionization"
                )
            ),
            "convolution_attachment_rate_coefficient_m3_s": float(
                sum(
                    r.mixture_weighted_rate_m3_s
                    for r in rates.rates
                    if r.process_type == "attachment"
                )
            ),
            "convolution_effective_rate_coefficient_m3_s": float(
                net_freq / number_density if number_density > 0.0 else np.nan
            ),
            "convolution_net_ionization_frequency_s-1": float(net_freq),
            "net_ionization_frequency_model": "convolution_effective_rate",
            "effective_townsend_1_m": float(
                net_freq / max(abs(flux.drift_velocity_m_s), 1.0e-300)
            ),
        }
        metadata.update(solution.metadata)
        metadata.update(ee_metadata)
        return SwarmCaseResult(
            solver=self.name,
            case_id=case_id,
            e_over_n_Td=case.e_over_n_Td,
            mean_energy_eV=mean_energy,
            drift_velocity_m_s=flux.drift_velocity_m_s,
            mobility_m2_V_s=flux.mobility_m2_V_s,
            reduced_mobility_m2_V_s_m3=flux.mobility_m2_V_s * number_density,
            diffusion_L_m2_s=diffusion_l,
            diffusion_T_m2_s=diffusion_t,
            reduced_diffusion_L_m2_s_m3=diffusion_l * number_density,
            reduced_diffusion_T_m2_s_m3=diffusion_t * number_density,
            net_ionization_frequency_s=net_freq,
            effective_townsend_m2=effective_townsend,
            energy_eV=solution.energy_eV,
            eedf=eedf,
            eepf=eepf,
            energy_widths_eV=solution.widths_eV,
            rates=list(rates.rates),
            metadata=metadata,
            transport=transport,
        )


__all__ = [
    "EnergyGrid",
    "LegendreBasis",
    "MultiTermSolver",
    "MultiTermCase",
    "MultiTermSolution",
    "build_multiterm_case",
    "electron_speed_m_s",
]
