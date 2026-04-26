"""Axisymmetric multi-term Boltzmann entry point integrated with electron_swarm."""

from __future__ import annotations

from time import perf_counter

import numpy as np

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.constants import BOLTZMANN_J_K, TOWNSEND
from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.solvers.base import SwarmSolver

from .angular import LegendreBasis
from .closure import (
    MomentClosureEngine,
    druyvesteyn_energy_pdf,
    maxwellian_energy_pdf,
)
from .grid import EnergyGrid, electron_speed_m_s, make_energy_grid
from .models import MultiTermCase, MultiTermSolution
from .operator import (
    DensityNormalizationConstraint,
    LegendreBlockLayout,
    MultiTermOperatorBackend,
    OperatorBackendUnavailable,
    OperatorAssemblyDiagnostics,
    OperatorSolveState,
    OperatorSystem,
    assemble_operator_system,
    build_density_normalization_constraint,
    build_operator_assembly_diagnostics,
    solve_operator_system,
)
from .projection import (
    ProjectedCollisionData,
    compute_rate_set,
    elastic_power_loss_eV_s,
    energy_loss_eV,
    project_collision_data,
    threshold_masked_sigma,
)

MULTITERM_SOLVER_NAME = "multiterm_boltzmann"


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
        grid=make_energy_grid(config),
        e_over_n_Td=float(e_over_n_Td),
        gas_number_density_m3=number_density,
        electric_field_V_m=float(e_over_n_Td) * TOWNSEND * number_density,
    )


def _project_collision_data(case: MultiTermCase) -> ProjectedCollisionData:
    return project_collision_data(case)


def _compute_rate_set(
    case: MultiTermCase,
    collisions: ProjectedCollisionData,
    energy_pdf_eV_inv: np.ndarray,
    case_id: str,
):
    return compute_rate_set(
        case,
        collisions,
        energy_pdf_eV_inv,
        case_id,
        MULTITERM_SOLVER_NAME,
    )


def _elastic_power_loss_eV_s(
    case: MultiTermCase,
    collisions: ProjectedCollisionData,
    F: np.ndarray,
    thermal_mean_eV: float,
) -> float:
    return elastic_power_loss_eV_s(case, collisions, F, thermal_mean_eV)


_threshold_masked_sigma = threshold_masked_sigma
_energy_loss_eV = energy_loss_eV


class MultiTermBoltzmannSolver(SwarmSolver):
    """Unified axisymmetric multi-term Boltzmann entry point.

    The solver exposes a conservative moment-closure estimator plus a scoped
    production lmax>1 sparse operator path. Operator hydrodynamic transport is
    emitted only when finite-k extraction succeeds.
    """

    name = MULTITERM_SOLVER_NAME

    def solve_case(self, e_over_n_Td: float, case_id: str) -> SwarmCaseResult:
        started = perf_counter()
        case = build_multiterm_case(self.config, self.cross_sections, e_over_n_Td)
        cfg = self.config.multiterm_boltzmann
        if cfg.method == "operator":
            solution = MultiTermOperatorBackend(self.name).solve(case, case_id)
            run_time_s = perf_counter() - started
            return self._to_case_result(
                case, solution, case_id, run_time_s, cfg.method
            )
        if cfg.method == "hybrid":
            raise NotImplementedError(
                "multiterm_boltzmann.method='hybrid' is reserved until the "
                "operator and moment-closure dispatch policy is implemented."
            )
        if cfg.method != "moment_closure":
            raise OperatorBackendUnavailable(
                f"Unsupported multiterm_boltzmann.method={cfg.method!r}"
            )
        solution = MomentClosureEngine(self.name).solve(case, case_id)
        run_time_s = perf_counter() - started
        return self._to_case_result(case, solution, case_id, run_time_s, cfg.method)

    def _to_case_result(
        self,
        case: MultiTermCase,
        solution: MultiTermSolution,
        case_id: str,
        run_time_s: float,
        requested_method: str,
    ) -> SwarmCaseResult:
        transport = solution.transport
        flux = transport.flux
        bulk = transport.bulk
        source = transport.source
        estimated_bulk = solution.estimated_bulk
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
        net_freq = solution.rates.effective_growth_frequency_s_inv
        effective_townsend = net_freq / max(
            abs(flux.drift_velocity_m_s) * number_density, 1.0e-300
        )
        eedf = solution.eedf_eV_inv
        eepf = eedf / np.sqrt(np.maximum(solution.energy_eV, 1.0e-30))
        estimated_gradient = (
            estimated_bulk.drift_velocity_m_s - flux.drift_velocity_m_s
            if estimated_bulk is not None
            else np.nan
        )
        is_operator_lmax1 = solution.method_used == "operator_lmax1_two_term"
        is_operator_flux = solution.method_used == "operator_flux"
        is_operator_hydrodynamic = solution.method_used == "operator_hydrodynamic"
        is_operator = is_operator_flux or is_operator_hydrodynamic
        if is_operator_lmax1:
            backend = "native_operator_lmax1_two_term"
            physical_validity = "operator_lmax1_two_term_reference"
        elif is_operator_hydrodynamic:
            backend = "native_operator_sparse"
            physical_validity = "operator_hydrodynamic_b0_dc_m0_integral_cross_sections"
        elif is_operator_flux:
            backend = "native_operator_sparse"
            physical_validity = "operator_flux_b0_dc_m0_integral_cross_sections"
        else:
            backend = "native_moment_closure"
            physical_validity = "moment_closure_estimate_not_multiterm_operator"
        metadata = {
            "backend": backend,
            "run_label": case_id,
            "sweep_param": "E_over_N_Td",
            "sweep_value": float(case.e_over_n_Td),
            "run_time_s": float(run_time_s),
            "gas_number_density_m-3": float(number_density),
            "electric_field_V_m": float(case.electric_field_V_m),
            "lmax": int(case.config.multiterm_boltzmann.lmax),
            "hydrodynamic": bool(bulk is not None),
            "hydrodynamic_requested": bool(case.config.multiterm_boltzmann.hydrodynamic),
            "transport_definition": "flux_bulk_source" if bulk is not None else "flux",
            "physical_validity": physical_validity,
            "cross_section_high_energy_extrapolation": (
                case.config.cross_sections.high_energy_extrapolation
            ),
            "multiterm_method_requested": requested_method,
            "multiterm_method_used": solution.method_used,
            "eedf_shape": (
                "two_term_reference"
                if is_operator_lmax1
                else "operator_solution"
                if is_operator
                else case.config.multiterm_boltzmann.eedf_shape
            ),
            "tail_fraction": solution.diagnostics.eedf_tail_fraction,
            "power_balance_relative_residual": (
                solution.diagnostics.power_balance_relative_residual
            ),
            "warnings": "; ".join(solution.diagnostics.warnings),
            "normalization_integral": float(np.sum(eedf * solution.widths_eV)),
            "convolution_ionization_rate_coefficient_m3_s": float(
                sum(
                    r.mixture_weighted_rate_m3_s
                    for r in solution.rates.rates
                    if r.process_type == "ionization"
                )
            ),
            "convolution_attachment_rate_coefficient_m3_s": float(
                sum(
                    r.mixture_weighted_rate_m3_s
                    for r in solution.rates.rates
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
            "source_gradient_velocity_m_s": source.gradient_velocity_m_s
            if source.gradient_velocity_m_s is not None
            else np.nan,
            "estimated_bulk_definition": (
                "moment_closure_estimate" if estimated_bulk is not None else ""
            ),
            "estimated_bulk_drift_velocity_m_s": (
                estimated_bulk.drift_velocity_m_s
                if estimated_bulk is not None
                else np.nan
            ),
            "estimated_bulk_reduced_diffusion_L_m-1_s-1": (
                estimated_bulk.diffusion_longitudinal_m2_s * number_density
                if estimated_bulk is not None
                and estimated_bulk.diffusion_longitudinal_m2_s is not None
                else np.nan
            ),
            "estimated_bulk_reduced_diffusion_T_m-1_s-1": (
                estimated_bulk.diffusion_transverse_m2_s * number_density
                if estimated_bulk is not None
                and estimated_bulk.diffusion_transverse_m2_s is not None
                else np.nan
            ),
            "estimated_source_gradient_velocity_m_s": float(estimated_gradient),
        }
        metadata.update(solution.metadata)
        if bulk is not None:
            metadata.update(
                {
                    "bulk_drift_velocity_m_s": bulk.drift_velocity_m_s,
                    "bulk_reduced_diffusion_L_m-1_s-1": (
                        bulk.diffusion_longitudinal_m2_s * number_density
                        if bulk.diffusion_longitudinal_m2_s is not None
                        else np.nan
                    ),
                    "bulk_reduced_diffusion_T_m-1_s-1": (
                        bulk.diffusion_transverse_m2_s * number_density
                        if bulk.diffusion_transverse_m2_s is not None
                        else np.nan
                    ),
                }
            )
        return SwarmCaseResult(
            solver=self.name,
            case_id=case_id,
            e_over_n_Td=case.e_over_n_Td,
            mean_energy_eV=solution.diagnostics.mean_energy_eV,
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
            rates=list(solution.rates.rates),
            metadata=metadata,
            transport=transport,
        )


__all__ = [
    "EnergyGrid",
    "DensityNormalizationConstraint",
    "LegendreBasis",
    "LegendreBlockLayout",
    "MultiTermBoltzmannSolver",
    "MultiTermCase",
    "MultiTermOperatorBackend",
    "MultiTermSolution",
    "OperatorBackendUnavailable",
    "OperatorAssemblyDiagnostics",
    "OperatorSolveState",
    "OperatorSystem",
    "assemble_operator_system",
    "build_density_normalization_constraint",
    "build_operator_assembly_diagnostics",
    "solve_operator_system",
    "build_multiterm_case",
    "druyvesteyn_energy_pdf",
    "electron_speed_m_s",
    "maxwellian_energy_pdf",
]
