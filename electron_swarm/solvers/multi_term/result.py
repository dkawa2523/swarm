"""Projection of an internal PN solution into the canonical product result."""

from __future__ import annotations

import numpy as np

from electron_swarm.core.result_metadata import (
    TRANSPORT_PN_F1_DRIFT_F0_DIFFUSION,
)
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.solvers.boltzmann_common.observables import weighted_integral

from .models import MultiTermCase, MultiTermSolution


def to_case_result(
    case: MultiTermCase,
    solution: MultiTermSolution,
    case_id: str,
    run_time_s: float,
    solver_name: str,
) -> SwarmCaseResult:
    transport = solution.transport
    number_density = case.gas_number_density_m3
    eedf = np.asarray(solution.eedf_eV_inv, dtype=float)
    rates = solution.rates
    net_frequency = rates.effective_growth_frequency_s_inv
    diagnostic = solution.diagnostics
    energy_grid = case.multi_term_config.energy_grid
    details = {
        "method": case.multi_term_config.method,
        "run_label": case_id,
        "sweep_param": "E_over_N_Td",
        "sweep_value": float(case.e_over_n_Td),
        "run_time_s": float(run_time_s),
        "gas_number_density_m-3": float(number_density),
        "electric_field_V_m": float(case.electric_field_V_m),
        "lmax": int(case.multi_term_config.lmax),
        "direct_pn_operator": True,
        "transport_definition": TRANSPORT_PN_F1_DRIFT_F0_DIFFUSION,
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
        "tail_fraction": diagnostic.eedf_tail_fraction,
        "pn_converged": True,
        "pn_iterations": diagnostic.iterations,
        "pn_growth_frequency_s_inv": diagnostic.growth_frequency_s_inv,
        "pn_full_relative_residual": diagnostic.full_pn_relative_residual,
        "pn_shape_change": diagnostic.shape_change,
        "pn_eigenvalue_change": diagnostic.eigenvalue_change,
        "normalization_integral": weighted_integral(
            eedf, solution.widths_eV
        ),
        "convolution_ionization_rate_coefficient_m3_s": float(
            sum(
                rate.mixture_weighted_rate_m3_s
                for rate in rates.rates
                if rate.process_type == "ionization"
            )
        ),
        "convolution_attachment_rate_coefficient_m3_s": float(
            sum(
                rate.mixture_weighted_rate_m3_s
                for rate in rates.rates
                if rate.process_type == "attachment"
            )
        ),
        "convolution_effective_rate_coefficient_m3_s": float(
            net_frequency / number_density if number_density > 0.0 else np.nan
        ),
        "convolution_net_ionization_frequency_s-1": float(net_frequency),
        "net_ionization_frequency_model": "convolution_effective_rate",
        "effective_townsend_1_m": float(
            net_frequency / max(abs(transport.drift_velocity_m_s), 1.0e-300)
        ),
    }
    details.update(solution.metadata)
    metadata = {
        key: solution.metadata[key]
        for key in (
            "solver_method",
            "angular_model",
            "angular_moment_source",
            "moment_table_provenance",
            "exact_dcs_based",
            "ordinary_integral_xs_closure",
            "lmax",
            "direct_pn_operator",
        )
        if key in solution.metadata
    }
    metadata["transport_definition"] = TRANSPORT_PN_F1_DRIFT_F0_DIFFUSION
    return SwarmCaseResult(
        solver=solver_name,
        case_id=case_id,
        e_over_n_Td=case.e_over_n_Td,
        mean_energy_eV=float(diagnostic.mean_energy_eV),
        net_ionization_frequency_s=net_frequency,
        effective_townsend_m2=(
            net_frequency
            / max(abs(transport.drift_velocity_m_s) * number_density, 1.0e-300)
        ),
        transport=transport,
        energy_eV=solution.energy_eV,
        eedf=eedf,
        energy_widths_eV=solution.widths_eV,
        rates=list(rates.rates),
        metadata=metadata,
        diagnostics={"multi_term": details},
    )


__all__ = ["to_case_result"]
