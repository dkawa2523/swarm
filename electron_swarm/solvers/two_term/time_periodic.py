"""Time-periodic two-term EEDF propagation."""

from __future__ import annotations

import math

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as spla

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.constants import TOWNSEND
from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.core.results import RfPhaseResult
from electron_swarm.core.solver_configs import TwoTermInternalConfig
from electron_swarm.solvers.boltzmann_common.observables import (
    compute_rates_from_eedf,
    normalize_eedf,
    weighted_integral,
)
from electron_swarm.solvers.boltzmann_common.operators import (
    assemble_energy_flux_operator,
)
from electron_swarm.solvers.two_term.models import TimePeriodicDistributionResult
from electron_swarm.solvers.two_term.steady import (
    mean_energy,
    solve_native_distribution,
    tail_probability,
)


_EPS = 1.0e-300


def solve_time_periodic_distribution(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    cfg: TwoTermInternalConfig,
    e_over_n_rms_Td: float,
    *,
    case_id: str = "rf_case",
    solver_name: str = "two_term",
) -> TimePeriodicDistributionResult:
    """Solve a periodic F0 with sinusoidal RMS E/N.

    The energy-space distribution is advanced with backward Euler over one RF
    cycle. The anisotropic two-term response is assumed to follow the electric
    field instantaneously, so the field-heating operator is evaluated at every
    phase. This is a controlled approximation, not a full time-dependent F0/F1
    solver.
    """

    field = config.physics.field.time_dependent
    if field.frequency_Hz is None or field.frequency_Hz <= 0.0:
        raise ValueError("time-dependent field frequency_Hz must be positive")
    if field.waveform != "sinusoidal":
        raise NotImplementedError(
            f"time-dependent waveform {field.waveform!r} is not implemented"
        )
    if field.amplitude_definition != "rms":
        raise NotImplementedError(
            "time-dependent two_term currently requires RMS E/N amplitudes"
        )
    if field.momentum_response != "instantaneous":
        raise NotImplementedError(
            "time-dependent two_term currently requires instantaneous F1 response"
        )

    # Size the adaptive grid at the sinusoidal peak, not at RMS. Using the RMS
    # grid under-resolves the high-energy tail during field peaks even when the
    # periodic propagation itself has converged.
    e_over_n_peak_Td = math.sqrt(2.0) * float(e_over_n_rms_Td)
    grid_seed = solve_native_distribution(
        config,
        cross_sections,
        cfg,
        e_over_n_peak_Td,
    )
    state_seed = solve_native_distribution(
        config,
        cross_sections,
        cfg,
        e_over_n_rms_Td,
    )
    energy = grid_seed.energy_eV
    widths = grid_seed.widths_eV
    operator_block = grid_seed.operator_block
    collisions = operator_block.collisions
    collision_matrix = operator_block.collision_matrix
    gas_density = operator_block.gas_number_density_m3

    phase_steps = field.phase_steps
    frequency_Hz = float(field.frequency_Hz)
    period_s = 1.0 / frequency_Hz
    dt_s = period_s / phase_steps
    phase_rad = (np.arange(phase_steps, dtype=float) + 0.5) * (
        2.0 * math.pi / phase_steps
    )
    instantaneous_en = math.sqrt(2.0) * float(e_over_n_rms_Td) * np.sin(phase_rad)
    identity = sparse.identity(len(energy), format="csc")
    implicit_steps: list[object] = []
    for e_over_n_Td in instantaneous_en:
        electric_field = float(e_over_n_Td) * TOWNSEND * gas_density
        energy_flux = assemble_energy_flux_operator(
            energy,
            widths,
            electric_field,
            collisions,
        )
        operator = (energy_flux + collision_matrix).tocsc()
        implicit_steps.append(spla.splu(identity - dt_s * operator))

    def propagate_period(values: np.ndarray) -> np.ndarray:
        propagated = np.asarray(values, dtype=float)
        for factorization in implicit_steps:
            propagated = factorization.solve(propagated)
        return propagated

    state = normalize_eedf(
        np.clip(
            np.interp(
                energy,
                state_seed.energy_eV,
                state_seed.eedf_eV_inv,
                left=state_seed.eedf_eV_inv[0],
                right=0.0,
            ),
            0.0,
            None,
        ),
        widths,
    )
    periodic_residual = math.inf
    cycle_multiplier = 1.0
    periods = 0
    converged = False
    for period_index in range(field.max_periods):
        periods = period_index + 1
        end_state = propagate_period(state)
        cycle_multiplier = weighted_integral(end_state, widths)
        if not np.isfinite(cycle_multiplier) or cycle_multiplier <= 0.0:
            raise FloatingPointError(
                "time-periodic Boltzmann propagation returned invalid mass"
            )
        end_state = normalize_eedf(
            np.clip(end_state, 0.0, None),
            widths,
        )
        periodic_residual = float(np.sum(np.abs(end_state - state) * widths))
        state = end_state
        if periodic_residual <= field.periodic_tolerance:
            converged = True
            break

    phase_eedf: list[np.ndarray] = []
    phase_mean_energy: list[float] = []
    phase_ionization: list[float] = []
    current = state
    for phase_index, factorization in enumerate(implicit_steps):
        next_state = factorization.solve(current)
        next_state = normalize_eedf(
            np.clip(next_state, 0.0, None),
            widths,
        )
        midpoint = normalize_eedf(
            np.clip(0.5 * (current + next_state), 0.0, None),
            widths,
        )
        phase_eedf.append(midpoint)
        phase_mean_energy.append(mean_energy(energy, midpoint, widths))
        convolution = compute_rates_from_eedf(
            config,
            cross_sections,
            energy,
            widths,
            midpoint,
            case_id=f"{case_id}_phase_{phase_index:04d}",
            e_over_n_Td=abs(float(instantaneous_en[phase_index])),
            solver_name=solver_name,
        )
        phase_ionization.append(convolution.ionization_rate_m3_s)
        current = next_state

    phase_matrix = np.asarray(phase_eedf, dtype=float)
    cycle_eedf = normalize_eedf(
        np.mean(phase_matrix, axis=0),
        widths,
    )
    phase_mean = np.asarray(phase_mean_energy, dtype=float)
    phase_l1 = np.sum(
        np.abs(phase_matrix - cycle_eedf[np.newaxis, :]) * widths[np.newaxis, :],
        axis=1,
    )
    hysteresis_l1 = 0.0
    for index in range(phase_steps):
        mirror = phase_steps - 1 - index
        hysteresis_l1 = max(
            hysteresis_l1,
            float(np.sum(np.abs(phase_matrix[index] - phase_matrix[mirror]) * widths)),
        )

    nu_eff = weighted_integral(collisions.nu_m * cycle_eedf, widths)
    omega = 2.0 * math.pi * frequency_Hz
    tail = tail_probability(cycle_eedf, widths, cfg)
    edge_to_peak = float(cycle_eedf[-1] / max(np.max(cycle_eedf), _EPS))
    diagnostics: dict[str, object] = {
        **grid_seed.metadata,
        "backend": "native_sg_time_periodic_f0",
        "field_type": "time_dependent",
        "field_waveform": "sinusoidal",
        "field_amplitude_definition": "rms",
        "adaptive_grid_seed": "sinusoidal_peak",
        "adaptive_grid_seed_E_over_N_Td": e_over_n_peak_Td,
        "periodic_state_seed": "stationary_rms",
        "rf_frequency_Hz": frequency_Hz,
        "reduced_angular_frequency_m3_s": omega / gas_density,
        "momentum_response": "instantaneous_f1",
        "time_integrator": "backward_euler",
        "phase_steps": phase_steps,
        "period_s": period_s,
        "periods": periods,
        "iterations": periods,
        "converged": converged,
        "periodic_residual_L1": periodic_residual,
        "residual_L1": periodic_residual,
        "periodic_tolerance": field.periodic_tolerance,
        "residual_tolerance": field.periodic_tolerance,
        "cycle_growth_frequency_s-1": math.log(cycle_multiplier) / period_s,
        "effective_momentum_frequency_s-1": nu_eff,
        "omega_over_effective_momentum_frequency": omega / max(nu_eff, _EPS),
        "phase_mean_energy_min_eV": float(np.min(phase_mean)),
        "phase_mean_energy_max_eV": float(np.max(phase_mean)),
        "phase_mean_energy_modulation_fraction": float(
            (np.max(phase_mean) - np.min(phase_mean)) / max(np.mean(phase_mean), _EPS)
        ),
        "max_phase_to_cycle_eedf_L1": float(np.max(phase_l1)),
        "same_field_hysteresis_eedf_L1": hysteresis_l1,
        "tail_probability": tail,
        "edge_to_peak": edge_to_peak,
        "transport_model": "cycle_averaged_f0_flux_integral",
    }
    phase = RfPhaseResult(
        phase_fraction=phase_rad / (2.0 * math.pi),
        phase_rad=phase_rad,
        instantaneous_e_over_n_Td=instantaneous_en,
        mean_energy_eV=phase_mean,
        ionization_rate_coefficient_m3_s=np.asarray(
            phase_ionization,
            dtype=float,
        ),
        eedf=phase_matrix,
    )
    return TimePeriodicDistributionResult(
        energy_eV=energy,
        widths_eV=widths,
        cycle_averaged_eedf_eV_inv=cycle_eedf,
        phase=phase,
        diagnostics=diagnostics,
        operator_block=operator_block,
    )
