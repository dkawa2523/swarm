"""Physical observables derived from a solved PN coefficient vector."""

from __future__ import annotations

import numpy as np

from electron_swarm.core.constants import TOWNSEND
from electron_swarm.core.result_metadata import (
    TRANSPORT_PN_F1_DRIFT_F0_DIFFUSION,
)
from electron_swarm.core.transport import ElectronTransport
from electron_swarm.solvers.boltzmann_common.grid import electron_speed_m_s
from electron_swarm.solvers.boltzmann_common.observables import (
    compute_rates_from_eedf,
    f0_reduced_transport_from_eedf,
    negative_mass_fraction,
)

from .models import (
    MultiTermCase,
    MultiTermSolution,
    PNOperator,
    RateSet,
    StationaryPNState,
)


def _rates_from_f0(
    case: MultiTermCase,
    operator: PNOperator,
    f0: np.ndarray,
    case_id: str,
    solver_name: str,
) -> RateSet:
    convolution = compute_rates_from_eedf(
        case.config,
        case.cross_sections,
        operator.energy_eV,
        operator.widths_eV,
        f0,
        case_id=case_id,
        e_over_n_Td=case.e_over_n_Td,
        solver_name=solver_name,
    )
    return RateSet(
        tuple(convolution.rates),
        case.gas_number_density_m3 * convolution.ionization_rate_m3_s,
        case.gas_number_density_m3 * convolution.attachment_rate_m3_s,
    )


def _transport_from_pn(
    case: MultiTermCase,
    operator: PNOperator,
    state: StationaryPNState,
) -> ElectronTransport:
    """Use F1 for flux drift and F0 only for the diffusion approximation."""

    f0 = state.coefficients[0]
    f1 = state.raw_coefficients[1]
    f1_energy = operator.moment_energy_eV[1]
    f1_weights = operator.moment_weights_eV[1]
    drift = float(
        np.sum(electron_speed_m_s(f1_energy) * f1 * f1_weights) / 3.0
    )
    reduced_field = case.e_over_n_Td * TOWNSEND
    if not np.isfinite(drift) or drift <= 0.0 or reduced_field <= 0.0:
        raise FloatingPointError(
            "PN F1 velocity moment did not produce positive finite drift"
        )
    reduced_mobility = drift / reduced_field
    _, reduced_diffusion, _, reduced_energy_diffusion = (
        f0_reduced_transport_from_eedf(
            operator.energy_eV,
            operator.widths_eV,
            f0,
            operator.sigma_m_m2,
        )
    )
    if (
        not np.isfinite(reduced_diffusion)
        or reduced_diffusion <= 0.0
        or not np.isfinite(reduced_energy_diffusion)
        or reduced_energy_diffusion <= 0.0
    ):
        raise FloatingPointError(
            "PN F0 diffusion moment is non-finite or non-positive"
        )
    return ElectronTransport(
        definition=TRANSPORT_PN_F1_DRIFT_F0_DIFFUSION,
        gas_number_density_m3=case.gas_number_density_m3,
        drift_velocity_m_s=drift,
        reduced_mobility_m2_V_s_m3=reduced_mobility,
        reduced_diffusion_L_m2_s_m3=float(reduced_diffusion),
        reduced_diffusion_T_m2_s_m3=float(reduced_diffusion),
        reduced_electron_energy_diffusion_m2_s_m3=(
            float(reduced_energy_diffusion)
        ),
    )


def build_solution(
    case: MultiTermCase,
    operator: PNOperator,
    state: StationaryPNState,
    case_id: str,
    solver_name: str,
) -> MultiTermSolution:
    coefficients = state.coefficients
    f0 = coefficients[0]
    rates = _rates_from_f0(case, operator, f0, case_id, solver_name)
    transport = _transport_from_pn(case, operator, state)
    method = case.multi_term_config.method
    angular_metadata = dict(operator.angular_metadata)
    exact_dcs = bool(angular_metadata.get("exact_dcs_based", False))
    metadata: dict[str, object] = {
        "solver_method": method,
        "physics_level": (
            "dcs_moment_direct_pn"
            if exact_dcs
            else (
                "moment_table_direct_pn"
                if method == "pn_dcs"
                else "ordinary_xs_angular_closure_direct_pn"
            )
        ),
        "direct_pn_operator": True,
        "exact_dcs_based": exact_dcs,
        "ordinary_integral_xs_closure": method == "pn_closure_direct",
        "lmax": operator.lmax,
        "angular_model": angular_metadata.get("angular_model", "unknown"),
        "angular_moment_source": angular_metadata.get(
            "angular_moment_source", "unknown"
        ),
        "pn_residual": state.diagnostics.full_pn_relative_residual,
        "pn_iterations": state.diagnostics.iterations,
        "pn_growth_frequency_s_inv": (
            state.diagnostics.growth_frequency_s_inv
        ),
        "field_coupling": "bidirectional_adjacent_legendre_moments",
        "drift_observable": "F1_velocity_moment",
        "diffusion_observable": "F0_gradient_reconstruction",
        "negative_mass_fraction": negative_mass_fraction(
            f0,
            operator.widths_eV,
        ),
    }
    if "moment_table_provenance" in angular_metadata:
        metadata["moment_table_provenance"] = angular_metadata[
            "moment_table_provenance"
        ]
    return MultiTermSolution(
        energy_eV=operator.energy_eV,
        widths_eV=operator.widths_eV,
        coefficients=coefficients,
        eedf_eV_inv=f0,
        rates=rates,
        transport=transport,
        diagnostics=state.diagnostics,
        metadata=metadata,
    )


__all__ = ["build_solution"]
