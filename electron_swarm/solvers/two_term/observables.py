"""Transport, reaction-rate, and case-result assembly for two-term solves."""

from __future__ import annotations

import numpy as np

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.core.solver_configs import TwoTermInternalConfig
from electron_swarm.solvers.boltzmann_common.observables import (
    compute_rates_from_eedf,
    mean_energy_from_eedf,
    weighted_integral,
)
from electron_swarm.solvers.boltzmann_common.operators import KineticOperatorBlock
from electron_swarm.solvers.two_term.transport import (
    discrete_elastic_energy_loss_rate_coefficient,
    transport_from_eedf,
)


def build_case_result(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    cfg: TwoTermInternalConfig,
    e_over_n_Td: float,
    case_id: str,
    operator_block: KineticOperatorBlock,
    eedf: np.ndarray,
    *,
    metadata: dict[str, object],
    solver_name: str,
    temporal_growth_frequency_s_inv: float | None = None,
) -> SwarmCaseResult:
    energy = operator_block.energy_eV
    widths = operator_block.widths_eV
    density = operator_block.gas_number_density_m3
    transport = transport_from_eedf(
        config,
        operator_block,
        eedf,
        e_over_n_Td,
        cfg,
        temporal_growth_frequency_s_inv=temporal_growth_frequency_s_inv,
    )
    mean_energy = mean_energy_from_eedf(energy, widths, eedf)
    elastic_energy_loss = discrete_elastic_energy_loss_rate_coefficient(
        operator_block,
        eedf,
    )
    rate_data = compute_rates_from_eedf(
        config,
        cross_sections,
        energy,
        widths,
        eedf,
        case_id=case_id,
        e_over_n_Td=e_over_n_Td,
        solver_name=solver_name,
    )
    rates = rate_data.rates
    ion_rate = rate_data.ionization_rate_m3_s
    attach_rate = rate_data.attachment_rate_m3_s
    net_ion_freq = rate_data.net_ionization_frequency_s
    effective_townsend = net_ion_freq / max(
        abs(transport.drift_velocity_m_s) * density,
        1.0e-300,
    )
    metadata = dict(metadata)
    metadata["characteristic_energy_eV"] = transport.characteristic_energy_L_eV
    metadata["normalization_integral"] = weighted_integral(eedf, widths)
    metadata["convolution_ionization_rate_coefficient_m3_s"] = ion_rate
    metadata["convolution_attachment_rate_coefficient_m3_s"] = attach_rate
    metadata["convolution_effective_rate_coefficient_m3_s"] = ion_rate - attach_rate
    metadata["convolution_net_ionization_frequency_s-1"] = net_ion_freq
    metadata["net_ionization_frequency_model"] = "convolution_effective_rate"
    metadata["effective_townsend_1_m"] = net_ion_freq / max(
        abs(transport.drift_velocity_m_s),
        1.0e-300,
    )
    metadata["elastic_energy_loss"] = {
        "schema": "swarm.elastic_energy_loss.v1",
        "status": "available",
        "symbol": "K_epsilon_el",
        "rate_coefficient_eV_m3_s": elastic_energy_loss,
        "estimator": "same_discrete_elastic_collision_operator_energy_moment",
        "operator": "native_finite_volume_scharfetter_gummel_elastic_A_D_zero_field",
        "operator_isolation": (
            "zero_field_reassembly_from_same_elastic_A_D_coefficients"
        ),
        "source_eedf": "same_solved_eedf",
        "gas_temperature_K": float(config.conditions.gas_temperature_K),
        "neutral_thermal_motion_model": "finite_temperature_fokker_planck",
        "gas_temperature_terms_included": True,
        "sign_convention": "positive_is_net_electron_energy_loss",
        "uncertainty_status": "deterministic_kinetic_solver",
    }
    diagnostics = {"two_term": metadata}
    return SwarmCaseResult(
        solver=solver_name,
        case_id=case_id,
        e_over_n_Td=e_over_n_Td,
        mean_energy_eV=mean_energy,
        net_ionization_frequency_s=net_ion_freq,
        effective_townsend_m2=effective_townsend,
        transport=transport,
        energy_eV=energy,
        eedf=eedf,
        energy_widths_eV=widths,
        rates=rates,
        metadata={},
        diagnostics=diagnostics,
    )
