"""Two-term temporal-growth closure and transport projection."""

from __future__ import annotations

import numpy as np

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.constants import TOWNSEND
from electron_swarm.core.solver_configs import TwoTermInternalConfig
from electron_swarm.core.transport import ElectronTransport
from electron_swarm.solvers.boltzmann_common.grid import electron_speed_m_s
from electron_swarm.solvers.boltzmann_common.observables import (
    f0_reduced_transport_from_eedf,
    f0_reduced_transport_from_inverse_sigma,
    weighted_integral,
)
from electron_swarm.solvers.boltzmann_common.operators import (
    KineticOperatorBlock,
    assemble_energy_flux_operator,
)


TWO_TERM_TEMPORAL_GROWTH_TRANSPORT_SCHEMA_VERSION = (
    "two_term_temporal_growth_transport.v1"
)
TWO_TERM_TEMPORAL_GROWTH_EFFECTIVE_MOMENTUM_MODEL = (
    "nu_m_eff(epsilon)=nu_m(epsilon)+growth_frequency"
)


def temporal_growth_effective_momentum_frequency(
    momentum_frequency_s_inv: np.ndarray,
    growth_frequency_s_inv: float,
) -> np.ndarray:
    """Return the PT effective momentum frequency ``nu_m + growth``."""

    base = np.asarray(momentum_frequency_s_inv, dtype=float)
    growth = float(growth_frequency_s_inv)
    if base.ndim != 1 or not np.all(np.isfinite(base)) or np.any(base <= 0.0):
        raise FloatingPointError(
            "PT momentum correction requires finite positive nu_m [s^-1]"
        )
    if not np.isfinite(growth):
        raise FloatingPointError(
            "PT momentum correction requires a finite growth frequency [s^-1]"
        )
    effective = base + growth
    if not np.all(np.isfinite(effective)) or np.any(effective <= 0.0):
        raise FloatingPointError(
            "PT effective momentum frequency nu_m + growth must remain "
            "finite and positive at every energy"
        )
    return effective


def discrete_elastic_energy_loss_rate_coefficient(
    block: KineticOperatorBlock,
    eedf: np.ndarray,
) -> float:
    """Return the elastic energy-loss coefficient from the solved operator."""

    values = np.asarray(eedf, dtype=float)
    if values.shape != block.energy_eV.shape or not np.all(np.isfinite(values)):
        raise ValueError("elastic energy-loss moment requires a finite grid EEDF")
    density = float(block.gas_number_density_m3)
    if not np.isfinite(density) or density <= 0.0:
        raise ValueError("elastic energy-loss moment requires positive gas density")
    elastic_operator = assemble_energy_flux_operator(
        block.energy_eV,
        block.widths_eV,
        0.0,
        block.collisions,
    )
    energy_change_eV_s = weighted_integral(
        block.energy_eV * (elastic_operator @ values),
        block.widths_eV,
    )
    coefficient = -energy_change_eV_s / density
    if not np.isfinite(coefficient):
        raise FloatingPointError("elastic energy-loss operator moment is non-finite")
    return float(coefficient)


def transport_from_eedf(
    config: SwarmConfig,
    block: KineticOperatorBlock,
    eedf: np.ndarray,
    e_over_n_Td: float,
    solver_config: TwoTermInternalConfig,
    *,
    temporal_growth_frequency_s_inv: float | None = None,
) -> ElectronTransport:
    """Project transport from an EEDF using its solved collision data."""

    energy = block.energy_eV
    widths = block.widths_eV
    gas_number_density_m3 = block.gas_number_density_m3
    coll = block.collisions
    if temporal_growth_frequency_s_inv is None:
        moments = f0_reduced_transport_from_eedf(energy, widths, eedf, coll.sigma_m)
    else:
        if (
            config.physics.field.type != "dc"
            or solver_config.nonconservative_model != "growth"
        ):
            raise ValueError(
                "PT momentum correction is limited to DC two_term growth mode"
            )
        effective_frequency = temporal_growth_effective_momentum_frequency(
            coll.nu_m,
            temporal_growth_frequency_s_inv,
        )
        inverse_sigma = (
            gas_number_density_m3 * electron_speed_m_s(energy) / effective_frequency
        )
        moments = f0_reduced_transport_from_inverse_sigma(
            np.asarray(energy, dtype=float),
            np.asarray(widths, dtype=float),
            np.asarray(eedf, dtype=float),
            inverse_sigma,
        )
    muN, diffN, energy_muN, energy_diffN = moments
    transport_moments = {
        "mobility": muN,
        "diffusion": diffN,
        "energy_mobility": energy_muN,
        "energy_diffusion": energy_diffN,
    }
    invalid = [
        name
        for name, value in transport_moments.items()
        if not np.isfinite(value) or value <= 0.0
    ]
    if invalid:
        raise FloatingPointError(
            "two-term EEDF transport moment is nonfinite or nonpositive: "
            + ", ".join(invalid)
            + "; refusing to replace kinetic output with a Drude fallback"
        )
    EN = e_over_n_Td * TOWNSEND
    return ElectronTransport(
        definition="flux",
        gas_number_density_m3=gas_number_density_m3,
        drift_velocity_m_s=float(muN) * EN,
        reduced_mobility_m2_V_s_m3=float(muN),
        reduced_diffusion_L_m2_s_m3=float(diffN),
        reduced_diffusion_T_m2_s_m3=float(diffN),
        reduced_electron_energy_mobility_m2_V_s_m3=float(energy_muN),
        reduced_electron_energy_diffusion_m2_s_m3=float(energy_diffN),
    )


__all__ = [
    "TWO_TERM_TEMPORAL_GROWTH_EFFECTIVE_MOMENTUM_MODEL",
    "TWO_TERM_TEMPORAL_GROWTH_TRANSPORT_SCHEMA_VERSION",
    "discrete_elastic_energy_loss_rate_coefficient",
    "temporal_growth_effective_momentum_frequency",
    "transport_from_eedf",
]
