"""Observable meanings; algorithm names never imply interchangeable quantities."""

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class ObservableDefinition:
    quantity: str
    method: str
    unit: str
    approximate: bool = False


_FLUX_DRIFT = ObservableDefinition("flux_drift", "velocity_moment", "m/s")
_FLUX_MOBILITY = ObservableDefinition(
    "flux_mobility", "flux_drift_over_field", "m2/(V s)"
)
_F0_L = ObservableDefinition(
    "isotropic_particle_diffusion_approximation", "f0_gradient", "m2/s", True
)
_F0_T = ObservableDefinition("transverse_flux_diffusion", "f0_gradient", "m2/s", True)
_TRAJECTORY_L = ObservableDefinition(
    "longitudinal_flux_diffusion", "trajectory_response", "m2/s"
)
_TRAJECTORY_T = ObservableDefinition(
    "transverse_flux_diffusion", "trajectory_response", "m2/s"
)
_RESTRICTED_L = ObservableDefinition(
    "restricted_density_packet_diffusion_L",
    "weighted_trajectory_response",
    "m2/s",
    True,
)
_RESTRICTED_T = ObservableDefinition(
    "restricted_density_packet_diffusion_T",
    "weighted_trajectory_response",
    "m2/s",
    True,
)

TRANSPORT_OBSERVABLES = {
    "flux": (_F0_L, _F0_T),
    "pn_f1_flux_drift_f0_gradient_diffusion": (_F0_L, _F0_T),
    "mc_flux_particle_tracking_fixed_population": (_TRAJECTORY_L, _TRAJECTORY_T),
    "mc_flux_particle_tracking_weighted_growth_population": (
        _RESTRICTED_L,
        _RESTRICTED_T,
    ),
    "propagator_flux_velocity_moment": (None, None),
}

_F0_ENERGY = {
    "electron_energy_mobility_m2_V_s": ObservableDefinition(
        "mean_energy_normalized_flux_energy_mobility", "f1_energy_moment", "m2/(V s)"
    ),
    "electron_energy_diffusion_m2_s": ObservableDefinition(
        "isotropic_energy_diffusion_approximation", "f0_energy_gradient", "m2/s", True
    ),
    "electron_energy_diffusion_L_m2_s": ObservableDefinition(
        "isotropic_energy_diffusion_approximation", "f0_energy_gradient", "m2/s", True
    ),
    "electron_energy_diffusion_T_m2_s": ObservableDefinition(
        "transverse_flux_energy_diffusion", "f0_energy_gradient", "m2/s", True
    ),
}
_TRACKING_ENERGY = {
    "electron_energy_mobility_m2_V_s": ObservableDefinition(
        "mean_energy_normalized_flux_energy_mobility",
        "trajectory_energy_current",
        "m2/(V s)",
    ),
    "electron_energy_diffusion_L_m2_s": ObservableDefinition(
        "longitudinal_flux_energy_diffusion", "trajectory_response", "m2/s"
    ),
    "electron_energy_diffusion_T_m2_s": ObservableDefinition(
        "transverse_flux_energy_diffusion", "trajectory_response", "m2/s"
    ),
}
_WEIGHTED_ENERGY = {
    "electron_energy_mobility_m2_V_s": _TRACKING_ENERGY[
        "electron_energy_mobility_m2_V_s"
    ],
    "electron_energy_diffusion_L_m2_s": ObservableDefinition(
        "restricted_density_packet_energy_diffusion_L",
        "weighted_trajectory_response",
        "m2/s",
        True,
    ),
    "electron_energy_diffusion_T_m2_s": ObservableDefinition(
        "restricted_density_packet_energy_diffusion_T",
        "weighted_trajectory_response",
        "m2/s",
        True,
    ),
}
ENERGY_TRANSPORT_OBSERVABLES = {
    "flux": _F0_ENERGY,
    "pn_f1_flux_drift_f0_gradient_diffusion": _F0_ENERGY,
    "mc_flux_particle_tracking_fixed_population": _TRACKING_ENERGY,
    "mc_flux_particle_tracking_weighted_growth_population": _WEIGHTED_ENERGY,
    "propagator_flux_velocity_moment": {},
}


def transport_observables(definition: str) -> dict[str, ObservableDefinition]:
    if definition not in TRANSPORT_OBSERVABLES:
        return {}
    longitudinal, transverse = TRANSPORT_OBSERVABLES[definition]
    values = {"drift_velocity_m_s": _FLUX_DRIFT, "mobility_m2_V_s": _FLUX_MOBILITY}
    values.update(ENERGY_TRANSPORT_OBSERVABLES[definition])
    if longitudinal is not None:
        values["diffusion_L_m2_s"] = longitudinal
    if transverse is not None:
        values["diffusion_T_m2_s"] = transverse
    return values


def observable_status(candidate, reference) -> str:
    if candidate is None or reference is None:
        return "definition_unknown"
    if candidate.quantity != reference.quantity or candidate.unit != reference.unit:
        return "definition_mismatch"
    return (
        "available_model_approximation"
        if candidate.approximate or reference.approximate
        else "available"
    )
