"""Shared product-result metadata names and compact value labels."""

from __future__ import annotations

from collections.abc import Iterable

from electron_swarm.core.results import SwarmCaseResult


SUMMARY_METADATA_KEYS = (
    "angular_model",
    "angular_moment_source",
    "angular_scattering_treatment",
    "angular_scattering_fidelity",
    "angular_scattering_assumption",
    "moment_table_provenance",
    "exact_dcs_based",
    "ordinary_integral_xs_closure",
    "lmax",
    "ionization_source_treatment",
    "electron_electron_treatment",
    "electron_electron_transport_stale",
    "magnetic_field_treatment",
    "rf_field_treatment",
    "rf_frequency_Hz",
    "rf_amplitude_definition",
    "tail_refinement_treatment",
    "monte_carlo_base_seed",
    "monte_carlo_case_seed",
    "transport_definition",
    "velocity_space_representation",
    "transport_components",
    "inelastic_angular_model",
    "inelastic_radial_transfer",
    "elastic_recoil_model",
    "elastic_collision_xs_role",
)

PRODUCT_CASE_METADATA_KEYS = frozenset(
    {
        "schema_version",
        "solver_method",
        "angular_model",
        "angular_moment_source",
        "angular_scattering_treatment",
        "angular_scattering_fidelity",
        "angular_scattering_assumption",
        "moment_table_provenance",
        "exact_dcs_based",
        "ordinary_integral_xs_closure",
        "lmax",
        "direct_pn_operator",
        "ionization_source_model",
        "ionization_source_treatment",
        "ionization_secondary_electron_energy_eV",
        "electron_electron_treatment",
        "electron_electron_transport_stale",
        "magnetic_field_treatment",
        "rf_field_treatment",
        "rf_frequency_Hz",
        "rf_amplitude_definition",
        "tail_refinement_treatment",
        "monte_carlo_base_seed",
        "monte_carlo_case_seed",
        "transport_definition",
        "velocity_space_representation",
        "transport_components",
        "inelastic_angular_model",
        "inelastic_radial_transfer",
        "elastic_recoil_model",
        "elastic_collision_xs_role",
        "elastic_total_xs_process_ids",
        "elastic_momentum_xs_process_ids",
    }
)

TRANSPORT_FLUX = "flux"
TRANSPORT_PN_F1_DRIFT_F0_DIFFUSION = "pn_f1_flux_drift_f0_gradient_diffusion"
TRANSPORT_MC_FIXED_POPULATION = "mc_flux_particle_tracking_fixed_population"
TRANSPORT_MC_WEIGHTED_GROWTH = "mc_flux_particle_tracking_weighted_growth_population"
TRANSPORT_PROPAGATOR_FLUX = "propagator_flux_velocity_moment"


def normalize_case_metadata(case: SwarmCaseResult) -> SwarmCaseResult:
    """Keep only product-contract metadata on a case result."""

    case.metadata = {
        key: value
        for key, value in case.metadata.items()
        if key in PRODUCT_CASE_METADATA_KEYS
    }
    return case


def normalize_product_metadata(
    cases: Iterable[SwarmCaseResult],
) -> list[SwarmCaseResult]:
    return [normalize_case_metadata(case) for case in cases]
