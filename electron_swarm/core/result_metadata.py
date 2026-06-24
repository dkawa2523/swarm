"""Shared product-result metadata names and compact value labels."""

from __future__ import annotations

from collections.abc import Iterable

from electron_swarm.core.results import SwarmCaseResult


SUMMARY_METADATA_KEYS = (
    "angular_model",
    "angular_moment_source",
    "moment_table_provenance",
    "exact_dcs_based",
    "ordinary_integral_xs_closure",
    "lmax",
    "ionization_source_treatment",
    "electron_electron_treatment",
    "electron_electron_transport_stale",
    "magnetic_field_treatment",
    "tail_refinement_treatment",
    "transport_definition",
)

PRODUCT_CASE_METADATA_KEYS = frozenset(
    {
        "schema_version",
        "solver_method",
        "angular_model",
        "angular_moment_source",
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
        "tail_refinement_treatment",
        "transport_definition",
    }
)

TRANSPORT_FLUX = "flux"
TRANSPORT_F0_GRADIENT_RECONSTRUCTION = "f0_gradient_reconstruction"
TRANSPORT_MC_FIXED_POPULATION = "mc_flux_particle_tracking_fixed_population"
TRANSPORT_MC_WEIGHTED_GROWTH = "mc_flux_particle_tracking_weighted_growth_population"


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
