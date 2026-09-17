from __future__ import annotations

from pathlib import Path

from swarm_workflow.comsol.models.gec_icp.mapping import load_gec_icp_mapping
from swarm_workflow.comsol.models.gec_icp.run_mapping import (
    load_gec_icp_run_mapping,
)


ROOT = Path(__file__).resolve().parents[1]
MAPS = ROOT / "comsol_modes" / "maps"
RUN_MAPS = {
    "two_term": "argon_gec_icp_two_term_function_eedf.yaml",
    "monte_carlo": ("argon_gec_icp_monte_carlo_function_eedf_restricted_lmea.yaml"),
    "propagator": "argon_gec_icp_propagator_function_eedf.yaml",
}


def test_production_model_map_matches_repository_icp_topology() -> None:
    mapping = load_gec_icp_mapping(
        MAPS / "argon_gec_icp_model.yaml", validate_files=True
    )

    assert mapping.model.input_mph == ROOT / "comsol_modes" / "argon_gec_icp.mph"
    assert mapping.model.component == "comp1"
    assert mapping.model.plasma_physics == "plas"
    assert mapping.model.magnetic_physics == "mf"
    assert mapping.model.study_feature == "ftrans"
    assert mapping.reactions.superelastic == "eir3"
    assert mapping.reactions.stepwise_ionization == "eir5"


def test_production_run_maps_share_one_restricted_lmea_contract() -> None:
    loaded = {
        source: load_gec_icp_run_mapping(MAPS / filename)
        for source, filename in RUN_MAPS.items()
    }

    assert set(loaded) == {
        mapping.bundle.expected_source for mapping in loaded.values()
    }
    assert len({mapping.run.output_mph for mapping in loaded.values()}) == 3
    assert len({mapping.output_directory for mapping in loaded.values()}) == 3
    assert len({mapping.log_path for mapping in loaded.values()}) == 3
    for mapping in loaded.values():
        assert mapping.closure.source_field == "steady_dc"
        assert mapping.closure.electron_transport == ("swarm_mobility_comsol_einstein")
        assert mapping.closure.reaction_model == "function_eedf"
        assert mapping.closure.elastic_energy_loss_model == (
            "comsol_cross_section_integral"
        )
        assert mapping.closure.chemistry_owner == "comsol_embedded_cross_sections"
        assert mapping.closure.rf_owner == "comsol_frequency_transient"
        assert mapping.run.power_W == 1500.0
        assert mapping.run.frequency_Hz == 13_560_000.0
        assert mapping.run.pressure_Pa == 2.66644
        assert mapping.run.gas_temperature_K == 300.0
        assert mapping.run.final_time_s == 0.01
        assert mapping.run.output_points_per_decade == 4
        assert mapping.run.clear_saved_solution is True
        assert mapping.run.convergence.mean_energy_support_relative_guard == 0.10

    assert loaded["two_term"].bundle.expected_transport_definition == "flux"
    assert loaded["monte_carlo"].bundle.expected_transport_definition == (
        "mc_flux_particle_tracking_weighted_growth_population"
    )
    assert loaded["monte_carlo"].bundle.expected_mc_qualification_profile == (
        "function_eedf_restricted_lmea"
    )
    assert loaded["propagator"].bundle.expected_transport_definition == (
        "propagator_flux_velocity_moment"
    )
