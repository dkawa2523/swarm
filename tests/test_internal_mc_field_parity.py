from __future__ import annotations

from copy import deepcopy

import numpy as np
import pytest

from electron_swarm import SwarmConfig, load_config, run
from electron_swarm.solvers.monte_carlo.evidence import (
    FIELD_PARITY_RESPONSE_ESTIMATOR_SCHEMA_VERSION,
)


def _small_weighted_config() -> SwarmConfig:
    config = load_config("examples/argon_gec_icp_monte_carlo_weighted_branching.yaml")
    config.run.e_over_n_Td = [20.0]
    config.physics.energy_grid_policy.threshold_refinement = False
    mc = config.solvers.monte_carlo
    mc.seed = 12
    mc.particles = 8
    mc.warmup_collisions = 8
    mc.max_collisions = 32
    mc.tail_max_collisions = 8
    mc.transport_correlation_lag_barriers = 4
    return config


@pytest.mark.mc
def test_field_parity_pairs_only_odd_transport_and_preserves_primary_path() -> None:
    single_config = _small_weighted_config()
    paired_config = deepcopy(single_config)
    paired_config.solvers.monte_carlo.transport_estimator = "paired_field_parity"

    single = run(single_config, write=False).cases[0]
    paired = run(paired_config, write=False).cases[0]
    repeated = run(paired_config, write=False).cases[0]

    # Cloning the generator before the primary leg consumes no random values:
    # the positive-field EEDF/rates and even transport remain bitwise stable.
    assert np.array_equal(paired.energy_eV, single.energy_eV)
    assert np.array_equal(paired.eedf, single.eedf)
    assert np.array_equal(paired.eedf_counts, single.eedf_counts)
    assert paired.mean_energy_eV == single.mean_energy_eV
    assert paired.transport.diffusion_L_m2_s == single.transport.diffusion_L_m2_s
    assert paired.transport.diffusion_T_m2_s == single.transport.diffusion_T_m2_s
    assert [rate.rate_coefficient_m3_s for rate in paired.rates] == [
        rate.rate_coefficient_m3_s for rate in single.rates
    ]

    diagnostics = paired.diagnostics["internal_monte_carlo_transport"]
    parity = diagnostics["field_parity_response"]
    assert diagnostics["component_estimators"] == {
        "drift_mobility": "normalized_weight_residence_velocity_moment",
        "particle_diffusion": (
            "synchronized_block_lag_position_velocity_flux_covariance"
        ),
        "energy_mobility": "normalized_weight_residence_energy_current",
        "energy_diffusion": (
            "synchronized_block_lag_restricted_density_packet_energy_current_covariance"
        ),
    }
    assert parity["estimator_schema_version"] == (
        FIELD_PARITY_RESPONSE_ESTIMATOR_SCHEMA_VERSION
    )
    for name, value in parity["paired_response"].items():
        assert value == pytest.approx(
            0.5
            * (
                parity["primary_positive_field"][name]
                + parity["mirror_negative_field"][name]
            )
        )
        assert diagnostics["production_estimates"][name] == pytest.approx(value)
    assert parity["eedf_rates_even_transport_source"] == ("primary_positive_field_leg")
    assert parity["primary_sampling_stream"] == (
        "bitwise_identical_to_single_field_same_seed_and_budget"
    )
    assert parity["mirror_tail_sampling"] is False
    assert paired.effective_townsend_m2 == pytest.approx(
        paired.net_ionization_frequency_s
        / (
            abs(paired.transport.drift_velocity_m_s)
            * paired.transport.gas_number_density_m3
        )
    )
    assert diagnostics["mc_run_provenance"]["batch_rng_advance"] == (
        "primary_positive_field_leg_only"
    )

    assert repeated.transport == paired.transport
    assert np.array_equal(repeated.eedf, paired.eedf)
    assert (
        repeated.diagnostics["internal_monte_carlo_transport"]["field_parity_response"]
        == parity
    )


def test_field_parity_rejects_a_magnetic_execution() -> None:
    config = _small_weighted_config()
    config.solvers.monte_carlo.transport_estimator = "paired_field_parity"
    magnetic = config.physics.field.magnetic_field
    magnetic.enabled = True
    magnetic.B_T = 0.01

    with pytest.raises(NotImplementedError, match="requires a DC electric field"):
        run(config, write=False)


@pytest.mark.mc
def test_field_parity_does_not_advance_the_primary_batch_stream() -> None:
    single_config = _small_weighted_config()
    single_config.run.e_over_n_Td = [20.0, 30.0]
    paired_config = deepcopy(single_config)
    paired_config.solvers.monte_carlo.transport_estimator = "paired_field_parity"

    single_cases = run(single_config, write=False).cases
    paired_cases = run(paired_config, write=False).cases

    for single, paired in zip(single_cases, paired_cases, strict=True):
        assert np.array_equal(paired.eedf, single.eedf)
        assert paired.mean_energy_eV == single.mean_energy_eV
        assert [rate.rate_coefficient_m3_s for rate in paired.rates] == [
            rate.rate_coefficient_m3_s for rate in single.rates
        ]
