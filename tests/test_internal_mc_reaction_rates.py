from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from electron_swarm import load_config, run
from electron_swarm.core.constants import ELECTRON_MASS_KG, EV_TO_J
from electron_swarm.core.cross_sections import CrossSectionProcess, ProcessType
from electron_swarm.physics.kinetics import cell_integrated_rate_coefficient
from electron_swarm.solvers.monte_carlo.reaction_rates import (
    MIN_REACTION_RATE_BLOCKS,
    TrajectoryReactionRateAccumulator,
    _correlated_batch_mean_statistics,
)
from electron_swarm.solvers.monte_carlo.orbit import (
    _dc_flight_quadrature,
    _magnetic_flight_energy_eV,
)
from product_helpers import base_product_config, write_config


def _constant_process(
    name: str,
    sigma_m2: float,
    *,
    process_type: ProcessType = ProcessType.EXCITATION,
    threshold_eV: float | None = 0.0,
    mass_amu: float | None = None,
) -> CrossSectionProcess:
    return CrossSectionProcess(
        species="Ar",
        process=name,
        process_type=process_type,
        threshold_eV=threshold_eV,
        mass_amu=mass_amu,
        energy_eV=np.array([0.0, 100.0]),
        cross_section_m2=np.array([sigma_m2, sigma_m2]),
    )


def test_cell_integrated_rate_convolution_resolves_an_internal_threshold() -> None:
    process = CrossSectionProcess(
        species="Ar",
        process="thresholded",
        process_type=ProcessType.EXCITATION,
        threshold_eV=0.75,
        energy_eV=np.array([0.75, 1.0]),
        cross_section_m2=np.array([2.0e-20, 2.0e-20]),
    )
    centers = np.array([0.5])
    widths = np.array([1.0])
    density = np.array([1.0])

    integrated = cell_integrated_rate_coefficient(
        centers,
        widths,
        density,
        process,
    )
    midpoint = float(process.sigma(centers)[0]) * math.sqrt(
        2.0 * centers[0] * EV_TO_J / ELECTRON_MASS_KG
    )

    assert midpoint == 0.0
    assert integrated > 0.0


def test_cell_integrated_rate_convolution_splits_near_edge_threshold() -> None:
    threshold_eV = 0.99
    sigma_m2 = 2.0e-20
    process = CrossSectionProcess(
        species="Ar",
        process="near-edge threshold",
        process_type=ProcessType.EXCITATION,
        threshold_eV=threshold_eV,
        energy_eV=np.array([0.0, 1.0]),
        cross_section_m2=np.array([sigma_m2, sigma_m2]),
    )

    integrated = cell_integrated_rate_coefficient(
        np.array([0.5]),
        np.array([1.0]),
        np.array([1.0]),
        process,
    )
    speed_factor = math.sqrt(2.0 * EV_TO_J / ELECTRON_MASS_KG)
    expected = (
        sigma_m2
        * speed_factor
        * (2.0 / 3.0)
        * (1.0 - threshold_eV**1.5)
    )

    assert integrated == pytest.approx(expected, rel=1.0e-12)


def test_cell_integrated_rate_convolution_splits_cross_section_knot() -> None:
    knot_eV = 0.99
    sigma_m2 = 3.0e-20
    process = CrossSectionProcess(
        species="Ar",
        process="near-edge interpolation knot",
        process_type=ProcessType.MOMENTUM,
        threshold_eV=None,
        energy_eV=np.array([0.0, knot_eV, 1.0]),
        cross_section_m2=np.array([0.0, 0.0, sigma_m2]),
    )

    integrated = cell_integrated_rate_coefficient(
        np.array([0.5]),
        np.array([1.0]),
        np.array([1.0]),
        process,
    )
    slope = sigma_m2 / (1.0 - knot_eV)
    speed_factor = math.sqrt(2.0 * EV_TO_J / ELECTRON_MASS_KG)
    expected = speed_factor * slope * (
        (2.0 / 5.0) * (1.0 - knot_eV**2.5)
        - knot_eV * (2.0 / 3.0) * (1.0 - knot_eV**1.5)
    )

    assert integrated == pytest.approx(expected, rel=1.0e-12)


def test_trajectory_rate_accumulator_integrates_sigma_v_and_marks_zero_events() -> None:
    energy_eV = 2.0
    sigma_m2 = 3.0e-20
    density_m3 = 2.0e22
    processes = (
        _constant_process("sampled", sigma_m2),
        _constant_process("postprocessed", 2.0 * sigma_m2),
    )
    accumulator = TrajectoryReactionRateAccumulator(
        processes=processes,
        fractions=(1.0, 1.0),
        event_sampled=(True, False),
        gas_number_density_m3=density_m3,
    )
    for step in range(20):
        for particle in range(2):
            accumulator.record_residence(
                block_index=accumulator.block_index_for_trial(
                    production_step=step,
                    production_steps=20,
                    particle_index=particle,
                    particles_this_round=2,
                ),
                energy_eV=energy_eV,
                weight=1.0,
                dt_s=0.25,
            )

    sampled, postprocessed = accumulator.estimates()
    speed = math.sqrt(2.0 * energy_eV * EV_TO_J / ELECTRON_MASS_KG)
    assert sampled.rate_coefficient_m3_s == pytest.approx(sigma_m2 * speed)
    assert sampled.block_count_requested == MIN_REACTION_RATE_BLOCKS
    assert sampled.block_count_used == MIN_REACTION_RATE_BLOCKS
    assert sampled.effective_sample_size == pytest.approx(20.0)
    assert sampled.standard_error_m3_s == pytest.approx(0.0)
    assert sampled.event_count == 0
    assert sampled.event_observation_status == "unobserved_zero_events_upper_bound_only"
    assert sampled.zero_event_rate_upper_95_m3_s is not None
    assert sampled.zero_event_rate_upper_95_m3_s > 0.0

    assert postprocessed.rate_coefficient_m3_s == pytest.approx(2.0 * sigma_m2 * speed)
    assert postprocessed.event_observation_status == "not_sampled_by_trajectory_model"
    assert postprocessed.zero_event_rate_upper_95_m3_s is None


def test_residence_pool_keeps_main_poisson_event_evidence() -> None:
    process = _constant_process("pooled", 2.0e-20)
    main = TrajectoryReactionRateAccumulator(
        processes=(process,),
        fractions=(1.0,),
        event_sampled=(True,),
        gas_number_density_m3=1.0e20,
    )
    tail = TrajectoryReactionRateAccumulator(
        processes=(process,),
        fractions=(1.0,),
        event_sampled=(True,),
        gas_number_density_m3=1.0e20,
        poisson_event_evidence=False,
    )
    for block in range(MIN_REACTION_RATE_BLOCKS):
        main.record_residence(
            block_index=block,
            energy_eV=1.0,
            weight=1.0,
            dt_s=1.0,
        )
        main.record_event(0, 1.0)
        tail.record_residence(
            block_index=block,
            energy_eV=4.0,
            weight=1.0,
            dt_s=1.0,
        )
        tail.record_event(0, 1.0)

    main.merge_residence_statistics(tail)
    estimate = main.estimates()[0]
    speed_1 = math.sqrt(2.0 * EV_TO_J / ELECTRON_MASS_KG)
    speed_4 = math.sqrt(8.0 * EV_TO_J / ELECTRON_MASS_KG)

    assert estimate.rate_coefficient_m3_s == pytest.approx(
        2.0e-20 * (speed_1 + speed_4) / 2.0
    )
    assert estimate.weighted_residence_time_s == pytest.approx(40.0)
    assert estimate.event_count == MIN_REACTION_RATE_BLOCKS
    assert estimate.event_observation_residence_time_s == pytest.approx(20.0)
    assert estimate.event_observation_status == "observed"


@pytest.mark.parametrize(
    ("process_type", "threshold_eV"),
    [
        (ProcessType.EXCITATION, 11.5),
        (ProcessType.IONIZATION, 15.8),
    ],
)
def test_direct_threshold_energy_loss_is_same_residence_integral_as_rate(
    process_type: ProcessType,
    threshold_eV: float,
) -> None:
    process = _constant_process(
        process_type.value,
        2.0e-20,
        process_type=process_type,
        threshold_eV=threshold_eV,
    )
    accumulator = TrajectoryReactionRateAccumulator(
        processes=(process,),
        fractions=(1.0,),
        event_sampled=(True,),
        gas_number_density_m3=1.0e21,
        angular_scattering_model="isotropic",
    )
    for block in range(MIN_REACTION_RATE_BLOCKS):
        accumulator.record_residence(
            block_index=block,
            energy_eV=20.0 + block,
            weight=1.0,
            dt_s=1.0,
        )

    estimate = accumulator.estimates()[0]
    loss = estimate.energy_loss

    assert loss.status == "direct_residence_integral"
    assert loss.model == "event_model_fixed_threshold_energy_change"
    assert loss.rate_coefficient_eV_m3_s == pytest.approx(
        estimate.rate_coefficient_m3_s * threshold_eV
    )
    assert loss.mean_energy_loss_per_collision_eV == pytest.approx(threshold_eV)
    assert loss.standard_error_eV_m3_s == pytest.approx(
        threshold_eV * estimate.standard_error_m3_s
    )
    assert loss.relative_standard_error == pytest.approx(
        estimate.relative_standard_error
    )
    metadata = estimate.metadata()["energy_loss"]
    assert metadata["estimator"] == ("trajectory_time_average_sigma_v_delta_energy")
    assert metadata["corresponding_rate_coefficient_m3_s"] == pytest.approx(
        estimate.rate_coefficient_m3_s
    )


def test_direct_elastic_energy_loss_uses_exact_thermal_event_estimator() -> None:
    energy_eV = 3.0
    mass_amu = 39.948
    process = _constant_process(
        "elastic",
        3.0e-20,
        process_type=ProcessType.ELASTIC,
        mass_amu=mass_amu,
    )
    accumulator = TrajectoryReactionRateAccumulator(
        processes=(process,),
        fractions=(1.0,),
        event_sampled=(True,),
        gas_number_density_m3=1.0e21,
        angular_scattering_model="isotropic",
    )
    for block in range(MIN_REACTION_RATE_BLOCKS):
        accumulator.record_residence(
            block_index=block,
            energy_eV=energy_eV,
            weight=1.0,
            dt_s=0.5,
        )
        accumulator.record_event(
            0,
            1.0,
            energy_loss_eV=2.0e-3,
            block_index=block,
        )

    estimate = accumulator.estimates()[0]
    assert estimate.energy_loss.status == "direct_event_estimator"
    assert estimate.energy_loss.model == (
        "maxwellian_target_exact_binary_collision_isotropic"
    )
    assert estimate.energy_loss.mean_energy_loss_per_collision_eV == pytest.approx(
        2.0e-3
    )
    assert estimate.energy_loss.rate_coefficient_eV_m3_s == pytest.approx(4.0e-24)
    metadata = estimate.metadata()["energy_loss"]
    assert metadata["estimator"] == (
        "trajectory_event_energy_change_per_target_density_residence_time"
    )
    assert metadata["neutral_thermal_motion_model"] == (
        "maxwellian_relative_speed_exact_binary_collision"
    )
    assert metadata["gas_temperature_terms_included"] is True


def test_direct_elastic_energy_loss_uses_maxent_p1_exact_thermal_events() -> None:
    energy_eV = 4.0
    mass_amu = 40.0
    total = _constant_process(
        "elastic total",
        4.0e-20,
        process_type=ProcessType.ELASTIC,
        mass_amu=mass_amu,
    )
    momentum = _constant_process(
        "momentum transfer",
        1.0e-20,
        process_type=ProcessType.MOMENTUM,
        mass_amu=mass_amu,
    )
    accumulator = TrajectoryReactionRateAccumulator(
        processes=(total, momentum),
        fractions=(1.0, 1.0),
        event_sampled=(True, False),
        gas_number_density_m3=1.0e21,
        angular_scattering_model="maxent_p1",
    )
    for block in range(MIN_REACTION_RATE_BLOCKS):
        accumulator.record_residence(
            block_index=block,
            energy_eV=energy_eV,
            weight=1.0,
            dt_s=1.0,
        )
        accumulator.record_event(
            0,
            1.0,
            energy_loss_eV=1.0e-3,
            block_index=block,
        )

    elastic, companion = accumulator.estimates()
    assert elastic.energy_loss.model == (
        "maxwellian_target_exact_binary_collision_maxent_p1"
    )
    assert elastic.energy_loss.mean_energy_loss_per_collision_eV == pytest.approx(
        1.0e-3
    )
    assert elastic.energy_loss.rate_coefficient_eV_m3_s == pytest.approx(1.0e-24)
    assert companion.energy_loss.status == (
        "unsupported_not_sampled_by_trajectory_model"
    )
    assert companion.energy_loss.rate_coefficient_eV_m3_s is None


def test_direct_elastic_energy_loss_fails_honestly_without_target_mass() -> None:
    process = _constant_process(
        "elastic",
        1.0e-20,
        process_type=ProcessType.ELASTIC,
        mass_amu=None,
    )
    accumulator = TrajectoryReactionRateAccumulator(
        processes=(process,),
        fractions=(1.0,),
        event_sampled=(True,),
        gas_number_density_m3=1.0e21,
        angular_scattering_model="isotropic",
    )
    accumulator.record_residence(
        block_index=0,
        energy_eV=2.0,
        weight=1.0,
        dt_s=1.0,
    )

    loss = accumulator.estimates()[0].energy_loss
    assert loss.status == "unsupported_missing_target_mass"
    assert loss.rate_coefficient_eV_m3_s is None
    assert loss.standard_error_eV_m3_s is None


def test_dc_flight_gauss_quadrature_integrates_energy_and_energy_flux() -> None:
    before = np.zeros(3, dtype=float)
    after = np.array([2.0e6, 0.0, 0.0], dtype=float)
    duration = 3.0e-9

    endpoint = 0.5 * ELECTRON_MASS_KG * float(np.dot(after, after)) / EV_TO_J
    quadrature = _dc_flight_quadrature(before, after, duration)
    integrated_energy = sum(energy * dt for energy, dt, _ in quadrature)
    integrated_energy_displacement = sum(
        energy * displacement for energy, _, displacement in quadrature
    )

    assert integrated_energy == pytest.approx(endpoint * duration / 3.0)
    assert integrated_energy_displacement == pytest.approx(
        endpoint * after * duration / 4.0
    )


def test_magnetic_flight_sample_does_not_create_spurious_energy_loss() -> None:
    speed = 2.0e6
    before = np.array([speed, 0.0, 0.0], dtype=float)
    after = np.array([0.0, speed, 0.0], dtype=float)

    midpoint = _magnetic_flight_energy_eV(
        before,
        after,
    )
    endpoint = 0.5 * ELECTRON_MASS_KG * speed * speed / EV_TO_J

    assert midpoint == pytest.approx(endpoint)


def test_zero_event_upper_bound_is_invariant_to_population_weight_scale() -> None:
    process = _constant_process("sampled", 1.0e-20)

    def estimate(weight: float):
        accumulator = TrajectoryReactionRateAccumulator(
            processes=(process,),
            fractions=(0.25,),
            event_sampled=(True,),
            gas_number_density_m3=2.0e22,
        )
        for block in range(MIN_REACTION_RATE_BLOCKS):
            accumulator.record_residence(
                block_index=block,
                energy_eV=2.0,
                weight=weight,
                dt_s=0.5,
            )
        return accumulator.estimates()[0]

    unit_weight = estimate(1.0)
    scaled_weight = estimate(10.0)

    assert scaled_weight.rate_coefficient_m3_s == pytest.approx(
        unit_weight.rate_coefficient_m3_s
    )
    assert scaled_weight.zero_event_rate_upper_95_m3_s == pytest.approx(
        unit_weight.zero_event_rate_upper_95_m3_s
    )
    assert scaled_weight.event_observation_residence_time_s == pytest.approx(
        unit_weight.event_observation_residence_time_s
    )
    assert scaled_weight.weighted_residence_time_s == pytest.approx(
        10.0 * unit_weight.weighted_residence_time_s
    )


@pytest.mark.parametrize("energy_eV", [0.75, 1.0, 1.5, 2.0, 2.5, 3.0])
def test_cached_cross_section_evaluation_preserves_process_endpoint_policies(
    energy_eV: float,
) -> None:
    processes = (
        CrossSectionProcess(
            species="Ar",
            process="zero",
            process_type=ProcessType.EXCITATION,
            energy_eV=np.array([1.0, 2.0]),
            cross_section_m2=np.array([2.0e-20, 4.0e-20]),
            metadata={"high_energy_extrapolation": "zero"},
        ),
        CrossSectionProcess(
            species="Ar",
            process="hold",
            process_type=ProcessType.EXCITATION,
            energy_eV=np.array([0.5, 3.0]),
            cross_section_m2=np.array([1.0e-20, 3.0e-20]),
            metadata={"high_energy_extrapolation": "hold"},
        ),
    )
    accumulator = TrajectoryReactionRateAccumulator(
        processes=processes,
        fractions=(1.0, 1.0),
        event_sampled=(True, True),
        gas_number_density_m3=1.0e20,
    )
    accumulator.record_residence(
        block_index=0,
        energy_eV=energy_eV,
        weight=1.0,
        dt_s=1.0,
    )

    estimates = accumulator.estimates()
    speed = math.sqrt(2.0 * energy_eV * EV_TO_J / ELECTRON_MASS_KG)
    for process, estimate in zip(processes, estimates, strict=True):
        expected = float(process.sigma(np.array([energy_eV]))[0]) * speed
        assert estimate.rate_coefficient_m3_s == pytest.approx(expected)


def test_cached_cross_section_evaluation_preserves_error_policy() -> None:
    process = CrossSectionProcess(
        species="Ar",
        process="strict",
        process_type=ProcessType.EXCITATION,
        energy_eV=np.array([0.0, 2.0]),
        cross_section_m2=np.array([1.0e-20, 2.0e-20]),
        metadata={"high_energy_extrapolation": "error"},
    )
    accumulator = TrajectoryReactionRateAccumulator(
        processes=(process,),
        fractions=(1.0,),
        event_sampled=(True,),
        gas_number_density_m3=1.0e20,
    )

    with pytest.raises(ValueError, match="above the tabulated range"):
        accumulator.record_residence(
            block_index=0,
            energy_eV=2.5,
            weight=1.0,
            dt_s=1.0,
        )


def test_batch_means_ess_accounts_for_positive_autocorrelation() -> None:
    correlated = np.linspace(1.0, 20.0, MIN_REACTION_RATE_BLOCKS)
    statistics = _correlated_batch_mean_statistics(
        correlated,
        reference_mean=float(np.mean(correlated)),
    )

    assert statistics.lag1_autocorrelation is not None
    assert statistics.lag1_autocorrelation > 0.0
    assert statistics.effective_sample_size is not None
    assert 1.0 <= statistics.effective_sample_size < MIN_REACTION_RATE_BLOCKS
    assert statistics.standard_error is not None
    assert statistics.standard_error > 0.0
    assert statistics.autocorrelation_time is not None
    assert statistics.autocorrelation_time > 1.0


def test_batch_means_ess_accounts_for_unequal_residence_exposure() -> None:
    values = np.tile(np.array([1.0, 2.0]), MIN_REACTION_RATE_BLOCKS // 2)
    equal = _correlated_batch_mean_statistics(
        values,
        reference_mean=float(np.mean(values)),
    )
    weights = np.ones(MIN_REACTION_RATE_BLOCKS, dtype=float)
    weights[0] = 20.0
    weighted_mean = float(np.average(values, weights=weights))
    unequal = _correlated_batch_mean_statistics(
        values,
        reference_mean=weighted_mean,
        weights=weights,
    )

    assert equal.effective_sample_size is not None
    assert unequal.effective_sample_size is not None
    assert unequal.effective_sample_size < equal.effective_sample_size
    assert unequal.standard_error is not None


@pytest.mark.mc
def test_internal_mc_uses_direct_rates_and_emits_sampling_metadata(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {
        "particles": 4,
        "warmup_collisions": 3,
        "max_collisions": 20,
        "seed": 123,
    }
    config = load_config(write_config(tmp_path, data))

    [case] = run(config, write=False, collect_diagnostics=True).cases

    diagnostics = case.diagnostics["internal_monte_carlo_reaction_rates"]
    assert diagnostics["estimator"] == "trajectory_time_average_sigma_v"
    rows = diagnostics["rates"]
    assert len(rows) == len(case.rates)
    for rate, row in zip(case.rates, rows, strict=True):
        assert row["estimator"] == "trajectory_time_average_sigma_v"
        assert row["rate_coefficient_m3_s"] == pytest.approx(rate.rate_coefficient_m3_s)
        assert row["batch_means_blocks_requested"] >= 20
        assert row["batch_means_blocks_used"] == 20
        assert row["batch_means_effective_sample_size"] is not None
        assert row["batch_means_standard_error_m3_s"] is not None
        loss = row["energy_loss"]
        if rate.process_type in {"elastic", "momentum", "effective"}:
            assert loss["status"] == "direct_event_estimator"
            assert loss["estimator"] == (
                "trajectory_event_energy_change_per_target_density_residence_time"
            )
            assert loss["neutral_thermal_motion_model"] == (
                "maxwellian_relative_speed_exact_binary_collision"
            )
            assert loss["gas_temperature_terms_included"] is True
        else:
            assert loss["status"] == "direct_residence_integral"
            assert loss["estimator"] == ("trajectory_time_average_sigma_v_delta_energy")
            assert loss["corresponding_rate_coefficient_m3_s"] == pytest.approx(
                rate.rate_coefficient_m3_s
            )
        assert loss["batch_means_standard_error_eV_m3_s"] is not None
        assert row["event_count"] >= 0
        if row["event_count"] == 0 and row["event_sampling_enabled"]:
            assert (
                row["event_observation_status"]
                == "unobserved_zero_events_upper_bound_only"
            )
            assert row["zero_event_rate_upper_95_m3_s"] is not None

    consistency = diagnostics["energy_audit_consistency"]
    assert consistency["status"] == "matched"
    assert consistency["direct_event_physical_energy_loss_eV"] > 0.0
    assert consistency["direct_event_physical_energy_loss_eV"] == pytest.approx(
        consistency["mc_energy_audit_physical_collision_loss_eV"],
        rel=1.0e-12,
        abs=1.0e-12,
    )
