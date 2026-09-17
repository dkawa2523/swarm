from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from electron_swarm import load_config
from electron_swarm.core.constants import (
    AMU_KG,
    BOLTZMANN_J_K,
    ELECTRON_MASS_KG,
    EV_TO_J,
)
from electron_swarm.core.cross_sections import CrossSectionProcess, ProcessType
from electron_swarm.physics.electron_neutral import (
    elastic_binary_collision_velocity_m_s,
    maxwell_velocity_component_std_m_s,
    speed_from_energy_m_s,
)
from electron_swarm.solvers.monte_carlo.collisions import (
    PreparedCollisionSampler,
    _ProjectedProcess,
    _maximum_speed_weighted_cross_section,
)

from product_helpers import base_product_config, write_config


def _thermal_config(tmp_path: Path, *, mass_amu: float = 39.948):
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["conditions"]["gas_mixture"][0]["mass_amu"] = mass_amu
    data["solvers"]["monte_carlo"] = {"numeric_kernel": "python"}
    return load_config(write_config(tmp_path, data, "thermal_mc.yaml"))


def _constant_elastic_process(
    mass_amu: float,
    *,
    extrapolation: str = "hold",
) -> CrossSectionProcess:
    return CrossSectionProcess(
        species="Ar",
        process="elastic",
        process_type=ProcessType.ELASTIC,
        energy_eV=np.asarray([0.0, 1.0], dtype=float),
        cross_section_m2=np.asarray([2.0e-20, 2.0e-20], dtype=float),
        mass_amu=mass_amu,
        metadata={"high_energy_extrapolation": extrapolation},
    )


def _mean_relative_speed(speed_m_s: float, component_std_m_s: float) -> float:
    speed = float(speed_m_s)
    std = float(component_std_m_s)
    if speed == 0.0:
        return 2.0 * std * math.sqrt(2.0 / math.pi)
    scaled = speed / (math.sqrt(2.0) * std)
    return (
        std * math.sqrt(2.0 / math.pi) * math.exp(-(scaled * scaled))
        + (speed + std * std / speed) * math.erf(scaled)
    )


def test_exact_binary_elastic_collision_conserves_pair_invariants() -> None:
    mass_amu = 39.948
    target_mass = mass_amu * AMU_KG
    electron = np.asarray([1.8e5, -2.5e5, 7.0e4], dtype=float)
    target = np.asarray([310.0, -120.0, 80.0], dtype=float)

    post_electron = elastic_binary_collision_velocity_m_s(
        electron,
        target,
        scattering_cosine=-0.37,
        azimuth_rad=1.2,
        target_mass_amu=mass_amu,
    )
    momentum = ELECTRON_MASS_KG * electron + target_mass * target
    post_target = (momentum - ELECTRON_MASS_KG * post_electron) / target_mass

    before_energy = 0.5 * (
        ELECTRON_MASS_KG * float(electron @ electron)
        + target_mass * float(target @ target)
    )
    after_energy = 0.5 * (
        ELECTRON_MASS_KG * float(post_electron @ post_electron)
        + target_mass * float(post_target @ post_target)
    )
    assert ELECTRON_MASS_KG * post_electron + target_mass * post_target == pytest.approx(
        momentum,
        rel=2.0e-14,
        abs=1.0e-40,
    )
    assert after_energy == pytest.approx(before_energy, rel=2.0e-14)
    assert np.linalg.norm(post_electron - post_target) == pytest.approx(
        np.linalg.norm(electron - target),
        rel=2.0e-14,
    )


def test_piecewise_linear_speed_cross_section_majorant_finds_interior_peak() -> None:
    energy = np.asarray([0.0, 4.0], dtype=float)
    sigma = np.asarray([4.0, 0.0], dtype=float)
    critical_energy = 4.0 / 3.0
    expected = math.sqrt(
        2.0 * EV_TO_J * critical_energy / ELECTRON_MASS_KG
    ) * (4.0 - critical_energy)

    assert _maximum_speed_weighted_cross_section(energy, sigma) == pytest.approx(
        expected,
        rel=2.0e-15,
    )


def test_hold_tail_marked_sampler_recovers_maxwell_relative_collision_rate(
    tmp_path: Path,
) -> None:
    mass_amu = 0.05
    config = _thermal_config(tmp_path, mass_amu=mass_amu)
    process = _constant_elastic_process(mass_amu)
    sampler = PreparedCollisionSampler.build(
        config,
        [_ProjectedProcess(process, 1.0)],
        angular_model_name="isotropic",
        max_energy_eV=0.2,
    )
    electron_speed = float(speed_from_energy_m_s(0.1))
    electron_velocity = np.asarray([0.0, 0.0, electron_speed], dtype=float)
    std = maxwell_velocity_component_std_m_s(300.0, mass_amu)
    expected = 2.0e-20 * _mean_relative_speed(electron_speed, std)

    rng = np.random.default_rng(9417)
    accepted = 0
    maximum_probability = 0.0
    trials = 50_000
    for _ in range(trials):
        trial = sampler.sample_trial(electron_velocity, rng)
        maximum_probability = max(
            maximum_probability,
            trial.acceptance_probability,
        )
        accepted += trial.process_index is not None
    estimated = sampler.trial_rate_coefficient_m3_s * accepted / trials

    assert maximum_probability <= 1.0
    assert estimated == pytest.approx(expected, rel=0.012)


def test_marked_exact_collision_chain_relaxes_to_gas_temperature(
    tmp_path: Path,
) -> None:
    # A light manufactured target accelerates relaxation without changing the
    # invariant Maxwell distribution checked by this test.
    mass_amu = 0.01
    temperature_K = 300.0
    config = _thermal_config(tmp_path, mass_amu=mass_amu)
    process = _constant_elastic_process(mass_amu)
    sampler = PreparedCollisionSampler.build(
        config,
        [_ProjectedProcess(process, 1.0)],
        angular_model_name="isotropic",
        max_energy_eV=0.5,
    )
    velocity = np.asarray(
        [0.0, 0.0, float(speed_from_energy_m_s(0.35))],
        dtype=float,
    )
    rng = np.random.default_rng(20260907)
    samples: list[float] = []
    for proposal in range(90_000):
        trial = sampler.sample_trial(velocity, rng)
        if trial.process_index is not None:
            if trial.target_velocity_m_s is None:
                raise AssertionError("accepted thermal collision has no target")
            velocity = elastic_binary_collision_velocity_m_s(
                velocity,
                trial.target_velocity_m_s,
                scattering_cosine=float(rng.uniform(-1.0, 1.0)),
                azimuth_rad=float(rng.uniform(0.0, 2.0 * math.pi)),
                target_mass_amu=mass_amu,
            )
        if proposal >= 10_000:
            samples.append(
                0.5 * ELECTRON_MASS_KG * float(velocity @ velocity) / EV_TO_J
            )

    expected_mean_eV = 1.5 * BOLTZMANN_J_K * temperature_K / EV_TO_J
    assert float(np.mean(samples)) == pytest.approx(expected_mean_eV, rel=0.04)
