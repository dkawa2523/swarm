"""Limited internal Monte Carlo backend with Boris magnetic push."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.constants import (
    AMU_KG,
    BOLTZMANN_J_K,
    E_CHARGE_C,
    ELECTRON_MASS_KG,
    EV_TO_J,
    TOWNSEND,
)
from electron_swarm.core.cross_sections import (
    CrossSectionProcess,
    CrossSectionSet,
    ProcessType,
    gas_mass_amu,
    mixture_fraction,
)
from electron_swarm.core.results import RateResult, SwarmCaseResult
from electron_swarm.physics.angular_scattering import (
    expected_angular_metadata,
    build_angular_model,
)


SCATTERING_TYPES = {
    ProcessType.MOMENTUM,
    ProcessType.ELASTIC,
    ProcessType.EFFECTIVE,
}
TRAJECTORY_REACTION_TYPES = {
    ProcessType.EXCITATION,
    ProcessType.IONIZATION,
    ProcessType.SUPERELASTIC,
}


def gas_number_density(config: SwarmConfig) -> float:
    cond = config.conditions
    if cond.gas_number_density_m3 is not None:
        return float(cond.gas_number_density_m3)
    assert cond.pressure_Pa is not None
    return float(cond.pressure_Pa / (BOLTZMANN_J_K * cond.gas_temperature_K))


def magnetic_field_vector(B_T: float, angle_EB_deg: float) -> np.ndarray:
    if B_T < 0.0 or not np.isfinite(B_T) or not np.isfinite(angle_EB_deg):
        raise ValueError("magnetic field magnitude and angle must be finite")
    if angle_EB_deg < 0.0 or angle_EB_deg > 180.0:
        raise ValueError("magnetic field angle_EB_deg must be in [0, 180]")
    angle = np.deg2rad(float(angle_EB_deg))
    return float(B_T) * np.array([np.sin(angle), 0.0, np.cos(angle)], dtype=float)


def boris_push(
    velocity_m_s: np.ndarray,
    electric_field_V_m: np.ndarray,
    magnetic_field_T: np.ndarray,
    dt_s: float,
) -> np.ndarray:
    """Advance electron velocity with the Boris algorithm."""

    v = np.asarray(velocity_m_s, dtype=float)
    E = np.asarray(electric_field_V_m, dtype=float)
    B = np.asarray(magnetic_field_T, dtype=float)
    qmdt2 = -E_CHARGE_C / ELECTRON_MASS_KG * float(dt_s) * 0.5
    v_minus = v + qmdt2 * E
    t = qmdt2 * B
    s = 2.0 * t / (1.0 + float(np.dot(t, t)))
    v_prime = v_minus + np.cross(v_minus, t)
    v_plus = v_minus + np.cross(v_prime, s)
    return v_plus + qmdt2 * E


@dataclass(frozen=True, slots=True)
class _ProjectedProcess:
    process: CrossSectionProcess
    fraction: float


def _project_processes(
    config: SwarmConfig, cross_sections: CrossSectionSet
) -> list[_ProjectedProcess]:
    out: list[_ProjectedProcess] = []
    for proc in cross_sections.processes:
        fraction = mixture_fraction(config.conditions, proc.species)
        if fraction > 0.0:
            out.append(_ProjectedProcess(proc, float(fraction)))
    if not out:
        raise ValueError("internal monte_carlo requires at least one active process")
    return out


def _energy_from_velocity(velocity: np.ndarray) -> float:
    return float(0.5 * ELECTRON_MASS_KG * np.dot(velocity, velocity) / EV_TO_J)


def _speed_from_energy(energy_eV: float) -> float:
    return float(np.sqrt(max(2.0 * energy_eV * EV_TO_J / ELECTRON_MASS_KG, 0.0)))


def _random_direction(rng: np.random.Generator) -> np.ndarray:
    mu = rng.uniform(-1.0, 1.0)
    phi = rng.uniform(0.0, 2.0 * np.pi)
    sint = np.sqrt(max(1.0 - mu * mu, 0.0))
    return np.array([sint * np.cos(phi), sint * np.sin(phi), mu], dtype=float)


def _sigma_arrays(
    processes: list[_ProjectedProcess], energy_eV: float
) -> tuple[np.ndarray, float, float, float]:
    raw_values = np.array(
        [
            item.fraction * float(item.process.sigma(np.array([energy_eV]))[0])
            for item in processes
        ],
        dtype=float,
    )
    raw_values = np.clip(
        np.nan_to_num(raw_values, nan=0.0, posinf=0.0, neginf=0.0), 0.0, None
    )
    trajectory_values = np.array(
        [
            value
            if item.process.process_type in SCATTERING_TYPES
            or item.process.process_type in TRAJECTORY_REACTION_TYPES
            else 0.0
            for value, item in zip(raw_values, processes, strict=True)
        ],
        dtype=float,
    )
    collision_total = float(np.sum(trajectory_values))
    scattering_total = float(
        np.sum(
            [
                value
                for value, item in zip(raw_values, processes, strict=True)
                if item.process.process_type in SCATTERING_TYPES
            ]
        )
    )
    sigma_momentum = scattering_total
    return trajectory_values, collision_total, scattering_total, sigma_momentum


def _trial_collision_frequency(
    processes: list[_ProjectedProcess], density_m3: float
) -> float:
    grids = [item.process.energy_eV for item in processes]
    energy = np.unique(np.concatenate(grids))
    energy = energy[np.isfinite(energy) & (energy >= 0.0)]
    if energy.size == 0:
        raise ValueError("internal monte_carlo requires finite cross-section energies")
    speed = np.maximum(
        np.sqrt(np.maximum(2.0 * energy * EV_TO_J / ELECTRON_MASS_KG, 0.0)),
        1.0,
    )
    total_sigma = np.zeros_like(energy, dtype=float)
    for item in processes:
        total_sigma += item.fraction * item.process.sigma(
            energy, right=float(item.process.cross_section_m2[-1])
        )
    nu = float(np.max(density_m3 * speed * total_sigma))
    if not np.isfinite(nu) or nu <= 0.0:
        raise ValueError("internal monte_carlo could not build a positive trial collision frequency")
    return 1.2 * nu


def _scatter_velocity(
    velocity: np.ndarray,
    mu: float,
    rng: np.random.Generator,
    energy_eV: float | None = None,
) -> np.ndarray:
    speed = (
        _speed_from_energy(float(energy_eV))
        if energy_eV is not None
        else float(np.linalg.norm(velocity))
    )
    if speed <= 0.0:
        return _random_direction(rng) * _speed_from_energy(1.0e-3)
    axis = velocity / speed
    ref = np.array([1.0, 0.0, 0.0]) if abs(axis[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = np.cross(axis, ref)
    e1 /= max(float(np.linalg.norm(e1)), 1.0e-300)
    e2 = np.cross(axis, e1)
    phi = rng.uniform(0.0, 2.0 * np.pi)
    sin_theta = np.sqrt(max(1.0 - mu * mu, 0.0))
    new_dir = mu * axis + sin_theta * (np.cos(phi) * e1 + np.sin(phi) * e2)
    return speed * new_dir


def _elastic_energy_after_collision(
    config: SwarmConfig,
    process: CrossSectionProcess,
    energy_eV: float,
    mu: float,
) -> float:
    mass_amu = process.mass_amu or gas_mass_amu(config.conditions, process.species)
    mass_ratio = 2.0 * ELECTRON_MASS_KG / max(mass_amu * AMU_KG, 1.0e-300)
    loss = float(energy_eV) * mass_ratio * max(1.0 - float(mu), 0.0)
    return max(float(energy_eV) - loss, 1.0e-4)


def _energy_loss_eV(process_type: ProcessType, threshold_eV: float | None) -> float:
    if process_type in {ProcessType.EXCITATION, ProcessType.IONIZATION}:
        return float(threshold_eV or 0.0)
    if process_type == ProcessType.SUPERELASTIC:
        return -abs(float(threshold_eV or 0.0))
    return 0.0


def _rate_results(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    energy: np.ndarray,
    eedf: np.ndarray,
    widths: np.ndarray,
    case_id: str,
    e_over_n_Td: float,
) -> tuple[list[RateResult], float, float]:
    speed = np.sqrt(np.maximum(2.0 * EV_TO_J * energy / ELECTRON_MASS_KG, 0.0))
    density = gas_number_density(config)
    rates: list[RateResult] = []
    ion_rate = 0.0
    attach_rate = 0.0
    for proc in cross_sections.processes:
        frac = mixture_fraction(config.conditions, proc.species)
        if frac <= 0.0:
            continue
        k = float(np.sum(proc.sigma(energy) * speed * eedf * widths))
        kmix = frac * k
        freq = density * kmix
        if proc.process_type == ProcessType.IONIZATION:
            ion_rate += kmix
        elif proc.process_type == ProcessType.ATTACHMENT:
            attach_rate += kmix
        rates.append(
            RateResult(
                solver="monte_carlo",
                case_id=case_id,
                e_over_n_Td=e_over_n_Td,
                species=proc.species,
                process=proc.process,
                process_type=proc.process_type.value,
                threshold_eV=proc.threshold_eV,
                rate_coefficient_m3_s=k,
                mixture_weighted_rate_m3_s=kmix,
                frequency_s_inv=freq,
                power_loss_eV_s=freq
                * _energy_loss_eV(proc.process_type, proc.threshold_eV),
            )
        )
    net_freq = density * (ion_rate - attach_rate)
    return rates, net_freq, ion_rate - attach_rate


def run_internal_monte_carlo(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
) -> list[SwarmCaseResult]:
    cfg = config.solvers.monte_carlo
    rng = np.random.default_rng(cfg.seed)
    particles = int(cfg.particles or 256)
    collisions = int(cfg.max_collisions or 100)
    density = gas_number_density(config)
    projected = _project_processes(config, cross_sections)
    trial_frequency = _trial_collision_frequency(projected, density)
    angular_model = build_angular_model(config)
    angular_name = config.physics.angular_scattering.model
    if angular_name not in {"isotropic", "maxent_p1"}:
        raise NotImplementedError(
            f"internal monte_carlo has no product sampler for angular model {angular_name!r}"
        )

    magnetic = config.physics.field.magnetic_field
    B = magnetic_field_vector(magnetic.B_T, magnetic.angle_EB_deg)
    magnetic_enabled = bool(magnetic.enabled and magnetic.B_T > 0.0)
    out: list[SwarmCaseResult] = []
    for i, e_over_n_Td in enumerate(config.run.e_over_n_Td):
        case_id = f"{config.run.case_prefix}_{i:04d}"
        E_scalar = float(e_over_n_Td) * TOWNSEND * density
        E = np.array([0.0, 0.0, E_scalar], dtype=float)
        positions = np.zeros((particles, 3), dtype=float)
        times = np.zeros(particles, dtype=float)
        energy_samples: list[float] = []
        initial_speed = _speed_from_energy(1.0)
        velocities = np.array(
            [_random_direction(rng) * initial_speed for _ in range(particles)],
            dtype=float,
        )
        for _ in range(collisions):
            for p in range(particles):
                dt = float(rng.exponential(1.0 / trial_frequency))
                if magnetic_enabled:
                    gyro = E_CHARGE_C * magnetic.B_T / ELECTRON_MASS_KG
                    dt = min(dt, 0.2 / max(gyro, 1.0))
                velocities[p] = boris_push(velocities[p], E, B, dt)
                positions[p] += velocities[p] * dt
                times[p] += dt
                energy_eV = min(max(_energy_from_velocity(velocities[p]), 1.0e-6), 1.0e4)
                energy_samples.append(energy_eV)
                speed = max(float(np.linalg.norm(velocities[p])), 1.0)
                sigmas, collision_total, scattering_total, sigma_momentum = _sigma_arrays(
                    projected, energy_eV
                )
                if collision_total <= 0.0:
                    continue
                collision_frequency = density * collision_total * speed
                if rng.random() > min(collision_frequency / trial_frequency, 1.0):
                    continue
                choice = int(rng.choice(len(projected), p=sigmas / collision_total))
                proc = projected[choice].process
                if proc.process_type in SCATTERING_TYPES:
                    if angular_name == "maxent_p1":
                        total = max(scattering_total, 1.0e-300)
                        momentum = max(sigma_momentum, 1.0e-300)
                        mu = float(
                            angular_model.sample_mu(
                                energy_eV,
                                rng,
                                sigma_total=total,
                                sigma_momentum=momentum,
                            )
                        )
                    else:
                        mu = float(angular_model.sample_mu(energy_eV, rng))
                    post_energy = _elastic_energy_after_collision(
                        config, proc, energy_eV, mu
                    )
                    velocities[p] = _scatter_velocity(
                        velocities[p], mu, rng, post_energy
                    )
                elif proc.process_type in TRAJECTORY_REACTION_TYPES:
                    new_energy = max(
                        energy_eV - _energy_loss_eV(proc.process_type, proc.threshold_eV),
                        1.0e-4,
                    )
                    velocities[p] = _random_direction(rng) * _speed_from_energy(new_energy)

        samples = np.asarray(energy_samples, dtype=float)
        max_energy = max(20.0, float(np.percentile(samples, 99.5)) * 1.2)
        edges = np.linspace(0.0, max_energy, 81)
        counts, edges = np.histogram(samples, bins=edges)
        widths = np.diff(edges)
        energy = 0.5 * (edges[:-1] + edges[1:])
        eedf = counts.astype(float) / max(float(np.sum(counts)), 1.0) / widths
        if not np.any(eedf > 0.0):
            eedf[0] = 1.0 / widths[0]
        eedf = eedf / max(float(np.sum(eedf * widths)), 1.0e-300)
        mean_energy = float(np.sum(energy * eedf * widths))
        mean_time = max(float(np.mean(times)), 1.0e-300)
        drift_velocity = -float(np.mean(positions[:, 2] / np.maximum(times, 1.0e-300)))
        mobility = drift_velocity / max(E_scalar, 1.0e-300)
        dz = -positions[:, 2] - drift_velocity * times
        diffusion_l = float(np.var(dz) / (2.0 * mean_time))
        transverse = positions[:, :2]
        diffusion_t = float(np.mean(np.sum(transverse * transverse, axis=1)) / (4.0 * mean_time))
        rates, net_freq, net_rate = _rate_results(
            config, cross_sections, energy, eedf, widths, case_id, float(e_over_n_Td)
        )
        effective_townsend = net_freq / max(abs(drift_velocity) * density, 1.0e-300)
        metadata = {
            "adapter": "internal_monte_carlo",
            "monte_carlo_backend": "internal",
            "particles": particles,
            "max_collisions": collisions,
            "field_integrator": "boris" if magnetic.enabled else "none",
            "magnetic_field_treatment": (
                "boris_lorentz_push" if magnetic.enabled else "none"
            ),
            "magnetic_field_B_T": float(magnetic.B_T),
            "magnetic_field_angle_EB_deg": float(magnetic.angle_EB_deg),
            "magnetic_field_Bx_T": float(B[0]),
            "magnetic_field_By_T": float(B[1]),
            "magnetic_field_Bz_T": float(B[2]),
            "gas_number_density_m-3": float(density),
            "electric_field_V_m": float(E_scalar),
            "convolution_effective_rate_coefficient_m3_s": float(net_rate),
            "reaction_rates_source": "eedf_convolution",
            "trajectory_reaction_model": "energy_loss_no_branching",
            "attachment_trajectory_treatment": "rate_convolution_only",
            **expected_angular_metadata(config),
        }
        out.append(
            SwarmCaseResult(
                solver="monte_carlo",
                case_id=case_id,
                e_over_n_Td=float(e_over_n_Td),
                mean_energy_eV=mean_energy,
                drift_velocity_m_s=drift_velocity,
                mobility_m2_V_s=mobility,
                reduced_mobility_m2_V_s_m3=mobility * density,
                diffusion_L_m2_s=diffusion_l,
                diffusion_T_m2_s=diffusion_t,
                reduced_diffusion_L_m2_s_m3=diffusion_l * density,
                reduced_diffusion_T_m2_s_m3=diffusion_t * density,
                net_ionization_frequency_s=net_freq,
                effective_townsend_m2=effective_townsend,
                energy_eV=energy,
                eedf=eedf,
                eepf=eedf / np.sqrt(np.maximum(energy, 1.0e-30)),
                rates=rates,
                metadata=metadata,
            )
        )
    return out
