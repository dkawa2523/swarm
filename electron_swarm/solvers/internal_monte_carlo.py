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


@dataclass(frozen=True, slots=True)
class _PostReactionOutcome:
    tracked_energy_eV: float
    inelastic_energy_loss_eV: float = 0.0
    ionization_threshold_loss_eV: float = 0.0
    ionization_untracked_secondary_energy_eV: float = 0.0


@dataclass(slots=True)
class _EnergyAudit:
    field_work_eV: float = 0.0
    elastic_energy_loss_eV: float = 0.0
    inelastic_energy_loss_eV: float = 0.0
    ionization_threshold_loss_eV: float = 0.0
    ionization_untracked_secondary_energy_eV: float = 0.0
    tracked_particle_initial_energy_eV: float = 0.0
    tracked_particle_final_energy_eV: float = 0.0
    null_collision_trials: int = 0
    accepted_collisions: int = 0
    max_collision_to_trial_ratio: float = 0.0

    def record_field_push(
        self, before_eV: float, after_eV: float, weight: float = 1.0
    ) -> None:
        self.field_work_eV += float(weight) * (float(after_eV) - float(before_eV))

    def record_trial_ratio(self, ratio: float) -> None:
        self.null_collision_trials += 1
        if np.isfinite(ratio):
            self.max_collision_to_trial_ratio = max(
                self.max_collision_to_trial_ratio, float(ratio)
            )

    def record_elastic_collision(
        self, before_eV: float, after_eV: float, weight: float = 1.0
    ) -> None:
        self.accepted_collisions += 1
        self.elastic_energy_loss_eV += float(weight) * (
            float(before_eV) - float(after_eV)
        )

    def record_reaction(
        self, outcome: _PostReactionOutcome, weight: float = 1.0
    ) -> None:
        self.accepted_collisions += 1
        self.inelastic_energy_loss_eV += float(weight) * float(
            outcome.inelastic_energy_loss_eV
        )
        self.ionization_threshold_loss_eV += float(weight) * float(
            outcome.ionization_threshold_loss_eV
        )
        self.ionization_untracked_secondary_energy_eV += float(weight) * float(
            outcome.ionization_untracked_secondary_energy_eV
        )

    def as_metadata(self) -> dict[str, float | str]:
        tracked_delta = (
            self.tracked_particle_final_energy_eV
            - self.tracked_particle_initial_energy_eV
        )
        expected_delta = (
            self.field_work_eV
            - self.elastic_energy_loss_eV
            - self.inelastic_energy_loss_eV
            - self.ionization_threshold_loss_eV
            - self.ionization_untracked_secondary_energy_eV
        )
        residual = tracked_delta - expected_delta
        scale = max(
            abs(tracked_delta),
            abs(self.field_work_eV)
            + abs(self.elastic_energy_loss_eV)
            + abs(self.inelastic_energy_loss_eV)
            + abs(self.ionization_threshold_loss_eV)
            + abs(self.ionization_untracked_secondary_energy_eV),
            1.0e-300,
        )
        residual_fraction = abs(residual) / scale
        if residual_fraction < 1.0e-3:
            status = "ok"
        elif residual_fraction < 1.0e-2:
            status = "warning"
        else:
            status = "fail"
        acceptance = (
            self.accepted_collisions / self.null_collision_trials
            if self.null_collision_trials
            else 0.0
        )
        return {
            "mc_field_work_eV": float(self.field_work_eV),
            "mc_elastic_energy_loss_eV": float(self.elastic_energy_loss_eV),
            "mc_inelastic_energy_loss_eV": float(self.inelastic_energy_loss_eV),
            "mc_ionization_threshold_loss_eV": float(
                self.ionization_threshold_loss_eV
            ),
            "mc_ionization_untracked_secondary_energy_eV": float(
                self.ionization_untracked_secondary_energy_eV
            ),
            "mc_tracked_particle_energy_delta_eV": float(tracked_delta),
            "mc_tracked_energy_balance_residual_eV": float(residual),
            "mc_tracked_energy_balance_residual_fraction": float(residual_fraction),
            "mc_physical_branching_gap_eV": float(
                self.ionization_untracked_secondary_energy_eV
            ),
            "mc_energy_balance_status": status,
            "mc_null_collision_acceptance_fraction": float(acceptance),
            "mc_max_collision_to_trial_ratio": float(
                self.max_collision_to_trial_ratio
            ),
        }


@dataclass(slots=True)
class _ParticleEnsemble:
    positions: np.ndarray
    velocities: np.ndarray
    times: np.ndarray
    weights: np.ndarray

    @classmethod
    def initialize(
        cls,
        particles: int,
        initial_speed: float,
        rng: np.random.Generator,
    ) -> "_ParticleEnsemble":
        return cls(
            positions=np.zeros((particles, 3), dtype=float),
            velocities=np.array(
                [_random_direction(rng) * initial_speed for _ in range(particles)],
                dtype=float,
            ),
            times=np.zeros(particles, dtype=float),
            weights=np.ones(particles, dtype=float),
        )

    def __len__(self) -> int:
        return int(self.weights.size)

    def total_weighted_energy_eV(self) -> float:
        return float(
            np.sum(
                self.weights
                * np.array([_energy_from_velocity(v) for v in self.velocities])
            )
        )

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
    processes: list[_ProjectedProcess],
    density_m3: float,
    *,
    max_energy_eV: float | None = None,
) -> float:
    grids = [item.process.energy_eV for item in processes]
    energy = np.unique(np.concatenate(grids))
    energy = energy[np.isfinite(energy) & (energy >= 0.0)]
    if energy.size == 0:
        raise ValueError("internal monte_carlo requires finite cross-section energies")
    if max_energy_eV is not None:
        max_energy = float(max_energy_eV)
        if not np.isfinite(max_energy) or max_energy <= 0.0:
            raise ValueError("physics.energy_grid_policy.max_eV_limit must be positive")
        if max_energy > float(energy[-1]):
            extension = np.linspace(float(energy[-1]), max_energy, 128)
            energy = np.unique(np.concatenate([energy, extension]))
    speed = np.maximum(
        np.sqrt(np.maximum(2.0 * energy * EV_TO_J / ELECTRON_MASS_KG, 0.0)),
        1.0,
    )
    total_sigma = np.zeros_like(energy, dtype=float)
    for item in processes:
        total_sigma += item.fraction * item.process.sigma(energy)
    nu = float(np.max(density_m3 * speed * total_sigma))
    if not np.isfinite(nu) or nu <= 0.0:
        raise ValueError("internal monte_carlo could not build a positive trial collision frequency")
    return 1.2 * nu


def _validate_trial_collision_frequency(
    collision_frequency_s_inv: float, trial_frequency_s_inv: float
) -> None:
    if collision_frequency_s_inv > trial_frequency_s_inv * (1.0 + 1.0e-12):
        raise RuntimeError(
            "internal_monte_carlo null-collision majorant was exceeded; "
            "increase physics.energy_grid_policy.max_eV_limit or revise cross sections"
        )


def _max_cross_section_energy(processes: list[_ProjectedProcess]) -> float:
    max_energy = max(float(np.nanmax(item.process.energy_eV)) for item in processes)
    if not np.isfinite(max_energy) or max_energy <= 0.0:
        raise ValueError("internal monte_carlo requires finite cross-section energies")
    return max_energy


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


def _post_reaction_energy(
    config: SwarmConfig,
    process_type: ProcessType,
    threshold_eV: float | None,
    energy_eV: float,
    rng: np.random.Generator,
) -> float:
    return _post_reaction_outcome(
        config, process_type, threshold_eV, energy_eV, rng
    ).tracked_energy_eV


def _post_reaction_outcome(
    config: SwarmConfig,
    process_type: ProcessType,
    threshold_eV: float | None,
    energy_eV: float,
    rng: np.random.Generator,
) -> _PostReactionOutcome:
    loss = _energy_loss_eV(process_type, threshold_eV)
    if process_type != ProcessType.IONIZATION:
        return _PostReactionOutcome(
            tracked_energy_eV=max(float(energy_eV) - loss, 1.0e-4),
            inelastic_energy_loss_eV=float(loss),
        )

    available = max(float(energy_eV) - float(threshold_eV or 0.0), 0.0)
    threshold_loss = float(threshold_eV or 0.0)
    sharing = config.physics.ionization.energy_sharing
    if sharing == "equal":
        tracked = 0.5 * available
        untracked = 0.5 * available
        return _PostReactionOutcome(
            tracked_energy_eV=max(tracked, 1.0e-4),
            ionization_threshold_loss_eV=threshold_loss,
            ionization_untracked_secondary_energy_eV=untracked,
        )
    if sharing == "primary_secondary":
        secondary = min(
            max(float(config.physics.ionization.secondary_electron_energy_eV), 0.0),
            available,
        )
        primary = max(available - secondary, 0.0)
        if rng.random() < 0.5:
            tracked, untracked = primary, secondary
        else:
            tracked, untracked = secondary, primary
        return _PostReactionOutcome(
            tracked_energy_eV=max(tracked, 1.0e-4),
            ionization_threshold_loss_eV=threshold_loss,
            ionization_untracked_secondary_energy_eV=untracked,
        )
    if sharing == "loss_only":
        return _PostReactionOutcome(
            tracked_energy_eV=max(available, 1.0e-4),
            ionization_threshold_loss_eV=threshold_loss,
        )
    raise ValueError(f"unsupported ionization energy sharing model: {sharing!r}")


def _ionization_branching_model(config: SwarmConfig) -> str:
    return (
        "none"
        if config.physics.ionization.energy_sharing == "loss_only"
        else "single_daughter_sampling"
    )


def _mc_energy_edges(max_energy_eV: float) -> np.ndarray:
    max_energy = float(max_energy_eV)
    if not np.isfinite(max_energy) or max_energy <= 0.0:
        raise ValueError("physics.energy_grid_policy.max_eV_limit must be positive")
    if max_energy <= 80.0:
        return np.linspace(0.0, max_energy, 81)
    parts = [np.linspace(0.0, 80.0, 161)]
    if max_energy > 80.0:
        parts.append(np.arange(90.0, min(max_energy, 260.0) + 1.0e-9, 10.0))
    if max_energy > 260.0:
        parts.append(np.arange(280.0, max_energy + 1.0e-9, 20.0))
    edges = np.unique(np.concatenate(parts))
    if edges[-1] < max_energy:
        edges = np.append(edges, max_energy)
    return edges


def _mc_bin_relative_standard_error(counts: np.ndarray) -> np.ndarray:
    values = np.asarray(counts, dtype=float)
    out = np.full_like(values, np.nan, dtype=float)
    mask = values > 0.0
    out[mask] = 1.0 / np.sqrt(values[mask])
    return out


def _mc_effective_bin_counts(
    weighted_sum: np.ndarray,
    weighted_square_sum: np.ndarray,
) -> np.ndarray:
    values = np.asarray(weighted_sum, dtype=float)
    squares = np.asarray(weighted_square_sum, dtype=float)
    out = np.zeros_like(values, dtype=float)
    mask = squares > 0.0
    out[mask] = values[mask] * values[mask] / squares[mask]
    return out


def _mc_tail_uncertainty_metadata(
    energy_eV: np.ndarray,
    counts: np.ndarray,
    threshold_eV: float,
    *,
    min_count: int = 20,
    max_weak_tail_fraction: float = 0.05,
) -> dict[str, float | int | str]:
    energy = np.asarray(energy_eV, dtype=float)
    values = np.asarray(counts, dtype=int)
    tail = energy >= float(threshold_eV)
    nonzero_tail = tail & (values > 0)
    if np.any(values >= min_count):
        max_resolved = float(np.max(energy[values >= min_count]))
    else:
        max_resolved = 0.0
    if not np.any(nonzero_tail):
        min_tail_count = 0
        status = "insufficient"
        weak_fraction = 1.0
    else:
        min_tail_count = int(np.min(values[nonzero_tail]))
        tail_total = float(np.sum(values[tail]))
        weak_total = float(np.sum(values[tail & (values < min_count)]))
        weak_fraction = weak_total / max(tail_total, 1.0e-300)
        status = (
            "ok"
            if min_tail_count >= min_count
            or weak_fraction <= max_weak_tail_fraction
            else "insufficient"
        )
    return {
        "mc_tail_uncertainty_status": status,
        "mc_min_tail_bin_count": min_tail_count,
        "mc_tail_effective_sample_count_min": float(min_tail_count),
        "mc_tail_weak_probability_fraction": float(weak_fraction),
        "mc_max_resolved_energy_eV": max_resolved,
    }


def _mc_tail_uncertainty_metadata_from_effective_counts(
    energy_eV: np.ndarray,
    effective_counts: np.ndarray,
    threshold_eV: float,
    *,
    min_count: int = 20,
    bin_probability: np.ndarray | None = None,
    max_weak_tail_fraction: float = 0.05,
) -> dict[str, float | int | str]:
    energy = np.asarray(energy_eV, dtype=float)
    values = np.asarray(effective_counts, dtype=float)
    tail = energy >= float(threshold_eV)
    nonzero_tail = tail & (values > 0.0)
    if np.any(values >= float(min_count)):
        max_resolved = float(np.max(energy[values >= float(min_count)]))
    else:
        max_resolved = 0.0
    if not np.any(nonzero_tail):
        min_tail_count = 0.0
        status = "insufficient"
        weak_fraction = 1.0
    else:
        min_tail_count = float(np.min(values[nonzero_tail]))
        weak = tail & (values < float(min_count))
        if bin_probability is None:
            mass = values
        else:
            mass = np.asarray(bin_probability, dtype=float)
        tail_total = float(np.sum(mass[tail]))
        weak_total = float(np.sum(mass[weak]))
        weak_fraction = weak_total / max(tail_total, 1.0e-300)
        status = (
            "ok"
            if min_tail_count >= float(min_count)
            or weak_fraction <= max_weak_tail_fraction
            else "insufficient"
        )
    return {
        "mc_tail_uncertainty_status": status,
        "mc_min_tail_bin_count": float(min_tail_count),
        "mc_tail_effective_sample_count_min": float(min_tail_count),
        "mc_tail_weak_probability_fraction": float(weak_fraction),
        "mc_max_resolved_energy_eV": max_resolved,
    }


def _mc_tail_comparison_status(
    *,
    energy_balance_status: str,
    tail_uncertainty_status: str,
) -> str:
    if energy_balance_status == "fail":
        return "energy_balance_fail"
    if energy_balance_status == "warning":
        return "energy_balance_warning"
    if tail_uncertainty_status != "ok":
        return "weak_tail_statistics"
    return "ok"


def _tail_threshold_eV(cross_sections: CrossSectionSet) -> float:
    thresholds = [
        float(proc.threshold_eV)
        for proc in cross_sections.processes
        if proc.process_type
        in {
            ProcessType.EXCITATION,
            ProcessType.IONIZATION,
            ProcessType.ATTACHMENT,
            ProcessType.SUPERELASTIC,
        }
        and proc.threshold_eV is not None
        and np.isfinite(proc.threshold_eV)
        and proc.threshold_eV > 0.0
    ]
    return max(thresholds) if thresholds else 20.0


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
    population_model = cfg.population_model
    collisions = int(cfg.max_collisions or 100)
    warmup_collisions = int(cfg.warmup_collisions or 0)
    density = gas_number_density(config)
    projected = _project_processes(config, cross_sections)
    max_energy_limit = float(config.physics.energy_grid_policy.max_eV_limit)
    trial_frequency = _trial_collision_frequency(
        projected, density, max_energy_eV=max_energy_limit
    )
    max_cross_section_energy = _max_cross_section_energy(projected)
    zero_high_energy_policy = config.cross_sections.high_energy_extrapolation == "zero"
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
        audit = _EnergyAudit()
        initial_speed = _speed_from_energy(1.0)
        ensemble = _ParticleEnsemble.initialize(particles, initial_speed, rng)
        audit.tracked_particle_initial_energy_eV = (
            ensemble.total_weighted_energy_eV()
        )
        edges = _mc_energy_edges(max_energy_limit)
        widths = np.diff(edges)
        energy = 0.5 * (edges[:-1] + edges[1:])
        counts = np.zeros_like(energy, dtype=int)
        weighted_hist = np.zeros_like(energy, dtype=float)
        weighted_square_hist = np.zeros_like(energy, dtype=float)
        total_collision_steps = warmup_collisions + collisions
        for step in range(total_collision_steps):
            sampling_enabled = step >= warmup_collisions
            if step == warmup_collisions and warmup_collisions:
                ensemble.positions[:] = 0.0
                ensemble.times[:] = 0.0
            particles_this_round = len(ensemble)
            for p in range(particles_this_round):
                if p >= len(ensemble):
                    break
                dt = float(rng.exponential(1.0 / trial_frequency))
                if magnetic_enabled:
                    gyro = E_CHARGE_C * magnetic.B_T / ELECTRON_MASS_KG
                    dt = min(dt, 0.2 / max(gyro, 1.0))
                weight = float(ensemble.weights[p])
                before_push_energy = _energy_from_velocity(ensemble.velocities[p])
                ensemble.velocities[p] = boris_push(ensemble.velocities[p], E, B, dt)
                after_push_energy = _energy_from_velocity(ensemble.velocities[p])
                audit.record_field_push(before_push_energy, after_push_energy, weight)
                ensemble.positions[p] += ensemble.velocities[p] * dt
                ensemble.times[p] += dt
                energy_eV = max(_energy_from_velocity(ensemble.velocities[p]), 1.0e-6)
                if zero_high_energy_policy and energy_eV > max_cross_section_energy:
                    raise RuntimeError(
                        "internal_monte_carlo particle exceeded the cross-section "
                        "energy range while high_energy_extrapolation=zero; use "
                        "hold extrapolation or extend the cross-section table"
                    )
                if energy_eV > max_energy_limit:
                    raise RuntimeError(
                        "internal_monte_carlo particle exceeded "
                        "physics.energy_grid_policy.max_eV_limit"
                    )
                bin_index = int(np.searchsorted(edges, energy_eV, side="right") - 1)
                if sampling_enabled and 0 <= bin_index < counts.size:
                    contribution = weight * dt
                    counts[bin_index] += 1
                    weighted_hist[bin_index] += contribution
                    weighted_square_hist[bin_index] += contribution * contribution
                speed = max(float(np.linalg.norm(ensemble.velocities[p])), 1.0)
                sigmas, collision_total, scattering_total, sigma_momentum = _sigma_arrays(
                    projected, energy_eV
                )
                if collision_total <= 0.0:
                    continue
                collision_frequency = density * collision_total * speed
                _validate_trial_collision_frequency(
                    collision_frequency, trial_frequency
                )
                collision_probability = collision_frequency / trial_frequency
                audit.record_trial_ratio(collision_probability)
                if rng.random() > collision_probability:
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
                    audit.record_elastic_collision(energy_eV, post_energy, weight)
                    ensemble.velocities[p] = _scatter_velocity(
                        ensemble.velocities[p], mu, rng, post_energy
                    )
                elif proc.process_type in TRAJECTORY_REACTION_TYPES:
                    outcome = _post_reaction_outcome(
                        config,
                        proc.process_type,
                        proc.threshold_eV,
                        energy_eV,
                        rng,
                    )
                    audit.record_reaction(outcome, weight)
                    ensemble.velocities[p] = _random_direction(
                        rng
                    ) * _speed_from_energy(outcome.tracked_energy_eV)

        audit.tracked_particle_final_energy_eV = ensemble.total_weighted_energy_eV()
        effective_counts = _mc_effective_bin_counts(weighted_hist, weighted_square_hist)
        total_hist_weight = float(np.sum(weighted_hist))
        eedf = weighted_hist / max(total_hist_weight, 1.0e-300) / widths
        if not np.any(eedf > 0.0):
            eedf[0] = 1.0 / widths[0]
        eedf = eedf / max(float(np.sum(eedf * widths)), 1.0e-300)
        nonzero_bins = counts > 0
        max_nonzero_energy = (
            float(np.max(energy[nonzero_bins])) if np.any(nonzero_bins) else 0.0
        )
        tail_uncertainty = _mc_tail_uncertainty_metadata_from_effective_counts(
            energy,
            effective_counts,
            _tail_threshold_eV(cross_sections),
            bin_probability=eedf * widths,
        )
        mean_energy = float(np.sum(energy * eedf * widths))
        total_weight = max(float(np.sum(ensemble.weights)), 1.0e-300)
        mean_time = max(
            float(np.sum(ensemble.weights * ensemble.times) / total_weight), 1.0e-300
        )
        drift_velocity = -float(
            np.sum(
                ensemble.weights
                * ensemble.positions[:, 2]
                / np.maximum(ensemble.times, 1.0e-300)
            )
            / total_weight
        )
        mobility = drift_velocity / max(E_scalar, 1.0e-300)
        dz = -ensemble.positions[:, 2] - drift_velocity * ensemble.times
        diffusion_l = float(
            np.sum(ensemble.weights * dz * dz) / total_weight / (2.0 * mean_time)
        )
        transverse = ensemble.positions[:, :2]
        diffusion_t = float(
            np.sum(ensemble.weights * np.sum(transverse * transverse, axis=1))
            / total_weight
            / (4.0 * mean_time)
        )
        rates, net_freq, net_rate = _rate_results(
            config, cross_sections, energy, eedf, widths, case_id, float(e_over_n_Td)
        )
        effective_townsend = net_freq / max(abs(drift_velocity) * density, 1.0e-300)
        audit_metadata = audit.as_metadata()
        tail_comparison_status = _mc_tail_comparison_status(
            energy_balance_status=str(audit_metadata["mc_energy_balance_status"]),
            tail_uncertainty_status=str(tail_uncertainty["mc_tail_uncertainty_status"]),
        )
        metadata = {
            "adapter": "internal_monte_carlo",
            "monte_carlo_backend": "internal",
            "particles": particles,
            "warmup_collisions": warmup_collisions,
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
            "trajectory_reaction_model": "energy_loss_single_daughter",
            "mc_population_model": population_model,
            "mc_warmup_collisions": warmup_collisions,
            "mc_production_collisions": collisions,
            "ionization_source_model": config.physics.ionization.energy_sharing,
            "ionization_source_treatment": config.physics.ionization.energy_sharing,
            "secondary_electron_tracking": False,
            "ionization_branching_model": _ionization_branching_model(config),
            "attachment_trajectory_treatment": "rate_convolution_only",
            "inelastic_angular_model": "isotropic_reset",
            "eedf_estimator": "time_sampled_null_clock",
            "mc_samples": int(np.sum(counts)),
            "mc_nonzero_bins": int(np.count_nonzero(nonzero_bins)),
            "mc_max_nonzero_energy_eV": max_nonzero_energy,
            "mc_energy_bin_max_eV": float(edges[-1]),
            "null_collision_trial_frequency_s_inv": float(trial_frequency),
            "mc_tail_comparison_status": tail_comparison_status,
            **audit_metadata,
            **tail_uncertainty,
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
                energy_widths_eV=widths,
                eedf_counts=counts.astype(int),
                eedf_effective_counts=effective_counts,
                rates=rates,
                metadata=metadata,
            )
        )
    return out
