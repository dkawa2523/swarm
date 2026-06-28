"""Limited internal Monte Carlo backend with Boris magnetic push."""

from __future__ import annotations

from dataclasses import dataclass
import logging

import numpy as np

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.constants import (
    BOLTZMANN_J_K,
    ELECTRON_MASS_KG,
    EV_TO_J,
    TOWNSEND,
)
from electron_swarm.core.cross_sections import (
    CrossSectionSet,
    ProcessType,
    mixture_fraction,
)
from electron_swarm.core.result_metadata import (
    TRANSPORT_MC_FIXED_POPULATION,
    TRANSPORT_MC_WEIGHTED_GROWTH,
)
from electron_swarm.core.results import RateResult, SwarmCaseResult
from electron_swarm.core.solver_configs import MonteCarloAdapterConfig
from electron_swarm.core.transport import ElectronTransport
from electron_swarm.solvers.base import SwarmSolver
from electron_swarm.physics.angular_scattering import (
    expected_angular_metadata,
    build_angular_model,
)
from electron_swarm.solvers._internal_mc.collisions import (
    SCATTERING_TYPES,
    TRAJECTORY_REACTION_TYPES,
    _PostReactionOutcome,
    _elastic_energy_after_collision,
    _energy_loss_eV,
    _ionization_daughters,
    _max_cross_section_energy,
    _moment_cross_sections,
    _post_reaction_outcome,
    _project_processes,
    _scatter_velocity,
    _sigma_arrays,
    _trial_collision_frequency,
    _validate_maxent_p1_cross_sections,
    _validate_trial_collision_frequency,
)
from electron_swarm.solvers._internal_mc.histogram import (
    _build_eedf_histogram,
    _mc_energy_edges,
    _mc_tail_comparison_status,
    _mc_tail_uncertainty_metadata_from_effective_counts,
    _tail_threshold_eV,
)
from electron_swarm.solvers._internal_mc.orbit import (
    _advance_to_trial_event,
    _energy_from_velocity,
    _random_direction,
    _speed_from_energy,
    magnetic_field_vector,
)


logger = logging.getLogger(__name__)



def gas_number_density(config: SwarmConfig) -> float:
    cond = config.conditions
    if cond.gas_number_density_m3 is not None:
        return float(cond.gas_number_density_m3)
    assert cond.pressure_Pa is not None
    return float(cond.pressure_Pa / (BOLTZMANN_J_K * cond.gas_temperature_K))


@dataclass(slots=True)
class _EnergyAudit:
    field_work_eV: float = 0.0
    elastic_energy_loss_eV: float = 0.0
    inelastic_energy_loss_eV: float = 0.0
    ionization_threshold_loss_eV: float = 0.0
    ionization_untracked_secondary_energy_eV: float = 0.0
    population_resampling_energy_adjustment_eV: float = 0.0
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

    def record_population_resampling(
        self, before_eV: float, after_eV: float
    ) -> None:
        self.population_resampling_energy_adjustment_eV += float(after_eV) - float(
            before_eV
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
            + self.population_resampling_energy_adjustment_eV
        )
        residual = tracked_delta - expected_delta
        scale = max(
            abs(tracked_delta),
            abs(self.field_work_eV)
            + abs(self.elastic_energy_loss_eV)
            + abs(self.inelastic_energy_loss_eV)
            + abs(self.ionization_threshold_loss_eV)
            + abs(self.ionization_untracked_secondary_energy_eV),
            abs(self.population_resampling_energy_adjustment_eV),
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
            "mc_population_resampling_energy_adjustment_eV": float(
                self.population_resampling_energy_adjustment_eV
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


class _NullEnergyAudit:
    tracked_particle_initial_energy_eV = 0.0
    tracked_particle_final_energy_eV = 0.0

    def record_field_push(
        self, before_eV: float, after_eV: float, weight: float = 1.0
    ) -> None:
        return None

    def record_trial_ratio(self, ratio: float) -> None:
        return None

    def record_elastic_collision(
        self, before_eV: float, after_eV: float, weight: float = 1.0
    ) -> None:
        return None

    def record_reaction(
        self, outcome: _PostReactionOutcome, weight: float = 1.0
    ) -> None:
        return None

    def record_population_resampling(
        self, before_eV: float, after_eV: float
    ) -> None:
        return None


@dataclass(slots=True)
class _MonteCarloRunAudit:
    seed: int | None
    particles: int
    population_model: str
    warmup_collisions: int
    production_collisions: int
    trial_collision_frequency_s_inv: float
    max_cross_section_energy_eV: float
    tail_threshold_eV: float
    max_sampled_energy_eV: float = 0.0
    energy_samples: int = 0
    energy_samples_above_xs_max: int = 0
    histogram_samples: int = 0
    tail_histogram_samples: int = 0
    collision_count_elastic: int = 0
    collision_count_excitation: int = 0
    collision_count_ionization: int = 0
    collision_count_attachment: int = 0
    collision_count_superelastic: int = 0
    collision_count_null: int = 0
    orbit_substeps: int = 0
    secondary_electron_count: int = 0
    branching_resample_count: int = 0
    population_total_weight_final: float = 0.0
    population_log_growth_estimate_s_inv: float = 0.0
    population_weight_cv: float = 0.0

    def record_energy_sample(self, energy_eV: float) -> None:
        energy = float(energy_eV)
        if np.isfinite(energy):
            self.max_sampled_energy_eV = max(self.max_sampled_energy_eV, energy)
        self.energy_samples += 1
        if energy > self.max_cross_section_energy_eV:
            self.energy_samples_above_xs_max += 1

    def record_histogram_sample(self, energy_eV: float) -> None:
        self.histogram_samples += 1
        if float(energy_eV) >= self.tail_threshold_eV:
            self.tail_histogram_samples += 1

    def record_null_collision(self) -> None:
        self.collision_count_null += 1

    def record_orbit_substep(self) -> None:
        self.orbit_substeps += 1

    def record_collision(self, process_type: ProcessType) -> None:
        if process_type in SCATTERING_TYPES:
            self.collision_count_elastic += 1
        elif process_type == ProcessType.EXCITATION:
            self.collision_count_excitation += 1
        elif process_type == ProcessType.IONIZATION:
            self.collision_count_ionization += 1
        elif process_type == ProcessType.ATTACHMENT:
            self.collision_count_attachment += 1
        elif process_type == ProcessType.SUPERELASTIC:
            self.collision_count_superelastic += 1

    def record_secondary_electron(self) -> None:
        self.secondary_electron_count += 1

    def record_branching_resample(self) -> None:
        self.branching_resample_count += 1

    def set_population_state(
        self, weights: np.ndarray, elapsed_time_s: float
    ) -> None:
        values = np.asarray(weights, dtype=float)
        total = float(np.sum(values))
        self.population_total_weight_final = total
        mean = float(np.mean(values)) if values.size else 0.0
        self.population_weight_cv = (
            float(np.std(values) / mean) if mean > 0.0 else 0.0
        )
        self.population_log_growth_estimate_s_inv = float(
            np.log(max(total, 1.0e-300) / max(float(self.particles), 1.0e-300))
            / max(float(elapsed_time_s), 1.0e-300)
        )

    def as_metadata(self) -> dict[str, float | int | str | None]:
        above_fraction = self.energy_samples_above_xs_max / max(
            self.energy_samples, 1
        )
        tail_fraction = self.tail_histogram_samples / max(self.histogram_samples, 1)
        return {
            "mc_seed": self.seed,
            "mc_particles": int(self.particles),
            "mc_population_model": self.population_model,
            "mc_warmup_collisions": int(self.warmup_collisions),
            "mc_production_collisions": int(self.production_collisions),
            "mc_trial_collision_frequency_s_inv": float(
                self.trial_collision_frequency_s_inv
            ),
            "mc_max_sampled_energy_eV": float(self.max_sampled_energy_eV),
            "mc_max_cross_section_energy_eV": float(self.max_cross_section_energy_eV),
            "mc_energy_samples_above_xs_max": int(
                self.energy_samples_above_xs_max
            ),
            "mc_energy_samples_above_xs_max_fraction": float(above_fraction),
            "mc_histogram_samples": int(self.histogram_samples),
            "mc_tail_histogram_samples": int(self.tail_histogram_samples),
            "mc_tail_histogram_sample_fraction": float(tail_fraction),
            "mc_collision_count_elastic": int(self.collision_count_elastic),
            "mc_collision_count_excitation": int(self.collision_count_excitation),
            "mc_collision_count_ionization": int(self.collision_count_ionization),
            "mc_collision_count_attachment": int(self.collision_count_attachment),
            "mc_collision_count_superelastic": int(self.collision_count_superelastic),
            "mc_collision_count_null": int(self.collision_count_null),
            "mc_orbit_substeps": int(self.orbit_substeps),
            "mc_secondary_electron_count": int(self.secondary_electron_count),
            "mc_branching_resample_count": int(self.branching_resample_count),
            "mc_population_total_weight_final": float(
                self.population_total_weight_final
            ),
            "mc_population_log_growth_estimate_s_inv": float(
                self.population_log_growth_estimate_s_inv
            ),
            "mc_population_weight_cv": float(self.population_weight_cv),
        }


class _NullMonteCarloRunAudit:
    def record_energy_sample(self, energy_eV: float) -> None:
        return None

    def record_histogram_sample(self, energy_eV: float) -> None:
        return None

    def record_null_collision(self) -> None:
        return None

    def record_orbit_substep(self) -> None:
        return None

    def record_collision(self, process_type: ProcessType) -> None:
        return None

    def record_secondary_electron(self) -> None:
        return None

    def record_branching_resample(self) -> None:
        return None

    def set_population_state(
        self, weights: np.ndarray, elapsed_time_s: float
    ) -> None:
        return None


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

    def append_particle(
        self,
        *,
        position: np.ndarray,
        velocity: np.ndarray,
        time_s: float,
        weight: float,
    ) -> None:
        self.positions = np.vstack([self.positions, np.asarray(position, dtype=float)])
        self.velocities = np.vstack(
            [self.velocities, np.asarray(velocity, dtype=float)]
        )
        self.times = np.append(self.times, float(time_s))
        self.weights = np.append(self.weights, float(weight))

    def normalize_total_weight(self, target_total_weight: float) -> None:
        total = float(np.sum(self.weights))
        if total > 0.0:
            self.weights *= float(target_total_weight) / total

    def systematic_resample(
        self,
        target_particles: int,
        rng: np.random.Generator,
    ) -> bool:
        target = int(target_particles)
        if target <= 0 or len(self) <= target:
            return False
        total = float(np.sum(self.weights))
        if not np.isfinite(total) or total <= 0.0:
            return False
        cumulative = np.cumsum(self.weights)
        positions = (float(rng.random()) + np.arange(target, dtype=float)) * (
            total / target
        )
        indices = np.searchsorted(cumulative, positions, side="right")
        indices = np.clip(indices, 0, len(self) - 1)
        self.positions = self.positions[indices].copy()
        self.velocities = self.velocities[indices].copy()
        self.times = self.times[indices].copy()
        self.weights = np.full(target, total / target, dtype=float)
        return True

def _weighted_branching_active(config: SwarmConfig, population_model: str) -> bool:
    return (
        population_model == "weighted_branching"
        and config.physics.ionization.energy_sharing != "loss_only"
    )


def _ionization_branching_model(config: SwarmConfig, population_model: str) -> str:
    if config.physics.ionization.energy_sharing == "loss_only":
        return "none"
    if population_model == "weighted_branching":
        return "two_daughter_weighted_resampling"
    return "single_daughter_sampling"


def _population_metadata(
    config: SwarmConfig, population_model: str
) -> dict[str, str | bool]:
    if _weighted_branching_active(config, population_model):
        return {
            "swarm_population_treatment": "weighted_branching_resampled",
            "nonconservative_growth_treatment": "explicit_weighted_branching",
            "secondary_electron_tracking": True,
            "ionization_branching_model": _ionization_branching_model(
                config, population_model
            ),
            "transport_definition": TRANSPORT_MC_WEIGHTED_GROWTH,
            "transport_has_bulk": False,
        }
    if population_model == "weighted_branching":
        return {
            "swarm_population_treatment": "fixed_population_loss_only_no_secondary",
            "nonconservative_growth_treatment": "not_tracked_loss_only",
            "secondary_electron_tracking": False,
            "ionization_branching_model": "none",
            "transport_definition": TRANSPORT_MC_FIXED_POPULATION,
            "transport_has_bulk": False,
        }
    return {
        "swarm_population_treatment": "fixed_population_no_growth",
        "nonconservative_growth_treatment": "not_tracked",
        "secondary_electron_tracking": False,
        "ionization_branching_model": _ionization_branching_model(
            config, population_model
        ),
        "transport_definition": TRANSPORT_MC_FIXED_POPULATION,
        "transport_has_bulk": False,
    }


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
                target_species_fraction=frac,
                energy_loss_eV=_energy_loss_eV(proc.process_type, proc.threshold_eV),
                gas_number_density_m3=density,
            )
        )
    net_freq = density * (ion_rate - attach_rate)
    return rates, net_freq, ion_rate - attach_rate


def _flux_transport_estimates(
    ensemble: _ParticleEnsemble,
    electric_field_V_m: float,
) -> tuple[float, float, float, float, float]:
    weights = np.asarray(ensemble.weights, dtype=float)
    times = np.asarray(ensemble.times, dtype=float)
    positions = np.asarray(ensemble.positions, dtype=float)
    total_weight = max(float(np.sum(weights)), 1.0e-300)
    total_time = max(float(np.sum(weights * times)), 1.0e-300)
    mean_time = total_time / total_weight
    drift_velocity = -float(np.sum(weights * positions[:, 2]) / total_time)
    mobility = drift_velocity / max(float(electric_field_V_m), 1.0e-300)
    dz = -positions[:, 2] - drift_velocity * times
    diffusion_l = float(np.sum(weights * dz * dz) / (2.0 * total_time))
    transverse = positions[:, :2]
    diffusion_t = float(
        np.sum(weights * np.sum(transverse * transverse, axis=1))
        / (4.0 * total_time)
    )
    return drift_velocity, mobility, diffusion_l, diffusion_t, mean_time


def run_internal_monte_carlo(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    solver_config: MonteCarloAdapterConfig,
    *,
    collect_audit: bool = False,
) -> list[SwarmCaseResult]:
    cfg = solver_config
    rng = np.random.default_rng(cfg.seed)
    particles = int(cfg.particles or 256)
    population_model = cfg.population_model
    branching_active = _weighted_branching_active(config, population_model)
    collisions = int(cfg.max_collisions or 100)
    warmup_collisions = int(cfg.warmup_collisions or 0)
    density = gas_number_density(config)
    projected = _project_processes(config, cross_sections)
    angular_model = build_angular_model(config)
    angular_name = config.physics.angular_scattering.model
    if angular_name not in {"isotropic", "maxent_p1"}:
        raise NotImplementedError(
            f"internal monte_carlo has no product sampler for angular model {angular_name!r}"
        )
    if angular_name == "maxent_p1":
        _validate_maxent_p1_cross_sections(projected)
    max_energy_limit = float(config.physics.energy_grid_policy.max_eV_limit)
    trial_frequency = _trial_collision_frequency(
        projected,
        density,
        angular_model_name=angular_name,
        max_energy_eV=max_energy_limit,
    )
    max_cross_section_energy = _max_cross_section_energy(projected)
    tail_threshold = _tail_threshold_eV(cross_sections)
    zero_high_energy_policy = config.cross_sections.high_energy_extrapolation == "zero"

    magnetic = config.physics.field.magnetic_field
    B = magnetic_field_vector(magnetic.B_T, magnetic.angle_EB_deg)
    magnetic_enabled = bool(magnetic.enabled and magnetic.B_T > 0.0)
    out: list[SwarmCaseResult] = []
    for i, e_over_n_Td in enumerate(config.run.e_over_n_Td):
        case_id = f"{config.run.case_prefix}_{i:04d}"
        E_scalar = float(e_over_n_Td) * TOWNSEND * density
        E = np.array([0.0, 0.0, E_scalar], dtype=float)
        audit = _EnergyAudit() if collect_audit else _NullEnergyAudit()
        run_audit = (
            _MonteCarloRunAudit(
                seed=cfg.seed,
                particles=particles,
                population_model=population_model,
                warmup_collisions=warmup_collisions,
                production_collisions=collisions,
                trial_collision_frequency_s_inv=trial_frequency,
                max_cross_section_energy_eV=max_cross_section_energy,
                tail_threshold_eV=tail_threshold,
            )
            if collect_audit
            else _NullMonteCarloRunAudit()
        )
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
                if branching_active:
                    ensemble.normalize_total_weight(float(particles))
                    audit = _EnergyAudit() if collect_audit else _NullEnergyAudit()
                    audit.tracked_particle_initial_energy_eV = (
                        ensemble.total_weighted_energy_eV()
                    )
            particles_this_round = len(ensemble)
            for p in range(particles_this_round):
                if p >= len(ensemble):
                    break
                trial_dt = float(rng.exponential(1.0 / trial_frequency))
                weight = float(ensemble.weights[p])
                energy_eV = _advance_to_trial_event(
                    ensemble=ensemble,
                    particle_index=p,
                    trial_dt_s=trial_dt,
                    electric_field_V_m=E,
                    magnetic_field_T=B,
                    magnetic_enabled=magnetic_enabled,
                    magnetic_B_T=magnetic.B_T,
                    audit=audit,
                    run_audit=run_audit,
                    edges=edges,
                    counts=counts,
                    weighted_hist=weighted_hist,
                    weighted_square_hist=weighted_square_hist,
                    sampling_enabled=sampling_enabled,
                    zero_high_energy_policy=zero_high_energy_policy,
                    max_cross_section_energy_eV=max_cross_section_energy,
                    max_energy_limit_eV=max_energy_limit,
                )
                speed = max(float(np.linalg.norm(ensemble.velocities[p])), 1.0)
                sigmas, collision_total = _sigma_arrays(
                    projected,
                    energy_eV,
                    angular_model_name=angular_name,
                )
                if collision_total <= 0.0:
                    if sampling_enabled:
                        run_audit.record_null_collision()
                    continue
                collision_frequency = density * collision_total * speed
                _validate_trial_collision_frequency(
                    collision_frequency, trial_frequency
                )
                collision_probability = collision_frequency / trial_frequency
                audit.record_trial_ratio(collision_probability)
                if rng.random() > collision_probability:
                    if sampling_enabled:
                        run_audit.record_null_collision()
                    continue
                choice = int(rng.choice(len(projected), p=sigmas / collision_total))
                proc = projected[choice].process
                if sampling_enabled:
                    run_audit.record_collision(proc.process_type)
                if proc.process_type in SCATTERING_TYPES:
                    if angular_name == "maxent_p1":
                        total, momentum = _moment_cross_sections(
                            projected,
                            proc.species,
                            energy_eV,
                        )
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
                    if proc.process_type == ProcessType.IONIZATION and branching_active:
                        daughters = _ionization_daughters(
                            config,
                            proc.threshold_eV,
                            energy_eV,
                        )
                        outcome = _PostReactionOutcome(
                            tracked_energy_eV=float(sum(daughters.energies_eV)),
                            ionization_threshold_loss_eV=daughters.threshold_loss_eV,
                        )
                        audit.record_reaction(outcome, weight)
                        energies = daughters.energies_eV
                        ensemble.velocities[p] = _random_direction(
                            rng
                        ) * _speed_from_energy(max(energies[0], 1.0e-4))
                        if len(energies) > 1:
                            ensemble.append_particle(
                                position=ensemble.positions[p].copy(),
                                velocity=_random_direction(rng)
                                * _speed_from_energy(max(energies[1], 1.0e-4)),
                                time_s=float(ensemble.times[p]),
                                weight=weight,
                            )
                            if sampling_enabled:
                                run_audit.record_secondary_electron()
                    else:
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
            if branching_active and len(ensemble) > 2 * particles:
                before_resample_energy = ensemble.total_weighted_energy_eV()
                if ensemble.systematic_resample(particles, rng):
                    after_resample_energy = ensemble.total_weighted_energy_eV()
                    audit.record_population_resampling(
                        before_resample_energy, after_resample_energy
                    )
                    if sampling_enabled:
                        run_audit.record_branching_resample()

        audit.tracked_particle_final_energy_eV = ensemble.total_weighted_energy_eV()
        energy, widths, eedf, counts, effective_counts = _build_eedf_histogram(
            edges,
            counts,
            weighted_hist,
            weighted_square_hist,
        )
        nonzero_bins = counts > 0
        max_nonzero_energy = (
            float(np.max(energy[nonzero_bins])) if np.any(nonzero_bins) else 0.0
        )
        mean_energy = float(np.sum(energy * eedf * widths))
        drift_velocity, mobility, diffusion_l, diffusion_t, mean_time = (
            _flux_transport_estimates(ensemble, E_scalar)
        )
        rates, net_freq, net_rate = _rate_results(
            config, cross_sections, energy, eedf, widths, case_id, float(e_over_n_Td)
        )
        effective_townsend = net_freq / max(abs(drift_velocity) * density, 1.0e-300)
        population_metadata = _population_metadata(config, population_model)
        mc_run_details = {
            "adapter": "internal_monte_carlo",
            "monte_carlo_backend": "internal",
            "mc_seed": cfg.seed,
            "mc_particles": int(particles),
            "mc_population_model": population_model,
            "magnetic_field_treatment": (
                "boris_lorentz_push" if magnetic.enabled else "none"
            ),
            "magnetic_field_B_T": float(magnetic.B_T),
            "magnetic_field_angle_EB_deg": float(magnetic.angle_EB_deg),
            "gas_number_density_m-3": float(density),
            "electric_field_V_m": float(E_scalar),
            "reaction_rates_source": "eedf_convolution",
            "ionization_source_model": config.physics.ionization.energy_sharing,
            "ionization_source_treatment": config.physics.ionization.energy_sharing,
            "eedf_estimator": "time_residence_midpoint",
            **population_metadata,
        }
        diagnostics = {}
        if collect_audit:
            run_audit.set_population_state(ensemble.weights, mean_time)
            audit_metadata = audit.as_metadata()
            run_audit_metadata = run_audit.as_metadata()
            tail_uncertainty = _mc_tail_uncertainty_metadata_from_effective_counts(
                energy,
                effective_counts,
                tail_threshold,
                bin_probability=eedf * widths,
            )
            tail_comparison_status = _mc_tail_comparison_status(
                energy_balance_status=str(audit_metadata["mc_energy_balance_status"]),
                tail_uncertainty_status=str(
                    tail_uncertainty["mc_tail_uncertainty_status"]
                ),
            )
            diagnostics["internal_monte_carlo_audit"] = {
                "mc_run": {
                    **run_audit_metadata,
                    **mc_run_details,
                    "mc_samples": int(np.sum(counts)),
                    "mc_nonzero_bins": int(np.count_nonzero(nonzero_bins)),
                    "mc_max_nonzero_energy_eV": max_nonzero_energy,
                    "mc_energy_bin_max_eV": float(edges[-1]),
                },
                "mc_energy_audit": audit_metadata,
                "mc_tail_audit": {
                    **tail_uncertainty,
                    "mc_tail_comparison_status": tail_comparison_status,
                },
            }
            logger.info(
                "internal MC case %s: seed=%s particles=%d samples=%d "
                "mean_energy_eV=%.6g max_sampled_energy_eV=%.6g "
                "xs_above_fraction=%.3g",
                case_id,
                cfg.seed,
                particles,
                int(run_audit_metadata["mc_histogram_samples"] or 0),
                mean_energy,
                float(run_audit_metadata["mc_max_sampled_energy_eV"] or 0.0),
                float(
                    run_audit_metadata[
                        "mc_energy_samples_above_xs_max_fraction"
                    ]
                    or 0.0
                ),
            )
        else:
            logger.info(
                "internal MC case %s: seed=%s particles=%d samples=%d "
                "mean_energy_eV=%.6g",
                case_id,
                cfg.seed,
                particles,
                int(np.sum(counts)),
                mean_energy,
            )
        metadata = {
            "magnetic_field_treatment": (
                "boris_lorentz_push" if magnetic.enabled else "none"
            ),
            "magnetic_field_B_T": float(magnetic.B_T),
            "magnetic_field_angle_EB_deg": float(magnetic.angle_EB_deg),
            "ionization_source_model": config.physics.ionization.energy_sharing,
            "ionization_source_treatment": config.physics.ionization.energy_sharing,
            "transport_definition": population_metadata["transport_definition"],
            **expected_angular_metadata(config),
        }
        transport = ElectronTransport.from_actual(
            definition=str(population_metadata["transport_definition"]),
            gas_number_density_m3=density,
            drift_velocity_m_s=drift_velocity,
            mobility_m2_V_s=mobility,
            diffusion_L_m2_s=diffusion_l,
            diffusion_T_m2_s=diffusion_t,
        )
        out.append(
            SwarmCaseResult(
                solver="monte_carlo",
                case_id=case_id,
                e_over_n_Td=float(e_over_n_Td),
                mean_energy_eV=mean_energy,
                net_ionization_frequency_s=net_freq,
                effective_townsend_m2=effective_townsend,
                transport=transport,
                energy_eV=energy,
                eedf=eedf,
                energy_widths_eV=widths,
                eedf_counts=counts.astype(int),
                eedf_effective_counts=effective_counts,
                rates=rates,
                metadata=metadata,
                diagnostics=diagnostics,
            )
        )
    return out


class MonteCarloSolver(SwarmSolver):
    name = "monte_carlo"

    def __init__(
        self,
        config: SwarmConfig,
        cross_sections: CrossSectionSet,
        solver_config: MonteCarloAdapterConfig,
    ) -> None:
        super().__init__(config, cross_sections)
        self.solver_config = solver_config

    def solve_all(self) -> list[SwarmCaseResult]:
        return run_internal_monte_carlo(
            self.config,
            self.cross_sections,
            self.solver_config,
            collect_audit=self.solver_config.collect_audit,
        )

    def solve_case(self, e_over_n_Td: float, case_id: str) -> SwarmCaseResult:
        raise NotImplementedError("monte_carlo executes through solve_all()")
