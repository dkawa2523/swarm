"""Energy-balance and run diagnostics for the internal Monte Carlo solver."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from electron_swarm.core.cross_sections import ProcessType
from electron_swarm.solvers.monte_carlo.collisions import (
    SCATTERING_TYPES,
    _PostReactionOutcome,
)
from electron_swarm.solvers.monte_carlo.weighted_ensemble import WeightedEnsemblePlan


@dataclass(slots=True)
class _EnergyAudit:
    enabled = True
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

    def record_population_resampling(self, before_eV: float, after_eV: float) -> None:
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
            "mc_ionization_threshold_loss_eV": float(self.ionization_threshold_loss_eV),
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
            "mc_max_collision_to_trial_ratio": float(self.max_collision_to_trial_ratio),
        }


class _NullEnergyAudit:
    enabled = False
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

    def record_population_resampling(self, before_eV: float, after_eV: float) -> None:
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
    collision_count_superelastic: int = 0
    collision_count_null: int = 0
    orbit_substeps: int = 0
    secondary_electron_count: int = 0
    branching_resample_count: int = 0
    tail_weighted_resample_count: int = 0
    tail_weighted_occupied_strata_max: int = 0
    tail_particles_before_resample_max: int = 0
    tail_particles_after_resample_max: int = 0
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
        elif process_type == ProcessType.SUPERELASTIC:
            self.collision_count_superelastic += 1

    def record_secondary_electron(self) -> None:
        self.secondary_electron_count += 1

    def record_branching_resample(self) -> None:
        self.branching_resample_count += 1

    def record_tail_weighted_resample(self, plan: WeightedEnsemblePlan) -> None:
        self.tail_weighted_resample_count += 1
        self.tail_weighted_occupied_strata_max = max(
            self.tail_weighted_occupied_strata_max,
            int(plan.occupied_strata),
        )
        self.tail_particles_before_resample_max = max(
            self.tail_particles_before_resample_max,
            int(plan.tail_particles_before),
        )
        self.tail_particles_after_resample_max = max(
            self.tail_particles_after_resample_max,
            int(plan.tail_particles_after),
        )

    def set_population_state(
        self,
        weights: np.ndarray,
        elapsed_time_s: float,
        cumulative_log_growth: float = 0.0,
    ) -> None:
        values = np.asarray(weights, dtype=float)
        total = float(np.sum(values))
        self.population_total_weight_final = total
        mean = float(np.mean(values)) if values.size else 0.0
        self.population_weight_cv = float(np.std(values) / mean) if mean > 0.0 else 0.0
        self.population_log_growth_estimate_s_inv = float(
            (
                float(cumulative_log_growth)
                + np.log(max(total, 1.0e-300) / max(float(self.particles), 1.0e-300))
            )
            / max(float(elapsed_time_s), 1.0e-300)
        )

    def as_metadata(self) -> dict[str, float | int | str | None]:
        above_fraction = self.energy_samples_above_xs_max / max(self.energy_samples, 1)
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
            "mc_energy_samples_above_xs_max": int(self.energy_samples_above_xs_max),
            "mc_energy_samples_above_xs_max_fraction": float(above_fraction),
            "mc_histogram_samples": int(self.histogram_samples),
            "mc_tail_histogram_samples": int(self.tail_histogram_samples),
            "mc_tail_histogram_sample_fraction": float(tail_fraction),
            "mc_collision_count_elastic": int(self.collision_count_elastic),
            "mc_collision_count_excitation": int(self.collision_count_excitation),
            "mc_collision_count_ionization": int(self.collision_count_ionization),
            "mc_collision_count_superelastic": int(self.collision_count_superelastic),
            "mc_collision_count_null": int(self.collision_count_null),
            "mc_orbit_substeps": int(self.orbit_substeps),
            "mc_secondary_electron_count": int(self.secondary_electron_count),
            "mc_branching_resample_count": int(self.branching_resample_count),
            "mc_tail_weighted_resample_count": int(self.tail_weighted_resample_count),
            "mc_tail_weighted_occupied_strata_max": int(
                self.tail_weighted_occupied_strata_max
            ),
            "mc_tail_particles_before_resample_max": int(
                self.tail_particles_before_resample_max
            ),
            "mc_tail_particles_after_resample_max": int(
                self.tail_particles_after_resample_max
            ),
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

    def record_tail_weighted_resample(self, plan: WeightedEnsemblePlan) -> None:
        return None

    def set_population_state(
        self,
        weights: np.ndarray,
        elapsed_time_s: float,
        cumulative_log_growth: float = 0.0,
    ) -> None:
        return None
