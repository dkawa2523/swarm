"""Trajectory-time reaction-rate estimators for the internal MC backend."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np

from electron_swarm.core.constants import ELECTRON_MASS_KG, EV_TO_J
from electron_swarm.core.cross_sections import (
    CrossSectionProcess,
    ProcessType,
    ScatteringRole,
)
from electron_swarm.solvers.monte_carlo.cross_section_table import (
    PreparedCrossSectionTable,
)
from electron_swarm.solvers.monte_carlo.flight_integration import (
    new_rate_integration_workspace,
)


MIN_REACTION_RATE_BLOCKS = 20
_ZERO_EVENT_CONFIDENCE = 0.95
_SCATTERING_ROLES = {
    ScatteringRole.ELASTIC_TOTAL,
    ScatteringRole.ELASTIC_MOMENTUM_TRANSFER,
    ScatteringRole.EFFECTIVE_MOMENTUM_TRANSFER,
}
_TOTAL_SCATTERING_ROLES = {
    ScatteringRole.ELASTIC_TOTAL,
}


@dataclass(frozen=True, slots=True)
class DirectEnergyLossEstimate:
    """Direct trajectory estimate of ``<sigma v Delta-epsilon>``."""

    status: str
    model: str
    rate_coefficient_eV_m3_s: float | None
    mean_energy_loss_per_collision_eV: float | None
    corresponding_rate_coefficient_m3_s: float
    standard_error_eV_m3_s: float | None
    relative_standard_error: float | None
    uncertainty_status: str
    event_loss_sample_count: int
    event_weighted_energy_loss_eV: float

    def metadata(self) -> dict[str, object]:
        thermal_elastic = self.model.startswith(
            "maxwellian_target_exact_binary_collision_"
        )
        return {
            "status": self.status,
            "estimator": (
                ("trajectory_event_energy_change_per_target_density_residence_time")
                if thermal_elastic and self.rate_coefficient_eV_m3_s is not None
                else "trajectory_time_average_sigma_v_delta_energy"
                if self.rate_coefficient_eV_m3_s is not None
                else None
            ),
            "model": self.model,
            "rate_coefficient_eV_m3_s": self.rate_coefficient_eV_m3_s,
            "mean_energy_loss_per_collision_eV": (
                self.mean_energy_loss_per_collision_eV
            ),
            "corresponding_rate_coefficient_m3_s": float(
                self.corresponding_rate_coefficient_m3_s
            ),
            "batch_means_standard_error_eV_m3_s": (self.standard_error_eV_m3_s),
            "batch_means_relative_standard_error": (self.relative_standard_error),
            "batch_means_uncertainty_status": self.uncertainty_status,
            "neutral_thermal_motion_model": (
                "maxwellian_relative_speed_exact_binary_collision"
                if thermal_elastic
                else None
            ),
            "gas_temperature_terms_included": (True if thermal_elastic else None),
        }


@dataclass(frozen=True, slots=True)
class ReactionRateEstimate:
    """One direct ``<sigma(E) v(E)>`` estimate and its sampling metadata."""

    rate_coefficient_m3_s: float
    block_count_requested: int
    block_count_used: int
    effective_sample_size: float | None
    standard_error_m3_s: float | None
    relative_standard_error: float | None
    lag1_autocorrelation: float | None
    autocorrelation_time_blocks: float | None
    uncertainty_status: str
    event_count: int
    weighted_event_count: float
    event_sampling_enabled: bool
    event_observation_status: str
    event_count_rate_coefficient_m3_s: float | None
    zero_event_rate_upper_95_m3_s: float | None
    event_observation_residence_time_s: float
    weighted_residence_time_s: float
    energy_loss: DirectEnergyLossEstimate

    def metadata(self) -> dict[str, object]:
        poisson_event_evidence = (
            self.event_observation_status
            != "weighted_ensemble_correlated_events_audit_only"
        )
        return {
            "estimator": "trajectory_time_average_sigma_v",
            "rate_coefficient_m3_s": float(self.rate_coefficient_m3_s),
            "batch_means_minimum_blocks": MIN_REACTION_RATE_BLOCKS,
            "batch_means_blocks_requested": int(self.block_count_requested),
            "batch_means_blocks_used": int(self.block_count_used),
            "batch_means_effective_sample_size": self.effective_sample_size,
            "batch_means_standard_error_m3_s": self.standard_error_m3_s,
            "batch_means_relative_standard_error": self.relative_standard_error,
            "batch_means_lag1_autocorrelation": self.lag1_autocorrelation,
            "batch_means_autocorrelation_time_blocks": (
                self.autocorrelation_time_blocks
            ),
            "batch_means_uncertainty_status": self.uncertainty_status,
            "event_count": int(self.event_count),
            "weighted_event_count": float(self.weighted_event_count),
            "event_sampling_enabled": bool(self.event_sampling_enabled),
            "event_observation_status": self.event_observation_status,
            "event_count_estimator": (
                "raw_events_per_target_density_unweighted_macro_trajectory_time"
                if poisson_event_evidence
                else "weighted_ensemble_correlated_raw_events_audit_only"
            ),
            "event_count_rate_coefficient_m3_s": (
                self.event_count_rate_coefficient_m3_s
            ),
            "zero_event_rate_upper_95_m3_s": (self.zero_event_rate_upper_95_m3_s),
            "zero_event_confidence": (
                _ZERO_EVENT_CONFIDENCE if poisson_event_evidence else None
            ),
            "zero_event_upper_bound_method": (
                "poisson_zero_count_unweighted_macro_trajectory_exposure"
                if poisson_event_evidence
                else None
            ),
            "event_observation_residence_time_s": float(
                self.event_observation_residence_time_s
            ),
            "weighted_residence_time_s": float(self.weighted_residence_time_s),
            "energy_loss": self.energy_loss.metadata(),
        }


@dataclass(frozen=True, slots=True)
class _CorrelatedMeanStatistics:
    effective_sample_size: float | None
    standard_error: float | None
    relative_standard_error: float | None
    lag1_autocorrelation: float | None
    autocorrelation_time: float | None


def _correlated_batch_mean_statistics(
    values: np.ndarray,
    *,
    reference_mean: float,
    weights: np.ndarray | None = None,
) -> _CorrelatedMeanStatistics:
    """Estimate SE using the initial-positive autocorrelation sequence.

    Negative autocorrelation is not used to claim extra information: the
    integrated autocorrelation time is floored at one block. This makes the
    reported ESS conservative for short MC runs.
    """

    samples = np.asarray(values, dtype=float)
    if weights is None:
        finite = np.isfinite(samples)
        sample_weights = np.ones(int(np.count_nonzero(finite)), dtype=float)
    else:
        raw_weights = np.asarray(weights, dtype=float)
        if raw_weights.shape != samples.shape:
            raise ValueError("reaction-rate block weights must match values")
        finite = np.isfinite(samples) & np.isfinite(raw_weights) & (raw_weights > 0.0)
        sample_weights = raw_weights[finite]
    samples = samples[finite]
    count = int(samples.size)
    if count < 2:
        return _CorrelatedMeanStatistics(None, None, None, None, None)

    centered = samples - float(np.mean(samples))
    variance_population = float(np.mean(centered * centered))
    weight_sum = float(np.sum(sample_weights))
    weight_square_sum = float(np.sum(sample_weights * sample_weights))
    independent_blocks = weight_sum * weight_sum / max(weight_square_sum, 1.0e-300)
    weighted_mean = float(np.sum(sample_weights * samples) / weight_sum)
    variance_denominator = weight_sum - weight_square_sum / weight_sum
    variance_sample = (
        float(
            np.sum(sample_weights * (samples - weighted_mean) ** 2)
            / variance_denominator
        )
        if variance_denominator > 0.0
        else 0.0
    )
    if variance_population <= 0.0 or not np.isfinite(variance_population):
        se = 0.0 if variance_sample == 0.0 else None
        relative = 0.0 if se == 0.0 and abs(reference_mean) > 0.0 else None
        return _CorrelatedMeanStatistics(
            float(independent_blocks), se, relative, 0.0, 1.0
        )

    positive_rhos: list[float] = []
    lag1: float | None = None
    # With twenty blocks, at most ten lags avoids the noisiest tail of the
    # empirical autocorrelation function.
    for lag in range(1, min(count // 2, count - 1) + 1):
        covariance = float(np.dot(centered[:-lag], centered[lag:]) / (count - lag))
        rho = float(np.clip(covariance / variance_population, -1.0, 1.0))
        if lag == 1:
            lag1 = rho
        if rho <= 0.0:
            break
        positive_rhos.append(rho)

    autocorrelation_time = max(1.0, 1.0 + 2.0 * sum(positive_rhos))
    effective_sample_size = float(
        np.clip(
            independent_blocks / autocorrelation_time,
            min(1.0, independent_blocks),
            independent_blocks,
        )
    )
    standard_error = math.sqrt(max(variance_sample, 0.0) / effective_sample_size)
    relative_standard_error = (
        standard_error / abs(reference_mean) if abs(reference_mean) > 0.0 else None
    )
    return _CorrelatedMeanStatistics(
        effective_sample_size,
        float(standard_error),
        (
            float(relative_standard_error)
            if relative_standard_error is not None
            else None
        ),
        float(lag1 if lag1 is not None else 0.0),
        float(autocorrelation_time),
    )


class TrajectoryReactionRateAccumulator:
    """Accumulate direct rate integrals in ordered trajectory-time blocks."""

    def __init__(
        self,
        *,
        processes: Sequence[CrossSectionProcess],
        fractions: Sequence[float],
        event_sampled: Sequence[bool],
        gas_number_density_m3: float,
        angular_scattering_model: str | None = None,
        block_count: int = MIN_REACTION_RATE_BLOCKS,
        poisson_event_evidence: bool = True,
    ) -> None:
        self.processes = tuple(processes)
        self.fractions = np.asarray(fractions, dtype=float)
        self.event_sampled = np.asarray(event_sampled, dtype=bool)
        process_count = len(self.processes)
        if process_count == 0:
            raise ValueError("reaction-rate accumulation requires at least one process")
        if self.fractions.shape != (process_count,):
            raise ValueError("reaction-rate fractions must match processes")
        if self.event_sampled.shape != (process_count,):
            raise ValueError("reaction-rate event flags must match processes")
        if np.any(~np.isfinite(self.fractions)) or np.any(self.fractions <= 0.0):
            raise ValueError("reaction-rate target fractions must be positive")
        density = float(gas_number_density_m3)
        if not np.isfinite(density) or density <= 0.0:
            raise ValueError("reaction-rate gas number density must be positive")
        blocks = int(block_count)
        if blocks < MIN_REACTION_RATE_BLOCKS:
            raise ValueError(
                f"reaction-rate batch means require at least "
                f"{MIN_REACTION_RATE_BLOCKS} blocks"
            )

        self.gas_number_density_m3 = density
        self.block_count = blocks
        self.poisson_event_evidence = bool(poisson_event_evidence)
        self._weighted_time = np.zeros(blocks, dtype=float)
        self._integrals = np.zeros((process_count, blocks), dtype=float)
        self._energy_loss_integrals = np.zeros((process_count, blocks), dtype=float)
        self._event_counts = np.zeros(process_count, dtype=np.int64)
        self._weighted_event_counts = np.zeros(process_count, dtype=float)
        self._event_loss_sample_counts = np.zeros(process_count, dtype=np.int64)
        self._weighted_event_energy_loss_eV = np.zeros(process_count, dtype=float)
        self._weighted_event_counts_by_block = np.zeros(
            (process_count, blocks),
            dtype=float,
        )
        self._weighted_event_energy_loss_eV_by_block = np.zeros(
            (process_count, blocks),
            dtype=float,
        )
        self._event_observation_time_s = 0.0
        self._prepared_cross_sections = PreparedCrossSectionTable.build(self.processes)
        thresholds = np.asarray(
            [
                float(process.threshold_eV)
                for process in self.processes
                if process.threshold_eV is not None
                and math.isfinite(float(process.threshold_eV))
            ],
            dtype=float,
        )
        self._flight_integration_breakpoints_eV = np.ascontiguousarray(
            np.unique(
                np.concatenate(
                    (self._prepared_cross_sections.grid_eV, thresholds)
                )
            )
        )
        (
            self._energy_loss_status,
            self._energy_loss_model,
            self._constant_energy_loss_eV,
            self._recoil_mass_ratio,
        ) = self._build_energy_loss_models(angular_scattering_model)
        self._angular_scattering_model = angular_scattering_model
        self._energy_loss_buffer = np.full(process_count, np.nan, dtype=float)
        self._flight_rate_workspace = new_rate_integration_workspace(process_count)

    def _build_energy_loss_models(
        self, angular_scattering_model: str | None
    ) -> tuple[list[str], list[str], np.ndarray, np.ndarray]:
        """Build event-model-consistent loss expectations on the shared grid."""

        process_count = len(self.processes)
        status = ["unsupported_process_energy_loss_not_defined"] * process_count
        model = ["none"] * process_count
        constant = np.full(process_count, np.nan, dtype=float)
        recoil_mass_ratio = np.full(process_count, np.nan, dtype=float)

        for index, process in enumerate(self.processes):
            process_type = process.process_type
            if process_type in {
                ProcessType.EXCITATION,
                ProcessType.IONIZATION,
                ProcessType.SUPERELASTIC,
            }:
                if process.threshold_eV is None:
                    status[index] = "unsupported_missing_process_threshold"
                    model[index] = "threshold_not_available"
                    continue
                threshold = float(process.threshold_eV)
                constant[index] = (
                    -abs(threshold)
                    if process_type == ProcessType.SUPERELASTIC
                    else threshold
                )
                status[index] = "direct_residence_integral"
                model[index] = "event_model_fixed_threshold_energy_change"
                continue

            if process.scattering_role not in _SCATTERING_ROLES:
                continue
            if not bool(self.event_sampled[index]):
                status[index] = "unsupported_not_sampled_by_trajectory_model"
                model[index] = "companion_cross_section_not_an_event_process"
                continue
            mass_amu = process.mass_amu
            if (
                mass_amu is None
                or not math.isfinite(float(mass_amu))
                or mass_amu <= 0.0
            ):
                status[index] = "unsupported_missing_target_mass"
                model[index] = "exact_binary_collision_requires_target_mass"
                continue
            if angular_scattering_model == "isotropic":
                model[index] = "maxwellian_target_exact_binary_collision_isotropic"
            elif angular_scattering_model == "maxent_p1":
                model[index] = "maxwellian_target_exact_binary_collision_maxent_p1"
            else:
                status[index] = "unsupported_angular_scattering_model"
                model[index] = str(angular_scattering_model or "not_provided")
                continue
            status[index] = "direct_event_estimator"

        return status, model, constant, recoil_mass_ratio

    def _evaluate_energy_loss_eV(
        self, energy_eV: float, sigma_v: np.ndarray
    ) -> np.ndarray:
        """Evaluate the physical event-model loss, without an EEDF closure."""

        energy = max(float(energy_eV), 0.0)
        out = self._energy_loss_buffer
        out.fill(np.nan)
        for index in range(len(self.processes)):
            if self._energy_loss_status[index] != "direct_residence_integral":
                continue
            fixed = self._constant_energy_loss_eV[index]
            if math.isfinite(float(fixed)):
                out[index] = float(fixed)
                continue
            mass_ratio = self._recoil_mass_ratio[index]
            if not math.isfinite(float(mass_ratio)):
                continue
            mean_one_minus_mu = 1.0
            if self._angular_scattering_model == "maxent_p1":
                species = self.processes[index].species
                total = sum(
                    float(sigma_v[other_index])
                    for other_index, other in enumerate(self.processes)
                    if other.species == species
                    and other.scattering_role in _TOTAL_SCATTERING_ROLES
                )
                momentum = sum(
                    float(sigma_v[other_index])
                    for other_index, other in enumerate(self.processes)
                    if other.species == species
                    and other.scattering_role
                    == ScatteringRole.ELASTIC_MOMENTUM_TRANSFER
                )
                mean_one_minus_mu = (
                    float(np.clip(momentum / total, 0.0, 2.0)) if total > 0.0 else 0.0
                )
            out[index] = energy * float(mass_ratio) * mean_one_minus_mu
        return out

    def _evaluate_sigma_v(self, energy_eV: float, speed_m_s: float) -> np.ndarray:
        """Evaluate all process cross sections with one shared interval lookup."""
        return self._prepared_cross_sections.evaluate(
            float(energy_eV),
            scale=float(speed_m_s),
        )

    def add_compiled_event_observation_time(self, residence_time_s: float) -> None:
        """Commit unweighted macro-trajectory exposure from a compiled block."""

        value = float(residence_time_s)
        if not math.isfinite(value) or value < 0.0:
            raise ValueError("compiled reaction-rate residence time is invalid")
        self._event_observation_time_s += value

    def block_index_for_trial(
        self,
        *,
        production_step: int,
        production_steps: int,
        particle_index: int,
        particles_this_round: int,
    ) -> int:
        """Map a trial event to an ordered, approximately equal-time batch."""

        steps = int(production_steps)
        particles = int(particles_this_round)
        step = int(production_step)
        particle = int(particle_index)
        if steps <= 0 or particles <= 0:
            raise ValueError("reaction-rate batch schedule must be positive")
        if step < 0 or step >= steps or particle < 0 or particle >= particles:
            raise IndexError("reaction-rate trial coordinate is out of range")
        progress = (step + (particle + 0.5) / particles) / steps
        return min(int(progress * self.block_count), self.block_count - 1)

    def record_residence(
        self,
        *,
        block_index: int,
        energy_eV: float,
        weight: float,
        dt_s: float,
    ) -> None:
        block = int(block_index)
        if block < 0 or block >= self.block_count:
            raise IndexError("reaction-rate block index is out of range")
        residence_time = float(dt_s)
        particle_weight = float(weight)
        if (
            not np.isfinite(residence_time)
            or residence_time <= 0.0
            or not np.isfinite(particle_weight)
            or particle_weight <= 0.0
        ):
            return
        contribution = particle_weight * residence_time
        energy = max(float(energy_eV), 0.0)
        speed = math.sqrt(max(2.0 * energy * EV_TO_J / ELECTRON_MASS_KG, 0.0))
        sigma_v = self._evaluate_sigma_v(energy, speed)
        energy_loss_eV = self._evaluate_energy_loss_eV(energy, sigma_v)
        self._weighted_time[block] += contribution
        self._integrals[:, block] += contribution * sigma_v
        supported = np.isfinite(energy_loss_eV)
        self._energy_loss_integrals[supported, block] += (
            contribution * sigma_v[supported] * energy_loss_eV[supported]
        )
        # A macro-particle's statistical event opportunity is one trajectory,
        # irrespective of the physical-population weight it represents.  A
        # global rescaling of weights must not make a zero-event upper bound
        # spuriously tighter.
        self._event_observation_time_s += residence_time

    def record_integrated_residence(
        self,
        *,
        block_index: int,
        weight: float,
        dt_s: float,
        rate_integrals_m3: np.ndarray,
        energy_loss_integrals_eV_m3: np.ndarray,
    ) -> None:
        """Commit error-controlled rate integrals for one complete flight."""

        block = int(block_index)
        if block < 0 or block >= self.block_count:
            raise IndexError("reaction-rate block index is out of range")
        residence_time = float(dt_s)
        particle_weight = float(weight)
        if residence_time <= 0.0 or particle_weight <= 0.0:
            return
        rates = np.asarray(rate_integrals_m3, dtype=float)
        losses = np.asarray(energy_loss_integrals_eV_m3, dtype=float)
        expected_shape = (len(self.processes),)
        if rates.shape != expected_shape or losses.shape != expected_shape:
            raise ValueError("integrated reaction-rate arrays do not match processes")
        self._weighted_time[block] += particle_weight * residence_time
        self._integrals[:, block] += particle_weight * rates
        self._energy_loss_integrals[:, block] += particle_weight * losses
        self._event_observation_time_s += residence_time

    def record_event(
        self,
        process_index: int,
        weight: float,
        *,
        energy_loss_eV: float | None = None,
        block_index: int | None = None,
    ) -> None:
        index = int(process_index)
        if index < 0 or index >= len(self.processes):
            raise IndexError("reaction event process index is out of range")
        block: int | None = None
        if energy_loss_eV is not None and math.isfinite(float(energy_loss_eV)):
            if block_index is None:
                if self._energy_loss_status[index] == "direct_event_estimator":
                    raise ValueError(
                        "event-based energy-loss estimation requires a block index"
                    )
            else:
                block = int(block_index)
                if block < 0 or block >= self.block_count:
                    raise IndexError("reaction event block index is out of range")
        self._event_counts[index] += 1
        self._weighted_event_counts[index] += float(weight)
        if energy_loss_eV is not None and math.isfinite(float(energy_loss_eV)):
            self._event_loss_sample_counts[index] += 1
            self._weighted_event_energy_loss_eV[index] += float(weight) * float(
                energy_loss_eV
            )
            if block is not None:
                self._weighted_event_counts_by_block[index, block] += float(weight)
                self._weighted_event_energy_loss_eV_by_block[index, block] += float(
                    weight
                ) * float(energy_loss_eV)

    def merge_residence_statistics(
        self,
        other: "TrajectoryReactionRateAccumulator",
    ) -> None:
        """Pool unbiased residence integrals while retaining this event audit.

        A threshold-weighted tail ensemble has correlated event counts, so its
        raw events cannot strengthen the ordinary Poisson evidence collected
        by the main trajectory ensemble. Its residence-time ``sigma v`` and
        energy-loss integrals remain unbiased and can be pooled block by block.
        """

        def process_identity(item: CrossSectionProcess) -> tuple[object, ...]:
            return (
                item.species,
                item.process,
                item.process_type,
                item.threshold_eV,
            )

        if (
            self.block_count != other.block_count
            or len(self.processes) != len(other.processes)
            or tuple(map(process_identity, self.processes))
            != tuple(map(process_identity, other.processes))
            or not np.array_equal(self.fractions, other.fractions)
            or not np.array_equal(self.event_sampled, other.event_sampled)
            or not math.isclose(
                self.gas_number_density_m3,
                other.gas_number_density_m3,
                rel_tol=0.0,
                abs_tol=0.0,
            )
            or self._energy_loss_status != other._energy_loss_status
            or self._energy_loss_model != other._energy_loss_model
        ):
            raise ValueError("reaction-rate residence accumulators are incompatible")
        self._weighted_time += other._weighted_time
        self._integrals += other._integrals
        self._energy_loss_integrals += other._energy_loss_integrals
        self._weighted_event_counts_by_block += other._weighted_event_counts_by_block
        self._weighted_event_energy_loss_eV_by_block += (
            other._weighted_event_energy_loss_eV_by_block
        )

    def estimates(self) -> list[ReactionRateEstimate]:
        total_time = float(np.sum(self._weighted_time))
        event_observation_time = float(self._event_observation_time_s)
        valid_blocks = self._weighted_time > 0.0
        blocks_used = int(np.count_nonzero(valid_blocks))
        out: list[ReactionRateEstimate] = []
        for index, _process in enumerate(self.processes):
            total_integral = float(np.sum(self._integrals[index]))
            direct_rate = total_integral / max(total_time, 1.0e-300)
            block_values = (
                self._integrals[index, valid_blocks] / self._weighted_time[valid_blocks]
            )
            statistics = _correlated_batch_mean_statistics(
                block_values,
                reference_mean=direct_rate,
                weights=self._weighted_time[valid_blocks],
            )
            if total_time <= 0.0:
                uncertainty_status = "no_exposure"
            elif blocks_used < MIN_REACTION_RATE_BLOCKS:
                uncertainty_status = "insufficient_blocks"
            else:
                uncertainty_status = "ok"

            count = int(self._event_counts[index])
            weighted_count = float(self._weighted_event_counts[index])
            sampled = bool(self.event_sampled[index])
            target_density = self.gas_number_density_m3 * float(self.fractions[index])
            event_rate: float | None = None
            upper: float | None = None
            if not sampled:
                event_status = "not_sampled_by_trajectory_model"
            elif not self.poisson_event_evidence:
                event_status = "weighted_ensemble_correlated_events_audit_only"
            elif event_observation_time <= 0.0:
                event_status = "no_exposure"
            elif count == 0:
                event_status = "unobserved_zero_events_upper_bound_only"
                upper = -math.log(1.0 - _ZERO_EVENT_CONFIDENCE) / (
                    target_density * event_observation_time
                )
            else:
                event_status = "observed"
                event_rate = count / (target_density * event_observation_time)

            loss_status = self._energy_loss_status[index]
            loss_coefficient: float | None = None
            mean_loss: float | None = None
            loss_statistics = _CorrelatedMeanStatistics(None, None, None, None, None)
            if loss_status == "direct_residence_integral":
                total_loss_integral = float(np.sum(self._energy_loss_integrals[index]))
                loss_coefficient = total_loss_integral / max(total_time, 1.0e-300)
                loss_block_values = (
                    self._energy_loss_integrals[index, valid_blocks]
                    / self._weighted_time[valid_blocks]
                )
                loss_statistics = _correlated_batch_mean_statistics(
                    loss_block_values,
                    reference_mean=loss_coefficient,
                    weights=self._weighted_time[valid_blocks],
                )
                if abs(direct_rate) > 0.0:
                    mean_loss = loss_coefficient / direct_rate
            elif loss_status == "direct_event_estimator":
                target_density = self.gas_number_density_m3 * float(
                    self.fractions[index]
                )
                total_event_loss = float(
                    np.sum(self._weighted_event_energy_loss_eV_by_block[index])
                )
                total_weighted_events = float(
                    np.sum(self._weighted_event_counts_by_block[index])
                )
                loss_coefficient = total_event_loss / max(
                    target_density * total_time,
                    1.0e-300,
                )
                loss_block_values = self._weighted_event_energy_loss_eV_by_block[
                    index, valid_blocks
                ] / (target_density * self._weighted_time[valid_blocks])
                loss_statistics = _correlated_batch_mean_statistics(
                    loss_block_values,
                    reference_mean=loss_coefficient,
                    weights=self._weighted_time[valid_blocks],
                )
                if total_weighted_events > 0.0:
                    mean_loss = total_event_loss / total_weighted_events
            if total_time <= 0.0:
                loss_uncertainty_status = "no_exposure"
            elif loss_status not in {
                "direct_residence_integral",
                "direct_event_estimator",
            }:
                loss_uncertainty_status = "not_identified"
            elif blocks_used < MIN_REACTION_RATE_BLOCKS:
                loss_uncertainty_status = "insufficient_blocks"
            else:
                loss_uncertainty_status = "ok"

            event_loss_samples = int(self._event_loss_sample_counts[index])
            event_weighted_loss = float(self._weighted_event_energy_loss_eV[index])

            out.append(
                ReactionRateEstimate(
                    rate_coefficient_m3_s=float(direct_rate),
                    block_count_requested=self.block_count,
                    block_count_used=blocks_used,
                    effective_sample_size=statistics.effective_sample_size,
                    standard_error_m3_s=statistics.standard_error,
                    relative_standard_error=statistics.relative_standard_error,
                    lag1_autocorrelation=statistics.lag1_autocorrelation,
                    autocorrelation_time_blocks=statistics.autocorrelation_time,
                    uncertainty_status=uncertainty_status,
                    event_count=count,
                    weighted_event_count=weighted_count,
                    event_sampling_enabled=sampled,
                    event_observation_status=event_status,
                    event_count_rate_coefficient_m3_s=event_rate,
                    zero_event_rate_upper_95_m3_s=upper,
                    event_observation_residence_time_s=event_observation_time,
                    weighted_residence_time_s=total_time,
                    energy_loss=DirectEnergyLossEstimate(
                        status=loss_status,
                        model=self._energy_loss_model[index],
                        rate_coefficient_eV_m3_s=loss_coefficient,
                        mean_energy_loss_per_collision_eV=mean_loss,
                        corresponding_rate_coefficient_m3_s=(
                            float(
                                np.sum(self._weighted_event_counts_by_block[index])
                                / max(target_density * total_time, 1.0e-300)
                            )
                            if loss_status == "direct_event_estimator"
                            else float(direct_rate)
                        ),
                        standard_error_eV_m3_s=loss_statistics.standard_error,
                        relative_standard_error=(
                            loss_statistics.relative_standard_error
                        ),
                        uncertainty_status=loss_uncertainty_status,
                        event_loss_sample_count=event_loss_samples,
                        event_weighted_energy_loss_eV=event_weighted_loss,
                    ),
                )
            )
        return out
