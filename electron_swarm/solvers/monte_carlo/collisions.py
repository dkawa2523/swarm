"""Collision sampling and reaction helpers for the internal Monte Carlo backend."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.constants import ELECTRON_MASS_KG, EV_TO_J
from electron_swarm.core.cross_sections import (
    ActiveMixtureInputs,
    CrossSectionProcess,
    ProcessType,
    ScatteringRole,
    gas_mass_amu,
)
from electron_swarm.core.scattering import resolve_particle_scattering
from electron_swarm.physics.electron_neutral import (
    IonizationDaughters,
    elastic_binary_collision_velocity_m_s,
    energy_from_speed_eV,
    ionization_daughters,
    maxwell_velocity_component_std_m_s,
)
from electron_swarm.solvers.monte_carlo.kinematics import (
    random_direction,
    speed_from_energy_m_s,
)
from electron_swarm.solvers.monte_carlo.cross_section_table import (
    PreparedCrossSectionTable,
)


SCATTERING_ROLES = {
    ScatteringRole.ELASTIC_TOTAL,
    ScatteringRole.ELASTIC_MOMENTUM_TRANSFER,
    ScatteringRole.EFFECTIVE_MOMENTUM_TRANSFER,
}
SCATTERING_TYPES = {
    ProcessType.MOMENTUM,
    ProcessType.ELASTIC,
    ProcessType.EFFECTIVE,
}
TOTAL_SCATTERING_ROLES = {ScatteringRole.ELASTIC_TOTAL}
MOMENTUM_TRANSFER_ROLES = {ScatteringRole.ELASTIC_MOMENTUM_TRANSFER}
TRAJECTORY_REACTION_TYPES = {
    ProcessType.EXCITATION,
    ProcessType.IONIZATION,
    ProcessType.SUPERELASTIC,
}

THERMAL_PROPOSAL_MAXWELL = 0
THERMAL_PROPOSAL_SPEED_ENVELOPE = 1
_SPEED_FACTOR = 2.0 * EV_TO_J / ELECTRON_MASS_KG
_TRIAL_BOUND_SAFETY = 1.0 + 1.0e-12


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


@dataclass(frozen=True, slots=True)
class SampledCollisionTrial:
    """One marked null-collision proposal and its accepted physical event."""

    process_index: int | None
    acceptance_probability: float
    target_velocity_m_s: np.ndarray | None = None
    relative_energy_eV: float | None = None


def _maximum_speed_weighted_cross_section(
    energy_eV: np.ndarray,
    cross_section_m2: np.ndarray,
) -> float:
    """Return the exact maximum of ``sqrt(2 e E/m) sigma(E)``.

    The prepared total cross section is piecewise linear in energy, so every
    interior extremum is available analytically in addition to the knots.
    """

    energy = np.asarray(energy_eV, dtype=float)
    sigma = np.asarray(cross_section_m2, dtype=float)
    if (
        energy.ndim != 1
        or sigma.shape != energy.shape
        or energy.size < 2
        or np.any(~np.isfinite(energy))
        or np.any(np.diff(energy) <= 0.0)
        or np.any(~np.isfinite(sigma))
        or np.any(sigma < 0.0)
    ):
        raise ValueError("invalid piecewise-linear cross section for majorant")

    candidates = list(map(float, energy))
    for left, right, sigma_left, sigma_right in zip(
        energy[:-1], energy[1:], sigma[:-1], sigma[1:], strict=True
    ):
        slope = float((sigma_right - sigma_left) / (right - left))
        intercept = float(sigma_left - slope * left)
        if slope < 0.0:
            critical = -intercept / (3.0 * slope)
            if float(left) < critical < float(right):
                candidates.append(float(critical))
    points = np.asarray(candidates, dtype=float)
    values = np.interp(points, energy, sigma)
    speed = np.sqrt(np.maximum(_SPEED_FACTOR * points, 0.0))
    maximum = float(np.max(speed * values, initial=0.0))
    if not math.isfinite(maximum) or maximum < 0.0:
        raise FloatingPointError("thermal collision majorant is invalid")
    return maximum


def _restricted_piecewise_values(
    energy_eV: np.ndarray,
    values: np.ndarray,
    upper_eV: float,
    right_value: float,
) -> tuple[np.ndarray, np.ndarray]:
    upper = float(upper_eV)
    if not math.isfinite(upper) or upper <= 0.0:
        raise ValueError("collision majorant energy limit must be positive")
    points = np.asarray(energy_eV, dtype=float)
    selected = points[points < upper]
    restricted = np.unique(np.concatenate(([0.0], selected, [upper])))
    sigma = np.interp(
        restricted,
        points,
        np.asarray(values, dtype=float),
        left=0.0,
        right=float(right_value),
    )
    return restricted, sigma


def _sample_maxwellian_target_velocity(
    component_std_m_s: float,
    rng: np.random.Generator,
) -> np.ndarray:
    scale = float(component_std_m_s)
    return np.asarray(
        [
            rng.normal(0.0, scale),
            rng.normal(0.0, scale),
            rng.normal(0.0, scale),
        ],
        dtype=float,
    )


def _sample_speed_weighted_maxwellian_target_velocity(
    component_std_m_s: float,
    rng: np.random.Generator,
) -> np.ndarray:
    # Multiplying the Maxwell speed density by speed changes
    # y=u^2/(2 s^2) from Gamma(3/2,1) to Gamma(2,1).  A shape-two gamma is the
    # sum of two unit exponentials and avoids a backend-specific gamma sampler.
    first = max(float(rng.random()), np.finfo(float).tiny)
    second = max(float(rng.random()), np.finfo(float).tiny)
    y = -math.log(first * second)
    speed = float(component_std_m_s) * math.sqrt(2.0 * y)
    return random_direction(rng) * speed


@dataclass(slots=True)
class PreparedCollisionSampler:
    """Exact thermal-elastic marked null-collision sampler.

    Inelastic channels retain the stationary-target ``v sigma(E)`` model used
    by the deterministic solvers.  Elastic channels sample a Maxwellian target
    and use ``g sigma(E_rel)``.  Zero-extrapolated elastic data use the sharp
    global ``max(g sigma)`` envelope.  Held tails use a provably dominating
    speed-weighted Maxwell proposal instead of truncating the neutral tail.
    """

    table: PreparedCrossSectionTable
    process_collision_group: np.ndarray
    cold_bound_m3_s: float
    thermal_group_bounds_m3_s: np.ndarray
    thermal_group_target_mass_amu: np.ndarray
    thermal_group_component_std_m_s: np.ndarray
    thermal_group_mean_speed_m_s: np.ndarray
    thermal_group_proposal_mode: np.ndarray
    thermal_group_sigma_max_m2: np.ndarray
    max_electron_speed_m_s: float
    max_energy_eV: float
    trial_rate_coefficient_m3_s: float

    @classmethod
    def build(
        cls,
        config: SwarmConfig,
        processes: list[_ProjectedProcess],
        *,
        angular_model_name: str,
        max_energy_eV: float,
    ) -> "PreparedCollisionSampler":
        multipliers = tuple(
            item.fraction
            if (
                _is_scattering_event_process(item, processes, angular_model_name)
                or item.process.process_type in TRAJECTORY_REACTION_TYPES
            )
            else 0.0
            for item in processes
        )
        table = PreparedCrossSectionTable.build(
            tuple(item.process for item in processes),
            multipliers=multipliers,
        )

        scattering = np.asarray(
            [
                _is_scattering_event_process(item, processes, angular_model_name)
                for item in processes
            ],
            dtype=bool,
        )
        cold = np.asarray(
            [
                item.process.process_type in TRAJECTORY_REACTION_TYPES
                for item in processes
            ],
            dtype=bool,
        )
        species = tuple(
            dict.fromkeys(
                item.process.species
                for item, active in zip(processes, scattering, strict=True)
                if active
            )
        )
        process_group = np.full(len(processes), -2, dtype=np.int64)
        process_group[cold] = -1
        for group, name in enumerate(species):
            process_group[
                np.asarray(
                    [
                        bool(active and item.process.species == name)
                        for item, active in zip(processes, scattering, strict=True)
                    ],
                    dtype=bool,
                )
            ] = group

        cold_values = (
            np.sum(table.values_m2[cold], axis=0)
            if np.any(cold)
            else np.zeros_like(table.grid_eV)
        )
        cold_right = float(np.sum(table.right_values_m2[cold])) if np.any(cold) else 0.0
        cold_energy, cold_sigma = _restricted_piecewise_values(
            table.grid_eV,
            cold_values,
            float(max_energy_eV),
            cold_right,
        )
        cold_bound = _maximum_speed_weighted_cross_section(
            cold_energy,
            cold_sigma,
        )

        max_speed = float(speed_from_energy_m_s(float(max_energy_eV)))
        temperature = float(config.conditions.gas_temperature_K)
        group_bounds: list[float] = []
        group_masses: list[float] = []
        group_stds: list[float] = []
        group_mean_speeds: list[float] = []
        group_modes: list[int] = []
        group_sigma_maxima: list[float] = []
        for group, name in enumerate(species):
            mask = process_group == group
            angular_companions = mask.copy()
            if angular_model_name == "maxent_p1":
                angular_companions = np.asarray(
                    [
                        item.process.species == name
                        and item.process.scattering_role in SCATTERING_ROLES
                        for item in processes
                    ],
                    dtype=bool,
                )
            if np.any(table.error_extrapolation & angular_companions):
                raise NotImplementedError(
                    "finite-temperature Monte Carlo elastic collisions require "
                    "zero or hold high-energy extrapolation because Maxwellian "
                    f"target velocities are unbounded; species={name}"
                )
            total_values = np.sum(table.values_m2[mask], axis=0)
            right_value = float(np.sum(table.right_values_m2[mask]))
            configured_mass = float(gas_mass_amu(config.conditions, name))
            declared_masses = {
                float(item.process.mass_amu or configured_mass)
                for item, active in zip(processes, mask, strict=True)
                if active
            }
            if len(declared_masses) != 1:
                raise ValueError(
                    "thermal elastic processes for one species must use one "
                    f"target mass; species={name}, masses={sorted(declared_masses)}"
                )
            target_mass = declared_masses.pop()
            component_std = maxwell_velocity_component_std_m_s(
                temperature,
                target_mass,
            )
            mean_speed = 2.0 * math.sqrt(2.0 / math.pi) * component_std
            sigma_max = max(
                float(np.max(total_values, initial=0.0)),
                right_value,
            )
            if right_value > 0.0:
                mode = THERMAL_PROPOSAL_SPEED_ENVELOPE
                bound = sigma_max * (max_speed + mean_speed)
            else:
                mode = THERMAL_PROPOSAL_MAXWELL
                bound = _maximum_speed_weighted_cross_section(
                    table.grid_eV,
                    total_values,
                )
            if not math.isfinite(bound) or bound <= 0.0:
                raise ValueError(
                    f"thermal elastic collision majorant is zero for {name}"
                )
            group_bounds.append(float(bound))
            group_masses.append(float(target_mass))
            group_stds.append(float(component_std))
            group_mean_speeds.append(float(mean_speed))
            group_modes.append(int(mode))
            group_sigma_maxima.append(float(sigma_max))

        total_bound = float(cold_bound + sum(group_bounds))
        if not math.isfinite(total_bound) or total_bound <= 0.0:
            raise ValueError(
                "internal monte_carlo could not build a positive collision majorant"
            )
        return cls(
            table=table,
            process_collision_group=process_group,
            cold_bound_m3_s=float(cold_bound),
            thermal_group_bounds_m3_s=np.asarray(group_bounds, dtype=float),
            thermal_group_target_mass_amu=np.asarray(group_masses, dtype=float),
            thermal_group_component_std_m_s=np.asarray(group_stds, dtype=float),
            thermal_group_mean_speed_m_s=np.asarray(group_mean_speeds, dtype=float),
            thermal_group_proposal_mode=np.asarray(group_modes, dtype=np.int64),
            thermal_group_sigma_max_m2=np.asarray(group_sigma_maxima, dtype=float),
            max_electron_speed_m_s=max_speed,
            max_energy_eV=float(max_energy_eV),
            trial_rate_coefficient_m3_s=float(
                np.nextafter(total_bound * _TRIAL_BOUND_SAFETY, math.inf)
            ),
        )

    def evaluate(self, energy_eV: float) -> tuple[np.ndarray, float]:
        values = self.table.evaluate(float(energy_eV))
        return values, float(np.sum(values))

    def trial_collision_frequency_s_inv(self, density_m3: float) -> float:
        density = float(density_m3)
        if not math.isfinite(density) or density <= 0.0:
            raise ValueError("collision majorant density must be positive")
        return density * self.trial_rate_coefficient_m3_s

    def sample_trial(
        self,
        electron_velocity_m_s: np.ndarray,
        rng: np.random.Generator,
    ) -> SampledCollisionTrial:
        velocity = np.asarray(electron_velocity_m_s, dtype=float)
        if velocity.shape != (3,) or np.any(~np.isfinite(velocity)):
            raise ValueError("collision trial electron velocity is invalid")
        speed = float(np.linalg.norm(velocity))
        if speed > self.max_electron_speed_m_s * (1.0 + 1.0e-12):
            raise RuntimeError(
                "internal_monte_carlo particle exceeded the collision-majorant "
                "energy support"
            )
        mark = float(rng.random()) * self.trial_rate_coefficient_m3_s
        group = -2
        group_bound = 0.0
        if mark < self.cold_bound_m3_s:
            group = -1
            group_bound = self.cold_bound_m3_s
        else:
            mark -= self.cold_bound_m3_s
            for index, bound in enumerate(self.thermal_group_bounds_m3_s):
                if mark < float(bound):
                    group = int(index)
                    group_bound = float(bound)
                    break
                mark -= float(bound)
        if group == -2:
            return SampledCollisionTrial(None, 0.0)

        target_velocity: np.ndarray | None = None
        relative_energy: float | None = None
        if group == -1:
            energy = float(energy_from_speed_eV(speed))
            values = self.table.evaluate(energy)
            mask = self.process_collision_group == -1
            total_sigma = float(np.sum(values[mask]))
            physical_rate_coefficient = speed * total_sigma
        else:
            component_std = float(self.thermal_group_component_std_m_s[group])
            mean_speed = float(self.thermal_group_mean_speed_m_s[group])
            mode = int(self.thermal_group_proposal_mode[group])
            if mode == THERMAL_PROPOSAL_MAXWELL:
                target_velocity = _sample_maxwellian_target_velocity(
                    component_std,
                    rng,
                )
            else:
                mixture_denominator = speed + mean_speed
                if float(rng.random()) < speed / mixture_denominator:
                    target_velocity = _sample_maxwellian_target_velocity(
                        component_std,
                        rng,
                    )
                else:
                    target_velocity = _sample_speed_weighted_maxwellian_target_velocity(
                        component_std,
                        rng,
                    )
            target_speed = float(np.linalg.norm(target_velocity))
            relative_speed = float(np.linalg.norm(velocity - target_velocity))
            relative_energy = float(energy_from_speed_eV(relative_speed))
            values = self.table.evaluate(relative_energy)
            mask = self.process_collision_group == group
            total_sigma = float(np.sum(values[mask]))
            if mode == THERMAL_PROPOSAL_MAXWELL:
                physical_rate_coefficient = relative_speed * total_sigma
            else:
                sigma_max = float(self.thermal_group_sigma_max_m2[group])
                denominator = max(speed + target_speed, np.finfo(float).tiny)
                physical_rate_coefficient = (
                    relative_speed * total_sigma * (speed + mean_speed) / denominator
                )
                group_bound = sigma_max * (self.max_electron_speed_m_s + mean_speed)

        probability = physical_rate_coefficient / max(group_bound, np.finfo(float).tiny)
        if probability > 1.0 + 2.0e-12:
            raise RuntimeError(
                "internal_monte_carlo thermal collision majorant was exceeded; "
                f"ratio={probability:.17g}"
            )
        probability = float(np.clip(probability, 0.0, 1.0))
        if float(rng.random()) > probability or total_sigma <= 0.0:
            return SampledCollisionTrial(None, probability)

        threshold = float(rng.random()) * total_sigma
        cumulative = 0.0
        choice: int | None = None
        for index in np.flatnonzero(mask):
            cumulative += float(values[index])
            if threshold < cumulative:
                choice = int(index)
                break
        if choice is None:
            choice = int(np.flatnonzero(mask)[-1])
        return SampledCollisionTrial(
            choice,
            probability,
            target_velocity_m_s=target_velocity,
            relative_energy_eV=relative_energy,
        )


def _project_processes(
    config: SwarmConfig, cross_sections: ActiveMixtureInputs
) -> list[_ProjectedProcess]:
    resolve_particle_scattering(
        cross_sections,
        angular_model=config.physics.angular_scattering.model,
        active_species=cross_sections.active_species,
    )
    fractions = {
        component.species: float(component.fraction)
        for component in cross_sections.components
    }
    out = [
        _ProjectedProcess(process, fractions[process.species])
        for process in cross_sections.processes
    ]
    if not out:
        raise ValueError("internal monte_carlo requires at least one active process")
    return out


def _species_has_total_scattering(
    processes: list[_ProjectedProcess],
    species: str,
) -> bool:
    return any(
        item.process.species == species
        and item.process.scattering_role in TOTAL_SCATTERING_ROLES
        for item in processes
    )


def _is_scattering_event_process(
    item: _ProjectedProcess,
    processes: list[_ProjectedProcess],
    angular_model_name: str,
) -> bool:
    role = item.process.scattering_role
    if role not in SCATTERING_ROLES:
        return False
    if angular_model_name == "maxent_p1":
        return role in TOTAL_SCATTERING_ROLES
    if _species_has_total_scattering(processes, item.process.species):
        return role in TOTAL_SCATTERING_ROLES
    return role in MOMENTUM_TRANSFER_ROLES


def _moment_cross_sections(
    processes: list[_ProjectedProcess],
    species: str,
    energy_eV: float,
) -> tuple[float, float]:
    energy = np.array([float(energy_eV)], dtype=float)
    total = 0.0
    momentum = 0.0
    for item in processes:
        if item.process.species != species:
            continue
        sigma = float(item.process.sigma(energy)[0])
        if item.process.scattering_role in TOTAL_SCATTERING_ROLES:
            total += sigma
        elif item.process.scattering_role in MOMENTUM_TRANSFER_ROLES:
            momentum += sigma
    if not np.isfinite(total) or total <= 0.0:
        raise ValueError(
            f"maxent_p1 total scattering cross section is zero for {species}"
        )
    if not np.isfinite(momentum) or momentum < 0.0:
        raise ValueError(
            f"maxent_p1 momentum-transfer cross section is invalid for {species}"
        )
    return total, momentum


def _max_cross_section_energy(processes: list[_ProjectedProcess]) -> float:
    max_energy = max(float(np.nanmax(item.process.energy_eV)) for item in processes)
    if not np.isfinite(max_energy) or max_energy <= 0.0:
        raise ValueError("internal monte_carlo requires finite cross-section energies")
    return max_energy


def _thermal_elastic_velocity_after_collision(
    config: SwarmConfig,
    process: CrossSectionProcess,
    electron_velocity_m_s: np.ndarray,
    target_velocity_m_s: np.ndarray,
    mu: float,
    rng: np.random.Generator,
) -> np.ndarray:
    mass_amu = process.mass_amu or gas_mass_amu(config.conditions, process.species)
    return elastic_binary_collision_velocity_m_s(
        electron_velocity_m_s,
        target_velocity_m_s,
        mu,
        float(rng.uniform(0.0, 2.0 * math.pi)),
        mass_amu,
    )


def _energy_loss_eV(process_type: ProcessType, threshold_eV: float | None) -> float:
    if process_type in {ProcessType.EXCITATION, ProcessType.IONIZATION}:
        return float(threshold_eV or 0.0)
    if process_type == ProcessType.SUPERELASTIC:
        return -abs(float(threshold_eV or 0.0))
    return 0.0


def _ionization_daughters(
    config: SwarmConfig,
    threshold_eV: float | None,
    energy_eV: float,
) -> IonizationDaughters:
    return ionization_daughters(
        energy_eV,
        float(threshold_eV or 0.0),
        model=config.physics.ionization.energy_sharing,
        secondary_electron_energy_eV=(
            config.physics.ionization.secondary_electron_energy_eV
        ),
    )


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

    daughters = _ionization_daughters(config, threshold_eV, energy_eV)
    energies = daughters.energies_eV
    if len(energies) == 1:
        tracked = energies[0]
        untracked = 0.0
    elif rng.random() < 0.5:
        tracked, untracked = energies[0], energies[1]
    else:
        tracked, untracked = energies[1], energies[0]
    return _PostReactionOutcome(
        tracked_energy_eV=max(tracked, 1.0e-4),
        ionization_threshold_loss_eV=daughters.threshold_loss_eV,
        ionization_untracked_secondary_energy_eV=max(untracked, 0.0),
    )
