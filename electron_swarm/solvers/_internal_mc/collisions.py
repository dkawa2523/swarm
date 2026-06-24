"""Collision sampling and reaction helpers for the internal Monte Carlo backend."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.constants import AMU_KG, ELECTRON_MASS_KG, EV_TO_J
from electron_swarm.core.cross_sections import (
    CrossSectionProcess,
    CrossSectionSet,
    ProcessType,
    gas_mass_amu,
    mixture_fraction,
)
from electron_swarm.solvers._internal_mc.orbit import (
    _random_direction,
    _speed_from_energy,
)


SCATTERING_TYPES = {
    ProcessType.MOMENTUM,
    ProcessType.ELASTIC,
    ProcessType.EFFECTIVE,
}
TOTAL_SCATTERING_TYPES = {
    ProcessType.ELASTIC,
    ProcessType.EFFECTIVE,
}
MOMENTUM_TRANSFER_TYPES = {
    ProcessType.MOMENTUM,
}
TRAJECTORY_REACTION_TYPES = {
    ProcessType.EXCITATION,
    ProcessType.IONIZATION,
    ProcessType.SUPERELASTIC,
}

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
class _IonizationDaughters:
    threshold_loss_eV: float
    primary_energy_eV: float
    secondary_energy_eV: float | None

    @property
    def energies_eV(self) -> tuple[float, ...]:
        if self.secondary_energy_eV is None:
            return (self.primary_energy_eV,)
        return (self.primary_energy_eV, self.secondary_energy_eV)

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


def _species_has_total_scattering(
    processes: list[_ProjectedProcess],
    species: str,
) -> bool:
    return any(
        item.process.species == species
        and item.process.process_type in TOTAL_SCATTERING_TYPES
        for item in processes
    )


def _is_scattering_event_process(
    item: _ProjectedProcess,
    processes: list[_ProjectedProcess],
    angular_model_name: str,
) -> bool:
    ptype = item.process.process_type
    if ptype not in SCATTERING_TYPES:
        return False
    if angular_model_name == "maxent_p1":
        return ptype in TOTAL_SCATTERING_TYPES
    if _species_has_total_scattering(processes, item.process.species):
        return ptype in TOTAL_SCATTERING_TYPES
    return True


def _active_scattering_species(processes: list[_ProjectedProcess]) -> set[str]:
    return {
        item.process.species
        for item in processes
        if item.process.process_type in SCATTERING_TYPES
    }


def _validate_maxent_p1_cross_sections(processes: list[_ProjectedProcess]) -> None:
    missing: list[str] = []
    for species in sorted(_active_scattering_species(processes)):
        has_total = any(
            item.process.species == species
            and item.process.process_type in TOTAL_SCATTERING_TYPES
            for item in processes
        )
        has_momentum = any(
            item.process.species == species
            and item.process.process_type in MOMENTUM_TRANSFER_TYPES
            for item in processes
        )
        if not (has_total and has_momentum):
            missing.append(species)
    if missing:
        raise NotImplementedError(
            "internal monte_carlo maxent_p1 requires both total/effective "
            "elastic and momentum-transfer cross sections for active species "
            f"{missing}"
        )


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
        if item.process.process_type in TOTAL_SCATTERING_TYPES:
            total += sigma
        elif item.process.process_type in MOMENTUM_TRANSFER_TYPES:
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


def _sigma_arrays(
    processes: list[_ProjectedProcess],
    energy_eV: float,
    *,
    angular_model_name: str = "isotropic",
) -> tuple[np.ndarray, float]:
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
            if (
                _is_scattering_event_process(item, processes, angular_model_name)
                or item.process.process_type in TRAJECTORY_REACTION_TYPES
            )
            else 0.0
            for value, item in zip(raw_values, processes, strict=True)
        ],
        dtype=float,
    )
    collision_total = float(np.sum(trajectory_values))
    return trajectory_values, collision_total

def _trial_collision_frequency(
    processes: list[_ProjectedProcess],
    density_m3: float,
    *,
    angular_model_name: str = "isotropic",
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
        if (
            _is_scattering_event_process(item, processes, angular_model_name)
            or item.process.process_type in TRAJECTORY_REACTION_TYPES
        ):
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

def _ionization_daughters(
    config: SwarmConfig,
    threshold_eV: float | None,
    energy_eV: float,
) -> _IonizationDaughters:
    available = max(float(energy_eV) - float(threshold_eV or 0.0), 0.0)
    threshold_loss = float(threshold_eV or 0.0)
    sharing = config.physics.ionization.energy_sharing
    if sharing == "equal":
        daughter = 0.5 * available
        return _IonizationDaughters(
            threshold_loss_eV=threshold_loss,
            primary_energy_eV=daughter,
            secondary_energy_eV=daughter,
        )
    if sharing == "primary_secondary":
        secondary = min(
            max(float(config.physics.ionization.secondary_electron_energy_eV), 0.0),
            available,
        )
        primary = max(available - secondary, 0.0)
        return _IonizationDaughters(
            threshold_loss_eV=threshold_loss,
            primary_energy_eV=primary,
            secondary_energy_eV=secondary,
        )
    if sharing == "loss_only":
        return _IonizationDaughters(
            threshold_loss_eV=threshold_loss,
            primary_energy_eV=available,
            secondary_energy_eV=None,
        )
    raise ValueError(f"unsupported ionization energy sharing model: {sharing!r}")

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
