"""State-resolved collision preprocessing.

This module keeps state-population handling outside individual solvers.  It can
augment a loaded :class:`CrossSectionSet` with generated superelastic processes
from optional YAML sections such as::

    state_resolved:
      enabled: true
      generate_superelastic: true
      transitions:
        - species: Ar
          process: excitation_4s
          initial_state: ground
          final_state: Ar_4s
    species_states:
      Ar:
        ground: {energy_eV: 0.0, population: 1.0}
        Ar_4s: {energy_eV: 11.55, population: 1.0e-4}

The generated processes are ordinary ``CrossSectionProcess`` objects with
``ProcessType.SUPERELASTIC`` so existing two-term, multi-term, and MC paths can
consume them without new solver contracts.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from electron_swarm.collisions.config import raw_config_from_source
from electron_swarm.core.cross_sections import (
    CrossSectionProcess,
    CrossSectionSet,
    ProcessType,
)


@dataclass(frozen=True, slots=True)
class SpeciesState:
    species: str
    name: str
    energy_eV: float
    population: float
    degeneracy: float | None = None


def _state_items(species: str, raw_states: object) -> list[SpeciesState]:
    states: list[SpeciesState] = []
    if isinstance(raw_states, dict):
        iterable = [dict(value or {}, name=name) for name, value in raw_states.items()]
    elif isinstance(raw_states, list):
        iterable = [item for item in raw_states if isinstance(item, dict)]
    else:
        iterable = []
    for item in iterable:
        name = str(item.get("name", "state"))
        states.append(
            SpeciesState(
                species=species,
                name=name,
                energy_eV=float(item.get("energy_eV", item.get("energy_ev", 0.0))),
                population=float(item.get("population", 0.0)),
                degeneracy=(
                    float(item["degeneracy"])
                    if item.get("degeneracy") is not None
                    else None
                ),
            )
        )
    return states


def parse_species_states(raw: dict[str, Any]) -> dict[tuple[str, str], SpeciesState]:
    """Parse explicit and molecular convenience state sections."""

    parsed: dict[tuple[str, str], SpeciesState] = {}
    for species, raw_states in (raw.get("species_states", {}) or {}).items():
        for state in _state_items(str(species), raw_states):
            parsed[(state.species, state.name)] = state

    # ``molecular_states`` is a user-friendly shorthand.  It expands to the
    # same internal SpeciesState map so solvers only see one state concept.
    for species, mol in (raw.get("molecular_states", {}) or {}).items():
        if not isinstance(mol, dict):
            continue
        vib = mol.get("vibrational", {}) or {}
        if isinstance(vib, dict):
            for item in vib.get("levels", []) or []:
                if not isinstance(item, dict):
                    continue
                label = item.get("name", f"v{item.get('v', len(parsed))}")
                state = SpeciesState(
                    species=str(species),
                    name=str(label),
                    energy_eV=float(item.get("energy_eV", item.get("energy_ev", 0.0))),
                    population=float(item.get("population", 0.0)),
                    degeneracy=(
                        float(item["degeneracy"])
                        if item.get("degeneracy") is not None
                        else None
                    ),
                )
                parsed[(state.species, state.name)] = state
    return parsed


def _transition_matches(process: CrossSectionProcess, transition: dict[str, Any]) -> bool:
    species = transition.get("species")
    if species is not None and str(species) != process.species:
        return False
    token = transition.get("process")
    if token is None:
        return process.process_type == ProcessType.EXCITATION
    token_s = str(token).lower()
    return token_s == process.process.lower() or token_s in process.process.lower()


def _superelastic_from_excitation(
    process: CrossSectionProcess,
    lower: SpeciesState,
    upper: SpeciesState,
    *,
    detailed_balance: str,
    energy_floor_eV: float | None,
) -> CrossSectionProcess | None:
    delta = float(upper.energy_eV - lower.energy_eV)
    if delta <= 0.0 or upper.population <= 0.0:
        return None
    energy = np.asarray(process.energy_eV, dtype=float)
    shifted = energy + delta
    sigma_exc = process.sigma(shifted, left=0.0, right=0.0)
    floor = (
        float(energy_floor_eV)
        if energy_floor_eV is not None
        else max(1.0e-3, 0.01 * delta)
    )
    degeneracy_ratio = 1.0
    if detailed_balance == "degeneracy" and lower.degeneracy and upper.degeneracy:
        degeneracy_ratio = float(lower.degeneracy) / float(upper.degeneracy)
    population_ratio = upper.population / max(lower.population, 1.0e-300)
    sigma = (
        sigma_exc
        * (shifted / np.maximum(energy, max(floor, 1.0e-12)))
        * degeneracy_ratio
        * population_ratio
    )
    sigma = np.nan_to_num(np.clip(sigma, 0.0, None), nan=0.0, posinf=0.0, neginf=0.0)
    if not np.any(sigma > 0.0):
        return None
    metadata = dict(process.metadata)
    metadata.update(
        {
            "generated_by": "state_resolved_superelastic",
            "source_process": process.process,
            "initial_state": upper.name,
            "final_state": lower.name,
            "state_energy_gap_eV": delta,
            "population_ratio": population_ratio,
            "detailed_balance": detailed_balance,
            "detailed_balance_energy_floor_eV": floor,
            "regularized_low_energy_singularity": True,
        }
    )
    return CrossSectionProcess(
        species=process.species,
        process=f"superelastic:{upper.name}->{lower.name}:{process.process}",
        process_type=ProcessType.SUPERELASTIC,
        energy_eV=energy,
        cross_section_m2=sigma,
        threshold_eV=delta,
        mass_amu=process.mass_amu,
        metadata=metadata,
    )


def generate_superelastic_processes(
    cross_sections: CrossSectionSet,
    raw: dict[str, Any],
) -> list[CrossSectionProcess]:
    state_cfg = raw.get("state_resolved", {}) or {}
    if not bool(state_cfg.get("enabled", False)):
        return []
    if not bool(state_cfg.get("generate_superelastic", True)):
        return []
    detailed_balance = str(state_cfg.get("detailed_balance", "simple"))
    floor_raw = state_cfg.get("superelastic_energy_floor_eV")
    energy_floor = float(floor_raw) if floor_raw is not None else None
    states = parse_species_states(raw)
    generated: list[CrossSectionProcess] = []
    for transition in state_cfg.get("transitions", []) or []:
        if not isinstance(transition, dict):
            continue
        species = str(transition.get("species", ""))
        lower_name = str(transition.get("initial_state", ""))
        upper_name = str(transition.get("final_state", ""))
        lower = states.get((species, lower_name))
        upper = states.get((species, upper_name))
        if lower is None or upper is None:
            continue
        for process in cross_sections.processes:
            if process.process_type != ProcessType.EXCITATION:
                continue
            if not _transition_matches(process, transition):
                continue
            new_process = _superelastic_from_excitation(
                process,
                lower,
                upper,
                detailed_balance=detailed_balance,
                energy_floor_eV=energy_floor,
            )
            if new_process is not None:
                generated.append(new_process)
    return generated


def augment_cross_sections_from_config(
    cross_sections: CrossSectionSet,
    config: object,
) -> CrossSectionSet:
    """Return a cross-section set augmented by optional state-resolved physics."""

    raw = raw_config_from_source(config)
    generated = generate_superelastic_processes(cross_sections, raw)
    if not generated:
        return cross_sections
    return CrossSectionSet([*cross_sections.processes, *generated])
