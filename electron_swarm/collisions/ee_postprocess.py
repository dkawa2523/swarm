"""Optional electron-electron relaxation post-processing.

This is the first, merge-safe e-e collision hook.  It deliberately avoids a full
Fokker-Planck operator and does not touch Monte Carlo cases.  When enabled in
YAML it relaxes Boltzmann-family EEDFs slightly toward a Maxwellian target that
preserves the mean energy unless configured otherwise.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Iterable

import numpy as np
import yaml

from electron_swarm.collisions.electron_electron import mean_energy_eV, relaxation_target
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.diagnostics.common import widths_from_centers


_BOLTZMANN_SOLVERS = {"boltzmann_two_term", "multiterm_boltzmann"}
_EPS = 1.0e-300


def _raw_config_from_source(config: object) -> dict[str, Any]:
    path = getattr(config, "source_path", None)
    if path is None:
        return {}
    try:
        data = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    except OSError:
        return {}
    return data if isinstance(data, dict) else {}


def _ee_section(raw: dict[str, Any]) -> dict[str, Any]:
    section = raw.get("electron_electron", raw.get("electron_electron_collision", {}))
    return section if isinstance(section, dict) else {}


def _normalise(energy_eV: np.ndarray, eedf: np.ndarray) -> np.ndarray:
    widths = widths_from_centers(energy_eV)
    total = float(np.sum(np.clip(eedf, 0.0, None) * widths))
    if total <= 0.0 or not np.isfinite(total):
        raise ValueError("Cannot normalize e-e relaxed EEDF")
    return np.clip(eedf, 0.0, None) / total


def _update_eepf(case: SwarmCaseResult) -> None:
    energy = np.asarray(case.energy_eV, dtype=float)
    denom = np.sqrt(np.maximum(energy, 1.0e-30))
    try:
        case.eepf = np.asarray(case.eedf, dtype=float) / denom
    except Exception:
        # Some result implementations may not expose mutable eepf.  The EEDF
        # and metadata are still updated; diagnostics will flag the case.
        case.metadata["ee_eepf_recomputed"] = False
    else:
        case.metadata["ee_eepf_recomputed"] = True


def _apply_to_case(case: SwarmCaseResult, section: dict[str, Any]) -> SwarmCaseResult:
    if case.solver not in _BOLTZMANN_SOLVERS:
        return case
    energy = np.asarray(case.energy_eV, dtype=float)
    old = np.asarray(case.eedf, dtype=float)
    widths = widths_from_centers(energy)
    if len(energy) == 0 or len(old) != len(energy):
        case.metadata["ee_relaxation_skipped"] = "invalid_eedf_shape"
        return case

    alpha = section.get("relaxation_fraction", section.get("strength_scale", 0.05))
    alpha = float(np.clip(float(alpha), 0.0, 1.0))
    conserve = bool(section.get("conserve_mean_energy", True))
    fallback_te = float(section.get("fallback_temperature_eV", section.get("temperature_eV", 2.0)))
    before = mean_energy_eV(energy, widths, old)
    target = relaxation_target(
        energy,
        widths,
        old,
        conserve_mean_energy=conserve,
        fallback_temperature_eV=fallback_te,
    )
    new = _normalise(energy, (1.0 - alpha) * old + alpha * target)
    after = mean_energy_eV(energy, widths, new)

    case.eedf = new
    _update_eepf(case)
    case.metadata.update(
        {
            "ee_enabled": True,
            "ee_model": "relaxation_postprocess",
            "ee_relaxation_fraction": alpha,
            "ee_conserve_mean_energy": conserve,
            "ee_mean_energy_before_eV": before,
            "ee_mean_energy_after_eV": after,
            "ee_rates_recomputed": False,
            "ee_transport_recomputed": False,
        }
    )
    return case


def apply_electron_electron_relaxation_from_config(
    cases: Iterable[SwarmCaseResult],
    config: object,
) -> list[SwarmCaseResult]:
    """Apply optional e-e relaxation to Boltzmann-family cases only."""

    raw = _raw_config_from_source(config)
    section = _ee_section(raw)
    if not bool(section.get("enabled", False)):
        return list(cases)
    if str(section.get("model", "relaxation")).lower() not in {"relaxation", "relaxation_postprocess"}:
        return list(cases)
    return [_apply_to_case(case, section) for case in cases]
