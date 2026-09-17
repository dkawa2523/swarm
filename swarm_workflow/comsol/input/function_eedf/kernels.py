"""Physical collision kernels used by adaptive Function-EEDF projection."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from electron_swarm.core.constants import E_CHARGE_C, ELECTRON_MASS_KG
from swarm_workflow.tables.contracts import (
    COLLISION_RATE_KERNEL_COLUMNS,
)

from .contracts import CollisionRateKernel, FunctionEedfError


_SPEED_PER_SQRT_EV = math.sqrt(2.0 * E_CHARGE_C / ELECTRON_MASS_KG)


def read_collision_rate_kernels(path: str | Path) -> tuple[CollisionRateKernel, ...]:
    """Read provenance-bound cross sections without depending on a solver."""

    grouped: dict[tuple[str, str, str, float, str], list[tuple[float, float]]] = {}
    try:
        with Path(path).open("r", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            if tuple(reader.fieldnames or ()) != COLLISION_RATE_KERNEL_COLUMNS:
                raise FunctionEedfError(
                    "collision-rate kernel columns are not canonical"
                )
            for row in reader:
                threshold = float(row["threshold_eV"])
                energy = float(row["electron_energy_eV"])
                cross_section = float(row["cross_section_m2"])
                extrapolation = row["high_energy_extrapolation"]
                if (
                    extrapolation not in {"zero", "hold", "error"}
                    or not math.isfinite(threshold)
                    or threshold < 0.0
                    or not math.isfinite(energy)
                    or energy < 0.0
                    or not math.isfinite(cross_section)
                    or cross_section < 0.0
                ):
                    raise FunctionEedfError(
                        "collision-rate kernel contains invalid data"
                    )
                key = (
                    row["species"],
                    row["process"],
                    row["process_type"],
                    threshold,
                    extrapolation,
                )
                grouped.setdefault(key, []).append((energy, cross_section))
    except (OSError, KeyError, ValueError) as exc:
        raise FunctionEedfError(f"invalid collision-rate kernel table: {exc}") from exc

    kernels: list[CollisionRateKernel] = []
    for key, points in grouped.items():
        ordered = sorted(points)
        energy = np.asarray([item[0] for item in ordered], dtype=float)
        sigma = np.asarray([item[1] for item in ordered], dtype=float)
        if len(energy) < 2 or np.any(np.diff(energy) <= 0.0):
            raise FunctionEedfError("collision-rate kernel grid must increase")
        species, process, process_type, threshold, extrapolation = key
        kernels.append(
            CollisionRateKernel(
                species=species,
                process=process,
                process_type=process_type,
                threshold_eV=threshold,
                electron_energies_eV=energy,
                cross_sections_m2=sigma,
                high_energy_extrapolation=extrapolation,
            )
        )
    return tuple(sorted(kernels, key=lambda item: item.identity))


def collision_rate_coefficient(
    energies_eV: np.ndarray,
    f0_eV_m32: np.ndarray,
    kernel: CollisionRateKernel,
) -> float:
    """Integrate ``<sigma v>`` for an EEPF normalized with ``sqrt(E)``."""

    energies = np.asarray(energies_eV, dtype=float)
    f0 = np.asarray(f0_eV_m32, dtype=float)
    if kernel.high_energy_extrapolation == "error" and float(energies[-1]) > float(
        kernel.electron_energies_eV[-1]
    ):
        raise FunctionEedfError(
            f"EEDF support exceeds cross-section support for {kernel.identity}"
        )
    right = (
        0.0
        if kernel.high_energy_extrapolation == "zero"
        else float(kernel.cross_sections_m2[-1])
    )
    sigma = np.interp(
        energies,
        kernel.electron_energies_eV,
        kernel.cross_sections_m2,
        left=0.0,
        right=right,
    )
    integrand = energies * sigma * f0
    return _SPEED_PER_SQRT_EV * float(np.trapezoid(integrand, energies))
