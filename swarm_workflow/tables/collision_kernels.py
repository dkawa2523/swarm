"""Materialize active cross sections for solver-independent EEDF projection."""

from __future__ import annotations

from dataclasses import replace
from hashlib import sha256
import math
from pathlib import Path
from typing import Any

from electron_swarm.core.config_parser import load_config
from electron_swarm.core.config import GasComponent
from electron_swarm.core.cross_sections import (
    CrossSectionProcess,
    load_active_mixture_inputs,
)

from .contracts import COLLISION_RATE_KERNEL_COLUMNS, TableBuildError


def active_collision_kernel_rows(
    metadata: dict[str, str],
    rates: list[dict[str, Any]],
    *,
    mixture: list[dict[str, str | float]],
) -> tuple[tuple[str, ...], list[dict[str, Any]]] | None:
    """Return the exact configured cross sections used by active rate rows."""

    raw_path = metadata.get("base_config_path")
    expected_hash = metadata.get("cross_sections_sha256")
    if not raw_path or not expected_hash:
        return None
    config_path = Path(raw_path)
    if not config_path.is_file():
        return None
    config = load_config(config_path)
    file_hashes = {
        str(file_config.path): _file_sha256(Path(file_config.path))
        for file_config in config.cross_sections.files
    }
    if _combined_hash(file_hashes) != expected_hash:
        raise TableBuildError(
            "configured cross sections changed after the workflow calculation"
        )
    conditions = replace(
        config.conditions,
        gas_mixture=[
            GasComponent(
                species=str(row["species"]),
                fraction=float(row["fraction"]),
                mass_amu=float(row["mass_amu"]),
            )
            for row in mixture
        ],
    )
    cross_sections = load_active_mixture_inputs(
        config.cross_sections,
        conditions,
    )
    active = {
        (
            str(row["species"]),
            str(row["process"]),
            str(row["process_type"]),
            _threshold(row.get("threshold_eV")),
        )
        for row in rates
    }
    selected: list[CrossSectionProcess] = []
    for identity in sorted(active):
        matches = [
            process
            for process in cross_sections.processes
            if _process_identity(process) == identity
        ]
        if len(matches) != 1:
            raise TableBuildError(
                "active rate process does not map to one configured cross section: "
                + ":".join(str(value) for value in identity[:3])
            )
        selected.append(matches[0])
    rows = [
        {
            "species": process.species,
            "process": process.process,
            "process_type": process.process_type.value,
            "threshold_eV": _threshold(process.threshold_eV),
            "electron_energy_eV": float(energy),
            "cross_section_m2": float(cross_section),
            "high_energy_extrapolation": str(
                process.metadata.get("high_energy_extrapolation", "hold")
            ).lower(),
        }
        for process in selected
        for energy, cross_section in zip(
            process.energy_eV,
            process.cross_section_m2,
            strict=True,
        )
    ]
    return COLLISION_RATE_KERNEL_COLUMNS, rows


def _process_identity(process: CrossSectionProcess) -> tuple[str, str, str, float]:
    return (
        process.species,
        process.process,
        process.process_type.value,
        _threshold(process.threshold_eV),
    )


def _threshold(value: object) -> float:
    if value in (None, ""):
        return 0.0
    result = float(value)
    if not math.isfinite(result) or result < 0.0:
        raise TableBuildError("collision threshold is invalid")
    return result


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _combined_hash(values: dict[str, str]) -> str:
    digest = sha256()
    for key, value in sorted(values.items()):
        digest.update(key.encode("utf-8"))
        digest.update(b"\0")
        digest.update(value.encode("ascii"))
        digest.update(b"\0")
    return digest.hexdigest()
