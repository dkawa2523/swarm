"""Canonical Argon cross-section evidence used by the GEC-ICP adapter."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from hashlib import sha256
import math
from pathlib import Path

from .contracts import GecIcpContractError


_REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
CANONICAL_ARGON_CROSS_SECTIONS = (
    _REPOSITORY_ROOT / "examples" / "cross_sections" / "argon_application_library.csv"
)
COMSOL_ARGON_IMPORT_REFERENCE = _REPOSITORY_ROOT / "Model" / "Ar_xsecs.txt"

_GROUND_SIGNATURES = {
    "elastic": ("Ar", "e+Ar=>e+Ar", "elastic", 0.0),
    "excitation": ("Ar", "e+Ar=>e+Ars", "excitation", 11.5),
    "ionization": ("Ar", "e+Ar=>2e+Ar+", "ionization", 15.8),
}
_COLUMNS = (
    "species",
    "process",
    "type",
    "threshold_eV",
    "mass_amu",
    "energy_eV",
    "cross_section_m2",
)


@dataclass(frozen=True, slots=True)
class CanonicalGroundCrossSection:
    species: str
    process: str
    process_type: str
    threshold_eV: float
    energy_eV: tuple[float, ...]
    cross_section_m2: tuple[float, ...]


def canonical_cross_sections_combined_sha256() -> str:
    """Return the exporter-compatible combined hash for the canonical CSV."""

    path = CANONICAL_ARGON_CROSS_SECTIONS
    if not path.is_file():
        raise GecIcpContractError(
            f"canonical Argon cross-section CSV is missing: {path}"
        )
    file_digest = sha256(path.read_bytes()).hexdigest()
    digest = sha256()
    digest.update(str(path.resolve()).encode("utf-8"))
    digest.update(b"\0")
    digest.update(file_digest.encode("ascii"))
    digest.update(b"\0")
    return digest.hexdigest()


def load_canonical_ground_cross_sections() -> dict[str, CanonicalGroundCrossSection]:
    """Load the three ground-state Argon processes used by the ICP MPH."""

    path = CANONICAL_ARGON_CROSS_SECTIONS
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            if tuple(reader.fieldnames or ()) != _COLUMNS:
                raise GecIcpContractError(
                    "canonical Argon cross-section CSV has an unexpected schema"
                )
            rows = list(reader)
    except OSError as exc:
        raise GecIcpContractError(
            f"cannot read canonical Argon cross-section CSV: {path}"
        ) from exc

    grouped: dict[str, list[tuple[float, float]]] = {
        role: [] for role in _GROUND_SIGNATURES
    }
    for row_number, row in enumerate(rows, start=2):
        try:
            threshold = float(row["threshold_eV"])
            mass = float(row["mass_amu"])
            energy = float(row["energy_eV"])
            cross_section = float(row["cross_section_m2"])
        except (KeyError, TypeError, ValueError) as exc:
            raise GecIcpContractError(
                f"canonical Argon cross-section CSV row {row_number} is malformed"
            ) from exc
        signature = (
            row.get("species"),
            row.get("process"),
            row.get("type"),
            threshold,
        )
        role = next(
            (
                candidate
                for candidate, expected in _GROUND_SIGNATURES.items()
                if signature[:3] == expected[:3]
                and math.isclose(threshold, expected[3], rel_tol=0.0, abs_tol=1.0e-12)
            ),
            None,
        )
        if role is None:
            raise GecIcpContractError(
                "canonical Argon cross-section CSV contains an unexpected "
                f"process at row {row_number}: {signature!r}"
            )
        if not math.isclose(mass, 39.948, rel_tol=0.0, abs_tol=1.0e-12):
            raise GecIcpContractError(
                f"canonical Argon cross-section CSV row {row_number} has mass {mass}"
            )
        if (
            not math.isfinite(energy)
            or not math.isfinite(cross_section)
            or energy < 0.0
            or cross_section < 0.0
        ):
            raise GecIcpContractError(
                f"canonical Argon cross-section CSV row {row_number} is nonphysical"
            )
        grouped[role].append((energy, cross_section))

    result: dict[str, CanonicalGroundCrossSection] = {}
    for role, signature in _GROUND_SIGNATURES.items():
        values = grouped[role]
        if len(values) < 2 or any(
            right[0] <= left[0] for left, right in zip(values, values[1:])
        ):
            raise GecIcpContractError(
                f"canonical Argon {role} cross-section grid is missing or unordered"
            )
        if not any(value > 0.0 for _, value in values):
            raise GecIcpContractError(
                f"canonical Argon {role} cross sections are all zero"
            )
        result[role] = CanonicalGroundCrossSection(
            species=signature[0],
            process=signature[1],
            process_type=signature[2],
            threshold_eV=signature[3],
            energy_eV=tuple(energy for energy, _ in values),
            cross_section_m2=tuple(value for _, value in values),
        )
    return result


__all__ = [
    "CANONICAL_ARGON_CROSS_SECTIONS",
    "COMSOL_ARGON_IMPORT_REFERENCE",
    "CanonicalGroundCrossSection",
    "canonical_cross_sections_combined_sha256",
    "load_canonical_ground_cross_sections",
]
