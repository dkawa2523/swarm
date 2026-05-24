"""Cross-section input parsing and interpolation utilities.

Canonical CSV format
--------------------
Each row represents one cross-section value:

    species,process,type,threshold_eV,mass_amu,energy_eV,cross_section_m2

Supported types are case-insensitive aliases of
``momentum``, ``elastic``, ``effective``, ``excitation``, ``ionization`` and
``attachment``. The canonical representation keeps all cross sections as
functions of electron energy in eV and cross section in m^2.

For LXCat/BOLSIG-style ``txt|lxcat|bolsig|dat`` inputs, the loader first uses
``lxcat_data_parser`` and then falls back to a small local block parser.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
import numpy as np
import pandas as pd
from lxcat_data_parser import CrossSectionReadingError
from lxcat_data_parser import CrossSectionSet as LxcatCrossSectionSet
from lxcat_data_parser import CrossSectionTypes as CST

from .config import ConditionsConfig, CrossSectionsConfig


class ProcessType(str, Enum):
    MOMENTUM = "momentum"
    ELASTIC = "elastic"
    EFFECTIVE = "effective"
    EXCITATION = "excitation"
    IONIZATION = "ionization"
    ATTACHMENT = "attachment"
    SUPERELASTIC = "superelastic"
    UNKNOWN = "unknown"


_TYPE_ALIASES = {
    "momentum": ProcessType.MOMENTUM,
    "momentum_transfer": ProcessType.MOMENTUM,
    "mt": ProcessType.MOMENTUM,
    "elastic": ProcessType.ELASTIC,
    "effective": ProcessType.EFFECTIVE,
    "total_elastic": ProcessType.EFFECTIVE,
    "excitation": ProcessType.EXCITATION,
    "exc": ProcessType.EXCITATION,
    "ionization": ProcessType.IONIZATION,
    "ionisation": ProcessType.IONIZATION,
    "ion": ProcessType.IONIZATION,
    "attachment": ProcessType.ATTACHMENT,
    "attach": ProcessType.ATTACHMENT,
    "superelastic": ProcessType.SUPERELASTIC,
    "super_elastic": ProcessType.SUPERELASTIC,
    "deexcitation": ProcessType.SUPERELASTIC,
    "de_excitation": ProcessType.SUPERELASTIC,
}


def normalize_process_type(value: str | None) -> ProcessType:
    if value is None:
        return ProcessType.UNKNOWN
    key = str(value).strip().lower().replace(" ", "_").replace("-", "_")
    return _TYPE_ALIASES.get(key, ProcessType.UNKNOWN)


def _process_type_from_lxcat(value: object) -> ProcessType:
    if value == CST.EFFECTIVE:
        return ProcessType.EFFECTIVE
    if value == CST.ELASTIC:
        return ProcessType.ELASTIC
    if value == CST.EXCITATION:
        return ProcessType.EXCITATION
    if value == CST.IONIZATION:
        return ProcessType.IONIZATION
    if value == CST.ATTACHMENT:
        return ProcessType.ATTACHMENT
    return normalize_process_type(getattr(value, "name", None))


def _extract_process_name(info: object, default: str) -> str:
    if isinstance(info, dict):
        process = info.get("PROCESS")
        if process:
            return str(process)
        nested = info.get("info")
        if isinstance(nested, dict) and nested.get("PROCESS"):
            return str(nested["PROCESS"])
    return default


@dataclass(slots=True)
class CrossSectionProcess:
    species: str
    process: str
    process_type: ProcessType
    energy_eV: np.ndarray
    cross_section_m2: np.ndarray
    threshold_eV: float | None = None
    mass_amu: float | None = None
    metadata: dict[str, str | float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        energy = np.asarray(self.energy_eV, dtype=float)
        sigma = np.asarray(self.cross_section_m2, dtype=float)
        if energy.ndim != 1 or sigma.ndim != 1 or len(energy) != len(sigma):
            raise ValueError(
                f"Invalid cross-section arrays for {self.species}:{self.process}"
            )
        if len(energy) < 2:
            raise ValueError(
                f"At least two energy points are required for {self.species}:{self.process}"
            )
        order = np.argsort(energy)
        energy = energy[order]
        sigma = sigma[order]
        keep = np.isfinite(energy) & np.isfinite(sigma) & (energy >= 0.0)
        energy = energy[keep]
        sigma = np.clip(sigma[keep], 0.0, None)
        uniq = np.unique(energy)
        if len(uniq) != len(energy):
            sigma2 = np.zeros_like(uniq)
            for j, value in enumerate(uniq):
                sigma2[j] = sigma[energy == value][-1]
            energy, sigma = uniq, sigma2
        self.energy_eV = energy
        self.cross_section_m2 = sigma

    def sigma(
        self, energy_eV: np.ndarray, *, left: float = 0.0, right: float | None = None
    ) -> np.ndarray:
        """Linearly interpolate cross section on an arbitrary energy grid."""

        if right is None:
            policy = str(
                self.metadata.get("high_energy_extrapolation", "hold")
            ).lower()
            if policy == "zero":
                right = 0.0
            elif policy == "hold":
                right = float(self.cross_section_m2[-1])
            elif policy == "error":
                energy = np.asarray(energy_eV, dtype=float)
                if np.any(energy > self.energy_eV[-1]):
                    raise ValueError(
                        "Cross-section interpolation requested above the "
                        f"tabulated range for {self.species}:{self.process}; "
                        "increase the energy grid or set "
                        "cross_sections.high_energy_extrapolation."
                    )
                right = float(self.cross_section_m2[-1])
            else:
                raise ValueError(
                    "Unsupported high-energy extrapolation policy for "
                    f"{self.species}:{self.process}: {policy!r}"
                )
        return np.interp(
            energy_eV, self.energy_eV, self.cross_section_m2, left=left, right=right
        )


@dataclass(slots=True)
class CrossSectionSet:
    processes: list[CrossSectionProcess]

    def by_type(self, *types: ProcessType) -> list[CrossSectionProcess]:
        return [p for p in self.processes if p.process_type in types]

    def by_species(self, species: str) -> list[CrossSectionProcess]:
        return [p for p in self.processes if p.species == species]

    @property
    def species(self) -> list[str]:
        return sorted({p.species for p in self.processes})


def _read_csv_long(path: Path, default_species: str | None) -> list[CrossSectionProcess]:
    df = pd.read_csv(path)
    columns = {c.lower(): c for c in df.columns}
    required = {"energy_ev", "cross_section_m2"}
    if not required.issubset(columns):
        raise ValueError("not long CSV")
    species_col = columns.get("species")
    process_col = columns.get("process")
    type_col = columns.get("type", columns.get("process_type"))
    threshold_col = columns.get("threshold_ev")
    mass_col = columns.get("mass_amu")
    energy_col = columns["energy_ev"]
    sigma_col = columns["cross_section_m2"]

    if process_col is None and type_col is None:
        raise ValueError(f"{path}: long CSV requires process or type column")

    group_cols = [
        c
        for c in [species_col, process_col, type_col, threshold_col, mass_col]
        if c is not None
    ]
    processes = []
    grouped = df.groupby(group_cols, dropna=False) if group_cols else [((), df)]
    for keys, group in grouped:
        if not isinstance(keys, tuple):
            keys = (keys,)
        key_map = dict(zip(group_cols, keys, strict=False))
        species = str(key_map.get(species_col, default_species or "gas"))
        ptype = normalize_process_type(
            str(key_map.get(type_col, "unknown")) if type_col is not None else None
        )
        process_name = str(
            key_map.get(
                process_col,
                ptype.value if ptype != ProcessType.UNKNOWN else "process",
            )
        )
        threshold = key_map.get(threshold_col) if threshold_col is not None else None
        mass = key_map.get(mass_col) if mass_col is not None else None
        processes.append(
            CrossSectionProcess(
                species=species,
                process=process_name,
                process_type=ptype,
                threshold_eV=None if pd.isna(threshold) else float(threshold),
                mass_amu=None if pd.isna(mass) else float(mass),
                energy_eV=group[energy_col].to_numpy(dtype=float),
                cross_section_m2=group[sigma_col].to_numpy(dtype=float),
            )
        )
    return processes


def _infer_type_from_wide_column(name: str) -> ProcessType:
    key = name.lower().replace("_m2", "")
    for token, ptype in [
        ("momentum", ProcessType.MOMENTUM),
        ("effective", ProcessType.EFFECTIVE),
        ("elastic", ProcessType.ELASTIC),
        ("excitation", ProcessType.EXCITATION),
        ("exc", ProcessType.EXCITATION),
        ("ionization", ProcessType.IONIZATION),
        ("ionisation", ProcessType.IONIZATION),
        ("attachment", ProcessType.ATTACHMENT),
        ("superelastic", ProcessType.SUPERELASTIC),
        ("deexcitation", ProcessType.SUPERELASTIC),
    ]:
        if token in key:
            return ptype
    return ProcessType.UNKNOWN


def _threshold_from_column(name: str) -> float | None:
    import re

    match = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*e[vV]", name)
    return float(match.group(1)) if match else None


def _read_csv_wide(
    path: Path, default_species: str | None, default_mass_amu: float | None
) -> list[CrossSectionProcess]:
    df = pd.read_csv(path)
    columns = {c.lower(): c for c in df.columns}
    if "energy_ev" not in columns:
        raise ValueError(f"{path}: CSV must contain energy_eV")
    energy = df[columns["energy_ev"]].to_numpy(dtype=float)
    processes = []
    for col in df.columns:
        if col == columns["energy_ev"]:
            continue
        if not col.lower().endswith("_m2"):
            continue
        ptype = _infer_type_from_wide_column(col)
        processes.append(
            CrossSectionProcess(
                species=default_species or "gas",
                process=col.replace("_m2", ""),
                process_type=ptype,
                threshold_eV=_threshold_from_column(col),
                mass_amu=default_mass_amu,
                energy_eV=energy,
                cross_section_m2=df[col].to_numpy(dtype=float),
            )
        )
    if not processes:
        raise ValueError(f"{path}: no *_m2 cross-section columns found")
    return processes


def _default_species(
    file_species: str | None, conditions: ConditionsConfig
) -> str | None:
    if file_species:
        return file_species
    if len(conditions.gas_mixture) == 1:
        return conditions.gas_mixture[0].species
    return None


def _read_lxcat_text(path: Path, default_species: str | None) -> list[CrossSectionProcess]:
    if default_species is None:
        raise ValueError(
            f"{path}: species must be set for txt/lxcat/bolsig inputs when the gas mixture has multiple species"
        )
    parsed = LxcatCrossSectionSet(str(path), default_species, None)
    processes: list[CrossSectionProcess] = []
    for section in parsed.cross_sections:
        processes.append(
            CrossSectionProcess(
                species=str(section.species),
                process=_extract_process_name(
                    getattr(section, "info", None),
                    f"{getattr(section.type, 'name', 'unknown')} {section.species}",
                ),
                process_type=_process_type_from_lxcat(section.type),
                threshold_eV=float(section.threshold) if section.threshold is not None else None,
                mass_amu=None,
                energy_eV=section.data["energy"].to_numpy(dtype=float),
                cross_section_m2=section.data["cross section"].to_numpy(dtype=float),
            )
        )
    return processes


def _read_bolsig_like(path: Path, default_species: str | None) -> list[CrossSectionProcess]:
    """Fallback parser for simple BOLSIG+/LXCat blocks."""

    processes: list[CrossSectionProcess] = []
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        ptype = normalize_process_type(line)
        if ptype == ProcessType.UNKNOWN:
            i += 1
            continue
        i += 1
        while i < len(lines) and not lines[i].strip():
            i += 1
        signature = lines[i].strip() if i < len(lines) else ptype.value
        species = default_species or signature.split()[0].split("->")[0]
        threshold = None
        i += 1
        data: list[tuple[float, float]] = []
        while i < len(lines):
            s = lines[i].strip()
            if not s:
                i += 1
                if data:
                    break
                continue
            parts = s.replace(",", " ").split()
            if len(parts) >= 2:
                try:
                    e = float(parts[0])
                    sig = float(parts[1])
                    data.append((e, sig))
                    i += 1
                    continue
                except ValueError:
                    pass
            if "threshold" in s.lower():
                for token in parts:
                    try:
                        threshold = float(token)
                        break
                    except ValueError:
                        continue
            if data:
                break
            i += 1
        if len(data) >= 2:
            arr = np.asarray(data, dtype=float)
            processes.append(
                CrossSectionProcess(
                    species=species,
                    process=signature,
                    process_type=ptype,
                    threshold_eV=threshold,
                    energy_eV=arr[:, 0],
                    cross_section_m2=arr[:, 1],
                )
            )
    if not processes:
        raise ValueError(f"{path}: no BOLSIG/LXCat-like processes parsed")
    return processes


def load_cross_sections(
    config: CrossSectionsConfig, conditions: ConditionsConfig
) -> CrossSectionSet:
    """Load all configured cross-section files."""

    mass_by_species = {g.species: g.mass_amu for g in conditions.gas_mixture}
    processes: list[CrossSectionProcess] = []
    for file_cfg in config.files:
        path = Path(file_cfg.path)
        fmt = (file_cfg.format or config.format).lower()
        default_species = _default_species(file_cfg.species, conditions)
        default_mass = mass_by_species.get(default_species or "")
        if fmt == "csv":
            try:
                processes.extend(_read_csv_long(path, default_species))
            except ValueError:
                processes.extend(_read_csv_wide(path, default_species, default_mass))
        elif fmt in {"bolsig", "lxcat", "txt", "dat"}:
            try:
                processes.extend(_read_lxcat_text(path, default_species))
            except (CrossSectionReadingError, ValueError, OSError):
                processes.extend(_read_bolsig_like(path, default_species))
        else:
            raise ValueError(f"Unsupported cross-section format: {fmt}")

    filled: list[CrossSectionProcess] = []
    for process in processes:
        if process.mass_amu is None and process.species in mass_by_species:
            process.mass_amu = mass_by_species[process.species]
        process.metadata.setdefault(
            "high_energy_extrapolation", config.high_energy_extrapolation
        )
        filled.append(process)

    missing_species = sorted(
        set(g.species for g in conditions.gas_mixture)
        - set(process.species for process in filled)
    )
    if missing_species:
        raise ValueError(f"No cross sections loaded for gas species: {missing_species}")
    return CrossSectionSet(filled)


def mixture_fraction(conditions: ConditionsConfig, species: str) -> float:
    for gas in conditions.gas_mixture:
        if gas.species == species:
            return gas.fraction
    return 0.0


def gas_mass_amu(conditions: ConditionsConfig, species: str) -> float:
    for gas in conditions.gas_mixture:
        if gas.species == species:
            return gas.mass_amu
    raise KeyError(species)
