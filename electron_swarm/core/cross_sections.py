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

from .config import ConditionsConfig, CrossSectionsConfig, GasComponent


class _NotLongCsv(ValueError):
    """Internal signal that a CSV uses the wide representation."""


class CrossSectionValidationError(ValueError):
    """A parsed cross-section process contains invalid physical data."""


class ProcessType(str, Enum):
    MOMENTUM = "momentum"
    ELASTIC = "elastic"
    EFFECTIVE = "effective"
    EXCITATION = "excitation"
    IONIZATION = "ionization"
    ATTACHMENT = "attachment"
    SUPERELASTIC = "superelastic"
    UNKNOWN = "unknown"


_INCIDENT_THRESHOLD_TYPES = frozenset(
    {
        ProcessType.ATTACHMENT,
        ProcessType.EXCITATION,
        ProcessType.IONIZATION,
    }
)


class ScatteringRole(str, Enum):
    """Physical role of an ordinary integral scattering cross section.

    ProcessType describes the collision family used by result tables.
    ScatteringRole records which integral the input actually supplies.
    Keeping the two separate prevents an effective or momentum-transfer
    cross section from silently becoming a particle collision frequency.
    """

    ELASTIC_TOTAL = "elastic_total"
    ELASTIC_MOMENTUM_TRANSFER = "elastic_momentum_transfer"
    EFFECTIVE_MOMENTUM_TRANSFER = "effective_momentum_transfer"


_TYPE_ALIASES = {
    "momentum": ProcessType.MOMENTUM,
    "momentum_transfer": ProcessType.MOMENTUM,
    "mt": ProcessType.MOMENTUM,
    "elastic": ProcessType.ELASTIC,
    "elastic_total": ProcessType.ELASTIC,
    "total_elastic": ProcessType.ELASTIC,
    "elastic_momentum_transfer": ProcessType.MOMENTUM,
    "effective": ProcessType.EFFECTIVE,
    "effective_momentum_transfer": ProcessType.EFFECTIVE,
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

_SCATTERING_ROLE_ALIASES = {
    "elastic_total": ScatteringRole.ELASTIC_TOTAL,
    "total_elastic": ScatteringRole.ELASTIC_TOTAL,
    "elastic_momentum_transfer": ScatteringRole.ELASTIC_MOMENTUM_TRANSFER,
    "momentum": ScatteringRole.ELASTIC_MOMENTUM_TRANSFER,
    "momentum_transfer": ScatteringRole.ELASTIC_MOMENTUM_TRANSFER,
    "mt": ScatteringRole.ELASTIC_MOMENTUM_TRANSFER,
    "effective": ScatteringRole.EFFECTIVE_MOMENTUM_TRANSFER,
    "effective_momentum_transfer": ScatteringRole.EFFECTIVE_MOMENTUM_TRANSFER,
}


def normalize_process_type(value: str | None) -> ProcessType:
    if value is None:
        return ProcessType.UNKNOWN
    key = str(value).strip().lower().replace(" ", "_").replace("-", "_")
    return _TYPE_ALIASES.get(key, ProcessType.UNKNOWN)


def normalize_scattering_role(
    value: str | None,
    *,
    source_format: str = "canonical",
) -> ScatteringRole | None:
    """Map an input label to its explicit scattering integral role.

    LXCat/BOLSIG ELASTIC is a momentum-transfer cross section. A canonical
    CSV must say elastic_total or total_elastic to identify a total elastic
    cross section. Bare elastic remains accepted in programmatic construction
    through CrossSectionProcess where its role defaults to total, but file
    loaders never infer total from that word.
    """

    if value is None:
        return None
    key = str(value).strip().lower().replace(" ", "_").replace("-", "_")
    if key == "elastic":
        if source_format.lower() in {"lxcat", "bolsig", "txt", "dat", "csv"}:
            return ScatteringRole.ELASTIC_MOMENTUM_TRANSFER
        return ScatteringRole.ELASTIC_TOTAL
    return _SCATTERING_ROLE_ALIASES.get(key)


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


def _scattering_role_from_lxcat(value: object) -> ScatteringRole | None:
    if value == CST.EFFECTIVE:
        return ScatteringRole.EFFECTIVE_MOMENTUM_TRANSFER
    if value == CST.ELASTIC:
        return ScatteringRole.ELASTIC_MOMENTUM_TRANSFER
    return None


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
    scattering_role: ScatteringRole | None = None
    threshold_eV: float | None = None
    mass_amu: float | None = None
    metadata: dict[str, str | float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        energy = np.asarray(self.energy_eV, dtype=float)
        sigma = np.asarray(self.cross_section_m2, dtype=float)
        if energy.ndim != 1 or sigma.ndim != 1 or len(energy) != len(sigma):
            raise CrossSectionValidationError(
                f"Invalid cross-section arrays for {self.species}:{self.process}"
            )
        if len(energy) < 2:
            raise CrossSectionValidationError(
                f"At least two energy points are required for {self.species}:{self.process}"
            )
        invalid_energy = np.flatnonzero(~np.isfinite(energy) | (energy < 0.0))
        if invalid_energy.size:
            raise CrossSectionValidationError(
                f"Invalid energy at row {int(invalid_energy[0])} for "
                f"{self.species}:{self.process}"
            )
        invalid_sigma = np.flatnonzero(~np.isfinite(sigma) | (sigma < 0.0))
        if invalid_sigma.size:
            raise CrossSectionValidationError(
                f"Invalid cross section at row {int(invalid_sigma[0])} for "
                f"{self.species}:{self.process}"
            )
        if self.threshold_eV is not None and (
            not np.isfinite(float(self.threshold_eV))
            or float(self.threshold_eV) < 0.0
        ):
            raise CrossSectionValidationError(
                f"Invalid threshold for {self.species}:{self.process}"
            )
        if self.mass_amu is not None and (
            not np.isfinite(float(self.mass_amu)) or float(self.mass_amu) <= 0.0
        ):
            raise CrossSectionValidationError(
                f"Invalid target mass for {self.species}:{self.process}"
            )
        if self.scattering_role is None:
            self.scattering_role = {
                ProcessType.ELASTIC: ScatteringRole.ELASTIC_TOTAL,
                ProcessType.MOMENTUM: ScatteringRole.ELASTIC_MOMENTUM_TRANSFER,
                ProcessType.EFFECTIVE: ScatteringRole.EFFECTIVE_MOMENTUM_TRANSFER,
            }.get(self.process_type)
        elif not isinstance(self.scattering_role, ScatteringRole):
            self.scattering_role = ScatteringRole(str(self.scattering_role))
        if self.process_type not in {
            ProcessType.ELASTIC,
            ProcessType.MOMENTUM,
            ProcessType.EFFECTIVE,
        } and self.scattering_role is not None:
            raise CrossSectionValidationError(
                f"Non-scattering process {self.species}:{self.process} "
                "cannot declare a scattering role"
            )

        order = np.argsort(energy, kind="stable")
        energy = energy[order]
        sigma = sigma[order]
        uniq = np.unique(energy)
        if len(uniq) != len(energy):
            sigma2 = np.zeros_like(uniq)
            for j, value in enumerate(uniq):
                sigma2[j] = sigma[energy == value][-1]
            energy, sigma = uniq, sigma2
        self.energy_eV = energy
        self.cross_section_m2 = sigma

    @property
    def incident_threshold_eV(self) -> float | None:
        """Return the lower incident-energy support, when physically applicable."""

        if (
            self.process_type in _INCIDENT_THRESHOLD_TYPES
            and self.threshold_eV is not None
        ):
            return float(self.threshold_eV)
        return None

    def sigma(
        self, energy_eV: np.ndarray, *, left: float = 0.0, right: float | None = None
    ) -> np.ndarray:
        """Linearly interpolate cross section on an arbitrary energy grid."""

        energy = np.asarray(energy_eV, dtype=float)
        if right is None:
            policy = str(
                self.metadata.get("high_energy_extrapolation", "hold")
            ).lower()
            if policy == "zero":
                right = 0.0
            elif policy == "hold":
                right = float(self.cross_section_m2[-1])
            elif policy == "error":
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
        values = np.interp(
            energy, self.energy_eV, self.cross_section_m2, left=left, right=right
        )
        threshold = self.incident_threshold_eV
        if threshold is not None:
            values = np.where(energy < threshold, 0.0, values)
        return values


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


@dataclass(slots=True)
class ActiveMixtureInputs(CrossSectionSet):
    """Solver-ready cross sections for positive-fraction gas components only."""

    components: tuple[GasComponent, ...]

    @property
    def active_species(self) -> tuple[str, ...]:
        return tuple(component.species for component in self.components)


_REQUIRED_THRESHOLD_TYPES = frozenset(
    {
        ProcessType.EXCITATION,
        ProcessType.IONIZATION,
        ProcessType.SUPERELASTIC,
    }
)


def prepare_active_mixture_inputs(
    cross_sections: CrossSectionSet,
    conditions: ConditionsConfig,
) -> ActiveMixtureInputs:
    """Filter and validate the one cross-section inventory used by every solver."""

    components = tuple(
        component
        for component in conditions.gas_mixture
        if component.fraction > 0.0
    )
    active_species = {component.species for component in components}
    if isinstance(cross_sections, ActiveMixtureInputs):
        if cross_sections.components != components:
            raise CrossSectionValidationError(
                "Active cross-section inventory does not match gas mixture"
            )
        _validate_active_processes(cross_sections.processes, active_species)
        return ActiveMixtureInputs(
            processes=list(cross_sections.processes),
            components=components,
        )

    processes = [
        process
        for process in cross_sections.processes
        if process.species in active_species
    ]

    _validate_active_processes(processes, active_species)
    return ActiveMixtureInputs(processes=processes, components=components)


def _validate_active_processes(
    processes: list[CrossSectionProcess], active_species: set[str]
) -> None:
    """Validate both newly filtered and already materialized active inputs."""

    unexpected_species = sorted(
        {process.species for process in processes} - active_species
    )
    if unexpected_species:
        raise CrossSectionValidationError(
            "Active cross-section inventory contains inactive gas species: "
            f"{unexpected_species}"
        )

    missing_species = sorted(
        active_species - {process.species for process in processes}
    )
    if missing_species:
        raise CrossSectionValidationError(
            f"No cross sections loaded for active gas species: {missing_species}"
        )

    for process in processes:
        identity = f"{process.species}:{process.process}"
        if process.process_type == ProcessType.UNKNOWN:
            raise CrossSectionValidationError(
                f"Unknown cross-section process type for active process {identity}"
            )
        if (
            process.process_type in _REQUIRED_THRESHOLD_TYPES
            and process.threshold_eV is None
        ):
            raise CrossSectionValidationError(
                "threshold_eV is required for active "
                f"{process.process_type.value} process {identity}"
            )


def _read_csv_long(path: Path, default_species: str | None) -> list[CrossSectionProcess]:
    df = pd.read_csv(path)
    columns = {c.lower(): c for c in df.columns}
    required = {"energy_ev", "cross_section_m2"}
    if not required.issubset(columns):
        raise _NotLongCsv("not long CSV")
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
        raw_type = (
            str(key_map.get(type_col, "unknown"))
            if type_col is not None
            else str(key_map.get(process_col, "unknown"))
        )
        ptype = normalize_process_type(raw_type)
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
                scattering_role=normalize_scattering_role(
                    raw_type,
                    source_format="csv",
                ),
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
                scattering_role=normalize_scattering_role(
                    col.removesuffix("_m2"),
                    source_format="csv",
                ),
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
    active = [
        component
        for component in conditions.gas_mixture
        if component.fraction > 0.0
    ]
    if len(active) == 1:
        return active[0].species
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
                scattering_role=_scattering_role_from_lxcat(section.type),
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
                    scattering_role=normalize_scattering_role(
                        line,
                        source_format="bolsig",
                    ),
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
    """Load all configured cross-section files without selecting a mixture."""

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
            except _NotLongCsv:
                processes.extend(_read_csv_wide(path, default_species, default_mass))
        elif fmt in {"bolsig", "lxcat", "txt", "dat"}:
            try:
                processes.extend(_read_lxcat_text(path, default_species))
            except CrossSectionValidationError:
                raise
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
        {
            component.species
            for component in conditions.gas_mixture
            if component.fraction > 0.0
        }
        - {process.species for process in filled}
    )
    if missing_species:
        raise ValueError(
            f"No cross sections loaded for active gas species: {missing_species}"
        )
    return CrossSectionSet(filled)


def load_active_mixture_inputs(
    config: CrossSectionsConfig,
    conditions: ConditionsConfig,
) -> ActiveMixtureInputs:
    """Load, select, and validate the inventory for one solver mixture."""

    return prepare_active_mixture_inputs(
        load_cross_sections(config, conditions),
        conditions,
    )


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
