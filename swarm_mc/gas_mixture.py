#  Copyright (c) 2020-2021 ETH Zurich

"""
Module for the GasMixture class, which aggregates data of different gases.

The GasMixture uses the cross_section module to read cross section data. It stores
linear interpolations of each cross section, as well as the mass ratio,
energy threshold and type of cross section, which are needed for the simulation.
The proportions of the gas mixture should sum up to 1. A common energy_vector for all
cross sections is calculated, as well as the total_cross_section.
"""

# Import Packages
from dataclasses import dataclass
import re
import numpy as np
from lxcat_data_parser import CrossSectionTypes as CST

# Import modules
from swarm_mc.cross_section import InterpolatedCrossSectionSet

np.seterr(all='raise')


class GasMixture:
    """
    Class representing a gas mixture.

    Attributes:
        cross_sections (ndarray): array of cross section interpolations (from all gases)
        types (ndarray): type of each cross section of each species
        thresholds (ndarray): threshold of each cross section of each species
        mass_ratios (ndarray): mass ratios (repeated for each cross section)
        is_attachment (ndarray): boolean mask for ATTACHMENT cross sections
        is_ionization (ndarray): boolean mask for IONIZATION cross sections
        energy_vector (ndarray): common energy vector, containing all distinct energies
            of all cross sections
        total_cross_section (ndarray): sum of all cross sections of all species,
            linearly interpolated at energy_vector
    """

    def __init__(self, gas_formulas: list, path_to_cross_section_files: list,
                 proportions: list, max_cross_section_energy: float):
        """
        Instantiates a GasMixture.

        Args:
            gas_formulas (list): list of the chemical formulas of the species
            path_to_cross_section_files (list): path to the cross section data of each
                species
            proportions (list): proportion of each species in the mixture
            max_cross_section_energy (float): maximum cross section energy (eV)

        Raises:
            ValueError: If the proportions do not sum up to 1
        """

        if not np.isclose(sum(proportions), 1, atol=1e-3):
            raise ValueError(
                "The sum of the proportions of the gases in the GasMixture must be 1 "
                f"(±1e-3), but it is {sum(proportions):.4f}"
            )

        gases = [InterpolatedCrossSectionSet(max_cross_section_energy, *args)
                 for args in zip(path_to_cross_section_files, gas_formulas)]

        cross_sections = [x for gas in gases for x in gas.cross_sections]
        self.cross_section_meta = _build_cross_section_meta(cross_sections)
        self.thresholds = np.ascontiguousarray(
            np.asarray([x.threshold for x in cross_sections], dtype=float)
        )
        self.mass_ratios = np.ascontiguousarray(
            np.asarray([x.mass_ratio for x in cross_sections], dtype=float)
        )
        self.types = [x.type for x in cross_sections]
        types_array = np.asarray(self.types)
        self.is_attachment = np.ascontiguousarray(types_array == CST.ATTACHMENT)
        self.is_ionization = np.ascontiguousarray(types_array == CST.IONIZATION)

        # combine the energy vectors into a total energy vector
        all_energies = np.concatenate([x.data['energy'] for x in cross_sections])
        self.energy_vector = np.unique(all_energies).astype(float)

        # pre-sample cross sections onto the common grid to avoid runtime interp1d
        scaled_tables = []
        for proportion, gas in zip(proportions, gases):
            for section in gas.cross_sections:
                values = np.interp(self.energy_vector,
                                   section.data['energy'],
                                   section.data['cross section'])
                scaled_tables.append(values * proportion)
        self.cross_section_table = np.ascontiguousarray(
            np.asarray(scaled_tables, dtype=float)
        )
        # keep callable cross sections for downstream rate calculations
        self.cross_sections = np.asarray([
            (lambda e, row=row: np.interp(e, self.energy_vector, row))
            for row in self.cross_section_table
        ], dtype=object)

        # calculate the total cross section
        self.total_cross_section = np.sum(self.cross_section_table, axis=0)

    @property
    def number_of_cross_sections(self) -> int:
        """
        int: Total number of cross sections, all gases considered.
        """

        return self.cross_section_table.shape[0]


@dataclass(frozen=True)
class CrossSectionMeta:
    index: int
    ctype: CST
    species: str
    process: str
    label: str


def _extract_process(info) -> str:
    if not isinstance(info, dict):
        return ""
    if "PROCESS" in info:
        return str(info.get("PROCESS") or "")
    nested = info.get("info")
    if isinstance(nested, dict):
        return str(nested.get("PROCESS") or "")
    return ""


def _clean_label(text: str) -> str:
    cleaned = re.sub(r"\\s+", " ", text.replace("\n", " ")).strip()
    cleaned = cleaned.rstrip(",;")
    return cleaned


def _build_cross_section_meta(cross_sections) -> list:
    metas = []
    seen = {}
    for i, x in enumerate(cross_sections):
        process = _extract_process(getattr(x, "info", None))
        base = process if process else f"{x.type.name} {x.species}"
        label = _clean_label(base)
        if not label:
            label = f"{x.type.name} {x.species} {i+1}"
        count = seen.get(label, 0) + 1
        seen[label] = count
        if count > 1:
            label = f"{label} ({count})"
        metas.append(CrossSectionMeta(
            index=i,
            ctype=x.type,
            species=x.species,
            process=process,
            label=label,
        ))
    return metas
