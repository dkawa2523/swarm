#  Copyright (c) 2020-2021 ETH Zurich

"""
Module for importing cross section data, based on the 'lxcat_data_parser' package.

Data points are added to each cross section at zero energy and at max_energy (which
is user-defined). The classes 'CrossSection' and 'CrossSectionSet' of the
lxcat_data_parser are extended to include 'interp1d' linear interpolations of each
cross section.
"""

# Import Packages
import warnings
from typing import Union
import pandas as pd
import molmass
import numpy as np
import scipy.interpolate
import scipy.constants as csts
import matplotlib.pyplot as plt
from lxcat_data_parser import CrossSectionTypes as CST
from lxcat_data_parser import CrossSection, CrossSectionSet, CrossSectionReadingError

np.seterr(all='raise')


class InterpolatedCrossSection(CrossSection):
    """
    Extension of the CrossSection class to add a linear interpolation.

    Attributes:
        type (str, CrossSectionType): type of collision
        species (str): chemical formula of the species, example: N2
        mass_ratio (float): ratio of electron mass to molecular mass
        threshold (float): energy threshold (eV) of the cross section
        database (str): name of the database
        data (DataFrame): pandas DataFrame with columns "energy" and "cross section"
        info: optional additional information on the cross section given via kwargs
        interpolation(interp1d): linear interpolation of the cross section
    """

    def __init__(self, cross_section_type: Union[str, CST], species: str,
                 data: pd.DataFrame, mass_ratio: float,
                 threshold: float, **kwargs):
        """
        Instantiate an InterpolatedCrossSection.

        Args:
            cross_section_type (str, CrossSectionType): type of collision
            species (str): chemical formula of the species, example: N2
            data (DataFrame): pandas DataFrame with columns "energy" and "cross section"
            mass_ratio (float): ratio of electron mass to molecular mass
            threshold (float): energy threshold (eV) of the cross section
        """

        super().__init__(cross_section_type, species, data, mass_ratio,
                         threshold, **kwargs)
        self.interpolation = scipy.interpolate.interp1d(
            self.data['energy'].to_numpy(), self.data['cross section'].to_numpy(),
            kind='linear')


class InterpolatedCrossSectionSet(CrossSectionSet):
    """
    Extension of the CrossSectionSet class to use InterpolatedCrossSections.

    Attributes:
        species (str): chemical formula of the species, example: N2
        database (str): name of the database
        cross_sections (list): list of CrossSection instances
    """

    def __init__(self, max_cross_section_energy: float, input_file: str,
                 imposed_species: str, imposed_database: str = None):
        """
        Instantiate an InterpolatedCrossSectionSet.

        Args:
            max_cross_section_energy (float): maximum cross section energy (eV), which
                is appended to the end of the cross section data
            input_file (str): path to the cross section data
            imposed_species (str): chemical formula of the species, example: N2
            imposed_database (str): (optional) name of the database

        Raises:
            CrossSectionReadingError: if the provided file does not contain a valid set.
        """

        super().__init__(input_file, imposed_species, imposed_database)
        raw_max_energy = max(x.data['energy'].iat[-1] for x in self.cross_sections)
        if raw_max_energy + 1e-9 < max_cross_section_energy:
            warnings.warn(
                "max_cross_section_energy exceeds available data; "
                f"capping to file maximum {raw_max_energy} eV for {input_file}"
            )
            max_cross_section_energy = raw_max_energy
        # find the mass ratio value attached to the ELASTIC or EFFECTIVE cross section
        x = [x for x in self.cross_sections if x.type in (CST.EFFECTIVE, CST.ELASTIC)]
        if not x:
            raise CrossSectionReadingError(
                "The cross section set should contain exactly one EFFECTIVE "
                "CrossSection or one ELASTIC CrossSection. None was found.")
        elif len(x) > 1:
            raise CrossSectionReadingError(
                "The cross section set should contain exactly one EFFECTIVE "
                "CrossSection or one ELASTIC CrossSection. Several were found.")
        mass_ratio = x[0].mass_ratio

        # check the mass ratio using the molmass package
        try:
            molmass_mass_ratio = csts.electron_mass \
                / (molmass.Formula(self.species).mass / (1000 * csts.Avogadro))
            if not np.isclose(mass_ratio, molmass_mass_ratio):
                warnings.warn(
                    "Incorrect mass ratio."
                    f" The mass ratio {mass_ratio} read from the file "
                    f"'{input_file}' for the species {imposed_species} does not match "
                    f"the value {molmass_mass_ratio:.4} calculated with the molmass "
                    f"package. You may consider fixing the input file."
                )
        except molmass.FormulaError:
            warnings.warn(
                "Could not check the input mass ratio using the molmass package."
                f" {imposed_species} is not a valid chemical formula."
            )

        # add a factor of 2 inside the mass ratio for future energy loss calculations
        mass_ratio *= 2

        for x in self.cross_sections:
            if x.mass_ratio is None:
                x.mass_ratio = 0
            else:
                x.mass_ratio = mass_ratio
            # set threshold to 0 when not defined
            if x.threshold is None:
                x.threshold = 0
            # add data points (if needed) at zero energy and at max energy
            if x.data['energy'].iat[0] > 0:
                zero_energy_line = pd.DataFrame(
                    {'energy': [0],
                     'cross section': [x.data['cross section'].iat[0]]})
                x.data = pd.concat([zero_energy_line, x.data], ignore_index=True)
            if x.data['energy'].iat[-1] < max_cross_section_energy:
                final_energy_line = pd.DataFrame(
                    {'energy': [max_cross_section_energy],
                     'cross section': [x.data['cross section'].iat[-1]]})
                x.data = pd.concat([x.data, final_energy_line], ignore_index=True)
            elif x.data['energy'].iat[-1] > max_cross_section_energy:
                interp = scipy.interpolate.interp1d(
                    x.data['energy'], x.data['cross section'], kind='linear')
                x.data = x.data.loc[x.data['energy'] < max_cross_section_energy]
                final_energy_line = pd.DataFrame(
                    {'energy': [max_cross_section_energy],
                     'cross section': [interp(max_cross_section_energy)]})
                x.data = pd.concat([x.data, final_energy_line], ignore_index=True)

        # replace CrossSections by InterpolatedCrossSection
        self.cross_sections = [
            InterpolatedCrossSection(
                x.type, x.species, x.data, x.mass_ratio, x.threshold, info=x.info
            ) for x in self.cross_sections
        ]

        # replace if needed the effective cross section by an elastic one
        self.effective_to_elastic()

    def effective_to_elastic(self) -> None:
        """
        If effective cross section is given instead of elastic, calculate elastic cross
        section and replace the effective cross section.
        """

        # check if there is an effective cross section
        x_eff = [x for x in self.cross_sections if x.type == CST.EFFECTIVE]
        if x_eff:

            # combine the energy vectors into a total energy vector
            all_energies = np.concatenate([x.data['energy']
                                           for x in self.cross_sections])
            energy_vector = np.unique(all_energies)

            # calculate the elastic cross section to replace the effective cross section
            x_ela = x_eff[0]
            x_ela.type = CST.ELASTIC
            x_ela.data = pd.DataFrame({
                'energy': energy_vector,
                'cross section': x_ela.interpolation(energy_vector)
            })

            # subtract all the excitation and ionization cross sections
            for x in self.cross_sections:
                if x.type in [CST.EXCITATION, CST.IONIZATION]:
                    x_ela.data['cross section'] -= x.interpolation(x_ela.data['energy'])
            # update the linear interpolation
            x_ela.interpolation = scipy.interpolate.interp1d(
                x_ela.data['energy'].to_numpy(), x_ela.data['cross section'].to_numpy(),
                kind='linear')

    def plot(self, block=True) -> None:
        """
        Plot the cross section data.

        Args:
            block (bool): block execution or not when showing plot
        """

        fig = plt.figure()
        ax = fig.add_subplot(111)

        colors = ['black', 'does_not_occur', 'grey', 'blue', 'red']
        labels = []
        for x in self.cross_sections:
            if x.type not in labels:
                labels.append(x.type)
                ax.loglog(x.data['energy'], x.data['cross section'],
                          color=colors[x.type.value], label=x.type.name)
            else:
                ax.loglog(x.data['energy'], x.data['cross section'],
                          color=colors[x.type.value])

        # some cosmetics
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),
                   title='type of cross section')
        plt.title('Cross Section Set of ' + self.species)
        plt.xlabel('energy [eV]')
        plt.ylabel('cross section [m^2]')

        # to make space for the legend on the right
        plt.subplots_adjust(right=0.72)
        plt.show(block=block)
