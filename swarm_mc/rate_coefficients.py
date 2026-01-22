#  Copyright (c) 2020-2021 ETH Zurich

"""
Module for the CountedRates and ConvolutedRates classes.

Both classes are based on the RateCoefficients abstract class. They implement two
different methods for calculating ionization/attachment rate coefficients. The
CountedRates class calculates the rate coefficients by averaging over time the number
of ionization and attachment collisions. The ConvolutedRates class calculates the
rate coefficients by calculating the convolution product of the cross section with the
time-averaged electron energy distribution of the electrons.
"""

# Import Packages
from dataclasses import dataclass
from abc import ABC
import numpy as np
import scipy.constants as csts
from lxcat_data_parser import CrossSectionTypes as CST

# Import modules
from swarm_mc.energy_distribution import TimeAveragedEnergyDistribution
from swarm_mc.gas_mixture import GasMixture
from swarm_mc.temporal_evolution import TimeSeries

np.seterr(all='raise')


class RateCoefficients(ABC):
    """
    Ionization, attachment and effective ionization rate coefficients

    Attributes:
        ionization: ionization reaction rate coefficient (m3.s-1)
        attachment: attachment reaction rate coefficient (m3.s-1)
        effective: effective ionization reaction rate coefficient (m3.s-1)
    """

    def __init__(self):

        self.ionization = None
        self.attachment = None
        self.effective = None


@dataclass(frozen=True)
class ReactionRate:
    label: str
    ctype: CST
    process: str
    rate: float


class ConvolutedRates(RateCoefficients):
    """
    Ionization, attachment and effective ionization rate coefficients calculated by
    convolution product of eedf and cross section.
    """

    def __init__(self):
        super().__init__()
        self.per_reaction = []

    def calculate_data(self, mix: GasMixture,
                       distri: TimeAveragedEnergyDistribution) -> None:
        """
        Calculates the effective, ionization and attachment rate coefficients.

        Args:
            mix (GasMixture): GasMixture object containing the cross section data
            distri (TimeAveragedEnergyDistribution): electron energy distribution (eV-1)
        """

        self.per_reaction = []
        ionization_rates = []
        attachment_rates = []
        for interp, meta in zip(mix.cross_sections, mix.cross_section_meta):
            if meta.ctype not in (CST.IONIZATION, CST.ATTACHMENT):
                continue
            rate = self.calculate_convolution(interp, distri)
            self.per_reaction.append(ReactionRate(
                label=meta.label,
                ctype=meta.ctype,
                process=meta.process,
                rate=rate,
            ))
            if meta.ctype == CST.IONIZATION:
                ionization_rates.append(rate)
            else:
                attachment_rates.append(rate)

        self.ionization = float(np.sum(ionization_rates)) if ionization_rates else 0.0
        self.attachment = float(np.sum(attachment_rates)) if attachment_rates else 0.0
        self.effective = self.ionization - self.attachment

    @staticmethod
    def calculate_convolution(interp,
                              distri: TimeAveragedEnergyDistribution) -> float:
        """
        Calculates the reaction rate of a collision type by calculating the convolution
        product of the eedf and the corresponding cross sections.

        Args:
            interp: linear interpolation of the cross section (m2)
            distri (TimeAveragedEnergyDistribution): energy distribution of the
                electrons (eV-1)

        Returns: reaction rate coefficient (m3.s-1)
        """

        cross_section = interp(distri.energy_bin_centers)

        A = distri.eedf * distri.energy_bin_centers ** 0.5 * np.diff(distri.energy_bins)
        rrate = np.sqrt(
            2 * csts.elementary_charge / csts.electron_mass) * np.dot(A, cross_section)
        return rrate


class CountedRates(RateCoefficients):
    """
    Ionization, attachment and effective ionization rate coefficients, calculated by
    counting ionization and attachment collisions.

    Attributes:
        gas_number_density (float): gas number density (m-3)
        conserve (bool): conservation of the electron number
    """

    def __init__(self, gas_number_density: float, conserve: bool):
        """
        Instantiation of CountedRates.

        Args:
            gas_number_density (float): gas number density (m-3) for normalization
            conserve (bool): bool indicating conservation of the electron number
                (affects the statistics when counting the collisions)
        """

        super().__init__()

        self.gas_number_density = gas_number_density
        self.conserve = conserve
        self.ionization_err = None
        self.attachment_err = None
        self.effective_err = None

    def calculate_data(self, time_series: TimeSeries) -> None:
        """
        Calculates the ionization, attachment and effective ionization rate coefficients
        values by counting collision events.

        Args:
            time_series (TimeSeries): temporal evolution of some quantities
        """

        steps_since_sst = time_series.time.size - 1 - time_series.ind_equ
        x = time_series.time[time_series.ind_equ:] \
            - time_series.time[time_series.ind_equ]

        if not self.conserve:
            # Ne(t-t_SST) = Ne(t_SST) * exp(nu_eff * (t-t_SST))
            # log(Ne(t-t_SST)) - log(Ne(t_SST)) = nu_eff * (t-t_SST)
            y = np.log(time_series.num_electrons[time_series.ind_equ:]) \
                - np.log(time_series.num_electrons[time_series.ind_equ])
            d = y[1:] / x[1:]
            nu_eff = np.mean(d)
            if nu_eff != 0:
                x = (np.exp(nu_eff * x) - 1) / nu_eff \
                    * time_series.num_electrons[time_series.ind_equ]
        else:
            y = time_series.num_electrons[time_series.ind_equ:] \
                - time_series.num_electrons[time_series.ind_equ]
            d = y[1:] / x[1:]

        # effective ionization
        self.effective = np.mean(d) / self.gas_number_density
        self.effective_err = np.std(d) / np.sqrt(steps_since_sst) \
            / self.gas_number_density

        # ionization (cations, positively charged)
        y = time_series.num_cations[time_series.ind_equ:] \
            - time_series.num_cations[time_series.ind_equ]
        d = y[1:] / x[1:]
        self.ionization = np.mean(d) / self.gas_number_density
        self.ionization_err = \
            np.std(d) / np.sqrt(steps_since_sst) / self.gas_number_density

        # attachment (anions, negatively charged)
        y = time_series.num_anions[time_series.ind_equ:] \
            - time_series.num_anions[time_series.ind_equ]
        d = y[1:] / x[1:]
        self.attachment = np.mean(d) / self.gas_number_density
        self.attachment_err = \
            np.std(d) / np.sqrt(steps_since_sst) / self.gas_number_density
