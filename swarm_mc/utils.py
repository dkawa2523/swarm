#  Copyright (c) 2020-2021 ETH Zurich

"""
Module with some utility methods.
"""

from typing import Optional, Union
import numpy as np
import scipy.constants as sc

num = Union[int, float, np.ndarray]


def velocity_from_energy(energy: num) -> num:
    """
    Calculate the velocity norm of electrons, based on their kinetic energy.

    Args:
        energy: Energy of the electron (eV)

    Returns: Velocity norm of the electrons (m.s-1)
    """

    return np.sqrt(2 * energy * sc.elementary_charge / sc.electron_mass)


def energy_from_velocity(velocity: num) -> num:
    """
    Calculate the kinetic energy of electrons, based on the norm of their velocity.

    Args:
        velocity: velocity norm of the electrons (m.s-1)

    Returns: Energy of the electrons (eV)
    """

    return 0.5 * sc.electron_mass * velocity ** 2 / sc.elementary_charge


def acceleration_from_electric_field(electric_field: num) -> num:
    """
    Calculates the acceleration of electrons, based on the local electric field.

    Args:
        electric_field: local electric field strength (V.m-1)

    Returns: acceleration (m.s-2)
    """

    return electric_field * (sc.elementary_charge / sc.electron_mass)


def maxwell_boltzmann_eedf(points: num, temperature: Union[int, float]) -> num:
    """
    Calculates the Maxwell-Boltzmann EEDF of the given points or point

    Args:
        points (Union[int, float, np.ndarray]): Scalar or vector of energies in eV
        temperature (Union[int, float]): Temperature at which to calculate the
            distribution
    """
    kB = sc.k
    me = sc.m_e
    q0 = sc.e
    T = temperature
    v = np.sqrt(2 * points * q0 / me)
    mb = np.sqrt(2/np.pi) * (me/(kB*T))**(3/2) * v**2 * np.exp(-me * v**2 / (kB*T))
    return mb * np.sqrt(points)


def maxwell_boltzmann_random(num: int, temperature: Union[float, int],
                             rng: Optional[np.random.Generator] = None) -> list:
    """
    Returns a list of `num` Maxwell Boltzmann distributed values (in eV) at
    temperature `temperature`

    Args:
        num (int): Number of electrons
        temperature (float): Temperature
    """

    rng = rng or np.random.default_rng()
    d = False
    i = 0
    last = -1
    di = 0.0001

    # determine the range for values
    while not d or not np.isclose(last, 0):
        now = maxwell_boltzmann_eedf(i, temperature)
        if now < last:
            d = True

        last = now
        i += di

    # divide the range into subintervals
    bins = np.arange(0, i, di)
    vals = []

    # for each subinterval associate a probability (the mean of the eedf value
    # at the right and left side of the subinterval)
    for b in bins:
        left = maxwell_boltzmann_eedf(b, temperature)
        right = maxwell_boltzmann_eedf(b+di, temperature)
        m = left + right / 2
        vals.append(m)

    # sum all probabilities together (imagine putting them all on a line)
    total = sum(vals)
    e = []

    # generate num values by choosing a random point the line and select a
    # random energy in the corresponding energy interval. Both random numbers
    # are chosen uniformly.
    for i in range(0, num):
        u = rng.uniform(0, total)

        s = 0
        for b, v in zip(bins, vals):
            s += v

            if u <= s:
                e.append(rng.uniform(b, b + di))
                break

    return e
