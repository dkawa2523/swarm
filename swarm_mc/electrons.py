#  Copyright (c) 2020-2021 ETH Zurich

"""
Module for the Electrons class.

The electrons class stores the position, velocity and acceleration of electrons
in (3,n)-sized ndarrays. The class provides numerous class properties to compute
electrons-related data such as the mean electron energy using cached class properties.
"""

import math
from typing import Optional, Tuple, Union
# Import Packages
import numpy as np
import numpy.linalg as linalg
import matplotlib.pyplot as plt

# Import modules
import swarm_mc.utils as utils
from swarm_mc.energy_distribution import InstantaneousEnergyDistribution

np.seterr(all='raise')


class Electrons:
    """
    Stores information on a electron swarm in motion.

    Attributes:
        position (ndarray): xyz-coordinates of each electron, dim=(num_e,3)
        velocity (ndarray): xyz-velocities of each electron, dim=(num_e,3)
        acceleration (ndarray): xyz-accelerations of each electron, dim=(num_e,3)
    """

    def __init__(self, num_e_initial: int, init_pos: list,
                 initial_std: list, e_field: float,
                 initial_energy_distribution: str = "zero",
                 initial_energy: float = None,
                 initial_direction: Union[Tuple[float, float, float], str] = None,
                 initial_temperature: float = None,
                 rng: Optional[np.random.Generator] = None):
        """
        Instantiates an electron swarm.

        Args:
            num_e_initial (int): Initial number of electrons
            init_pos (list): initial xyz position (m) of electron swarm
            initial_std (list): initial xyz std (m) of electron swarm
            e_field (float): electric field strength (V/m) along z-direction
            initial_energy_distribution (str, optional): The initial energy
                distribution. See
                :py:attr:`swarm_mc.config.Config.initial_energy_distribution`.
                ``"maxwell-boltzmann"`` requires temperature to be set.
                Defaults to ``"zero"``
            initial_energy (float, optional): Initial electron energy in eV.
                Needed if initial_energy_distribution is not ``"zero"``.
            initial_direction (Union[Tuple[float, float, float], str], optional):
                See :py:attr:`swarm_mc.config.Config.initial_direction`.
            temperature (float, optional): Starting temperature in K. Required
                for Maxwell-Boltzmann distribution.
        """

        self.position = np.zeros([3, num_e_initial])
        self.rng = rng or np.random.default_rng()
        for i in range(3):
            self.position[i, :] = \
                init_pos[i] + self.rng.normal(loc=0, scale=initial_std[i],
                                              size=num_e_initial)

        if initial_energy_distribution == "maxwell-boltzmann":
            if initial_temperature is None:
                raise TypeError("argument temperature must be given for "
                                "Maxwell-Boltzmann distribution")

            if initial_temperature < 0:
                raise ValueError("temperature cannot be negative")

            velocities = np.zeros([3, num_e_initial])
            energies = utils.maxwell_boltzmann_random(num_e_initial,
                                                      initial_temperature,
                                                      rng=self.rng)
            for i, r in zip(range(num_e_initial), energies):
                theta = self.rng.random() * 2 * math.pi
                phi = self.rng.random() * math.pi
                x = r * math.sin(theta) * math.cos(phi)
                y = r * math.sin(theta) * math.sin(phi)
                z = r * math.cos(theta)
                velocities[0, i] = x
                velocities[1, i] = y
                velocities[2, i] = z
            self.velocity = velocities
        elif initial_energy_distribution == "zero" or initial_energy == 0:
            self.velocity = np.zeros([3, num_e_initial])
        elif initial_energy_distribution == "fixed":
            if initial_energy < 0:
                raise ValueError("initial_energy cannot be negative")
            if initial_direction == "random":
                velocities = np.zeros([3, num_e_initial])
                r = utils.velocity_from_energy(initial_energy)
                for i in range(num_e_initial):
                    theta = self.rng.random() * 2 * math.pi
                    phi = self.rng.random() * math.pi
                    x = r * math.sin(theta) * math.cos(phi)
                    y = r * math.sin(theta) * math.sin(phi)
                    z = r * math.cos(theta)
                    velocities[0, i] = x
                    velocities[1, i] = y
                    velocities[2, i] = z
                self.velocity = velocities
            else:
                if initial_direction == [0, 0, 0]:
                    raise ValueError("initial_direction cannot be all zero")

                normed = initial_direction / linalg.norm(initial_direction)
                initial_velocity = normed * utils.velocity_from_energy(initial_energy)
                self.velocity = \
                    np.full([num_e_initial, 3], initial_velocity).transpose()
        else:
            raise ValueError("Invalid value for initial_energy_distribution")

        self.acceleration = None
        self.accelerate(e_field)

        self.reset_cache()

    def accelerate(self, e_field: float) -> None:
        """
        Calculates the electrons acceleration in the electric field.

        Args:
            e_field (float): electric field strength (V/m) along z direction
        """

        self.acceleration = np.zeros([3, self.num_e])
        self.acceleration[2, :] = \
            utils.acceleration_from_electric_field(e_field) * np.ones(self.num_e)

    def free_flight(self, duration: float) -> None:
        """
        Update the attributes of Electrons after a free-flight in the electric field.

        Args:
            duration (float): duration (s) of the free flight
        """

        self.position += self.velocity * duration \
            + 0.5 * self.acceleration * duration ** 2
        self.velocity += self.acceleration * duration
        self.reset_cache()

    def apply_scatter(self, pos: np.ndarray, vel: np.ndarray, e_field: float) -> None:
        """
        Updates the attributes of Electrons after scattering by the gas.

        Args:
            pos (ndarray): new positions of the electrons (attachment and ionization)
            vel (ndarray): new velocities of the electrons (scattering)
            e_field (float): electric field strength (V/m) along z direction
        """

        self.position = pos
        self.velocity = vel
        self.accelerate(e_field)
        self.reset_cache()

    @property
    def num_e(self) -> int:
        return self.position.shape[1]

    @property
    def mean_position(self) -> np.ndarray:
        if self._mean_position is None:
            self._mean_position = np.mean(self.position, axis=1, keepdims=True).T
        return self._mean_position

    @property
    def var_position(self) -> np.ndarray:
        if self._var_position is None:
            self._var_position = np.var(self.position, axis=1, keepdims=True).T
        return self._var_position

    @property
    def mean_velocity(self) -> np.ndarray:
        if self._mean_velocity is None:
            self._mean_velocity = np.mean(self.velocity, axis=1, keepdims=True).T
        return self._mean_velocity

    @property
    def mean_velocity_moment(self) -> np.ndarray:
        if self._mean_velocity_moment is None:
            self._mean_velocity_moment = \
                np.mean(self.position * self.velocity, axis=1, keepdims=True).T
        return self._mean_velocity_moment

    @property
    def velocity_norm(self) -> np.ndarray:
        if self._velocity_norm is None:
            self._velocity_norm = linalg.norm(self.velocity, axis=0)
        return self._velocity_norm

    @property
    def max_velocity_norm(self) -> float:
        if self._max_velocity_norm is None:
            self._max_velocity_norm = np.max(self.velocity_norm)
        return self._max_velocity_norm

    @property
    def energy(self) -> np.ndarray:
        if self._energy is None:
            self._energy = utils.energy_from_velocity(self.velocity_norm)
        return self._energy

    @property
    def mean_energy(self) -> float:
        if self._mean_energy is None:
            self._mean_energy = np.float64(np.mean(self.energy))
        return self._mean_energy

    @property
    def std_energy(self) -> float:
        if self._std_energy is None:
            self._std_energy = np.float64(np.std(self.energy))
        return self._std_energy

    @property
    def max_energy(self) -> float:
        if self._max_energy is None:
            self._max_energy = np.max(self.energy)
        return self._max_energy

    @property
    def acceleration_norm(self) -> np.ndarray:
        if self._acceleration_norm is None:
            self._acceleration_norm = linalg.norm(self.acceleration, axis=0)
        return self._acceleration_norm

    @property
    def max_acceleration_norm(self) -> float:
        if self._max_acceleration_norm is None:
            self._max_acceleration_norm = np.max(self.acceleration_norm)
        return self._max_acceleration_norm

    @property
    def energy_distribution(self) -> InstantaneousEnergyDistribution:
        if self._energy_distribution.eedf is None:
            self._energy_distribution.calculate_distribution(self.energy)
        return self._energy_distribution

    def reset_cache(self) -> None:
        """
        Resets the cache to ensure it is updated.
        """

        self._mean_position = None
        self._var_position = None
        self._mean_velocity = None
        self._mean_velocity_moment = None
        self._energy = None
        self._mean_energy = None
        self._std_energy = None
        self._max_energy = None
        self._position_norm = None
        self._velocity_norm = None
        self._max_velocity_norm = None
        self._acceleration_norm = None
        self._max_acceleration_norm = None
        self._energy_distribution = InstantaneousEnergyDistribution()

    def plot_all(self, show: bool = True, block: bool = True) -> None:
        """
        Produces three Matplotlib Figures to visualize the position, velocity, and
        energy distribution of the electrons.

        Args:
            show (bool): calls plt.show() if True, else does nothing
            block (bool): if show is True, plt.show(block) is called
        """

        self.plot_position(show=False, block=False)
        self.plot_velocity(show=False, block=False)
        self.plot_energy(show=show, block=block)

    def plot_position(self, show: bool = True, block: bool = True) -> plt.Figure:
        """
        Produces a figure showing the electron positions. The figure contains four
        subplots showing the 3D scatter, and histograms along x, y, and z directions of
        the electron positions.

        Args:
            show (bool): calls plt.show() if True, else does nothing
            block (bool): if show is True, plt.show(block) is called
        Returns: Matplotlib figure object
        """

        fig = plt.figure()
        plt.suptitle("Electron positions")
        ax = fig.add_subplot(221, projection='3d')
        ax.scatter(self.position[0, :], self.position[1, :], self.position[2, :])
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.set_zlabel('z (m)')
        ax = fig.add_subplot(222)
        (_, bins) = np.histogram(self.position[0, :], bins='auto')
        ax.hist(self.position[0, :], bins=bins)
        ax.set_xlabel('x (m)')
        ax.set_ylabel('number of electrons')
        ax = fig.add_subplot(223)
        (_, bins) = np.histogram(self.position[1, :], bins='auto')
        ax.hist(self.position[1, :], bins=bins)
        ax.set_xlabel('y (m)')
        ax.set_ylabel('number of electrons')
        ax = fig.add_subplot(224)
        (_, bins) = np.histogram(self.position[2, :], bins='auto')
        ax.hist(self.position[2, :], bins=bins)
        ax.set_xlabel('z (m)')
        ax.set_ylabel('number of electrons')

        if show:
            plt.show(block=block)
        return fig

    def plot_velocity(self, show: bool = True, block: bool = True) -> plt.Figure:
        """
        Produces a figure showing the electron positions in velocity space. The figure
        contains four subplots showing the 3D scatter, and histograms along v_x, v_y,
        and v_z directions.

        Args:
            show (bool): calls plt.show() if True, else does nothing
            block (bool): if show is True, plt.show(block) is called
        Returns: Matplotlib figure object
        """

        fig = plt.figure()
        plt.suptitle("Electron velocities")
        ax = fig.add_subplot(221, projection='3d')
        ax.scatter(self.velocity[0, :], self.velocity[1, :], self.velocity[2, :])
        ax.set_xlabel('v$_x$ (m s$^{-1}$)')
        ax.set_ylabel('v$_y$ (m s$^{-1}$)')
        ax.set_zlabel('v$_z$ (m s$^{-1}$)')
        ax = fig.add_subplot(222)
        (_, bins) = np.histogram(self.velocity[0, :], bins='auto')
        ax.hist(self.velocity[0, :], bins=bins)
        ax.set_xlabel('v$_x$ (m s$^{-1}$)')
        ax.set_ylabel('number of electrons')
        ax = fig.add_subplot(223)
        (_, bins) = np.histogram(self.velocity[1, :], bins='auto')
        ax.hist(self.velocity[1, :], bins=bins)
        ax.set_xlabel('v$_y$ (m s$^{-1}$)')
        ax.set_ylabel('number of electrons')
        ax = fig.add_subplot(224)
        (_, bins) = np.histogram(self.velocity[2, :], bins='auto')
        ax.hist(self.velocity[2, :], bins=bins)
        ax.set_xlabel('v$_z$ (m s$^{-1}$)')
        ax.set_ylabel('number of electrons')

        if show:
            plt.show(block=block)
        return fig

    def plot_energy(self, show: bool = True, block: bool = True) -> plt.Figure:
        """
        Produces a figure showing the electron energy distribution function and the
        electron energy probability function.

        Args:
            show (bool): calls plt.show() if True, else does nothing
            block (bool): if show is True, plt.show(block) is called
        Returns: Matplotlib figure object
        """

        distri = self.energy_distribution

        fig, ax1 = plt.subplots()
        plt.title("Electron Energy Distribution Function (EEDF) \n"
                  "and Electron Energy Probability Function (EEPF)")
        color = 'tab:blue'
        ax1.hist(distri.energy_bin_centers, bins=distri.energy_bins,
                 weights=distri.eedf, color=color, histtype='step')
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.set_ylabel('EEDF (eV$^{-1}$)', color=color)
        ax1.set_xlabel('energy (eV)')
        ax2 = ax1.twinx()
        color = 'tab:red'
        ax2.hist(distri.energy_bin_centers, bins=distri.energy_bins,
                 weights=distri.eepf, color=color, histtype='step')
        ax2.tick_params(axis='y', labelcolor=color)
        ax2.set_ylabel('EEPF (eV$^{-3/2}$)', color=color)
        fig.tight_layout()

        if show:
            plt.show(block=block)
        return fig
