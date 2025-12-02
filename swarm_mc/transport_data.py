#  Copyright (c) 2020-2021 ETH Zurich

"""
Module for the FluxData and BulkData classes.

Both classes are based on the TransportData abstract class, and calculate the drift
velocity and diffusion coefficient of electrons based on the temporal evolution of their
position and velocity. The BulkData class calculates the bulk transport parameters,
which describe the transport of the center-of-mass of the swarm (time derivative of
spatial averages). The FluxData class calculates the flux transport parameters
(spatial averages of time derivatives).
"""

# Import Packages
from abc import ABC
import numpy as np
from scipy import stats

# Import modules
from swarm_mc.temporal_evolution import TimeSeries

np.seterr(all='raise')


class TransportData(ABC):
    """
    Stores the drift velocity and diffusion coefficient.

    Attributes:
        w (ndarray): drift velocity along x, y and z directions (m.s-1)
        w_err (ndarray): drift velocity error along x, y and z directions (m.s-1)
        DN (ndarray): diffusion coeff. along x, y and z directions (m-1.s-1)
        DN_err (ndarray): diffusion coeff. error along x, y and z directions (m-1.s-1)
    """

    def __init__(self, gas_number_density: float):
        self.gas_number_density = gas_number_density
        self.w = np.empty(3) + np.nan
        self.w_err = np.empty(3) + np.nan
        self.DN = np.empty(3) + np.nan
        self.DN_err = np.empty(3) + np.nan


class BulkData(TransportData):
    """
    Stores the bulk drift velocity and bulk diffusion coefficient.

    Attributes:
        w (ndarray): bulk drift velocity along x, y and z directions (m.s-1)
        w_err (ndarray): bulk drift velocity error along x, y and z directions (m.s-1)
        DN (ndarray): bulk diffusion coeff. along x, y and z directions (m-1.s-1)
        DN_err (ndarray): bulk diffusion coeff. error along x, y and z directions
            (m-1.s-1)
    """

    def calculate_data(self, time_series: TimeSeries) -> None:
        """
        Calculates the bulk drift velocity and bulk diffusion coefficient values
        (time derivative of spatial averages).

        Args:
            time_series (TimeSeries): temporal evolution of some quantities
        """

        x = time_series.time[time_series.ind_equ:]
        # for x,y,z directions
        for i in range(3):
            # drift coefficients
            y = time_series.mean_position[time_series.ind_equ:, i]
            slope, _, _, _, std_err = stats.linregress(x, y)
            self.w[i] = slope
            self.w_err[i] = 1.96 * std_err

            # diffusion coefficients
            y = time_series.var_position[time_series.ind_equ:, i]
            slope, _, _, _, std_err = stats.linregress(x, y)
            self.DN[i] = 0.5 * self.gas_number_density * slope
            self.DN_err[i] = 0.5 * self.gas_number_density * 1.96 * std_err


class FluxData(TransportData):
    """
    Flux drift velocity and flux diffusion coefficient.

    Attributes:
        w (ndarray): flux drift velocity along x, y and z directions (m.s-1)
        w_err (ndarray): flux drift velocity error along x, y and z directions (m.s-1)
        DN (ndarray): flux diffusion coeff. along x, y and z directions (m-1.s-1)
        DN_err (ndarray): flux diffusion coeff. error along x, y and z directions
            (m-1.s-1)
    """

    def __init__(self, gas_number_density: float):
        super().__init__(gas_number_density)

    def calculate_data(self, time_series: TimeSeries) -> None:
        """
        Calculates the flux drift velocity and the flux diffusion coefficient
        (spatial averages of time derivatives).

        Args:
            time_series (TimeSeries): temporal evolution of some quantities
        """

        steps_since_sst = time_series.time.size - 1 - time_series.ind_equ

        # calculate the drift velocity
        # calculate the drift velocity
        self.w = np.mean(time_series.mean_velocity[time_series.ind_equ:, :], axis=0)

        # calculate the diffusion coefficient
        y = time_series.mean_velocity_moment \
            - time_series.mean_position * time_series.mean_velocity
        self.DN = self.gas_number_density * np.mean(y[time_series.ind_equ:, :], axis=0)

        # At least 10 steps are required to estimate the errors
        if steps_since_sst > 9:
            # calculate the drift velocity error
            self.w_err = \
                np.std(time_series.mean_velocity[time_series.ind_equ:, :], axis=0) \
                / np.sqrt(steps_since_sst)

            # calculate the diffusion coefficient error
            y = time_series.mean_velocity_moment \
                - time_series.mean_position * time_series.mean_velocity
            self.DN_err = self.gas_number_density \
                * np.std(y[time_series.ind_equ:, :], axis=0) / np.sqrt(steps_since_sst)
