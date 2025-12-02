#  Copyright (c) 2020-2021 ETH Zurich

"""
Module for the Output class.
"""

# Import Packages
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Import modules
from swarm_mc.config import Config
from swarm_mc.electrons import Electrons
from swarm_mc.energy_distribution import TimeAveragedEnergyDistribution
from swarm_mc.transport_data import BulkData, FluxData
from swarm_mc.rate_coefficients import ConvolutedRates, CountedRates
from swarm_mc.temporal_evolution import TimeSeries

np.seterr(all='raise')


class Output:
    """
    The Output class instantiates all output-related classes: TimeSeries,
    TimeAveragedEnergyDistribution, FluxData, BulkData, ConvolutedRates, CountedRates.
    It also provides methods to save or to plot the output data. Finally, it provides
    the check_sst method which test if the swarm is at equilibrium based on the
    evolution of the mean electron energy.

    Attributes:
        config (Config): configuration of the simulation
        version (str): version of swarm_mc as a string
        time_series (TimeSeries): temporal evolution data
        energy_distribution (TimeAveragedEnergyDistribution): energy distribution data
        flux (FluxData): flux transport data
        bulk (BulkData): flux transport data
        rates_conv (ConvolutedRates): convoluted rate coefficients
        rates_count (CountedRates): counted rate coefficients
    """

    def __init__(self, cfg: Config, ver: str, electrons: Electrons):
        """
        Instantiates the Output class.

        Args:
            cfg (Config): configuration of the simulation
            ver (str): version of swarm_mc as a string
            electrons (Electrons): electron data
        """

        self.config = cfg
        self.version = ver

        # temporal evolution of some quantities
        self.time_series = TimeSeries(electrons, log_every_k=cfg.log_every_k)

        # energy-related data:
        self.energy_distribution = TimeAveragedEnergyDistribution()

        # flux data
        self.flux = FluxData(self.config.gas_number_density)

        # bulk data
        self.bulk = BulkData(self.config.gas_number_density)

        # rate coefficients calculated by convolution
        self.rates_conv = ConvolutedRates()

        # rate coefficients calculated by counting events
        self.rates_count = CountedRates(cfg.gas_number_density, cfg.conserve)

    # ------------------------------------------------------------------
    # Helpers
    def _wants_csv(self) -> bool:
        return self.config.output_format in ("csv", "both")

    def _wants_json(self) -> bool:
        return self.config.output_format in ("json", "both")

    def _path(self, name: str, extension: str) -> Path:
        return Path(self.config.output_directory) / f"{name}.{extension}"

    @staticmethod
    def _safe_float(value):
        return float(value) if value is not None else np.nan

    def check_sst(self) -> bool:
        """
        Checks if the swarm energy is at equilibrium (steady-state). This is done by
        evaluating a linear trend and residual variance over a sliding window with a
        conservative fallback to the legacy percentile check.

        Returns: True if equilibrium was reached, False otherwise
        """

        n_points = self.time_series.mean_energy.size
        window = min(self.config.sst_max_window, n_points)
        if window < self.config.sst_min_window:
            return False

        energies = self.time_series.mean_energy[-window:]
        times = self.time_series.time[-window:]
        mean_energy = np.mean(energies)
        if mean_energy == 0:
            return False

        slope, intercept = np.polyfit(times, energies, 1)
        fitted = slope * times + intercept
        residual = energies - fitted
        trend_ok = abs(slope) / mean_energy < 0.02 \
            and np.std(residual) / mean_energy < 0.08

        # Legacy fallback: last 10% not greater than previous 10%
        legacy_ok = False
        n = max(round(n_points / 10), 1)
        if n >= 10:
            legacy_ok = np.mean(self.time_series.mean_energy[-2 * n:-n]) >= \
                        np.mean(self.time_series.mean_energy[-n:])

        if trend_ok or legacy_ok:
            self.time_series.ind_equ = self.time_series.time.size
            return True
        return False

    def save_temporal_evolution(self, name: str = None) -> None:
        """
        Saves the temporal evolution of the swarm data to a csv and/or json file
        (depending on :py:attr:`swarm_mc.config.Config.output_format`).

        Args:
            name: name of the json file created
        """

        if name is None:
            name = '_'.join([self.config.base_name, "temporal_evolution"])

        dataframe = self.time_series.to_dataframe()

        if self._wants_json():
            with open(self._path(name, "json"), "w") as output_file:
                json.dump(dataframe.to_dict('list'), output_file, indent=2)

        if self._wants_csv():
            dataframe.to_csv(self._path(name, "csv"), index=False)

    def save_swarm_parameters(self, name: str = None) -> None:
        """
        Saves the final swarm parameters to a csv and/or json file.

        Args:
            name: name of the json file created
        """

        coefficients = self.swarm_parameters_dict()
        if not coefficients:
            return

        if name is None:
            name = '_'.join([self.config.base_name, "swarm_parameters"])

        if self._wants_json():
            with open(self._path(name, "json"), "w") as output_file:
                json.dump(coefficients, output_file, indent=2)

        if self._wants_csv():
            pd.DataFrame([coefficients]).to_csv(self._path(name, "csv"),
                                                index=False)

    def swarm_parameters_dict(self) -> dict:
        """
        Return swarm parameters as a plain dict for saving or summarising.
        """
        if self.flux.w_err[2] is None:
            return {}

        return {
            'E/N (Td)': float(self.config.EN),
            'mean energy (eV)': self._safe_float(self.energy_distribution.energy_mean),
            'mean energy error (eV)':
                self._safe_float(self.energy_distribution.energy_mean_err),
            'bulk drift velocity (m.s-1)': self._safe_float(self.bulk.w[2]),
            'bulk drift velocity error (m.s-1)':
                self._safe_float(self.bulk.w_err[2]),
            'bulk L diffusion coeff. * N (m-1.s-1)': self._safe_float(self.bulk.DN[2]),
            'bulk L diffusion coeff. error * N (m-1.s-1)':
                self._safe_float(self.bulk.DN_err[2]),
            'bulk T diffusion coeff. * N (m-1.s-1)': self._safe_float(self.bulk.DN[0]),
            'bulk T diffusion coeff. error * N (m-1.s-1)':
                self._safe_float(self.bulk.DN_err[0]),
            'flux drift velocity (m.s-1)': self._safe_float(self.flux.w[2]),
            'flux drift velocity error (m.s-1)': self._safe_float(self.flux.w_err[2]),
            'flux L diffusion coeff. * N (m-1.s-1)': self._safe_float(self.flux.DN[2]),
            'flux L diffusion coeff. error * N (m-1.s-1)':
                self._safe_float(self.flux.DN_err[2]),
            'flux T diffusion coeff. * N (m-1.s-1)': self._safe_float(self.flux.DN[0]),
            'flux T diffusion coeff. error * N (m-1.s-1)':
                self._safe_float(self.flux.DN_err[0]),
            'effective ionization rate coeff. (counted) (m3.s-1)':
                self._safe_float(self.rates_count.effective),
            'effective ionization rate coeff. (counted) error (m3.s-1)':
                self._safe_float(self.rates_count.effective_err),
            'ionization rate coeff. (counted) (m3.s-1)':
                self._safe_float(self.rates_count.ionization),
            'ionization rate coeff. (counted) error (m3.s-1)':
                self._safe_float(self.rates_count.ionization_err),
            'attachment rate coeff. (counted) (m3.s-1)':
                self._safe_float(self.rates_count.attachment),
            'attachment rate coeff. (counted) error (m3.s-1)':
                self._safe_float(self.rates_count.attachment_err),
            'effective ionization rate coeff. (convolution) (m3.s-1)':
                self._safe_float(self.rates_conv.effective),
            'ionization rate coeff. (convolution) (m3.s-1)':
                self._safe_float(self.rates_conv.ionization),
            'attachment rate coeff. (convolution) (m3.s-1)':
                self._safe_float(self.rates_conv.attachment),
        }

    def summary_row(self) -> dict:
        """
        Compact summary of transport coefficients for sweep aggregation.
        """
        return self.swarm_parameters_dict()

    def save_energy_distribution(self, name: str = None) -> None:
        """
        Saves the final time-averaged electron energy distribution function and electron
        energy probability function to a csv and/or json file.

        Args:
            name: name of the json file created
        """

        if self.energy_distribution.eedf is None:
            return

        energy_distributions = pd.DataFrame({
            'energy (eV)': self.energy_distribution.energy_bin_centers,
            'eedf (eV-1)': self.energy_distribution.eedf,
            'eepf (eV-3/2)': self.energy_distribution.eepf
        })

        if name is None:
            name = '_'.join([self.config.base_name, "energy_distribution"])

        if self._wants_json():
            with open(self._path(name, "json"), "w") as output_file:
                json.dump(energy_distributions.to_dict('list'), output_file, indent=2)

        if self._wants_csv():
            energy_distributions.to_csv(self._path(name, "csv"), index=False)

    def save_perf_profile(self, perf_log: list, name: str = None) -> None:
        """
        Save per-step performance breakdown (profiling opt-in).

        Args:
            perf_log (list): list of dicts containing timing data
            name (str): optional base filename
        """

        if not perf_log:
            return
        if name is None:
            name = '_'.join([self.config.base_name, "perf_breakdown"])
        pd.DataFrame(perf_log).to_csv(self._path(name, "csv"), index=False)

    @staticmethod
    def plot_sst_line(axes, t_sst) -> None:
        """
        Plots a vertical line at the time instant where equilibrium is reached.

        Args:
            axes: axes to plot the line on
            t_sst: equilibration time
        """

        y_limits = axes.get_ylim()
        axes.plot([t_sst, t_sst], y_limits, 'k-', label='equilibration time')
        axes.set_ylim(y_limits)

    def plot_temporal_evolution(self, show: bool = True,
                                block: bool = True) -> plt.Figure:
        """
        Produces a figure showing the temporal evolution of the swarm. The figure
        contains five subplots showing the number of particles, the mean electron
        energy, the number of collisions, the mean electron position and the
        variance of electron positions.

        Args:
            show (bool): calls plt.show() if True, else does nothing
            block (bool): if show is True, plt.show(block) is called

        Returns: Matplotlib figure object
        """

        data = self.time_series
        simu_time = data.time
        pos_mean = data.mean_position
        pos_variance = data.var_position

        # plot a line to indicate the time when equilibrium is reached
        t_sst = None
        if self.time_series.ind_equ is not None:
            t_sst = data.time[self.time_series.ind_equ]

        fig = plt.figure()
        ax = fig.add_subplot(231)
        ax.plot(simu_time, data.num_electrons - data.num_electrons[0],
                label='n_electrons - n_0', linewidth=3)
        ax.plot(simu_time, data.num_cations, label='n_cations')
        ax.plot(simu_time, data.num_anions, label='n_anions')
        Output.plot_sst_line(ax, t_sst)
        ax.legend()
        ax.set_title('Number of Particles')
        ax.set_ylabel('number of particles')
        ax.set_xlabel('time (s)')
        ax.grid(True)

        ax = fig.add_subplot(232)
        ax.plot(simu_time, data.mean_energy)
        Output.plot_sst_line(ax, t_sst)
        ax.legend()
        ax.set_title("Mean Energy of Electrons")
        ax.set_ylabel('mean electron energy (eV)')
        ax.set_xlabel('time (s)')
        ax.grid(True)

        ax = fig.add_subplot(233)
        ax.plot(simu_time[1:], np.diff(data.num_collisions), label='collisions')
        ax.plot(simu_time[1:], data.num_electrons[1:] - np.diff(data.num_collisions),
                label='null collisions')
        Output.plot_sst_line(ax, t_sst)
        ax.legend()
        ax.set_title('Number of Collisions Per Simulation Step')
        ax.set_ylabel('number of collisions')
        ax.set_xlabel('time (s)')
        ax.grid(True)

        ax = fig.add_subplot(234)
        ax.plot(simu_time, pos_mean[:, 0], label="x")
        ax.plot(simu_time, pos_mean[:, 1], label="y")
        ax.plot(simu_time, pos_mean[:, 2], label="z")
        Output.plot_sst_line(ax, t_sst)
        ax.legend()
        ax.set_title("Mean Position of Electrons")
        ax.set_ylabel('mean position (m)')
        ax.set_xlabel('time (s)')
        ax.grid(True)

        ax = fig.add_subplot(235)
        ax.plot(simu_time, pos_variance[:, 0], label="x")
        ax.plot(simu_time, pos_variance[:, 1], label="y")
        ax.plot(simu_time, pos_variance[:, 2], label="z")
        Output.plot_sst_line(ax, t_sst)
        ax.legend()
        ax.set_title("Variance of Electrons Positions")
        ax.set_ylabel('variance (m$^2$)')
        ax.set_xlabel('time (s)')
        ax.grid(True)

        if show:
            plt.show(block=block)
        return fig

    def plot_energy_distribution(self, show: bool = True,
                                 block: bool = True) -> plt.Figure:
        """
        Produces a figure showing the time-averaged electron energy distribution
        function and electron energy probability function.

        Args:
            show (bool): calls plt.show() if True, else does nothing
            block (bool): if show is True, plt.show(block) is called

        Returns: Matplotlib figure object
        """

        if self.energy_distribution.eedf is not None:
            energy = self.energy_distribution.energy_bin_centers
            eedf = np.clip(self.energy_distribution.eedf,
                           np.finfo(float).tiny, None)
            eepf = np.clip(self.energy_distribution.eepf,
                           np.finfo(float).tiny, None)

            def _log_limits(arr, floor_ratio=1e-2, padding=1.5):
                pos = np.asarray(arr)
                pos = pos[np.isfinite(pos) & (pos > 0)]
                if pos.size == 0:
                    return (1e-30, 1)
                lo = pos.min()
                hi = pos.max()
                lo = max(lo, hi * floor_ratio)
                return (lo / padding, hi * padding)

            fig, ax1 = plt.subplots()
            plt.title("Electron Energy Distribution Function (EEDF) \n"
                      "and Electron Energy Probability Function (EEPF)")
            color = 'tab:blue'
            ax1.plot(energy, eedf, color=color)
            ax1.set_yscale('log')
            ax1.set_ylim(_log_limits(self.energy_distribution.eedf))
            ax1.tick_params(axis='y', labelcolor=color)
            ax1.set_ylabel('EEDF (eV$^{-1}$)', color=color)
            ax1.set_xlabel('energy (eV)')
            ax1.grid(True)

            ax2 = ax1.twinx()
            color = 'tab:red'
            ax2.plot(energy, eepf, color=color)
            ax2.set_yscale('log')
            ax2.set_ylim(_log_limits(self.energy_distribution.eepf))
            ax2.tick_params(axis='y', labelcolor=color)
            ax2.set_ylabel('EEPF (eV$^{-3/2}$)', color=color)

            if show:
                plt.show(block=block)
            return fig
