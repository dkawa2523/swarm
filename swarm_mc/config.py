#  Copyright (c) 2020-2021 ETH Zurich

"""
Module for the Config class of the simulation.
"""

# Import Packages
import os
from typing import Tuple, Union
import warnings
from collections.abc import Callable
import json
import json5
import yaml
import numpy as np
import scipy.constants as csts

num = Union[int, float]


class Config:
    """
    Configuration class for the Monte-Carlo simulation.

    The configuration can be loaded from, or saved to, json, json5 or yaml files.
    Alternatively, it can be provided as, or exported to, a (nested) dictionary.
    The gas number density is not a configuration parameter, but a
    cached property of the Config class, which is computed from the pressure and
    temperature.

    Attributes:
        paths_to_cross_section_files (list): paths to the cross section files in
            txt format
        gases (list): sum formulae of gases
        fractions (list): proportions of the gases in the gas mixture
        max_cross_section_energy (float): maximum cross section energy (eV)

        output_directory (str): path to the output directory
        base_name (str): prefix of the output filename
        save_simulation_pickle (bool): save the simulation as pickle file
        save_temporal_evolution (bool): save temporal evolution
        save_swarm_parameters (bool): save swarm parameters
        save_energy_distribution (bool): save energy distribution
        output_format (str): persistence format: ``csv``, ``json`` or ``both``
            (case-insensitive). Defaults to ``csv``.

        EN (float): E/N ratio in (Td)
        _pressure (float): gas pressure in Pa
        _temperature (float): gas temperature in K
        _gas_number_density (float): gas number density in m-3

        num_e_initial (int): initial number of electrons
        initial_pos_electrons (list): initial position [x, y, z] of the electrons'
            center of mass
        initial_std_electrons (list): initial broadening of gaussian distributed
            electrons in x, y and z direction
        initial_energy_distribution (str): The initial energy distribution of
            the electrons. Can be either ``"zero"`` (all electrons have zero
            kinetic energy), ``"fixed"`` (all electrons have the same energy)
            or ``"maxwell-boltzmann"`` (at temperature :py:attr:`initial_temperature`).
            Maxwell-Boltzmann support is experimental, check
            :py:func:`swarm_mc.utils.maxwell_boltzmann_random` and its test case
            to see if the required precision is achieved.
            The default is ``"zero"``. If the initial distribution is `"fixed"`,
            :py:attr:`initial_energy` and :py:attr:`initial_direction` must be set.
        initial_energy (float): The initial energy of the electrons in eV.
        initial_direction (Union[Tuple[float, float, float], str]): The initial
            direction of the electrons. Either the string ``"random"`` to give
            each electron a random direction or a tuple with three elements x,
            y, z specifying a single direction for all electrons.
        initial_temperature (Union[Float, int]): The initial temperature in K. Used
            for the Maxwell-Boltzmann distribution.

        num_energy_bins (int): number of energy bins to group the electrons for the
            energy distribution
        energy_sharing_factor (float): energy sharing factor for ionization collisions
        energy_sharing_floor (float): optional minimum sharing factor applied at low
            incident energies to slow chain ionization; enabled when set
        energy_sharing_taper_energy (float): characteristic energy (eV) over which the
            sharing factor blends from floor to nominal value
        isotropic_scattering (bool): scattering: isotropic (true), non-isotropic
            according to Vahedi et al. (false)
        use_velocity_lut (bool): if True, precompute velocity_from_energy on the
            cross-section grid and use LUT interpolation to avoid repeated sqrt
        use_scatter_lut (bool): if True, use a pre-tabulated inverse CDF for cos_chi
            to reduce per-step trig cost in isotropic scattering
        use_aniso_scatter_lut (bool): if True, pre-tabulate the anisotropic cos_chi
            inverse CDF over an energy grid and sample via LUT (fallback to analytic)
        use_max_coll_lut (bool): if True, compress max-collision-frequency tables to a
            smaller LUT for faster timestep interpolation
        max_coll_lut_size (int): target size of compressed LUT (<= original grid size)
        use_jit_energy_loss (bool): if True, use numba energy loss kernel; can disable
            to match NumPy path exactly
        use_jit_collision_classify (bool): if True, classify attachment/ionization via
            numba; disable to match NumPy path
        jit_parallel_threshold (int): minimum particle count to use prange-based JIT
            kernels; below this falls back to NumPy for lower overhead
        conserve (bool): conservation of the number of electrons
        num_e_max (int): maximum allowed electron number (when it is reached, the
            number of electrons is then conserved until simulation ends)
        max_ionization_growth_factor (float): optional cap on new electrons created
            per step, expressed as a multiple of the current population
        num_e_throttle_ratio (float): fraction of num_e_max at which ionization
            throttling begins (linear fade to zero at num_e_max)
        num_e_throttle_min_factor (float): lower bound on throttling factor
        timestep_min (float): optional lower bound on dt (s)
        timestep_max (float): optional upper bound on dt (s)
        scattering_transition_low (float): lower energy bound (eV) to start blending
            toward anisotropic scattering
        scattering_transition_high (float): upper energy bound (eV) to finish
            blending toward anisotropic scattering
        sst_min_window (int): minimum samples for steady-state regression window
        sst_max_window (int): maximum samples for steady-state regression window
        sst_energy_slope_tol (float): relative slope tolerance for energy trend
        sst_energy_resid_tol (float): relative residual std tolerance for energy trend
        sst_transport_slope_tol (float): relative slope tolerance for transport trends
        sst_transport_resid_tol (float): relative residual std tolerance for transport
        dn_block_size (int): block size for DN/w estimation and jackknife error
        dn_step_skip (int): optional sampling stride for velocity moments (>=1)
        dn_use_jackknife (bool): if True, enable block jackknife for DN/w errors
        log_every_k (int): optional stride to downsample TimeSeries logging (>=1)
        seed (int, str): optional. If set to an integer it is used to seed
            the Simulation. If set to `"random"` entropy is drawn from the OS.
            If set to `"run_id"` a deterministic seed derived from the run label
            is used. Default value is `"random"`.

        end_condition_type (str): Specifies the end condition. Can be
            ``"steady-state"``, ``"num_col_max"``, ``"w_tol+ND_tol"`` or
            ``"custom"``.  The ``"custom"`` end condition requires
            :py:attr:`is_done` to be set as well.  Defaults to ``"w_tol+ND_tol"``
        w_tol (float): tolerance on the flux drift velocity. simulation ends
            when w_err/w < w_tol
        DN_tol (float): tolerance on the flux diffusion coefficient. simulation ends
            when DN_err/w < DN_tol
        num_col_max (int): maximum number of collisions during the simulation,
            simulation ends when it is reached
        is_done (Callable): This function gets called to determine whether to end the
            simulation or not. Gets passed the simulation object as argument.
            Return ``True`` to stop the simulation, ``False`` otherwise.
        timeout (int): End the simulation after ``timeout`` seconds. Zero means no
            timeout. Defaults to zero.
    """

    def __init__(self, config: Union[str, dict]):
        """
        Instantiate the config.

        Args:
            config (str, dict): path to a json or json5 config file, or dictionary.
        """

        def _as_int(value, field_name: str) -> int:
            try:
                return int(value)
            except (TypeError, ValueError):
                try:
                    return int(float(value))
                except (TypeError, ValueError):
                    raise ValueError(f"{field_name} must be an integer-like value")

        if isinstance(config, str):
            if config.endswith('.json5'):
                with open(config, "r") as json_file:
                    config = json5.load(json_file)
            elif config.endswith('.json'):
                with open(config, "r") as json_file:
                    config = json.load(json_file)
            elif config.endswith('.yaml') or config.endswith('.yml'):
                with open(config, "r") as yaml_file:
                    config = yaml.safe_load(yaml_file)
            else:
                raise ValueError(f"Configuration file '{config}' has invalid extension."
                                 " Extensions '.json', '.json5', '.yaml' or '.yml' are"
                                 " expected.")

        if isinstance(config, dict) and 'base_config' in config:
            # allow passing full experiment YAML/JSON and extracting base_config
            config = config['base_config']

        # gases
        input_gases = config['input_gases']
        self.paths_to_cross_section_files: list = \
            input_gases['paths_to_cross_section_files']
        self.gases: list = input_gases['gases']
        self.fractions: list = input_gases['fractions']
        self.max_cross_section_energy: float = \
            float(input_gases['max_cross_section_energy'])
        if self.max_cross_section_energy <= 0:
            raise ValueError("max_cross_section_energy must be positive")
        if not np.isclose(sum(self.fractions), 1, atol=1e-3):
            raise ValueError("Fractions must sum to 1 within ±1e-3.")
        if any(f < 0 for f in self.fractions):
            raise ValueError("Fractions cannot be negative.")

        # output
        output = config['output']
        self.output_directory: str = output['output_directory']
        if not self.output_directory.endswith(os.sep):
            self.output_directory += os.sep
        self.base_name: str = output['base_name']
        self.save_simulation_pickle: bool = output['save_simulation_pickle']
        self.save_temporal_evolution: bool = output['save_temporal_evolution']
        self.save_swarm_parameters: bool = output['save_swarm_parameters']
        self.save_energy_distribution: bool = output['save_energy_distribution']
        self.output_format: str = output.get('output_format', 'csv').lower()
        if self.output_format not in ('csv', 'json', 'both'):
            raise ValueError("output_format must be one of: csv, json, both")

        # physical conditions
        physical_conditions = config['physical_conditions']
        self.EN: float = float(physical_conditions['EN'])
        self._pressure: float = float(physical_conditions['pressure'])
        self._temperature: float = float(physical_conditions['temperature'])
        self._gas_number_density: float = None
        if self.EN < 0:
            raise ValueError("EN must be non-negative")
        if self._pressure <= 0:
            raise ValueError("pressure must be > 0")
        if self._temperature <= 0:
            raise ValueError("temperature must be > 0")

        # initial state
        initial_state = config['initial_state']
        self.num_e_initial: int = _as_int(initial_state['num_e_initial'],
                                          "num_e_initial")
        self.initial_pos_electrons: list = initial_state['initial_pos_electrons']
        self.initial_std_electrons: list = initial_state['initial_std_electrons']
        self.initial_energy_distribution: str = "zero"
        self.initial_energy: float = None
        self.initial_direction: Union[Tuple[float, float, float], str] = None
        self.initial_temperature: Union[int, float] = None

        if 'initial_energy_distribution' in initial_state:
            val = initial_state['initial_energy_distribution']

            if val not in ("zero", "fixed", "maxwell-boltzmann"):
                raise ValueError(
                    "initial_energy_distribution must be zero, fixed or "
                    "maxwell-boltzmann")

            self.initial_energy_distribution = val

        if 'initial_energy' in initial_state and \
                initial_state['initial_energy'] is not None:
            self.initial_energy = float(initial_state['initial_energy'])
            if self.initial_energy < 0:
                raise ValueError("initial_energy cannot be negative")

        if 'initial_direction' in initial_state and \
                initial_state['initial_direction'] is not None:
            val = initial_state['initial_direction']

            if isinstance(val, str) and val == "random":
                self.initial_direction = val
            elif isinstance(val, (list, tuple)):
                if len(val) != 3:
                    raise ValueError("initial_direction must be \"random\" "
                                     "or list of three floats")

                self.initial_direction = [float(item) for item in val]

                if self.initial_direction == [0, 0, 0]:
                    raise ValueError("initial_direction cannot be all zero")
            else:
                raise ValueError("Invalid value for initial_direction")

        if 'initial_temperature' in initial_state:
            self.initial_temperature = float(initial_state['initial_temperature'])

            if self.initial_temperature < 0:
                raise ValueError("initial_temperature must be positive")

        if self.initial_energy_distribution in ("zero", "maxwell-boltzmann"):
            if self.initial_energy is not None:
                warnings.warn("initial_energy setting useless with "
                              f"{self.initial_energy_distribution} distribution")
            if self.initial_direction is not None:
                warnings.warn("initial_direction setting useless with "
                              f"{self.initial_energy_distribution} distribution")
        else:
            assert self.initial_energy_distribution == "fixed"
            if self.initial_energy is None:
                raise ValueError("Must set initital_energy")
            if self.initial_direction is None:
                raise ValueError("Must set initial_direction")

        if self.initial_energy_distribution == "maxwell-boltzmann" and \
                self.initial_temperature is None:
            raise ValueError("Must set initial_temperature for Maxwell-Boltzmann")

        # simulation settings
        simulation = config['simulation_settings']
        self.num_energy_bins: int = simulation['num_energy_bins']
        self.energy_sharing_factor: float = float(simulation['energy_sharing_factor'])
        self.energy_sharing_floor = simulation.get('energy_sharing_floor')
        if self.energy_sharing_floor is not None:
            self.energy_sharing_floor = float(self.energy_sharing_floor)
            if not 0 < self.energy_sharing_floor <= 1:
                raise ValueError("energy_sharing_floor must be in (0, 1]")
            if self.energy_sharing_floor > self.energy_sharing_factor:
                warnings.warn("energy_sharing_floor exceeds energy_sharing_factor; "
                              "sharing will clamp to the higher floor.")
        self.energy_sharing_taper_energy: float = float(
            simulation.get('energy_sharing_taper_energy', 10.0))
        if self.energy_sharing_taper_energy <= 0:
            raise ValueError("energy_sharing_taper_energy must be positive")
        self.isotropic_scattering: bool = simulation['isotropic_scattering']
        self.use_velocity_lut: bool = bool(simulation.get('use_velocity_lut', True))
        self.use_scatter_lut: bool = bool(simulation.get('use_scatter_lut', True))
        self.use_aniso_scatter_lut: bool = bool(
            simulation.get('use_aniso_scatter_lut', False)
        )
        self.use_max_coll_lut: bool = bool(simulation.get('use_max_coll_lut', True))
        self.max_coll_lut_size: int = int(simulation.get('max_coll_lut_size', 256))
        if self.max_coll_lut_size < 16:
            raise ValueError("max_coll_lut_size must be >= 16")
        self.use_jit_energy_loss: bool = bool(simulation.get('use_jit_energy_loss', True))
        self.use_jit_collision_classify: bool = bool(simulation.get('use_jit_collision_classify', True))
        self.jit_parallel_threshold: int = int(simulation.get('jit_parallel_threshold', 2000))
        if self.jit_parallel_threshold < 1:
            raise ValueError("jit_parallel_threshold must be >= 1")
        self.conserve: bool = simulation['conserve']
        self.num_e_max: int = _as_int(simulation['num_e_max'], "num_e_max")
        self.max_ionization_growth_factor = \
            simulation.get('max_ionization_growth_factor', None)
        if self.max_ionization_growth_factor is not None:
            self.max_ionization_growth_factor = \
                float(self.max_ionization_growth_factor)
            if self.max_ionization_growth_factor <= 0:
                raise ValueError("max_ionization_growth_factor must be > 0 "
                                 "or null to disable")
        self.num_e_throttle_ratio = simulation.get('num_e_throttle_ratio', None)
        if self.num_e_throttle_ratio is not None:
            self.num_e_throttle_ratio = float(self.num_e_throttle_ratio)
            if not 0 < self.num_e_throttle_ratio < 1:
                raise ValueError("num_e_throttle_ratio must be in (0, 1) or null")
        self.num_e_throttle_min_factor = float(
            simulation.get('num_e_throttle_min_factor', 0.0))
        if self.num_e_throttle_min_factor < 0 or self.num_e_throttle_min_factor > 1:
            raise ValueError("num_e_throttle_min_factor must be in [0, 1]")
        self.timestep_min = simulation.get('timestep_min')
        if self.timestep_min is not None:
            self.timestep_min = float(self.timestep_min)
            if self.timestep_min <= 0:
                raise ValueError("timestep_min must be > 0 or null")
        self.timestep_max = simulation.get('timestep_max')
        if self.timestep_max is not None:
            self.timestep_max = float(self.timestep_max)
            if self.timestep_max <= 0:
                raise ValueError("timestep_max must be > 0 or null")
        if self.timestep_min and self.timestep_max \
                and self.timestep_min > self.timestep_max:
            raise ValueError("timestep_min cannot exceed timestep_max")
        self.scattering_transition_low = float(
            simulation.get('scattering_transition_low', 1.0))
        self.scattering_transition_high = float(
            simulation.get('scattering_transition_high', 5.0))
        if self.scattering_transition_low <= 0 or \
                self.scattering_transition_high <= self.scattering_transition_low:
            raise ValueError("scattering_transition_* must satisfy 0<low<high")
        self.sst_min_window = int(simulation.get('sst_min_window', 30))
        self.sst_max_window = int(simulation.get('sst_max_window', 50))
        if self.sst_min_window <= 5 or self.sst_max_window < self.sst_min_window:
            raise ValueError("sst_min_window must be >5 and <= sst_max_window")
        self.sst_energy_slope_tol = float(
            simulation.get('sst_energy_slope_tol', 0.02))
        self.sst_energy_resid_tol = float(
            simulation.get('sst_energy_resid_tol', 0.08))
        self.sst_transport_slope_tol = float(
            simulation.get('sst_transport_slope_tol', 0.02))
        self.sst_transport_resid_tol = float(
            simulation.get('sst_transport_resid_tol', 0.1))
        self.dn_block_size = int(simulation.get('dn_block_size', 20))
        if self.dn_block_size < 5:
            raise ValueError("dn_block_size must be >= 5")
        self.dn_step_skip = int(simulation.get('dn_step_skip', 1))
        if self.dn_step_skip < 1:
            raise ValueError("dn_step_skip must be >= 1")
        self.dn_use_jackknife = bool(simulation.get('dn_use_jackknife', False))
        self.log_every_k = int(simulation.get('log_every_k', 1))
        if self.log_every_k < 1:
            raise ValueError("log_every_k must be >= 1")
        self.use_jit = bool(simulation.get('use_jit', False))
        self.profile_per_step = bool(simulation.get('profile_per_step', False))
        self.seed: Union[int, str] = "random"
        if 'seed' in simulation:
            seed = simulation['seed']
            if isinstance(seed, int) \
                    or (isinstance(seed, str)
                        and seed in ("random", "run_id")):
                self.seed = seed
            else:
                raise ValueError("seed must be an integer or one of: \"random\", "
                                 "\"run_id\"")

        # end conditions
        end_conditions = config['end_conditions']
        self.end_condition_type: str = "w_tol+ND_tol"

        if 'end_condition_type' in end_conditions:
            val = str(end_conditions['end_condition_type'])

            if val not in ("steady-state", "num_col_max", "w_tol+ND_tol", "custom"):
                raise ValueError("end_condition_type must be \"steady-state\", "
                                 "\"num_col_max\", \"w_tol+ND_tol\" or \"custom\"")

            self.end_condition_type = val

        self.w_tol: float = float(end_conditions['w_tol'])
        self.DN_tol: float = float(end_conditions['DN_tol'])
        self.num_col_max: int = _as_int(end_conditions['num_col_max'], "num_col_max")
        self.is_done: Callable = None
        self.timeout: int = 0

        if 'is_done' in end_conditions and end_conditions['is_done'] is not None:
            self.is_done = end_conditions['is_done']

            if not isinstance(self.is_done, Callable):
                raise TypeError("is_done must be a function")

        if self.end_condition_type == "custom" and self.is_done is None:
            raise ValueError("custom requires is_done to be set to a callback")
        elif self.end_condition_type != "custom" and self.is_done is not None:
            warnings.warn("Setting is_done is useless without end_condition_type "
                          "custom")

        if 'timeout' in end_conditions:
            self.timeout = int(end_conditions['timeout'])

            if self.timeout < 0:
                raise ValueError("timeout must be >= 0")

    @property
    def gas_number_density(self) -> float:
        if self._gas_number_density is None:
            self._gas_number_density = \
                self.pressure / (csts.Boltzmann * self.temperature)
        return self._gas_number_density

    @property
    def pressure(self) -> float:
        return self._pressure

    @pressure.setter
    def pressure(self, value: float):
        """
        Pressure setter. If a new value is set, resets the cache for the gas
        number density.

        Args:
            value: pressure in Pascal
        """
        self._pressure = value
        self._gas_number_density = None

    @property
    def temperature(self) -> float:
        return self._temperature

    @temperature.setter
    def temperature(self, value: float) -> None:
        """
        Temperature setter. If a new value is set, resets the cache for the gas
        number density.

        Args:
            value: temperature in Kelvin
        """
        self._temperature = value
        self._gas_number_density = None

    def to_dict(self) -> dict:
        """
        Returns the current configuration as a dictionary.

        Returns: dict of configuration
        """

        return {
            'input_gases': {
                'gases': self.gases,
                'paths_to_cross_section_files': self.paths_to_cross_section_files,
                'fractions': self.fractions,
                'max_cross_section_energy': self.max_cross_section_energy,
            },
            'output': {
                'output_directory': self.output_directory,
                'base_name': self.base_name,
                'save_simulation_pickle': self.save_simulation_pickle,
                'save_temporal_evolution': self.save_temporal_evolution,
                'save_swarm_parameters': self.save_swarm_parameters,
                'save_energy_distribution': self.save_energy_distribution,
                'output_format': self.output_format,
            },
            'physical_conditions': {
                'EN': self.EN,
                'pressure': self.pressure,
                'temperature': self.temperature,
            },
            'initial_state': {
                'num_e_initial': self.num_e_initial,
                'initial_pos_electrons': self.initial_pos_electrons,
                'initial_std_electrons': self.initial_std_electrons,
                'initial_energy_distribution': self.initial_energy_distribution,
                'initial_energy': self.initial_energy,
                'initial_direction': self.initial_direction,
                'initial_temperature': self.initial_temperature,
            },
            'simulation_settings': {
                'num_energy_bins': self.num_energy_bins,
                'energy_sharing_factor': self.energy_sharing_factor,
                'energy_sharing_floor': self.energy_sharing_floor,
                'energy_sharing_taper_energy': self.energy_sharing_taper_energy,
                'isotropic_scattering': self.isotropic_scattering,
                'use_velocity_lut': self.use_velocity_lut,
                'use_scatter_lut': self.use_scatter_lut,
                'use_aniso_scatter_lut': self.use_aniso_scatter_lut,
                'use_max_coll_lut': self.use_max_coll_lut,
                'max_coll_lut_size': self.max_coll_lut_size,
                'use_jit_energy_loss': self.use_jit_energy_loss,
                'use_jit_collision_classify': self.use_jit_collision_classify,
                'jit_parallel_threshold': self.jit_parallel_threshold,
                'conserve': self.conserve,
                'num_e_max': self.num_e_max,
                'max_ionization_growth_factor': self.max_ionization_growth_factor,
                'num_e_throttle_ratio': self.num_e_throttle_ratio,
                'num_e_throttle_min_factor': self.num_e_throttle_min_factor,
                'timestep_min': self.timestep_min,
                'timestep_max': self.timestep_max,
                'scattering_transition_low': self.scattering_transition_low,
                'scattering_transition_high': self.scattering_transition_high,
                'sst_min_window': self.sst_min_window,
                'sst_max_window': self.sst_max_window,
                'sst_energy_slope_tol': self.sst_energy_slope_tol,
                'sst_energy_resid_tol': self.sst_energy_resid_tol,
                'sst_transport_slope_tol': self.sst_transport_slope_tol,
                'sst_transport_resid_tol': self.sst_transport_resid_tol,
                'dn_block_size': self.dn_block_size,
                'dn_step_skip': self.dn_step_skip,
                'dn_use_jackknife': self.dn_use_jackknife,
                'log_every_k': self.log_every_k,
                'use_jit': self.use_jit,
                'profile_per_step': self.profile_per_step,
                'seed': self.seed,
            },
            'end_conditions': {
                'end_condition_type': self.end_condition_type,
                'w_tol': self.w_tol,
                'DN_tol': self.DN_tol,
                'num_col_max': self.num_col_max,
                'is_done': self.is_done,
                'timeout': self.timeout,
            }
        }

    def save_json5(self, path: str = 'config.json5') -> None:
        """
        Saves the current configuration to a json5 file.

        Args:
            path (str): path including the file name and extension,
                example: 'data/config.json5'
        """

        d = self.to_dict()

        if d['end_conditions']['is_done'] is not None:
            del d['end_conditions']['is_done']
            warnings.warn("Cannot save custom callback to json5!")

        with open(path, "w") as config_file:
            json5.dump(d, config_file, indent=2)

    def save_json(self, path: str = 'config.json') -> None:
        """
        Saves the current configuration to a json file.

        Args:
            path (str): path including the file name and extension,
                example: 'data/config.json'
        """

        d = self.to_dict()

        if d['end_conditions']['is_done'] is not None:
            del d['end_conditions']['is_done']
            warnings.warn("Cannot save custom callback to json!")

        with open(path, "w") as config_file:
            json.dump(d, config_file, indent=2)
