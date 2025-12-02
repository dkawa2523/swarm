#  Copyright (c) 2020-2021 ETH Zurich

"""
Module for the TimeSeries class.
"""

# Import Packages
import numpy as np
import pandas as pd

# Import modules
from swarm_mc.electrons import Electrons


class TimeSeries:
    """
    Temporal evolution of some quantities during the simulation.

    The TimeSeries class stores the temporal evolution of some quantities in arrays of
    ever-increasing length. The TimeSeries can be exported to a pandas DataFrame for
    further processing.

    Attributes:
        ind_equ (int): index of equilibration time
        time (ndarray): simulated time
        num_collisions (ndarray): cumulative number of collisions
        num_electrons (ndarray): number of electrons
        num_cations (ndarray): cumulative number of cations
        num_anions (ndarray): cumulative number of anions
        mean_energy (ndarray): mean energy of electrons
        mean_position (ndarray): mean position of electrons
        var_position (ndarray): variance of electrons positions
        mean_velocity (ndarray): mean velocity of electrons
        mean_velocity_moment (ndarray): mean velocity moment of electrons
    """

    def __init__(self, electrons: Electrons, log_every_k: int = 1):
        """
        Instantiates a TimeSeries.

        Args:
            electrons (Electrons): electron related data
            log_every_k (int): optional stride to downsample logging
        """

        # number of time steps until equilibration
        self.ind_equ = None
        self.log_every_k = max(1, int(log_every_k))
        self._step_counter = 0

        self._time = [0.0]
        self._num_collisions = [0]
        self._num_electrons = [electrons.num_e]
        self._num_cations = [0]
        self._num_anions = [0]
        self._mean_energy = [electrons.mean_energy]
        self._mean_position = [electrons.mean_position.flatten()]
        self._var_position = [electrons.var_position.flatten()]
        self._mean_velocity = [electrons.mean_velocity.flatten()]
        self._mean_velocity_moment = [electrons.mean_velocity_moment.flatten()]
        self._cache = {}
        self.dn_block_size = None

    def _invalidate_cache(self) -> None:
        self._cache = {}

    def _array(self, name: str, data) -> np.ndarray:
        if name in self._cache:
            return self._cache[name]
        arr = np.asarray(data, dtype=float)
        self._cache[name] = arr
        return arr

    @property
    def time(self) -> np.ndarray:
        return self._array("time", self._time)

    @property
    def num_collisions(self) -> np.ndarray:
        return self._array("num_collisions", self._num_collisions)

    @property
    def num_electrons(self) -> np.ndarray:
        return self._array("num_electrons", self._num_electrons)

    @property
    def num_cations(self) -> np.ndarray:
        return self._array("num_cations", self._num_cations)

    @property
    def num_anions(self) -> np.ndarray:
        return self._array("num_anions", self._num_anions)

    @property
    def mean_energy(self) -> np.ndarray:
        return self._array("mean_energy", self._mean_energy)

    @property
    def mean_position(self) -> np.ndarray:
        return self._array("mean_position", self._mean_position)

    @property
    def var_position(self) -> np.ndarray:
        return self._array("var_position", self._var_position)

    @property
    def mean_velocity(self) -> np.ndarray:
        return self._array("mean_velocity", self._mean_velocity)

    @property
    def mean_velocity_moment(self) -> np.ndarray:
        return self._array("mean_velocity_moment", self._mean_velocity_moment)

    def append_data(self, electrons: Electrons, dt: float, n_coll: int, n_cations: int,
                    n_anions: int) -> None:
        """
        Appends latest data to the arrays storing the temporal evolution of some
        quantities.

        Args:
            electrons (Electrons): current data related to the electrons
            dt (float): duration of the current time-step (s)
            n_coll (int): number of collisions during the current time-step
            n_cations (int): number of cations produced during the current time-step
            n_anions (int): number of anions produced during the current time-step
        """

        self._step_counter += 1
        if (self._step_counter - 1) % self.log_every_k != 0:
            return

        self._time.append(self._time[-1] + dt)
        self._num_collisions.append(self._num_collisions[-1] + n_coll)
        self._num_electrons.append(electrons.num_e)
        self._num_cations.append(self._num_cations[-1] + n_cations)
        self._num_anions.append(self._num_anions[-1] + n_anions)
        self._mean_energy.append(electrons.mean_energy)
        self._mean_position.append(electrons.mean_position.flatten())
        self._var_position.append(electrons.var_position.flatten())
        self._mean_velocity.append(electrons.mean_velocity.flatten())
        self._mean_velocity_moment.append(electrons.mean_velocity_moment.flatten())
        self._invalidate_cache()

    def to_dataframe(self) -> pd.DataFrame:
        """
        Creates a pandas DataFrame with the data containe in the TimeSeries and
        returns it.

        Returns: pandas DataFrame containing the TimeSeries data.
        """

        n = self.time.size
        num_coll = np.diff(self.num_collisions, prepend=0)
        num_null_coll = self.num_electrons - num_coll

        df = pd.DataFrame(data=np.concatenate([
            self.time.reshape((n, 1)),
            num_coll.reshape((n, 1)),
            num_null_coll.reshape((n, 1)),
            self.num_electrons.reshape((n, 1)),
            self.num_cations.reshape((n, 1)),
            self.num_anions.reshape((n, 1)),
            self.mean_energy.reshape((n, 1)),
            self.mean_position,
            self.var_position,
            self.mean_velocity,
            self.mean_velocity_moment
        ], axis=1), columns=[
            'time',
            'num_collisions',
            'num_null_collisions',
            'num_electrons',
            'num_cations',
            'num_anions',
            'mean_energy',
            'mean_position_x',
            'mean_position_y',
            'mean_position_z',
            'var_position_x',
            'var_position_y',
            'var_position_z',
            'mean_velocity_x',
            'mean_velocity_y',
            'mean_velocity_z',
            'mean_velocity_moment_x',
            'mean_velocity_moment_y',
            'mean_velocity_moment_z',
        ])
        df.attrs = {'index_equilibration': self.ind_equ}
        return df
