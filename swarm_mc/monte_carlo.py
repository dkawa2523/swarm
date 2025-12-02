#  Copyright (c) 2020-2021 ETH Zurich

"""
Module for the MonteCarlo class.

The MonteCarlo class implements all random-number based methods to simulate the motion
of electrons. The simulation time-step is determined with the null collision technique
(calculate_max_coll_freq creates a lookup table for the choice of the
trial collision frequency, determine_timestep calculates the time-step based on the
current maximum electron energy and acceleration). The collision processes are randomly
chosen based on the collision frequency of each process. The scattering angles are
randomly chosen with either an isotropic or an anisotropic model.
"""

# Import Packages
from typing import Tuple
import numpy as np
import numpy.linalg as linalg
import scipy.constants as csts
import scipy.interpolate

# Import modules
import swarm_mc.utils as utils
from swarm_mc.__about__ import __version__
from swarm_mc.config import Config
from swarm_mc.gas_mixture import GasMixture
from swarm_mc import jit_kernels

np.seterr(all='raise')


class MonteCarlo:
    """
    Class implementing all random-number based simulation methods

    Attributes:
        config (Config): configuration of the simulation
        trial_coll_freq (float): trial collision frequency for the null collision
            technique
        max_coll_freq (interp1d): cumulative maximum of the collision frequency as a
            function of the electron energy
        max_coll_period (ndarray): array inversely proportional to the cumulative
            maximum of the collision frequency
        max_coll_period_squared (ndarray): array inversely proportional to the square
            of the cumulative maximum of the collision frequency
        collision_by_electron (ndarray): array of collision index for each electron
            (starting at 0), an index equal to the number of cross sections indicates
            a null collision
    """

    version = f"swarm_mc version {__version__}\n"

    def __init__(self, cfg: Config, rng=None):
        """
        Instantiates the MonteCarlo class.

        Args:
            cfg (Config): configuration of the simulation
            rng: numpy.random.Generator for reproducibility
        """

        self.config = cfg
        self.rng = rng or np.random.default_rng()

        # trial collision frequency
        self.trial_coll_freq = None
        # data for the calculation of the trial collision frequency
        self.max_coll_freq = None
        self.max_coll_period = None
        self.max_coll_period_squared = None
        self._max_coll_energy_grid = None
        self._fast_energy_grid = None
        self._fast_coll_freq = None
        self._fast_coll_period = None
        self._fast_coll_period_sq = None
        self._jit_enabled = cfg.use_jit and jit_kernels.jit_supported()
        self._vel_lut_grid = None
        self._vel_lut_values = None
        if self._jit_enabled:
            self._velocity_from_energy = self._jit_velocity
        elif cfg.use_velocity_lut:
            self._velocity_from_energy = self._lut_velocity
        else:
            self._velocity_from_energy = utils.velocity_from_energy
        self._rng_state = None
        if self._jit_enabled:
            seed_val = np.uint64(self.rng.integers(1, np.iinfo(np.uint64).max))
            self._rng_state = np.array([seed_val], dtype=np.uint64)
        self._cos_chi_lut = None
        if self.config.use_scatter_lut and self.config.isotropic_scattering:
            self._cos_chi_lut = jit_kernels.cos_chi_inverse_cdf_lut()
        self._aniso_lut = None
        self._aniso_e_grid = None
        if self.config.use_aniso_scatter_lut and not self.config.isotropic_scattering:
            # reuse max collision energy grid for LUT
            self._aniso_e_grid = None  # will be built when max_coll_freq is ready
        self._jit_threshold = max(1, int(self.config.jit_parallel_threshold))

        # vector of which collision which electron undergoes
        self.collision_by_electron = None
        self._velocity_from_energy = self._jit_velocity if self._jit_enabled else utils.velocity_from_energy

    def calculate_max_coll_freq(self, gas_mixture: GasMixture):
        """
        Calculates the maximum collision frequency in the given gas mixture.

        Args:
            gas_mixture (GasMixture): gas mixture
        """

        gas_density = self.config.gas_number_density
        energy_grid = gas_mixture.energy_vector
        velocity = utils.velocity_from_energy(energy_grid)
        freq = gas_density * velocity * gas_mixture.total_cross_section
        if freq.size > 1 and freq[0] == 0:
            freq[0] = freq[1]
        freq = np.maximum.accumulate(freq)
        self._max_coll_energy_grid = energy_grid
        if self.config.use_velocity_lut:
            self._vel_lut_grid = energy_grid
            self._vel_lut_values = velocity
        if self.config.use_aniso_scatter_lut and not self.config.isotropic_scattering:
            self._aniso_e_grid = energy_grid
            self._aniso_lut = jit_kernels.build_aniso_cos_chi_lut(energy_grid)
        self.max_coll_freq = freq
        self.max_coll_period = 0.5 * csts.electron_mass / csts.elementary_charge / freq
        self.max_coll_period_squared = 0.5 * csts.electron_mass \
            / csts.elementary_charge / freq ** 2
        if self.config.use_max_coll_lut:
            target = min(int(self.config.max_coll_lut_size), freq.size)
            if target < freq.size:
                e_min = float(energy_grid[0])
                e_max = float(energy_grid[-1])
                compressed_e = np.linspace(e_min, e_max, target)
                compressed_freq = np.interp(compressed_e, energy_grid, freq)
                self._fast_energy_grid = compressed_e
                self._fast_coll_freq = compressed_freq
                self._fast_coll_period = 0.5 * csts.electron_mass / csts.elementary_charge / compressed_freq
                self._fast_coll_period_sq = 0.5 * csts.electron_mass \
                    / csts.elementary_charge / compressed_freq ** 2
            else:
                self._fast_energy_grid = energy_grid
                self._fast_coll_freq = freq
                self._fast_coll_period = self.max_coll_period
                self._fast_coll_period_sq = self.max_coll_period_squared

    def determine_timestep(self, max_velocity: float, max_acceleration: float) -> float:
        """
        Determine the duration of the next time-step in the simulation with the
        null-collision technique.

        Args:
            max_velocity (float): current maximum electron velocity
            max_acceleration (float): current maximum electron acceleration

        Returns: time-step duration (s)
        """

        rand = float(self._rand_scalar())
        freq = self._fast_coll_freq if (self.config.use_max_coll_lut and self._fast_coll_freq is not None) else self.max_coll_freq
        period = self._fast_coll_period if (self.config.use_max_coll_lut and self._fast_coll_period is not None) else self.max_coll_period
        period_sq = self._fast_coll_period_sq if (self.config.use_max_coll_lut and self._fast_coll_period_sq is not None) else self.max_coll_period_squared
        energy_grid = self._fast_energy_grid if (self.config.use_max_coll_lut and self._fast_energy_grid is not None) else self._max_coll_energy_grid
        if self._jit_enabled:
            dt, trial = jit_kernels.determine_timestep_jit(
                rand, max_velocity, max_acceleration,
                freq, period,
                period_sq, energy_grid,
                self.config.max_cross_section_energy
            )
            self.trial_coll_freq = trial
            return dt

        # fallback numpy path
        s = - np.log(rand)
        max_energy = utils.energy_from_velocity(max_velocity)
        de = 2 * max_velocity * max_acceleration * s * period \
            + max_acceleration ** 2 * s ** 2 * period_sq
        e_end = np.clip(max_energy + de,
                        a_min=None,
                        a_max=self.config.max_cross_section_energy)
        freq_end = np.interp(e_end, energy_grid, freq)
        freq_needed = float(np.max(freq_end))
        idx = int(np.searchsorted(freq, freq_needed, side="left"))
        idx = int(np.clip(idx, 0, freq.size - 1))
        self.trial_coll_freq = 1.01 * freq[idx]
        return s / float(self.trial_coll_freq)

    def determine_collisions(self, gas_mixture: GasMixture,
                             velocity_norm: np.ndarray,
                             energy: np.ndarray) -> None:
        """
        Calculates the collision frequencies for all electrons with all cross sections,
        and chooses a collision type via a random number.

        Args:
            gas_mixture (GasMixture): cross section data
            velocity_norm (ndarray): norm of the velocity of each electron
            energy (ndarray): energy of each electron
        """

        energies = np.clip(energy, 0, self._max_coll_energy_grid[-1])
        idx = np.searchsorted(self._max_coll_energy_grid, energies, side="right")
        idx = np.clip(idx, 1, self._max_coll_energy_grid.size - 1)
        e_low = self._max_coll_energy_grid[idx - 1]
        e_high = self._max_coll_energy_grid[idx]
        denom = np.where(e_high == e_low, 1, e_high - e_low)
        weights = (energies - e_low) / denom
        cs_low = np.take_along_axis(gas_mixture.cross_section_table, (idx - 1)[None, :],
                                    axis=1)
        cs_high = np.take_along_axis(gas_mixture.cross_section_table, idx[None, :],
                                     axis=1)
        cross_sections = cs_low + (cs_high - cs_low) * weights
        collision_matrix = (self.config.gas_number_density * cross_sections
                            * velocity_norm[np.newaxis, :] / self.trial_coll_freq)

        # cumulative probabilities per electron; avoids large np.repeat
        cumprob = np.cumsum(collision_matrix, axis=0, dtype=float)
        rand = self._rand_array(energy.size)
        self.collision_by_electron = np.sum(rand > cumprob, axis=0)

    def perform_collisions(
            self, gas_mixture: GasMixture, position: np.ndarray, velocity: np.ndarray,
            energy: np.ndarray) -> Tuple[np.ndarray, np.ndarray, int, int, int]:
        """
        Calculates the electrons positions (created/removed) and velocities (scattered)
        after the collisions listed in the collision_by_electron array.

        Args:
            gas_mixture (GasMixture): cross section data
            position (ndarray): coordinates (x,y,z) of each electron (m)
            velocity (ndarray): velocity of each electron in (x,y,z) directions (m.s-1)
            energy (ndarray): energy of each electron (eV)

        Returns: the new position and velocity of electrons, the number of collisions,
            cations and anions produced
        """

        null_collisions = \
            self.collision_by_electron >= gas_mixture.number_of_cross_sections
        null_coll_p = position[:, null_collisions]
        null_coll_v = velocity[:, null_collisions]

        collided_mask = ~null_collisions
        collisions = self.collision_by_electron[collided_mask]
        n_coll = collisions.size
        pos = position[:, collided_mask]
        vel = velocity[:, collided_mask]
        ener = energy[collided_mask]

        use_parallel = self._jit_enabled and ener.size >= self._jit_threshold
        if use_parallel and self.config.use_jit_collision_classify:
            attachment_mask, ionization_mask_full = jit_kernels.classify_collisions(
                collisions, gas_mixture.is_attachment, gas_mixture.is_ionization
            )
        else:
            attachment_mask = gas_mixture.is_attachment[collisions]
            ionization_mask_full = gas_mixture.is_ionization[collisions]
        n_att = int(np.count_nonzero(attachment_mask))
        survive_mask = ~attachment_mask
        collisions = collisions[survive_mask]
        pos = pos[:, survive_mask]
        vel = vel[:, survive_mask]
        ener = ener[survive_mask]
        ionization_mask_full = ionization_mask_full[survive_mask]

        v_scattered_hat, cos_chi = self.unit_scattered_velocity(
            ener, vel, self.config.isotropic_scattering)

        mass_ratios = gas_mixture.mass_ratios[collisions]
        ionization_collisions = ionization_mask_full
        thresholds = gas_mixture.thresholds[collisions]
        sharing = np.ones(thresholds.shape)
        ion_share = None
        if np.any(ionization_collisions):
            ion_share = self._ionization_share_factor(ener[ionization_collisions])
            sharing[ionization_collisions] = ion_share
        if use_parallel and self.config.use_jit_energy_loss:
            energy_after_coll = jit_kernels.energy_after_collision_jit(
                ener, thresholds, mass_ratios, cos_chi
            )
        else:
            losses = thresholds + ener * mass_ratios * (1 - cos_chi)
            energy_after_coll = np.maximum(ener - losses, 0)
        v_scattered = v_scattered_hat * self._velocity_from_energy(
            energy_after_coll * sharing)

        new_e_pos = pos[:, ionization_collisions]
        if ion_share is None:
            ion_share = np.array([], dtype=float)
        new_e_energy = energy_after_coll[ionization_collisions] * (1 - ion_share)
        if new_e_energy.size:
            new_dirs = self.unit_scattered_velocity(
                ener[ionization_collisions],
                vel[:, ionization_collisions],
                True
            )[0]
            new_e_velocity = new_dirs * self._velocity_from_energy(new_e_energy)
        else:
            new_e_velocity = np.empty((3, 0), dtype=velocity.dtype)

        new_e_pos, new_e_velocity, new_e_energy = self._apply_ionization_limits(
            position.shape[1], new_e_pos, new_e_velocity, new_e_energy)
        n_ion = new_e_energy.size

        total_count = null_coll_p.shape[1] + pos.shape[1] + new_e_pos.shape[1]
        pos_end = np.empty((3, total_count), dtype=position.dtype)
        vel_end = np.empty_like(pos_end)
        if self._jit_enabled:
            jit_kernels.assemble_particles(
                null_coll_p, null_coll_v,
                pos, v_scattered,
                new_e_pos, new_e_velocity,
                pos_end, vel_end
            )
        else:
            cursor = 0
            if null_coll_p.size:
                span = null_coll_p.shape[1]
                pos_end[:, cursor:cursor + span] = null_coll_p
                vel_end[:, cursor:cursor + span] = null_coll_v
                cursor += span
            if pos.size:
                span = pos.shape[1]
                pos_end[:, cursor:cursor + span] = pos
                vel_end[:, cursor:cursor + span] = v_scattered
                cursor += span
            if new_e_pos.size:
                span = new_e_pos.shape[1]
                pos_end[:, cursor:cursor + span] = new_e_pos
                vel_end[:, cursor:cursor + span] = new_e_velocity

        if self.config.conserve:
            ne_ini = position.shape[1]
            ne_end = pos_end.shape[1]
            if ne_end == 0:
                pos_end = position[:, :1]
                vel_end = velocity[:, :1]
            elif ne_end > ne_ini:
                keep = self.rng.choice(ne_end, size=ne_ini, replace=False)
                pos_end = pos_end[:, keep]
                vel_end = vel_end[:, keep]
            elif ne_end < ne_ini:
                need = ne_ini - ne_end
                dup_idx = self.rng.choice(ne_end, size=need, replace=True)
                pos_end = np.hstack([pos_end, pos_end[:, dup_idx]])
                vel_end = np.hstack([vel_end, vel_end[:, dup_idx]])

        return pos_end, vel_end, n_coll, n_ion, n_att

    def _ionization_share_factor(self, energies: np.ndarray) -> np.ndarray:
        """
        Blend the ionization energy sharing factor toward a low-energy floor.
        """
        base = float(self.config.energy_sharing_factor)
        floor = getattr(self.config, "energy_sharing_floor", None)
        if floor is None:
            return np.full_like(energies, base, dtype=float)
        taper = float(self.config.energy_sharing_taper_energy)
        blend = np.clip(energies / taper, 0.0, 1.0)
        share = floor + (base - floor) * blend
        lo = min(base, floor)
        hi = max(base, floor)
        return np.clip(share, lo, hi)

    def _apply_ionization_limits(self, n_before: int,
                                 new_pos: np.ndarray,
                                 new_vel: np.ndarray,
                                 new_energy: np.ndarray
                                 ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Apply per-step growth caps and throttling near num_e_max to newly created
        electrons. Returns possibly down-sampled new electron arrays.
        """
        target = new_energy.size
        # hard cap relative to current population
        if self.config.max_ionization_growth_factor is not None:
            cap = int(np.floor(n_before * self.config.max_ionization_growth_factor))
            target = min(target, cap)

        # throttle as we approach num_e_max
        throttle = 1.0
        if self.config.num_e_throttle_ratio is not None and self.config.num_e_max > 0:
            start = self.config.num_e_throttle_ratio * self.config.num_e_max
            if n_before >= start:
                denom = max(self.config.num_e_max - start, 1e-9)
                throttle = max(0.0, 1.0 - (n_before - start) / denom)
                throttle = max(throttle, self.config.num_e_throttle_min_factor)
        if throttle < 1.0:
            target = min(target, int(np.floor(new_energy.size * throttle)))

        target = max(0, target)
        if target < new_energy.size:
            if target > 0:
                idx = self.rng.choice(new_energy.size, size=target, replace=False)
                new_pos = new_pos[:, idx]
                new_vel = new_vel[:, idx]
                new_energy = new_energy[idx]
            else:
                new_pos = new_pos[:, :0]
                new_vel = new_vel[:, :0]
                new_energy = new_energy[:0]
        return new_pos, new_vel, new_energy

    def unit_scattered_velocity(self, energy: np.ndarray, velocity: np.ndarray,
                                iso: bool) -> Tuple[np.ndarray, np.ndarray]:
        """
        Calculates the new direction of the velocity vector after scattering.

        Args:
            energy (ndarray): energy of each electron (eV)
            velocity (ndarray): velocity of each electron in (x,y,z) directions before
                the collision (m.s-1)
            iso (bool): isotropic scattering (True) or anisotropic scattering (False)

        Returns: normed velocities after collisions, cosine of polar scattering angle
        """

        if self._jit_enabled and energy.size >= self._jit_threshold:
            rand_phi = self._rand_array(energy.size)
            rand_chi = self._rand_array(energy.size)
            e_low = float(self.config.scattering_transition_low)
            e_high = float(self.config.scattering_transition_high)
            return jit_kernels.unit_scattered_velocity_jit(
                energy, velocity, iso, rand_phi, rand_chi, e_low, e_high
            )

        cos_chi, sin_chi, cos_phi, sin_phi = self.scattering_angles(energy, iso)

        v_hat = velocity / linalg.norm(velocity, axis=0, keepdims=True)

        e_x = np.zeros_like(velocity)
        e_x[0, :] = 1

        theta = np.arccos(v_hat[0, :])
        sin_theta = np.sin(theta)
        near_axis = np.abs(sin_theta) < 1e-12
        if np.any(near_axis):
            alt = np.zeros_like(velocity[:, near_axis])
            alt[1, :] = 1
            e_x[:, near_axis] = alt
            sin_theta[near_axis] = 1.0

        v_new_dir = \
            v_hat * cos_chi \
            + np.cross(v_hat, e_x, axis=0) * sin_chi * sin_phi / sin_theta \
            + np.cross(v_hat, np.cross(e_x, v_hat, axis=0), axis=0) \
            * sin_chi * cos_phi / sin_theta

        return v_new_dir, cos_chi

    def scattering_angles(self, energy: np.ndarray, iso: bool
                          ) -> (np.ndarray, np.ndarray, np.ndarray, np.ndarray):
        """
        Generates values for the polar (chi) and azimuthal (phi)
        isotropic or anisotropic scattering angles according to Vahedi (1995).

        Args:
            energy: array of electron energies
            iso: isotropic scattering or not

        Returns: 4 arrays cos(chi), sin(chi), cos(phi), sin(phi)
        """

        # choose scattering angle phi
        phi = 2 * np.pi * self.rng.random(energy.size)
        sin_phi = np.sin(phi)
        cos_phi = np.cos(phi)

        # choose scattering angle chi
        if iso:
            if self._cos_chi_lut is not None:
                cos_chi = jit_kernels.sample_cos_chi_from_lut(
                    self._cos_chi_lut, self._rand_array(energy.size)
                )
            else:
                cos_chi = 1 - 2 * self._rand_array(energy.size)
        else:
            if self.config.use_aniso_scatter_lut and self._aniso_lut is not None:
                cos_chi = jit_kernels.sample_aniso_cos_chi_from_lut(
                    self._aniso_lut, self._aniso_e_grid,
                    energy, self._rand_array(energy.size)
                )
            else:
                e_safe = np.maximum(energy, 1e-9)
                e_low, e_high = 1.0, 5.0
                anisotropic = (2 + e_safe
                               - 2 * (1 + e_safe) ** self._rand_array(energy.size)) \
                    / e_safe
                isotropic = 1 - 2 * self._rand_array(energy.size)
                blend = np.clip((energy - e_low) / (e_high - e_low), 0, 1)
                cos_chi = (1 - blend) * isotropic + blend * anisotropic
        cos_chi = np.clip(cos_chi, -1, 1)

        sin_chi = np.sqrt(1 - cos_chi ** 2)

        return cos_chi, sin_chi, cos_phi, sin_phi

    @staticmethod
    def _jit_velocity(energy: np.ndarray) -> np.ndarray:
        """
        JIT-backed velocity_from_energy helper.
        """
        return jit_kernels.velocity_from_energy_jit(np.asarray(energy, dtype=float))

    def _lut_velocity(self, energy: np.ndarray) -> np.ndarray:
        """
        LUT-backed velocity_from_energy helper using precomputed values on the
        max-collision energy grid.
        """
        if self._vel_lut_grid is None or self._vel_lut_values is None:
            return utils.velocity_from_energy(energy)
        return np.interp(energy, self._vel_lut_grid, self._vel_lut_values)

    def _rand_scalar(self) -> float:
        if self._jit_enabled and self._rng_state is not None:
            return float(jit_kernels.rng_uniform_scalar(self._rng_state))
        return float(self.rng.random())

    def _rand_array(self, n: int) -> np.ndarray:
        if self._jit_enabled and self._rng_state is not None:
            return jit_kernels.rng_uniform(self._rng_state, int(n))
        return self.rng.random(n)
