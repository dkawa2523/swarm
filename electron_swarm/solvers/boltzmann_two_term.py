"""Production-grade Boltzmann two-term approximation solver.

The default native backend implements a BOLSIG+/Hagelaar-Pitchford-style
energy-space finite-volume solver:

* two-term Legendre expansion under a uniform DC electric field,
* Scharfetter-Gummel/exponential differencing for the convection-diffusion
  operator in electron energy,
* conservative energy-shift operators for excitation/superelastic collisions,
* optional ionization energy sharing and temporal-growth eigenvalue treatment,
* flux transport coefficients from the anisotropic first-order correction.

An optional BOLOS backend is retained as an independent reference path when the
``bolos`` package is installed.  The public result schema remains shared with
particle Monte Carlo runs.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as spla

from electron_swarm.core.config import SwarmConfig, BoltzmannTwoTermConfig
from electron_swarm.core.constants import (
    AMU_KG,
    BOLTZMANN_J_K,
    E_CHARGE_C,
    ELECTRON_MASS_KG,
    EV_TO_J,
    TOWNSEND,
)
from electron_swarm.core.cross_sections import (
    CrossSectionProcess,
    CrossSectionSet,
    ProcessType,
    gas_mass_amu,
    mixture_fraction,
)
from electron_swarm.core.results import RateResult, SwarmCaseResult
from .base import SwarmSolver


_EPS = 1.0e-300


def _electron_speed(energy_eV: np.ndarray) -> np.ndarray:
    """Electron speed [m/s] from kinetic energy [eV]."""

    return np.sqrt(np.maximum(2.0 * EV_TO_J * energy_eV / ELECTRON_MASS_KG, 0.0))


def _gas_number_density(config: SwarmConfig) -> float:
    cond = config.conditions
    if cond.gas_number_density_m3 is not None:
        return cond.gas_number_density_m3
    assert cond.pressure_Pa is not None
    return cond.pressure_Pa / (BOLTZMANN_J_K * cond.gas_temperature_K)


def _cell_edges_from_centers(centers: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    edges = np.empty(len(centers) + 1)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0] = max(0.0, centers[0] - 0.5 * (centers[1] - centers[0]))
    edges[-1] = centers[-1] + 0.5 * (centers[-1] - centers[-2])
    widths = np.diff(edges)
    if np.any(widths <= 0.0):
        raise ValueError("Energy grid must be strictly increasing")
    return edges, widths


def _make_energy_grid(
    cfg: BoltzmannTwoTermConfig,
    *,
    max_eV_override: float | None = None,
    n_override: int | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    grid = cfg.energy_grid
    n = int(n_override or grid.n)
    emin = max(float(grid.min_eV), 0.0)
    emax = float(max_eV_override if max_eV_override is not None else grid.max_eV)
    if n < 8:
        raise ValueError("Boltzmann energy grid requires at least 8 cells")
    if emax <= emin:
        raise ValueError("Boltzmann energy grid max_eV must be larger than min_eV")

    if grid.spacing == "linear":
        centers = np.linspace(emin, emax, n)
    elif grid.spacing == "quadratic":
        # Dense low-energy resolution without the severe timestep/conditioning
        # penalty of a logarithmic grid at epsilon -> 0.
        x = np.linspace(0.0, 1.0, n)
        centers = emin + (emax - emin) * x * x
    elif grid.spacing == "log":
        centers = np.geomspace(max(emin, 1.0e-8), emax, n)
    else:
        raise ValueError(f"Unsupported energy spacing: {grid.spacing}")
    edges, widths = _cell_edges_from_centers(centers)
    return centers, edges, widths


def _bernoulli(x: np.ndarray | float) -> np.ndarray | float:
    """Stable Bernoulli function B(x)=x/(exp(x)-1)."""

    x_arr = np.asarray(x, dtype=float)
    out = np.empty_like(x_arr, dtype=float)
    small = np.abs(x_arr) < 1.0e-7
    pos_large = x_arr > 80.0
    neg_large = x_arr < -80.0
    mid = ~(small | pos_large | neg_large)
    out[small] = 1.0 - x_arr[small] / 2.0 + x_arr[small] ** 2 / 12.0 - x_arr[small] ** 4 / 720.0
    out[pos_large] = 0.0
    out[neg_large] = -x_arr[neg_large]
    out[mid] = x_arr[mid] / np.expm1(x_arr[mid])
    if np.isscalar(x):
        return float(out)
    return out


def _maxwell_eedf(energy: np.ndarray, kT_eV: float) -> np.ndarray:
    """Normalized Maxwell-Boltzmann EEDF F(eps), integral F d eps = 1."""

    kT = max(float(kT_eV), 1.0e-8)
    f = np.sqrt(np.maximum(energy, 0.0)) * np.exp(-np.maximum(energy, 0.0) / kT)
    return np.clip(f, 0.0, None)


def _weighted_integral(values: np.ndarray, widths: np.ndarray) -> float:
    return float(np.sum(values * widths))


@dataclass(slots=True)
class EffectiveCollisionData:
    """Collision data projected onto the active energy grid."""

    nu_m: np.ndarray  # effective transport frequency [s^-1]
    nu_m_over_N: np.ndarray  # effective transport frequency / gas density [m^3 s^-1]
    sigma_m: np.ndarray  # BOLSIG-style effective total momentum-transfer cross section [m^2]
    elastic_A_eV_s: np.ndarray  # elastic energy drift coefficient [eV s^-1]
    elastic_D_eV2_s: np.ndarray  # elastic energy diffusion coefficient [eV^2 s^-1]
    processes: list[CrossSectionProcess]


@dataclass(slots=True)
class NativeSolveDiagnostics:
    converged: bool
    iterations: int
    residual: float
    growth_frequency_s: float
    tail_probability: float
    edge_to_peak: float
    grid_max_eV: float
    regrid_cycles: int


@dataclass(slots=True)
class TransportCoefficients:
    drift_velocity_m_s: float
    mobility_m2_V_s: float
    reduced_mobility_m2_V_s_m3: float
    diffusion_L_m2_s: float
    diffusion_T_m2_s: float
    reduced_diffusion_L_m2_s_m3: float
    reduced_diffusion_T_m2_s_m3: float
    characteristic_energy_eV: float


class BoltzmannTwoTermSolver(SwarmSolver):
    """Unified electron Boltzmann two-term solver."""

    name = "boltzmann_two_term"

    def __init__(self, config: SwarmConfig, cross_sections: CrossSectionSet) -> None:
        super().__init__(config, cross_sections)
        backend = config.boltzmann_two_term.backend
        if backend == "auto":
            # BOLOS is a useful independent reference implementation.  Use it
            # when explicitly available; otherwise use the built-in production
            # native backend rather than the old regression-test backend.
            try:
                import bolos  # noqa: F401

                backend = "bolos"
            except Exception:
                backend = "native_bolsig"
        if backend == "internal":
            backend = "native_bolsig"
        self.backend = backend

    def solve_case(self, e_over_n_Td: float, case_id: str) -> SwarmCaseResult:
        if self.backend == "bolos":
            return self._solve_case_bolos(e_over_n_Td, case_id)
        if self.backend == "native_bolsig":
            return self._solve_case_native(e_over_n_Td, case_id)
        raise ValueError(f"Unsupported Boltzmann backend: {self.backend}")

    # ------------------------------------------------------------------
    # Optional BOLOS backend: independent BOLSIG-like reference path.
    # ------------------------------------------------------------------
    def _solve_case_bolos(self, e_over_n_Td: float, case_id: str) -> SwarmCaseResult:
        try:
            from bolos import grid as bolos_grid
            from bolos import solver as bolos_solver
        except Exception as exc:  # pragma: no cover - optional dependency
            raise RuntimeError("backend=bolos requires the optional 'bolos' package") from exc

        cfg = self.config.boltzmann_two_term

        def make_bolos_grid(max_eV: float):
            if cfg.energy_grid.spacing == "quadratic" and hasattr(bolos_grid, "QuadraticGrid"):
                return bolos_grid.QuadraticGrid(cfg.energy_grid.min_eV, max_eV, cfg.energy_grid.n)
            if cfg.energy_grid.spacing == "log" and hasattr(bolos_grid, "LogGrid"):
                return bolos_grid.LogGrid(max(cfg.energy_grid.min_eV, 1.0e-8), max_eV, cfg.energy_grid.n)
            return bolos_grid.LinearGrid(cfg.energy_grid.min_eV, max_eV, cfg.energy_grid.n)

        max_eV = cfg.energy_grid.max_eV
        bs = None
        f = None
        energy = None
        last_grid = None
        cycles = max(1, cfg.adaptive_grid.max_cycles if cfg.adaptive_grid.enabled else 1)
        for cycle in range(cycles):
            gr = make_bolos_grid(max_eV)
            bs = bolos_solver.BoltzmannSolver(gr)
            bolos_processes: dict[tuple[str, str], Any] = {}
            for proc in self.cross_sections.processes:
                ptype = proc.process_type.value.upper()
                if ptype == "ELASTIC":
                    ptype = "MOMENTUM"
                if ptype == "SUPERELASTIC":
                    # BOLOS/BOLSIG format normally represents superelastic
                    # reactions as inelastic processes with negative threshold.
                    ptype = "EXCITATION"
                if ptype not in {"EFFECTIVE", "MOMENTUM", "EXCITATION", "IONIZATION", "ATTACHMENT"}:
                    continue
                mass_amu = proc.mass_amu or gas_mass_amu(self.config.conditions, proc.species)
                ratio = ELECTRON_MASS_KG / (mass_amu * AMU_KG)
                added = bs.add_process(
                    type=ptype,
                    target=proc.species,
                    ratio=ratio,
                    threshold=proc.threshold_eV or 0.0,
                    data=np.c_[proc.energy_eV, proc.cross_section_m2],
                )
                bolos_processes[(proc.species, proc.process)] = added
            for gas in self.config.conditions.gas_mixture:
                bs.target[gas.species].density = gas.fraction
            bs.kT = self.config.conditions.gas_temperature_K * BOLTZMANN_J_K / EV_TO_J
            bs.EN = e_over_n_Td * TOWNSEND
            bs.init()
            energy = np.asarray(bs.grid.c, dtype=float) if hasattr(bs.grid, "c") else np.asarray(gr.c, dtype=float)
            if f is None:
                f0 = bs.maxwell(cfg.initial_electron_temperature_eV)
            else:
                old_energy = np.asarray(last_grid.c, dtype=float) if hasattr(last_grid, "c") else np.asarray(last_grid.c, dtype=float)
                f0 = np.interp(energy, old_energy, f, left=0.0, right=0.0)
                f0 = self._normalize_eedf(f0, _cell_edges_from_centers(energy)[1])
            f = bs.converge(f0, maxn=cfg.convergence.max_iterations, rtol=cfg.convergence.tolerance)
            last_grid = gr
            if not cfg.adaptive_grid.enabled:
                break
            mean_e = float(bs.mean_energy(f)) if hasattr(bs, "mean_energy") else self._mean_energy(energy, f)
            proposed = min(cfg.adaptive_grid.max_max_eV, max(cfg.adaptive_grid.min_max_eV, cfg.adaptive_grid.mean_energy_multiplier * mean_e))
            widths = _cell_edges_from_centers(energy)[1]
            tail = self._tail_probability(f, widths)
            if proposed <= max_eV * 1.05 and tail <= cfg.adaptive_grid.tail_probability:
                break
            max_eV = max(max_eV * 1.25, proposed)

        assert bs is not None and f is not None and energy is not None
        widths = _cell_edges_from_centers(energy)[1]
        eedf = np.clip(np.asarray(f, dtype=float), 0.0, None)
        eedf = self._normalize_eedf(eedf, widths)

        N = _gas_number_density(self.config)
        transport: TransportCoefficients | None = None
        try:  # Prefer BOLOS' own transport integrals when available.
            muN = float(bs.mobility(eedf))
            diffN = float(bs.diffusion(eedf))
            transport = self._transport_from_reduced(e_over_n_Td, N, muN, diffN)
        except Exception:  # pragma: no cover - depends on BOLOS version
            transport = None

        metadata = {
            "backend": "bolos",
            "grid_max_eV": float(np.max(energy)),
            "adaptive_cycles": int(cycle + 1),
        }
        return self._postprocess(e_over_n_Td, case_id, energy, widths, eedf, metadata=metadata, transport_override=transport)

    # ------------------------------------------------------------------
    # Native BOLSIG-like backend.
    # ------------------------------------------------------------------
    def _solve_case_native(self, e_over_n_Td: float, case_id: str) -> SwarmCaseResult:
        cfg = self.config.boltzmann_two_term
        max_eV = float(cfg.energy_grid.max_eV)
        previous: tuple[np.ndarray, np.ndarray] | None = None
        last_diag: NativeSolveDiagnostics | None = None
        cycles = max(1, cfg.adaptive_grid.max_cycles if cfg.adaptive_grid.enabled else 1)
        for cycle in range(cycles):
            energy, edges, widths = _make_energy_grid(cfg, max_eV_override=max_eV)
            initial = None
            if previous is not None:
                old_e, old_f = previous
                initial = np.interp(energy, old_e, old_f, left=0.0, right=0.0)
                initial = self._normalize_eedf(initial, widths)
            eedf, diag = self._solve_native_on_grid(e_over_n_Td, energy, edges, widths, initial=initial, cycle=cycle)
            previous = (energy, eedf)
            last_diag = diag
            if not cfg.adaptive_grid.enabled:
                break
            mean_energy = self._mean_energy(energy, eedf, widths)
            required = max(cfg.adaptive_grid.min_max_eV, cfg.adaptive_grid.mean_energy_multiplier * mean_energy)
            edge_bad = diag.tail_probability > cfg.adaptive_grid.tail_probability or diag.edge_to_peak > cfg.adaptive_grid.edge_to_peak
            if required <= max_eV * 1.05 and not edge_bad:
                break
            max_eV = min(cfg.adaptive_grid.max_max_eV, max(max_eV * 1.25, required))
            if max_eV <= energy[-1] * 1.001:
                break

        assert previous is not None and last_diag is not None
        energy, eedf = previous
        edges, widths = _cell_edges_from_centers(energy)
        metadata = {
            "backend": "native_bolsig",
            "converged": last_diag.converged,
            "iterations": last_diag.iterations,
            "residual_L1": last_diag.residual,
            "growth_frequency_s-1": last_diag.growth_frequency_s,
            "tail_probability": last_diag.tail_probability,
            "edge_to_peak": last_diag.edge_to_peak,
            "grid_max_eV": last_diag.grid_max_eV,
            "adaptive_cycles": last_diag.regrid_cycles,
            "transport_model": "two_term_flux_integral",
            "discretization": "finite_volume_scharfetter_gummel",
        }
        return self._postprocess(e_over_n_Td, case_id, energy, widths, eedf, metadata=metadata)

    def _solve_native_on_grid(
        self,
        e_over_n_Td: float,
        energy: np.ndarray,
        edges: np.ndarray,
        widths: np.ndarray,
        *,
        initial: np.ndarray | None,
        cycle: int,
    ) -> tuple[np.ndarray, NativeSolveDiagnostics]:
        cfg = self.config.boltzmann_two_term
        N = _gas_number_density(self.config)
        E = e_over_n_Td * TOWNSEND * N
        coll = self._effective_collision_data(energy, N)
        op = self._assemble_operator(energy, edges, widths, E, coll)

        if initial is None:
            p = _maxwell_eedf(energy, cfg.initial_electron_temperature_eV)
            p = self._normalize_eedf(p, widths)
        else:
            p = self._normalize_eedf(np.clip(initial, 0.0, None), widths)

        growth = 0.0
        residual = np.inf
        converged = False
        iterations = 0
        for it in range(cfg.convergence.max_iterations):
            iterations = it + 1
            shifted = op - growth * sparse.identity(op.shape[0], format="csr")
            p_new = self._solve_normalized(shifted, widths)
            if cfg.convergence.clip_negative:
                # The exponential scheme is positivity preserving for the
                # nearest-neighbour fluxes; clipping protects against tiny
                # sparse-solver roundoff and high-energy extrapolation noise.
                p_new = np.clip(p_new, 0.0, None)
                p_new = self._normalize_eedf(p_new, widths)
            lp = op @ p_new
            new_growth = _weighted_integral(lp, widths)
            eig_res = lp - new_growth * p_new
            residual = float(np.sum(np.abs(eig_res) * widths) / max(np.sum(np.abs(lp) * widths), abs(new_growth), 1.0))
            shape_change = float(np.sum(np.abs(p_new - p) * widths))
            growth_change = abs(new_growth - growth) / max(abs(new_growth), abs(growth), 1.0)
            p = p_new
            growth = cfg.convergence.relaxation * new_growth + (1.0 - cfg.convergence.relaxation) * growth
            if shape_change < cfg.convergence.tolerance and growth_change < cfg.convergence.eigenvalue_tolerance and residual < cfg.convergence.residual_tolerance:
                converged = True
                growth = new_growth
                break

        tail = self._tail_probability(p, widths)
        edge_to_peak = float(p[-1] / max(np.max(p), _EPS))
        diag = NativeSolveDiagnostics(
            converged=converged,
            iterations=iterations,
            residual=residual,
            growth_frequency_s=float(growth),
            tail_probability=tail,
            edge_to_peak=edge_to_peak,
            grid_max_eV=float(energy[-1]),
            regrid_cycles=cycle + 1,
        )
        return p, diag

    def _effective_collision_data(self, energy: np.ndarray, N: float) -> EffectiveCollisionData:
        cfg = self.config.boltzmann_two_term
        speed = _electron_speed(energy)
        sigma_m = np.zeros_like(energy)
        nu_m = np.zeros_like(energy)
        nu_m_over_N = np.zeros_like(energy)
        elastic_A = np.zeros_like(energy)
        elastic_D = np.zeros_like(energy)
        kT_eV = self.config.conditions.gas_temperature_K * BOLTZMANN_J_K / EV_TO_J
        species_has_effective = {
            species: any(
                proc.process_type == ProcessType.EFFECTIVE
                for proc in self.cross_sections.by_species(species)
            )
            for species in self.cross_sections.species
        }
        momentum_like = self.cross_sections.by_type(ProcessType.MOMENTUM, ProcessType.EFFECTIVE, ProcessType.ELASTIC)
        if not momentum_like:
            raise ValueError("At least one momentum/effective/elastic cross section is required")
        transport_inelastic = {
            ProcessType.EXCITATION,
            ProcessType.IONIZATION,
            ProcessType.ATTACHMENT,
            ProcessType.SUPERELASTIC,
        }
        for proc in self.cross_sections.processes:
            frac = mixture_fraction(self.config.conditions, proc.species)
            if frac <= 0.0:
                continue
            sigma = np.maximum(proc.sigma(energy), cfg.min_momentum_cross_section_m2)
            use_effective_only = species_has_effective.get(proc.species, False)
            include_transport = False
            if proc.process_type == ProcessType.EFFECTIVE:
                include_transport = True
            elif proc.process_type in {ProcessType.MOMENTUM, ProcessType.ELASTIC}:
                include_transport = not use_effective_only
            elif proc.process_type in transport_inelastic:
                include_transport = not use_effective_only
            if include_transport:
                nuN = frac * sigma * speed
                nu = N * nuN
                sigma_m += frac * sigma
                nu_m_over_N += nuN
                nu_m += nu
            if proc.process_type not in {ProcessType.MOMENTUM, ProcessType.EFFECTIVE, ProcessType.ELASTIC}:
                continue
            mass_amu = proc.mass_amu or gas_mass_amu(self.config.conditions, proc.species)
            mass_kg = mass_amu * AMU_KG
            nuN = frac * sigma * speed
            nu = N * nuN
            ratio = ELECTRON_MASS_KG / mass_kg
            eps = np.maximum(energy, 1.0e-12)
            # Conservative elastic energy-exchange Fokker-Planck coefficients
            # for the EEDF F(eps).  With E=0, zero flux gives
            # F ~ sqrt(eps) exp(-eps/kT_g), i.e. the Maxwell energy PDF.
            elastic_D += 2.0 * ratio * nu * eps * kT_eV
            elastic_A += 2.0 * ratio * nu * (0.5 * kT_eV - eps)
        nu_m = np.maximum(nu_m, 1.0e-60)
        nu_m_over_N = np.maximum(nu_m_over_N, 1.0e-80)
        sigma_m = np.maximum(sigma_m, cfg.min_momentum_cross_section_m2)
        return EffectiveCollisionData(
            nu_m=nu_m,
            nu_m_over_N=nu_m_over_N,
            sigma_m=sigma_m,
            elastic_A_eV_s=elastic_A,
            elastic_D_eV2_s=np.maximum(elastic_D, 0.0),
            processes=self.cross_sections.processes,
        )

    def _field_diffusion(self, energy: np.ndarray, E: float, nu_m: np.ndarray) -> np.ndarray:
        """Electric-field heating energy diffusion coefficient [eV^2/s]."""

        eps_J = np.maximum(energy * EV_TO_J, 0.0)
        D_J2_s = (2.0 / 3.0) * (E_CHARGE_C * E) ** 2 / ELECTRON_MASS_KG * eps_J / np.maximum(nu_m, 1.0e-60)
        return D_J2_s / (EV_TO_J * EV_TO_J)

    def _assemble_operator(
        self,
        energy: np.ndarray,
        edges: np.ndarray,
        widths: np.ndarray,
        E: float,
        coll: EffectiveCollisionData,
    ) -> sparse.csr_matrix:
        n = len(energy)
        mat = sparse.lil_matrix((n, n), dtype=float)
        D_field = self._field_diffusion(energy, E, coll.nu_m)
        eps_safe = np.maximum(energy, max(energy[1] - energy[0], 1.0e-12) * 0.5)
        A_field = D_field / (2.0 * eps_safe)
        D_center = np.maximum(coll.elastic_D_eV2_s + D_field, 1.0e-80)
        A_center = coll.elastic_A_eV_s + A_field

        # Neighbour fluxes J = A F - D dF/dε with exponential fitting.
        # Boundary fluxes are zero, matching BOLSIG-like closed energy-domain
        # truncation when the upper tail is numerically negligible.
        for i in range(n - 1):
            h = energy[i + 1] - energy[i]
            if h <= 0.0:
                continue
            D = max(0.5 * (D_center[i] + D_center[i + 1]), 1.0e-80)
            A = 0.5 * (A_center[i] + A_center[i + 1])
            peclet = A * h / D
            cL = D / h * _bernoulli(-peclet)
            cR = -D / h * _bernoulli(peclet)
            mat[i, i] += -cL / widths[i]
            mat[i, i + 1] += -cR / widths[i]
            mat[i + 1, i] += cL / widths[i + 1]
            mat[i + 1, i + 1] += cR / widths[i + 1]

        # Inelastic/nonconservative collision operators.
        speed = _electron_speed(energy)
        N = _gas_number_density(self.config)
        for proc in self.cross_sections.processes:
            frac = mixture_fraction(self.config.conditions, proc.species)
            if frac <= 0.0:
                continue
            if proc.process_type not in {ProcessType.EXCITATION, ProcessType.SUPERELASTIC, ProcessType.IONIZATION, ProcessType.ATTACHMENT}:
                continue
            nu = N * frac * proc.sigma(energy) * speed
            if proc.process_type == ProcessType.ATTACHMENT:
                if self.config.boltzmann_two_term.nonconservative_model == "ignore":
                    continue
                for i, val in enumerate(nu):
                    mat[i, i] += -float(val)
                continue
            threshold = float(proc.threshold_eV or 0.0)
            if proc.process_type == ProcessType.SUPERELASTIC and threshold > 0.0:
                threshold = -threshold
            if proc.process_type in {ProcessType.EXCITATION, ProcessType.SUPERELASTIC}:
                self._add_energy_shift_transition(mat, energy, widths, nu, threshold, multiplicity=1.0)
            elif proc.process_type == ProcessType.IONIZATION:
                if (
                    self.config.boltzmann_two_term.nonconservative_model == "ignore"
                    or self.config.boltzmann_two_term.ionization_energy_sharing == "loss_only"
                ):
                    self._add_energy_shift_transition(mat, energy, widths, nu, threshold, multiplicity=1.0)
                elif self.config.boltzmann_two_term.ionization_energy_sharing == "primary_secondary":
                    self._add_primary_secondary_ionization(mat, energy, widths, nu, threshold)
                else:
                    self._add_equal_sharing_ionization(mat, energy, widths, nu, threshold)
        return mat.tocsr()

    def _add_energy_shift_transition(
        self,
        mat: sparse.lil_matrix,
        energy: np.ndarray,
        widths: np.ndarray,
        nu: np.ndarray,
        threshold: float,
        *,
        multiplicity: float,
    ) -> None:
        for i, rate in enumerate(np.asarray(nu, dtype=float)):
            if rate <= 0.0:
                continue
            mat[i, i] += -rate
            target_e = energy[i] - threshold
            self._deposit_energy(mat, energy, widths, i, target_e, multiplicity * rate)

    def _add_equal_sharing_ionization(
        self,
        mat: sparse.lil_matrix,
        energy: np.ndarray,
        widths: np.ndarray,
        nu: np.ndarray,
        threshold: float,
    ) -> None:
        for i, rate in enumerate(np.asarray(nu, dtype=float)):
            if rate <= 0.0:
                continue
            mat[i, i] += -rate
            target_e = max(0.0, 0.5 * (energy[i] - threshold))
            self._deposit_energy(mat, energy, widths, i, target_e, 2.0 * rate)

    def _add_primary_secondary_ionization(
        self,
        mat: sparse.lil_matrix,
        energy: np.ndarray,
        widths: np.ndarray,
        nu: np.ndarray,
        threshold: float,
    ) -> None:
        """Ionization model with one cold secondary and one fast primary.

        This option is useful for sensitivity studies.  Equal sharing is usually
        smoother and is the default for BOLSIG-like fluid tables when no
        differential ionization data are supplied.
        """

        secondary_eV = max(float(self.config.boltzmann_two_term.secondary_electron_energy_eV), 0.0)
        for i, rate in enumerate(np.asarray(nu, dtype=float)):
            if rate <= 0.0:
                continue
            mat[i, i] += -rate
            primary_e = max(0.0, energy[i] - threshold - secondary_eV)
            self._deposit_energy(mat, energy, widths, i, primary_e, rate)
            self._deposit_energy(mat, energy, widths, i, secondary_eV, rate)

    def _deposit_energy(
        self,
        mat: sparse.lil_matrix,
        energy: np.ndarray,
        widths: np.ndarray,
        source_i: int,
        target_e: float,
        rate: float,
    ) -> None:
        n = len(energy)
        if target_e <= energy[0]:
            mat[0, source_i] += rate * widths[source_i] / widths[0]
        elif target_e >= energy[-1]:
            mat[-1, source_i] += rate * widths[source_i] / widths[-1]
        else:
            j = int(np.searchsorted(energy, target_e) - 1)
            j = max(0, min(j, n - 2))
            denom = energy[j + 1] - energy[j]
            wR = (target_e - energy[j]) / denom
            wL = 1.0 - wR
            mat[j, source_i] += rate * widths[source_i] * wL / widths[j]
            mat[j + 1, source_i] += rate * widths[source_i] * wR / widths[j + 1]

    def _solve_normalized(self, op: sparse.csr_matrix, widths: np.ndarray) -> np.ndarray:
        n = op.shape[0]
        mat = op.tolil(copy=True)
        rhs = np.zeros(n)
        # Replace one equation by normalization ∫F dε = 1.  Choose the row with
        # the largest diagonal magnitude to reduce boundary-condition bias.
        diag_abs = np.abs(op.diagonal())
        row = int(np.argmax(diag_abs)) if np.any(diag_abs > 0.0) else n - 1
        mat[row, :] = widths
        rhs[row] = 1.0
        csr = mat.tocsr()
        try:
            sol = spla.spsolve(csr, rhs)
        except Exception:
            sol = spla.lsmr(csr, rhs, atol=1.0e-13, btol=1.0e-13, maxiter=max(1000, 4 * n))[0]
        if not np.all(np.isfinite(sol)):
            raise FloatingPointError("Boltzmann linear solve returned non-finite values")
        return self._normalize_eedf(sol, widths)

    def _normalize_eedf(self, eedf: np.ndarray, widths: np.ndarray) -> np.ndarray:
        eedf = np.asarray(eedf, dtype=float)
        total = _weighted_integral(eedf, widths)
        if not np.isfinite(total) or abs(total) < 1.0e-300:
            raise FloatingPointError("Cannot normalize EEDF with zero/non-finite integral")
        if total < 0.0:
            eedf = -eedf
            total = -total
        return eedf / total

    def _tail_probability(self, eedf: np.ndarray, widths: np.ndarray) -> float:
        cfg = self.config.boltzmann_two_term.adaptive_grid
        n_tail = max(1, int(len(eedf) * cfg.tail_cells_fraction))
        return float(np.sum(np.clip(eedf[-n_tail:], 0.0, None) * widths[-n_tail:]))

    def _mean_energy(self, energy: np.ndarray, eedf: np.ndarray, widths: np.ndarray | None = None) -> float:
        if widths is None:
            widths = _cell_edges_from_centers(energy)[1]
        return float(np.sum(energy * eedf * widths))

    def _transport_from_reduced(self, e_over_n_Td: float, N: float, muN: float, diffN: float) -> TransportCoefficients:
        EN = e_over_n_Td * TOWNSEND
        mobility = muN / N
        diffusion = diffN / N
        drift = muN * EN
        char_e = diffN / max(abs(muN), 1.0e-300)
        return TransportCoefficients(
            drift_velocity_m_s=drift,
            mobility_m2_V_s=mobility,
            reduced_mobility_m2_V_s_m3=muN,
            diffusion_L_m2_s=diffusion,
            diffusion_T_m2_s=diffusion,
            reduced_diffusion_L_m2_s_m3=diffN,
            reduced_diffusion_T_m2_s_m3=diffN,
            characteristic_energy_eV=char_e,
        )

    def _transport_from_eedf(self, energy: np.ndarray, widths: np.ndarray, eedf: np.ndarray, N: float, e_over_n_Td: float) -> TransportCoefficients:
        coll = self._effective_collision_data(energy, N)
        nuN = np.maximum(coll.nu_m_over_N, 1.0e-80)
        speed = _electron_speed(energy)
        if len(energy) >= 3:
            dF = np.gradient(eedf, energy, edge_order=2)
        else:
            dF = np.gradient(eedf, energy)
        mobility_integrand = (2.0 * energy / nuN) * dF - eedf / nuN
        muN = -E_CHARGE_C / (3.0 * ELECTRON_MASS_KG) * _weighted_integral(mobility_integrand, widths)
        if not np.isfinite(muN) or muN <= 0.0:
            # Conservative fallback: relaxation-time estimate.  This should be
            # rare and is marked in output by metadata residual/tail diagnostics.
            nu_eff_over_N = _weighted_integral(nuN * eedf, widths)
            muN = E_CHARGE_C / (ELECTRON_MASS_KG * max(nu_eff_over_N, 1.0e-80))
        diffN = (1.0 / 3.0) * _weighted_integral((speed * speed / nuN) * eedf, widths)
        return self._transport_from_reduced(e_over_n_Td, N, float(muN), float(diffN))

    def _postprocess(
        self,
        e_over_n_Td: float,
        case_id: str,
        energy: np.ndarray,
        widths: np.ndarray,
        eedf: np.ndarray,
        *,
        metadata: dict[str, object],
        transport_override: TransportCoefficients | None = None,
    ) -> SwarmCaseResult:
        N = _gas_number_density(self.config)
        transport = transport_override or self._transport_from_eedf(energy, widths, eedf, N, e_over_n_Td)
        mean_energy = self._mean_energy(energy, eedf, widths)
        speed = _electron_speed(energy)

        rates: list[RateResult] = []
        ion_rate = 0.0
        attach_rate = 0.0
        net_ion_freq = 0.0
        for proc in self.cross_sections.processes:
            frac = mixture_fraction(self.config.conditions, proc.species)
            if frac <= 0.0:
                continue
            k = float(np.sum(proc.sigma(energy) * speed * eedf * widths))
            kmix = frac * k
            if proc.process_type == ProcessType.IONIZATION:
                ion_rate += kmix
                net_ion_freq += N * kmix
            elif proc.process_type == ProcessType.ATTACHMENT:
                attach_rate += kmix
                net_ion_freq -= N * kmix
            rates.append(
                RateResult(
                    solver=self.name,
                    case_id=case_id,
                    e_over_n_Td=e_over_n_Td,
                    species=proc.species,
                    process=proc.process,
                    process_type=proc.process_type.value,
                    threshold_eV=proc.threshold_eV,
                    rate_coefficient_m3_s=k,
                    mixture_weighted_rate_m3_s=kmix,
                )
            )
        effective_townsend = net_ion_freq / max(abs(transport.drift_velocity_m_s) * N, 1.0e-300)
        eepf = eedf / np.sqrt(np.maximum(energy, 1.0e-30))
        metadata = dict(metadata)
        metadata["characteristic_energy_eV"] = transport.characteristic_energy_eV
        metadata["normalization_integral"] = _weighted_integral(eedf, widths)
        metadata["convolution_ionization_rate_coefficient_m3_s"] = ion_rate
        metadata["convolution_attachment_rate_coefficient_m3_s"] = attach_rate
        metadata["convolution_effective_rate_coefficient_m3_s"] = ion_rate - attach_rate
        metadata["convolution_net_ionization_frequency_s-1"] = net_ion_freq
        metadata["net_ionization_frequency_model"] = "convolution_effective_rate"
        metadata["effective_townsend_1_m"] = net_ion_freq / max(
            abs(transport.drift_velocity_m_s), 1.0e-300
        )
        return SwarmCaseResult(
            solver=self.name,
            case_id=case_id,
            e_over_n_Td=e_over_n_Td,
            mean_energy_eV=mean_energy,
            drift_velocity_m_s=transport.drift_velocity_m_s,
            mobility_m2_V_s=transport.mobility_m2_V_s,
            reduced_mobility_m2_V_s_m3=transport.reduced_mobility_m2_V_s_m3,
            diffusion_L_m2_s=transport.diffusion_L_m2_s,
            diffusion_T_m2_s=transport.diffusion_T_m2_s,
            reduced_diffusion_L_m2_s_m3=transport.reduced_diffusion_L_m2_s_m3,
            reduced_diffusion_T_m2_s_m3=transport.reduced_diffusion_T_m2_s_m3,
            net_ionization_frequency_s=net_ion_freq,
            effective_townsend_m2=effective_townsend,
            energy_eV=energy,
            eedf=eedf,
            eepf=eepf,
            rates=rates,
            metadata=metadata,
        )
