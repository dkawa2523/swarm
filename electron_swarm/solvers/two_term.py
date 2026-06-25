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

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.constants import (
    AMU_KG,
    BOLTZMANN_J_K,
    ELECTRON_MASS_KG,
    EV_TO_J,
    TOWNSEND,
)
from electron_swarm.core.cross_sections import (
    CrossSectionSet,
    gas_mass_amu,
)
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.core.solver_configs import TwoTermInternalConfig
from electron_swarm.core.transport import ElectronTransport
from electron_swarm.solvers.kinetic import (
    EffectiveCollisionData,
    KineticGrid,
    KineticOperatorBlock,
    assemble_collision_operator,
    assemble_energy_flux_operator,
    assemble_native_operator_blocks,
    build_effective_collision_data,
    cell_edges_from_centers,
    compute_rates_from_eedf,
    electron_speed_m_s,
    gas_number_density,
    make_two_term_energy_grid,
    mean_energy_from_eedf,
    normalize_eedf,
    transport_from_eedf,
    transport_from_reduced,
    weighted_integral,
)
from .base import SwarmSolver


_EPS = 1.0e-300


def _electron_speed(energy_eV: np.ndarray) -> np.ndarray:
    """Electron speed [m/s] from kinetic energy [eV]."""

    return electron_speed_m_s(energy_eV)


def _gas_number_density(config: SwarmConfig) -> float:
    return gas_number_density(config)


def _cell_edges_from_centers(centers: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return cell_edges_from_centers(centers)


def _make_energy_grid(
    cfg: TwoTermInternalConfig,
    *,
    max_eV_override: float | None = None,
    n_override: int | None = None,
    cross_sections: CrossSectionSet | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    grid = make_two_term_energy_grid(
        cfg,
        max_eV_override=max_eV_override,
        n_override=n_override,
        cross_sections=cross_sections,
    )
    return grid.energy_eV, grid.edges_eV, grid.widths_eV


def _grid_metadata(
    cfg: TwoTermInternalConfig, energy: np.ndarray
) -> dict[str, bool | float | int | str]:
    grid = cfg.energy_grid
    return {
        "grid_n_cells": int(len(energy)),
        "grid_min_eV": float(np.min(energy)) if len(energy) else float("nan"),
        "grid_max_eV": float(np.max(energy)) if len(energy) else float("nan"),
        "grid_spacing": grid.spacing,
        "threshold_refined": bool(grid.refine.enabled and len(energy) > grid.n),
    }


def _maxwell_eedf(energy: np.ndarray, kT_eV: float) -> np.ndarray:
    """Normalized Maxwell-Boltzmann EEDF F(eps), integral F d eps = 1."""

    kT = max(float(kT_eV), 1.0e-8)
    f = np.sqrt(np.maximum(energy, 0.0)) * np.exp(-np.maximum(energy, 0.0) / kT)
    return np.clip(f, 0.0, None)


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
class NativeDistributionResult:
    """Reusable native two-term distribution block for future direct-PN work."""

    energy_eV: np.ndarray
    edges_eV: np.ndarray
    widths_eV: np.ndarray
    eedf_eV_inv: np.ndarray
    diagnostics: NativeSolveDiagnostics
    metadata: dict[str, object]


NativeOperatorBlock = KineticOperatorBlock


class TwoTermSolver(SwarmSolver):
    """Unified electron Boltzmann two-term solver."""

    name = "two_term"

    def __init__(
        self,
        config: SwarmConfig,
        cross_sections: CrossSectionSet,
        solver_config: TwoTermInternalConfig,
    ) -> None:
        super().__init__(config, cross_sections)
        self.solver_config = solver_config
        backend = solver_config.backend
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

        cfg = self.solver_config

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
        transport: ElectronTransport | None = None
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
            "grid_n_cells": int(len(energy)),
            "grid_min_eV": float(np.min(energy)),
            "grid_spacing": self.solver_config.energy_grid.spacing,
            "threshold_refined": False,
        }
        return self._postprocess(
            e_over_n_Td,
            case_id,
            energy,
            widths,
            eedf,
            metadata=metadata,
            transport_override=transport,
        )

    # ------------------------------------------------------------------
    # Native BOLSIG-like backend.
    # ------------------------------------------------------------------
    def _solve_case_native(self, e_over_n_Td: float, case_id: str) -> SwarmCaseResult:
        native = self.solve_native_distribution(e_over_n_Td)
        return self._postprocess(
            e_over_n_Td,
            case_id,
            native.energy_eV,
            native.widths_eV,
            native.eedf_eV_inv,
            metadata=native.metadata,
        )

    def solve_native_reference_case(
        self, e_over_n_Td: float, case_id: str
    ) -> SwarmCaseResult:
        """Run the native Scharfetter-Gummel path regardless of backend setting."""

        return self._solve_case_native(e_over_n_Td, case_id)

    def solve_native_distribution(self, e_over_n_Td: float) -> NativeDistributionResult:
        """Solve only the native EEDF block and return reusable arrays/metadata."""

        cfg = self.solver_config
        max_eV = float(cfg.energy_grid.max_eV)
        previous: tuple[np.ndarray, np.ndarray] | None = None
        last_diag: NativeSolveDiagnostics | None = None
        cycles = max(1, cfg.adaptive_grid.max_cycles if cfg.adaptive_grid.enabled else 1)
        for cycle in range(cycles):
            energy, edges, widths = _make_energy_grid(
                cfg,
                max_eV_override=max_eV,
                cross_sections=self.cross_sections,
            )
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
            "cross_section_high_energy_extrapolation": (
                self.config.cross_sections.high_energy_extrapolation
            ),
        }
        metadata.update(_grid_metadata(cfg, energy))
        return NativeDistributionResult(
            energy_eV=energy,
            edges_eV=edges,
            widths_eV=widths,
            eedf_eV_inv=eedf,
            diagnostics=last_diag,
            metadata=metadata,
        )

    def make_native_energy_grid(
        self,
        *,
        max_eV_override: float | None = None,
        n_override: int | None = None,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return the native two-term energy grid used by the SG operator."""

        return _make_energy_grid(
            self.solver_config,
            max_eV_override=max_eV_override,
            n_override=n_override,
            cross_sections=self.cross_sections,
        )

    def assemble_native_operator_block(
        self,
        e_over_n_Td: float,
        energy: np.ndarray,
        edges: np.ndarray,
        widths: np.ndarray,
    ) -> NativeOperatorBlock:
        """Assemble the reusable native two-term energy-space operator.

        Future direct-PN work may reuse this validated lmax=1
        Scharfetter-Gummel block, but schema v2 does not route multi_term
        product runs through it.
        """

        energy = np.asarray(energy, dtype=float)
        edges = np.asarray(edges, dtype=float)
        widths = np.asarray(widths, dtype=float)
        if energy.ndim != 1 or edges.ndim != 1 or widths.ndim != 1:
            raise ValueError("Energy grid arrays must be one-dimensional")
        if len(edges) != len(energy) + 1 or len(widths) != len(energy):
            raise ValueError("Energy grid edges/widths do not match centers")
        if np.any(np.diff(energy) <= 0.0) or np.any(widths <= 0.0):
            raise ValueError("Energy grid must be strictly increasing")

        return assemble_native_operator_blocks(
            self.config,
            self.cross_sections,
            e_over_n_Td,
            KineticGrid(energy, edges, widths, electron_speed_m_s(energy)),
            self.solver_config,
        )

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
        cfg = self.solver_config
        block = self.assemble_native_operator_block(e_over_n_Td, energy, edges, widths)
        op = block.matrix

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
            new_growth = weighted_integral(lp, widths)
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
        return build_effective_collision_data(
            self.config,
            self.cross_sections,
            energy,
            N,
            self.solver_config,
        )

    def _field_diffusion(self, energy: np.ndarray, E: float, nu_m: np.ndarray) -> np.ndarray:
        """Electric-field heating energy diffusion coefficient [eV^2/s]."""

        from electron_swarm.solvers.kinetic import field_diffusion_eV2_s

        return field_diffusion_eV2_s(energy, E, nu_m)

    def _assemble_energy_flux_operator(
        self,
        energy: np.ndarray,
        edges: np.ndarray,
        widths: np.ndarray,
        E: float,
        coll: EffectiveCollisionData,
    ) -> sparse.csr_matrix:
        return assemble_energy_flux_operator(energy, widths, E, coll)

    def _assemble_collision_operator(
        self,
        energy: np.ndarray,
        widths: np.ndarray,
    ) -> sparse.csr_matrix:
        return assemble_collision_operator(
            self.config,
            self.cross_sections,
            energy,
            widths,
            _gas_number_density(self.config),
            self.solver_config,
        )

    def _assemble_operator(
        self,
        energy: np.ndarray,
        edges: np.ndarray,
        widths: np.ndarray,
        E: float,
        coll: EffectiveCollisionData,
    ) -> sparse.csr_matrix:
        return (
            self._assemble_energy_flux_operator(energy, edges, widths, E, coll)
            + self._assemble_collision_operator(energy, widths)
        ).tocsr()

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
        return normalize_eedf(eedf, widths)

    def _tail_probability(self, eedf: np.ndarray, widths: np.ndarray) -> float:
        cfg = self.solver_config.adaptive_grid
        n_tail = max(1, int(len(eedf) * cfg.tail_cells_fraction))
        return float(np.sum(np.clip(eedf[-n_tail:], 0.0, None) * widths[-n_tail:]))

    def _mean_energy(self, energy: np.ndarray, eedf: np.ndarray, widths: np.ndarray | None = None) -> float:
        if widths is None:
            widths = _cell_edges_from_centers(energy)[1]
        return mean_energy_from_eedf(energy, widths, eedf)

    def _transport_from_reduced(
        self, e_over_n_Td: float, N: float, muN: float, diffN: float
    ) -> ElectronTransport:
        return transport_from_reduced(e_over_n_Td, N, muN, diffN)

    def _transport_from_eedf(
        self,
        energy: np.ndarray,
        widths: np.ndarray,
        eedf: np.ndarray,
        N: float,
        e_over_n_Td: float,
    ) -> ElectronTransport:
        return transport_from_eedf(
            self.config,
            self.cross_sections,
            energy,
            widths,
            eedf,
            N,
            e_over_n_Td,
            self.solver_config,
        )

    def transport_from_eedf(
        self,
        energy: np.ndarray,
        widths: np.ndarray,
        eedf: np.ndarray,
        gas_number_density_m3: float,
        e_over_n_Td: float,
    ) -> ElectronTransport:
        """Compute native flux transport from a normalized EEDF.

        This is the product-facing helper used by direct-PN reduction tests so
        multi-term code does not need to reach into private two-term methods.
        """

        return self._transport_from_eedf(
            energy,
            widths,
            eedf,
            gas_number_density_m3,
            e_over_n_Td,
        )

    def _postprocess(
        self,
        e_over_n_Td: float,
        case_id: str,
        energy: np.ndarray,
        widths: np.ndarray,
        eedf: np.ndarray,
        *,
        metadata: dict[str, object],
        transport_override: ElectronTransport | None = None,
    ) -> SwarmCaseResult:
        N = _gas_number_density(self.config)
        transport = transport_override or self._transport_from_eedf(
            energy, widths, eedf, N, e_over_n_Td
        )
        if transport.reduced_electron_energy_mobility_eV_m2_V_s_m3 is None:
            energy_transport = self._transport_from_eedf(
                energy, widths, eedf, N, e_over_n_Td
            )
            transport = ElectronTransport(
                definition=transport.definition,
                gas_number_density_m3=transport.gas_number_density_m3,
                drift_velocity_m_s=transport.drift_velocity_m_s,
                reduced_mobility_m2_V_s_m3=transport.reduced_mobility_m2_V_s_m3,
                reduced_diffusion_L_m2_s_m3=transport.reduced_diffusion_L_m2_s_m3,
                reduced_diffusion_T_m2_s_m3=transport.reduced_diffusion_T_m2_s_m3,
                reduced_electron_energy_mobility_eV_m2_V_s_m3=(
                    energy_transport.reduced_electron_energy_mobility_eV_m2_V_s_m3
                ),
                reduced_electron_energy_diffusion_eV_m2_s_m3=(
                    energy_transport.reduced_electron_energy_diffusion_eV_m2_s_m3
                ),
            )
        mean_energy = self._mean_energy(energy, eedf, widths)
        rate_data = compute_rates_from_eedf(
            self.config,
            self.cross_sections,
            energy,
            widths,
            eedf,
            case_id=case_id,
            e_over_n_Td=e_over_n_Td,
            solver_name=self.name,
        )
        rates = rate_data.rates
        ion_rate = rate_data.ionization_rate_m3_s
        attach_rate = rate_data.attachment_rate_m3_s
        net_ion_freq = rate_data.net_ionization_frequency_s
        effective_townsend = net_ion_freq / max(abs(transport.drift_velocity_m_s) * N, 1.0e-300)
        metadata = dict(metadata)
        metadata["characteristic_energy_eV"] = transport.characteristic_energy_L_eV
        metadata["normalization_integral"] = weighted_integral(eedf, widths)
        metadata["convolution_ionization_rate_coefficient_m3_s"] = ion_rate
        metadata["convolution_attachment_rate_coefficient_m3_s"] = attach_rate
        metadata["convolution_effective_rate_coefficient_m3_s"] = ion_rate - attach_rate
        metadata["convolution_net_ionization_frequency_s-1"] = net_ion_freq
        metadata["net_ionization_frequency_model"] = "convolution_effective_rate"
        metadata["effective_townsend_1_m"] = net_ion_freq / max(
            abs(transport.drift_velocity_m_s), 1.0e-300
        )
        diagnostics = {"two_term": metadata}
        return SwarmCaseResult(
            solver=self.name,
            case_id=case_id,
            e_over_n_Td=e_over_n_Td,
            mean_energy_eV=mean_energy,
            net_ionization_frequency_s=net_ion_freq,
            effective_townsend_m2=effective_townsend,
            transport=transport,
            energy_eV=energy,
            eedf=eedf,
            energy_widths_eV=widths,
            rates=rates,
            metadata={},
            diagnostics=diagnostics,
        )
