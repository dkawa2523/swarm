"""Shared kinetic projection utilities for Boltzmann-family solvers."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse

from electron_swarm.core.config import SwarmConfig, TwoTermInternalConfig
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
from electron_swarm.core.results import RateResult
from electron_swarm.grids.energy import build_energy_grid


@dataclass(frozen=True, slots=True)
class KineticGrid:
    energy_eV: np.ndarray
    edges_eV: np.ndarray
    widths_eV: np.ndarray
    speed_m_s: np.ndarray


@dataclass(slots=True)
class EffectiveCollisionData:
    nu_m: np.ndarray
    nu_m_over_N: np.ndarray
    sigma_m: np.ndarray
    sigma_total_like: np.ndarray
    inelastic_loss_frequency_s_inv: np.ndarray
    elastic_A_eV_s: np.ndarray
    elastic_D_eV2_s: np.ndarray
    processes: list[CrossSectionProcess]
    effective_momentum_model: str
    sigma_total_like_model: str
    suppressed_momentum_processes: tuple[str, ...]


@dataclass(slots=True)
class KineticOperatorBlock:
    grid: KineticGrid
    gas_number_density_m3: float
    electric_field_V_m: float
    collisions: EffectiveCollisionData
    energy_flux_matrix: sparse.csr_matrix
    collision_matrix: sparse.csr_matrix
    matrix: sparse.csr_matrix
    discretization: str = "finite_volume_scharfetter_gummel"

    @property
    def energy_eV(self) -> np.ndarray:
        return self.grid.energy_eV

    @property
    def edges_eV(self) -> np.ndarray:
        return self.grid.edges_eV

    @property
    def widths_eV(self) -> np.ndarray:
        return self.grid.widths_eV


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


@dataclass(frozen=True, slots=True)
class RateConvolution:
    rates: list[RateResult]
    ionization_rate_m3_s: float
    attachment_rate_m3_s: float
    net_ionization_frequency_s: float


def electron_speed_m_s(energy_eV: np.ndarray) -> np.ndarray:
    return np.sqrt(np.maximum(2.0 * EV_TO_J * energy_eV / ELECTRON_MASS_KG, 0.0))


def gas_number_density(config: SwarmConfig) -> float:
    cond = config.conditions
    if cond.gas_number_density_m3 is not None:
        return float(cond.gas_number_density_m3)
    assert cond.pressure_Pa is not None
    return float(cond.pressure_Pa / (BOLTZMANN_J_K * cond.gas_temperature_K))


def cell_edges_from_centers(centers: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    edges = np.empty(len(centers) + 1)
    edges[1:-1] = 0.5 * (centers[:-1] + centers[1:])
    edges[0] = max(0.0, centers[0] - 0.5 * (centers[1] - centers[0]))
    edges[-1] = centers[-1] + 0.5 * (centers[-1] - centers[-2])
    widths = np.diff(edges)
    if np.any(widths <= 0.0):
        raise ValueError("Energy grid must be strictly increasing")
    return edges, widths


def make_two_term_energy_grid(
    cfg: TwoTermInternalConfig,
    *,
    max_eV_override: float | None = None,
    n_override: int | None = None,
    cross_sections: CrossSectionSet | None = None,
) -> KineticGrid:
    grid = cfg.energy_grid
    n = int(n_override or grid.n)
    emin = max(float(grid.min_eV), 0.0)
    emax = float(max_eV_override if max_eV_override is not None else grid.max_eV)
    if n < 8:
        raise ValueError("Boltzmann energy grid requires at least 8 cells")
    if emax <= emin:
        raise ValueError("Boltzmann energy grid max_eV must be larger than min_eV")

    if grid.refine.enabled:
        refined = build_energy_grid(
            min_eV=emin,
            max_eV=emax,
            n=n,
            spacing=grid.spacing,
            cross_sections=cross_sections,
            refine=True,
            threshold_padding_eV=grid.refine.threshold_padding_eV,
            points_per_threshold=grid.refine.points_per_threshold,
            max_extra_points=grid.refine.max_extra_points,
        )
        energy = refined.centers_eV
        edges = refined.edges_eV
        widths = refined.widths_eV
    elif grid.spacing == "linear":
        energy = np.linspace(emin, emax, n)
        edges, widths = cell_edges_from_centers(energy)
    elif grid.spacing == "quadratic":
        x = np.linspace(0.0, 1.0, n)
        energy = emin + (emax - emin) * x * x
        edges, widths = cell_edges_from_centers(energy)
    elif grid.spacing == "log":
        energy = np.geomspace(max(emin, 1.0e-8), emax, n)
        edges, widths = cell_edges_from_centers(energy)
    else:
        raise ValueError(f"Unsupported energy spacing: {grid.spacing}")
    return KineticGrid(energy, edges, widths, electron_speed_m_s(energy))


def _bernoulli(x: np.ndarray | float) -> np.ndarray | float:
    x_arr = np.asarray(x, dtype=float)
    out = np.empty_like(x_arr, dtype=float)
    small = np.abs(x_arr) < 1.0e-7
    pos_large = x_arr > 80.0
    neg_large = x_arr < -80.0
    mid = ~(small | pos_large | neg_large)
    out[small] = (
        1.0
        - x_arr[small] / 2.0
        + x_arr[small] ** 2 / 12.0
        - x_arr[small] ** 4 / 720.0
    )
    out[pos_large] = 0.0
    out[neg_large] = -x_arr[neg_large]
    out[mid] = x_arr[mid] / np.expm1(x_arr[mid])
    if np.isscalar(x):
        return float(out)
    return out


def build_effective_collision_data(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    energy: np.ndarray,
    gas_number_density_m3: float,
) -> EffectiveCollisionData:
    cfg = config.internal.two_term
    speed = electron_speed_m_s(energy)
    sigma_m = np.zeros_like(energy)
    sigma_total_like = np.zeros_like(energy)
    nu_m = np.zeros_like(energy)
    nu_m_over_N = np.zeros_like(energy)
    inelastic_loss = np.zeros_like(energy)
    elastic_A = np.zeros_like(energy)
    elastic_D = np.zeros_like(energy)
    kT_eV = config.conditions.gas_temperature_K * BOLTZMANN_J_K / EV_TO_J
    species_has_effective = {
        species: any(
            proc.process_type == ProcessType.EFFECTIVE
            for proc in cross_sections.by_species(species)
        )
        for species in cross_sections.species
    }
    if not cross_sections.by_type(
        ProcessType.MOMENTUM, ProcessType.EFFECTIVE, ProcessType.ELASTIC
    ):
        raise ValueError("At least one momentum/effective/elastic cross section is required")
    transport_inelastic = {
        ProcessType.EXCITATION,
        ProcessType.IONIZATION,
        ProcessType.ATTACHMENT,
        ProcessType.SUPERELASTIC,
    }
    suppressed: list[str] = []
    for proc in cross_sections.processes:
        frac = mixture_fraction(config.conditions, proc.species)
        if frac <= 0.0:
            continue
        sigma = np.maximum(proc.sigma(energy), cfg.min_momentum_cross_section_m2)
        use_effective_only = species_has_effective.get(proc.species, False)
        include_transport = False
        if proc.process_type == ProcessType.EFFECTIVE:
            include_transport = True
        elif proc.process_type in {ProcessType.MOMENTUM, ProcessType.ELASTIC}:
            include_transport = not use_effective_only
            if use_effective_only:
                suppressed.append(f"{proc.species}:{proc.process_type.value}:{proc.process}")
        elif proc.process_type in transport_inelastic:
            include_transport = not use_effective_only
            inelastic_loss += gas_number_density_m3 * frac * sigma * speed
        if include_transport:
            nuN = frac * sigma * speed
            nu = gas_number_density_m3 * nuN
            sigma_m += frac * sigma
            nu_m_over_N += nuN
            nu_m += nu
        if proc.process_type == ProcessType.EFFECTIVE:
            sigma_total_like += frac * sigma
        elif proc.process_type in {ProcessType.MOMENTUM, ProcessType.ELASTIC}:
            if not use_effective_only:
                sigma_total_like += frac * sigma
        if proc.process_type not in {
            ProcessType.MOMENTUM,
            ProcessType.EFFECTIVE,
            ProcessType.ELASTIC,
        }:
            continue
        mass_amu = proc.mass_amu or gas_mass_amu(config.conditions, proc.species)
        mass_kg = mass_amu * AMU_KG
        nu = gas_number_density_m3 * frac * sigma * speed
        ratio = ELECTRON_MASS_KG / mass_kg
        eps = np.maximum(energy, 1.0e-12)
        elastic_D += 2.0 * ratio * nu * eps * kT_eV
        elastic_A += 2.0 * ratio * nu * (0.5 * kT_eV - eps)
    nu_m = np.maximum(nu_m, 1.0e-60)
    nu_m_over_N = np.maximum(nu_m_over_N, 1.0e-80)
    sigma_m = np.maximum(sigma_m, cfg.min_momentum_cross_section_m2)
    sigma_total_like = np.maximum(sigma_total_like, cfg.min_momentum_cross_section_m2)
    model = "effective_suppresses_elastic_momentum" if suppressed else "direct_momentum"
    return EffectiveCollisionData(
        nu_m=nu_m,
        nu_m_over_N=nu_m_over_N,
        sigma_m=sigma_m,
        sigma_total_like=sigma_total_like,
        inelastic_loss_frequency_s_inv=inelastic_loss,
        elastic_A_eV_s=elastic_A,
        elastic_D_eV2_s=np.maximum(elastic_D, 0.0),
        processes=cross_sections.processes,
        effective_momentum_model=model,
        sigma_total_like_model=model,
        suppressed_momentum_processes=tuple(suppressed),
    )


def field_diffusion_eV2_s(
    energy: np.ndarray,
    electric_field_V_m: float,
    momentum_frequency_s_inv: np.ndarray,
) -> np.ndarray:
    eps_J = np.maximum(energy * EV_TO_J, 0.0)
    D_J2_s = (
        (2.0 / 3.0)
        * (E_CHARGE_C * electric_field_V_m) ** 2
        / ELECTRON_MASS_KG
        * eps_J
        / np.maximum(momentum_frequency_s_inv, 1.0e-60)
    )
    return D_J2_s / (EV_TO_J * EV_TO_J)


def assemble_energy_flux_operator(
    energy: np.ndarray,
    widths: np.ndarray,
    electric_field_V_m: float,
    collisions: EffectiveCollisionData,
) -> sparse.csr_matrix:
    n = len(energy)
    mat = sparse.lil_matrix((n, n), dtype=float)
    D_field = field_diffusion_eV2_s(energy, electric_field_V_m, collisions.nu_m)
    eps_safe = np.maximum(energy, max(energy[1] - energy[0], 1.0e-12) * 0.5)
    A_field = D_field / (2.0 * eps_safe)
    D_center = np.maximum(collisions.elastic_D_eV2_s + D_field, 1.0e-80)
    A_center = collisions.elastic_A_eV_s + A_field
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
    return mat.tocsr()


def _deposit_energy(
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


def _add_energy_shift_transition(
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
        _deposit_energy(mat, energy, widths, i, target_e, multiplicity * rate)


def _add_equal_sharing_ionization(
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
        _deposit_energy(mat, energy, widths, i, target_e, 2.0 * rate)


def _add_primary_secondary_ionization(
    mat: sparse.lil_matrix,
    energy: np.ndarray,
    widths: np.ndarray,
    nu: np.ndarray,
    threshold: float,
    secondary_electron_energy_eV: float,
) -> None:
    secondary_eV = max(float(secondary_electron_energy_eV), 0.0)
    for i, rate in enumerate(np.asarray(nu, dtype=float)):
        if rate <= 0.0:
            continue
        mat[i, i] += -rate
        primary_e = max(0.0, energy[i] - threshold - secondary_eV)
        _deposit_energy(mat, energy, widths, i, primary_e, rate)
        _deposit_energy(mat, energy, widths, i, secondary_eV, rate)


def assemble_collision_operator(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    energy: np.ndarray,
    widths: np.ndarray,
    gas_number_density_m3: float,
) -> sparse.csr_matrix:
    n = len(energy)
    mat = sparse.lil_matrix((n, n), dtype=float)
    speed = electron_speed_m_s(energy)
    cfg = config.internal.two_term
    for proc in cross_sections.processes:
        frac = mixture_fraction(config.conditions, proc.species)
        if frac <= 0.0:
            continue
        if proc.process_type not in {
            ProcessType.EXCITATION,
            ProcessType.SUPERELASTIC,
            ProcessType.IONIZATION,
            ProcessType.ATTACHMENT,
        }:
            continue
        nu = gas_number_density_m3 * frac * proc.sigma(energy) * speed
        if proc.process_type == ProcessType.ATTACHMENT:
            if cfg.nonconservative_model == "ignore":
                continue
            for i, val in enumerate(nu):
                mat[i, i] += -float(val)
            continue
        threshold = float(proc.threshold_eV or 0.0)
        if proc.process_type == ProcessType.SUPERELASTIC and threshold > 0.0:
            threshold = -threshold
        if proc.process_type in {ProcessType.EXCITATION, ProcessType.SUPERELASTIC}:
            _add_energy_shift_transition(
                mat, energy, widths, nu, threshold, multiplicity=1.0
            )
        elif proc.process_type == ProcessType.IONIZATION:
            if (
                cfg.nonconservative_model == "ignore"
                or cfg.ionization_energy_sharing == "loss_only"
            ):
                _add_energy_shift_transition(
                    mat, energy, widths, nu, threshold, multiplicity=1.0
                )
            elif cfg.ionization_energy_sharing == "primary_secondary":
                _add_primary_secondary_ionization(
                    mat,
                    energy,
                    widths,
                    nu,
                    threshold,
                    cfg.secondary_electron_energy_eV,
                )
            else:
                _add_equal_sharing_ionization(mat, energy, widths, nu, threshold)
    return mat.tocsr()


def assemble_native_operator_blocks(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    e_over_n_Td: float,
    grid: KineticGrid,
) -> KineticOperatorBlock:
    N = gas_number_density(config)
    electric_field = e_over_n_Td * TOWNSEND * N
    collisions = build_effective_collision_data(
        config,
        cross_sections,
        grid.energy_eV,
        N,
    )
    energy_flux = assemble_energy_flux_operator(
        grid.energy_eV,
        grid.widths_eV,
        electric_field,
        collisions,
    )
    collision = assemble_collision_operator(
        config,
        cross_sections,
        grid.energy_eV,
        grid.widths_eV,
        N,
    )
    return KineticOperatorBlock(
        grid=grid,
        gas_number_density_m3=float(N),
        electric_field_V_m=float(electric_field),
        collisions=collisions,
        energy_flux_matrix=energy_flux,
        collision_matrix=collision,
        matrix=(energy_flux + collision).tocsr(),
    )


def weighted_integral(values: np.ndarray, widths: np.ndarray) -> float:
    return float(np.sum(np.asarray(values, dtype=float) * widths))


def normalize_eedf(eedf: np.ndarray, widths: np.ndarray) -> np.ndarray:
    values = np.asarray(eedf, dtype=float)
    widths = np.asarray(widths, dtype=float)
    if (
        values.shape != widths.shape
        or not np.all(np.isfinite(values))
        or not np.all(np.isfinite(widths))
        or np.any(widths <= 0.0)
    ):
        raise FloatingPointError("Cannot normalize non-finite EEDF")
    total = weighted_integral(values, widths)
    if not np.isfinite(total) or abs(total) < 1.0e-300:
        raise FloatingPointError("Cannot normalize EEDF with zero/non-finite integral")
    if total < 0.0:
        values = -values
        total = -total
    return values / total


def mean_energy_from_eedf(
    energy_eV: np.ndarray,
    widths_eV: np.ndarray,
    eedf: np.ndarray,
) -> float:
    energy = np.asarray(energy_eV, dtype=float)
    widths = np.asarray(widths_eV, dtype=float)
    values = np.asarray(eedf, dtype=float)
    if (
        energy.shape != widths.shape
        or values.shape != widths.shape
        or not np.all(np.isfinite(energy))
        or not np.all(np.isfinite(widths))
        or not np.all(np.isfinite(values))
        or np.any(widths <= 0.0)
    ):
        raise ValueError("Energy, widths, and EEDF arrays must have matching shape")
    return weighted_integral(energy * values, widths)


def eepf_from_eedf(energy_eV: np.ndarray, eedf: np.ndarray) -> np.ndarray:
    return np.asarray(eedf, dtype=float) / np.sqrt(
        np.maximum(np.asarray(energy_eV, dtype=float), 1.0e-30)
    )


def negative_mass_fraction(eedf: np.ndarray, widths: np.ndarray) -> float:
    values = np.asarray(eedf, dtype=float)
    widths = np.asarray(widths, dtype=float)
    if (
        values.shape != widths.shape
        or not np.all(np.isfinite(values))
        or not np.all(np.isfinite(widths))
        or np.any(widths <= 0.0)
    ):
        raise FloatingPointError("Cannot evaluate negative mass for non-finite EEDF")
    negative = float(np.sum(np.maximum(-values, 0.0) * widths))
    mass = float(np.sum(np.abs(values) * widths))
    return negative / max(mass, 1.0e-300)


def transport_from_reduced(
    e_over_n_Td: float,
    gas_number_density_m3: float,
    muN: float,
    diffN: float,
) -> TransportCoefficients:
    EN = e_over_n_Td * TOWNSEND
    mobility = muN / gas_number_density_m3
    diffusion = diffN / gas_number_density_m3
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


def transport_from_eedf(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    energy: np.ndarray,
    widths: np.ndarray,
    eedf: np.ndarray,
    gas_number_density_m3: float,
    e_over_n_Td: float,
) -> TransportCoefficients:
    coll = build_effective_collision_data(
        config, cross_sections, energy, gas_number_density_m3
    )
    nuN = np.maximum(coll.nu_m_over_N, 1.0e-80)
    speed = electron_speed_m_s(energy)
    if len(energy) >= 3:
        dF = np.gradient(eedf, energy, edge_order=2)
    else:
        dF = np.gradient(eedf, energy)
    mobility_integrand = (2.0 * energy / nuN) * dF - eedf / nuN
    muN = -E_CHARGE_C / (3.0 * ELECTRON_MASS_KG) * weighted_integral(
        mobility_integrand, widths
    )
    if not np.isfinite(muN) or muN <= 0.0:
        nu_eff_over_N = weighted_integral(nuN * eedf, widths)
        muN = E_CHARGE_C / (ELECTRON_MASS_KG * max(nu_eff_over_N, 1.0e-80))
    diffN = (1.0 / 3.0) * weighted_integral((speed * speed / nuN) * eedf, widths)
    return transport_from_reduced(
        e_over_n_Td, gas_number_density_m3, float(muN), float(diffN)
    )


def compute_rates_from_eedf(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    energy: np.ndarray,
    widths: np.ndarray,
    eedf: np.ndarray,
    *,
    case_id: str,
    e_over_n_Td: float,
    solver_name: str,
) -> RateConvolution:
    speed = electron_speed_m_s(energy)
    N = gas_number_density(config)
    rates: list[RateResult] = []
    ion_rate = 0.0
    attach_rate = 0.0
    net_freq = 0.0
    for proc in cross_sections.processes:
        frac = mixture_fraction(config.conditions, proc.species)
        if frac <= 0.0:
            continue
        k = float(np.sum(proc.sigma(energy) * speed * eedf * widths))
        kmix = frac * k
        if proc.process_type == ProcessType.IONIZATION:
            ion_rate += kmix
            net_freq += N * kmix
        elif proc.process_type == ProcessType.ATTACHMENT:
            attach_rate += kmix
            net_freq -= N * kmix
        loss = 0.0
        if proc.process_type in {ProcessType.EXCITATION, ProcessType.IONIZATION}:
            loss = float(proc.threshold_eV or 0.0)
        elif proc.process_type == ProcessType.SUPERELASTIC:
            loss = -abs(float(proc.threshold_eV or 0.0))
        rates.append(
            RateResult(
                solver=solver_name,
                case_id=case_id,
                e_over_n_Td=e_over_n_Td,
                species=proc.species,
                process=proc.process,
                process_type=proc.process_type.value,
                threshold_eV=proc.threshold_eV,
                rate_coefficient_m3_s=k,
                mixture_weighted_rate_m3_s=kmix,
                frequency_s_inv=N * kmix,
                power_loss_eV_s=N * kmix * loss,
            )
        )
    return RateConvolution(rates, ion_rate, attach_rate, net_freq)
