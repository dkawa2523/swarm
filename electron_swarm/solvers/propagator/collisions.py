"""Conservative collision blocks for the stationary response propagator.

Elastic scattering is split into a reciprocal angular collision and a
finite-temperature, first-mass-ratio energy Fokker--Planck block.  The latter
is an equilibrium-fitted Scharfetter--Gummel generator; it replaces the old
one-way mean-recoil remap.  Inelastic transfers are integrated over donor
energy cells and weakly deposited into destination cells.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.special import erfcx

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.constants import (
    AMU_KG,
    BOLTZMANN_J_K,
    ELECTRON_MASS_KG,
    EV_TO_J,
)
from electron_swarm.core.cross_sections import (
    CrossSectionProcess,
    CrossSectionSet,
    ProcessType,
    gas_mass_amu,
    mixture_fraction,
)
from electron_swarm.core.scattering import resolve_particle_scattering
from electron_swarm.physics.electron_neutral import (
    ionization_daughter_energy_arrays,
    speed_from_energy_m_s,
)
from electron_swarm.physics.kinetics import gas_number_density
from electron_swarm.solvers.propagator.angular_kernel import (
    isotropic_destination_weights,
    maxent_p1_transition_stack,
)
from electron_swarm.solvers.propagator.models import (
    CollisionOperatorData,
    ElasticTransfer,
    InelasticTransfer,
    PropagatorGrid,
)


_CELL_QUADRATURE_ORDER = 8


def _process_id(species: str, process: str) -> str:
    return f"{species}:{process}"


def _bernoulli(value: float) -> float:
    if abs(value) < 1.0e-7:
        return float(
            1.0
            - value / 2.0
            + value * value / 12.0
            - value**4 / 720.0
        )
    if value > 80.0:
        return 0.0
    if value < -80.0:
        return -value
    return float(value / np.expm1(value))


def _maxwell_log_cell_mass(edges_eV: np.ndarray, kT_eV: float) -> np.ndarray:
    """Return log cell integrals of ``sqrt(E) exp(-E/kT)`` without underflow."""

    scaled = np.asarray(edges_eV, dtype=float) / float(kT_eV)
    root = np.sqrt(np.maximum(scaled, 0.0))
    log_upper_gamma = -scaled + np.log(
        root + 0.5 * np.sqrt(np.pi) * erfcx(root)
    )
    log_ratio = log_upper_gamma[1:] - log_upper_gamma[:-1]
    if np.any(log_ratio >= 0.0):
        raise FloatingPointError("Maxwell cell measure is not strictly decreasing")
    return (
        1.5 * np.log(float(kT_eV))
        + log_upper_gamma[:-1]
        + np.log(-np.expm1(log_ratio))
    )


def cell_average_collision_rate(
    processes: tuple[CrossSectionProcess, ...],
    density_m3: float,
    fraction: float,
    grid: PropagatorGrid,
) -> np.ndarray:
    result = np.zeros(grid.energy_cells, dtype=float)
    for index, (lower, upper) in enumerate(
        zip(grid.energy_edges_eV[:-1], grid.energy_edges_eV[1:], strict=True)
    ):
        energy, quadrature_weights = _cell_average_quadrature(
            float(lower),
            float(upper),
            processes,
        )
        sigma = np.zeros_like(energy)
        for process in processes:
            process_sigma = process.sigma(energy)
            threshold = process.threshold_eV
            if threshold is not None:
                process_sigma = np.where(
                    energy >= float(threshold),
                    process_sigma,
                    0.0,
                )
            sigma += process_sigma
        rate = density_m3 * fraction * sigma * speed_from_energy_m_s(energy)
        result[index] = float(np.dot(quadrature_weights, rate))
    return result


def _cell_average_quadrature(
    lower_eV: float,
    upper_eV: float,
    processes: tuple[CrossSectionProcess, ...],
    *,
    kinematic_breaks_eV: tuple[float, ...] = (),
) -> tuple[np.ndarray, np.ndarray]:
    """Return a cell-average quadrature aligned with every physical kink.

    Cross sections are piecewise linear, while reaction thresholds and some
    ionization sharing laws introduce additional derivative discontinuities.
    Integrating across those points with one fixed Gaussian rule makes the
    weak collision operator depend on their accidental position inside a
    donor cell.  Split first, then apply the same positive Gaussian rule on
    each smooth interval.  The returned weights include division by the full
    donor-cell width and therefore sum to one.
    """

    lower = float(lower_eV)
    upper = float(upper_eV)
    if not np.isfinite(lower) or not np.isfinite(upper) or upper <= lower:
        raise ValueError("collision quadrature requires an increasing cell")
    breaks: list[float] = [lower, upper]
    for process in processes:
        table = np.asarray(process.energy_eV, dtype=float)
        first = int(np.searchsorted(table, lower, side="right"))
        last = int(np.searchsorted(table, upper, side="left"))
        if last > first:
            breaks.extend(map(float, table[first:last]))
        if process.threshold_eV is not None:
            threshold = float(process.threshold_eV)
            if lower < threshold < upper:
                breaks.append(threshold)
    breaks.extend(
        value
        for value in map(float, kinematic_breaks_eV)
        if lower < value < upper
    )
    partition = np.unique(np.asarray(breaks, dtype=float))
    nodes, weights = np.polynomial.legendre.leggauss(_CELL_QUADRATURE_ORDER)
    energies: list[np.ndarray] = []
    average_weights: list[np.ndarray] = []
    width = upper - lower
    for left, right in zip(partition[:-1], partition[1:], strict=True):
        half_width = 0.5 * float(right - left)
        energies.append(half_width * nodes + 0.5 * float(right + left))
        average_weights.append(weights * half_width / width)
    return np.concatenate(energies), np.concatenate(average_weights)


def _deposit_target(
    rows: list[int],
    columns: list[int],
    entries: list[float],
    centers: np.ndarray,
    source: int,
    target_eV: float,
    amount_s_inv: float,
) -> None:
    if amount_s_inv == 0.0:
        return
    if target_eV <= centers[0]:
        rows.append(0)
        columns.append(source)
        entries.append(amount_s_inv)
        return
    if target_eV >= centers[-1]:
        rows.append(len(centers) - 1)
        columns.append(source)
        entries.append(amount_s_inv)
        return
    upper = int(np.searchsorted(centers, target_eV, side="right"))
    lower = upper - 1
    upper_weight = (target_eV - centers[lower]) / (
        centers[upper] - centers[lower]
    )
    rows.extend((lower, upper))
    columns.extend((source, source))
    entries.extend(
        ((1.0 - upper_weight) * amount_s_inv, upper_weight * amount_s_inv)
    )


def _weak_inelastic_transfer(
    config: SwarmConfig,
    process: CrossSectionProcess,
    density_m3: float,
    fraction: float,
    grid: PropagatorGrid,
) -> InelasticTransfer:
    centers = grid.energy_centers_eV
    frequency = np.zeros(grid.energy_cells, dtype=float)
    incident_energy_rate = np.zeros(grid.energy_cells, dtype=float)
    daughter_energy_rate = np.zeros(grid.energy_cells, dtype=float)
    rows: list[int] = []
    columns: list[int] = []
    entries: list[float] = []
    threshold = (
        0.0
        if process.process_type == ProcessType.ATTACHMENT
        and process.threshold_eV is None
        else float(process.threshold_eV)
    )
    daughter_count = {
        ProcessType.ATTACHMENT: 0,
        ProcessType.EXCITATION: 1,
        ProcessType.IONIZATION: (
            1
            if config.physics.ionization.energy_sharing == "loss_only"
            else 2
        ),
    }[process.process_type]
    sharing_breaks: tuple[float, ...] = ()
    if (
        process.process_type == ProcessType.IONIZATION
        and config.physics.ionization.energy_sharing == "primary_secondary"
    ):
        sharing_breaks = (
            threshold
            + float(
                config.physics.ionization.secondary_electron_energy_eV
            ),
        )

    for source, (lower, upper) in enumerate(
        zip(grid.energy_edges_eV[:-1], grid.energy_edges_eV[1:], strict=True)
    ):
        energy, quadrature_weights = _cell_average_quadrature(
            float(lower),
            float(upper),
            (process,),
            kinematic_breaks_eV=sharing_breaks,
        )
        sigma = np.where(
            energy >= threshold,
            process.sigma(energy),
            0.0,
        )
        rate = (
            density_m3
            * fraction
            * sigma
            * speed_from_energy_m_s(energy)
        )
        quadrature_rate = quadrature_weights * rate
        frequency[source] = float(np.sum(quadrature_rate))
        incident_energy_rate[source] = float(
            np.dot(quadrature_rate, energy)
        )
        if daughter_count == 0:
            continue
        if process.process_type == ProcessType.EXCITATION:
            daughters = (np.maximum(energy - threshold, 0.0),)
        else:
            daughters = ionization_daughter_energy_arrays(
                energy,
                threshold,
                model=config.physics.ionization.energy_sharing,
                secondary_electron_energy_eV=(
                    config.physics.ionization.secondary_electron_energy_eV
                ),
            )
        for node_index, event_rate in enumerate(quadrature_rate):
            for daughter in daughters:
                daughter_energy_rate[source] += (
                    float(event_rate) * float(daughter[node_index])
                )
                _deposit_target(
                    rows,
                    columns,
                    entries,
                    centers,
                    source,
                    float(daughter[node_index]),
                    float(event_rate),
                )

    gain = sparse.coo_matrix(
        (entries, (rows, columns)),
        shape=(grid.energy_cells, grid.energy_cells),
        dtype=float,
    ).tocsr()
    gain.sum_duplicates()
    expected = daughter_count * frequency
    actual = np.asarray(gain.sum(axis=0)).reshape(-1)
    if not np.allclose(actual, expected, rtol=2.0e-13, atol=1.0e-12):
        raise FloatingPointError(
            f"{process.species}:{process.process} weak transfer lost daughters"
        )
    if process.process_type in {
        ProcessType.EXCITATION,
        ProcessType.IONIZATION,
    }:
        energy_loss_rate = incident_energy_rate - daughter_energy_rate
        expected_loss_rate = threshold * frequency
        energy_scale = np.maximum(
            np.maximum(
                np.abs(incident_energy_rate),
                np.abs(daughter_energy_rate),
            ),
            np.abs(expected_loss_rate),
        )
        if np.any(
            np.abs(energy_loss_rate - expected_loss_rate)
            > 3.0e-13 * energy_scale + 1.0e-12
        ):
            raise FloatingPointError(
                f"{process.species}:{process.process} weak transfer violated "
                "the reaction energy identity"
            )
    deposited_energy_rate = np.asarray(centers @ gain).reshape(-1)
    boundary_energy_radius = max(
        float(centers[0]),
        float(grid.energy_edges_eV[-1] - centers[-1]),
    )
    projection_bound = expected * boundary_energy_radius
    if np.any(
        np.abs(deposited_energy_rate - daughter_energy_rate)
        > projection_bound * (1.0 + 2.0e-13) + 1.0e-12
    ):
        raise FloatingPointError(
            f"{process.species}:{process.process} weak energy projection "
            "exceeded the finite-volume boundary bound"
        )
    return InelasticTransfer(
        species=process.species,
        process=process.process,
        process_type=process.process_type.value,
        frequency_s_inv=frequency,
        incident_energy_rate_eV_s_inv=incident_energy_rate,
        daughter_energy_rate_eV_s_inv=daughter_energy_rate,
        gain_matrix_s_inv=gain,
        daughter_count=daughter_count,
        number_change_per_event=daughter_count - 1,
        threshold_eV=threshold,
    )


def _elastic_energy_generator(
    config: SwarmConfig,
    grid: PropagatorGrid,
    elastic: tuple[ElasticTransfer, ...],
) -> tuple[sparse.csr_matrix, np.ndarray, np.ndarray, np.ndarray]:
    """Build a momentum-transfer-driven reversible SG energy generator."""

    energy = grid.energy_centers_eV
    widths = grid.energy_widths_eV
    kT_eV = config.conditions.gas_temperature_K * BOLTZMANN_J_K / EV_TO_J
    if not np.isfinite(kT_eV) or kT_eV <= 0.0:
        raise ValueError("propagator requires a finite positive gas temperature")
    drift = np.zeros_like(energy)
    diffusion = np.zeros_like(energy)
    for transfer in elastic:
        ratio = ELECTRON_MASS_KG / (transfer.target_mass_amu * AMU_KG)
        drift += (
            2.0
            * ratio
            * transfer.momentum_frequency_s_inv
            * (0.5 * kT_eV - energy)
        )
        diffusion += (
            2.0
            * ratio
            * transfer.momentum_frequency_s_inv
            * energy
            * kT_eV
        )

    density_operator = sparse.lil_matrix(
        (grid.energy_cells, grid.energy_cells), dtype=float
    )
    log_equilibrium_mass = _maxwell_log_cell_mass(
        grid.energy_edges_eV,
        kT_eV,
    )
    log_equilibrium = log_equilibrium_mass - np.log(widths)
    for index in range(grid.energy_cells - 1):
        spacing = float(energy[index + 1] - energy[index])
        edge_diffusion = 0.5 * (
            float(diffusion[index]) + float(diffusion[index + 1])
        )
        if spacing <= 0.0 or edge_diffusion <= 0.0:
            continue
        peclet = float(log_equilibrium[index + 1] - log_equilibrium[index])
        left = edge_diffusion / spacing * _bernoulli(-peclet)
        right = -edge_diffusion / spacing * _bernoulli(peclet)
        density_operator[index, index] += -left / widths[index]
        density_operator[index, index + 1] += -right / widths[index]
        density_operator[index + 1, index] += left / widths[index + 1]
        density_operator[index + 1, index + 1] += right / widths[index + 1]

    density_csr = density_operator.tocsr()
    mass_operator = (
        sparse.diags(widths)
        @ density_csr
        @ sparse.diags(1.0 / widths)
    ).tocsr()
    mass_operator.eliminate_zeros()
    column_error = float(
        np.max(np.abs(np.asarray(mass_operator.sum(axis=0)).reshape(-1)))
    )
    scale = max(float(np.max(np.abs(mass_operator.data), initial=0.0)), 1.0)
    if column_error > 2.0e-12 * scale:
        raise FloatingPointError("elastic SG energy generator is not conservative")
    diagonal = mass_operator.diagonal()
    off_diagonal = mass_operator - sparse.diags(diagonal)
    if off_diagonal.nnz and float(np.min(off_diagonal.data)) < -2.0e-14 * scale:
        raise FloatingPointError("elastic SG energy generator is not positive")
    equilibrium_mass = np.exp(
        log_equilibrium_mass - float(np.max(log_equilibrium_mass))
    )
    residual = mass_operator @ equilibrium_mass
    operator_one_norm = float(
        np.max(
            np.asarray(np.abs(mass_operator).sum(axis=0)).reshape(-1),
            initial=0.0,
        )
    )
    equilibrium_error = float(np.linalg.norm(residual, ord=1)) / max(
        operator_one_norm * float(np.linalg.norm(equilibrium_mass, ord=1)),
        1.0e-300,
    )
    if equilibrium_error > 2.0e-12:
        raise FloatingPointError(
            "elastic SG energy generator does not preserve Maxwell equilibrium"
        )
    equilibrium_mass /= float(np.sum(equilibrium_mass))
    return mass_operator, drift, diffusion, equilibrium_mass


def build_collision_operator(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    grid: PropagatorGrid,
) -> CollisionOperatorData:
    angular_model = config.physics.angular_scattering.model
    if angular_model not in {"isotropic", "maxent_p1"}:
        raise NotImplementedError(
            f"propagator angular model {angular_model!r} is not implemented"
        )
    density = gas_number_density(config)
    isotropic_weights = isotropic_destination_weights(grid)
    elastic_maps: list[ElasticTransfer] = []
    inelastic_maps: list[InelasticTransfer] = []
    roles: list[str] = []
    total_ids: list[str] = []
    momentum_ids: list[str] = []

    for channels in resolve_particle_scattering(
        cross_sections,
        angular_model=angular_model,
        active_species=(
            component.species
            for component in config.conditions.gas_mixture
            if component.fraction > 0.0
        ),
    ):
        fraction = mixture_fraction(config.conditions, channels.species)
        if fraction <= 0.0:
            continue
        collision_processes = tuple(channels.collision_total)
        collision_frequency = cell_average_collision_rate(
            collision_processes,
            density,
            fraction,
            grid,
        )
        for process in collision_processes:
            total_ids.append(_process_id(process.species, process.process))
        target_mass = gas_mass_amu(config.conditions, channels.species)
        explicit_momentum_processes = tuple(channels.momentum_transfer)
        momentum_processes = (
            explicit_momentum_processes
            if explicit_momentum_processes
            else collision_processes
        )
        momentum_frequency = cell_average_collision_rate(
            momentum_processes,
            density,
            fraction,
            grid,
        )
        for process in explicit_momentum_processes:
            momentum_ids.append(_process_id(process.species, process.process))
        if angular_model == "maxent_p1":
            mean_cosine = 1.0 - momentum_frequency / np.maximum(
                collision_frequency, 1.0e-300
            )
            active = collision_frequency > 0.0
            if np.any(
                active
                & (
                    (mean_cosine < -1.0 - 1.0e-10)
                    | (mean_cosine > 1.0 + 1.0e-10)
                )
            ):
                raise ValueError(
                    f"{channels.species}: cell-integrated total and momentum "
                    "cross sections imply an angular moment outside [-1, 1]"
                )
            mean_cosine = np.where(
                active,
                np.clip(mean_cosine, -1.0, 1.0),
                0.0,
            )
            angular_kernel = maxent_p1_transition_stack(
                grid,
                mean_cosine,
            )
            isotropic = False
        else:
            if explicit_momentum_processes:
                relative_mismatch = np.abs(
                    collision_frequency - momentum_frequency
                ) / np.maximum(
                    np.maximum(collision_frequency, momentum_frequency),
                    1.0e-300,
                )
                if float(np.max(relative_mismatch)) > 1.0e-8:
                    raise ValueError(
                        f"{channels.species}: isotropic scattering requires "
                        "equal elastic total and momentum-transfer cross "
                        "sections; use maxent_p1 for anisotropic inputs"
                    )
                collision_frequency = momentum_frequency.copy()
            else:
                momentum_frequency = collision_frequency.copy()
            mean_cosine = np.zeros(grid.energy_cells, dtype=float)
            angular_kernel = None
            isotropic = True
        elastic_maps.append(
            ElasticTransfer(
                species=channels.species,
                collision_frequency_s_inv=collision_frequency,
                momentum_frequency_s_inv=momentum_frequency,
                mean_cosine=mean_cosine,
                angular_kernel=angular_kernel,
                isotropic=isotropic,
                target_mass_amu=target_mass,
                xs_role=channels.collision_total_model,
                collision_processes=collision_processes,
                momentum_processes=momentum_processes,
                density_scale_m3=density * fraction,
            )
        )
        roles.append(channels.collision_total_model)

    supported_inelastic = {
        ProcessType.EXCITATION,
        ProcessType.IONIZATION,
        ProcessType.ATTACHMENT,
    }
    for process in cross_sections.processes:
        if process.process_type in {
            ProcessType.ELASTIC,
            ProcessType.MOMENTUM,
            ProcessType.EFFECTIVE,
        }:
            continue
        fraction = mixture_fraction(config.conditions, process.species)
        if fraction <= 0.0:
            continue
        if process.process_type == ProcessType.SUPERELASTIC:
            raise NotImplementedError(
                "propagator P1 does not implement superelastic collisions"
            )
        if process.process_type not in supported_inelastic:
            raise NotImplementedError(
                f"propagator does not implement process type "
                f"{process.process_type.value!r}"
            )
        inelastic_maps.append(
            _weak_inelastic_transfer(
                config,
                process,
                density,
                fraction,
                grid,
            )
        )

    if not elastic_maps:
        raise ValueError("propagator requires an elastic scattering channel")
    elastic_tuple = tuple(elastic_maps)
    elastic_energy, elastic_A, elastic_D, elastic_equilibrium = _elastic_energy_generator(
        config,
        grid,
        elastic_tuple,
    )
    elastic_diagonal = elastic_energy.diagonal()
    elastic_gain = (
        elastic_energy - sparse.diags(elastic_diagonal)
    ).tocsr()
    elastic_outflow = -elastic_diagonal
    inelastic_outflow = np.zeros(grid.energy_cells, dtype=float)
    for item in inelastic_maps:
        inelastic_outflow += item.frequency_s_inv
    nonlocal_outflow = elastic_outflow + inelastic_outflow

    memory = (
        isotropic_weights.nbytes
        + elastic_A.nbytes
        + elastic_D.nbytes
        + elastic_equilibrium.nbytes
        + elastic_outflow.nbytes
        + inelastic_outflow.nbytes
        + nonlocal_outflow.nbytes
        + elastic_energy.data.nbytes
        + elastic_energy.indices.nbytes
        + elastic_energy.indptr.nbytes
        + sum(
            item.collision_frequency_s_inv.nbytes
            + item.momentum_frequency_s_inv.nbytes
            + item.mean_cosine.nbytes
            + (0 if item.angular_kernel is None else item.angular_kernel.nbytes)
            for item in elastic_tuple
        )
        + sum(
            item.frequency_s_inv.nbytes
            + item.incident_energy_rate_eV_s_inv.nbytes
            + item.daughter_energy_rate_eV_s_inv.nbytes
            + item.gain_matrix_s_inv.data.nbytes
            + item.gain_matrix_s_inv.indices.nbytes
            + item.gain_matrix_s_inv.indptr.nbytes
            for item in inelastic_maps
        )
    )
    return CollisionOperatorData(
        elastic=elastic_tuple,
        inelastic=tuple(inelastic_maps),
        elastic_energy_generator_s_inv=elastic_energy,
        elastic_energy_gain_s_inv=elastic_gain,
        elastic_energy_outflow_s_inv=elastic_outflow,
        inelastic_outflow_s_inv=inelastic_outflow,
        nonlocal_outflow_s_inv=nonlocal_outflow,
        elastic_A_eV_s=elastic_A,
        elastic_D_eV2_s=elastic_D,
        elastic_equilibrium_mass=elastic_equilibrium,
        isotropic_weights=isotropic_weights,
        memory_bytes=int(memory),
        elastic_xs_roles=tuple(roles),
        elastic_total_process_ids=tuple(sorted(set(total_ids))),
        elastic_momentum_process_ids=tuple(sorted(set(momentum_ids))),
    )


def angular_collision_loss_matrix(
    operator: CollisionOperatorData,
    energy_index: int,
    mu_widths: np.ndarray,
) -> np.ndarray:
    """Return loss-minus-gain acting on angular density in one shell."""

    widths = np.asarray(mu_widths, dtype=float)
    result = np.zeros((len(widths), len(widths)), dtype=float)
    for transfer in operator.elastic:
        frequency = float(transfer.collision_frequency_s_inv[energy_index])
        if transfer.isotropic:
            probability = np.repeat(
                operator.isotropic_weights[:, None], len(widths), axis=1
            )
        else:
            assert transfer.angular_kernel is not None
            probability = transfer.angular_kernel[energy_index]
        density_transition = (
            probability
            * widths[None, :]
            / widths[:, None]
        )
        result += frequency * (np.eye(len(widths)) - density_transition)
    weighted_error = float(np.max(np.abs(widths @ result)))
    scale = max(float(np.linalg.norm(result, ord=1)), 1.0)
    off_diagonal = result.copy()
    np.fill_diagonal(off_diagonal, 0.0)
    if weighted_error > 2.0e-12 * scale or np.max(off_diagonal) > 2.0e-14 * scale:
        raise FloatingPointError("angular collision block is not a conservative M-matrix")
    return result


def nonlocal_collision_inflow(
    operator: CollisionOperatorData,
    population: np.ndarray,
) -> np.ndarray:
    values = np.asarray(population, dtype=float)
    if values.shape != (
        operator.nonlocal_outflow_s_inv.size,
        operator.isotropic_weights.size,
    ):
        raise ValueError("population shape does not match collision operator")
    out = np.asarray(operator.elastic_energy_gain_s_inv @ values, dtype=float)
    shell_population = np.sum(values, axis=1)
    for transfer in operator.inelastic:
        daughter_rate = transfer.gain_matrix_s_inv @ shell_population
        out += daughter_rate[:, None] * operator.isotropic_weights[None, :]
    return out


def collision_inflow(
    operator: CollisionOperatorData,
    population: np.ndarray,
    output: np.ndarray | None = None,
) -> np.ndarray:
    values = np.asarray(population, dtype=float)
    out = nonlocal_collision_inflow(operator, values)
    for transfer in operator.elastic:
        source = transfer.collision_frequency_s_inv[:, None] * values
        if transfer.isotropic:
            redirected = (
                np.sum(source, axis=1)[:, None]
                * operator.isotropic_weights[None, :]
            )
        else:
            assert transfer.angular_kernel is not None
            redirected = np.einsum(
                "eij,ej->ei", transfer.angular_kernel, source, optimize=True
            )
        out += redirected
    if output is not None:
        output[:] = out
        return output
    return out


def collision_generator_action(
    operator: CollisionOperatorData,
    population: np.ndarray,
) -> np.ndarray:
    values = np.asarray(population, dtype=float)
    angular_outflow = np.zeros_like(values)
    for transfer in operator.elastic:
        angular_outflow += transfer.collision_frequency_s_inv[:, None] * values
    return collision_inflow(operator, values) - (
        angular_outflow
        + operator.nonlocal_outflow_s_inv[:, None] * values
    )


def reaction_number_rate_s_inv(
    operator: CollisionOperatorData,
    population: np.ndarray,
) -> float:
    shell_population = np.sum(np.asarray(population, dtype=float), axis=1)
    return float(
        sum(
            transfer.number_change_per_event
            * float(np.dot(transfer.frequency_s_inv, shell_population))
            for transfer in operator.inelastic
        )
    )
