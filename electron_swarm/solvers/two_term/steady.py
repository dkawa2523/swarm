"""Stationary native EEDF solve and adaptive-grid refinement."""

from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as spla

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.constants import TOWNSEND
from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.core.solver_configs import TwoTermInternalConfig
from electron_swarm.physics.kinetics import gas_number_density
from electron_swarm.solvers.boltzmann_common.collisions import (
    assemble_collision_operator,
    build_effective_collision_data,
)
from electron_swarm.solvers.boltzmann_common.grid import (
    KineticGrid,
    cell_edges_from_centers,
    electron_speed_m_s,
)
from electron_swarm.solvers.boltzmann_common.observables import (
    mean_energy_from_eedf,
    normalize_eedf,
    weighted_integral,
)
from electron_swarm.solvers.boltzmann_common.operators import (
    KineticOperatorBlock,
    assemble_energy_flux_operator,
)
from electron_swarm.solvers.two_term.grid import (
    grid_metadata,
    make_energy_grid,
    maxwell_eedf,
)
from electron_swarm.solvers.two_term.models import (
    NativeDistributionResult,
    NativeSolveDiagnostics,
)
from electron_swarm.solvers.two_term.transport import (
    temporal_growth_effective_momentum_frequency,
)


_EPS = 1.0e-300


def solve_native_distribution(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    cfg: TwoTermInternalConfig,
    e_over_n_Td: float,
) -> NativeDistributionResult:
    """Solve the stationary EEDF and retain its final case assembly."""

    max_eV = float(cfg.energy_grid.max_eV)
    previous: tuple[np.ndarray, np.ndarray] | None = None
    last_diag: NativeSolveDiagnostics | None = None
    last_block: KineticOperatorBlock | None = None
    cycles = max(1, cfg.adaptive_grid.max_cycles if cfg.adaptive_grid.enabled else 1)
    for cycle in range(cycles):
        energy, edges, widths = make_energy_grid(
            cfg,
            max_eV_override=max_eV,
            cross_sections=cross_sections,
        )
        initial = None
        if previous is not None:
            old_e, old_f = previous
            initial = np.interp(energy, old_e, old_f, left=0.0, right=0.0)
            initial = normalize_eedf(initial, widths)
        eedf, diag, block = solve_native_on_grid(
            config,
            cross_sections,
            cfg,
            e_over_n_Td,
            energy,
            edges,
            widths,
            initial=initial,
            cycle=cycle,
        )
        previous = (energy, eedf)
        last_diag = diag
        last_block = block
        if not cfg.adaptive_grid.enabled:
            break
        mean_energy = mean_energy_from_eedf(energy, widths, eedf)
        required = max(
            cfg.adaptive_grid.min_max_eV,
            cfg.adaptive_grid.mean_energy_multiplier * mean_energy,
        )
        edge_bad = (
            diag.tail_probability > cfg.adaptive_grid.tail_probability
            or diag.edge_to_peak > cfg.adaptive_grid.edge_to_peak
        )
        if required <= max_eV * 1.05 and not edge_bad:
            break
        max_eV = min(
            cfg.adaptive_grid.max_max_eV,
            max(max_eV * 1.25, required),
        )
        if max_eV <= energy[-1] * 1.001:
            break

    assert previous is not None and last_diag is not None and last_block is not None
    _, eedf = previous
    energy = last_block.energy_eV
    edges = last_block.edges_eV
    widths = last_block.widths_eV
    metadata = {
        "backend": "native_sg",
        "converged": last_diag.converged,
        "iterations": last_diag.iterations,
        "residual_L1": last_diag.residual,
        "residual_requested_tolerance": last_diag.residual_requested_tolerance,
        "residual_roundoff_bound": last_diag.residual_roundoff_bound,
        "residual_backward_error": last_diag.residual_backward_error,
        "residual_roundoff_limited": last_diag.residual_roundoff_limited,
        "growth_frequency_s-1": last_diag.growth_frequency_s,
        "tail_probability": last_diag.tail_probability,
        "tail_probability_target": cfg.adaptive_grid.tail_probability,
        "edge_to_peak": last_diag.edge_to_peak,
        "edge_to_peak_target": cfg.adaptive_grid.edge_to_peak,
        "grid_max_eV": last_diag.grid_max_eV,
        "grid_max_limit_eV": cfg.adaptive_grid.max_max_eV,
        "adaptive_cycles": last_diag.regrid_cycles,
        "residual_tolerance": last_diag.residual_effective_tolerance,
        "residual_acceptance_model": (
            "max_requested_and_componentwise_sparse_matvec_roundoff_bound"
        ),
        "transport_model": "two_term_flux_integral",
        "discretization": "finite_volume_scharfetter_gummel",
        "cross_section_high_energy_extrapolation": (
            config.cross_sections.high_energy_extrapolation
        ),
        "temporal_growth_momentum_correction": {
            "applied": bool(
                config.physics.field.type == "dc"
                and cfg.nonconservative_model == "growth"
            ),
            "model": "nu_m_eff(epsilon)=nu_m(epsilon)+growth_frequency",
            "growth_frequency_s-1": last_diag.growth_frequency_s,
            "frequency_source": (
                "converged_temporal_growth_eigenvalue"
                if last_diag.converged
                else "last_temporal_growth_fixed_point_iterate"
            ),
            "sign_convention": "positive_net_creation_negative_net_loss",
        },
    }
    metadata.update(grid_metadata(cfg, energy))
    return NativeDistributionResult(
        energy_eV=energy,
        edges_eV=edges,
        widths_eV=widths,
        eedf_eV_inv=eedf,
        diagnostics=last_diag,
        metadata=metadata,
        operator_block=last_block,
    )


def assemble_native_operator_block(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    cfg: TwoTermInternalConfig,
    e_over_n_Td: float,
    energy: np.ndarray,
    edges: np.ndarray,
    widths: np.ndarray,
) -> KineticOperatorBlock:
    """Assemble the reusable native two-term energy-space operator."""

    energy = np.asarray(energy, dtype=float)
    edges = np.asarray(edges, dtype=float)
    widths = np.asarray(widths, dtype=float)
    if energy.ndim != 1 or edges.ndim != 1 or widths.ndim != 1:
        raise ValueError("Energy grid arrays must be one-dimensional")
    if len(edges) != len(energy) + 1 or len(widths) != len(energy):
        raise ValueError("Energy grid edges/widths do not match centers")
    if np.any(np.diff(energy) <= 0.0) or np.any(widths <= 0.0):
        raise ValueError("Energy grid must be strictly increasing")

    grid = KineticGrid(energy, edges, widths, electron_speed_m_s(energy))
    density = gas_number_density(config)
    electric_field = e_over_n_Td * TOWNSEND * density
    collisions = build_effective_collision_data(
        config,
        cross_sections,
        grid.energy_eV,
        density,
        cfg,
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
        density,
        cfg,
    )
    return KineticOperatorBlock(
        grid=grid,
        gas_number_density_m3=float(density),
        electric_field_V_m=float(electric_field),
        collisions=collisions,
        energy_flux_matrix=energy_flux,
        collision_matrix=collision,
        matrix=(energy_flux + collision).tocsr(),
    )


def solve_native_on_grid(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    cfg: TwoTermInternalConfig,
    e_over_n_Td: float,
    energy: np.ndarray,
    edges: np.ndarray,
    widths: np.ndarray,
    *,
    initial: np.ndarray | None,
    cycle: int,
) -> tuple[np.ndarray, NativeSolveDiagnostics, KineticOperatorBlock]:
    block = assemble_native_operator_block(
        config,
        cross_sections,
        cfg,
        e_over_n_Td,
        energy,
        edges,
        widths,
    )
    use_pt_momentum = (
        config.physics.field.type == "dc" and cfg.nonconservative_model == "growth"
    )

    if initial is None:
        p = maxwell_eedf(energy, cfg.initial_electron_temperature_eV)
        p = normalize_eedf(p, widths)
    else:
        p = normalize_eedf(np.clip(initial, 0.0, None), widths)

    growth = 0.0
    residual = np.inf
    residual_effective_tolerance = cfg.convergence.residual_tolerance
    residual_roundoff_bound = 0.0
    residual_backward_error = np.inf
    residual_roundoff_limited = False
    converged = False
    iterations = 0
    for it in range(cfg.convergence.max_iterations):
        iterations = it + 1
        if use_pt_momentum:
            effective_momentum = temporal_growth_effective_momentum_frequency(
                block.collisions.nu_m,
                growth,
            )
            energy_flux = assemble_energy_flux_operator(
                energy,
                widths,
                block.electric_field_V_m,
                block.collisions,
                momentum_frequency_s_inv=effective_momentum,
            )
            op = (energy_flux + block.collision_matrix).tocsr()
        else:
            op = block.matrix
        shifted = op - growth * sparse.identity(op.shape[0], format="csr")
        p_new = solve_normalized(shifted, widths)
        if cfg.convergence.clip_negative:
            # The exponential scheme is positivity preserving for the
            # nearest-neighbour fluxes; clipping protects against tiny
            # sparse-solver roundoff and high-energy extrapolation noise.
            p_new = np.clip(p_new, 0.0, None)
            p_new = normalize_eedf(p_new, widths)
        lp = op @ p_new
        new_growth = weighted_integral(lp, widths)
        eig_res = lp - new_growth * p_new
        residual_numerator = float(np.sum(np.abs(eig_res) * widths))
        residual_denominator = max(
            float(np.sum(np.abs(lp) * widths)), abs(new_growth), 1.0
        )
        residual = residual_numerator / residual_denominator
        action_scale = weighted_integral(
            np.asarray(np.abs(op) @ np.abs(p_new)).reshape(-1), widths
        ) + abs(new_growth) * weighted_integral(np.abs(p_new), widths)
        terms_per_row = int(np.max(np.diff(op.indptr))) + 1
        machine_epsilon = np.finfo(float).eps
        gamma = (terms_per_row * machine_epsilon) / (
            1.0 - terms_per_row * machine_epsilon
        )
        residual_roundoff_bound = gamma * action_scale / residual_denominator
        residual_backward_error = residual_numerator / max(action_scale, _EPS)
        residual_effective_tolerance = max(
            cfg.convergence.residual_tolerance,
            residual_roundoff_bound,
        )
        residual_roundoff_limited = bool(
            residual > cfg.convergence.residual_tolerance
            and residual <= residual_effective_tolerance
        )
        shape_change = float(np.sum(np.abs(p_new - p) * widths))
        growth_change = abs(new_growth - growth) / max(
            abs(new_growth), abs(growth), 1.0
        )
        p = p_new
        growth = (
            cfg.convergence.relaxation * new_growth
            + (1.0 - cfg.convergence.relaxation) * growth
        )
        if (
            shape_change < cfg.convergence.tolerance
            and growth_change < cfg.convergence.eigenvalue_tolerance
            and residual <= residual_effective_tolerance
        ):
            converged = True
            growth = new_growth
            break

    tail = tail_probability(p, widths, cfg)
    edge_to_peak = float(p[-1] / max(np.max(p), _EPS))
    diag = NativeSolveDiagnostics(
        converged=converged,
        iterations=iterations,
        residual=residual,
        residual_requested_tolerance=cfg.convergence.residual_tolerance,
        residual_effective_tolerance=residual_effective_tolerance,
        residual_roundoff_bound=residual_roundoff_bound,
        residual_backward_error=residual_backward_error,
        residual_roundoff_limited=residual_roundoff_limited,
        growth_frequency_s=float(growth),
        tail_probability=tail,
        edge_to_peak=edge_to_peak,
        grid_max_eV=float(energy[-1]),
        regrid_cycles=cycle + 1,
    )
    return p, diag, block


def solve_normalized(op: sparse.csr_matrix, widths: np.ndarray) -> np.ndarray:
    n = op.shape[0]
    mat = op.tolil(copy=True)
    rhs = np.zeros(n)
    # Replace one equation by normalization integral F dE = 1. Choose the row
    # with the largest diagonal magnitude to reduce boundary-condition bias.
    diag_abs = np.abs(op.diagonal())
    row = int(np.argmax(diag_abs)) if np.any(diag_abs > 0.0) else n - 1
    mat[row, :] = widths
    rhs[row] = 1.0
    csr = mat.tocsr()
    try:
        sol = spla.spsolve(csr, rhs)
    except Exception:
        sol = spla.lsmr(
            csr,
            rhs,
            atol=1.0e-13,
            btol=1.0e-13,
            maxiter=max(1000, 4 * n),
        )[0]
    if not np.all(np.isfinite(sol)):
        raise FloatingPointError("Boltzmann linear solve returned non-finite values")
    return normalize_eedf(sol, widths)


def tail_probability(
    eedf: np.ndarray,
    widths: np.ndarray,
    cfg: TwoTermInternalConfig,
) -> float:
    n_tail = max(1, int(len(eedf) * cfg.adaptive_grid.tail_cells_fraction))
    return float(np.sum(np.clip(eedf[-n_tail:], 0.0, None) * widths[-n_tail:]))


def mean_energy(
    energy: np.ndarray,
    eedf: np.ndarray,
    widths: np.ndarray | None = None,
) -> float:
    if widths is None:
        widths = cell_edges_from_centers(energy)[1]
    return mean_energy_from_eedf(energy, widths, eedf)
