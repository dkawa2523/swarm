"""Sparse axisymmetric PN operator for homogeneous DC electron swarms.

Even Legendre moments live at energy-cell centers and odd moments at internal
cell faces.  This staggered layout removes the odd-even null mode of a
collocated first derivative and makes the l=0 field coupling conservative.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse

from electron_swarm.core.constants import ELECTRON_MASS_KG, EV_TO_J
from electron_swarm.core.cross_sections import ScatteringRole
from electron_swarm.physics.angular_scattering import build_angular_model
from electron_swarm.solvers.boltzmann_common.collisions import (
    build_effective_collision_data,
    build_inelastic_collision_data,
)
from electron_swarm.solvers.boltzmann_common.grid import electron_speed_m_s
from electron_swarm.solvers.boltzmann_common.operators import (
    assemble_energy_flux_operator,
)

from .models import MultiTermCase, PNOperator


_SPEED_PER_SQRT_EV = np.sqrt(2.0 * EV_TO_J / ELECTRON_MASS_KG)


def _angular_inputs(
    case: MultiTermCase,
    energy_eV: np.ndarray,
    lmax: int,
    angular: object,
) -> tuple[np.ndarray, np.ndarray]:
    """Close each species before adding its linear collision generator."""

    total = np.zeros_like(energy_eV)
    moment_integrals = np.zeros((lmax + 1, len(energy_eV)))
    for component in case.config.conditions.gas_mixture:
        if component.fraction <= 0.0:
            continue
        processes = case.cross_sections.by_species(component.species)
        totals = [
            p for p in processes if p.scattering_role == ScatteringRole.ELASTIC_TOTAL
        ]
        momenta = [
            p
            for p in processes
            if p.scattering_role == ScatteringRole.ELASTIC_MOMENTUM_TRANSFER
        ]
        sigma_total = sum(
            (p.sigma(energy_eV) for p in totals), np.zeros_like(energy_eV)
        )
        sigma_momentum = sum(
            (p.sigma(energy_eV) for p in momenta), np.zeros_like(energy_eV)
        )
        if angular.name in {"momentum_power", "maxent_p1"}:
            if not totals or not momenta:
                raise ValueError(
                    f"{component.species}: {angular.name} requires total and momentum cross sections"
                )
            active = sigma_total > 0.0
            if np.any((~active) & (sigma_momentum > 0.0)):
                raise ValueError(
                    f"{component.species}: momentum cross section exceeds zero total"
                )
            moments = np.zeros_like(moment_integrals)
            moments[0] = 1.0
            moments[:, active] = angular.moments(
                energy_eV[active],
                lmax,
                sigma_total=sigma_total[active],
                sigma_momentum=sigma_momentum[active],
            )
        else:
            moments = angular.moments(energy_eV, lmax)
            if not totals:
                if not momenta:
                    raise ValueError(
                        f"{component.species}: no elastic scattering integral"
                    )
                denominator = 1.0 - moments[1]
                if np.any(denominator <= 1e-8):
                    raise ValueError(
                        "angular moments too close to m1=1 to reconstruct total cross section"
                    )
                sigma_total = sigma_momentum / denominator
            elif momenta and not np.allclose(
                sigma_momentum,
                sigma_total * (1.0 - moments[1]),
                rtol=1e-6,
                atol=1e-30,
            ):
                raise ValueError(
                    f"{component.species}: elastic total and momentum integrals disagree with angular moments; use maxent_p1 for unequal isotropic integrals"
                )
        weighted_total = component.fraction * sigma_total
        total += weighted_total
        moment_integrals += moments * weighted_total

    moments = np.divide(
        moment_integrals,
        total[None, :],
        out=np.zeros_like(moment_integrals),
        where=total[None, :] > 0.0,
    )
    moments[0] = 1.0
    return moments, total


def _center_to_face_derivative(centers: np.ndarray) -> sparse.csr_matrix:
    n = len(centers)
    matrix = sparse.lil_matrix((n - 1, n), dtype=float)
    separation = np.diff(centers)
    rows = np.arange(n - 1)
    matrix[rows, rows] = -1.0 / separation
    matrix[rows, rows + 1] = 1.0 / separation
    return matrix.tocsr()


def _face_to_center_derivative(widths: np.ndarray) -> sparse.csr_matrix:
    """Divergence of internal-face data with zero values at both boundaries."""

    n = len(widths)
    matrix = sparse.lil_matrix((n, n - 1), dtype=float)
    for face in range(n - 1):
        matrix[face, face] += 1.0 / widths[face]
        matrix[face + 1, face] -= 1.0 / widths[face + 1]
    return matrix.tocsr()


def _center_projection(
    centers: np.ndarray,
    edges: np.ndarray,
    ell: int,
) -> sparse.csr_matrix:
    """Map a staggered coefficient to centers for output, not for solving."""

    n = len(centers)
    if ell % 2 == 0:
        return sparse.identity(n, format="csr")
    projection = sparse.lil_matrix((n, n - 1), dtype=float)
    for cell in range(n):
        span = edges[cell + 1] - edges[cell]
        right_weight = (centers[cell] - edges[cell]) / span
        left_weight = 1.0 - right_weight
        if cell > 0:
            projection[cell, cell - 1] = left_weight
        if cell < n - 1:
            projection[cell, cell] = right_weight
    return projection.tocsr()


def build_field_coupling(
    energy_eV: np.ndarray,
    edges_eV: np.ndarray,
    electric_field_V_m: float,
    receiver_l: int,
    source_l: int,
) -> sparse.csr_matrix:
    """Return one physical first-order ``l <-> l +/- 1`` field block."""

    if receiver_l < 0 or abs(receiver_l - source_l) != 1:
        raise ValueError("PN field coupling requires adjacent nonnegative moments")
    centers = np.asarray(energy_eV, dtype=float)
    edges = np.asarray(edges_eV, dtype=float)
    if np.any(centers <= 0.0) or np.any(edges[1:-1] <= 0.0):
        raise ValueError("PN field coupling requires positive energy coordinates")
    widths = np.diff(edges)
    faces = edges[1:-1]
    ell = receiver_l
    scale = -float(electric_field_V_m) * _SPEED_PER_SQRT_EV

    if ell % 2 == 0:
        derivative = _face_to_center_derivative(widths)
        receiver_energy = centers
        source_energy = faces
    else:
        derivative = _center_to_face_derivative(centers)
        receiver_energy = faces
        source_energy = centers

    if source_l == ell - 1:
        left_power = (ell + 1.0) / 2.0
        right_power = -ell / 2.0
    else:
        left_power = -ell / 2.0
        right_power = (ell + 1.0) / 2.0
    block = (
        scale
        * sparse.diags(receiver_energy**left_power, format="csr")
        @ derivative
        @ sparse.diags(source_energy**right_power, format="csr")
    ).tocsr()
    if block.data.size and not np.all(np.isfinite(block.data)):
        raise FloatingPointError("PN field coupling contains non-finite entries")
    return block


def build_angular_damping(
    case: MultiTermCase,
    energy_eV: np.ndarray,
    speed_m_s: np.ndarray,
    sigma_total_m2: np.ndarray,
    inelastic_loss_s_inv: np.ndarray,
    moments: np.ndarray,
) -> np.ndarray:
    """Build the collision eigenvalue ``nu_l`` for every angular moment."""

    lmax = moments.shape[0] - 1
    shape = np.asarray(energy_eV).shape
    if lmax < 1 or moments.shape[1:] != shape:
        raise ValueError("angular moments must have shape (lmax + 1, n_energy)")
    if not np.all(np.isfinite(moments)) or np.any(np.abs(moments) > 1.0):
        raise ValueError("angular moments must be finite and within [-1, 1]")
    if not np.allclose(moments[0], 1.0, rtol=1.0e-10, atol=1.0e-10):
        raise ValueError("angular moments must have m0=1")
    sigma_total = np.asarray(sigma_total_m2, dtype=float)
    if sigma_total.shape != shape or np.any(~np.isfinite(sigma_total)):
        raise ValueError("elastic total cross section must match the PN energy grid")
    if np.any(sigma_total < 0.0):
        raise ValueError(
            "elastic total cross section must be nonnegative for PN damping"
        )

    damping = np.zeros_like(moments)
    base = case.gas_number_density_m3 * np.asarray(speed_m_s, dtype=float) * sigma_total
    loss = np.asarray(inelastic_loss_s_inv, dtype=float)
    for ell in range(1, lmax + 1):
        damping[ell] = np.maximum(
            base * (1.0 - moments[ell]) + loss,
            case.gas_number_density_m3
            * speed_m_s
            * case.multi_term_config.min_momentum_cross_section_m2,
        )
    if not np.all(np.isfinite(damping[1:])) or np.any(damping[1:] <= 0.0):
        raise ValueError("PN angular damping must be finite and positive")
    return damping


def assemble_pn_operator(case: MultiTermCase) -> PNOperator:
    """Assemble the complete ``F_0 .. F_L`` stationary generator."""

    if case.config.physics.field.type != "dc":
        raise NotImplementedError("axisymmetric multi_term supports dc fields only")
    magnetic = case.config.physics.field.magnetic_field
    if magnetic.enabled and magnetic.B_T != 0.0:
        raise NotImplementedError(
            "axisymmetric m=0 multi_term does not support a magnetic field"
        )
    lmax = int(case.multi_term_config.lmax)
    if lmax < 1:
        raise ValueError("multi_term lmax must be at least 1")

    energy = np.asarray(case.grid.centers_eV, dtype=float)
    edges = np.asarray(case.grid.edges_eV, dtype=float)
    widths = np.asarray(case.grid.widths_eV, dtype=float)
    speed = electron_speed_m_s(energy)
    collisions = build_effective_collision_data(
        case.config,
        case.cross_sections,
        energy,
        case.gas_number_density_m3,
        case.multi_term_config,
    )
    inelastic = build_inelastic_collision_data(
        case.config,
        case.cross_sections,
        energy,
        widths,
        case.gas_number_density_m3,
        case.multi_term_config,
    )
    # Only neutral thermal recoil/elastic energy relaxation is reused here.
    # The eliminated two-term electric-field diffusion is deliberately absent.
    elastic_energy = assemble_energy_flux_operator(
        energy,
        widths,
        0.0,
        collisions,
    )
    isotropic_collision = (inelastic.matrix + elastic_energy).tocsr()

    method = case.multi_term_config.method
    angular = build_angular_model(case.config)
    if method == "pn_closure_direct" and angular.name == "moment_table":
        raise NotImplementedError(
            "pn_closure_direct uses ordinary integral cross sections; "
            "use pn_dcs for a moment table"
        )
    if method == "pn_dcs" and angular.name != "moment_table":
        raise NotImplementedError(
            "pn_dcs requires physics.angular_scattering.model=moment_table"
        )
    moments, elastic_total = _angular_inputs(
        case,
        energy,
        lmax,
        angular,
    )
    damping = build_angular_damping(
        case,
        energy,
        speed,
        elastic_total,
        collisions.inelastic_loss_frequency_s_inv,
        moments,
    )

    sizes = tuple(
        len(energy) if ell % 2 == 0 else len(energy) - 1 for ell in range(lmax + 1)
    )
    starts = np.cumsum((0, *sizes))
    slices = tuple(
        slice(int(starts[ell]), int(starts[ell + 1])) for ell in range(lmax + 1)
    )
    moment_energy = tuple(
        energy if ell % 2 == 0 else edges[1:-1] for ell in range(lmax + 1)
    )
    moment_weights = tuple(
        widths if ell % 2 == 0 else 0.5 * (edges[2:] - edges[:-2])
        for ell in range(lmax + 1)
    )
    projections = tuple(
        _center_projection(energy, edges, ell) for ell in range(lmax + 1)
    )
    blocks: list[list[sparse.csr_matrix]] = [
        [
            sparse.csr_matrix((sizes[row], sizes[column]), dtype=float)
            for column in range(lmax + 1)
        ]
        for row in range(lmax + 1)
    ]
    blocks[0][0] = isotropic_collision
    for ell in range(1, lmax + 1):
        damping_values = (
            damping[ell]
            if ell % 2 == 0
            else np.interp(edges[1:-1], energy, damping[ell])
        )
        blocks[ell][ell] = sparse.diags(-damping_values, format="csr")
    for ell in range(lmax + 1):
        if ell > 0:
            lower = build_field_coupling(
                energy, edges, case.electric_field_V_m, ell, ell - 1
            )
            blocks[ell][ell - 1] = (ell / (2.0 * ell - 1.0)) * lower
        if ell < lmax:
            upper = build_field_coupling(
                energy, edges, case.electric_field_V_m, ell, ell + 1
            )
            blocks[ell][ell + 1] = ((ell + 1.0) / (2.0 * ell + 3.0)) * upper

    matrix = sparse.bmat(blocks, format="csr")
    expected_dimension = int(sum(sizes))
    if matrix.shape != (expected_dimension, expected_dimension) or (
        matrix.data.size and not np.all(np.isfinite(matrix.data))
    ):
        raise FloatingPointError("assembled PN operator is invalid")
    angular_metadata = dict(angular.metadata())
    angular_metadata.update(
        {
            "ordinary_integral_xs_closure": method == "pn_closure_direct",
            "exact_dcs_based": bool(
                method == "pn_dcs" and angular_metadata.get("exact_dcs_based", False)
            ),
        }
    )
    return PNOperator(
        matrix=matrix,
        energy_eV=energy,
        edges_eV=edges,
        widths_eV=widths,
        speed_m_s=speed,
        angular_moments=moments,
        angular_damping_s_inv=damping,
        sigma_m_m2=(
            damping[1] / np.maximum(case.gas_number_density_m3 * speed, 1.0e-300)
        ),
        angular_metadata=angular_metadata,
        lmax=lmax,
        moment_slices=slices,
        moment_energy_eV=moment_energy,
        moment_weights_eV=moment_weights,
        center_projections=projections,
    )


__all__ = [
    "assemble_pn_operator",
    "build_angular_damping",
    "build_field_coupling",
]
