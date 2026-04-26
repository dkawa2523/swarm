"""Sparse operator core for production multi-term Boltzmann runs."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import sparse
from scipy.sparse import linalg as spla

from electron_swarm.core.cross_sections import ProcessType
from electron_swarm.solvers.boltzmann_two_term import BoltzmannTwoTermSolver

from .angular import mu_coupling_coefficients
from .models import MultiTermCase
from .projection import ProjectedCollisionData, project_collision_data


@dataclass(frozen=True, slots=True)
class LegendreBlockLayout:
    """Index layout for coefficients ordered by Legendre term, then energy."""

    lmax: int
    n_energy_cells: int

    def __post_init__(self) -> None:
        if self.lmax < 0:
            raise ValueError("lmax must be non-negative")
        if self.n_energy_cells <= 0:
            raise ValueError("n_energy_cells must be positive")

    @property
    def n_legendre_terms(self) -> int:
        return self.lmax + 1

    @property
    def n_unknowns(self) -> int:
        return self.n_legendre_terms * self.n_energy_cells

    def offset(self, ell: int) -> int:
        self._validate_ell(ell)
        return ell * self.n_energy_cells

    def slice_for_l(self, ell: int) -> slice:
        start = self.offset(ell)
        return slice(start, start + self.n_energy_cells)

    def empty_matrix(self, *, format: str = "lil") -> sparse.spmatrix:
        return sparse.lil_matrix(
            (self.n_unknowns, self.n_unknowns), dtype=float
        ).asformat(format)

    def _validate_ell(self, ell: int) -> None:
        if ell < 0 or ell > self.lmax:
            raise IndexError(f"Legendre term ell={ell} outside 0..{self.lmax}")


class _SparseBlockBuilder:
    """Small internal helper for sparse Legendre block assembly."""

    def __init__(self, layout: LegendreBlockLayout) -> None:
        self.layout = layout
        self._matrix = layout.empty_matrix(format="lil")

    def add_block(
        self,
        row_l: int,
        col_l: int,
        block: sparse.spmatrix | np.ndarray,
    ) -> None:
        row = self.layout.slice_for_l(row_l)
        col = self.layout.slice_for_l(col_l)
        block_csr = sparse.csr_matrix(block)
        expected = (self.layout.n_energy_cells, self.layout.n_energy_cells)
        if block_csr.shape != expected:
            raise ValueError(
                f"Block shape {block_csr.shape} does not match energy grid {expected}"
            )
        self._matrix[row, col] = block_csr

    def add_diagonal(self, ell: int, values: np.ndarray) -> None:
        values = np.asarray(values, dtype=float)
        if values.shape != (self.layout.n_energy_cells,):
            raise ValueError("Diagonal values must match the energy grid length")
        self.add_block(ell, ell, sparse.diags(values, format="csr"))

    def to_csr(self) -> sparse.csr_matrix:
        return self._matrix.tocsr()


@dataclass(frozen=True, slots=True)
class _CouplingBlock:
    row_l: int
    col_l: int
    coefficient: float


@dataclass(frozen=True, slots=True)
class DensityNormalizationConstraint:
    """Linear constraint enforcing integral f_0(epsilon) d epsilon = 1."""

    layout: LegendreBlockLayout
    weights: np.ndarray
    active_l: int = 0
    definition: str = "integral_f0_dE_equals_1"
    coefficient_order: str = "legendre_major_energy_minor"

    @property
    def nonzero_count(self) -> int:
        return int(np.count_nonzero(self.weights))

    def apply(self, coefficients_flat: np.ndarray) -> float:
        coeff = np.asarray(coefficients_flat, dtype=float)
        if coeff.shape != (self.layout.n_unknowns,):
            raise ValueError("Coefficient vector length does not match layout")
        return float(np.dot(self.weights, coeff))


@dataclass(frozen=True, slots=True)
class OperatorAssemblyDiagnostics:
    """Developer-facing checks for sparse operator assembly contracts.

    The ``scaffold_*`` metadata names are retained for CSV compatibility with
    earlier builds, but the values describe the executable operator base matrix
    rather than a separate non-executing scaffold object.
    """

    layout: LegendreBlockLayout
    density_constraint: DensityNormalizationConstraint
    scaffold_status: str
    scaffold_nnz: int
    field_coupling_topology_nnz: int
    n_field_coupling_blocks: int
    l0_block_matches_native: bool
    density_constraint_matches_grid: bool
    native_lmax1_reference_ready: bool
    nonphysical_field_scaling: bool
    coefficient_order: str = "legendre_major_energy_minor"

    def as_metadata(self, prefix: str = "operator_assembly") -> dict[str, object]:
        return {
            f"{prefix}_coefficient_order": self.coefficient_order,
            f"{prefix}_n_unknowns": self.layout.n_unknowns,
            f"{prefix}_n_legendre_terms": self.layout.n_legendre_terms,
            f"{prefix}_n_energy_cells": self.layout.n_energy_cells,
            f"{prefix}_normalization": self.density_constraint.definition,
            f"{prefix}_normalization_nonzeros": (
                self.density_constraint.nonzero_count
            ),
            f"{prefix}_scaffold_status": self.scaffold_status,
            f"{prefix}_scaffold_nnz": self.scaffold_nnz,
            f"{prefix}_field_coupling_topology_nnz": (
                self.field_coupling_topology_nnz
            ),
            f"{prefix}_field_coupling_blocks": self.n_field_coupling_blocks,
            f"{prefix}_l0_block_matches_native": self.l0_block_matches_native,
            f"{prefix}_density_constraint_matches_grid": (
                self.density_constraint_matches_grid
            ),
            f"{prefix}_native_lmax1_reference_ready": (
                self.native_lmax1_reference_ready
            ),
            f"{prefix}_nonphysical_field_scaling": self.nonphysical_field_scaling,
        }


@dataclass(frozen=True, slots=True)
class OperatorSystem:
    """Executable sparse multi-term operator system."""

    layout: LegendreBlockLayout
    matrix: sparse.csr_matrix
    normalization: DensityNormalizationConstraint
    collisions: ProjectedCollisionData
    diagnostics: OperatorAssemblyDiagnostics
    base_matrix: sparse.csr_matrix
    field_coupling_matrix: sparse.csr_matrix
    streaming_matrix: sparse.csr_matrix
    inelastic_sink_frequency_s_inv: np.ndarray
    status: str = "operator"


@dataclass(frozen=True, slots=True)
class OperatorSolveState:
    coefficients_flat: np.ndarray
    growth_frequency_s_inv: float
    residual_L1: float
    iterations: int
    converged: bool
    warnings: tuple[str, ...]


def build_density_normalization_constraint(
    layout: LegendreBlockLayout, widths_eV: np.ndarray
) -> DensityNormalizationConstraint:
    """Build the normalization row for the flattened coefficient vector."""

    widths = np.asarray(widths_eV, dtype=float)
    if widths.shape != (layout.n_energy_cells,):
        raise ValueError("Energy widths must match layout.n_energy_cells")
    if np.any(widths <= 0.0) or not np.all(np.isfinite(widths)):
        raise ValueError("Energy widths must be positive finite values")
    weights = np.zeros(layout.n_unknowns, dtype=float)
    weights[layout.slice_for_l(0)] = widths
    return DensityNormalizationConstraint(layout=layout, weights=weights)


def _legendre_coupling_blocks(layout: LegendreBlockLayout) -> tuple[_CouplingBlock, ...]:
    blocks: list[_CouplingBlock] = []
    for ell in range(layout.n_legendre_terms):
        for target_l, coefficient in mu_coupling_coefficients(ell):
            if target_l is None or target_l > layout.lmax or coefficient == 0.0:
                continue
            blocks.append(
                _CouplingBlock(
                    row_l=int(target_l),
                    col_l=int(ell),
                    coefficient=float(coefficient),
                )
            )
    return tuple(blocks)


def _upwind_energy_advection_matrix(
    widths_eV: np.ndarray, face_advection_eV_s: np.ndarray
) -> sparse.csr_matrix:
    """Conservative energy-cell advection with zero boundary flux."""

    widths = np.asarray(widths_eV, dtype=float)
    face = np.asarray(face_advection_eV_s, dtype=float)
    n = len(widths)
    if n < 2 or np.any(widths <= 0.0) or not np.all(np.isfinite(widths)):
        raise ValueError("Energy widths must be positive finite values")
    if face.shape != (n + 1,) or not np.all(np.isfinite(face)):
        raise ValueError("Face advection must have n_energy_cells + 1 values")
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    for i in range(n):
        inv_w = 1.0 / widths[i]
        left = face[i]
        right = face[i + 1]
        if i > 0 and left != 0.0:
            upwind_col = i - 1 if left >= 0.0 else i
            rows.append(i)
            cols.append(upwind_col)
            vals.append(-left * inv_w)
        if i < n - 1 and right != 0.0:
            upwind_col = i if right >= 0.0 else i + 1
            rows.append(i)
            cols.append(upwind_col)
            vals.append(right * inv_w)
    return sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))


def _field_coupling_matrix(
    case: MultiTermCase,
    layout: LegendreBlockLayout,
    blocks: tuple[_CouplingBlock, ...],
) -> sparse.csr_matrix:
    scale = float(case.config.multiterm_boltzmann.field_coupling_scale)
    if scale == 0.0 or case.electric_field_V_m == 0.0 or not blocks:
        return sparse.csr_matrix((layout.n_unknowns, layout.n_unknowns))
    cell_advection = (
        -scale * float(case.electric_field_V_m) * np.maximum(case.grid.speeds_m_s, 0.0)
    )
    face_advection = np.zeros(layout.n_energy_cells + 1, dtype=float)
    face_advection[1:-1] = 0.5 * (cell_advection[:-1] + cell_advection[1:])
    energy_block = _upwind_energy_advection_matrix(
        case.grid.widths_eV, face_advection
    )
    builder = _SparseBlockBuilder(layout)
    for block in blocks:
        builder.add_block(block.row_l, block.col_l, block.coefficient * energy_block)
    return builder.to_csr()


def _operator_relaxation_frequencies(
    collisions: ProjectedCollisionData,
) -> tuple[np.ndarray, np.ndarray]:
    """Return angular momentum relaxation and isotropic inelastic loss.

    Integral cross sections do not define full differential anisotropy. The
    operator therefore damps higher Legendre terms with the available momentum
    relaxation, while inelastic/source processes remove anisotropic population
    once as a sink. Post-collision source remains in the isotropic l=0 block.
    """

    momentum = np.zeros_like(collisions.momentum_frequency_s_inv)
    inelastic_sink = np.zeros_like(collisions.momentum_frequency_s_inv)
    for item in collisions.processes:
        ptype = item.process.process_type
        if ptype in {ProcessType.MOMENTUM, ProcessType.EFFECTIVE, ProcessType.ELASTIC}:
            if item.contributes_to_momentum:
                momentum += item.frequency_s_inv
        elif ptype in {
            ProcessType.EXCITATION,
            ProcessType.IONIZATION,
            ProcessType.ATTACHMENT,
            ProcessType.SUPERELASTIC,
        }:
            inelastic_sink += item.frequency_s_inv
    return np.maximum(momentum, 1.0), inelastic_sink


def _diagnostics_from_blocks(
    case: MultiTermCase,
    layout: LegendreBlockLayout,
    normalization: DensityNormalizationConstraint,
    base_matrix: sparse.csr_matrix,
    native_l0_matrix: sparse.csr_matrix,
    n_coupling_blocks: int,
) -> OperatorAssemblyDiagnostics:
    l0_slice = layout.slice_for_l(0)
    l0_block = base_matrix[l0_slice, l0_slice]
    delta = l0_block - native_l0_matrix
    if delta.nnz == 0 or delta.data.size == 0:
        l0_matches = True
    else:
        l0_matches = bool(np.max(np.abs(delta.data)) < 1.0e-30)
    density_matches = bool(
        np.allclose(
            normalization.weights[l0_slice],
            case.grid.widths_eV,
            rtol=0.0,
            atol=0.0,
        )
        and (
            layout.n_legendre_terms == 1
            or np.count_nonzero(normalization.weights[layout.slice_for_l(1)]) == 0
        )
    )
    return OperatorAssemblyDiagnostics(
        layout=layout,
        density_constraint=normalization,
        scaffold_status="operator",
        scaffold_nnz=int(base_matrix.nnz),
        field_coupling_topology_nnz=int(n_coupling_blocks * layout.n_energy_cells),
        n_field_coupling_blocks=int(n_coupling_blocks),
        l0_block_matches_native=l0_matches,
        density_constraint_matches_grid=density_matches,
        native_lmax1_reference_ready=bool(l0_matches and density_matches),
        nonphysical_field_scaling=not np.isclose(
            float(case.config.multiterm_boltzmann.field_coupling_scale),
            1.0,
            rtol=0.0,
            atol=1.0e-12,
        ),
    )


def _streaming_matrix(
    case: MultiTermCase,
    layout: LegendreBlockLayout,
    blocks: tuple[_CouplingBlock, ...],
) -> sparse.csr_matrix:
    """Build the z-streaming operator v mu for finite-k hydrodynamics."""

    builder = _SparseBlockBuilder(layout)
    speed_block = sparse.diags(case.grid.speeds_m_s, format="csr")
    for block in blocks:
        builder.add_block(block.row_l, block.col_l, block.coefficient * speed_block)
    return builder.to_csr()


def assemble_operator_system(
    case: MultiTermCase,
    collisions: ProjectedCollisionData | None = None,
) -> OperatorSystem:
    """Assemble the executable multi-term operator system."""

    collisions = collisions or project_collision_data(case)
    layout = LegendreBlockLayout(
        lmax=int(case.config.multiterm_boltzmann.lmax),
        n_energy_cells=case.grid.n_cells,
    )
    normalization = build_density_normalization_constraint(layout, case.grid.widths_eV)
    native_e0 = BoltzmannTwoTermSolver(
        case.config, case.cross_sections
    ).assemble_native_operator_block(
        0.0,
        case.grid.centers_eV,
        case.grid.edges_eV,
        case.grid.widths_eV,
    )
    builder = _SparseBlockBuilder(layout)
    builder.add_block(0, 0, native_e0.matrix)
    momentum_relaxation, inelastic_sink = _operator_relaxation_frequencies(collisions)
    for ell in range(1, layout.n_legendre_terms):
        relaxation = float(ell) * momentum_relaxation + inelastic_sink
        builder.add_diagonal(ell, -relaxation)
    base_matrix = builder.to_csr()
    coupling_blocks = _legendre_coupling_blocks(layout)
    field_matrix = _field_coupling_matrix(case, layout, coupling_blocks)
    streaming_matrix = _streaming_matrix(case, layout, coupling_blocks)
    diagnostics = _diagnostics_from_blocks(
        case,
        layout,
        normalization,
        base_matrix,
        native_e0.matrix,
        len(coupling_blocks),
    )
    return OperatorSystem(
        layout=layout,
        matrix=(base_matrix + field_matrix).tocsr(),
        normalization=normalization,
        collisions=collisions,
        diagnostics=diagnostics,
        base_matrix=base_matrix,
        field_coupling_matrix=field_matrix,
        streaming_matrix=streaming_matrix,
        inelastic_sink_frequency_s_inv=inelastic_sink,
    )


def build_operator_assembly_diagnostics(
    case: MultiTermCase,
) -> OperatorAssemblyDiagnostics:
    """Return assembly diagnostics for the executable operator."""

    return assemble_operator_system(case).diagnostics


def solve_operator_system(
    case: MultiTermCase,
    system: OperatorSystem,
) -> OperatorSolveState:
    """Solve the multi-term operator dominant growth mode."""

    cfg = case.config.boltzmann_two_term.convergence
    layout = system.layout
    widths = case.grid.widths_eV
    weights_all = np.tile(widths, layout.n_legendre_terms)
    warnings: list[str] = []
    n = layout.n_unknowns
    k_eigs = min(6, n - 2)
    try:
        values, vectors = spla.eigs(
            system.matrix,
            k=k_eigs,
            sigma=0.0,
            which="LM",
            tol=cfg.residual_tolerance,
            maxiter=max(cfg.max_iterations * n, 1000),
        )
    except Exception:
        try:
            values, vectors = spla.eigs(
                system.matrix,
                k=k_eigs,
                which="LR",
                tol=cfg.residual_tolerance,
                maxiter=max(cfg.max_iterations * n, 1000),
            )
        except Exception:
            if n > 600:
                raise
            values, vectors = np.linalg.eig(system.matrix.toarray())
    order = np.argsort(np.abs(values))
    best_coeff: np.ndarray | None = None
    best_growth = 0.0
    best_residual = np.inf
    for idx in order:
        value = values[idx]
        vector = np.asarray(vectors[:, idx])
        if np.max(np.abs(np.imag(vector))) > 1.0e-7 * max(
            np.max(np.abs(np.real(vector))), 1.0
        ):
            continue
        coeff = np.real(vector).astype(float)
        norm = system.normalization.apply(coeff)
        if abs(norm) < 1.0e-200:
            continue
        coeff = coeff / norm
        f0 = coeff[layout.slice_for_l(0)]
        if float(np.sum(f0 * widths)) < 0.0:
            coeff = -coeff
            f0 = -f0
        max_f0 = max(float(np.max(np.abs(f0))), 1.0e-300)
        negative_mass = float(np.sum(np.clip(-f0, 0.0, None) * widths))
        if negative_mass > 1.0e-5 or float(np.min(f0)) < -1.0e-4 * max_f0:
            continue
        if np.any(f0 < 0.0):
            warnings.append("operator_f0_has_small_negative_roundoff")
        growth = float(np.real(value))
        lp = system.matrix @ coeff
        eig_res = lp - growth * coeff
        residual = float(
            np.sum(np.abs(eig_res) * weights_all)
            / max(np.sum(np.abs(lp) * weights_all), abs(growth), 1.0)
        )
        if residual < best_residual:
            best_coeff = coeff
            best_growth = growth
            best_residual = residual
    if best_coeff is None:
        raise RuntimeError("operator solve found no physical f0 mode")
    if not np.isfinite(best_residual) or best_residual > max(
        2.0e-5, 200.0 * cfg.residual_tolerance
    ):
        raise RuntimeError(
            "operator eigen residual is too large "
            f"({best_residual:.3e})"
        )
    return OperatorSolveState(
        coefficients_flat=best_coeff,
        growth_frequency_s_inv=float(best_growth),
        residual_L1=float(best_residual),
        iterations=1,
        converged=True,
        warnings=tuple(dict.fromkeys(warnings)),
    )
