"""Stationary collision-coupled response operator for homogeneous DC swarms."""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import SuperLU, splu
from threadpoolctl import threadpool_limits

from electron_swarm.solvers.propagator.collisions import (
    build_collision_operator,
    collision_generator_action,
    nonlocal_collision_inflow,
    reaction_number_rate_s_inv,
)
from electron_swarm.solvers.propagator.models import (
    CollisionOperatorData,
    PropagatorGrid,
)
from electron_swarm.solvers.propagator.origin_response import (
    AffineOriginResponse,
    OriginOperator,
    build_origin_operator,
    build_origin_response,
)
from electron_swarm.solvers.propagator.shell_response import (
    AffineShellResponse,
    build_shell_response,
)


# Shell responses consist of very many dense systems no larger than the polar
# grid.  Letting a BLAS runtime fan each of those systems across the whole host
# is slower and makes Arnoldi roundoff depend on unrelated process settings.
# Keep the scope local to response assembly; the global sparse solve and
# eigensolver retain the caller's normal execution policy.
RESPONSE_BUILD_BLAS_THREADS = 1


@dataclass(frozen=True, slots=True)
class ResponseApplication:
    population: np.ndarray
    outer_escape_s_inv: float
    boundary_residual_L1: float


@dataclass(frozen=True, slots=True)
class GrowthBounds:
    lower_s_inv: float
    upper_s_inv: float
    shift_s_inv: float


@dataclass(slots=True)
class StationaryResponseSystem:
    operator: "PropagatorOperator"
    acceleration_m_s2: float
    growth_frequency_s_inv: float
    shift_s_inv: float
    origin: AffineOriginResponse
    shells: tuple[AffineShellResponse, ...]
    boundary_factor: SuperLU
    build_time_s: float
    response_applications: int = 0

    @property
    def memory_bytes(self) -> int:
        factor = self.boundary_factor
        factor_bytes = int(
            factor.L.data.nbytes
            + factor.L.indices.nbytes
            + factor.L.indptr.nbytes
            + factor.U.data.nbytes
            + factor.U.indices.nbytes
            + factor.U.indptr.nbytes
        )
        return int(
            self.origin.memory_bytes
            + sum(item.memory_bytes for item in self.shells)
            + factor_bytes
        )

    @property
    def maximum_shell_rational_error_estimate(self) -> float:
        return max(
            (item.rational_error_estimate for item in self.shells),
            default=0.0,
        )

    @property
    def maximum_shell_coefficient_error_estimate(self) -> float:
        return max(
            (item.coefficient_error_estimate for item in self.shells),
            default=0.0,
        )

    @property
    def maximum_shell_coefficient_segments(self) -> int:
        return max(
            (item.coefficient_segments for item in self.shells),
            default=0,
        )

    def apply(self, population: np.ndarray) -> ResponseApplication:
        values = np.asarray(population, dtype=float)
        grid = self.operator.grid
        if values.shape != (grid.energy_cells, grid.polar_cells):
            raise ValueError("population shape does not match response system")
        self.response_applications += 1
        source = nonlocal_collision_inflow(
            self.operator.collisions,
            values,
        ) + self.shift_s_inv * values
        regular_count = len(self.shells)
        polar_cells = grid.polar_cells
        half = polar_cells // 2
        source_outflow = [
            shell.source_to_outflow @ source[index + 1]
            for index, shell in enumerate(self.shells)
        ]
        right_hand = np.zeros(regular_count * polar_cells, dtype=float)
        right_hand[:half] = (
            self.origin.incoming_to_outflow @ source_outflow[0][half:]
            + self.origin.source_to_outflow @ source[0]
        )
        for regular in range(1, regular_count):
            row = regular * polar_cells
            right_hand[row : row + half] = source_outflow[regular - 1][:half]
        for regular in range(regular_count - 1):
            row = regular * polar_cells + half
            right_hand[row : row + half] = source_outflow[regular + 1][half:]
        incoming_flat = self.boundary_factor.solve(right_hand)
        incoming = incoming_flat.reshape(regular_count, polar_cells)
        outgoing = np.empty_like(incoming)
        result = np.empty_like(values)
        for regular, shell in enumerate(self.shells):
            outgoing[regular] = (
                shell.scattering @ incoming[regular]
                + source_outflow[regular]
            )
            result[regular + 1] = (
                shell.incoming_to_population @ incoming[regular]
                + shell.source_to_population @ source[regular + 1]
            )
        origin_incoming = outgoing[0, half:]
        result[0] = (
            self.origin.incoming_to_population @ origin_incoming
            + self.origin.source_to_population @ source[0]
        )
        outer_escape = float(np.sum(outgoing[-1, :half]))
        origin_connection = (
            incoming[0, :half]
            - self.origin.incoming_to_outflow @ origin_incoming
            - self.origin.source_to_outflow @ source[0]
        )
        internal_connection = 0.0
        if regular_count > 1:
            internal_connection = max(
                float(
                    np.max(
                        np.abs(incoming[1:, :half] - outgoing[:-1, :half])
                    )
                ),
                float(
                    np.max(
                        np.abs(incoming[:-1, half:] - outgoing[1:, half:])
                    )
                ),
            )
        boundary_scale = max(
            float(np.sum(np.abs(incoming))),
            float(np.sum(np.abs(outgoing))),
            float(np.sum(np.abs(source))),
            1.0e-300,
        )
        boundary_residual = (
            float(np.sum(np.abs(origin_connection))) + internal_connection
        ) / boundary_scale
        return ResponseApplication(result, outer_escape, boundary_residual)


@dataclass(frozen=True, slots=True)
class PropagatorOperator:
    grid: PropagatorGrid
    collisions: CollisionOperatorData
    origin: OriginOperator

    @classmethod
    def build(
        cls,
        grid: PropagatorGrid,
        collisions: CollisionOperatorData,
    ) -> "PropagatorOperator":
        return cls(
            grid=grid,
            collisions=collisions,
            origin=build_origin_operator(grid, collisions),
        )

    @property
    def memory_bytes(self) -> int:
        return int(
            self.grid.cell_volumes_v3.nbytes
            + self.collisions.memory_bytes
            + self.origin.memory_bytes
        )

    def collision_generator_action(self, population: np.ndarray) -> np.ndarray:
        return collision_generator_action(self.collisions, population)

    def reaction_number_rate_s_inv(self, population: np.ndarray) -> float:
        return reaction_number_rate_s_inv(self.collisions, population)

    def growth_bounds(self, acceleration_m_s2: float) -> GrowthBounds:
        reaction_columns = np.zeros(self.grid.energy_cells, dtype=float)
        for transfer in self.collisions.inelastic:
            reaction_columns += (
                transfer.number_change_per_event * transfer.frequency_s_inv
            )
        outer_width = float(
            self.grid.speed_edges_m_s[-1] - self.grid.speed_edges_m_s[-2]
        )
        maximum_escape = (
            float(acceleration_m_s2)
            * float(np.max(np.maximum(self.grid.mu_centers, 0.0)))
            / outer_width
        )
        lower = float(np.min(reaction_columns) - maximum_escape)
        upper = float(np.max(reaction_columns))
        return GrowthBounds(lower, upper, 1.0)

    def minimum_response_shift(self, growth_frequency_s_inv: float) -> float:
        """Return only the shift required to keep local loss positive.

        A large global safety shift is algebraically redundant in the
        continuum but creates a P0 source/reconstruction error.  Choosing the
        local minimum keeps the positive inverse without changing the core
        transport scale.
        """

        growth = float(growth_frequency_s_inv)
        minimum_physical_loss = float(
            np.min(self.collisions.nonlocal_outflow_s_inv)
        )
        scale = max(
            abs(growth),
            float(np.max(self.collisions.nonlocal_outflow_s_inv)),
            1.0,
        )
        return max(
            1.0,
            -growth - minimum_physical_loss + 1.0e-12 * scale,
        )

    def build_response_system(
        self,
        acceleration_m_s2: float,
        growth_frequency_s_inv: float,
        shift_s_inv: float,
    ) -> StationaryResponseSystem:
        started = perf_counter()
        growth = float(growth_frequency_s_inv)
        shift = float(shift_s_inv)
        scalar = self.collisions.nonlocal_outflow_s_inv + growth + shift
        if np.any(~np.isfinite(scalar)) or np.any(scalar <= 0.0):
            raise ValueError(
                "response loss must remain positive over the growth bracket"
            )
        with threadpool_limits(
            limits=RESPONSE_BUILD_BLAS_THREADS,
            user_api="blas",
        ):
            origin = build_origin_response(
                self.origin,
                float(self.grid.speed_edges_m_s[1]),
                acceleration_m_s2,
                float(scalar[0]),
            )
            shells = tuple(
                build_shell_response(
                    self.grid,
                    self.collisions,
                    energy_index,
                    acceleration_m_s2,
                    float(scalar[energy_index]),
                )
                for energy_index in range(1, self.grid.energy_cells)
            )
            boundary = _build_boundary_matrix(origin, shells)
            factor = splu(boundary.tocsc())
        return StationaryResponseSystem(
            operator=self,
            acceleration_m_s2=float(acceleration_m_s2),
            growth_frequency_s_inv=growth,
            shift_s_inv=shift,
            origin=origin,
            shells=shells,
            boundary_factor=factor,
            build_time_s=perf_counter() - started,
        )

    def explicit_response_matrix(
        self,
        system: StationaryResponseSystem,
        *,
        maximum_cells: int = 512,
    ) -> sparse.csr_matrix:
        if self.grid.cells > maximum_cells:
            raise MemoryError(
                "explicit propagator response matrix is limited to "
                f"{maximum_cells} cells"
            )
        columns: list[np.ndarray] = []
        shape = (self.grid.energy_cells, self.grid.polar_cells)
        for index in range(self.grid.cells):
            basis = np.zeros(self.grid.cells, dtype=float)
            basis[index] = 1.0
            columns.append(system.apply(basis.reshape(shape)).population.reshape(-1))
        return sparse.csc_matrix(np.column_stack(columns)).tocsr()


def _append_dense_block(
    rows: list[np.ndarray],
    columns: list[np.ndarray],
    entries: list[np.ndarray],
    row_start: int,
    column_start: int,
    block: np.ndarray,
) -> None:
    values = np.asarray(block, dtype=float)
    row_count, column_count = values.shape
    rows.append(
        np.repeat(np.arange(row_start, row_start + row_count), column_count)
    )
    columns.append(
        np.tile(np.arange(column_start, column_start + column_count), row_count)
    )
    entries.append(values.reshape(-1))


def _build_boundary_matrix(
    origin: AffineOriginResponse,
    shells: tuple[AffineShellResponse, ...],
) -> sparse.csr_matrix:
    if not shells:
        raise ValueError("response operator requires at least one regular shell")
    polar_cells = shells[0].scattering.shape[0]
    half = polar_cells // 2
    regular_count = len(shells)
    size = regular_count * polar_cells
    rows: list[np.ndarray] = [np.arange(size)]
    columns: list[np.ndarray] = [np.arange(size)]
    entries: list[np.ndarray] = [np.ones(size, dtype=float)]
    _append_dense_block(
        rows,
        columns,
        entries,
        0,
        0,
        -(origin.incoming_to_outflow @ shells[0].scattering[half:, :]),
    )
    for regular in range(1, regular_count):
        _append_dense_block(
            rows,
            columns,
            entries,
            regular * polar_cells,
            (regular - 1) * polar_cells,
            -shells[regular - 1].scattering[:half, :],
        )
    for regular in range(regular_count - 1):
        _append_dense_block(
            rows,
            columns,
            entries,
            regular * polar_cells + half,
            (regular + 1) * polar_cells,
            -shells[regular + 1].scattering[half:, :],
        )
    matrix = sparse.coo_matrix(
        (
            np.concatenate(entries),
            (np.concatenate(rows), np.concatenate(columns)),
        ),
        shape=(size, size),
        dtype=float,
    ).tocsr()
    matrix.sum_duplicates()
    return matrix


__all__ = [
    "GrowthBounds",
    "PropagatorOperator",
    "ResponseApplication",
    "StationaryResponseSystem",
    "build_collision_operator",
]
