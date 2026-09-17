"""Characteristic-aligned response of the first velocity-space ball.

The origin is represented in ``u=v_parallel/R`` and
``rho=v_perpendicular**2/R**2``.  Acceleration is a one-way Cartesian flux in
``u`` and the axisymmetric measure is ``pi du d rho``.  No artificial inner
radial boundary or angular interpolation is introduced.
"""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu

from electron_swarm.physics.electron_neutral import speed_from_energy_m_s
from electron_swarm.solvers.propagator.angular_kernel import (
    maxent_p1_transition_matrices,
)
from electron_swarm.solvers.propagator.models import (
    CollisionOperatorData,
    ElasticTransfer,
    PropagatorGrid,
)


_GAUSS_NODES, _GAUSS_WEIGHTS = np.polynomial.legendre.leggauss(8)


def _radial_band_average_rates(
    transfer: ElasticTransfer,
    outer_energy_eV: float,
    radial_bands: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Volume-average total and momentum rates inside the origin ball."""

    edges = np.linspace(0.0, 1.0, int(radial_bands) + 1)
    collision = np.zeros(radial_bands, dtype=float)
    momentum = np.zeros(radial_bands, dtype=float)
    for index, (lower, upper) in enumerate(
        zip(edges[:-1], edges[1:], strict=True)
    ):
        radius = (
            0.5 * (upper - lower) * _GAUSS_NODES
            + 0.5 * (upper + lower)
        )
        energy = float(outer_energy_eV) * radius * radius
        speed = speed_from_energy_m_s(energy)
        total_sigma = np.zeros_like(radius)
        momentum_sigma = np.zeros_like(radius)
        for process in transfer.collision_processes:
            total_sigma += process.sigma(energy)
        for process in transfer.momentum_processes:
            momentum_sigma += process.sigma(energy)
        normalization = (upper**3 - lower**3) / 3.0
        measure = 0.5 * (upper - lower) * _GAUSS_WEIGHTS * radius * radius
        collision[index] = (
            transfer.density_scale_m3
            * float(np.dot(measure, total_sigma * speed))
            / normalization
        )
        momentum[index] = (
            transfer.density_scale_m3
            * float(np.dot(measure, momentum_sigma * speed))
            / normalization
        )
    return collision, momentum


def _positive_part_integral(
    lower_u: float,
    upper_u: float,
    lower_rho: float,
    upper_rho: float,
    lower_radius_squared: float,
    upper_radius_squared: float,
) -> float:
    points = [lower_u, upper_u]
    for squared in (
        lower_radius_squared - lower_rho,
        lower_radius_squared - upper_rho,
        upper_radius_squared - lower_rho,
        upper_radius_squared - upper_rho,
        1.0 - lower_rho,
        1.0 - upper_rho,
    ):
        if squared <= 0.0:
            continue
        root = float(np.sqrt(squared))
        if lower_u < root < upper_u:
            points.append(root)
        if lower_u < -root < upper_u:
            points.append(-root)
    points = sorted(set(points))
    integral = 0.0
    for left, right in zip(points[:-1], points[1:], strict=True):
        if right <= left:
            continue
        nodes = 0.5 * (right - left) * _GAUSS_NODES + 0.5 * (right + left)
        lower = np.maximum(lower_rho, lower_radius_squared - nodes * nodes)
        upper = np.minimum.reduce(
            (
                np.full_like(nodes, upper_rho),
                upper_radius_squared - nodes * nodes,
                1.0 - nodes * nodes,
            )
        )
        integral += 0.5 * (right - left) * float(
            np.dot(_GAUSS_WEIGHTS, np.maximum(upper - lower, 0.0))
        )
    return integral


def _boundary_rho_measure(
    lower_u: float,
    upper_u: float,
    lower_rho: float,
    upper_rho: float,
    *,
    incoming: bool,
) -> float:
    if incoming:
        left = max(lower_u, -1.0)
        right = min(upper_u, 0.0)
        if right <= left:
            return 0.0
        mapped_lower = 1.0 - left * left
        mapped_upper = 1.0 - right * right
    else:
        left = max(lower_u, 0.0)
        right = min(upper_u, 1.0)
        if right <= left:
            return 0.0
        mapped_lower = 1.0 - right * right
        mapped_upper = 1.0 - left * left
    return max(
        0.0,
        min(upper_rho, mapped_upper) - max(lower_rho, mapped_lower),
    )


def _vertical_rho_measure(
    u_face: float,
    lower_rho: float,
    upper_rho: float,
) -> float:
    if abs(u_face) >= 1.0:
        return 0.0
    return max(0.0, min(upper_rho, 1.0 - u_face * u_face) - lower_rho)


@dataclass(frozen=True, slots=True)
class OriginGeometry:
    volumes: np.ndarray
    streaming_generator: sparse.csr_matrix
    boundary_inflow: np.ndarray
    boundary_outflow: sparse.csr_matrix
    collision_overlap: sparse.csr_matrix
    collision_cell_volumes: np.ndarray
    source_basis: sparse.csr_matrix
    population_aggregation: sparse.csr_matrix
    radial_bands: int
    polar_cells: int
    build_time_s: float
    partition_error: float

    @property
    def memory_bytes(self) -> int:
        matrices = (
            self.streaming_generator,
            self.boundary_outflow,
            self.collision_overlap,
            self.source_basis,
            self.population_aggregation,
        )
        return int(
            self.volumes.nbytes
            + self.boundary_inflow.nbytes
            + self.collision_cell_volumes.nbytes
            + sum(
                item.data.nbytes + item.indices.nbytes + item.indptr.nbytes
                for item in matrices
            )
        )


@dataclass(frozen=True, slots=True)
class OriginOperator:
    geometry: OriginGeometry
    collision_generator_s_inv: sparse.csr_matrix
    detailed_balance_error: float
    column_balance_error: float

    @property
    def memory_bytes(self) -> int:
        matrix = self.collision_generator_s_inv
        return int(
            self.geometry.memory_bytes
            + matrix.data.nbytes
            + matrix.indices.nbytes
            + matrix.indptr.nbytes
        )


@dataclass(frozen=True, slots=True)
class AffineOriginResponse:
    incoming_to_outflow: np.ndarray
    source_to_outflow: np.ndarray
    incoming_to_population: np.ndarray
    source_to_population: np.ndarray
    minimum_entry: float
    factor_time_s: float

    @property
    def memory_bytes(self) -> int:
        return int(
            self.incoming_to_outflow.nbytes
            + self.source_to_outflow.nbytes
            + self.incoming_to_population.nbytes
            + self.source_to_population.nbytes
        )


def _collision_partition(
    cells: list[tuple[int, int]],
    u_edges: np.ndarray,
    rho_edges: np.ndarray,
    radial_bands: int,
    theta_edges: np.ndarray,
    volumes: np.ndarray,
) -> tuple[sparse.csr_matrix, np.ndarray, float]:
    """Integrate cut cells against radial/angular collision cells."""

    radial_edges = np.linspace(0.0, 1.0, radial_bands + 1)
    mu_faces = np.cos(theta_edges)
    polar_cells = len(theta_edges) - 1
    rows: list[int] = []
    columns: list[int] = []
    entries: list[float] = []
    maximum_partition_error = 0.0
    quadrature_nodes, quadrature_weights = np.polynomial.legendre.leggauss(8)

    for state, (u_index, rho_index) in enumerate(cells):
        lower_u = float(u_edges[u_index])
        upper_u = float(u_edges[u_index + 1])
        lower_rho = float(rho_edges[rho_index])
        upper_rho = float(rho_edges[rho_index + 1])
        u_breakpoints = [lower_u, upper_u]
        for radius_face in radial_edges:
            radius_squared = float(radius_face * radius_face)
            for rho_face in (lower_rho, upper_rho):
                squared = radius_squared - rho_face
                if squared <= 0.0:
                    continue
                root = float(np.sqrt(squared))
                if lower_u < root < upper_u:
                    u_breakpoints.append(root)
                if lower_u < -root < upper_u:
                    u_breakpoints.append(-root)
        for face in np.abs(mu_faces[1:-1]):
            if face <= 1.0e-15 or face >= 1.0 - 1.0e-15:
                continue
            coefficient = 1.0 / float(face * face) - 1.0
            for rho_face in (lower_rho, upper_rho):
                if rho_face <= 0.0:
                    continue
                root = float(np.sqrt(rho_face / coefficient))
                if lower_u < root < upper_u:
                    u_breakpoints.append(root)
                if lower_u < -root < upper_u:
                    u_breakpoints.append(-root)
        nodes_parts: list[np.ndarray] = []
        weights_parts: list[np.ndarray] = []
        for u_left, u_right in zip(
            sorted(set(u_breakpoints))[:-1],
            sorted(set(u_breakpoints))[1:],
            strict=True,
        ):
            if u_right <= u_left:
                continue
            nodes_parts.append(
                0.5 * (u_right - u_left) * quadrature_nodes
                + 0.5 * (u_right + u_left)
            )
            weights_parts.append(
                0.5 * (u_right - u_left) * quadrature_weights
            )
        nodes = np.concatenate(nodes_parts)
        weights = np.concatenate(weights_parts)
        local: dict[int, float] = {}
        for u, weight in zip(nodes, weights, strict=True):
            rho_low = lower_rho
            rho_high = min(upper_rho, 1.0 - float(u) ** 2)
            if rho_high <= rho_low:
                continue
            breakpoints = [rho_low, rho_high]
            u2 = float(u) ** 2
            for radial_face in radial_edges[1:-1]:
                value = float(radial_face * radial_face - u2)
                if rho_low < value < rho_high:
                    breakpoints.append(value)
            absolute_u = abs(float(u))
            if absolute_u > 0.0:
                for face in np.abs(mu_faces[1:-1]):
                    if face <= 1.0e-15:
                        continue
                    value = u2 * (1.0 / float(face * face) - 1.0)
                    if rho_low < value < rho_high:
                        breakpoints.append(value)
            breakpoints = sorted(set(breakpoints))
            for rho_left, rho_right in zip(
                breakpoints[:-1], breakpoints[1:], strict=True
            ):
                if rho_right <= rho_left:
                    continue
                rho_mid = 0.5 * (rho_left + rho_right)
                radius = np.sqrt(u2 + rho_mid)
                if radius <= 0.0:
                    continue
                radial = min(
                    int(np.searchsorted(radial_edges, radius, side="right") - 1),
                    radial_bands - 1,
                )
                mu = float(u) / radius
                theta = float(np.arccos(np.clip(mu, -1.0, 1.0)))
                polar = min(
                    int(np.searchsorted(theta_edges, theta, side="right") - 1),
                    polar_cells - 1,
                )
                collision_cell = radial * polar_cells + polar
                local[collision_cell] = local.get(collision_cell, 0.0) + (
                    np.pi * float(weight) * (rho_right - rho_left)
                )
        total = float(sum(local.values()))
        if total <= 0.0:
            raise FloatingPointError("origin collision partition lost a cut cell")
        maximum_partition_error = max(
            maximum_partition_error,
            abs(total - volumes[state]) / volumes[state],
        )
        scale = volumes[state] / total
        for collision_cell, value in local.items():
            rows.append(collision_cell)
            columns.append(state)
            entries.append(value * scale)

    overlap = sparse.coo_matrix(
        (entries, (rows, columns)),
        shape=(radial_bands * polar_cells, len(cells)),
        dtype=float,
    ).tocsr()
    overlap.sum_duplicates()
    angular_weights = 0.5 * (
        np.cos(theta_edges[:-1]) - np.cos(theta_edges[1:])
    )
    radial_volumes = (4.0 * np.pi / 3.0) * np.diff(radial_edges**3)
    desired_rows = np.kron(radial_volumes, angular_weights)
    desired_rows *= float(np.sum(volumes)) / float(np.sum(desired_rows))
    if np.any(desired_rows <= 0.0):
        raise FloatingPointError("origin collision partition has an empty cell")
    # The one-dimensional quadrature determines the sparse intersection
    # pattern.  Balance its two analytically known finite-volume margins so
    # restriction and prolongation are exact adjoints.  This is construction
    # of the common-refinement mass matrix, not a correction of a solution.
    base_overlap = overlap.tocsr()
    row_scale = np.ones(base_overlap.shape[0], dtype=float)
    column_scale = np.ones(base_overlap.shape[1], dtype=float)
    for _ in range(20000):
        row_denominator = np.asarray(
            base_overlap @ column_scale
        ).reshape(-1)
        if np.any(row_denominator <= 0.0):
            raise FloatingPointError(
                "origin common-refinement contains an empty row support"
            )
        row_scale = desired_rows / row_denominator
        column_denominator = np.asarray(
            base_overlap.T @ row_scale
        ).reshape(-1)
        if np.any(column_denominator <= 0.0):
            raise FloatingPointError(
                "origin common-refinement contains an empty column support"
            )
        column_scale = volumes / column_denominator
        balanced_rows = row_scale * np.asarray(
            base_overlap @ column_scale
        ).reshape(-1)
        balanced_columns = column_scale * np.asarray(
            base_overlap.T @ row_scale
        ).reshape(-1)
        row_error = float(
            np.max(np.abs(balanced_rows - desired_rows))
        ) / max(float(np.max(desired_rows)), 1.0e-300)
        column_error = float(
            np.max(np.abs(balanced_columns - volumes))
        ) / max(float(np.max(volumes)), 1.0e-300)
        if max(row_error, column_error) <= 2.0e-10:
            break
    else:
        raise FloatingPointError(
            "origin common-refinement margins did not converge: "
            f"row_error={row_error:.3e}, column_error={column_error:.3e}"
        )
    overlap = (
        sparse.diags(row_scale)
        @ base_overlap
        @ sparse.diags(column_scale)
    ).tocsr()
    collision_volumes = np.asarray(overlap.sum(axis=1)).reshape(-1)
    return overlap, collision_volumes, maximum_partition_error


def build_origin_geometry(
    grid: PropagatorGrid,
    *,
    u_cells: int | None = None,
    radial_bands: int | None = None,
) -> OriginGeometry:
    started = perf_counter()
    if grid.polar_cells % 2:
        raise ValueError("origin response requires an even polar grid")
    polar_cells = grid.polar_cells
    half = polar_cells // 2
    u_count = max(32, 4 * polar_cells) if u_cells is None else int(u_cells)
    band_count = max(16, u_count // 2) if radial_bands is None else int(radial_bands)
    if u_count < 8 or u_count % 2 or band_count < 4:
        raise ValueError("invalid origin cut-cell resolution")
    external_rho_edges = np.sin(grid.theta_edges_rad[: half + 1]) ** 2
    rho_edges = external_rho_edges.copy()
    u_edges = np.linspace(-1.0, 1.0, u_count + 1)
    cells: list[tuple[int, int]] = []
    lookup: dict[tuple[int, int], int] = {}
    volumes: list[float] = []
    external_bins: list[int] = []
    incoming_areas: list[float] = []
    outgoing_areas: list[float] = []
    for rho_index, (lower_rho, upper_rho) in enumerate(
        zip(rho_edges[:-1], rho_edges[1:], strict=True)
    ):
        for u_index, (lower_u, upper_u) in enumerate(
            zip(u_edges[:-1], u_edges[1:], strict=True)
        ):
            volume = np.pi * _positive_part_integral(
                float(lower_u),
                float(upper_u),
                float(lower_rho),
                float(upper_rho),
                0.0,
                1.0,
            )
            if volume <= 1.0e-14:
                continue
            lookup[(u_index, rho_index)] = len(cells)
            cells.append((u_index, rho_index))
            volumes.append(volume)
            external_bins.append(rho_index)
            incoming_areas.append(
                np.pi
                * _boundary_rho_measure(
                    float(lower_u),
                    float(upper_u),
                    float(lower_rho),
                    float(upper_rho),
                    incoming=True,
                )
            )
            outgoing_areas.append(
                np.pi
                * _boundary_rho_measure(
                    float(lower_u),
                    float(upper_u),
                    float(lower_rho),
                    float(upper_rho),
                    incoming=False,
                )
            )

    volume_array = np.asarray(volumes, dtype=float)
    incoming_area = np.asarray(incoming_areas, dtype=float)
    outgoing_area = np.asarray(outgoing_areas, dtype=float)
    external_bin = np.asarray(external_bins, dtype=int)
    states = len(cells)
    rows: list[int] = []
    columns: list[int] = []
    entries: list[float] = []
    streaming_outflow = np.zeros(states, dtype=float)
    for source, (u_index, rho_index) in enumerate(cells):
        upper_u = float(u_edges[u_index + 1])
        interface_area = np.pi * _vertical_rho_measure(
            upper_u,
            float(rho_edges[rho_index]),
            float(rho_edges[rho_index + 1]),
        )
        if interface_area > 1.0e-14:
            destination = lookup.get((u_index + 1, rho_index))
            if destination is None:
                raise FloatingPointError("origin characteristic lost an interface")
            coefficient = interface_area / volume_array[source]
            rows.append(destination)
            columns.append(source)
            entries.append(coefficient)
            streaming_outflow[source] += coefficient
        streaming_outflow[source] += outgoing_area[source] / volume_array[source]
    streaming_gain = sparse.coo_matrix(
        (entries, (rows, columns)), shape=(states, states), dtype=float
    ).tocsr()
    streaming = streaming_gain - sparse.diags(streaming_outflow)

    projected_area = np.pi * np.diff(external_rho_edges)
    boundary_inflow = np.zeros((states, half), dtype=float)
    for state in range(states):
        boundary_inflow[state, external_bin[state]] += incoming_area[state]
    boundary_inflow /= projected_area[None, :]
    out_rows: list[int] = []
    out_columns: list[int] = []
    out_entries: list[float] = []
    for state in range(states):
        if outgoing_area[state] <= 0.0:
            continue
        out_rows.append(int(external_bin[state]))
        out_columns.append(state)
        out_entries.append(float(outgoing_area[state] / volume_array[state]))
    boundary_outflow = sparse.coo_matrix(
        (out_entries, (out_rows, out_columns)),
        shape=(half, states),
        dtype=float,
    ).tocsr()

    overlap, collision_volumes, partition_error = _collision_partition(
        cells,
        u_edges,
        rho_edges,
        band_count,
        grid.theta_edges_rad,
        volume_array,
    )
    angular_overlap = sparse.vstack(
        [
            sparse.csr_matrix(
                overlap[index::polar_cells].sum(axis=0)
            )
            for index in range(polar_cells)
        ],
        format="csr",
    )
    angular_volumes = np.asarray(angular_overlap.sum(axis=1)).reshape(-1)
    source_basis = (
        angular_overlap.T @ sparse.diags(1.0 / angular_volumes)
    ).tocsr()
    population_aggregation = (
        angular_overlap @ sparse.diags(1.0 / volume_array)
    ).tocsr()
    return OriginGeometry(
        volumes=volume_array,
        streaming_generator=streaming.tocsr(),
        boundary_inflow=boundary_inflow,
        boundary_outflow=boundary_outflow,
        collision_overlap=overlap,
        collision_cell_volumes=collision_volumes,
        source_basis=source_basis,
        population_aggregation=population_aggregation,
        radial_bands=band_count,
        polar_cells=polar_cells,
        build_time_s=perf_counter() - started,
        partition_error=partition_error,
    )


def build_origin_operator(
    grid: PropagatorGrid,
    collisions: CollisionOperatorData,
    *,
    u_cells: int | None = None,
    radial_bands: int | None = None,
) -> OriginOperator:
    geometry = build_origin_geometry(
        grid,
        u_cells=u_cells,
        radial_bands=radial_bands,
    )
    overlap = geometry.collision_overlap
    volumes = geometry.volumes
    collision_volumes = geometry.collision_cell_volumes
    radial_bands = geometry.radial_bands
    polar_cells = geometry.polar_cells
    restriction = overlap @ sparse.diags(1.0 / volumes)
    prolongation = overlap.T @ sparse.diags(1.0 / collision_volumes)
    block = sparse.csr_matrix(
        (radial_bands * polar_cells, radial_bands * polar_cells), dtype=float
    )
    loss_rate = np.zeros(radial_bands * polar_cells, dtype=float)
    for transfer in collisions.elastic:
        collision_frequency, momentum_frequency = _radial_band_average_rates(
            transfer,
            float(grid.energy_edges_eV[1]),
            radial_bands,
        )
        if transfer.isotropic:
            relative_mismatch = np.abs(
                collision_frequency - momentum_frequency
            ) / np.maximum(
                np.maximum(collision_frequency, momentum_frequency),
                1.0e-300,
            )
            if float(np.max(relative_mismatch)) > 1.0e-8:
                raise ValueError(
                    f"{transfer.species}: origin isotropic closure has "
                    "inconsistent total and momentum-transfer rates"
                )
            collision_frequency = momentum_frequency
            band_kernels = None
        else:
            active = collision_frequency > 0.0
            mean_cosine = np.where(
                active,
                1.0
                - momentum_frequency
                / np.maximum(collision_frequency, 1.0e-300),
                0.0,
            )
            if np.any(
                active
                & (
                    (mean_cosine < -1.0 - 1.0e-10)
                    | (mean_cosine > 1.0 + 1.0e-10)
                )
            ):
                raise ValueError(
                    f"{transfer.species}: origin total and momentum rates "
                    "imply an angular moment outside [-1, 1]"
                )
            band_kernels = maxent_p1_transition_matrices(
                grid,
                np.clip(mean_cosine, -1.0, 1.0),
            )
        transition_blocks: list[sparse.csr_matrix] = []
        for radial in range(radial_bands):
            frequency = float(collision_frequency[radial])
            section = slice(
                radial * polar_cells,
                (radial + 1) * polar_cells,
            )
            local_weights = collision_volumes[section]
            local_weights = local_weights / float(np.sum(local_weights))
            if transfer.isotropic:
                transition = np.repeat(
                    local_weights[:, None], polar_cells, axis=1
                )
            else:
                assert band_kernels is not None
                base = band_kernels[radial]
                base_joint = (
                    base * collisions.isotropic_weights[None, :]
                )
                scaling = np.ones(polar_cells, dtype=float)
                for _ in range(100):
                    margin = scaling * (base_joint @ scaling)
                    ratio = local_weights / np.maximum(margin, 1.0e-300)
                    scaling *= np.sqrt(ratio)
                    if float(np.max(np.abs(ratio - 1.0))) <= 2.0e-13:
                        break
                joint = (
                    scaling[:, None] * base_joint * scaling[None, :]
                )
                transition = joint / local_weights[None, :]
            transition_blocks.append(
                frequency * sparse.csr_matrix(transition)
            )
        block += sparse.block_diag(transition_blocks, format="csr")
        loss_rate += np.repeat(collision_frequency, polar_cells)
    gain = (prolongation @ block @ restriction).tocsr()
    state_loss = np.asarray(restriction.T @ loss_rate).reshape(-1)
    generator = (gain - sparse.diags(state_loss)).tocsr()
    column_error = float(
        np.max(np.abs(np.asarray(generator.sum(axis=0)).reshape(-1)))
    )
    weighted = generator @ sparse.diags(volumes)
    defect = weighted - weighted.T
    detailed_error = float(np.max(np.abs(defect.data))) if defect.nnz else 0.0
    scale = max(float(np.max(np.abs(generator.data), initial=0.0)), 1.0)
    if column_error > 2.0e-11 * scale:
        raise FloatingPointError("origin collision generator is not conservative")
    if detailed_error > 2.0e-11 * scale:
        raise FloatingPointError("origin collision generator is not reciprocal")
    return OriginOperator(
        geometry=geometry,
        collision_generator_s_inv=generator,
        detailed_balance_error=detailed_error,
        column_balance_error=column_error,
    )


def build_origin_response(
    operator: OriginOperator,
    radius_m_s: float,
    acceleration_m_s2: float,
    local_scalar_loss_s_inv: float,
) -> AffineOriginResponse:
    geometry = operator.geometry
    radius = float(radius_m_s)
    acceleration = float(acceleration_m_s2)
    loss = float(local_scalar_loss_s_inv)
    if radius <= 0.0 or acceleration <= 0.0 or loss <= 0.0:
        raise ValueError("origin response requires positive radius, field, and loss")
    streaming_rate = acceleration / radius
    generator = (
        streaming_rate * geometry.streaming_generator
        + operator.collision_generator_s_inv
        - sparse.eye(len(geometry.volumes), format="csr") * loss
    ).tocsc()
    factor_started = perf_counter()
    factor = splu(-generator)
    factor_time = perf_counter() - factor_started
    right = np.column_stack(
        (
            geometry.boundary_inflow,
            geometry.source_basis.toarray(),
        )
    )
    interior = factor.solve(right)
    half = geometry.boundary_inflow.shape[1]
    boundary_state = interior[:, :half]
    source_state = interior[:, half:]
    outgoing = streaming_rate * geometry.boundary_outflow
    incoming_to_outflow = np.asarray(outgoing @ boundary_state)
    source_to_outflow = np.asarray(outgoing @ source_state)
    incoming_to_population = np.asarray(
        geometry.population_aggregation @ boundary_state
    )
    source_to_population = np.asarray(
        geometry.population_aggregation @ source_state
    )
    minimum = float(
        min(
            np.min(incoming_to_outflow),
            np.min(source_to_outflow),
            np.min(incoming_to_population),
            np.min(source_to_population),
        )
    )
    scale = max(
        float(np.max(incoming_to_outflow)),
        float(np.max(source_to_outflow)),
        float(np.max(incoming_to_population)),
        float(np.max(source_to_population)),
        1.0,
    )
    if minimum < -2.0e-11 * scale:
        raise FloatingPointError("origin response is not positive")
    return AffineOriginResponse(
        incoming_to_outflow=incoming_to_outflow,
        source_to_outflow=source_to_outflow,
        incoming_to_population=incoming_to_population,
        source_to_population=source_to_population,
        minimum_entry=minimum,
        factor_time_s=factor_time,
    )
