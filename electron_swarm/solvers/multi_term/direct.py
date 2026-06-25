"""Direct PN reduction path for the product multi-term solver."""

from __future__ import annotations

import numpy as np
from scipy import optimize, sparse
from scipy.sparse import linalg as spla

from electron_swarm.core.cross_sections import ProcessType
from electron_swarm.core.result_metadata import TRANSPORT_F0_GRADIENT_RECONSTRUCTION
from electron_swarm.core.transport import ElectronTransport
from electron_swarm.physics.angular_scattering import build_angular_model
from electron_swarm.solvers.kinetic import (
    EffectiveCollisionData,
    KineticOperatorBlock,
    assemble_native_operator_blocks,
    assemble_energy_flux_operator,
    compute_rates_from_eedf,
    make_two_term_energy_grid,
    mean_energy_from_eedf,
    negative_mass_fraction,
    normalize_eedf,
    transport_from_eedf,
    weighted_integral,
)

from .diagnostics import SolverDiagnostics
from .models import MultiTermCase, MultiTermSolution, RateSet


_NEGATIVE_MASS_LIMIT = 1.0e-8
_HIGHER_L_NEGATIVE_MASS_LIMIT = 1.0e-7


def build_higher_l_collision_damping(
    native: KineticOperatorBlock,
    moments: np.ndarray,
) -> np.ndarray:
    """Build validated angular-closure damping arrays for coefficient-space PN."""

    lmax = moments.shape[0] - 1
    energy_shape = native.grid.energy_eV.shape
    if lmax < 1 or moments.shape[1:] != energy_shape:
        raise ValueError("angular moments must have shape (lmax + 1, n_energy)")
    if not np.all(np.isfinite(moments)) or not np.all(np.abs(moments) <= 1.0):
        raise ValueError("angular moments must be finite and within [-1, 1]")
    if not np.allclose(moments[0], 1.0, rtol=1.0e-10, atol=1.0e-10):
        raise ValueError("angular moments must have m0=1")
    sigma_total = np.asarray(native.collisions.sigma_total_like, dtype=float)
    if sigma_total.shape != energy_shape or np.any(~np.isfinite(sigma_total)):
        raise ValueError("sigma_total_like is required for higher-l damping")
    if np.any(sigma_total <= 0.0):
        raise ValueError("sigma_total_like must be positive for higher-l damping")

    damping = np.zeros_like(moments, dtype=float)
    damping[1] = native.collisions.nu_m
    for ell in range(2, lmax + 1):
        damping[ell] = (
            native.gas_number_density_m3
            * native.grid.speed_m_s
            * sigma_total
            * (1.0 - moments[ell])
            + native.collisions.inelastic_loss_frequency_s_inv
        )
    if not np.all(np.isfinite(damping[1:])) or np.any(damping[1:] < 0.0):
        raise ValueError("higher-l collision damping must be finite and nonnegative")
    return damping


def _native_with_moment_table_l1_damping(
    native: KineticOperatorBlock,
    moments: np.ndarray,
) -> KineticOperatorBlock:
    if moments.shape[0] < 2 or moments.shape[1:] != native.grid.energy_eV.shape:
        raise ValueError("moment_table must provide m0 and m1 on the solver grid")
    if not np.allclose(moments[0], 1.0, rtol=1.0e-10, atol=1.0e-10):
        raise ValueError("moment_table m0 values must be 1")
    m1 = np.asarray(moments[1], dtype=float)
    sigma_total = np.asarray(native.collisions.sigma_total_like, dtype=float)
    sigma_m = sigma_total * (1.0 - m1)
    if (
        sigma_m.shape != native.grid.energy_eV.shape
        or not np.all(np.isfinite(sigma_m))
        or np.any(sigma_m <= 0.0)
    ):
        raise ValueError("moment_table m1 produced nonpositive momentum damping")
    nu_m = native.gas_number_density_m3 * native.grid.speed_m_s * sigma_m
    nu_m_over_N = native.grid.speed_m_s * sigma_m
    collisions = EffectiveCollisionData(
        nu_m=np.maximum(nu_m, 1.0e-60),
        nu_m_over_N=np.maximum(nu_m_over_N, 1.0e-80),
        sigma_m=sigma_m,
        sigma_total_like=sigma_total,
        inelastic_loss_frequency_s_inv=native.collisions.inelastic_loss_frequency_s_inv,
        elastic_A_eV_s=native.collisions.elastic_A_eV_s,
        elastic_D_eV2_s=native.collisions.elastic_D_eV2_s,
        processes=native.collisions.processes,
        effective_momentum_model="moment_table_l1",
        sigma_total_like_model=native.collisions.sigma_total_like_model,
        suppressed_momentum_processes=native.collisions.suppressed_momentum_processes,
    )
    energy_flux = assemble_energy_flux_operator(
        native.grid.energy_eV,
        native.grid.widths_eV,
        native.electric_field_V_m,
        collisions,
    )
    return KineticOperatorBlock(
        grid=native.grid,
        gas_number_density_m3=native.gas_number_density_m3,
        electric_field_V_m=native.electric_field_V_m,
        collisions=collisions,
        energy_flux_matrix=energy_flux,
        collision_matrix=native.collision_matrix,
        matrix=(energy_flux + native.collision_matrix).tocsr(),
    )


def assemble_coefficient_pn_operator(
    case: MultiTermCase,
    native: KineticOperatorBlock,
    lmax: int,
    *,
    allow_moment_table: bool = False,
) -> sparse.csr_matrix:
    """Assemble the limited ordinary-XS higher-l angular-closure PN block."""

    if lmax < 2:
        raise ValueError("coefficient-space PN assembly is only for lmax > 1")
    if case.config.physics.field.type != "dc":
        raise NotImplementedError("pn_closure_direct lmax>1 supports dc fields only")
    magnetic = case.config.physics.field.magnetic_field
    if magnetic.enabled and magnetic.B_T > 0.0:
        raise NotImplementedError(
            "pn_closure_direct lmax>1 does not support B-field"
        )
    if any(proc.process_type == ProcessType.SUPERELASTIC for proc in case.cross_sections.processes):
        raise NotImplementedError(
            "pn_closure_direct lmax>1 does not support superelastic l>0 treatment"
        )
    angular = build_angular_model(case.config)
    if angular.name == "moment_table" and not allow_moment_table:
        raise NotImplementedError(
            "pn_closure_direct lmax>1 ordinary-XS scope does not use moment_table"
        )
    moments = angular.moments(
        native.grid.energy_eV,
        lmax,
        sigma_total=native.collisions.sigma_total_like,
        sigma_momentum=native.collisions.sigma_m,
    )
    damping = build_higher_l_collision_damping(native, moments)
    n = len(native.grid.energy_eV)
    zero = sparse.csr_matrix((n, n))
    identity = sparse.identity(n, format="csr")
    blocks: list[list[sparse.csr_matrix]] = [
        [zero.copy() for _ in range(lmax + 1)] for _ in range(lmax + 1)
    ]

    blocks[0][0] = native.collision_matrix
    blocks[0][1] = identity

    inv_damp = [
        sparse.diags(1.0 / np.maximum(damping[ell], 1.0), format="csr")
        for ell in range(lmax + 1)
    ]
    field = native.energy_flux_matrix
    blocks[1][0] = -field
    blocks[1][1] = identity
    for ell in range(2, lmax + 1):
        lower = ell / (2.0 * ell - 1.0)
        upper = (ell + 1.0) / (2.0 * ell + 3.0)
        blocks[ell][ell] = identity
        blocks[ell][ell - 1] = lower * (inv_damp[ell] @ field)
        if ell < lmax:
            blocks[ell][ell + 1] = upper * (inv_damp[ell] @ field)
    return sparse.bmat(blocks, format="csr")


def _normalize_signed(
    f0: np.ndarray,
    f1: np.ndarray,
    widths: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    total = weighted_integral(f0, widths)
    if not np.isfinite(total) or abs(total) < 1.0e-300:
        raise FloatingPointError("direct PN solve produced a non-normalizable f0")
    if total < 0.0:
        f0 = -f0
        f1 = -f1
        total = -total
    return f0 / total, f1 / total


def _solve_block_system(
    block_matrix: sparse.csr_matrix,
    widths: np.ndarray,
    growth: float,
) -> tuple[np.ndarray, np.ndarray]:
    n = len(widths)
    identity = sparse.identity(n, format="csr")
    l0_left = block_matrix[:n, :n] - growth * identity
    l0_right = identity
    l1_left = -block_matrix[n:, :n]
    l1_right = identity
    mat = sparse.bmat(
        [[l0_left, l0_right], [l1_left, l1_right]],
        format="lil",
    )
    rhs = np.zeros(2 * n, dtype=float)
    diag_abs = np.abs((l0_left + block_matrix[n:, :n]).diagonal())
    row = int(np.argmax(diag_abs)) if np.any(diag_abs > 0.0) else n - 1
    mat[row, :] = 0.0
    mat[row, :n] = widths
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
            maxiter=max(1000, 8 * n),
        )[0]
    if not np.all(np.isfinite(sol)):
        raise FloatingPointError("direct PN block solve returned non-finite values")
    return sol[:n], sol[n:]


def _solve_extended_block_system(
    block_matrix: sparse.csr_matrix,
    widths: np.ndarray,
    lmax: int,
    growth: float,
) -> np.ndarray:
    n = len(widths)
    mat = (block_matrix - growth * sparse.block_diag(
        [sparse.identity(n, format="csr")]
        + [sparse.csr_matrix((n, n)) for _ in range(lmax)],
        format="csr",
    )).tolil()
    rhs = np.zeros((lmax + 1) * n, dtype=float)
    diag_abs = np.abs(block_matrix[:n, :n].diagonal())
    row = int(np.argmax(diag_abs)) if np.any(diag_abs > 0.0) else n - 1
    mat[row, :] = 0.0
    mat[row, :n] = widths
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
            maxiter=max(1000, 8 * (lmax + 1) * n),
        )[0]
    if not np.all(np.isfinite(sol)):
        raise FloatingPointError("direct PN lmax>1 solve returned non-finite values")
    if negative_mass_fraction(sol[:n], widths) > _HIGHER_L_NEGATIVE_MASS_LIMIT:
        lower = np.full_like(sol, -np.inf)
        lower[:n] = 0.0
        upper = np.full_like(sol, np.inf)
        bounded = optimize.lsq_linear(
            csr,
            rhs,
            bounds=(lower, upper),
            lsmr_tol="auto",
            max_iter=200,
        )
        if not bounded.success or not np.all(np.isfinite(bounded.x)):
            raise RuntimeError("direct PN lmax>1 bounded sparse solve failed")
        sol = bounded.x
    return sol.reshape((lmax + 1, n))


def _transport_set_from_eedf(
    case: MultiTermCase,
    native: KineticOperatorBlock,
    rates: RateSet,
    f0: np.ndarray,
    solver_name: str,
) -> ElectronTransport:
    transport = transport_from_eedf(
        case.config,
        case.cross_sections,
        native.grid.energy_eV,
        native.grid.widths_eV,
        f0,
        native.gas_number_density_m3,
        case.e_over_n_Td,
        case.two_term_config,
    )
    return ElectronTransport(
        definition=TRANSPORT_F0_GRADIENT_RECONSTRUCTION,
        gas_number_density_m3=transport.gas_number_density_m3,
        drift_velocity_m_s=transport.drift_velocity_m_s,
        reduced_mobility_m2_V_s_m3=transport.reduced_mobility_m2_V_s_m3,
        reduced_diffusion_L_m2_s_m3=transport.reduced_diffusion_L_m2_s_m3,
        reduced_diffusion_T_m2_s_m3=transport.reduced_diffusion_T_m2_s_m3,
    )


def _rate_set_from_f0(
    case: MultiTermCase,
    native: KineticOperatorBlock,
    f0: np.ndarray,
    case_id: str,
    solver_name: str,
) -> RateSet:
    rate_data = compute_rates_from_eedf(
        case.config,
        case.cross_sections,
        native.grid.energy_eV,
        native.grid.widths_eV,
        f0,
        case_id=case_id,
        e_over_n_Td=case.e_over_n_Td,
        solver_name=solver_name,
    )
    return RateSet(
        tuple(rate_data.rates),
        case.gas_number_density_m3 * rate_data.ionization_rate_m3_s,
        case.gas_number_density_m3 * rate_data.attachment_rate_m3_s,
    )


def _solve_direct_higher_l(
    case: MultiTermCase,
    case_id: str,
    solver_name: str,
    native: KineticOperatorBlock,
    lmax: int,
    *,
    method_used: str = "pn_closure_direct",
    metadata_override: dict[str, object] | None = None,
    angular_moments: np.ndarray | None = None,
    allow_moment_table: bool = False,
) -> MultiTermSolution:
    block = assemble_coefficient_pn_operator(
        case,
        native,
        lmax,
        allow_moment_table=allow_moment_table,
    )
    widths = native.grid.widths_eV
    energy = native.grid.energy_eV
    edges = native.grid.edges_eV
    conv = case.two_term_config.convergence
    growth = 0.0
    previous_f0: np.ndarray | None = None
    residual = np.inf
    coefficients = np.zeros((lmax + 1, len(energy)), dtype=float)
    negative_mass = np.inf
    for _ in range(1, conv.max_iterations + 1):
        coefficients = _solve_extended_block_system(block, widths, lmax, growth)
        f0 = coefficients[0]
        total = weighted_integral(f0, widths)
        if not np.isfinite(total) or abs(total) < 1.0e-300:
            raise FloatingPointError("direct PN lmax>1 produced non-normalizable F0")
        if total < 0.0:
            coefficients *= -1.0
            total = -total
        coefficients /= total
        f0 = coefficients[0]
        negative_mass = negative_mass_fraction(f0, widths)
        if negative_mass > _HIGHER_L_NEGATIVE_MASS_LIMIT:
            raise RuntimeError(
                "direct PN lmax>1 produced nonphysical negative F0 mass "
                f"({negative_mass:.3e})"
            )
        op = block @ coefficients.reshape(-1)
        new_growth = weighted_integral(op[: len(energy)], widths)
        residual_vector = op.copy()
        residual_vector[: len(energy)] -= new_growth * f0
        residual = float(
            np.sum(np.abs(residual_vector[: len(energy)]) * widths)
            / max(np.sum(np.abs(op[: len(energy)]) * widths), abs(new_growth), 1.0)
        )
        shape_change = (
            float(np.sum(np.abs(f0 - previous_f0) * widths))
            if previous_f0 is not None
            else np.inf
        )
        growth_change = abs(new_growth - growth) / max(
            abs(new_growth), abs(growth), 1.0
        )
        previous_f0 = f0.copy()
        growth = conv.relaxation * new_growth + (1.0 - conv.relaxation) * growth
        if (
            shape_change < conv.tolerance
            and growth_change < conv.eigenvalue_tolerance
            and residual < max(conv.residual_tolerance, 1.0e-8)
        ):
            growth = new_growth
            break

    f0 = coefficients[0]
    rates = _rate_set_from_f0(case, native, f0, case_id, solver_name)
    transport_set = _transport_set_from_eedf(case, native, rates, f0, solver_name)
    tail_mask = energy >= 0.9 * edges[-1]
    diagnostics = SolverDiagnostics(
        warnings=(),
        mean_energy_eV=mean_energy_from_eedf(energy, widths, f0),
        eedf_tail_fraction=float(np.sum(f0[tail_mask] * widths[tail_mask])),
        power_balance_relative_residual=residual,
    )
    metadata = {
        "solver_method": method_used,
        "physics_level": "angular_closure",
        "direct_pn_operator": True,
        "exact_dcs_based": False,
        "ordinary_integral_xs_closure": True,
        "lmax": int(lmax),
        "angular_model": case.config.physics.angular_scattering.model,
        "pn_residual": float(residual),
        "negative_mass_fraction": float(negative_mass),
    }
    if metadata_override:
        metadata.update(metadata_override)
    return MultiTermSolution(
        energy,
        widths,
        coefficients,
        f0,
        rates,
        transport_set,
        diagnostics,
        method_used=method_used,
        metadata=metadata,
        angular_moments=angular_moments,
    )


def solve_direct_lmax1(
    case: MultiTermCase,
    case_id: str,
    solver_name: str,
    *,
    native_override: KineticOperatorBlock | None = None,
    method_used: str = "pn_closure_direct",
    metadata_override: dict[str, object] | None = None,
    angular_moments: np.ndarray | None = None,
    allow_moment_table: bool = False,
) -> MultiTermSolution:
    """Solve the lmax=1 direct PN reduction without using two-term f0/rates."""

    lmax = int(case.config.solvers.multi_term.lmax)
    native = native_override
    if native is None:
        kinetic_grid = make_two_term_energy_grid(
            case.two_term_config,
            cross_sections=case.cross_sections,
        )
        native = assemble_native_operator_blocks(
            case.config,
            case.cross_sections,
            case.e_over_n_Td,
            kinetic_grid,
            case.two_term_config,
        )
    if lmax != 1:
        return _solve_direct_higher_l(
            case,
            case_id,
            solver_name,
            native,
            lmax,
            method_used=method_used,
            metadata_override=metadata_override,
            angular_moments=angular_moments,
            allow_moment_table=allow_moment_table,
        )

    energy = native.grid.energy_eV
    edges = native.grid.edges_eV
    widths = native.grid.widths_eV
    n = len(energy)
    system_core = sparse.bmat(
        [
            [native.collision_matrix, sparse.identity(n, format="csr")],
            [native.energy_flux_matrix, sparse.csr_matrix((n, n))],
        ],
        format="csr",
    )

    conv = case.two_term_config.convergence
    growth = 0.0
    previous_f0: np.ndarray | None = None
    residual = np.inf
    negative_mass = np.inf
    f0 = np.zeros(n, dtype=float)
    f1 = np.zeros(n, dtype=float)
    for iteration in range(1, conv.max_iterations + 1):
        f0, f1 = _solve_block_system(system_core, widths, growth)
        f0, f1 = _normalize_signed(f0, f1, widths)
        negative_mass = negative_mass_fraction(f0, widths)
        if negative_mass > _NEGATIVE_MASS_LIMIT:
            raise RuntimeError(
                "direct PN solve produced nonphysical negative f0 mass "
                f"({negative_mass:.3e})"
            )
        if negative_mass > 0.0:
            f0 = np.clip(f0, 0.0, None)
            f0 = normalize_eedf(f0, widths)
            f1 = native.energy_flux_matrix @ f0

        op_f0 = native.matrix @ f0
        new_growth = weighted_integral(op_f0, widths)
        eig_res = op_f0 - new_growth * f0
        residual = float(
            np.sum(np.abs(eig_res) * widths)
            / max(np.sum(np.abs(op_f0) * widths), abs(new_growth), 1.0)
        )
        shape_change = (
            float(np.sum(np.abs(f0 - previous_f0) * widths))
            if previous_f0 is not None
            else np.inf
        )
        growth_change = abs(new_growth - growth) / max(
            abs(new_growth), abs(growth), 1.0
        )
        previous_f0 = f0.copy()
        growth = conv.relaxation * new_growth + (1.0 - conv.relaxation) * growth
        if (
            shape_change < conv.tolerance
            and growth_change < conv.eigenvalue_tolerance
            and residual < conv.residual_tolerance
        ):
            growth = new_growth
            break

    rates = _rate_set_from_f0(case, native, f0, case_id, solver_name)
    transport_set = _transport_set_from_eedf(case, native, rates, f0, solver_name)
    coefficients = np.zeros((2, n), dtype=float)
    coefficients[0] = f0
    coefficients[1] = f1
    tail_mask = energy >= 0.9 * edges[-1]
    tail = float(np.sum(f0[tail_mask] * widths[tail_mask]))
    mean_energy = mean_energy_from_eedf(energy, widths, f0)
    diagnostics = SolverDiagnostics(
        warnings=(),
        mean_energy_eV=mean_energy,
        eedf_tail_fraction=tail,
        power_balance_relative_residual=residual,
    )
    metadata = {
        "solver_method": method_used,
        "physics_level": "angular_closure",
        "direct_pn_operator": True,
        "exact_dcs_based": False,
        "ordinary_integral_xs_closure": True,
        "lmax": 1,
        "angular_model": case.config.physics.angular_scattering.model,
        "pn_residual": float(residual),
        "negative_mass_fraction": float(negative_mass),
        "grid_spacing": case.two_term_config.energy_grid.spacing,
        "threshold_refined": bool(
            case.two_term_config.energy_grid.refine.enabled
            and len(energy) > case.two_term_config.energy_grid.n
        ),
    }
    if method_used == "pn_closure_direct":
        metadata["lmax1_regression_target"] = "two_term"
    if metadata_override:
        metadata.update(metadata_override)
    return MultiTermSolution(
        energy,
        widths,
        coefficients,
        f0,
        rates,
        transport_set,
        diagnostics,
        method_used=method_used,
        metadata=metadata,
        angular_moments=angular_moments,
    )


def solve_pn_dcs(
    case: MultiTermCase,
    case_id: str,
    solver_name: str,
) -> MultiTermSolution:
    """Run pn_dcs with validated normalized Legendre moment-table input."""

    if case.config.physics.angular_scattering.model != "moment_table":
        raise NotImplementedError(
            "multi_term method 'pn_dcs' requires physics.angular_scattering.model=moment_table"
        )
    lmax = int(case.config.solvers.multi_term.lmax)
    kinetic_grid = make_two_term_energy_grid(
        case.two_term_config,
        cross_sections=case.cross_sections,
    )
    native = assemble_native_operator_blocks(
        case.config,
        case.cross_sections,
        case.e_over_n_Td,
        kinetic_grid,
        case.two_term_config,
    )
    angular = build_angular_model(case.config)
    moments = angular.moments(native.grid.energy_eV, max(lmax, 1))
    native = _native_with_moment_table_l1_damping(native, moments)
    angular_metadata = angular.metadata()
    exact = bool(angular_metadata.get("exact_dcs_based", False))
    metadata_override = {
        "solver_method": "pn_dcs",
        "physics_level": "dcs_moment_based" if exact else "table_moment_based",
        "direct_pn_operator": True,
        "angular_model": "moment_table",
        "angular_moment_source": "moment_table",
        "moment_table_provenance": angular_metadata.get(
            "moment_table_provenance", "unknown"
        ),
        "exact_dcs_based": exact,
        "ordinary_integral_xs_closure": False,
    }
    return solve_direct_lmax1(
        case,
        case_id,
        solver_name,
        native_override=native,
        method_used="pn_dcs",
        metadata_override=metadata_override,
        angular_moments=moments[: lmax + 1],
        allow_moment_table=True,
    )
