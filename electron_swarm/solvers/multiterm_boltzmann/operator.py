"""Sparse multi-term operator backend.

This module keeps the public backend, lmax=1 two-term compatibility path, and
result extraction. The executable lmax>1 assembly/solve core lives in
:mod:`operator_core` so the extension surface stays small.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.cross_sections import CrossSectionSet
from electron_swarm.core.results import RateResult, SwarmCaseResult
from electron_swarm.core.transport import FluxTransport, TransportMetadata, TransportSet
from electron_swarm.solvers.boltzmann_two_term import (
    BoltzmannTwoTermSolver,
    NativeOperatorBlock,
)

from .diagnostics import SolverDiagnostics
from .grid import electron_speed_m_s
from .hydro import HydroConfig, HydrodynamicModeSolver
from .models import MultiTermCase, MultiTermSolution
from .models import RateSet
from .operator_core import (
    DensityNormalizationConstraint,
    LegendreBlockLayout,
    OperatorAssemblyDiagnostics,
    OperatorSolveState,
    OperatorSystem,
    assemble_operator_system,
    build_density_normalization_constraint,
    build_operator_assembly_diagnostics,
    solve_operator_system,
)
from .projection import (
    compute_rate_set,
    validity_warnings,
)

DEFAULT_LMAX1_TOLERANCES = {
    "mean_energy_eV": 0.35,
    "drift_velocity_m_s": 0.75,
    "diffusion_L_m2_s": 0.75,
    "net_ionization_frequency_s": 1.50,
}


class OperatorBackendUnavailable(NotImplementedError):
    """Raised when an unsupported operator method is requested."""


@dataclass(frozen=True, slots=True)
class MetricComparison:
    metric: str
    candidate: float
    reference: float
    relative_difference: float
    tolerance: float

    @property
    def passed(self) -> bool:
        if not np.isfinite(self.relative_difference):
            return False
        return abs(self.relative_difference) <= self.tolerance


@dataclass(frozen=True, slots=True)
class LmaxOneReferenceReport:
    candidate_solver: str
    reference_solver: str
    comparisons: tuple[MetricComparison, ...]

    @property
    def passed(self) -> bool:
        return all(item.passed for item in self.comparisons)

    def as_metadata(self, prefix: str = "lmax1_reference") -> dict[str, float | bool]:
        metadata: dict[str, float | bool] = {f"{prefix}_passed": self.passed}
        for item in self.comparisons:
            metadata[f"{prefix}_{item.metric}_relative_difference"] = (
                item.relative_difference
            )
        return metadata


class MultiTermOperatorBackend:
    """Operator backend for lmax=1 reference and production lmax>1 runs."""

    method_name = "operator"

    def __init__(self, solver_name: str = "multiterm_boltzmann") -> None:
        self.solver_name = solver_name

    def solve(
        self,
        case: MultiTermCase,
        case_id: str,
        previous: MultiTermSolution | None = None,
    ) -> MultiTermSolution:
        lmax = int(case.config.multiterm_boltzmann.lmax)
        if lmax == 1:
            return solve_lmax1_operator_compat(case, case_id, self.solver_name)
        return solve_operator_lmax_gt1(case, case_id, self.solver_name)


def solve_lmax1_two_term_reference(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    e_over_n_Td: float,
    case_id: str = "lmax1_reference",
) -> SwarmCaseResult:
    """Run the native two-term block as the required lmax=1 reference."""

    return BoltzmannTwoTermSolver(config, cross_sections).solve_native_reference_case(
        e_over_n_Td, case_id
    )


def assemble_lmax1_native_operator_block(
    config: SwarmConfig,
    cross_sections: CrossSectionSet,
    e_over_n_Td: float,
    *,
    max_eV_override: float | None = None,
    n_override: int | None = None,
) -> NativeOperatorBlock:
    """Assemble the reusable lmax=1 Scharfetter-Gummel operator block.

    Future multi-term operator work should treat this as the regression anchor:
    the general Legendre block system must reduce to this native two-term block
    when only the l=0/1 terms are retained.
    """

    solver = BoltzmannTwoTermSolver(config, cross_sections)
    energy, edges, widths = solver.make_native_energy_grid(
        max_eV_override=max_eV_override,
        n_override=n_override,
    )
    return solver.assemble_native_operator_block(e_over_n_Td, energy, edges, widths)


def _velocity_moment_vector(case: MultiTermCase, layout: LegendreBlockLayout) -> np.ndarray:
    moment = np.zeros(layout.n_unknowns, dtype=float)
    if layout.n_legendre_terms > 1:
        moment[layout.slice_for_l(1)] = case.grid.speeds_m_s * case.grid.widths_eV / 3.0
    return moment


def _hydrodynamic_k_values(
    case: MultiTermCase, system: OperatorSystem, f0: np.ndarray
) -> tuple[float, ...]:
    length = case.config.conditions.length_scale_m
    k0 = 1.0e-4
    if length is not None and length > 0.0:
        k0 = min(k0, 1.0e-3 / length)
    else:
        mean_speed = float(np.sum(case.grid.speeds_m_s * f0 * case.grid.widths_eV))
        mean_nu = float(
            np.sum(
                np.maximum(system.collisions.total_frequency_s_inv, 0.0)
                * f0
                * case.grid.widths_eV
            )
        )
        if mean_speed > 0.0 and mean_nu > 0.0:
            k0 = min(k0, 1.0e-3 * mean_nu / mean_speed)
    k0 = max(k0, 1.0e-8)
    return (-2.0 * k0, -k0, 0.0, k0, 2.0 * k0)


def _coefficients_from_state(
    case: MultiTermCase, system: OperatorSystem, state: OperatorSolveState
) -> tuple[np.ndarray, np.ndarray]:
    layout = system.layout
    coeff = state.coefficients_flat.reshape(
        layout.n_legendre_terms, layout.n_energy_cells
    )
    f0 = case.grid.normalize_energy_pdf(np.clip(coeff[0], 0.0, None))
    coeff = coeff.copy()
    coeff[0] = f0
    return coeff, f0


def _operator_quality_metadata(
    case: MultiTermCase,
    system: OperatorSystem,
    coeff: np.ndarray,
    f0: np.ndarray,
) -> tuple[dict[str, float | bool], tuple[str, ...]]:
    """Compute compact solution-quality indicators without new public contracts."""

    widths = case.grid.widths_eV
    tail_mask = case.grid.centers_eV >= 0.9 * case.grid.edges_eV[-1]
    tail_probability = float(np.sum(f0[tail_mask] * widths[tail_mask]))
    with np.errstate(under="ignore"):
        total_rate = float(
            np.sum(system.collisions.total_frequency_s_inv * f0 * widths)
        )
        tail_rate = float(
            np.sum(
                system.collisions.total_frequency_s_inv[tail_mask]
                * f0[tail_mask]
                * widths[tail_mask]
            )
        )
    tail_rate_fraction = tail_rate / total_rate if total_rate > 0.0 else 0.0
    f0_l1 = max(float(np.sum(np.abs(f0) * widths)), 1.0e-300)
    f0_l2 = max(float(np.sum(f0 * f0 * widths)), 1.0e-300)
    highest = coeff[-1]
    highest_l1 = float(np.sum(np.abs(highest) * widths)) / f0_l1
    highest_l2 = float(np.sqrt(np.sum(highest * highest * widths) / f0_l2))
    tolerance = float(case.config.multiterm_boltzmann.lmax_convergence_tolerance)
    warnings: list[str] = []
    if highest_l1 > tolerance:
        warnings.append("operator_lmax_may_be_underresolved")
    metadata = {
        "operator_tail_probability": tail_probability,
        "operator_tail_rate_fraction": float(tail_rate_fraction),
        "operator_highest_l_relative_l1": highest_l1,
        "operator_highest_l_relative_l2": highest_l2,
        "operator_lmax_convergence_tolerance": tolerance,
        "operator_lmax_convergence_ok": highest_l1 <= tolerance,
    }
    return metadata, tuple(warnings)


def _ionization_source_model(config: SwarmConfig) -> str:
    model = config.boltzmann_two_term.ionization_energy_sharing
    if model == "equal":
        return "two_term_native_equal_sharing"
    if model == "primary_secondary":
        return "two_term_native_primary_secondary"
    return "two_term_native_loss_only"


def solve_operator_lmax_gt1(
    case: MultiTermCase, case_id: str, solver_name: str = "multiterm_boltzmann"
) -> MultiTermSolution:
    """Run the production lmax>1 sparse-operator backend."""

    system = assemble_operator_system(case)
    if (
        case.config.multiterm_boltzmann.hydrodynamic
        and system.diagnostics.nonphysical_field_scaling
    ):
        raise RuntimeError(
            "multiterm_boltzmann.hydrodynamic requires field_coupling_scale=1.0"
        )
    state = solve_operator_system(case, system)
    layout = system.layout
    coeff, f0 = _coefficients_from_state(case, system, state)
    drift = 0.0
    if layout.n_legendre_terms > 1:
        drift = float(np.sum(case.grid.speeds_m_s * coeff[1] * case.grid.widths_eV / 3.0))
    mobility = (
        drift / case.electric_field_V_m
        if case.electric_field_V_m != 0.0
        else float("nan")
    )
    rates = compute_rate_set(case, system.collisions, f0, case_id, solver_name)
    flux = FluxTransport.with_characteristic_energies(
        drift,
        mobility,
        diffusion_longitudinal_m2_s=None,
        diffusion_transverse_m2_s=None,
    )
    transport = TransportSet.from_flux_only(
        flux,
        rates.ionization_frequency_s_inv,
        rates.attachment_frequency_s_inv,
        TransportMetadata(
            solver=solver_name,
            coefficient_definition="flux",
            swarm_condition="local_flux",
            notes=("operator_flux_b0_dc_m0_integral_cross_sections",),
        ),
    )
    method_used = "operator_flux"
    if case.config.multiterm_boltzmann.hydrodynamic:
        hydro_config = HydroConfig(
            k_values=_hydrodynamic_k_values(case, system, f0),
            dense_threshold=int(case.config.multiterm_boltzmann.dense_threshold),
        )
        hydro_fit, transport = HydrodynamicModeSolver().solve(
            system.matrix,
            system.streaming_matrix,
            _velocity_moment_vector(case, layout),
            case.electric_field_V_m,
            system.normalization.weights,
            rates.ionization_frequency_s_inv,
            rates.attachment_frequency_s_inv,
            hydro_config,
        )
        method_used = "operator_hydrodynamic"
    else:
        hydro_fit = None
    tail_mask = case.grid.centers_eV >= 0.9 * case.grid.edges_eV[-1]
    tail = float(np.sum(f0[tail_mask] * case.grid.widths_eV[tail_mask]))
    quality_metadata, quality_warnings = _operator_quality_metadata(
        case, system, coeff, f0
    )
    warnings = [
        "operator_scope_b0_dc_m0_integral_cross_sections",
        "operator_l_gt_0_inelastic_source_isotropic_l0_only",
    ]
    if method_used == "operator_flux":
        warnings.extend(
            [
                "operator_diffusion_not_computed",
                "operator_bulk_source_gradient_not_computed",
            ]
        )
    if system.diagnostics.nonphysical_field_scaling:
        warnings.append("operator_nonphysical_field_scaling")
    warnings.extend(validity_warnings(case, system.collisions, f0))
    warnings.extend(quality_warnings)
    warnings.extend(state.warnings)
    diagnostics = SolverDiagnostics(
        warnings=tuple(dict.fromkeys(warnings)),
        mean_energy_eV=float(np.sum(case.grid.centers_eV * f0 * case.grid.widths_eV)),
        eedf_tail_fraction=tail,
        power_balance_relative_residual=state.residual_L1,
    )
    metadata = {
        "operator_residual_L1": state.residual_L1,
        "operator_growth_frequency_s-1": state.growth_frequency_s_inv,
        "operator_iterations": state.iterations,
        "operator_converged": state.converged,
        "operator_lmax": layout.lmax,
        "operator_coefficient_order": system.normalization.coefficient_order,
        "operator_warnings": "; ".join(tuple(dict.fromkeys(warnings))),
        "operator_matrix_nnz": int(system.matrix.nnz),
        "operator_base_matrix_nnz": int(system.base_matrix.nnz),
        "operator_field_coupling_matrix_nnz": int(system.field_coupling_matrix.nnz),
        "operator_streaming_matrix_nnz": int(system.streaming_matrix.nnz),
        "operator_normalization": system.normalization.definition,
        "operator_ionization_source_model": _ionization_source_model(case.config),
        "operator_l_gt_0_inelastic_model": "sink_only_isotropic_l0_source",
    }
    metadata.update(quality_metadata)
    if hydro_fit is not None:
        metadata.update(
            {
                "operator_hydro_fit_residual": hydro_fit.fit_residual,
                "operator_hydro_symmetry_error": hydro_fit.symmetry_error,
                "operator_hydro_mode_continuity_error": (
                    hydro_fit.mode_continuity_error
                ),
                "operator_hydro_k_min_1_m": float(np.min(hydro_fit.k_values)),
                "operator_hydro_k_max_1_m": float(np.max(hydro_fit.k_values)),
                "operator_hydro_growth_frequency_s-1": hydro_fit.nu_eff_s_inv,
            }
        )
    metadata.update(system.diagnostics.as_metadata())
    return MultiTermSolution(
        case.grid.centers_eV,
        case.grid.widths_eV,
        coeff,
        f0,
        rates,
        transport,
        None,
        diagnostics,
        method_used=method_used,
        metadata=metadata,
    )


def _widths_from_centers(energy_eV: np.ndarray) -> np.ndarray:
    energy = np.asarray(energy_eV, dtype=float)
    if len(energy) < 2:
        return np.ones_like(energy)
    edges = np.empty(len(energy) + 1, dtype=float)
    edges[1:-1] = 0.5 * (energy[:-1] + energy[1:])
    edges[0] = max(0.0, energy[0] - 0.5 * (energy[1] - energy[0]))
    edges[-1] = energy[-1] + 0.5 * (energy[-1] - energy[-2])
    return np.diff(edges)


def _l1_coefficient_from_reference(
    energy_eV: np.ndarray,
    widths_eV: np.ndarray,
    eedf_eV_inv: np.ndarray,
    drift_velocity_m_s: float,
) -> np.ndarray:
    if len(energy_eV) >= 3:
        derivative = np.gradient(eedf_eV_inv, energy_eV, edge_order=2)
    else:
        derivative = np.gradient(eedf_eV_inv, energy_eV)
    anis = -derivative
    norm = np.max(np.abs(anis)) or 1.0
    anis = anis / norm * (np.max(eedf_eV_inv) or 1.0)
    speeds = electron_speed_m_s(energy_eV)
    moment = float(np.sum(speeds * anis * widths_eV / 3.0))
    scale = drift_velocity_m_s / moment if abs(moment) > 1.0e-300 else 0.0
    return anis * scale


def _convert_reference_rates(
    reference: SwarmCaseResult, solver_name: str, case_id: str, number_density: float
) -> RateSet:
    rates: list[RateResult] = []
    ion_freq = 0.0
    att_freq = 0.0
    for rate in reference.rates:
        freq = number_density * rate.mixture_weighted_rate_m3_s
        if rate.process_type == "ionization":
            ion_freq += freq
        elif rate.process_type == "attachment":
            att_freq += freq
        rates.append(
            RateResult(
                solver=solver_name,
                case_id=case_id,
                e_over_n_Td=reference.e_over_n_Td,
                species=rate.species,
                process=rate.process,
                process_type=rate.process_type,
                threshold_eV=rate.threshold_eV,
                rate_coefficient_m3_s=rate.rate_coefficient_m3_s,
                mixture_weighted_rate_m3_s=rate.mixture_weighted_rate_m3_s,
                frequency_s_inv=freq,
                power_loss_eV_s=rate.power_loss_eV_s,
            )
        )
    return RateSet(tuple(rates), float(ion_freq), float(att_freq))


def solve_lmax1_operator_compat(
    case: MultiTermCase, case_id: str, solver_name: str = "multiterm_boltzmann"
) -> MultiTermSolution:
    """Execute the current operator method through the validated two-term path.

    This is the first executable operator milestone: lmax=1 uses the existing
    Scharfetter-Gummel two-term implementation as the reference block. The
    general l>1 sparse block system is handled by the production path.
    """

    reference = solve_lmax1_two_term_reference(
        case.config, case.cross_sections, case.e_over_n_Td, case_id
    )
    energy = np.asarray(reference.energy_eV, dtype=float)
    widths = _widths_from_centers(energy)
    eedf = np.asarray(reference.eedf, dtype=float)
    coeff = np.zeros((2, len(energy)), dtype=float)
    coeff[0] = eedf
    coeff[1] = _l1_coefficient_from_reference(
        energy, widths, eedf, reference.drift_velocity_m_s
    )
    flux = FluxTransport.with_characteristic_energies(
        reference.drift_velocity_m_s,
        reference.mobility_m2_V_s,
        diffusion_longitudinal_m2_s=reference.diffusion_L_m2_s,
        diffusion_transverse_m2_s=reference.diffusion_T_m2_s,
    )
    rate_set = _convert_reference_rates(
        reference, solver_name, case_id, case.gas_number_density_m3
    )
    transport = TransportSet.from_flux_only(
        flux,
        rate_set.ionization_frequency_s_inv,
        rate_set.attachment_frequency_s_inv,
        TransportMetadata(
            solver=solver_name,
            coefficient_definition="flux",
            swarm_condition="local_flux",
            notes=("operator_lmax1_two_term_reference",),
        ),
    )
    tail_mask = energy >= 0.9 * float(np.max(energy))
    tail = float(np.sum(eedf[tail_mask] * widths[tail_mask]))
    diagnostics = SolverDiagnostics(
        warnings=("operator_lmax1_uses_two_term_reference_backend",),
        mean_energy_eV=float(reference.mean_energy_eV),
        eedf_tail_fraction=tail,
    )
    return MultiTermSolution(
        energy,
        widths,
        coeff,
        eedf,
        rate_set,
        transport,
        None,
        diagnostics,
        method_used="operator_lmax1_two_term",
    )


def compare_lmax1_transport_to_reference(
    candidate: SwarmCaseResult,
    reference: SwarmCaseResult,
    tolerances: Mapping[str, float] | None = None,
) -> LmaxOneReferenceReport:
    """Compare a candidate multi-term lmax=1 case to the two-term reference."""

    active_tolerances = dict(DEFAULT_LMAX1_TOLERANCES)
    if tolerances:
        active_tolerances.update({str(k): float(v) for k, v in tolerances.items()})

    comparisons = []
    for metric, tolerance in active_tolerances.items():
        candidate_value = float(getattr(candidate, metric))
        reference_value = float(getattr(reference, metric))
        denom = max(abs(reference_value), 1.0e-300)
        relative = (candidate_value - reference_value) / denom
        comparisons.append(
            MetricComparison(
                metric=metric,
                candidate=candidate_value,
                reference=reference_value,
                relative_difference=float(relative),
                tolerance=float(tolerance),
            )
        )
    return LmaxOneReferenceReport(
        candidate_solver=candidate.solver,
        reference_solver=reference.solver,
        comparisons=tuple(comparisons),
    )
