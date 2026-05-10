"""Reference-gated multi-term operator backend.

This module keeps the public backend, lmax=1 two-term compatibility path, and
result extraction. For lmax>1 with integral cross sections only, public results
are anchored to the validated two-term solve and expose bounded higher-order
coefficients as diagnostics.
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
from .models import MultiTermCase, MultiTermSolution
from .models import RateSet
from .projection import (
    project_collision_data,
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
    """Operator backend for lmax=1 reference and experimental lmax>1 runs."""

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


def _interpolate_collision_frequency(
    case: MultiTermCase, energy_eV: np.ndarray, field: str
) -> np.ndarray:
    collisions = project_collision_data(case)
    values = np.asarray(getattr(collisions, field), dtype=float)
    return np.interp(
        energy_eV,
        case.grid.centers_eV,
        values,
        left=float(values[0]),
        right=float(values[-1]),
    )


def _anchored_lmax_coefficients(
    case: MultiTermCase,
    energy_eV: np.ndarray,
    widths_eV: np.ndarray,
    eedf_eV_inv: np.ndarray,
    drift_velocity_m_s: float,
) -> np.ndarray:
    """Build bounded l>1 diagnostic coefficients anchored to the two-term EEDF."""

    lmax = int(case.config.multiterm_boltzmann.lmax)
    coeff = np.zeros((lmax + 1, len(energy_eV)), dtype=float)
    coeff[0] = eedf_eV_inv
    if lmax >= 1:
        coeff[1] = _l1_coefficient_from_reference(
            energy_eV, widths_eV, eedf_eV_inv, drift_velocity_m_s
        )
    if lmax < 2:
        return coeff

    momentum = np.maximum(
        _interpolate_collision_frequency(
            case, energy_eV, "momentum_frequency_s_inv"
        ),
        1.0,
    )
    field_advection = (
        float(case.config.multiterm_boltzmann.field_coupling_scale)
        * case.electric_field_V_m
        * electron_speed_m_s(energy_eV)
    )
    for ell in range(2, lmax + 1):
        source_coupling = ell / (2 * ell - 1)
        derivative = np.gradient(coeff[ell - 1], energy_eV, edge_order=2)
        raw = -source_coupling * field_advection * derivative / (ell * momentum)
        raw_l1 = float(np.sum(np.abs(raw) * widths_eV))
        previous_l1 = float(np.sum(np.abs(coeff[ell - 1]) * widths_eV))
        max_l1 = 0.35 * previous_l1
        if raw_l1 > max_l1 > 0.0:
            raw *= max_l1 / raw_l1
        coeff[ell] = np.nan_to_num(raw, nan=0.0, posinf=0.0, neginf=0.0)
    return coeff


def _anchored_quality_metadata(
    case: MultiTermCase,
    energy_eV: np.ndarray,
    widths_eV: np.ndarray,
    coeff: np.ndarray,
    f0: np.ndarray,
) -> tuple[dict[str, float | bool], tuple[str, ...]]:
    total_frequency = np.maximum(
        _interpolate_collision_frequency(case, energy_eV, "total_frequency_s_inv"),
        0.0,
    )
    tail_mask = energy_eV >= 0.9 * float(np.max(energy_eV))
    tail_probability = float(np.sum(f0[tail_mask] * widths_eV[tail_mask]))
    total_rate = float(np.sum(total_frequency * f0 * widths_eV))
    tail_rate = float(
        np.sum(total_frequency[tail_mask] * f0[tail_mask] * widths_eV[tail_mask])
    )
    tail_rate_fraction = tail_rate / total_rate if total_rate > 0.0 else 0.0
    f0_l1 = max(float(np.sum(np.abs(f0) * widths_eV)), 1.0e-300)
    f0_l2 = max(float(np.sum(f0 * f0 * widths_eV)), 1.0e-300)
    highest = coeff[-1]
    highest_l1 = float(np.sum(np.abs(highest) * widths_eV)) / f0_l1
    highest_l2 = float(np.sqrt(np.sum(highest * highest * widths_eV) / f0_l2))
    tolerance = float(case.config.multiterm_boltzmann.lmax_convergence_tolerance)
    warnings: list[str] = []
    if highest_l1 > tolerance:
        warnings.append("operator_lmax_may_be_underresolved")
    return (
        {
            "operator_tail_probability": tail_probability,
            "operator_tail_rate_fraction": float(tail_rate_fraction),
            "operator_highest_l_relative_l1": highest_l1,
            "operator_highest_l_relative_l2": highest_l2,
            "operator_lmax_convergence_tolerance": tolerance,
            "operator_lmax_convergence_ok": highest_l1 <= tolerance,
        },
        tuple(warnings),
    )


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
    """Run the validation-gated lmax>1 integral-cross-section closure.

    True lmax>1 transport needs higher-order differential scattering moments.
    With integral/effective cross sections only, the stable and honest path is
    to anchor f0, rates, and flux transport to the validated two-term
    Scharfetter-Gummel solve, then emit bounded higher-Legendre coefficients as
    diagnostics for future differential-collision work.
    """

    reference = solve_lmax1_two_term_reference(
        case.config, case.cross_sections, case.e_over_n_Td, case_id
    )
    energy = np.asarray(reference.energy_eV, dtype=float)
    widths = _widths_from_centers(energy)
    f0 = np.asarray(reference.eedf, dtype=float)
    f0 = f0 / max(float(np.sum(f0 * widths)), 1.0e-300)
    coeff = _anchored_lmax_coefficients(
        case, energy, widths, f0, reference.drift_velocity_m_s
    )
    rates = _convert_reference_rates(
        reference, solver_name, case_id, case.gas_number_density_m3
    )
    flux = FluxTransport.with_characteristic_energies(
        reference.drift_velocity_m_s,
        reference.mobility_m2_V_s,
        diffusion_longitudinal_m2_s=reference.diffusion_L_m2_s,
        diffusion_transverse_m2_s=reference.diffusion_T_m2_s,
    )
    transport = TransportSet.from_flux_only(
        flux,
        rates.ionization_frequency_s_inv,
        rates.attachment_frequency_s_inv,
        TransportMetadata(
            solver=solver_name,
            coefficient_definition="flux",
            swarm_condition="local_flux",
            notes=(
                "operator_lmax_gt1_reference_anchored_integral_closure",
            ),
        ),
    )
    method_used = "operator_reference_anchored_lmax_gt1"
    tail_mask = energy >= 0.9 * float(np.max(energy))
    tail = float(np.sum(f0[tail_mask] * widths[tail_mask]))
    quality_metadata, quality_warnings = _anchored_quality_metadata(
        case, energy, widths, coeff, f0
    )
    warnings = [
        "operator_lmax_gt1_reference_anchored_integral_closure",
        "operator_l_gt_0_elastic_integral_momentum_extension",
        "operator_l_gt_0_inelastic_source_isotropic_l0_only",
    ]
    if case.config.multiterm_boltzmann.hydrodynamic:
        warnings.append("operator_hydrodynamic_requested_but_not_computed")
    if not np.isclose(
        float(case.config.multiterm_boltzmann.field_coupling_scale),
        1.0,
        rtol=0.0,
        atol=1.0e-12,
    ):
        warnings.append("operator_nonphysical_field_scaling")
    warnings.extend(
        validity_warnings(
            case,
            project_collision_data(case),
            np.interp(
                case.grid.centers_eV,
                energy,
                f0,
                left=float(f0[0]),
                right=0.0,
            ),
        )
    )
    warnings.extend(quality_warnings)
    diagnostics = SolverDiagnostics(
        warnings=tuple(dict.fromkeys(warnings)),
        mean_energy_eV=float(reference.mean_energy_eV),
        eedf_tail_fraction=tail,
        power_balance_relative_residual=float(
            reference.metadata.get("residual_L1", np.nan)
        ),
    )
    metadata = {
        "operator_residual_L1": reference.metadata.get("residual_L1", np.nan),
        "operator_growth_frequency_s-1": reference.metadata.get(
            "growth_frequency_s-1", np.nan
        ),
        "operator_iterations": reference.metadata.get("iterations", 0),
        "operator_converged": reference.metadata.get("converged", False),
        "operator_lmax": int(case.config.multiterm_boltzmann.lmax),
        "operator_coefficient_order": "legendre_major_energy_minor",
        "operator_warnings": "; ".join(tuple(dict.fromkeys(warnings))),
        "operator_matrix_nnz": 0,
        "operator_base_matrix_nnz": 0,
        "operator_field_coupling_matrix_nnz": 0,
        "operator_streaming_matrix_nnz": 0,
        "operator_normalization": "integral_f0_dE_equals_1",
        "operator_ionization_source_model": _ionization_source_model(case.config),
        "operator_l_gt_0_elastic_model": "reference_anchored_integral_momentum_closure",
        "operator_l_gt_1_model": "bounded_field_relaxation_diagnostic",
        "operator_l_gt_0_inelastic_model": "sink_only_isotropic_l0_source",
        "operator_reference_gate_status": "anchored_to_native_two_term",
        "operator_reference_relerr_max": 0.0,
        "operator_transport_reused_from": "native_two_term_reference",
        "operator_rates_reused_from": "native_two_term_reference",
        "operator_hydrodynamic_computed": False,
    }
    metadata.update(quality_metadata)
    return MultiTermSolution(
        energy,
        widths,
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
    future l>1 block system must reduce to this reference before it can replace
    the current reference-anchored closure.
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
