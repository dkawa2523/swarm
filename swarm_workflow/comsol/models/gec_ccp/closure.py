"""Shared closure decisions and input contracts for GEC-CCP workflows."""

from __future__ import annotations

from typing import Any

import numpy as np

from swarm_workflow.comsol.models.gec_ccp.contracts import (
    FUNCTION_EEDF_PREINTEGRATED_INELASTIC,
    GEC_EXTERNAL_ELASTIC_ENERGY_LOSS,
    GEC_LOOKUP_INTERPOLATION,
    GEC_RESTRICTED_TRANSPORT_CLOSURE,
    GEC_TWO_TERM_HYBRID_TRANSPORT_CLOSURE,
    GecCcpMapping,
    GecCcpWorkflowError,
    GecClosureSpec,
    GecReactionSpec,
)
from swarm_workflow.comsol.models.gec_ccp.data import _required_float


GEC_SCALAR_DIFFUSION_RELATIVE_TOLERANCE = 1.0e-10
GEC_TRANSPORT_PROPERTIES = (
    (
        "sw_muN_e",
        "muN",
        "reduced_mobility_m2_V_s_m3",
        "1/(V*m*s)",
        "ptp.muerr",
    ),
    (
        "sw_DeN_e",
        "DeN",
        "reduced_diffusion_L_m2_s_m3",
        "1/(m*s)",
        "ptp.Derr",
    ),
    (
        "sw_muenN_e",
        "muenN",
        "reduced_electron_energy_mobility_m2_V_s_m3",
        "1/(V*m*s)",
        "ptp.muenrr",
    ),
    (
        "sw_DenN_e",
        "DenN",
        "reduced_electron_energy_diffusion_m2_s_m3",
        "1/(m*s)",
        "ptp.Denrr",
    ),
)
GEC_RESTRICTED_TRANSPORT_FUNCTIONS = (
    (
        "sw_muN_e",
        "muN",
        "reduced_mobility_m2_V_s_m3",
        "1/(V*m*s)",
    ),
    (
        "sw_DeN_L_e",
        "DeN_L",
        "reduced_diffusion_L_m2_s_m3",
        "1/(m*s)",
    ),
    (
        "sw_DeN_T_e",
        "DeN_T",
        "reduced_diffusion_T_m2_s_m3",
        "1/(m*s)",
    ),
    (
        "sw_muenN_e",
        "muenN",
        "reduced_electron_energy_mobility_m2_V_s_m3",
        "1/(V*m*s)",
    ),
    (
        "sw_DenN_L_e",
        "DenN_L",
        "reduced_electron_energy_diffusion_L_m2_s_m3",
        "1/(m*s)",
    ),
    (
        "sw_DenN_T_e",
        "DenN_T",
        "reduced_electron_energy_diffusion_T_m2_s_m3",
        "1/(m*s)",
    ),
)
GEC_TWO_TERM_HYBRID_TRANSPORT_FUNCTIONS = (
    (
        "sw_muN_e",
        "muN",
        "reduced_mobility_m2_V_s_m3",
        "1/(V*m*s)",
    ),
    (
        "sw_muenN_e",
        "muenN",
        "reduced_electron_energy_mobility_m2_V_s_m3",
        "1/(V*m*s)",
    ),
    (
        "sw_DenN_e",
        "DenN",
        "reduced_electron_energy_diffusion_m2_s_m3",
        "1/(m*s)",
    ),
)
GEC_MONTE_CARLO_HYBRID_TRANSPORT_FUNCTIONS = (
    GEC_TWO_TERM_HYBRID_TRANSPORT_FUNCTIONS[0],
    GEC_TWO_TERM_HYBRID_TRANSPORT_FUNCTIONS[1],
    (
        "sw_DenN_L_e",
        "DenN_L",
        "reduced_electron_energy_diffusion_L_m2_s_m3",
        "1/(m*s)",
    ),
    (
        "sw_DenN_T_e",
        "DenN_T",
        "reduced_electron_energy_diffusion_T_m2_s_m3",
        "1/(m*s)",
    ),
)
GEC_TWO_TERM_RESTRICTED_TRANSPORT_FUNCTIONS = tuple(
    item[:4] for item in GEC_TRANSPORT_PROPERTIES
)
# COMSOL flattens axisymmetric tensors column-major in the r, phi, z basis.
# The field-aligned tensors are symmetric, so the two r-z entries are equal.
GEC_TENSOR_COMPONENTS = (
    "rr", "phir", "zr", "rphi", "phiphi", "zphi", "rz", "phiz", "zz",
)
# Electron density, energy density, and electrostatic potential are scalar
# axisymmetric fields.  Their gradients have only r and z components, and the
# divergence of their fluxes contains only the r and z flux components.  COMSOL
# still stores the input properties as 3-by-3 tensors, but only this r-z block
# is equation-active in the present 2D-axisymmetric model.
GEC_AXISYMMETRIC_ACTIVE_TRANSPORT_COMPONENTS = ("rr", "zr", "rz", "zz")
GEC_AXISYMMETRIC_ACTIVE_TRANSPORT_INDICES = tuple(
    GEC_TENSOR_COMPONENTS.index(component)
    for component in GEC_AXISYMMETRIC_ACTIVE_TRANSPORT_COMPONENTS
)
GEC_ACTUAL_TRANSPORT_COMPONENTS = {
    "muN": tuple(f"ptp.mue{item}" for item in GEC_TENSOR_COMPONENTS),
    "DeN": tuple(f"ptp.De{item}" for item in GEC_TENSOR_COMPONENTS),
    "muenN": tuple(f"ptp.muen{item}" for item in GEC_TENSOR_COMPONENTS),
    "DenN": tuple(f"ptp.Den{item}" for item in GEC_TENSOR_COMPONENTS),
}


def _uses_transport(closure: GecClosureSpec) -> bool:
    return closure.electron_transport != "comsol"


def _uses_hybrid_einstein_transport(closure: GecClosureSpec) -> bool:
    return (
        closure.electron_transport
        == GEC_TWO_TERM_HYBRID_TRANSPORT_CLOSURE
    )


def _uses_external_elastic_energy_loss(closure: GecClosureSpec) -> bool:
    return (
        closure.elastic_energy_loss_model
        == GEC_EXTERNAL_ELASTIC_ENERGY_LOSS
    )


def _thermal_diffusion_enabled(closure: GecClosureSpec) -> bool:
    return closure.thermal_diffusion_model == "comsol_grad_diffusivity"


def _scalar_diffusion_audit(
    rows: list[dict[str, Any]],
    *,
    required_for: str,
) -> dict[str, Any]:
    pairs = {
        "particle_diffusion": (
            "reduced_diffusion_L_m2_s_m3",
            "reduced_diffusion_T_m2_s_m3",
        ),
        "energy_diffusion": (
            "reduced_electron_energy_diffusion_L_m2_s_m3",
            "reduced_electron_energy_diffusion_T_m2_s_m3",
        ),
    }
    results: dict[str, Any] = {}
    passed = True
    for name, (longitudinal_column, transverse_column) in pairs.items():
        longitudinal = np.asarray(
            [_required_float(row, longitudinal_column) for row in rows],
            dtype=float,
        )
        transverse = np.asarray(
            [_required_float(row, transverse_column) for row in rows],
            dtype=float,
        )
        relative_error = np.abs(longitudinal - transverse) / np.maximum(
            np.maximum(np.abs(longitudinal), np.abs(transverse)),
            1.0e-300,
        )
        pair_passed = bool(
            np.all(np.isfinite(relative_error))
            and np.all(longitudinal > 0.0)
            and np.all(transverse > 0.0)
            and np.max(relative_error)
            <= GEC_SCALAR_DIFFUSION_RELATIVE_TOLERANCE
        )
        passed = passed and pair_passed
        results[name] = {
            "passed": pair_passed,
            "longitudinal_column": longitudinal_column,
            "transverse_column": transverse_column,
            "maximum_relative_difference": float(np.max(relative_error)),
            "role": (
                "active_scalar_energy_diffusion_preflight"
                if name == "energy_diffusion"
                else "external_particle_diffusion_evidence_only_preflight"
            ),
        }
    scalar_energy = np.asarray(
        [
            _required_float(
                row, "reduced_electron_energy_diffusion_m2_s_m3"
            )
            for row in rows
        ],
        dtype=float,
    )
    energy_longitudinal = np.asarray(
        [
            _required_float(
                row,
                "reduced_electron_energy_diffusion_L_m2_s_m3",
            )
            for row in rows
        ],
        dtype=float,
    )
    energy_transverse = np.asarray(
        [
            _required_float(
                row,
                "reduced_electron_energy_diffusion_T_m2_s_m3",
            )
            for row in rows
        ],
        dtype=float,
    )
    scalar_relative_error = np.maximum(
        np.abs(scalar_energy - energy_longitudinal),
        np.abs(scalar_energy - energy_transverse),
    ) / np.maximum.reduce(
        (
            np.abs(scalar_energy),
            np.abs(energy_longitudinal),
            np.abs(energy_transverse),
            np.full_like(scalar_energy, 1.0e-300),
        )
    )
    scalar_binding_passed = bool(
        np.all(np.isfinite(scalar_relative_error))
        and np.all(scalar_energy > 0.0)
        and np.all(energy_longitudinal > 0.0)
        and np.all(energy_transverse > 0.0)
        and np.max(scalar_relative_error)
        <= GEC_SCALAR_DIFFUSION_RELATIVE_TOLERANCE
    )
    passed = passed and scalar_binding_passed
    results["active_energy_scalar_vs_LT"] = {
        "passed": scalar_binding_passed,
        "scalar_column": "reduced_electron_energy_diffusion_m2_s_m3",
        "longitudinal_column": (
            "reduced_electron_energy_diffusion_L_m2_s_m3"
        ),
        "transverse_column": (
            "reduced_electron_energy_diffusion_T_m2_s_m3"
        ),
        "maximum_relative_difference": float(
            np.max(scalar_relative_error)
        ),
        "role": "active_scalar_energy_diffusion_binding_preflight",
    }
    return {
        "passed": passed,
        "required_for": required_for,
        "relative_tolerance": GEC_SCALAR_DIFFUSION_RELATIVE_TOLERANCE,
        "pairs": results,
    }


def _validate_restricted_closure_policy(
    mapping: GecCcpMapping,
    *,
    source: str,
    transport_rows: list[dict[str, Any]] | None = None,
) -> dict[str, Any] | None:
    closure = mapping.closure
    if _uses_hybrid_einstein_transport(closure):
        if source not in {"two_term", "monte_carlo"}:
            raise GecCcpWorkflowError(
                "swarm_hybrid_einstein_de is qualified only for two_term "
                "or monte_carlo input"
            )
        if _thermal_diffusion_enabled(closure):
            raise GecCcpWorkflowError(
                "swarm_hybrid_einstein_de requires "
                "thermal_diffusion_model=off_restricted_diagonal"
            )
        if source == "two_term":
            if transport_rows is None:
                raise GecCcpWorkflowError(
                    "two-term hybrid transport requires scalar diffusion "
                    "verification"
                )
            scalar_audit = _scalar_diffusion_audit(
                transport_rows,
                required_for=(
                    "two_term_hybrid_active_energy_scalar_and_particle_"
                    "diffusion_evidence"
                ),
            )
            if not scalar_audit["passed"]:
                raise GecCcpWorkflowError(
                    "two-term hybrid transport requires its active scalar "
                    "energy diffusion to agree with L/T evidence"
                )
            return scalar_audit
        return None
    if closure.electron_transport != GEC_RESTRICTED_TRANSPORT_CLOSURE:
        if _thermal_diffusion_enabled(closure):
            raise GecCcpWorkflowError(
                "comsol_grad_diffusivity requires "
                "electron_transport=comsol_specify_all_restricted"
            )
        return None
    if source not in {"two_term", "monte_carlo"}:
        raise GecCcpWorkflowError(
            "comsol_specify_all_restricted is qualified only for two_term "
            "or monte_carlo input"
        )
    if _thermal_diffusion_enabled(closure) and source != "two_term":
        raise GecCcpWorkflowError(
            "comsol_grad_diffusivity is eligible only for the two_term "
            "restricted local-f0 closure"
        )
    scalar_audit: dict[str, Any] | None = None
    requires_two_term_scalar_binding = source == "two_term"
    if _thermal_diffusion_enabled(closure) or requires_two_term_scalar_binding:
        if transport_rows is None:
            raise GecCcpWorkflowError(
                "two-term restricted transport requires transport-table "
                "scalar diffusion verification"
            )
        scalar_audit = _scalar_diffusion_audit(
            transport_rows,
            required_for=(
                "two_term_isotropic_standard_local_energy_binding"
                if requires_two_term_scalar_binding
                else "comsol_grad_diffusivity"
            ),
        )
        if not scalar_audit["passed"]:
            raise GecCcpWorkflowError(
                "two-term restricted transport requires longitudinal and "
                "transverse particle/energy diffusion coefficients to agree"
            )
    return scalar_audit


def _restricted_gradient_response_metadata(
    mapping: GecCcpMapping,
    *,
    scalar_diffusion_audit: dict[str, Any] | None,
) -> dict[str, Any] | None:
    closure = mapping.closure
    if closure.electron_transport == "swarm_mobility_einstein":
        return {
            "model": "swarm_mobility_einstein",
            "closure_scope": (
                "external_swarm_mobility_with_COMSOL_SpecifyMueOnly_"
                "local_mean_energy_closure"
            ),
            "gradient_response_policy": closure.gradient_response_policy,
            "standard_local_energy_selected": True,
            "full_gradient_response_identified": False,
            "independent_2x2_gradient_response_identified": False,
            "comsol_executable_eligible": True,
            "restricted_closure_eligible": True,
            "full_physical_eligible": False,
            "thermal_diffusion_model": closure.thermal_diffusion_model,
            "transport_ownership": {
                "electron_mobility": "external_swarm_muN",
                "particle_diffusion": (
                    "COMSOL_mean_energy_dependent_Einstein_from_external_muN"
                ),
                "electron_energy_mobility": "COMSOL_SpecifyMueOnly_closure",
                "electron_energy_diffusion": "COMSOL_SpecifyMueOnly_closure",
            },
            "unconsumed_swarm_transport": [
                "particle_diffusion_L_T",
                "electron_energy_mobility",
                "electron_energy_diffusion_L_T",
            ],
            "limitation": (
                "external_swarm_data_do_not_identify_the_independent_"
                "density_and_mean_energy_gradient_response_matrix"
            ),
            "scalar_diffusion_audit": scalar_diffusion_audit,
        }
    if _uses_hybrid_einstein_transport(closure):
        source = mapping.bundle.expected_source
        return {
            "model": GEC_TWO_TERM_HYBRID_TRANSPORT_CLOSURE,
            "closure_scope": (
                "standard_COMSOL_electron_density_and_energy_density_flux"
            ),
            "gradient_response_policy": closure.gradient_response_policy,
            "standard_local_energy_selected": True,
            "full_gradient_response_identified": False,
            "independent_2x2_gradient_response_identified": False,
            "full_physical_eligible": False,
            "thermal_diffusion_model": closure.thermal_diffusion_model,
            "response_matrix_assumptions": {
                "particle_diffusion": "De=mu*Te_COMSOL_Einstein",
                "particle_cross_response": "0_restricted_assumption",
                "energy_diffusion": (
                    "external_MC_Den_L_T_field_aligned"
                    if source == "monte_carlo"
                    else "external_two_term_scalar_Den"
                ),
                "energy_cross_response": (
                    "same_Den_applied_by_standard_COMSOL_local_energy_closure"
                ),
            },
            "external_particle_diffusion": "evidence_only_not_consumed",
            "scalar_diffusion_audit": scalar_diffusion_audit,
        }
    if closure.electron_transport != GEC_RESTRICTED_TRANSPORT_CLOSURE:
        return None
    source = mapping.bundle.expected_source
    mc_density_packet = source == "monte_carlo"
    executable = source in {"two_term", "monte_carlo"}
    standard_selected = (
        closure.gradient_response_policy == "standard_local_energy"
    )
    particle_cross_response = (
        "dDe/dln(mean_energy)"
        if _thermal_diffusion_enabled(closure)
        else "0"
    )
    return {
        "model": GEC_RESTRICTED_TRANSPORT_CLOSURE,
        "closure_scope": (
            "standard_COMSOL_electron_density_and_energy_density_flux"
        ),
        "thermal_diffusion_model": closure.thermal_diffusion_model,
        "gradient_response_policy": closure.gradient_response_policy,
        "include_thermal_diffusion": _thermal_diffusion_enabled(closure),
        "full_gradient_response_identified": False,
        "independent_2x2_gradient_response_identified": False,
        "comsol_executable_eligible": executable,
        "restricted_closure_eligible": executable,
        "full_physical_eligible": False,
        "acceptance": {
            "executable": {
                "eligible": executable,
                "requires_cold_direct_solve_and_runtime_audits": True,
            },
            "standard_local_energy": {
                "selected": standard_selected,
                "eligible": executable and standard_selected,
                "model": (
                    "COMSOL_four_property_density_and_energy_density_flux"
                ),
            },
            "full_gradient_response": {
                "selected": not standard_selected,
                "eligible": False,
                "reason": "independent_2x2_gradient_response_not_identified",
            },
        },
        "response_matrix_normalization": (
            "Gamma_n/n=-mu_n*E-A_nn*grad(ln(n))-A_nE*"
            "grad(ln(mean_energy)); Gamma_u/u=-mu_u*E-A_un*"
            "grad(ln(n))-A_uE*grad(ln(mean_energy))"
        ),
        "response_matrix_assumptions": {
            "A_nn": "De",
            "A_nE": particle_cross_response,
            "A_un": "Den",
            "A_uE": "Den",
        },
        "response_matrix_provenance": {
            "A_nn": (
                "Monte_Carlo_density_packet_particle_diffusion"
                if mc_density_packet
                else "two_term_local_f0_particle_diffusion_moment"
            ),
            "A_nE": (
                "COMSOL_derivative_of_supplied_De"
                if _thermal_diffusion_enabled(closure)
                else "restricted_zero_assumption"
            ),
            "A_un": (
                "Monte_Carlo_density_packet_energy_diffusion_for_"
                "COMSOL_energy_density_flux"
                if mc_density_packet
                else "two_term_local_f0_energy_diffusion_moment"
            ),
            "A_uE": (
                "same_MC_Den_applied_to_grad_log_energy_density_by_"
                "standard_COMSOL_closure"
                if mc_density_packet
                else "two_term_local_f0_restricted_energy_flux_closure"
            ),
        },
        "source_interpretation": (
            "Monte_Carlo_direct_particle_and_energy_density_packet_"
            "coefficients_under_standard_COMSOL_restricted_closure;_"
            "not_an_independent_2x2_gradient_response"
            if mc_density_packet
            else "two_term_local_f0_restricted_moment_closure"
        ),
        "scalar_diffusion_audit": scalar_diffusion_audit,
    }


def _uses_external_rates(closure: GecClosureSpec) -> bool:
    return closure.reaction_model in {
        "external_rates",
        FUNCTION_EEDF_PREINTEGRATED_INELASTIC,
    }


def _uses_function_eedf(closure: GecClosureSpec) -> bool:
    return closure.reaction_model in {
        "function_eedf",
        FUNCTION_EEDF_PREINTEGRATED_INELASTIC,
    }


def _reaction_uses_preintegrated_rate(
    closure: GecClosureSpec,
    reaction: GecReactionSpec,
) -> bool:
    if (
        reaction.process_type == "elastic"
        and _uses_external_elastic_energy_loss(closure)
    ):
        return False
    return (
        closure.reaction_model == "external_rates"
        and reaction.process_type in closure.external_rate_processes
    ) or (
        closure.reaction_model == FUNCTION_EEDF_PREINTEGRATED_INELASTIC
        and reaction.process_type in {"excitation", "ionization"}
    )


def _function_eedf_rate_reactions(
    mapping: GecCcpMapping,
) -> tuple[tuple[int, GecReactionSpec], ...]:
    """Return rates derived from this EEDF, including its offline quadrature."""

    return tuple(
        (index, reaction)
        for index, reaction in enumerate(mapping.reactions, start=1)
        if (
            not _reaction_uses_preintegrated_rate(mapping.closure, reaction)
            or mapping.closure.reaction_model == FUNCTION_EEDF_PREINTEGRATED_INELASTIC
        )
        and not (
            reaction.process_type == "elastic"
            and _uses_external_elastic_energy_loss(mapping.closure)
        )
    )


def _reaction_handling_metadata(
    mapping: GecCcpMapping,
) -> list[dict[str, str]]:
    function_tag = (
        mapping.closure.function_eedf.function_tag
        if mapping.closure.function_eedf is not None
        else "physics_interface_eedf"
    )
    handling: list[dict[str, str]] = []
    for reaction in mapping.reactions:
        preintegrated = _reaction_uses_preintegrated_rate(
            mapping.closure, reaction
        )
        dense_function_rate = (
            mapping.closure.reaction_model
            == FUNCTION_EEDF_PREINTEGRATED_INELASTIC
            and preintegrated
        )
        if dense_function_rate:
            rate_source = "active_function_eedf_dense_reintegration"
            eedf_source = "same_swarm_eedf_offline_preintegration"
        elif (
            reaction.process_type == "elastic"
            and _uses_external_elastic_energy_loss(mapping.closure)
        ):
            rate_source = "disabled_COMSOL_elastic_reaction_rate"
            eedf_source = "same_solver_native_elastic_energy_loss_estimator"
        elif preintegrated and mapping.bundle.expected_source == "monte_carlo":
            rate_source = "monte_carlo_trajectory_time_average_sigma_v"
            eedf_source = (
                "independent_monte_carlo_eedf_consistency_evidence_not_rate_owner"
            )
        elif preintegrated and mapping.bundle.expected_source == "two_term":
            rate_source = "two_term_eedf_integrated_rate_table"
            eedf_source = "same_two_term_eedf_rate_integration"
        else:
            rate_source = "comsol_cross_section_integral"
            eedf_source = function_tag
        handling.append(
            {
                "feature": reaction.feature,
                "process_type": reaction.process_type,
                "binding": (
                    "DisabledElasticReaction"
                    if reaction.process_type == "elastic"
                    and _uses_external_elastic_energy_loss(mapping.closure)
                    else (
                        "RateConstant"
                        if preintegrated
                        else "UseCrossSectionData"
                    )
                ),
                "rate_source": rate_source,
                "eedf_source": eedf_source,
                "energy_loss": (
                    "external_solver_native_GeneralPowerDeposition"
                    if reaction.process_type == "elastic"
                    and _uses_external_elastic_energy_loss(mapping.closure)
                    else "preserved_existing_ElectronImpactReaction_de"
                ),
            }
        )
    return handling


def _positive_rate_suffix(
    rows: list[dict[str, Any]],
    *,
    process_type: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Select the final contiguous positive support without adding a floor."""

    selected = [
        row for row in rows if row.get("process_type") == process_type
    ]
    ordered = sorted(
        selected,
        key=lambda row: _required_float(row, "mean_energy_eV"),
    )
    if not ordered:
        raise GecCcpWorkflowError(
            f"missing rate rows for process {process_type}"
        )
    rates = np.asarray(
        [_required_float(row, "rate_coefficient_m3_s") for row in ordered],
        dtype=float,
    )
    if np.any(rates < 0.0):
        raise GecCcpWorkflowError(
            f"negative rate is not valid for log-PCHIP: {process_type}"
        )
    nonpositive = np.flatnonzero(rates <= 0.0)
    start = int(nonpositive[-1] + 1) if nonpositive.size else 0
    active = ordered[start:]
    active_energies = sorted(
        {
            _required_float(row, "mean_energy_eV")
            for row in active
        }
    )
    if len(active_energies) < 2:
        raise GecCcpWorkflowError(
            f"rate {process_type} needs at least two points in its final "
            "contiguous strictly positive support"
        )
    return active, {
        "process_type": process_type,
        "raw_rows": len(ordered),
        "zero_rows_observed": int(np.count_nonzero(rates == 0.0)),
        "rows_omitted_before_positive_suffix": start,
        "positive_rows_discarded_before_support": int(
            np.count_nonzero(rates[:start] > 0.0)
        ),
        "active_rows": len(active),
        "active_unique_mean_energy_points": len(active_energies),
        "positive_support_mean_energy_eV": [
            active_energies[0], active_energies[-1]
        ],
        "zero_policy": (
            "preserve_raw_zero_or_sampling_censored_observations_as_evidence;_"
            "bind_only_the_final_"
            "contiguous_strictly_positive_suffix; no_artificial_floor"
        ),
    }


def _mc_nonnegative_rate_table(
    rows: list[dict[str, Any]],
    *,
    process_type: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Retain every direct-MC rate anchor, including censored zero estimates."""

    ordered = sorted(
        (
            row
            for row in rows
            if row.get("process_type") == process_type
        ),
        key=lambda row: _required_float(row, "mean_energy_eV"),
    )
    if not ordered:
        raise GecCcpWorkflowError(
            f"missing rate rows for process {process_type}"
        )
    energies = np.asarray(
        [_required_float(row, "mean_energy_eV") for row in ordered],
        dtype=float,
    )
    rates = np.asarray(
        [_required_float(row, "rate_coefficient_m3_s") for row in ordered],
        dtype=float,
    )
    if (
        np.any(~np.isfinite(energies))
        or np.any(energies <= 0.0)
        or np.any(np.diff(energies) <= 0.0)
        or np.any(~np.isfinite(rates))
        or np.any(rates < 0.0)
    ):
        raise GecCcpWorkflowError(
            "Monte Carlo rate binding requires finite, strictly increasing "
            f"mean energy and finite nonnegative rates: {process_type}"
        )
    if len(ordered) < 2 or not np.any(rates > 0.0):
        raise GecCcpWorkflowError(
            "Monte Carlo rate binding needs at least two anchors and one "
            f"positive estimate: {process_type}"
        )
    scale = float(np.max(rates))
    return ordered, {
        "process_type": process_type,
        "raw_rows": len(ordered),
        "active_rows": len(ordered),
        "active_unique_mean_energy_points": len(ordered),
        "zero_rows_observed": int(np.count_nonzero(rates == 0.0)),
        "rows_omitted_before_positive_suffix": 0,
        "positive_rows_discarded_before_support": 0,
        "active_support_mean_energy_eV": [
            float(energies[0]), float(energies[-1])
        ],
        "normalization_rate_m3_s": scale,
        "interpolation": (
            "raw_nonnegative_rate_over_process_max_PCHIP_on_log_mean_energy"
        ),
        "zero_policy": (
            "bind_direct_MC_zero_point_estimates_exactly; no_floor; "
            "censoring_uncertainty_is_a_separate_physics_qualification"
        ),
    }


def _active_rate_table(
    mapping: GecCcpMapping,
    rows: list[dict[str, Any]],
    *,
    process_type: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    if _uses_direct_mc_rates(mapping):
        return _mc_nonnegative_rate_table(
            rows, process_type=process_type
        )
    return _positive_rate_suffix(rows, process_type=process_type)


def _uses_direct_mc_rates(mapping: GecCcpMapping) -> bool:
    return (
        mapping.bundle.expected_source == "monte_carlo"
        and mapping.closure.reaction_model == "external_rates"
    )


def _active_rate_support_metadata(
    mapping: GecCcpMapping,
    rows: list[dict[str, Any]],
) -> dict[str, Any] | None:
    if not _uses_external_rates(mapping.closure):
        return None
    processes: dict[str, Any] = {}
    for reaction in mapping.reactions:
        if not _reaction_uses_preintegrated_rate(
            mapping.closure, reaction
        ):
            continue
        _, metadata = _active_rate_table(
            mapping, rows, process_type=reaction.process_type
        )
        processes[reaction.process_type] = metadata
    return {
        "interpolation": GEC_LOOKUP_INTERPOLATION,
        "strictly_positive_active_support": (
            not _uses_direct_mc_rates(mapping)
        ),
        "nonnegative_zero_preserving_active_support": (
            _uses_direct_mc_rates(mapping)
        ),
        "artificial_floor": False,
        "processes": processes,
    }


def _active_transport_functions(
    closure: GecClosureSpec,
    *,
    source: str | None = None,
) -> tuple[tuple[str, str, str, str], ...]:
    if closure.electron_transport == "swarm_mobility_einstein":
        tag, name, column, unit, _ = GEC_TRANSPORT_PROPERTIES[0]
        return ((tag, name, column, unit),)
    if _uses_hybrid_einstein_transport(closure):
        return (
            GEC_MONTE_CARLO_HYBRID_TRANSPORT_FUNCTIONS
            if source == "monte_carlo"
            else GEC_TWO_TERM_HYBRID_TRANSPORT_FUNCTIONS
        )
    if closure.electron_transport == GEC_RESTRICTED_TRANSPORT_CLOSURE:
        if source == "two_term":
            return GEC_TWO_TERM_RESTRICTED_TRANSPORT_FUNCTIONS
        return GEC_RESTRICTED_TRANSPORT_FUNCTIONS
    return ()


def _transport_input_contract(
    closure: GecClosureSpec,
    *,
    source: str | None = None,
) -> dict[str, dict[str, Any]]:
    """Describe every COMSOL electron-transport property independently."""

    built_in = {
        name: {
            "source": "comsol_local_energy_closure",
            "external_table_column": None,
            "formula": None,
        }
        for name in ("muN", "DeN", "muenN", "DenN")
    }
    if closure.electron_transport == "comsol":
        return built_in
    built_in["muN"] = {
        "source": "external_swarm_table",
        "external_table_column": "reduced_mobility_m2_V_s_m3",
        "formula": "muN(mean_energy)",
    }
    if closure.electron_transport == "swarm_mobility_einstein":
        return built_in
    if _uses_hybrid_einstein_transport(closure):
        energy_diffusion_column: str | list[str] = (
            [
                "reduced_electron_energy_diffusion_L_m2_s_m3",
                "reduced_electron_energy_diffusion_T_m2_s_m3",
            ]
            if source == "monte_carlo"
            else "reduced_electron_energy_diffusion_m2_s_m3"
        )
        built_in.update(
            {
                "DeN": {
                    "source": "comsol_einstein_from_external_mobility",
                    "external_table_column": None,
                    "formula": "DeN=muN*ptp.Te",
                },
                "muenN": {
                    "source": "external_swarm_table",
                    "external_table_column": (
                        "reduced_electron_energy_mobility_m2_V_s_m3"
                    ),
                    "formula": "muenN(mean_energy)",
                },
                "DenN": {
                    "source": (
                        "external_swarm_table_field_aligned"
                        if source == "monte_carlo"
                        else "external_swarm_table"
                    ),
                    "external_table_column": energy_diffusion_column,
                    "formula": (
                        "DenN_L/T(mean_energy)"
                        if source == "monte_carlo"
                        else "DenN(mean_energy)"
                    ),
                },
            }
        )
        return built_in
    if source == "two_term":
        built_in.update(
            {
                "DeN": {
                    "source": "external_swarm_table_isotropic",
                    "external_table_column": (
                        "reduced_diffusion_L_m2_s_m3"
                    ),
                    "formula": "DeN(mean_energy); L=T preflight required",
                },
                "muenN": {
                    "source": "external_swarm_table",
                    "external_table_column": (
                        "reduced_electron_energy_mobility_m2_V_s_m3"
                    ),
                    "formula": "muenN(mean_energy)",
                },
                "DenN": {
                    "source": "external_swarm_table_isotropic",
                    "external_table_column": (
                        "reduced_electron_energy_diffusion_m2_s_m3"
                    ),
                    "formula": "DenN(mean_energy); L=T preflight required",
                },
            }
        )
        return built_in
    built_in.update(
        {
            "DeN": {
                "source": "external_swarm_table_field_aligned",
                "external_table_column": [
                    "reduced_diffusion_L_m2_s_m3",
                    "reduced_diffusion_T_m2_s_m3",
                ],
                "formula": "DeN_L/T(mean_energy)",
            },
            "muenN": {
                "source": "external_swarm_table",
                "external_table_column": (
                    "reduced_electron_energy_mobility_m2_V_s_m3"
                ),
                "formula": "muenN(mean_energy)",
            },
            "DenN": {
                "source": "external_swarm_table_field_aligned",
                "external_table_column": [
                    "reduced_electron_energy_diffusion_L_m2_s_m3",
                    "reduced_electron_energy_diffusion_T_m2_s_m3",
                ],
                "formula": "DenN_L/T(mean_energy)",
            },
        }
    )
    return built_in


def _transport_tensor_metadata(
    mapping: GecCcpMapping,
    *,
    source: str,
) -> dict[str, Any] | None:
    closure = mapping.closure
    if closure.electron_transport == "swarm_mobility_einstein":
        return {
            "basis": "isotropic_unmagnetized_external_mobility",
            "scope": (
                "external_muN_only_with_COMSOL_SpecifyMueOnly_"
                "local_mean_energy_closure"
            ),
            "particle_diffusion": (
                "COMSOL_mean_energy_dependent_Einstein_from_external_muN"
            ),
            "energy_transport": "COMSOL_SpecifyMueOnly_closure",
            "thermal_diffusion": False,
            "full_eedf_consistent_transport": False,
            "full_gradient_response_identified": False,
            "external_particle_diffusion_consumed": False,
            "external_energy_transport_consumed": False,
            "comsol_magnetized_tensor_computation": False,
        }
    if _uses_hybrid_einstein_transport(closure):
        if source == "monte_carlo":
            return {
                "basis": "instantaneous_local_electric_field_r_phi_z",
                "scope": (
                    "external_MC_muN_muenN_DenN_L_T_with_COMSOL_"
                    "Einstein_DeN"
                ),
                "particle_diffusion": "DeN=muN*ptp.Te",
                "energy_diffusion": "external_MC_field_aligned_DenN_L_T",
                "thermal_diffusion": False,
                "full_eedf_consistent_transport": False,
                "full_gradient_response_identified": False,
                "external_particle_diffusion_consumed": False,
                "external_MC_particle_diffusion_role": "evidence_only",
                "zero_field_limit": "trace_preserving_isotropic",
                "zero_field_isotropization_Td": (
                    closure.zero_field_isotropization_Td
                ),
                "comsol_magnetized_tensor_computation": False,
            }
        return {
            "basis": "isotropic_unmagnetized_two_term",
            "scope": (
                "external_muN_muenN_DenN_with_COMSOL_Einstein_DeN"
            ),
            "particle_diffusion": "DeN=muN*ptp.Te",
            "energy_diffusion": "external_two_term_scalar_DenN",
            "thermal_diffusion": False,
            "full_eedf_consistent_transport": False,
            "full_gradient_response_identified": False,
            "external_particle_diffusion_consumed": False,
            "comsol_magnetized_tensor_computation": False,
        }
    if closure.electron_transport != GEC_RESTRICTED_TRANSPORT_CLOSURE:
        return None
    if source == "two_term":
        return {
            "basis": "isotropic_unmagnetized_two_term",
            "particle_diffusion": "external_two_term_scalar_DeN",
            "energy_diffusion": "external_two_term_scalar_DenN",
            "longitudinal_transverse_preflight_required": True,
            "electric_field_orientation_used": False,
            "zero_field_isotropization_used": False,
            "thermal_diffusion": _thermal_diffusion_enabled(closure),
            "scope": (
                "four_COMSOL_SpecifyAll_transport_properties; "
                "standard_local_energy_restricted_not_full_gradient_response"
            ),
            "gradient_closure": GEC_RESTRICTED_TRANSPORT_CLOSURE,
            "full_gradient_response_identified": False,
            "energy_diffusivity_interpretation": (
                "two_term_scalar_f0_moment_closure"
            ),
            "thermal_diffusion_interpretation": (
                "COMSOL_expands_grad(De*ne),_giving_"
                "A_nE=dDe/dln(mean_energy)"
                if _thermal_diffusion_enabled(closure)
                else "A_nE=0_restricted_assumption;_smallness_not_established"
            ),
            "comsol_magnetized_tensor_computation": False,
        }
    return {
        "basis": "instantaneous_local_electric_field_r_phi_z",
        "strong_field_limit": "D_T I + (D_L-D_T) ee",
        "zero_field_limit": "trace_preserving_isotropic",
        "thermal_diffusion": _thermal_diffusion_enabled(closure),
        "scope": (
            "four_COMSOL_SpecifyAll_transport_properties; "
            "restricted_not_full_kinetic_gradient_response"
        ),
        "gradient_closure": GEC_RESTRICTED_TRANSPORT_CLOSURE,
        "full_gradient_response_identified": False,
        "energy_diffusivity_interpretation": (
            "direct_MC_density_packet_effective_coefficient_with_"
            "unmeasured_energy_gradient_column"
            if source == "monte_carlo"
            else "two_term_scalar_f0_moment_closure"
        ),
        "thermal_diffusion_interpretation": (
            "COMSOL_expands_grad(De*ne),_giving_"
            "A_nE=dDe/dln(mean_energy)"
            if _thermal_diffusion_enabled(closure)
            else "A_nE=0_restricted_assumption;_smallness_not_established"
        ),
        "comsol_magnetized_tensor_computation": False,
    }


def _transport_reduced_expression(
    function_tag: str,
    reduced_unit: str,
    *,
    argument: str = "sw_logeps",
) -> str:
    log_tag = f"sw_log_{function_tag.removeprefix('sw_')}"
    return f"exp({log_tag}({argument}))*1[{reduced_unit}]"


def _isotropic_tensor(expression: str, unit: str) -> tuple[str, ...]:
    zero = f"0[{unit}]"
    return tuple(
        expression if index in {0, 4, 8} else zero for index in range(9)
    )


def _field_aligned_tensor(
    longitudinal: str,
    transverse: str,
    *,
    zero_field_isotropization_Td: float,
) -> tuple[str, ...]:
    """Return a positive, trace-preserving r/phi/z diffusion tensor.

    The strong-field limit is ``D_T I + (D_L-D_T) ee``.  At an RF field
    reversal the direction ``e`` is undefined, so the traceless anisotropic
    part tends smoothly to zero and the tensor tends to
    ``(D_L+2 D_T) I/3``.  The transition scale is explicit mapping data, not
    a nonlinear-solver damping parameter.
    """

    delta = f"(({longitudinal})-({transverse}))"
    isotropic = f"(({longitudinal})+2*({transverse}))/3"
    field_squared = "(ptp.Er^2+ptp.Ez^2)"
    regularization = (
        f"({zero_field_isotropization_Td:.17g}*1e-21[V*m^2]*ptp.Nn)"
    )
    denominator = f"({field_squared}+({regularization})^2)"
    trace_part = f"{field_squared}/(3*{denominator})"
    rr = (
        f"({isotropic})+({delta})*(ptp.Er^2/{denominator}-{trace_part})"
    )
    phiphi = f"({isotropic})-({delta})*({trace_part})"
    zz = (
        f"({isotropic})+({delta})*(ptp.Ez^2/{denominator}-{trace_part})"
    )
    rz = f"({delta})*ptp.Er*ptp.Ez/{denominator}"
    zero = "0[1/(m*s)]"
    return (rr, zero, rz, zero, phiphi, zero, rz, zero, zz)


def _transport_property_tensors(
    closure: GecClosureSpec,
    *,
    source: str | None = None,
) -> tuple[tuple[str, str, tuple[str, ...], str], ...]:
    if closure.electron_transport == "comsol":
        return ()
    argument = (
        "sw_logeps_transport" if source == "monte_carlo" else "sw_logeps"
    )
    mobility = _transport_reduced_expression(
        "sw_muN_e", "1/(V*m*s)", argument=argument
    )
    bindings: list[tuple[str, str, tuple[str, ...], str]] = [
        (
            "muN",
            "muN",
            _isotropic_tensor(mobility, "1/(V*m*s)"),
            "m^2/(V*s)",
        )
    ]
    if closure.electron_transport == "swarm_mobility_einstein":
        return tuple(bindings)
    if _uses_hybrid_einstein_transport(closure):
        particle_diffusion = f"({mobility})*ptp.Te"
        energy_mobility = _transport_reduced_expression(
            "sw_muenN_e", "1/(V*m*s)", argument=argument
        )
        if source == "monte_carlo":
            scale = closure.zero_field_isotropization_Td
            if scale is None:  # Parser validation makes this unreachable.
                raise GecCcpWorkflowError("missing zero-field tensor scale")
            energy_diffusion = _field_aligned_tensor(
                _transport_reduced_expression(
                    "sw_DenN_L_e", "1/(m*s)", argument=argument
                ),
                _transport_reduced_expression(
                    "sw_DenN_T_e", "1/(m*s)", argument=argument
                ),
                zero_field_isotropization_Td=scale,
            )
        else:
            scalar_energy_diffusion = _transport_reduced_expression(
                "sw_DenN_e", "1/(m*s)", argument=argument
            )
            energy_diffusion = _isotropic_tensor(
                scalar_energy_diffusion, "1/(m*s)"
            )
        bindings.extend(
            (
                (
                    "DeN",
                    "DeN",
                    _isotropic_tensor(particle_diffusion, "1/(m*s)"),
                    "m^2/s",
                ),
                (
                    "muenN",
                    "muenN",
                    _isotropic_tensor(energy_mobility, "1/(V*m*s)"),
                    "m^2/(V*s)",
                ),
                (
                    "DenN",
                    "DenN",
                    energy_diffusion,
                    "m^2/s",
                ),
            )
        )
        return tuple(bindings)
    if source == "two_term":
        particle_diffusion = _transport_reduced_expression(
            "sw_DeN_e", "1/(m*s)", argument=argument
        )
        energy_mobility = _transport_reduced_expression(
            "sw_muenN_e", "1/(V*m*s)", argument=argument
        )
        energy_diffusion = _transport_reduced_expression(
            "sw_DenN_e", "1/(m*s)", argument=argument
        )
        bindings.extend(
            (
                (
                    "DeN",
                    "DeN",
                    _isotropic_tensor(particle_diffusion, "1/(m*s)"),
                    "m^2/s",
                ),
                (
                    "muenN",
                    "muenN",
                    _isotropic_tensor(energy_mobility, "1/(V*m*s)"),
                    "m^2/(V*s)",
                ),
                (
                    "DenN",
                    "DenN",
                    _isotropic_tensor(energy_diffusion, "1/(m*s)"),
                    "m^2/s",
                ),
            )
        )
        return tuple(bindings)
    scale = closure.zero_field_isotropization_Td
    if scale is None:  # Parser validation makes this unreachable.
        raise GecCcpWorkflowError("missing zero-field tensor scale")
    diffusion = _field_aligned_tensor(
        _transport_reduced_expression(
            "sw_DeN_L_e", "1/(m*s)", argument=argument
        ),
        _transport_reduced_expression(
            "sw_DeN_T_e", "1/(m*s)", argument=argument
        ),
        zero_field_isotropization_Td=scale,
    )
    energy_mobility = _transport_reduced_expression(
        "sw_muenN_e", "1/(V*m*s)", argument=argument
    )
    energy_diffusion = _field_aligned_tensor(
        _transport_reduced_expression(
            "sw_DenN_L_e", "1/(m*s)", argument=argument
        ),
        _transport_reduced_expression(
            "sw_DenN_T_e", "1/(m*s)", argument=argument
        ),
        zero_field_isotropization_Td=scale,
    )
    bindings.extend(
        (
            ("DeN", "DeN", diffusion, "m^2/s"),
            (
                "muenN",
                "muenN",
                _isotropic_tensor(energy_mobility, "1/(V*m*s)"),
                "m^2/(V*s)",
            ),
            ("DenN", "DenN", energy_diffusion, "m^2/s"),
        )
    )
    return tuple(bindings)


def _transport_audit_components(
    closure: GecClosureSpec,
    *,
    source: str | None = None,
) -> tuple[tuple[str, str, str, str, str], ...]:
    components: list[tuple[str, str, str, str, str]] = []
    for _, quantity, tensor, local_unit in _transport_property_tensors(
        closure, source=source
    ):
        actual = GEC_ACTUAL_TRANSPORT_COMPONENTS[quantity]
        for component, actual_expression, reduced_expression in zip(
            GEC_TENSOR_COMPONENTS, actual, tensor, strict=True
        ):
            if component not in GEC_AXISYMMETRIC_ACTIVE_TRANSPORT_COMPONENTS:
                continue
            expected = f"({reduced_expression})/(ptp.Nn)"
            components.append(
                (quantity, component, actual_expression, expected, local_unit)
            )
    return tuple(components)
