"""Canonical native result-table projections, shared by publication and analysis."""

from dataclasses import asdict
import pandas as pd
from electron_swarm.core.results import SwarmCaseResult, SwarmRunResult
from electron_swarm.core.result_metadata import SUMMARY_METADATA_KEYS
from electron_swarm.core.numerics import eepf_from_eedf, widths_from_centers

SUMMARY_COLUMNS = [
    "solver",
    "solver_method",
    "case_id",
    "E_over_N_Td",
    "mean_energy_eV",
    "drift_velocity_m_s",
    "mobility_m2_V_s",
    "reduced_mobility_m2_V_s_m3",
    "diffusion_L_m2_s",
    "diffusion_T_m2_s",
    "reduced_diffusion_L_m2_s_m3",
    "reduced_diffusion_T_m2_s_m3",
    "reduced_electron_energy_mobility_m2_V_s_m3",
    "reduced_electron_energy_diffusion_m2_s_m3",
    "reduced_electron_energy_diffusion_L_m2_s_m3",
    "reduced_electron_energy_diffusion_T_m2_s_m3",
    "net_ionization_frequency_s",
    "effective_townsend_m2",
]

SUMMARY_METADATA_COLUMNS = list(SUMMARY_METADATA_KEYS)

RATES_COLUMNS = [
    "solver",
    "case_id",
    "E_over_N_Td",
    "species",
    "process",
    "process_type",
    "threshold_eV",
    "rate_coefficient_m3_s",
    "target_species_fraction",
    "energy_loss_eV",
    "energy_loss_rate_coefficient_eV_m3_s",
    "mixture_weighted_rate_m3_s",
    "frequency_s_inv",
    "power_loss_eV_s",
    "tail_fraction",
]

EEDF_COLUMNS = [
    "solver",
    "case_id",
    "E_over_N_Td",
    "energy_eV",
    "energy_width_eV",
    "eedf",
    "eepf",
    "sample_count",
    "weight_effective_count",
]

ENERGY_ANGLE_COLUMNS = [
    "solver",
    "case_id",
    "E_over_N_Td",
    "energy_eV",
    "energy_width_eV",
    "mu",
    "mu_width",
    "energy_angle_density_eV_inv",
]

RF_PHASE_COLUMNS = [
    "solver",
    "case_id",
    "E_over_N_rms_Td",
    "phase_index",
    "phase_fraction",
    "phase_rad",
    "instantaneous_E_over_N_Td",
    "mean_energy_eV",
    "ionization_rate_coefficient_m3_s",
]

RF_PHASE_EEDF_COLUMNS = [
    "solver",
    "case_id",
    "E_over_N_rms_Td",
    "phase_index",
    "phase_fraction",
    "phase_rad",
    "instantaneous_E_over_N_Td",
    "mean_energy_eV",
    "energy_eV",
    "energy_width_eV",
    "eedf",
    "eepf",
]

SOLVER_PLAN_COLUMNS = [
    "solver",
    "state",
    "requested",
    "runnable",
    "skipped",
    "degraded",
    "skip_reason",
    "warnings",
    "angular_scattering_support",
    "angular_scattering_fidelity",
    "angular_scattering_assumption",
    "effective_angular_scattering",
    "effective_collision_processes",
    "effective_ionization_source",
    "effective_electron_electron",
    "effective_magnetic_field",
    "effective_rf_field",
    "effective_tail_refinement",
    "effective_finite_k",
]

COMPARISON_SUMMARY_COLUMNS = [
    "case_id",
    "E_over_N_Td",
    "reference_solver",
    "candidate_solver",
    "angular_model_status",
    "physics_status",
    "physics_differences",
    "comparison_status",
    "reference_numerical_status",
    "candidate_numerical_status",
    "reference_transport_definition",
    "candidate_transport_definition",
    "mean_energy_eV_relative_difference",
    "mean_energy_eV_status",
    "drift_velocity_relative_difference",
    "drift_velocity_m_s_status",
    "mobility_relative_difference",
    "mobility_m2_V_s_status",
    "diffusion_L_relative_difference",
    "diffusion_L_m2_s_status",
    "diffusion_T_relative_difference",
    "diffusion_T_m2_s_status",
    "net_ionization_frequency_relative_difference",
    "net_ionization_frequency_s_status",
    "eedf_l1_error",
]


def _summary_columns() -> list[str]:
    return SUMMARY_COLUMNS + [f"meta_{key}" for key in SUMMARY_METADATA_COLUMNS]


def _canonical_frame(
    rows: list[dict[str, object]],
    columns: list[str],
) -> pd.DataFrame:
    return pd.DataFrame(rows, dtype=object).reindex(columns=columns)


def _summary_frame(cases: list[SwarmCaseResult]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for case in cases:
        row = {
            "solver": case.solver,
            "case_id": case.case_id,
            "E_over_N_Td": case.e_over_n_Td,
            "mean_energy_eV": case.mean_energy_eV,
            "drift_velocity_m_s": case.drift_velocity_m_s,
            "mobility_m2_V_s": case.mobility_m2_V_s,
            "reduced_mobility_m2_V_s_m3": case.reduced_mobility_m2_V_s_m3,
            "diffusion_L_m2_s": case.diffusion_L_m2_s,
            "diffusion_T_m2_s": case.diffusion_T_m2_s,
            "reduced_diffusion_L_m2_s_m3": case.reduced_diffusion_L_m2_s_m3,
            "reduced_diffusion_T_m2_s_m3": case.reduced_diffusion_T_m2_s_m3,
            "reduced_electron_energy_mobility_m2_V_s_m3": (
                case.reduced_electron_energy_mobility_m2_V_s_m3
            ),
            "reduced_electron_energy_diffusion_m2_s_m3": (
                case.reduced_electron_energy_diffusion_m2_s_m3
            ),
            "reduced_electron_energy_diffusion_L_m2_s_m3": (
                case.reduced_electron_energy_diffusion_L_m2_s_m3
            ),
            "reduced_electron_energy_diffusion_T_m2_s_m3": (
                case.reduced_electron_energy_diffusion_T_m2_s_m3
            ),
            "net_ionization_frequency_s": case.net_ionization_frequency_s,
            "effective_townsend_m2": case.effective_townsend_m2,
        }
        row["solver_method"] = case.metadata.get("solver_method", "")
        for key in SUMMARY_METADATA_COLUMNS:
            row[f"meta_{key}"] = case.metadata.get(key)
        rows.append(row)
    return _canonical_frame(rows, _summary_columns())


def _eedf_frame(cases: list[SwarmCaseResult]) -> pd.DataFrame:
    rows = []
    for case in cases:
        counts = case.eedf_counts
        effective_counts = case.eedf_effective_counts
        widths = case.energy_widths_eV
        if widths is None or len(widths) != len(case.energy_eV):
            widths = widths_from_centers(case.energy_eV)
        eepf = eepf_from_eedf(case.energy_eV, case.eedf)
        for index, (energy, eedf, eepf_value) in enumerate(
            zip(case.energy_eV, case.eedf, eepf, strict=False)
        ):
            width = (
                float(widths[index])
                if widths is not None
                and index < len(widths)
                and pd.notna(widths[index])
                else None
            )
            count = (
                counts[index] if counts is not None and index < len(counts) else None
            )
            effective_count = (
                float(effective_counts[index])
                if effective_counts is not None and index < len(effective_counts)
                else None
            )
            rows.append(
                {
                    "solver": case.solver,
                    "case_id": case.case_id,
                    "E_over_N_Td": case.e_over_n_Td,
                    "energy_eV": energy,
                    "energy_width_eV": width,
                    "eedf": eedf,
                    "eepf": eepf_value,
                    "sample_count": count,
                    "weight_effective_count": effective_count,
                }
            )
    return _canonical_frame(rows, EEDF_COLUMNS)


def _rates_frame(cases: list[SwarmCaseResult]) -> pd.DataFrame:
    rows = []
    for case in cases:
        for rate in case.rates:
            rows.append(
                {
                    "solver": rate.solver,
                    "case_id": rate.case_id,
                    "E_over_N_Td": rate.e_over_n_Td,
                    "species": rate.species,
                    "process": rate.process,
                    "process_type": rate.process_type,
                    "threshold_eV": rate.threshold_eV,
                    "rate_coefficient_m3_s": rate.rate_coefficient_m3_s,
                    "target_species_fraction": rate.target_species_fraction,
                    "energy_loss_eV": rate.energy_loss_eV,
                    "energy_loss_rate_coefficient_eV_m3_s": (
                        rate.energy_loss_rate_coefficient_eV_m3_s
                    ),
                    "mixture_weighted_rate_m3_s": rate.mixture_weighted_rate_m3_s,
                    "frequency_s_inv": rate.frequency_s_inv,
                    "power_loss_eV_s": rate.power_loss_eV_s,
                    "tail_fraction": rate.tail_fraction,
                }
            )
    return _canonical_frame(rows, RATES_COLUMNS)


def _energy_angle_frame(cases: list[SwarmCaseResult]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for case in cases:
        distribution = case.energy_angle_distribution
        if distribution is None:
            continue
        expected_shape = (
            len(distribution.energy_eV),
            len(distribution.mu),
        )
        if distribution.density_eV_inv.shape != expected_shape:
            raise ValueError(
                "energy-angle density shape does not match its coordinates"
            )
        for energy_index, energy in enumerate(distribution.energy_eV):
            for polar_index, mu in enumerate(distribution.mu):
                rows.append(
                    {
                        "solver": case.solver,
                        "case_id": case.case_id,
                        "E_over_N_Td": case.e_over_n_Td,
                        "energy_eV": energy,
                        "energy_width_eV": (
                            distribution.energy_widths_eV[energy_index]
                        ),
                        "mu": mu,
                        "mu_width": distribution.mu_widths[polar_index],
                        "energy_angle_density_eV_inv": (
                            distribution.density_eV_inv[
                                energy_index,
                                polar_index,
                            ]
                        ),
                    }
                )
    return _canonical_frame(rows, ENERGY_ANGLE_COLUMNS)


def _rf_phase_frames(
    cases: list[SwarmCaseResult],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    phase_rows: list[dict[str, object]] = []
    eedf_rows: list[dict[str, object]] = []
    for case in cases:
        phase = case.rf_phase
        if phase is None:
            continue
        widths = case.energy_widths_eV
        if widths is None or len(widths) != len(case.energy_eV):
            widths = widths_from_centers(case.energy_eV)
        for phase_index in range(len(phase.phase_fraction)):
            base = {
                "solver": case.solver,
                "case_id": case.case_id,
                "E_over_N_rms_Td": case.e_over_n_Td,
                "phase_index": phase_index,
                "phase_fraction": phase.phase_fraction[phase_index],
                "phase_rad": phase.phase_rad[phase_index],
                "instantaneous_E_over_N_Td": (
                    phase.instantaneous_e_over_n_Td[phase_index]
                ),
                "mean_energy_eV": phase.mean_energy_eV[phase_index],
                "ionization_rate_coefficient_m3_s": (
                    phase.ionization_rate_coefficient_m3_s[phase_index]
                ),
            }
            phase_rows.append(base)
            eepf = eepf_from_eedf(case.energy_eV, phase.eedf[phase_index])
            for energy_index, energy in enumerate(case.energy_eV):
                eedf_rows.append(
                    {
                        **base,
                        "energy_eV": energy,
                        "energy_width_eV": widths[energy_index],
                        "eedf": phase.eedf[phase_index, energy_index],
                        "eepf": eepf[energy_index],
                    }
                )
    return (
        _canonical_frame(phase_rows, RF_PHASE_COLUMNS),
        _canonical_frame(eedf_rows, RF_PHASE_EEDF_COLUMNS),
    )


def _solver_plan_frame(result: SwarmRunResult) -> pd.DataFrame:
    rows = result.metadata.get("solver_plan", [])
    return _canonical_frame(rows, SOLVER_PLAN_COLUMNS)


_STRING_COLUMNS = {
    "solver",
    "solver_method",
    "case_id",
    "species",
    "process",
    "process_type",
    "reference_solver",
    "candidate_solver",
    "reference_transport_definition",
    "candidate_transport_definition",
    "physics_differences",
}
_BOOLEAN_COLUMNS = {
    "requested",
    "runnable",
    "skipped",
    "degraded",
    "meta_exact_dcs_based",
    "meta_ordinary_integral_xs_closure",
    "meta_electron_electron_transport_stale",
}
_INTEGER_COLUMNS = {
    "phase_index",
    "sample_count",
    "meta_lmax",
    "meta_monte_carlo_base_seed",
    "meta_monte_carlo_case_seed",
    "meta_ionization_daughter_count",
}


def table_kinds(key: str, columns) -> dict[str, str]:
    kinds = {}
    for name in columns:
        if name in _BOOLEAN_COLUMNS:
            kind = "boolean"
        elif name in _INTEGER_COLUMNS:
            kind = "integer"
        elif (
            name in _STRING_COLUMNS
            or name.endswith("_status")
            or key == "solver_plan_csv"
            or (name.startswith("meta_") and name != "meta_rf_frequency_Hz")
        ):
            kind = "string"
        else:
            kind = "number"
        kinds[name] = kind
    return kinds


def result_tables(result: SwarmRunResult) -> dict[str, pd.DataFrame]:
    tables = (
        {"solver_plan_csv": _solver_plan_frame(result)}
        if result.metadata.get("solver_plan")
        else {}
    )
    if result.cases:
        tables.update(
            {
                "summary_csv": _summary_frame(result.cases),
                "eedf_csv": _eedf_frame(result.cases),
                "rates_csv": _rates_frame(result.cases),
            }
        )
        angular = _energy_angle_frame(result.cases)
        if not angular.empty:
            tables["energy_angle_distribution_csv"] = angular
        phase, eedf = _rf_phase_frames(result.cases)
        if not phase.empty:
            tables.update({"rf_phase_csv": phase, "rf_phase_eedf_csv": eedf})
    if result.metadata.get("comparison", {}).get("status", "disabled") not in {
        "disabled",
        "not_evaluated",
    }:
        tables["comparison_summary_csv"] = pd.DataFrame(
            result.metadata.get("comparison_summary_rows", []),
            columns=COMPARISON_SUMMARY_COLUMNS,
            dtype=object,
        )
    return tables


def case_descriptions(cases):
    return [
        {
            "solver": case.solver,
            "case_id": case.case_id,
            "E_over_N_Td": case.e_over_n_Td,
            "physical_context": asdict(case.physical_context)
            if case.physical_context
            else None,
            "transport_observables": {
                name: asdict(definition)
                for name, definition in (
                    case.transport_observables or case.transport.observable_definitions
                ).items()
            },
        }
        for case in cases
    ]
