"""Canonical product output writers."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from electron_swarm.core.config import ComparisonConfig, OutputConfig
from electron_swarm.core.results import (
    SUMMARY_METADATA_KEYS,
    SwarmCaseResult,
    SwarmRunResult,
)


SUMMARY_COLUMNS = [
    "solver",
    "solver_method",
    "case_id",
    "E_over_N_Td",
    "mean_energy_eV",
    "drift_velocity_m_s",
    "mobility_m2_V_s",
    "diffusion_L_m2_s",
    "diffusion_T_m2_s",
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
    "eedf",
    "eepf",
]

SOLVER_PLAN_COLUMNS = [
    "solver",
    "requested",
    "runnable",
    "skipped",
    "degraded",
    "skip_reason",
    "warnings",
    "effective_angular_scattering",
    "effective_ionization_source",
    "effective_electron_electron",
    "effective_magnetic_field",
    "effective_tail_refinement",
    "effective_bulk_transport",
    "effective_finite_k",
    "capability_electron_neutral",
    "capability_angular_scattering",
    "capability_ionization_source",
    "capability_electron_electron",
    "capability_magnetic_field",
    "capability_tail_refinement",
    "capability_bulk_transport",
]

COMPARISON_SUMMARY_COLUMNS = [
    "case_id",
    "E_over_N_Td",
    "reference_solver",
    "candidate_solver",
    "angular_model",
    "same_angular_model",
    "angular_model_warning",
    "angular_model_status",
    "reference_angular_model",
    "candidate_angular_model",
    "reference_angular_moment_source",
    "candidate_angular_moment_source",
    "mean_energy_eV_relative_difference",
    "drift_velocity_relative_difference",
    "mobility_relative_difference",
    "diffusion_L_relative_difference",
    "net_ionization_frequency_relative_difference",
    "eedf_l1_error",
]


def _summary_frame(cases: list[SwarmCaseResult]) -> pd.DataFrame:
    rows = []
    for case in cases:
        row = {
            "solver": case.solver,
            "schema_version": case.schema_version,
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
            "net_ionization_frequency_s": case.net_ionization_frequency_s,
            "effective_townsend_m2": case.effective_townsend_m2,
        }
        row["solver_method"] = case.metadata.get("solver_method", "")
        for key in SUMMARY_METADATA_COLUMNS:
            row[f"meta_{key}"] = case.metadata.get(key, "")
        rows.append(row)
    frame = pd.DataFrame(rows)
    if frame.empty:
        return pd.DataFrame(columns=SUMMARY_COLUMNS)
    ordered = [col for col in SUMMARY_COLUMNS if col in frame.columns]
    ordered.extend(
        f"meta_{col}" for col in SUMMARY_METADATA_COLUMNS if f"meta_{col}" in frame.columns
    )
    return frame[ordered]


def _eedf_frame(cases: list[SwarmCaseResult]) -> pd.DataFrame:
    rows = []
    for case in cases:
        for energy, eedf, eepf in zip(
            case.energy_eV, case.eedf, case.eepf, strict=False
        ):
            rows.append(
                {
                    "solver": case.solver,
                    "case_id": case.case_id,
                    "E_over_N_Td": case.e_over_n_Td,
                    "energy_eV": energy,
                    "eedf": eedf,
                    "eepf": eepf,
                }
            )
    return pd.DataFrame(rows, columns=EEDF_COLUMNS)


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
                    "mixture_weighted_rate_m3_s": rate.mixture_weighted_rate_m3_s,
                    "frequency_s_inv": rate.frequency_s_inv,
                    "power_loss_eV_s": rate.power_loss_eV_s,
                    "tail_fraction": rate.tail_fraction,
                }
            )
    return pd.DataFrame(rows, columns=RATES_COLUMNS)


def _solver_plan_frame(result: SwarmRunResult) -> pd.DataFrame:
    rows = result.metadata.get("solver_plan", [])
    frame = pd.DataFrame(rows)
    if frame.empty:
        return pd.DataFrame(columns=SOLVER_PLAN_COLUMNS)
    ordered = [col for col in SOLVER_PLAN_COLUMNS if col in frame.columns]
    return frame[ordered]


def _write_csv(frame: pd.DataFrame, path: Path, float_format: str) -> None:
    frame.to_csv(path, index=False, float_format=float_format)

def _write_comparison_outputs(
    result: SwarmRunResult,
    output: OutputConfig,
    comparison: ComparisonConfig,
    paths: dict[str, Path],
) -> None:
    if not comparison.enabled:
        return
    summary_rows = result.metadata.get("comparison_summary_rows", [])
    summary_path = output.directory / f"{output.base_name}_comparison_summary.csv"
    _write_csv(
        pd.DataFrame(summary_rows, columns=COMPARISON_SUMMARY_COLUMNS),
        summary_path,
        output.float_format,
    )
    paths["comparison_summary_csv"] = summary_path


def write_outputs(
    result: SwarmRunResult,
    output: OutputConfig,
    *,
    comparison: ComparisonConfig | None = None,
) -> dict[str, Path]:
    output.directory.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}

    summary_path = output.directory / f"{output.base_name}_summary.csv"
    _write_csv(_summary_frame(result.cases), summary_path, output.float_format)
    paths["summary_csv"] = summary_path

    eedf_path = output.directory / f"{output.base_name}_eedf.csv"
    _write_csv(_eedf_frame(result.cases), eedf_path, output.float_format)
    paths["eedf_csv"] = eedf_path

    rates_path = output.directory / f"{output.base_name}_rates.csv"
    _write_csv(_rates_frame(result.cases), rates_path, output.float_format)
    paths["rates_csv"] = rates_path

    plan_path = output.directory / f"{output.base_name}_solver_plan.csv"
    _write_csv(_solver_plan_frame(result), plan_path, output.float_format)
    paths["solver_plan_csv"] = plan_path

    if comparison is not None:
        _write_comparison_outputs(result, output, comparison, paths)

    return paths
