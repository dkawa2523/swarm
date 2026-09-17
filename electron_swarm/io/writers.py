"""Canonical product output writers."""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import tempfile

import pandas as pd

from electron_swarm.core.config import ComparisonConfig, OutputConfig
from electron_swarm.core.numerics import eepf_from_eedf, widths_from_centers
from electron_swarm.core.result_metadata import SUMMARY_METADATA_KEYS
from electron_swarm.core.results import (
    SwarmCaseResult,
    SwarmRunResult,
    validate_case_result,
)


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
    "effective_sample_count",
    "relative_standard_error",
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
    "mean_energy_eV_relative_difference",
    "mean_energy_eV_status",
    "drift_velocity_relative_difference",
    "drift_velocity_m_s_status",
    "mobility_relative_difference",
    "mobility_m2_V_s_status",
    "diffusion_L_relative_difference",
    "diffusion_L_m2_s_status",
    "net_ionization_frequency_relative_difference",
    "net_ionization_frequency_s_status",
    "eedf_l1_error",
]


def _summary_columns() -> list[str]:
    return SUMMARY_COLUMNS + [
        f"meta_{key}" for key in SUMMARY_METADATA_COLUMNS
    ]


def _canonical_frame(
    rows: list[dict[str, object]],
    columns: list[str],
) -> pd.DataFrame:
    return pd.DataFrame(rows).reindex(columns=columns)


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
            row[f"meta_{key}"] = case.metadata.get(key, "")
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
                int(counts[index])
                if counts is not None and index < len(counts)
                else None
            )
            effective_count = (
                float(effective_counts[index])
                if effective_counts is not None and index < len(effective_counts)
                else (float(count) if count is not None else None)
            )
            relative_error = (
                float(1.0 / (effective_count ** 0.5))
                if effective_count is not None and effective_count > 0.0
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
                    "effective_sample_count": effective_count,
                    "relative_standard_error": relative_error,
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
    return pd.DataFrame(rows, columns=RATES_COLUMNS)


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
    return pd.DataFrame(rows, columns=ENERGY_ANGLE_COLUMNS)


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
        pd.DataFrame(phase_rows, columns=RF_PHASE_COLUMNS),
        pd.DataFrame(eedf_rows, columns=RF_PHASE_EEDF_COLUMNS),
    )


def _solver_plan_frame(result: SwarmRunResult) -> pd.DataFrame:
    rows = result.metadata.get("solver_plan", [])
    return _canonical_frame(rows, SOLVER_PLAN_COLUMNS)


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


def _write_generation(
    result: SwarmRunResult,
    output: OutputConfig,
    *,
    comparison: ComparisonConfig | None = None,
) -> dict[str, Path]:
    paths: dict[str, Path] = {}

    summary_path = output.directory / f"{output.base_name}_summary.csv"
    _write_csv(_summary_frame(result.cases), summary_path, output.float_format)
    paths["summary_csv"] = summary_path

    eedf_path = output.directory / f"{output.base_name}_eedf.csv"
    _write_csv(_eedf_frame(result.cases), eedf_path, output.float_format)
    paths["eedf_csv"] = eedf_path

    energy_angle = _energy_angle_frame(result.cases)
    if not energy_angle.empty:
        energy_angle_path = (
            output.directory
            / f"{output.base_name}_energy_angle_distribution.csv"
        )
        _write_csv(energy_angle, energy_angle_path, output.float_format)
        paths["energy_angle_distribution_csv"] = energy_angle_path

    rates_path = output.directory / f"{output.base_name}_rates.csv"
    _write_csv(_rates_frame(result.cases), rates_path, output.float_format)
    paths["rates_csv"] = rates_path

    plan_path = output.directory / f"{output.base_name}_solver_plan.csv"
    _write_csv(_solver_plan_frame(result), plan_path, output.float_format)
    paths["solver_plan_csv"] = plan_path

    rf_phase, rf_phase_eedf = _rf_phase_frames(result.cases)
    if not rf_phase.empty:
        phase_path = output.directory / f"{output.base_name}_rf_phase.csv"
        _write_csv(rf_phase, phase_path, output.float_format)
        paths["rf_phase_csv"] = phase_path
        phase_eedf_path = (
            output.directory / f"{output.base_name}_rf_phase_eedf.csv"
        )
        _write_csv(rf_phase_eedf, phase_eedf_path, output.float_format)
        paths["rf_phase_eedf_csv"] = phase_eedf_path

    if comparison is not None:
        _write_comparison_outputs(result, output, comparison, paths)

    return paths


_PRODUCT_OUTPUT_SUFFIXES = (
    "summary.csv",
    "eedf.csv",
    "energy_angle_distribution.csv",
    "rates.csv",
    "solver_plan.csv",
    "rf_phase.csv",
    "rf_phase_eedf.csv",
    "comparison_summary.csv",
)


def write_outputs(
    result: SwarmRunResult,
    output: OutputConfig,
    *,
    comparison: ComparisonConfig | None = None,
) -> dict[str, Path]:
    """Publish one complete product-output generation.

    All frames are built and written in a same-filesystem staging directory.
    Existing files are replaced only after the full generation succeeds, and
    optional files omitted by the new run are removed.
    """

    for case in result.cases:
        validate_case_result(case)

    target_directory = output.directory.resolve()
    target_directory.mkdir(parents=True, exist_ok=True)
    stage_directory = Path(
        tempfile.mkdtemp(
            prefix=f".{output.base_name}.staging-",
            dir=target_directory,
        )
    ).resolve()
    if stage_directory.parent != target_directory:
        raise RuntimeError("product output staging escaped its target directory")

    staged_output = OutputConfig(
        directory=stage_directory,
        base_name=output.base_name,
        float_format=output.float_format,
    )
    try:
        staged_paths = _write_generation(
            result,
            staged_output,
            comparison=comparison,
        )
        final_paths: dict[str, Path] = {}
        published_names: set[str] = set()
        for key, staged_path in staged_paths.items():
            target = target_directory / staged_path.name
            os.replace(staged_path, target)
            final_paths[key] = target
            published_names.add(target.name)

        for suffix in _PRODUCT_OUTPUT_SUFFIXES:
            stale = target_directory / f"{output.base_name}_{suffix}"
            if stale.name not in published_names and stale.exists():
                stale.unlink()
        return final_paths
    finally:
        shutil.rmtree(stage_directory, ignore_errors=True)
