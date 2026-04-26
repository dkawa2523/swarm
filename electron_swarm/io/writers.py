"""Unified output writers.

The writer produces one shared schema for Monte Carlo, Boltzmann, and combined
runs. When both solvers are executed, all rows are written to the same unified
files and legacy compatibility tables are emitted per solver so existing tools
can keep reading ``summary.csv`` / ``eedf_table.csv`` / ``energy_table.csv``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from electron_swarm.core.config import OutputConfig
from electron_swarm.core.results import SwarmCaseResult, SwarmRunResult
from electron_swarm.core.solver_registry import solver_legacy_tag
from electron_swarm.io.transport_writer import write_transport_outputs


def _summary_frame(cases: list[SwarmCaseResult]) -> pd.DataFrame:
    return pd.DataFrame([case.summary_dict() for case in cases])


def _eedf_frame(cases: list[SwarmCaseResult]) -> pd.DataFrame:
    rows = []
    for case in cases:
        for e, f, g in zip(case.energy_eV, case.eedf, case.eepf, strict=False):
            rows.append(
                {
                    "solver": case.solver,
                    "case_id": case.case_id,
                    "E_over_N_Td": case.e_over_n_Td,
                    "mean_energy_eV": case.mean_energy_eV,
                    "energy_eV": e,
                    "eedf_eV-1": f,
                    "eepf_eV-3/2": g,
                }
            )
    return pd.DataFrame(rows)


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
                    "type": rate.process_type,
                    "threshold_eV": rate.threshold_eV,
                    "rate_coefficient_m3_s": rate.rate_coefficient_m3_s,
                    "mixture_weighted_rate_m3_s": rate.mixture_weighted_rate_m3_s,
                    "frequency_s_inv": rate.frequency_s_inv,
                    "power_loss_eV_s": rate.power_loss_eV_s,
                }
            )
    return pd.DataFrame(rows)


def _scalar_metadata(case: SwarmCaseResult, key: str, default=np.nan):
    return case.metadata.get(key, default)


def _bulk_metadata(case: SwarmCaseResult, key: str, default=np.nan):
    if case.transport is not None and case.transport.bulk is None:
        return np.nan
    return _scalar_metadata(case, key, default)


def _legacy_convolution_columns(case: SwarmCaseResult) -> dict[str, float]:
    cols: dict[str, float] = {}
    for rate in case.rates:
        base = f"k_conv_{rate.process_type}_{rate.process}"
        key = base
        suffix = 2
        while key in cols:
            key = f"{base}_{suffix}"
            suffix += 1
        cols[key] = float(rate.mixture_weighted_rate_m3_s)
    return cols


def _rate_sum(case: SwarmCaseResult, process_type: str) -> float:
    return float(
        sum(
            float(rate.mixture_weighted_rate_m3_s)
            for rate in case.rates
            if rate.process_type == process_type
        )
    )


def _legacy_summary_row(case: SwarmCaseResult) -> dict[str, object]:
    ion_rate = _rate_sum(case, "ionization")
    attach_rate = _rate_sum(case, "attachment")
    metadata_cols = {
        f"meta_{key}": value
        for key, value in case.metadata.items()
        if isinstance(value, (str, int, float, bool))
    }
    row: dict[str, object] = {
        "solver": case.solver,
        "run_label": _scalar_metadata(case, "run_label", case.case_id),
        "sweep_param": _scalar_metadata(case, "sweep_param", "E_over_N_Td"),
        "sweep_value": _scalar_metadata(case, "sweep_value", case.e_over_n_Td),
        "run_time (s)": _scalar_metadata(case, "run_time_s", np.nan),
        "E/N (Td)": case.e_over_n_Td,
        "mean energy (eV)": case.mean_energy_eV,
        "mean energy error (eV)": _scalar_metadata(
            case, "mean_energy_error_eV", np.nan
        ),
        "bulk drift velocity (m.s-1)": _bulk_metadata(
            case, "bulk_drift_velocity_m_s", case.drift_velocity_m_s
        ),
        "bulk drift velocity error (m.s-1)": _bulk_metadata(
            case, "bulk_drift_velocity_error_m_s", np.nan
        ),
        "bulk L diffusion coeff. * N (m-1.s-1)": _bulk_metadata(
            case,
            "bulk_reduced_diffusion_L_m-1_s-1",
            case.reduced_diffusion_L_m2_s_m3,
        ),
        "bulk L diffusion coeff. error * N (m-1.s-1)": _bulk_metadata(
            case, "bulk_reduced_diffusion_L_error_m-1_s-1", np.nan
        ),
        "bulk T diffusion coeff. * N (m-1.s-1)": _bulk_metadata(
            case,
            "bulk_reduced_diffusion_T_m-1_s-1",
            case.reduced_diffusion_T_m2_s_m3,
        ),
        "bulk T diffusion coeff. error * N (m-1.s-1)": _bulk_metadata(
            case, "bulk_reduced_diffusion_T_error_m-1_s-1", np.nan
        ),
        "flux drift velocity (m.s-1)": case.drift_velocity_m_s,
        "flux drift velocity error (m.s-1)": _scalar_metadata(
            case, "flux_drift_velocity_error_m_s", np.nan
        ),
        "flux L diffusion coeff. * N (m-1.s-1)": case.reduced_diffusion_L_m2_s_m3,
        "flux L diffusion coeff. error * N (m-1.s-1)": _scalar_metadata(
            case, "flux_reduced_diffusion_L_error_m-1_s-1", np.nan
        ),
        "flux T diffusion coeff. * N (m-1.s-1)": case.reduced_diffusion_T_m2_s_m3,
        "flux T diffusion coeff. error * N (m-1.s-1)": _scalar_metadata(
            case, "flux_reduced_diffusion_T_error_m-1_s-1", np.nan
        ),
        "effective ionization rate coeff. (counted) (m3.s-1)": _scalar_metadata(
            case, "counted_effective_rate_coefficient_m3_s", np.nan
        ),
        "effective ionization rate coeff. (counted) error (m3.s-1)": _scalar_metadata(
            case, "counted_effective_rate_error_m3_s", np.nan
        ),
        "ionization rate coeff. (counted) (m3.s-1)": _scalar_metadata(
            case, "counted_ionization_rate_coefficient_m3_s", np.nan
        ),
        "ionization rate coeff. (counted) error (m3.s-1)": _scalar_metadata(
            case, "counted_ionization_rate_error_m3_s", np.nan
        ),
        "attachment rate coeff. (counted) (m3.s-1)": _scalar_metadata(
            case, "counted_attachment_rate_coefficient_m3_s", np.nan
        ),
        "attachment rate coeff. (counted) error (m3.s-1)": _scalar_metadata(
            case, "counted_attachment_rate_error_m3_s", np.nan
        ),
        "effective ionization rate coeff. (convolution) (m3.s-1)": _scalar_metadata(
            case, "convolution_effective_rate_coefficient_m3_s", ion_rate - attach_rate
        ),
        "ionization rate coeff. (convolution) (m3.s-1)": _scalar_metadata(
            case, "convolution_ionization_rate_coefficient_m3_s", ion_rate
        ),
        "attachment rate coeff. (convolution) (m3.s-1)": _scalar_metadata(
            case, "convolution_attachment_rate_coefficient_m3_s", attach_rate
        ),
    }
    row.update(_legacy_convolution_columns(case))
    row.update(metadata_cols)
    return row


def _legacy_summary_frame(cases: list[SwarmCaseResult]) -> pd.DataFrame:
    rows = [_legacy_summary_row(case) for case in cases]
    frame = pd.DataFrame(rows)
    if "E/N (Td)" in frame.columns:
        frame = frame.sort_values(["E/N (Td)", "run_label"]).reset_index(drop=True)
    return frame


def _legacy_distribution_rows(
    case: SwarmCaseResult, *, include_sweep_metadata: bool
) -> list[dict[str, object]]:
    run_label = _scalar_metadata(case, "run_label", case.case_id)
    base_row = {
        "solver": case.solver,
        "run_label": run_label,
        "E/N (Td)": case.e_over_n_Td,
        "mean energy (eV)": case.mean_energy_eV,
    }
    if include_sweep_metadata:
        base_row["sweep_param"] = _scalar_metadata(
            case, "sweep_param", "E_over_N_Td"
        )
        base_row["sweep_value"] = _scalar_metadata(
            case, "sweep_value", case.e_over_n_Td
        )

    rows = []
    for energy_eV, eedf in zip(case.energy_eV, case.eedf, strict=False):
        row = dict(base_row)
        row["energy (eV)"] = energy_eV
        row["eedf (eV-1)"] = eedf
        rows.append(row)
    return rows


def _legacy_energy_frame(cases: list[SwarmCaseResult]) -> pd.DataFrame:
    rows = []
    for case in cases:
        rows.extend(_legacy_distribution_rows(case, include_sweep_metadata=True))
    frame = pd.DataFrame(rows)
    if "mean energy (eV)" in frame.columns:
        frame = frame.sort_values(
            ["mean energy (eV)", "energy (eV)"]
        ).reset_index(drop=True)
    return frame


def _legacy_eedf_table_frame(cases: list[SwarmCaseResult]) -> pd.DataFrame:
    rows = []
    for case in cases:
        rows.extend(_legacy_distribution_rows(case, include_sweep_metadata=False))
    frame = pd.DataFrame(rows)
    if "mean energy (eV)" in frame.columns:
        frame = frame.sort_values(
            ["mean energy (eV)", "energy (eV)"]
        ).reset_index(drop=True)
    return frame


def _solver_groups(
    cases: list[SwarmCaseResult],
) -> dict[str, list[SwarmCaseResult]]:
    groups: dict[str, list[SwarmCaseResult]] = {}
    for case in cases:
        groups.setdefault(case.solver, []).append(case)
    return groups


def _choose_primary_solver(
    output: OutputConfig, groups: dict[str, list[SwarmCaseResult]]
) -> str | None:
    preferred = output.compatibility.primary_solver
    if preferred in groups:
        return preferred
    if groups:
        return next(iter(groups))
    return None


def _write_csv(frame: pd.DataFrame, path: Path, float_format: str) -> None:
    frame.to_csv(path, index=False, float_format=float_format)


def write_outputs(result: SwarmRunResult, output: OutputConfig) -> dict[str, Path]:
    output.directory.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}

    summary = _summary_frame(result.cases)
    summary_path = output.directory / f"{output.base_name}_summary.csv"
    _write_csv(summary, summary_path, output.float_format)
    paths["summary_csv"] = summary_path

    if output.write_eedf:
        eedf = _eedf_frame(result.cases)
        eedf_path = output.directory / f"{output.base_name}_eedf.csv"
        _write_csv(eedf, eedf_path, output.float_format)
        paths["eedf_csv"] = eedf_path

    if output.write_rates:
        rates = _rates_frame(result.cases)
        rates_path = output.directory / f"{output.base_name}_rates.csv"
        _write_csv(rates, rates_path, output.float_format)
        paths["rates_csv"] = rates_path

    transport_path = write_transport_outputs(result, output)
    if transport_path is not None:
        paths["transport_csv"] = transport_path

    if not output.compatibility.write_legacy_tables:
        return paths

    groups = _solver_groups(result.cases)
    primary_solver = _choose_primary_solver(output, groups)
    primary_summary: pd.DataFrame | None = None
    primary_energy: pd.DataFrame | None = None
    primary_eedf: pd.DataFrame | None = None

    for solver, cases in groups.items():
        tag = solver_legacy_tag(solver)
        solver_summary = _legacy_summary_frame(cases)
        solver_energy = _legacy_energy_frame(cases)
        solver_eedf = _legacy_eedf_table_frame(cases)

        summary_path = output.directory / f"summary_{tag}.csv"
        energy_path = output.directory / f"energy_table_{tag}.csv"
        eedf_path = output.directory / f"eedf_table_{tag}.csv"
        _write_csv(solver_summary, summary_path, output.float_format)
        _write_csv(solver_energy, energy_path, output.float_format)
        _write_csv(solver_eedf, eedf_path, output.float_format)

        paths[f"legacy_summary_{tag}_csv"] = summary_path
        paths[f"legacy_energy_table_{tag}_csv"] = energy_path
        paths[f"legacy_eedf_table_{tag}_csv"] = eedf_path

        if solver == primary_solver:
            primary_summary = solver_summary
            primary_energy = solver_energy
            primary_eedf = solver_eedf

    if primary_summary is not None:
        summary_alias = output.directory / "summary.csv"
        energy_alias = output.directory / "energy_table.csv"
        eedf_alias = output.directory / "eedf_table.csv"
        _write_csv(primary_summary, summary_alias, output.float_format)
        _write_csv(primary_energy, energy_alias, output.float_format)
        _write_csv(primary_eedf, eedf_alias, output.float_format)
        paths["legacy_summary_csv"] = summary_alias
        paths["legacy_energy_table_csv"] = energy_alias
        paths["legacy_eedf_table_csv"] = eedf_alias

    return paths
