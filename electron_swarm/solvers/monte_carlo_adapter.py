"""Adapter for the existing particle Monte Carlo solver.

The current repository already contains the MC implementation. This adapter is
an integration seam: it can call the existing code through either a command line
or a Python function and then normalize the resulting files/objects into the
shared ``SwarmCaseResult`` schema.
"""

from __future__ import annotations

import importlib
import os
import shlex
import subprocess
from typing import Any

import numpy as np
import pandas as pd

from electron_swarm.core.constants import BOLTZMANN_J_K
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.physics.angular_scattering import (
    ANGULAR_METADATA_KEYS,
    expected_angular_metadata,
)
from .base import SwarmSolver
from .internal_monte_carlo import run_internal_monte_carlo

MAGNETIC_METADATA_KEYS = (
    "magnetic_field_treatment",
    "magnetic_field_B_T",
    "magnetic_field_angle_EB_deg",
    "field_integrator",
)


def _first_present(row: pd.Series, *names: str) -> float:
    for name in names:
        if name in row and pd.notna(row[name]):
            return float(row[name])
    return float("nan")


def _first_metadata_value(row: pd.Series, key: str) -> Any:
    for name in (key, f"meta_{key}"):
        if name in row and pd.notna(row[name]) and row[name] != "":
            return row[name]
    return None


def _coerce_metadata_value(key: str, value: Any) -> object:
    if key in {"exact_dcs_based", "ordinary_integral_xs_closure"}:
        if isinstance(value, bool):
            return value
        if isinstance(value, (int, float)) and np.isfinite(float(value)):
            return bool(value)
        normalized = str(value).strip().lower()
        if normalized in {"true", "1", "yes"}:
            return True
        if normalized in {"false", "0", "no"}:
            return False
    return str(value)


def _metadata_float(metadata: dict[str, Any], key: str) -> float:
    value = metadata[key]
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"monte_carlo magnetic_field metadata {key} must be numeric") from exc
    if not np.isfinite(out):
        raise ValueError(f"monte_carlo magnetic_field metadata {key} must be finite")
    return out


def _infer_gas_number_density(config, row: pd.Series) -> float:
    for name in (
        "meta_gas_number_density_m-3",
        "gas_number_density_m-3",
        "gas_number_density_m3",
    ):
        if name in row and pd.notna(row[name]):
            return float(row[name])
    cond = config.conditions
    if cond.gas_number_density_m3 is not None:
        return float(cond.gas_number_density_m3)
    if cond.pressure_Pa is not None:
        return float(cond.pressure_Pa / (BOLTZMANN_J_K * cond.gas_temperature_K))
    return float("nan")


def _infer_effective_frequency(row: pd.Series, gas_number_density: float) -> float:
    direct = _first_present(row, "net_ionization_frequency_s")
    if np.isfinite(direct):
        return direct
    for rate_name in (
        "meta_convolution_effective_rate_coefficient_m3_s",
        "effective ionization rate coeff. (convolution) (m3.s-1)",
        "meta_counted_effective_rate_coefficient_m3_s",
        "effective ionization rate coeff. (counted) (m3.s-1)",
    ):
        rate = _first_present(row, rate_name)
        if np.isfinite(rate) and np.isfinite(gas_number_density):
            return float(rate * gas_number_density)
    return float("nan")


def _infer_effective_townsend(
    row: pd.Series, frequency_s: float, drift_velocity_m_s: float, gas_number_density: float
) -> float:
    direct = _first_present(row, "effective_townsend_m2")
    if np.isfinite(direct):
        return direct
    townsend_1_m = _first_present(row, "meta_effective_townsend_1_m")
    if np.isfinite(townsend_1_m) and np.isfinite(gas_number_density):
        return float(townsend_1_m / gas_number_density)
    if (
        np.isfinite(frequency_s)
        and np.isfinite(drift_velocity_m_s)
        and np.isfinite(gas_number_density)
        and drift_velocity_m_s != 0.0
    ):
        return float(frequency_s / (abs(drift_velocity_m_s) * gas_number_density))
    return float("nan")


class MonteCarloAdapter(SwarmSolver):
    name = "monte_carlo"

    def solve_all(self) -> list[SwarmCaseResult]:
        cfg = self.config.internal.monte_carlo
        if cfg.backend == "internal":
            cases = run_internal_monte_carlo(self.config, self.cross_sections)
            return self._apply_magnetic_contract(self._apply_angular_contract(cases))
        if not cfg.command and not cfg.python_api:
            raise RuntimeError(
                "run.solvers includes monte_carlo but solvers.monte_carlo.command or "
                "solvers.monte_carlo.python_api is not configured. Set one of these to "
                "delegate to the repository's existing MC solver."
            )
        if cfg.python_api:
            cases = self._run_python_api()
        else:
            cases = self._run_command()
        return self._apply_magnetic_contract(self._apply_angular_contract(cases))

    def solve_case(self, e_over_n_Td: float, case_id: str) -> SwarmCaseResult:
        raise NotImplementedError("MonteCarloAdapter executes through solve_all()")

    def _run_python_api(self) -> list[SwarmCaseResult]:
        cfg = self.config.internal.monte_carlo
        assert cfg.python_api is not None
        module_name, func_name = cfg.python_api.split(":", 1)
        func = getattr(importlib.import_module(module_name), func_name)
        obj = func(self.config, self.cross_sections, **cfg.passthrough)
        if isinstance(obj, list) and all(isinstance(x, SwarmCaseResult) for x in obj):
            return obj
        if isinstance(obj, dict):
            return self._results_from_mapping(obj)
        if cfg.output_summary_csv:
            return self._read_existing_outputs()
        raise TypeError(
            "monte_carlo.python_api must return list[SwarmCaseResult], a supported mapping, "
            "or write output_summary_csv/output_eedf_csv configured in YAML."
        )

    def _run_command(self) -> list[SwarmCaseResult]:
        cfg = self.config.internal.monte_carlo
        assert cfg.command is not None
        env = os.environ.copy()
        env.update(cfg.environment)
        command = cfg.command.format(
            config=self.config.source_path or "",
            output_dir=self.config.output.directory,
        )
        subprocess.run(
            shlex.split(command),
            cwd=str(cfg.working_directory) if cfg.working_directory else None,
            env=env,
            check=True,
            timeout=cfg.timeout_s,
        )
        return self._read_existing_outputs()

    def _read_existing_outputs(self) -> list[SwarmCaseResult]:
        cfg = self.config.internal.monte_carlo
        if cfg.output_summary_csv is None:
            raise RuntimeError(
                "monte_carlo.output_summary_csv is required to parse command output"
            )
        summary = pd.read_csv(cfg.output_summary_csv)
        eedf_df = pd.read_csv(cfg.output_eedf_csv) if cfg.output_eedf_csv else None
        return self._from_dataframes(summary, eedf_df)

    def _results_from_mapping(self, obj: dict[str, Any]) -> list[SwarmCaseResult]:
        summary = pd.DataFrame(obj.get("summary", obj.get("cases", [])))
        eedf = pd.DataFrame(obj.get("eedf", [])) if obj.get("eedf") is not None else None
        return self._from_dataframes(summary, eedf)

    def _find_en_column(self, summary: pd.DataFrame) -> str | None:
        cols = {c.lower(): c for c in summary.columns}
        for candidate in ("e_over_n_td", "e/n (td)", "e/n_td", "e_n_td"):
            if candidate in cols:
                return cols[candidate]
        for column in summary.columns:
            normalized = (
                str(column)
                .lower()
                .replace("/", "_")
                .replace("(", "")
                .replace(")", "")
                .replace(" ", "_")
            )
            if normalized in {"e_over_n_td", "e_n_td"}:
                return column
        return None

    def _from_dataframes(
        self, summary: pd.DataFrame, eedf: pd.DataFrame | None
    ) -> list[SwarmCaseResult]:
        if summary.empty:
            return []
        e_col = self._find_en_column(summary)
        if e_col is None:
            raise ValueError("MC summary CSV must include E_over_N_Td or E/N (Td)")

        out: list[SwarmCaseResult] = []
        for i, row in summary.iterrows():
            e_over_n = float(row[e_col])
            case_id = str(
                row.get(
                    "case_id",
                    row.get("run_label", f"{self.config.run.case_prefix}_{i:04d}"),
                )
            )
            energy = np.array([0.0, 1.0], dtype=float)
            dist = np.array([1.0, 0.0], dtype=float)
            if eedf is not None and not eedf.empty:
                rows = eedf
                if "case_id" in rows.columns:
                    rows = rows[rows["case_id"].astype(str) == case_id]
                elif "run_label" in rows.columns:
                    rows = rows[rows["run_label"].astype(str) == case_id]
                elif "E_over_N_Td" in rows.columns:
                    rows = rows[np.isclose(rows["E_over_N_Td"].astype(float), e_over_n)]
                elif "E/N (Td)" in rows.columns:
                    rows = rows[np.isclose(rows["E/N (Td)"].astype(float), e_over_n)]
                if not rows.empty:
                    energy_col = (
                        "energy_eV"
                        if "energy_eV" in rows.columns
                        else "energy (eV)"
                        if "energy (eV)" in rows.columns
                        else rows.columns[0]
                    )
                    eedf_col = (
                        "eedf"
                        if "eedf" in rows.columns
                        else "eedf (eV-1)"
                        if "eedf (eV-1)" in rows.columns
                        else rows.columns[-1]
                    )
                    energy = rows[energy_col].to_numpy(dtype=float)
                    dist = rows[eedf_col].to_numpy(dtype=float)
                    integral = np.trapezoid(dist, energy) if len(energy) > 1 else 1.0
                    if integral > 0:
                        dist = dist / integral

            drift = _first_present(
                row, "drift_velocity_m_s", "flux drift velocity (m.s-1)"
            )
            gas_number_density = _infer_gas_number_density(self.config, row)
            mobility = _first_present(row, "mobility_m2_V_s")
            reduced_mobility = _first_present(
                row,
                "reduced_mobility_m2_V_s_m3",
            )
            if np.isnan(reduced_mobility) and np.isfinite(mobility) and np.isfinite(
                gas_number_density
            ):
                reduced_mobility = mobility * gas_number_density
            diffusion_l = _first_present(
                row, "diffusion_L_m2_s", "diffusion_m2_s"
            )
            diffusion_t = _first_present(
                row, "diffusion_T_m2_s", "diffusion_m2_s"
            )
            reduced_diffusion_l = _first_present(
                row,
                "reduced_diffusion_L_m2_s_m3",
                "flux L diffusion coeff. * N (m-1.s-1)",
            )
            reduced_diffusion_t = _first_present(
                row,
                "reduced_diffusion_T_m2_s_m3",
                "flux T diffusion coeff. * N (m-1.s-1)",
            )
            if np.isnan(reduced_diffusion_l) and np.isfinite(diffusion_l) and np.isfinite(
                gas_number_density
            ):
                reduced_diffusion_l = diffusion_l * gas_number_density
            if np.isnan(reduced_diffusion_t) and np.isfinite(diffusion_t) and np.isfinite(
                gas_number_density
            ):
                reduced_diffusion_t = diffusion_t * gas_number_density
            net_ionization_frequency = _infer_effective_frequency(
                row, gas_number_density
            )
            effective_townsend = _infer_effective_townsend(
                row, net_ionization_frequency, drift, gas_number_density
            )
            metadata: dict[str, Any] = {"adapter": "command_or_python_api"}
            for key in ANGULAR_METADATA_KEYS:
                value = _first_metadata_value(row, key)
                if value is not None:
                    metadata[key] = _coerce_metadata_value(key, value)
            for key in MAGNETIC_METADATA_KEYS:
                value = _first_metadata_value(row, key)
                if value is not None:
                    metadata[key] = (
                        float(value)
                        if key in {"magnetic_field_B_T", "magnetic_field_angle_EB_deg"}
                        else str(value)
                    )
            out.append(
                SwarmCaseResult(
                    solver=self.name,
                    case_id=case_id,
                    e_over_n_Td=e_over_n,
                    mean_energy_eV=_first_present(
                        row, "mean_energy_eV", "mean energy (eV)"
                    ),
                    drift_velocity_m_s=drift,
                    mobility_m2_V_s=mobility,
                    reduced_mobility_m2_V_s_m3=reduced_mobility,
                    diffusion_L_m2_s=diffusion_l,
                    diffusion_T_m2_s=diffusion_t,
                    reduced_diffusion_L_m2_s_m3=reduced_diffusion_l,
                    reduced_diffusion_T_m2_s_m3=reduced_diffusion_t,
                    net_ionization_frequency_s=net_ionization_frequency,
                    effective_townsend_m2=effective_townsend,
                    energy_eV=energy,
                    eedf=dist,
                    eepf=dist / np.sqrt(np.maximum(energy, 1.0e-30)),
                    rates=[],
                    metadata=metadata,
                )
            )
        return out

    def _apply_angular_contract(
        self, cases: list[SwarmCaseResult]
    ) -> list[SwarmCaseResult]:
        mode = self.config.solvers.monte_carlo.angular_scattering
        expected = expected_angular_metadata(self.config)
        for case in cases:
            reported = {
                key: _coerce_metadata_value(key, case.metadata[key])
                for key in ANGULAR_METADATA_KEYS
                if key in case.metadata and case.metadata[key] not in {"", None}
            }
            if mode == "same_as_physics":
                missing = [key for key in ANGULAR_METADATA_KEYS if key not in reported]
                if missing:
                    raise ValueError(
                        "monte_carlo same_as_physics requires angular metadata "
                        f"{missing} in MC output"
                    )
                mismatched = [
                    key
                    for key in ANGULAR_METADATA_KEYS
                    if reported[key] != expected[key]
                ]
                if mismatched:
                    raise ValueError(
                        "monte_carlo angular metadata mismatch for "
                        f"{mismatched}; same_as_physics requires MC output to "
                        "match physics.angular_scattering"
                    )
                case.metadata.update(expected)
                case.metadata["monte_carlo_angular_scattering"] = "same_as_physics"
            else:
                if not reported:
                    case.metadata.update(
                        {
                            "angular_model": "unknown",
                            "angular_moment_source": "external_adapter",
                            "exact_dcs_based": False,
                            "ordinary_integral_xs_closure": False,
                        }
                    )
                else:
                    case.metadata.update(reported)
                case.metadata["monte_carlo_angular_scattering"] = "external"
        return cases

    def _apply_magnetic_contract(
        self, cases: list[SwarmCaseResult]
    ) -> list[SwarmCaseResult]:
        magnetic = self.config.physics.field.magnetic_field
        if not magnetic.enabled:
            return cases
        required = {
            "magnetic_field_treatment",
            "magnetic_field_B_T",
            "magnetic_field_angle_EB_deg",
        }
        for case in cases:
            missing = [
                key
                for key in required
                if key not in case.metadata or case.metadata[key] in {"", None}
            ]
            if missing:
                raise ValueError(
                    "monte_carlo magnetic_field requires MC output metadata "
                    f"{missing}"
                )
            treatment = str(case.metadata["magnetic_field_treatment"])
            if treatment in {"none", "unsupported"}:
                raise ValueError(
                    "monte_carlo magnetic_field output must report an active "
                    "magnetic_field_treatment"
                )
            reported_B_T = _metadata_float(case.metadata, "magnetic_field_B_T")
            reported_angle = _metadata_float(
                case.metadata, "magnetic_field_angle_EB_deg"
            )
            if not np.isclose(reported_B_T, float(magnetic.B_T), rtol=1.0e-9, atol=1.0e-12):
                raise ValueError(
                    "monte_carlo magnetic_field metadata mismatch for magnetic_field_B_T"
                )
            if not np.isclose(
                reported_angle,
                float(magnetic.angle_EB_deg),
                rtol=1.0e-9,
                atol=1.0e-9,
            ):
                raise ValueError(
                    "monte_carlo magnetic_field metadata mismatch for "
                    "magnetic_field_angle_EB_deg"
                )
        return cases
