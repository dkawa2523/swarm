"""Audit GEC-CCP conservation exports and domain phase state."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any

import numpy as np

from swarm_workflow._io import write_json
from swarm_workflow.comsol.models.gec_ccp.contracts import GecCcpPlan, GecCcpWorkflowError


def audit_gec_ccp_conservation_run(plan: GecCcpPlan) -> Path:
    """Evaluate the fixed GEC weak-balance export contract.

    A failed balance is a quality result, not a false solver failure.  The
    caller records it separately from COMSOL execution completion.
    """

    balance_limit = 0.01
    power_limit = 1.0e-3
    cases: dict[str, Any] = {}
    result_directories = {
        "external": plan.output_directory / "swarm_tables"
    }
    if plan.mapping.run.include_builtin_reference:
        result_directories["baseline"] = (
            plan.output_directory / "builtin_druyvesteyn"
        )
    for name, directory in result_directories.items():
        volume = _read_fixed_comsol_table(
            directory / "conservation_volume.csv", 6
        )
        wall = _read_fixed_comsol_table(
            directory / "conservation_wall.csv", 8
        )
        terminal = _read_fixed_comsol_table(
            directory / "conservation_terminal_power.csv", 4
        )
        electron = _conservation_balance(
            volume[0], wall[0], sign=1.0, limit=balance_limit,
            unit="1/s",
        )
        argon_ion = _conservation_balance(
            volume[1], wall[1], sign=1.0, limit=balance_limit,
            unit="1/s",
        )
        electron_energy = _conservation_balance(
            volume[2] + volume[3] + volume[4], wall[2], sign=1.0,
            limit=balance_limit,
            unit="W",
        )
        weak_energy_identity_residual = wall[2] - (
            -wall[3] + wall[4] + wall[5]
        )
        gross_energy = {
            "absorbed_power_W": volume[2],
            "collisional_signed_power_W": volume[3],
            "external_elastic_signed_power_W": volume[4],
            "random_motion_wall_loss_W": wall[3],
            "secondary_wall_input_W": wall[4],
            "thermionic_wall_input_W": wall[5],
            "gross_input_W": volume[2] + wall[4] + wall[5],
            "gross_output_W": -volume[3] - volume[4] + wall[3],
            "gross_residual_W": (
                volume[2] + wall[4] + wall[5]
                + volume[3] + volume[4] - wall[3]
            ),
            "weak_term_identity_residual_W": weak_energy_identity_residual,
            "weak_term_identity_passed": abs(weak_energy_identity_residual)
            <= 1.0e-10 * max(abs(wall[2]), abs(wall[3]), 1.0),
        }
        gross_energy["gross_relative_residual"] = abs(
            gross_energy["gross_residual_W"]
        ) / max(
            abs(gross_energy["gross_input_W"]),
            abs(gross_energy["gross_output_W"]),
            1.0e-300,
        )
        gross_energy["gross_relative_residual_limit"] = balance_limit
        gross_energy["gross_balance_passed"] = (
            gross_energy["gross_relative_residual"] <= balance_limit
        )
        terminal_power = _conservation_balance(
            terminal[0], terminal[1], sign=-1.0, limit=power_limit,
            unit="W",
        )
        terminal_power["exported_residual"] = terminal[2]
        terminal_power["exported_residual_consistent"] = bool(
            math.isclose(
                terminal_power["residual"],
                terminal[2],
                rel_tol=1.0e-8,
                abs_tol=1.0e-12,
            )
        )
        terminal_power["period_averaged_current_A"] = terminal[3]
        geometry_error = abs(
            volume[5] - 0.00260191657233603
        ) / 0.00260191657233603
        geometry = {
            "axisymmetric_plasma_volume_m3": volume[5],
            "expected_m3": 0.00260191657233603,
            "relative_error": geometry_error,
            "limit": 1.0e-9,
            "passed": geometry_error <= 1.0e-9,
        }
        state = _audit_domain_phase_state(
            directory / "domain_phase_closure.csv"
        )
        balances = {
            "electron_particles": electron,
            "argon_ions": argon_ion,
            "electron_energy": electron_energy,
            "terminal_power": terminal_power,
        }
        strong_trace_diagnostics = {
            "electron_particles": _conservation_balance(
                volume[0], wall[6], sign=-1.0, limit=balance_limit,
                unit="1/s",
            ),
            "electron_energy": _conservation_balance(
                volume[2] + volume[3] + volume[4], wall[7], sign=-1.0,
                limit=balance_limit,
                unit="W",
            ),
            "quality_gate": False,
            "interpretation": (
                "down(gflux) is an element strong trace, not the P1 log-FEM "
                "natural weak boundary flux"
            ),
        }
        passed = all(item["passed"] for item in balances.values()) and (
            geometry["passed"]
        ) and state["passed"] and gross_energy["weak_term_identity_passed"] and gross_energy[
            "gross_balance_passed"
        ]
        cases[name] = {
            "status": "passed" if passed else "failed",
            "balances": balances,
            "geometry": geometry,
            "domain_phase_state": state,
            "strong_trace_mesh_diagnostics": strong_trace_diagnostics,
            "electron_energy_gross_channels": gross_energy,
        }
    all_cases_status = (
        "passed"
        if all(case["status"] == "passed" for case in cases.values())
        else "failed"
    )
    status = cases["external"]["status"]
    output = plan.output_directory / "conservation_audit.json"
    write_json(
        output,
        {
            "stage": "audit-gec-ccp-conservation",
            "status": status,
            "all_cases_status": all_cases_status,
            "quality_gate_case": "external",
            "quality_gate_only": True,
            "limits": {
                "particle_and_electron_energy_relative_residual": (
                    balance_limit
                ),
                "terminal_power_relative_residual": power_limit,
            },
            "axisymmetric_measure": {
                "volume": "IntSurface intvolume=true",
                "wall": "IntLine intsurface=true",
                "manual_2_pi_r": False,
            },
            "cases": cases,
            "interpretation": (
                "Primary particle and electron-energy balances use the exact "
                "weak boundary terms. down(gflux) is retained only as a "
                "non-gating mesh diagnostic."
            ),
        },
    )
    return output


def _audit_domain_phase_state(path: Path) -> dict[str, Any]:
    headers, values = _read_comsol_numeric_csv(path)
    specifications = {
        "electron_density_m3": ("ptp.ne ", 0.0, None),
        "mean_electron_energy_eV": ("ptp.ebar ", 0.0, None),
        "argon_ion_mass_fraction": ("ptp.wAr_1p ", 0.0, 1.0),
    }
    quantities: dict[str, Any] = {}
    passed = True
    for name, (prefix, lower, upper) in specifications.items():
        indices = [
            index for index, header in enumerate(headers)
            if header.startswith(prefix)
        ]
        if not indices:
            quantities[name] = {
                "passed": False,
                "reason": "phase-resolved columns are absent",
            }
            passed = False
            continue
        selected = values[:, indices]
        finite = bool(np.all(np.isfinite(selected)))
        minimum = float(np.min(selected)) if finite else math.nan
        maximum = float(np.max(selected)) if finite else math.nan
        range_valid = finite and minimum > lower and (
            upper is None or maximum <= upper * (1.0 + 1.0e-9)
        )
        passed = passed and range_valid
        quantities[name] = {
            "passed": range_valid,
            "finite": finite,
            "strictly_positive": finite and minimum > lower,
            "minimum": minimum,
            "maximum": maximum,
            "upper_limit": upper,
            "phase_columns": len(indices),
            "spatial_rows": int(values.shape[0]),
        }
    return {
        "passed": passed,
        "coverage": "all exported plasma-domain nodes and all RF phases",
        "quantities": quantities,
    }


def _conservation_balance(
    volume: float,
    wall: float,
    *,
    sign: float,
    limit: float,
    unit: str,
) -> dict[str, Any]:
    residual = volume + sign * wall
    scale = max(abs(volume), abs(wall), 1.0e-300)
    relative = abs(residual) / scale
    return {
        "volume_component": volume,
        "wall_component": wall,
        "wall_sign_in_residual": sign,
        "residual": residual,
        "relative_residual": relative,
        "limit": limit,
        "unit": unit,
        "passed": relative <= limit,
    }


def _read_comsol_numeric_csv(path: Path) -> tuple[list[str], np.ndarray]:
    if not path.exists():
        raise GecCcpWorkflowError(f"missing COMSOL numeric export: {path}")
    lines = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    header_line: str | None = None
    data_lines: list[str] = []
    for line in lines:
        if line.startswith("%"):
            candidate = line[1:].lstrip()
            if candidate.startswith("R,") or candidate.startswith("Z,"):
                header_line = candidate
            continue
        if line.strip():
            data_lines.append(line)
    if header_line is None or not data_lines:
        raise GecCcpWorkflowError(
            f"COMSOL numeric export has no readable header/data: {path}"
        )
    headers = next(csv.reader([header_line]))
    rows = list(csv.reader(data_lines))
    if any(len(row) != len(headers) for row in rows):
        raise GecCcpWorkflowError(
            f"COMSOL numeric export is not rectangular: {path}"
        )
    try:
        values = np.asarray(rows, dtype=float)
    except ValueError as exc:
        raise GecCcpWorkflowError(
            f"COMSOL numeric export contains a nonnumeric value: {path}"
        ) from exc
    return headers, values


def _read_fixed_comsol_table(path: Path, value_count: int) -> list[float]:
    """Read the final row of this runner's fixed Derived Values table."""

    if not path.exists():
        raise GecCcpWorkflowError(f"missing COMSOL table export: {path}")
    rows: list[list[float]] = []
    for line in path.read_text(
        encoding="utf-8-sig", errors="replace"
    ).splitlines():
        if not line.strip() or line.lstrip().startswith("%"):
            continue
        try:
            numeric = [float(value) for value in next(csv.reader([line]))]
        except ValueError:
            continue
        if len(numeric) >= value_count:
            rows.append(numeric[-value_count:])
    if not rows or any(not math.isfinite(value) for value in rows[-1]):
        raise GecCcpWorkflowError(
            f"COMSOL table export lacks {value_count} finite values: {path}"
        )
    return rows[-1]
