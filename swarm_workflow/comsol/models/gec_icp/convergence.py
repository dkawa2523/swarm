"""Physical-time convergence assessment for GEC-ICP runs."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any

from swarm_workflow._io import write_json

from .java import TRANSPORT_TABLE, output_times, read_mobility_table
from .run_contracts import GecIcpPlan, GecIcpWorkflowError


def assess_gec_icp_convergence(plan: GecIcpPlan) -> dict[str, Any]:
    """Assess a fixed late-time window and the lookup-table support."""

    output = plan.output_directory
    volume = _numeric_rows(output / "convergence_volume.csv", columns=7)
    minimum = _numeric_rows(output / "convergence_mean_energy_min.csv", columns=2)
    maximum = _numeric_rows(output / "convergence_mean_energy_max.csv", columns=2)
    coil = _numeric_rows(output / "convergence_coil_power.csv", columns=2)
    _require_same_times(volume, minimum, maximum, coil)
    expected_times = output_times(plan.mapping)
    _require_expected_times(volume, expected_times)
    if len(volume) < 4:
        raise GecIcpWorkflowError(
            "GEC-ICP convergence requires at least four planned time rows"
        )
    late_volume = volume[-4:]
    late_means = [
        _ratio(row[3], row[2], "electron-weighted mean energy") for row in late_volume
    ]
    final = volume[-1]
    final_time = final[0]
    expected_time = plan.mapping.run.final_time_s
    time_tolerance = max(1.0e-12, 1.0e-9 * expected_time)

    final_mean = late_means[-1]
    observed_minimum = min(row[1] for row in minimum)
    observed_maximum = max(row[1] for row in maximum)
    support_axis, _ = read_mobility_table(plan.mapping.bundle.path / TRANSPORT_TABLE)
    support_minimum, support_maximum = support_axis[0], support_axis[-1]
    tolerance = plan.mapping.run.convergence
    support_guard = tolerance.mean_energy_support_relative_guard
    required_support_maximum = observed_maximum * (1.0 + support_guard)

    metrics = {
        "electron_inventory_relative_change": _maximum_adjacent_relative_change(
            [row[2] for row in late_volume]
        ),
        "ion_inventory_relative_change": _maximum_adjacent_relative_change(
            [row[5] for row in late_volume]
        ),
        "metastable_inventory_relative_change": (
            _maximum_adjacent_relative_change([row[4] for row in late_volume])
        ),
        "mean_energy_relative_change": _maximum_adjacent_relative_change(late_means),
        "absorbed_power_relative_change": _maximum_adjacent_relative_change(
            [row[6] for row in late_volume]
        ),
        "coil_power_relative_error": _relative_change(
            plan.mapping.run.power_W, coil[-1][1]
        ),
    }
    thresholds = {name: float(getattr(tolerance, name)) for name in metrics}
    gates = {
        name: {
            "passed": metrics[name] <= thresholds[name],
            "value": metrics[name],
            "limit": thresholds[name],
        }
        for name in metrics
    }
    gates.update(
        {
            "final_time_reached": {
                "passed": abs(final_time - expected_time) <= time_tolerance,
                "value_s": final_time,
                "expected_s": expected_time,
            },
            "planned_time_axis": {
                "passed": True,
                "points": len(volume),
                "expected_points": len(expected_times),
            },
            "mean_energy_support": {
                "passed": (
                    observed_minimum >= support_minimum * (1.0 - 1.0e-9)
                    and required_support_maximum
                    <= support_maximum * (1.0 + 1.0e-9)
                ),
                "observed_eV": [observed_minimum, observed_maximum],
                "table_eV": [support_minimum, support_maximum],
                "relative_guard": support_guard,
                "required_upper_eV": required_support_maximum,
            },
            "physical_values": {
                "passed": all(
                    value > 0.0
                    for row in volume
                    for value in (row[1], row[2], row[4], row[5], row[6])
                )
                and all(value > 0.0 for value in late_means)
                and all(row[1] > 0.0 for row in minimum + maximum + coil),
            },
            "passive_absorbed_power": {
                "passed": final[6]
                <= coil[-1][1]
                * (1.0 + plan.mapping.run.convergence.coil_power_relative_error),
                "absorbed_power_W": final[6],
                "coil_power_W": coil[-1][1],
                "relative_allowance": (
                    plan.mapping.run.convergence.coil_power_relative_error
                ),
            },
        }
    )
    passed = all(bool(gate["passed"]) for gate in gates.values())
    result = {
        "schema": "swarm.gec_icp_convergence.v1",
        "status": "passed" if passed else "failed",
        "passed": passed,
        "source": plan.mapping.bundle.expected_source,
        "stationarity_contract": (
            "maximum adjacent relative change over the final four output times"
        ),
        "comparison_times_s": [row[0] for row in late_volume],
        "final_values": {
            "axisymmetric_volume_m3": final[1],
            "electron_inventory": final[2],
            "electron_weighted_mean_energy_eV": final_mean,
            "metastable_inventory": final[4],
            "ion_inventory": final[5],
            "absorbed_power_W": final[6],
            "coil_power_W": coil[-1][1],
            "minimum_mean_energy_eV": minimum[-1][1],
            "maximum_mean_energy_eV": maximum[-1][1],
        },
        "metrics": metrics,
        "thresholds": thresholds,
        "gates": gates,
    }
    write_json(output / "convergence.json", result)
    return result


def _numeric_rows(path: Path, *, columns: int) -> list[list[float]]:
    if not path.is_file():
        raise GecIcpWorkflowError(f"missing COMSOL convergence export: {path}")
    try:
        lines = [
            line
            for line in path.read_text(encoding="utf-8-sig").splitlines()
            if line.strip() and not line.lstrip().startswith("%")
        ]
        parsed = [
            next(csv.reader([line])) if "," in line else line.split() for line in lines
        ]
        rows = [[float(value) for value in row] for row in parsed]
    except (OSError, UnicodeError, csv.Error, ValueError) as exc:
        raise GecIcpWorkflowError(
            f"invalid numeric COMSOL convergence export: {path}"
        ) from exc
    if not rows or any(len(row) != columns for row in rows):
        raise GecIcpWorkflowError(f"incomplete COMSOL convergence export: {path}")
    if any(not math.isfinite(value) for row in rows for value in row):
        raise GecIcpWorkflowError(f"nonfinite COMSOL convergence value: {path}")
    if any(right[0] <= left[0] for left, right in zip(rows, rows[1:])):
        raise GecIcpWorkflowError(
            f"COMSOL convergence times are not increasing: {path}"
        )
    return rows


def _require_same_times(*tables: list[list[float]]) -> None:
    reference = [row[0] for row in tables[0]]
    for rows in tables[1:]:
        candidate = [row[0] for row in rows]
        if len(candidate) != len(reference) or any(
            not math.isclose(left, right, rel_tol=1.0e-10, abs_tol=1.0e-15)
            for left, right in zip(reference, candidate, strict=True)
        ):
            raise GecIcpWorkflowError(
                "COMSOL convergence exports do not share one time axis"
            )


def _require_expected_times(
    rows: list[list[float]], expected: tuple[float, ...]
) -> None:
    observed = [row[0] for row in rows]
    if len(observed) != len(expected) or any(
        not math.isclose(left, right, rel_tol=1.0e-10, abs_tol=1.0e-15)
        for left, right in zip(observed, expected, strict=True)
    ):
        raise GecIcpWorkflowError(
            "COMSOL convergence export does not match the planned time axis"
        )


def _ratio(numerator: float, denominator: float, name: str) -> float:
    if denominator <= 0.0:
        raise GecIcpWorkflowError(f"cannot evaluate {name} from zero inventory")
    value = numerator / denominator
    if not math.isfinite(value) or value <= 0.0:
        raise GecIcpWorkflowError(f"invalid {name}: {value}")
    return value


def _relative_change(reference: float, candidate: float) -> float:
    return abs(candidate - reference) / max(abs(reference), abs(candidate), 1.0e-300)


def _maximum_adjacent_relative_change(values: list[float]) -> float:
    if len(values) < 2:
        raise GecIcpWorkflowError("stationarity window requires at least two values")
    return max(_relative_change(left, right) for left, right in zip(values, values[1:]))


__all__ = ["assess_gec_icp_convergence"]
