"""Post-solve closure, operating-condition, and canonical summary audits."""

from __future__ import annotations

from hashlib import sha256
from dataclasses import asdict
import json
import math
from pathlib import Path
import re
from typing import Any

from .lookup import read_lookup
from .results import RESULT_SCHEMA, PositiveColumnResults
from .config import CLOSURE_QUANTITIES, closure_quantity_activity
from .contracts import (
    PositiveColumnModelConditions,
    PositiveColumnPlan,
    PositiveColumnWorkflowError,
)
from .java import (
    GAS_TEMPERATURE_LOG_PREFIX,
    MEAN_ENERGY_MODEL_LOG_PREFIX,
    MESH_ELEMENTS_LOG_PREFIX,
    PRESSURE_LOG_PREFIX,
    VOLTAGE_LOG_PREFIX,
)

RUN_SUMMARY_JSON = "positive_column_run_summary.json"
REQUIRED_CLOSURE_ITEMS = CLOSURE_QUANTITIES


def write_run_status(
    path: Path,
    *,
    status: str,
    stage: str,
    mapping_path: str | Path,
    bundle_path: str | Path | None,
    error: Exception | None = None,
) -> None:
    """Replace prior acceptance as soon as a new execution is requested."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "stage": "run-positive-column",
        "status": status,
        "execution_stage": stage,
        "quality_passed": False,
        "mapping": str(Path(mapping_path).resolve()),
        "bundle": str(Path(bundle_path).resolve()) if bundle_path is not None else None,
    }
    if error is not None:
        payload["failure_reason"] = str(error)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def write_positive_column_run_summary(
    plan: PositiveColumnPlan,
    *,
    result: PositiveColumnResults,
    closure_checks: dict,
    closure_verification: dict,
    model_condition_verification: dict,
    apply_summary: Any | None = None,
    verify_summary: Any | None = None,
    run_summary: Any | None = None,
    export_summary: Any | None = None,
) -> Path:
    """Write the result, closure checks, and measured workflow timing."""

    missing = {field.name for field in plan.mapping.results.fields} - set(result.units)
    if missing:
        raise PositiveColumnWorkflowError(
            f"result CSV omitted configured output fields: {sorted(missing)}"
        )
    if (
        sha256(plan.result_csv_path.read_bytes()).hexdigest()
        != result.metadata["csv_sha256"]
    ):
        raise PositiveColumnWorkflowError("result CSV changed during acceptance")
    result_stats = _result_stats(result)
    timing = {
        "apply": _maybe_number(getattr(apply_summary, "total_time_s", None)),
        "verify": _maybe_number(
            getattr(getattr(verify_summary, "execution", None), "total_time_s", None)
        ),
        "run": _maybe_number(getattr(run_summary, "total_time_s", None)),
        "export": _maybe_number(getattr(export_summary, "total_time_s", None)),
    }
    measured = [value for value in timing.values() if value is not None]
    timing["total"] = sum(measured) if measured else None
    verify_execution = getattr(verify_summary, "execution", None)
    executed_voltages = _executed_voltage_sequence(run_summary)
    formulation = plan.mapping.closure.mean_energy_formulation
    formulation_activity = closure_quantity_activity(formulation)
    formulation_readback = None
    if isinstance(model_condition_verification, dict):
        value = model_condition_verification.get("mean_energy_formulation")
        if isinstance(value, dict):
            formulation_readback = value.get("readback_comsol_value")
    summary = {
        "stage": "run-positive-column",
        "artifact_schema": RESULT_SCHEMA,
        "status": "completed",
        "physical_context": asdict(plan.mapping.physical_context),
        "mapping": str(plan.original_mapping_path),
        "effective_mapping": str(plan.effective_mapping_path),
        "bundle": str(plan.bundle.bundle_path),
        "bundle_manifest_sha256": plan.bundle.manifest_sha256,
        "input_mph": str(plan.mapping.model.input_mph),
        "output_mph": str(plan.mapping.model.output_mph),
        "study": plan.mapping.model.study,
        "voltage_parameter": plan.mapping.run.voltage_parameter,
        "voltages_V": list(plan.mapping.run.voltages_V),
        "executed_voltage_sequence_V": (
            list(executed_voltages) if executed_voltages is not None else None
        ),
        "final_voltage_V": plan.mapping.run.voltages_V[-1],
        "result_selection": {
            "dataset": plan.mapping.results.dataset,
            "solutions": plan.mapping.results.solutions,
            "field_units": {
                field.name: field.unit for field in plan.mapping.results.fields
            },
            "time_unit": "s",
            "solution_parameter_unit": "parameter_unit column in result CSV",
        },
        "result_csv": plan.result_csv_path.name,
        "result_csv_sha256": result.metadata["csv_sha256"],
        "valid_E_over_N_Td": plan.bundle.valid_e_over_n_Td,
        "quality_passed": plan.bundle.quality_passed,
        "mean_energy_formulation": {
            "mode": formulation.mode,
            "property_path": (
                f"{plan.mapping.model.component}/{plan.mapping.model.physics}/"
                f"{formulation.property_group}/{formulation.property}"
            ),
            "expected_comsol_value": formulation.comsol_value,
            "readback_comsol_value": formulation_readback,
            "quantity_activity": formulation_activity,
            "active_items": [
                quantity
                for quantity, item in formulation_activity.items()
                if item["status"] == "active"
            ],
            "inactive_items": [
                quantity
                for quantity, item in formulation_activity.items()
                if item["status"] == "inactive"
            ],
        },
        "result_stats": result_stats,
        "closure_checks": closure_checks,
        "closure_verification": closure_verification,
        "model_conditions": model_condition_verification,
        "logs": {
            "apply": _maybe_path(getattr(apply_summary, "log_dir", None)),
            "verify_summary_csv": _maybe_path(
                getattr(verify_summary, "summary_csv", None)
            ),
            "run": _maybe_path(getattr(run_summary, "log_dir", None)),
            "export": _maybe_path(getattr(export_summary, "log_dir", None)),
        },
        "provenance": {
            "apply": _maybe_path(getattr(apply_summary, "provenance_json", None)),
            "verify": _maybe_path(getattr(verify_execution, "provenance_json", None)),
            "run": _maybe_path(getattr(run_summary, "provenance_json", None)),
            "export": _maybe_path(getattr(export_summary, "provenance_json", None)),
        },
        "comsol": {
            "version": getattr(run_summary, "comsol_version", None),
            "build": getattr(run_summary, "comsol_build", None),
        },
        "timing_s": timing,
        "note": "External closure passed; COMSOL agreement is reported separately.",
    }
    path = plan.output_dir / RUN_SUMMARY_JSON
    path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return path


def _result_stats(result: PositiveColumnResults) -> dict[str, Any]:
    rows = result.rows
    stats: dict[str, Any] = {
        "rows": len(rows),
        "columns": list(rows[0].keys()) if rows else [],
        "units": result.units,
    }
    for column in (
        "x",
        "E_over_N",
        "electron_density",
        "mean_electron_energy",
        "electric_potential",
        "electron_current_density",
        "ion_current_density",
        "total_current_density",
        "excitation_source",
        "ionization_source",
        "applied_voltage",
        "gas_pressure",
    ):
        if column in result.units:
            values = [row[column] for row in rows]
            stats[column] = {
                "min": min(values),
                "max": max(values),
            }
    if "E_over_N" in stats:
        stats["E_over_N_Td"] = {
            "min": stats["E_over_N"]["min"] / 1.0e-21,
            "max": stats["E_over_N"]["max"] / 1.0e-21,
        }
    return stats


def validate_positive_column_apply_conditions(
    plan: PositiveColumnPlan,
    apply_summary: Any,
) -> dict[str, Any]:
    """Require formulation and optional model-condition read-back markers."""

    requested = _model_conditions_dict(plan.model_conditions)
    stdout_paths = getattr(apply_summary, "stdout_paths", None)
    if not isinstance(stdout_paths, (tuple, list)) or not stdout_paths:
        raise PositiveColumnWorkflowError(
            "COMSOL apply summary has no stdout for formulation/condition verification"
        )
    texts: list[str] = []
    for value in stdout_paths:
        path = Path(value)
        if path.is_file():
            texts.append(path.read_text(encoding="utf-8", errors="replace"))
    text = "\n".join(texts)
    formulation = plan.mapping.closure.mean_energy_formulation
    configured_formulation = _parse_text_marker(
        text,
        MEAN_ENERGY_MODEL_LOG_PREFIX,
    )
    if configured_formulation != formulation.comsol_value:
        raise PositiveColumnWorkflowError(
            "COMSOL mean-energy formulation read-back does not match mapping: "
            f"expected {formulation.comsol_value!r}, got "
            f"{configured_formulation!r}"
        )
    formulation_check = {
        "status": "passed",
        "mode": formulation.mode,
        "property_path": (
            f"{plan.mapping.model.component}/{plan.mapping.model.physics}/"
            f"{formulation.property_group}/{formulation.property}"
        ),
        "expected_comsol_value": formulation.comsol_value,
        "readback_comsol_value": configured_formulation,
    }
    applied: dict[str, float | int] = {}
    marker_contract = (
        ("pressure_Pa", PRESSURE_LOG_PREFIX),
        ("gas_temperature_K", GAS_TEMPERATURE_LOG_PREFIX),
        ("mesh_elements", MESH_ELEMENTS_LOG_PREFIX),
    )
    for key, prefix in marker_contract:
        expected = requested[key]
        if expected is None:
            continue
        actual = _parse_condition_marker(text, prefix)
        if key == "mesh_elements":
            rounded = int(round(actual))
            if not math.isclose(actual, rounded, rel_tol=0.0, abs_tol=1.0e-9):
                raise PositiveColumnWorkflowError(
                    f"COMSOL mesh read-back is not an integer: {actual}"
                )
            applied[key] = rounded
            matches = rounded == int(expected)
        else:
            applied[key] = actual
            matches = math.isclose(
                actual,
                float(expected),
                rel_tol=1.0e-12,
                abs_tol=1.0e-12,
            )
        if not matches:
            raise PositiveColumnWorkflowError(
                f"COMSOL condition read-back {key}={actual} does not match "
                f"requested {expected}"
            )
    return {
        "status": "passed",
        "requested": requested,
        "applied": applied,
        "mean_energy_formulation": formulation_check,
    }


def validate_positive_column_model_conditions(
    plan: PositiveColumnPlan,
    *,
    condition_markers: dict[str, Any],
    result: PositiveColumnResults,
) -> dict[str, Any]:
    """Validate configured conditions against read-back and spatial output."""

    requested = _model_conditions_dict(plan.model_conditions)
    active = {key: value for key, value in requested.items() if value is not None}
    if condition_markers.get("status") != "passed":
        raise PositiveColumnWorkflowError(
            "COMSOL formulation/model-condition apply read-back did not pass"
        )
    formulation_check = condition_markers.get("mean_energy_formulation")
    if not isinstance(formulation_check, dict) or (
        formulation_check.get("status") != "passed"
    ):
        raise PositiveColumnWorkflowError(
            "COMSOL mean-energy formulation apply read-back did not pass"
        )
    if not active:
        return {
            "status": "not_requested",
            "requested": requested,
            "checks": {},
            "mean_energy_formulation": formulation_check,
        }
    rows = result.rows
    if not rows:
        raise PositiveColumnWorkflowError(
            "positive-column result CSV is empty during condition validation"
        )
    checks: dict[str, Any] = {
        "apply_feature_readback": {
            "status": "passed",
            "scope": "COMSOL apply-stage property read-back markers",
            "values": condition_markers["applied"],
        }
    }
    pressure = plan.model_conditions.pressure_Pa
    if pressure is not None:
        if "gas_pressure" not in result.units:
            raise PositiveColumnWorkflowError(
                "gas_pressure is missing in the spatial result"
            )
        values = [row["gas_pressure"] for row in rows]
        max_error = max(abs(value - pressure) for value in values)
        tolerance = max(1.0e-8, abs(pressure) * 1.0e-7)
        if max_error > tolerance:
            raise PositiveColumnWorkflowError(
                "spatial gas_pressure does not match requested pressure: "
                f"max_abs_error={max_error:.6g} Pa"
            )
        checks["pressure_spatial_profile"] = {
            "status": "passed",
            "scope": "all_exported_spatial_points",
            "control_points": len(values),
            "min_Pa": min(values),
            "max_Pa": max(values),
            "max_absolute_error_Pa": max_error,
        }
    if plan.model_conditions.gas_temperature_K is not None:
        checks["gas_temperature"] = {
            "status": "passed",
            "scope": (
                "COMSOL plas/pes1 T property read back immediately after "
                "setting; no stable spatial temperature probe is exported"
            ),
        }
    if plan.model_conditions.mesh_elements is not None:
        checks["mesh_elements"] = {
            "status": "passed",
            "scope": (
                "COMSOL comp1/mesh1/edg1/dis1 elemcount property read back "
                "after mesh rebuild"
            ),
        }
    return {
        "status": "passed",
        "requested": requested,
        "applied": condition_markers["applied"],
        "mean_energy_formulation": formulation_check,
        "model_tag_contract": {
            "pressure_parameter": "p0",
            "plasma_feature": (
                f"{plan.mapping.model.component}/{plan.mapping.model.physics}/"
                f"{plan.mapping.closure.feature}"
            ),
            "gas_temperature_property": "T",
            "mesh": (f"{plan.mapping.model.component}/mesh1/edg1/dis1:elemcount"),
        },
        "checks": checks,
    }


def validate_positive_column_result(
    plan: PositiveColumnPlan, result: PositiveColumnResults
) -> dict[str, float]:
    """Reject converged-looking results that did not use the external closure."""

    rows = result.rows
    if not rows:
        raise PositiveColumnWorkflowError("positive-column result CSV is empty")
    required = {
        "E_over_N",
        "mean_electron_energy",
        "reduced_mobility",
        "reduced_diffusion_L",
        "excitation_townsend",
        "ionization_townsend",
    }
    missing = sorted(required - set(rows[0]))
    if missing:
        raise PositiveColumnWorkflowError(
            "positive-column result is missing closure checks: " + ", ".join(missing)
        )
    e_values = [row["E_over_N"] / 1.0e-21 for row in rows]
    lower, upper = plan.bundle.valid_e_over_n_Td
    solved_lower = min(e_values)
    solved_upper = max(e_values)
    if (
        solved_upper > upper
        or solved_lower < 0.0
        or (solved_lower < lower and lower > 0.05)
    ):
        raise PositiveColumnWorkflowError(
            "solved E/N is outside the Swarm table range: "
            f"solved=[{solved_lower:.6g}, {solved_upper:.6g}] Td, "
            f"table=[{lower:.6g}, {upper:.6g}] Td"
        )
    transport_data = read_lookup(
        plan.bundle.bundle_path / "transport_vs_mean_energy.csv",
        "mean_energy_eV",
        "reduced_mobility_m2_V_s_m3",
    )
    transport_x, transport_mu = transport_data.x, transport_data.y
    mean_values = [row["mean_electron_energy"] for row in rows]
    if min(mean_values) < transport_x[0] or max(mean_values) > transport_x[-1]:
        raise PositiveColumnWorkflowError(
            "solved mean energy is outside the Swarm transport table range: "
            f"solved=[{min(mean_values):.6g}, {max(mean_values):.6g}] eV, "
            f"table=[{transport_x[0]:.6g}, {transport_x[-1]:.6g}] eV"
        )
    low_field_points = sum(value < lower for value in e_values)
    errors = {
        "low_field_clamp_from_Td": solved_lower if solved_lower < lower else 0.0,
        "low_field_clamp_points": float(low_field_points),
        "low_field_clamp_point_fraction": low_field_points / len(rows),
        "mobility_max_relative_error": _max_lookup_error(
            rows,
            "mean_electron_energy",
            "reduced_mobility",
            transport_x,
            transport_mu,
        ),
    }
    if plan.mapping.closure.mean_energy_formulation.mode == "local_field":
        mean_energy_data = read_lookup(
            plan.bundle.bundle_path / "mean_energy_vs_en.csv",
            "E_over_N_Td",
            "mean_energy_eV",
        )
        mean_energy_x, mean_energy_y = mean_energy_data.x, mean_energy_data.y
        errors["mean_energy_max_relative_error"] = _max_lookup_error(
            rows,
            "E_over_N",
            "mean_electron_energy",
            [value * 1.0e-21 for value in mean_energy_x],
            mean_energy_y,
        )
    diffusion_data = read_lookup(
        plan.bundle.bundle_path / "transport_vs_mean_energy.csv",
        "mean_energy_eV",
        "reduced_diffusion_L_m2_s_m3",
    )
    transport_x, transport_diffusion = diffusion_data.x, diffusion_data.y
    errors["longitudinal_diffusion_max_relative_error"] = _max_lookup_error(
        rows,
        "mean_electron_energy",
        "reduced_diffusion_L",
        transport_x,
        transport_diffusion,
    )
    for process, column in (
        ("excitation", "excitation_townsend"),
        ("ionization", "ionization_townsend"),
    ):
        lookup = next(
            item
            for item in plan.mapping.reaction_lookups
            if item.process_type == process
        )
        rate_data = read_lookup(
            lookup.path, lookup.x_column, lookup.y_column, reaction=lookup.identity
        )
        rate_x, rate_y = rate_data.x, rate_data.y
        errors[f"{process}_townsend_max_relative_error"] = _max_lookup_error(
            rows,
            "mean_electron_energy",
            column,
            rate_x,
            rate_y,
        )
    failed = {
        name: value
        for name, value in errors.items()
        if name.endswith("_relative_error") and value > 1.0e-2
    }
    if failed:
        detail = ", ".join(f"{name}={value:.3g}" for name, value in failed.items())
        raise PositiveColumnWorkflowError(
            "COMSOL solution did not follow the external closure: " + detail
        )
    return errors


def validate_positive_column_closure_verification(
    plan: PositiveColumnPlan,
    verify_summary: Any,
    *,
    result: PositiveColumnResults,
    closure_checks: dict,
) -> dict[str, Any]:
    """Audit seven written properties and report their mode-dependent activity."""

    summary_value = getattr(verify_summary, "summary_json", None)
    if summary_value is None:
        raise PositiveColumnWorkflowError(
            "seven-item COMSOL closure verification summary path is unset"
        )
    summary_path = Path(summary_value)
    if not summary_path.is_file():
        raise PositiveColumnWorkflowError(
            f"seven-item COMSOL closure verification summary is missing: {summary_path}"
        )
    payload = _read_json(summary_path)
    rows = result.rows
    checks = closure_checks
    raw_items = payload.get("closure_items")
    if not isinstance(raw_items, list):
        raise PositiveColumnWorkflowError(
            "COMSOL verification summary has no closure_items metadata"
        )
    by_quantity = {
        str(item.get("quantity")): item
        for item in raw_items
        if isinstance(item, dict) and item.get("quantity")
    }
    missing = sorted(set(REQUIRED_CLOSURE_ITEMS) - set(by_quantity))
    if missing:
        raise PositiveColumnWorkflowError(
            "COMSOL verification cannot verify required closure quantity/quantities: "
            + ", ".join(missing)
        )
    failed = [
        quantity
        for quantity in REQUIRED_CLOSURE_ITEMS
        if by_quantity[quantity].get("status") != "passed"
        or int(by_quantity[quantity].get("control_points", 0)) <= 0
    ]
    if failed:
        raise PositiveColumnWorkflowError(
            "COMSOL feature-table read-back failed for: " + ", ".join(failed)
        )
    formulation = plan.mapping.closure.mean_energy_formulation
    activity = closure_quantity_activity(formulation)
    profile_errors = {
        "reduced_mobility": "mobility_max_relative_error",
        "reduced_longitudinal_diffusion": ("longitudinal_diffusion_max_relative_error"),
        "excitation_townsend": "excitation_townsend_max_relative_error",
        "ionization_townsend": "ionization_townsend_max_relative_error",
    }
    if formulation.mode == "local_field":
        profile_errors["mean_energy"] = "mean_energy_max_relative_error"
    items: dict[str, Any] = {}
    for quantity in REQUIRED_CLOSURE_ITEMS:
        feature_item = by_quantity[quantity]
        quantity_activity = activity[quantity]
        item: dict[str, Any] = {
            "status": "passed",
            "activity": quantity_activity["status"],
            "activity_reason": quantity_activity["reason"],
            "feature_table_readback": {
                "scope": "all_written_table_rows",
                "control_points": int(feature_item["control_points"]),
                "max_relative_error": float(
                    feature_item.get("max_relative_error", 0.0)
                ),
            },
        }
        error_key = profile_errors.get(quantity)
        if error_key is not None:
            item["spatial_profile_check"] = {
                "scope": "all_exported_spatial_points",
                "control_points": len(rows),
                "max_relative_error": checks[error_key],
            }
        elif quantity_activity["status"] == "inactive":
            item["spatial_profile_check"] = {
                "scope": "not_applicable",
                "reason": quantity_activity["reason"],
            }
        else:
            item["spatial_profile_check"] = {
                "scope": "not_available",
                "reason": (
                    "COMSOL exposes no stable postprocessing expression for this "
                    "coefficient; the model feature table is read back directly."
                ),
            }
        items[quantity] = item
    active_items = [
        quantity
        for quantity in REQUIRED_CLOSURE_ITEMS
        if activity[quantity]["status"] == "active"
    ]
    inactive_items = [
        quantity
        for quantity in REQUIRED_CLOSURE_ITEMS
        if activity[quantity]["status"] == "inactive"
    ]
    return {
        "status": "passed",
        "mean_energy_formulation": {
            "mode": formulation.mode,
            "property_group": formulation.property_group,
            "property": formulation.property,
            "comsol_value": formulation.comsol_value,
        },
        "written_property_audit": {
            "status": "passed",
            "required_items": list(REQUIRED_CLOSURE_ITEMS),
        },
        "active_items": active_items,
        "inactive_items": inactive_items,
        "required_items": list(REQUIRED_CLOSURE_ITEMS),
        "items": items,
    }


def _max_lookup_error(
    rows: list[dict[str, float]],
    argument: str,
    actual: str,
    x_values: list[float],
    y_values: list[float],
) -> float:
    errors = []
    for row in rows:
        expected = _interpolate(x_values, y_values, row[argument])
        scale = max(abs(expected), 1.0e-30)
        errors.append(abs(row[actual] - expected) / scale)
    return max(errors)


def _interpolate(x_values: list[float], y_values: list[float], value: float) -> float:
    if value <= x_values[0]:
        return y_values[0]
    if value >= x_values[-1]:
        return y_values[-1]
    for index in range(1, len(x_values)):
        if value <= x_values[index]:
            fraction = (value - x_values[index - 1]) / (
                x_values[index] - x_values[index - 1]
            )
            return y_values[index - 1] + fraction * (
                y_values[index] - y_values[index - 1]
            )
    raise AssertionError("unreachable interpolation interval")


def _model_conditions_dict(
    conditions: PositiveColumnModelConditions,
) -> dict[str, float | int | None]:
    return {
        "pressure_Pa": conditions.pressure_Pa,
        "gas_temperature_K": conditions.gas_temperature_K,
        "mesh_elements": conditions.mesh_elements,
    }


def _parse_condition_marker(text: str, prefix: str) -> float:
    match = re.search(
        re.escape(prefix) + r"\s*([0-9eE.+-]+)",
        text,
    )
    if match is None:
        raise PositiveColumnWorkflowError(
            f"required COMSOL model-condition marker is missing: {prefix}"
        )
    try:
        value = float(match.group(1))
    except ValueError as exc:
        raise PositiveColumnWorkflowError(
            f"COMSOL model-condition marker is not numeric: {prefix}"
        ) from exc
    if not math.isfinite(value):
        raise PositiveColumnWorkflowError(
            f"COMSOL model-condition marker is nonfinite: {prefix}"
        )
    return value


def _parse_text_marker(text: str, prefix: str) -> str:
    match = re.search(re.escape(prefix) + r"\s*([^\r\n]+)", text)
    if match is None:
        raise PositiveColumnWorkflowError(
            f"required COMSOL formulation marker is missing: {prefix}"
        )
    value = match.group(1).strip()
    if not value:
        raise PositiveColumnWorkflowError(
            f"COMSOL formulation marker is empty: {prefix}"
        )
    return value


def _read_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise PositiveColumnWorkflowError(f"JSON root must be a mapping: {path}")
    return data


def _maybe_path(value: Any) -> str | None:
    return str(value) if value is not None else None


def _maybe_number(value: Any) -> float | None:
    if isinstance(value, (float, int)) and math.isfinite(value):
        return float(value)
    return None


def _executed_voltage_sequence(summary: Any | None) -> tuple[float, ...] | None:
    """Read positive-column continuation markers from a run log."""

    stdout_paths = getattr(summary, "stdout_paths", ())
    if not stdout_paths:
        return None
    path = Path(stdout_paths[-1])
    if not path.is_file():
        return None
    stdout = path.read_text(encoding="utf-8", errors="replace")
    executed = tuple(
        float(value)
        for value in re.findall(
            re.escape(VOLTAGE_LOG_PREFIX)
            + r"\s*([+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)"
            + r"(?:[eE][+-]?[0-9]+)?)\s*\[V\]",
            stdout,
        )
    )
    return executed or None
