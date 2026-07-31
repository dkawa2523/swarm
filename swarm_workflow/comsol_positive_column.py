"""End-to-end positive-column COMSOL workflow with external Swarm tables."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
import re
from typing import Any, Iterable

import yaml

from ._paths import discover_repo_root
from .comsol_adapter import (
    ComsolExecutionSummary,
    execute_apply_comsol,
    execute_generated_comsol_java,
)
from .comsol_java import (
    GAS_TEMPERATURE_LOG_PREFIX,
    MEAN_ENERGY_MODEL_LOG_PREFIX,
    MESH_ELEMENTS_LOG_PREFIX,
    PRESSURE_LOG_PREFIX,
    extract_probe_csv_from_stdout,
    generate_apply_java_source,
    generate_continuation_java_source,
    generate_probe_table_export_java_source,
)
from .comsol_mapping import (
    CLOSURE_QUANTITIES,
    ComsolMappingError,
    ComsolModelMapping,
    closure_quantity_activity,
    load_comsol_mapping,
    validate_comsol_mapping_files,
)
from .comsol_verify import (
    VerifyComsolSummary,
    execute_verify_comsol_functions,
    generate_verify_java_source,
    plan_verify_comsol_functions,
)


RESULTS_CSV = "positive_column_results.csv"
RUN_SUMMARY_JSON = "positive_column_run_summary.json"
RESULT_EXPORT_CLASS = "SwarmPositiveColumnResultsExport"
RUN_CLASS = "SwarmPositiveColumnRun"
DEFAULT_OUTPUT_DIR = "outputs/comsol_positive_column"
REQUIRED_TABLES = (
    "mean_energy_vs_en.csv",
    "transport_vs_en.csv",
    "transport_vs_mean_energy.csv",
    "rates_vs_mean_energy.csv",
    "quality.csv",
)
REQUIRED_TRANSPORT_COLUMNS = (
    "E_over_N_Td",
    "E_over_N_V_m2",
    "mean_energy_eV",
    "reduced_mobility_m2_V_s_m3",
    "reduced_diffusion_L_m2_s_m3",
    "reduced_electron_energy_mobility_m2_V_s_m3",
    "reduced_electron_energy_diffusion_m2_s_m3",
)
REQUIRED_RATE_COLUMNS = (
    "mean_energy_eV",
    "process_type",
    "reduced_townsend_m2",
)
REQUIRED_RESULT_PROBES = (
    "x",
    "electron_density",
    "mean_electron_energy",
    "electric_potential",
    "electron_current_density",
    "ion_current_density",
    "total_current_density",
    "excitation_source",
    "ionization_source",
    "E_over_N",
    "reduced_mobility",
    "reduced_diffusion_L",
    "excitation_townsend",
    "ionization_townsend",
    "applied_voltage",
    "gas_pressure",
)
REQUIRED_CLOSURE_ITEMS = CLOSURE_QUANTITIES


class PositiveColumnWorkflowError(RuntimeError):
    """Raised when the positive-column E2E workflow cannot proceed."""


@dataclass(frozen=True, slots=True)
class BundleValidationSummary:
    bundle_path: Path
    valid_e_over_n_Td: tuple[float, float] | None
    quality_passed: bool | None
    tables_checked: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class PositiveColumnModelConditions:
    pressure_Pa: float | None
    gas_temperature_K: float | None
    mesh_elements: int | None


@dataclass(frozen=True, slots=True)
class PositiveColumnPlan:
    mapping: ComsolModelMapping
    original_mapping_path: Path
    effective_mapping_path: Path
    bundle: BundleValidationSummary
    model_conditions: PositiveColumnModelConditions
    output_dir: Path
    apply_java_path: Path
    verify_java_path: Path
    run_java_path: Path
    result_java_path: Path
    result_csv_path: Path
    steps: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class PositiveColumnRunSummary:
    plan: PositiveColumnPlan
    apply: ComsolExecutionSummary
    verify: VerifyComsolSummary
    run: ComsolExecutionSummary
    export: ComsolExecutionSummary
    run_summary_json: Path


def prepare_positive_column_run(
    mapping_path: str | Path,
    *,
    bundle_path: str | Path,
    write_java: bool = True,
    pressure_Pa: float | None = None,
    gas_temperature_K: float | None = None,
    mesh_elements: int | None = None,
) -> PositiveColumnPlan:
    original = Path(mapping_path).resolve()
    output_dir = _default_output_dir(original)
    output_dir.mkdir(parents=True, exist_ok=True)
    effective_mapping_path = output_dir / "positive_column_effective_mapping.yaml"
    _write_effective_mapping(original, bundle_path, effective_mapping_path)
    mapping = load_comsol_mapping(effective_mapping_path)
    _validate_mapping_for_positive_column(mapping)
    bundle = validate_positive_column_bundle(mapping.bundle.path)
    model_conditions = _validate_model_conditions(
        pressure_Pa=pressure_Pa,
        gas_temperature_K=gas_temperature_K,
        mesh_elements=mesh_elements,
    )
    validate_comsol_mapping_files(mapping, require_unit_metadata=True)
    apply_java = output_dir / "SwarmComsolApply.java"
    verify_plan = plan_verify_comsol_functions(effective_mapping_path)
    _validate_verify_plan_for_positive_column(verify_plan)
    verify_java = output_dir / "SwarmComsolVerify.java"
    run_java = output_dir / f"{RUN_CLASS}.java"
    result_java = output_dir / f"{RESULT_EXPORT_CLASS}.java"
    result_csv = output_dir / RESULTS_CSV
    plan = PositiveColumnPlan(
        mapping=mapping,
        original_mapping_path=original,
        effective_mapping_path=effective_mapping_path,
        bundle=bundle,
        model_conditions=model_conditions,
        output_dir=output_dir,
        apply_java_path=apply_java,
        verify_java_path=verify_java,
        run_java_path=run_java,
        result_java_path=result_java,
        result_csv_path=result_csv,
        steps=(
            "validate bundle manifest and mapping",
            "apply COMSOL interpolation tables",
            "verify imported COMSOL interpolation functions",
            f"run target {mapping.run.voltages_V[-1]:g} V; fallback continuation "
            + " -> ".join(f"{value:g} V" for value in mapping.run.voltages_V),
            "export positive-column result fields",
        ),
    )
    if write_java:
        apply_java.write_text(
            generate_apply_java_source(
                mapping,
                pressure_Pa=model_conditions.pressure_Pa,
                gas_temperature_K=model_conditions.gas_temperature_K,
                mesh_elements=model_conditions.mesh_elements,
            ),
            encoding="utf-8",
        )
        verify_java.write_text(generate_verify_java_source(verify_plan), encoding="utf-8")
        run_java.write_text(
            generate_continuation_java_source(mapping, class_name=RUN_CLASS),
            encoding="utf-8",
        )
        result_java.write_text(generate_result_export_java(plan), encoding="utf-8")
    return plan


def execute_positive_column_run(
    mapping_path: str | Path,
    *,
    bundle_path: str | Path,
    comsol_executable: str | Path | None = None,
    pressure_Pa: float | None = None,
    gas_temperature_K: float | None = None,
    mesh_elements: int | None = None,
) -> PositiveColumnRunSummary:
    plan = prepare_positive_column_run(
        mapping_path,
        bundle_path=bundle_path,
        write_java=True,
        pressure_Pa=pressure_Pa,
        gas_temperature_K=gas_temperature_K,
        mesh_elements=mesh_elements,
    )
    apply_kwargs: dict[str, Any] = {}
    if plan.model_conditions.pressure_Pa is not None:
        apply_kwargs["pressure_Pa"] = plan.model_conditions.pressure_Pa
    if plan.model_conditions.gas_temperature_K is not None:
        apply_kwargs["gas_temperature_K"] = (
            plan.model_conditions.gas_temperature_K
        )
    if plan.model_conditions.mesh_elements is not None:
        apply_kwargs["mesh_elements"] = plan.model_conditions.mesh_elements
    apply_summary = execute_apply_comsol(
        plan.effective_mapping_path,
        comsol_executable=comsol_executable,
        **apply_kwargs,
    )
    condition_markers = validate_positive_column_apply_conditions(
        plan,
        apply_summary,
    )
    verify_summary = execute_verify_comsol_functions(
        plan.effective_mapping_path,
        comsol_executable=comsol_executable,
    )
    run_summary = execute_generated_comsol_java(
        plan.mapping,
        plan.run_java_path,
        operation="positive_column_run",
        comsol_executable=comsol_executable,
    )
    plan.result_csv_path.unlink(missing_ok=True)
    export_summary = execute_generated_comsol_java(
        plan.mapping,
        plan.result_java_path,
        operation="positive_column_export",
        comsol_executable=comsol_executable,
    )
    if not plan.result_csv_path.exists():
        extract_probe_csv_from_stdout(export_summary.stdout_paths[-1], plan.result_csv_path)
    if not plan.result_csv_path.exists():
        raise PositiveColumnWorkflowError(
            f"positive-column result CSV was not produced: {plan.result_csv_path}; "
            f"see log_dir: {export_summary.log_dir}"
        )
    validate_positive_column_result(plan)
    closure_verification = validate_positive_column_closure_verification(
        plan,
        verify_summary,
    )
    model_condition_verification = validate_positive_column_model_conditions(
        plan,
        condition_markers=condition_markers,
    )
    run_summary_json = write_positive_column_run_summary(
        plan,
        apply_summary=apply_summary,
        verify_summary=verify_summary,
        run_summary=run_summary,
        export_summary=export_summary,
        closure_verification=closure_verification,
        model_condition_verification=model_condition_verification,
    )
    return PositiveColumnRunSummary(
        plan=plan,
        apply=apply_summary,
        verify=verify_summary,
        run=run_summary,
        export=export_summary,
        run_summary_json=run_summary_json,
    )


def write_positive_column_run_summary(
    plan: PositiveColumnPlan,
    *,
    apply_summary: Any | None = None,
    verify_summary: Any | None = None,
    run_summary: Any | None = None,
    export_summary: Any | None = None,
    closure_verification: dict[str, Any] | None = None,
    model_condition_verification: dict[str, Any] | None = None,
) -> Path:
    """Write the result, closure checks, and measured workflow timing."""

    result_stats = summarize_positive_column_result_csv(plan.result_csv_path)
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
    executed_voltages = getattr(
        run_summary,
        "executed_voltage_sequence_V",
        None,
    )
    formulation = plan.mapping.closure.mean_energy_formulation
    formulation_activity = closure_quantity_activity(formulation)
    formulation_readback = None
    if isinstance(model_condition_verification, dict):
        value = model_condition_verification.get("mean_energy_formulation")
        if isinstance(value, dict):
            formulation_readback = value.get("readback_comsol_value")
    summary = {
        "stage": "run-comsol",
        "status": "completed",
        "mapping": str(plan.original_mapping_path),
        "effective_mapping": str(plan.effective_mapping_path),
        "bundle": str(plan.bundle.bundle_path),
        "input_mph": str(plan.mapping.model.input_mph),
        "output_mph": str(plan.mapping.model.output_mph),
        "study": plan.mapping.model.study,
        "voltages_V": list(plan.mapping.run.voltages_V),
        "executed_voltage_sequence_V": (
            list(executed_voltages) if executed_voltages is not None else None
        ),
        "final_voltage_V": plan.mapping.run.voltages_V[-1],
        "result_csv": str(plan.result_csv_path),
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
        "closure_checks": validate_positive_column_result(plan),
        "closure_verification": (
            closure_verification
            if closure_verification is not None
            else validate_positive_column_closure_verification(
                plan,
                verify_summary,
            )
        ),
        "model_conditions": (
            model_condition_verification
            if model_condition_verification is not None
            else {
                "status": (
                    "not_requested"
                    if all(
                        value is None
                        for value in _model_conditions_dict(
                            plan.model_conditions
                        ).values()
                    )
                    else "not_evaluated"
                ),
                "requested": _model_conditions_dict(plan.model_conditions),
            }
        ),
        "logs": {
            "apply": _maybe_path(getattr(apply_summary, "log_dir", None)),
            "verify_summary_csv": _maybe_path(getattr(verify_summary, "summary_csv", None)),
            "run": _maybe_path(getattr(run_summary, "log_dir", None)),
            "export": _maybe_path(getattr(export_summary, "log_dir", None)),
        },
        "provenance": {
            "apply": _maybe_path(getattr(apply_summary, "provenance_json", None)),
            "verify": _maybe_path(
                getattr(verify_execution, "provenance_json", None)
            ),
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
    path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return path


def summarize_positive_column_result_csv(path: str | Path) -> dict[str, Any]:
    result_path = Path(path)
    if not result_path.exists():
        raise PositiveColumnWorkflowError(f"positive-column result CSV does not exist: {result_path}")
    rows = _read_numeric_rows(result_path)
    stats: dict[str, Any] = {
        "rows": len(rows),
        "columns": list(rows[0].keys()) if rows else [],
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
        values = _finite_values(rows, column)
        if values:
            stats[column] = {
                "min": min(values),
                "max": max(values),
            }
    stats["negative_counts"] = {
        column: sum(1 for value in _finite_values(rows, column) if value < 0.0)
        for column in (
            "electron_density",
            "excitation_source",
            "ionization_source",
        )
        if column in stats["columns"]
    }
    stats["nan_or_nonfinite_count"] = _nonfinite_count(rows)
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
    rows = _read_numeric_rows(plan.result_csv_path)
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
        values = _finite_values(rows, "gas_pressure")
        if len(values) != len(rows):
            raise PositiveColumnWorkflowError(
                "gas_pressure is missing or nonfinite in the spatial result"
            )
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
            "mesh": (
                f"{plan.mapping.model.component}/mesh1/edg1/dis1:elemcount"
            ),
        },
        "checks": checks,
    }


def validate_positive_column_result(plan: PositiveColumnPlan) -> dict[str, float]:
    """Reject converged-looking results that did not use the external closure."""

    rows = _read_numeric_rows(plan.result_csv_path)
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
    invalid = [
        column
        for column in sorted(required)
        if len(_finite_values(rows, column)) != len(rows)
    ]
    if invalid:
        raise PositiveColumnWorkflowError(
            "positive-column result has missing or nonfinite closure values: "
            + ", ".join(invalid)
        )
    e_values = [row["E_over_N"] / 1.0e-21 for row in rows]
    lower, upper = plan.bundle.valid_e_over_n_Td
    solved_lower = min(e_values)
    solved_upper = max(e_values)
    if solved_upper > upper or solved_lower < 0.0 or (
        solved_lower < lower and lower > 0.05
    ):
        raise PositiveColumnWorkflowError(
            "solved E/N is outside the Swarm table range: "
            f"solved=[{solved_lower:.6g}, {solved_upper:.6g}] Td, "
            f"table=[{lower:.6g}, {upper:.6g}] Td"
        )
    transport_x, transport_mu = _read_lookup(
        plan.bundle.bundle_path / "transport_vs_mean_energy.csv",
        "mean_energy_eV",
        "reduced_mobility_m2_V_s_m3",
    )
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
        mean_energy_x, mean_energy_y = _read_lookup(
            plan.bundle.bundle_path / "mean_energy_vs_en.csv",
            "E_over_N_Td",
            "mean_energy_eV",
        )
        errors["mean_energy_max_relative_error"] = _max_lookup_error(
            rows,
            "E_over_N",
            "mean_electron_energy",
            [value * 1.0e-21 for value in mean_energy_x],
            mean_energy_y,
        )
    transport_x, transport_diffusion = _read_lookup(
        plan.bundle.bundle_path / "transport_vs_mean_energy.csv",
        "mean_energy_eV",
        "reduced_diffusion_L_m2_s_m3",
    )
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
        rate_x, rate_y = _read_lookup(
            plan.bundle.bundle_path / "rates_vs_mean_energy.csv",
            "mean_energy_eV",
            "reduced_townsend_m2",
            process_type=process,
        )
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
            "seven-item COMSOL closure verification summary is missing: "
            f"{summary_path}"
        )
    payload = _read_json(summary_path)
    rows = _read_numeric_rows(plan.result_csv_path)
    checks = validate_positive_column_result(plan)
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
        "reduced_longitudinal_diffusion": (
            "longitudinal_diffusion_max_relative_error"
        ),
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


def _read_lookup(
    path: Path,
    x_column: str,
    y_column: str,
    *,
    process_type: str | None = None,
) -> tuple[list[float], list[float]]:
    values: dict[float, float] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            if process_type is not None and row.get("process_type") != process_type:
                continue
            values[float(row[x_column])] = float(row[y_column])
    x_values = sorted(values)
    return x_values, [values[value] for value in x_values]


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


def validate_positive_column_bundle(bundle_path: str | Path) -> BundleValidationSummary:
    bundle = Path(bundle_path).resolve()
    manifest_path = bundle / "manifest.json"
    if not manifest_path.exists():
        raise PositiveColumnWorkflowError(
            f"COMSOL bundle manifest.json does not exist: {manifest_path}"
        )
    manifest = _read_json(manifest_path)
    if manifest.get("status") == "failed":
        raise PositiveColumnWorkflowError(f"COMSOL bundle status is failed: {manifest_path}")
    tables = manifest.get("tables", {})
    if not isinstance(tables, dict):
        raise PositiveColumnWorkflowError("COMSOL bundle manifest.tables must be a mapping")
    for table_name in REQUIRED_TABLES:
        if table_name not in tables:
            raise PositiveColumnWorkflowError(
                f"COMSOL bundle missing required table metadata: {table_name}"
            )
        if not (bundle / table_name).exists():
            raise PositiveColumnWorkflowError(
                f"COMSOL bundle missing required CSV: {bundle / table_name}"
            )
    _require_columns(tables, "transport_vs_en.csv", REQUIRED_TRANSPORT_COLUMNS)
    _require_columns(
        tables,
        "transport_vs_mean_energy.csv",
        REQUIRED_TRANSPORT_COLUMNS,
    )
    _require_columns(tables, "rates_vs_mean_energy.csv", REQUIRED_RATE_COLUMNS)
    quality_summary = manifest.get("quality_summary", {})
    quality_passed = None
    if isinstance(quality_summary, dict) and "passed" in quality_summary:
        quality_passed = bool(quality_summary["passed"])
        if not quality_passed:
            raise PositiveColumnWorkflowError(
                f"COMSOL bundle quality_summary.passed is false: {manifest_path}"
            )
    ranges = manifest.get("valid_ranges", {})
    valid_e = None
    if isinstance(ranges, dict) and isinstance(ranges.get("E_over_N_Td"), list):
        values = ranges["E_over_N_Td"]
        if len(values) == 2:
            valid_e = (float(values[0]), float(values[1]))
    if valid_e is None:
        raise PositiveColumnWorkflowError(
            f"COMSOL bundle valid_ranges.E_over_N_Td is missing: {manifest_path}"
        )
    return BundleValidationSummary(
        bundle_path=bundle,
        valid_e_over_n_Td=valid_e,
        quality_passed=quality_passed,
        tables_checked=REQUIRED_TABLES,
    )


def generate_result_export_java(plan: PositiveColumnPlan) -> str:
    return generate_probe_table_export_java_source(
        class_name=RESULT_EXPORT_CLASS,
        model_name="swarmPositiveColumnResults",
        input_mph=plan.mapping.model.output_mph,
        output_csv=plan.result_csv_path,
        eval_tag="swpc_eval",
        probes=(
            (probe.name, probe.expression, probe.unit)
            for probe in plan.mapping.probes
        ),
    )


def format_positive_column_plan(plan: PositiveColumnPlan) -> str:
    formulation = plan.mapping.closure.mean_energy_formulation
    activity = closure_quantity_activity(formulation)
    lines = [
        "COMSOL positive-column dry-run",
        f"mapping: {plan.original_mapping_path}",
        f"effective_mapping: {plan.effective_mapping_path}",
        f"bundle: {plan.bundle.bundle_path}",
        f"valid_E_over_N_Td: {plan.bundle.valid_e_over_n_Td}",
        f"quality_passed: {plan.bundle.quality_passed}",
        f"input_mph: {plan.mapping.model.input_mph}",
        f"output_mph: {plan.mapping.model.output_mph}",
        (
            "mean_energy_formulation: "
            f"{formulation.mode} ({formulation.property_group}/"
            f"{formulation.property}={formulation.comsol_value})"
        ),
        "active_closure_items: "
        + ", ".join(
            quantity
            for quantity, item in activity.items()
            if item["status"] == "active"
        ),
        "inactive_written_items: "
        + ", ".join(
            quantity
            for quantity, item in activity.items()
            if item["status"] == "inactive"
        ),
        "model_conditions: "
        + json.dumps(
            _model_conditions_dict(plan.model_conditions),
            ensure_ascii=False,
            sort_keys=True,
        ),
        f"logs: {plan.mapping.logs.path}",
        f"apply_java: {plan.apply_java_path}",
        f"verify_java: {plan.verify_java_path}",
        f"run_java: {plan.run_java_path}",
        f"result_java: {plan.result_java_path}",
        f"result_csv: {plan.result_csv_path}",
        "steps:",
    ]
    lines.extend(f"  - {step}" for step in plan.steps)
    lines.append("COMSOL command: not executed")
    return "\n".join(lines)


def format_positive_column_summary(summary: PositiveColumnRunSummary) -> str:
    return "\n".join(
        [
            "COMSOL positive-column workflow completed",
            f"output_mph: {summary.plan.mapping.model.output_mph}",
            f"result_csv: {summary.plan.result_csv_path}",
            f"run_summary: {summary.run_summary_json}",
            f"effective_mapping: {summary.plan.effective_mapping_path}",
            f"apply_log: {summary.apply.log_dir}",
            f"verify_summary: {summary.verify.summary_csv}",
            f"run_log: {summary.run.log_dir}",
            f"export_log: {summary.export.log_dir}",
            f"generated_java: {summary.plan.result_java_path}",
        ]
    )


def _validate_mapping_for_positive_column(mapping: ComsolModelMapping) -> None:
    if mapping.model.input_mph == mapping.model.output_mph:
        raise PositiveColumnWorkflowError(
            "model.output_mph must differ from model.input_mph; refusing to "
            f"overwrite template .mph: {mapping.model.input_mph}"
        )
    if not mapping.probes:
        raise PositiveColumnWorkflowError("mapping results.probes must not be empty")
    names = {probe.name for probe in mapping.probes}
    missing = [name for name in REQUIRED_RESULT_PROBES if name not in names]
    if missing:
        raise PositiveColumnWorkflowError(
            "mapping results.probes missing required positive-column outputs: "
            + ", ".join(missing)
        )
    function_tags = {function.tag for function in mapping.functions}
    if function_tags != {"sw_meanE", "sw_muN"}:
        raise PositiveColumnWorkflowError(
            "positive-column functions must be exactly sw_meanE and sw_muN"
        )
    process_types = {lookup.process_type for lookup in mapping.reaction_lookups}
    if process_types != {"excitation", "ionization"}:
        raise PositiveColumnWorkflowError(
            "positive-column reactions must be exactly excitation and ionization"
        )


def _validate_model_conditions(
    *,
    pressure_Pa: float | None,
    gas_temperature_K: float | None,
    mesh_elements: int | None,
) -> PositiveColumnModelConditions:
    for name, value in (
        ("pressure_Pa", pressure_Pa),
        ("gas_temperature_K", gas_temperature_K),
    ):
        if value is None:
            continue
        number = float(value)
        if not math.isfinite(number) or number <= 0.0:
            raise PositiveColumnWorkflowError(
                f"{name} must be finite and positive"
            )
    if mesh_elements is not None and mesh_elements not in (200, 400):
        raise PositiveColumnWorkflowError("mesh_elements must be 200 or 400")
    return PositiveColumnModelConditions(
        pressure_Pa=float(pressure_Pa) if pressure_Pa is not None else None,
        gas_temperature_K=(
            float(gas_temperature_K)
            if gas_temperature_K is not None
            else None
        ),
        mesh_elements=int(mesh_elements) if mesh_elements is not None else None,
    )


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


def _validate_verify_plan_for_positive_column(verify_plan: Any) -> None:
    quantities = {point.quantity for point in verify_plan.feature_points}
    missing = sorted(set(REQUIRED_CLOSURE_ITEMS) - quantities)
    if missing:
        raise PositiveColumnWorkflowError(
            "positive-column mapping cannot audit all seven written closure "
            "properties; missing: " + ", ".join(missing)
        )


def _write_effective_mapping(
    mapping_path: Path,
    bundle_path: str | Path,
    output_path: Path,
) -> None:
    raw = yaml.safe_load(mapping_path.read_text(encoding="utf-8")) or {}
    if not isinstance(raw, dict):
        raise ComsolMappingError("mapping YAML root must be a mapping")
    raw.setdefault("bundle", {})
    if not isinstance(raw["bundle"], dict):
        raise ComsolMappingError("bundle must be a mapping")
    raw["bundle"]["path"] = str(Path(bundle_path).resolve())
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")


def _default_output_dir(mapping_path: Path) -> Path:
    root = discover_repo_root(mapping_path)
    return (root / DEFAULT_OUTPUT_DIR).resolve()


def _require_columns(
    tables: dict[str, Any],
    table_name: str,
    columns: tuple[str, ...],
) -> None:
    table = tables.get(table_name, {})
    listed = table.get("columns") if isinstance(table, dict) else None
    if not isinstance(listed, list):
        raise PositiveColumnWorkflowError(
            f"COMSOL bundle table {table_name} has no manifest columns"
        )
    missing = [column for column in columns if column not in listed]
    if missing:
        raise PositiveColumnWorkflowError(
            f"COMSOL bundle table {table_name} missing columns: "
            + ", ".join(missing)
        )


def _read_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise PositiveColumnWorkflowError(f"JSON root must be a mapping: {path}")
    return data


def _read_numeric_rows(path: Path) -> list[dict[str, float | str | None]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    parsed: list[dict[str, float | str | None]] = []
    for row in rows:
        parsed_row: dict[str, float | str | None] = {}
        for key, value in row.items():
            if value is None or value == "":
                parsed_row[key] = None
                continue
            try:
                parsed_row[key] = float(value)
            except ValueError:
                parsed_row[key] = value
        parsed.append(parsed_row)
    return parsed


def _finite_values(rows: Iterable[dict[str, Any]], column: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        value = row.get(column)
        if isinstance(value, (float, int)) and math.isfinite(value):
            values.append(float(value))
    return values


def _nonfinite_count(rows: Iterable[dict[str, Any]]) -> int:
    count = 0
    for row in rows:
        for value in row.values():
            if isinstance(value, (float, int)) and not math.isfinite(value):
                count += 1
    return count


def _maybe_path(value: Any) -> str | None:
    return str(value) if value is not None else None


def _maybe_number(value: Any) -> float | None:
    if isinstance(value, (float, int)) and math.isfinite(value):
        return float(value)
    return None
