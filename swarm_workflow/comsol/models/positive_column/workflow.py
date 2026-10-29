"""End-to-end positive-column COMSOL workflow with external Swarm tables."""

from __future__ import annotations

from electron_swarm.core.yaml_input import read_yaml_mapping

from dataclasses import asdict, replace
import json
import math
from pathlib import Path
from typing import Any

import yaml

from ...runtime import (
    ComsolAdapterError,
    ComsolExecutionSummary,
    execute_comsol_java_source,
    execute_generated_comsol_java,
)
from .config import (
    CLOSURE_QUANTITIES,
    ComsolMappingError,
    ComsolModelMapping,
    closure_quantity_activity,
    comsol_run_context,
    load_comsol_mapping,
    validate_comsol_mapping_files,
)
from .java import (
    DEFAULT_CLASS_NAME,
    generate_apply_java_source,
    generate_continuation_java_source,
    generate_field_export_java_source,
)
from .verify import (
    execute_verify_comsol_functions,
    generate_verify_java_source,
    plan_verify_comsol_functions,
)
from .contracts import (
    PositiveColumnPlan,
    PositiveColumnRunSummary,
    PositiveColumnWorkflowError,
    PositiveColumnModelConditions,
)
from . import bundle as _bundle
from . import audit as _audit
from .results import read_mapped_results
from .result_units import canonical_unit


RESULTS_CSV = "positive_column_results.csv"
RESULT_EXPORT_CLASS = "SwarmPositiveColumnResultsExport"
RUN_CLASS = "SwarmPositiveColumnRun"
REQUIRED_RESULT_FIELDS = (
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


def execute_apply_comsol(
    mapping_path: str | Path,
    comsol_executable: str | Path | None = None,
    *,
    pressure_Pa: float | None = None,
    gas_temperature_K: float | None = None,
    mesh_elements: int | None = None,
) -> ComsolExecutionSummary:
    """Apply positive-column lookup tables to a fresh model copy."""

    mapping = load_comsol_mapping(mapping_path)
    if gas_temperature_K is None:
        gas_temperature_K = mapping.physical_context.gas_temperature_K
    else:
        mapping = replace(
            mapping,
            physical_context=replace(
                mapping.physical_context, gas_temperature_K=float(gas_temperature_K)
            ),
        )
    _bundle.validate_positive_column_bundle(
        mapping.bundle.path, mapping=mapping, gas_temperature_K=gas_temperature_K
    )
    validate_comsol_mapping_files(mapping, require_unit_metadata=True)
    if mapping.model.input_mph == mapping.model.output_mph:
        raise ComsolAdapterError(
            "mapping model.output_mph must differ from model.input_mph; "
            f"refusing to overwrite input .mph: {mapping.model.input_mph}"
        )
    source = generate_apply_java_source(
        mapping,
        pressure_Pa=pressure_Pa,
        gas_temperature_K=gas_temperature_K,
        mesh_elements=mesh_elements,
    )
    return execute_comsol_java_source(
        comsol_run_context(mapping),
        source,
        class_name=DEFAULT_CLASS_NAME,
        operation="apply",
        comsol_executable=comsol_executable,
        study=mapping.model.study,
        require_fresh_output=True,
    )


def prepare_positive_column_run(
    mapping_path: str | Path,
    *,
    bundle_path: str | Path | None = None,
    write_java: bool = True,
    pressure_Pa: float | None = None,
    gas_temperature_K: float | None = None,
    mesh_elements: int | None = None,
) -> PositiveColumnPlan:
    original = Path(mapping_path).resolve()
    original_mapping = load_comsol_mapping(original)
    _validate_mapping_for_positive_column(original_mapping)
    output_dir = original_mapping.results.output_directory
    output_dir.mkdir(parents=True, exist_ok=True)
    effective_mapping_path = output_dir / "positive_column_effective_mapping.yaml"
    _write_effective_mapping(
        original,
        bundle_path,
        effective_mapping_path,
        gas_temperature_K=gas_temperature_K,
    )
    mapping = load_comsol_mapping(effective_mapping_path)
    _validate_mapping_for_positive_column(mapping)
    model_conditions = _validate_model_conditions(
        pressure_Pa=pressure_Pa,
        gas_temperature_K=gas_temperature_K
        if gas_temperature_K is not None
        else mapping.physical_context.gas_temperature_K,
        mesh_elements=mesh_elements,
    )
    bundle = _bundle.validate_positive_column_bundle(
        mapping.bundle.path,
        mapping=mapping,
        gas_temperature_K=model_conditions.gas_temperature_K,
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
        verify_java.write_text(
            generate_verify_java_source(verify_plan), encoding="utf-8"
        )
        run_java.write_text(
            generate_continuation_java_source(mapping, class_name=RUN_CLASS),
            encoding="utf-8",
        )
        result_java.write_text(
            generate_field_export_java_source(
                mapping, class_name=RESULT_EXPORT_CLASS, output_csv=result_csv
            ),
            encoding="utf-8",
        )
    return plan


def execute_positive_column_run(
    mapping_path: str | Path,
    *,
    bundle_path: str | Path | None = None,
    comsol_executable: str | Path | None = None,
    pressure_Pa: float | None = None,
    gas_temperature_K: float | None = None,
    mesh_elements: int | None = None,
) -> PositiveColumnRunSummary:
    mapping = load_comsol_mapping(mapping_path)
    _validate_mapping_for_positive_column(mapping)
    status_path = mapping.results.output_directory / _audit.RUN_SUMMARY_JSON
    stage = "preflight"
    _audit.write_run_status(
        status_path,
        status="running",
        stage=stage,
        mapping_path=mapping_path,
        bundle_path=bundle_path,
    )
    try:
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
            apply_kwargs["gas_temperature_K"] = plan.model_conditions.gas_temperature_K
        if plan.model_conditions.mesh_elements is not None:
            apply_kwargs["mesh_elements"] = plan.model_conditions.mesh_elements
        _bundle.assert_bundle_unchanged(plan.bundle)
        stage = "apply"
        apply_summary = execute_apply_comsol(
            plan.effective_mapping_path,
            comsol_executable=comsol_executable,
            **apply_kwargs,
        )
        condition_markers = _audit.validate_positive_column_apply_conditions(
            plan,
            apply_summary,
        )
        _bundle.assert_bundle_unchanged(plan.bundle)
        stage = "verify"
        verify_summary = execute_verify_comsol_functions(
            plan.effective_mapping_path,
            comsol_executable=comsol_executable,
        )
        _bundle.assert_bundle_unchanged(plan.bundle)
        stage = "solve"
        run_summary = execute_generated_comsol_java(
            comsol_run_context(plan.mapping),
            plan.run_java_path,
            operation="positive_column_run",
            comsol_executable=comsol_executable,
            study=plan.mapping.model.study,
        )
        _bundle.assert_bundle_unchanged(plan.bundle)
        stage = "export"
        plan.result_csv_path.unlink(missing_ok=True)
        export_summary = execute_generated_comsol_java(
            comsol_run_context(plan.mapping),
            plan.result_java_path,
            operation="positive_column_export",
            comsol_executable=comsol_executable,
            study=plan.mapping.model.study,
        )
        if not plan.result_csv_path.exists():
            raise PositiveColumnWorkflowError(
                f"positive-column result CSV was not produced: {plan.result_csv_path}; "
                f"see log_dir: {export_summary.log_dir}"
            )
        _bundle.assert_bundle_unchanged(plan.bundle)
        stage = "acceptance"
        result = read_mapped_results(plan.result_csv_path, plan.mapping.results)
        closure_checks = _audit.validate_positive_column_result(plan, result)
        closure_verification = _audit.validate_positive_column_closure_verification(
            plan,
            verify_summary,
            result=result,
            closure_checks=closure_checks,
        )
        model_condition_verification = _audit.validate_positive_column_model_conditions(
            plan,
            condition_markers=condition_markers,
            result=result,
        )
        run_summary_json = _audit.write_positive_column_run_summary(
            plan,
            result=result,
            closure_checks=closure_checks,
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
    except Exception as exc:
        _audit.write_run_status(
            status_path,
            status="blocked" if stage == "preflight" else "failed",
            stage=stage,
            mapping_path=mapping_path,
            bundle_path=bundle_path,
            error=exc,
        )
        raise


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
            asdict(plan.model_conditions),
            ensure_ascii=False,
            sort_keys=True,
        ),
        f"logs: {plan.mapping.logs.path}",
        f"apply_java: {plan.apply_java_path}",
        f"verify_java: {plan.verify_java_path}",
        f"run_java: {plan.run_java_path}",
        f"result_java: {plan.result_java_path}",
        f"result_csv: {plan.result_csv_path}",
        "result_units: "
        + json.dumps(
            {
                field.name: {
                    "export": field.unit,
                    "canonical": canonical_unit(field.name, field.unit)[0],
                }
                for field in plan.mapping.results.fields
            }
        ),
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
    if mapping.results is None:
        raise PositiveColumnWorkflowError(
            "positive-column run requires explicit results settings"
        )
    names = {field.name for field in mapping.results.fields}
    missing = [name for name in REQUIRED_RESULT_FIELDS if name not in names]
    if missing:
        raise PositiveColumnWorkflowError(
            "mapping results.fields missing required positive-column outputs: "
            + ", ".join(missing)
        )
    function_tags = {function.tag for function in mapping.functions}
    if function_tags != {"sw_meanE", "sw_muN"}:
        raise PositiveColumnWorkflowError(
            "positive-column functions must be exactly sw_meanE and sw_muN"
        )
    process_types = {lookup.process_type for lookup in mapping.reaction_lookups}
    if (
        process_types != {"excitation", "ionization"}
        or len(mapping.reaction_lookups) != 2
    ):
        raise PositiveColumnWorkflowError(
            "positive-column reactions must be exactly one excitation and one ionization channel"
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
            raise PositiveColumnWorkflowError(f"{name} must be finite and positive")
    if mesh_elements is not None and mesh_elements not in (200, 400):
        raise PositiveColumnWorkflowError("mesh_elements must be 200 or 400")
    return PositiveColumnModelConditions(
        pressure_Pa=float(pressure_Pa) if pressure_Pa is not None else None,
        gas_temperature_K=(
            float(gas_temperature_K) if gas_temperature_K is not None else None
        ),
        mesh_elements=int(mesh_elements) if mesh_elements is not None else None,
    )


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
    bundle_path: str | Path | None,
    output_path: Path,
    *,
    gas_temperature_K: float | None = None,
) -> None:
    raw = read_yaml_mapping(mapping_path)
    if not isinstance(raw, dict):
        raise ComsolMappingError("mapping YAML root must be a mapping")
    raw.setdefault("bundle", {})
    if not isinstance(raw["bundle"], dict):
        raise ComsolMappingError("bundle must be a mapping")
    mapping = load_comsol_mapping(mapping_path)
    raw["model"]["input_mph"] = str(mapping.model.input_mph)
    raw["model"]["output_mph"] = str(mapping.model.output_mph)
    raw["bundle"]["path"] = str(
        Path(bundle_path).resolve() if bundle_path is not None else mapping.bundle.path
    )
    raw.setdefault("logs", {})["path"] = str(mapping.logs.path)
    raw.setdefault("verify", {})["output_path"] = str(mapping.verify.output_path)
    if mapping.results is not None:
        raw["results"]["output_directory"] = str(mapping.results.output_directory)
    if gas_temperature_K is not None:
        if not isinstance(raw.get("physical_context"), dict):
            raise ComsolMappingError(
                "schema v2 migration required: physical_context must be a mapping"
            )
        raw["physical_context"]["gas_temperature_K"] = gas_temperature_K
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(yaml.safe_dump(raw, sort_keys=False), encoding="utf-8")
