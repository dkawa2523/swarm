"""Materialization stages for a validated GEC-CCP COMSOL run plan."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
from pathlib import Path
from typing import Any

from swarm_workflow._io import write_json
from swarm_workflow.comsol.input.eedf_audit import (
    ComsolEedfAuditPlan,
    prepare_comsol_eedf_audit,
)
from swarm_workflow.comsol.models.gec_ccp.closure import (
    _reaction_handling_metadata,
    _transport_input_contract,
    _transport_tensor_metadata,
    _uses_function_eedf,
    _uses_transport,
)
from swarm_workflow.comsol.models.gec_ccp.contracts import (
    GecCcpMapping,
    GecCcpPlan,
    GecCcpWorkflowError,
    MphContract,
)
from swarm_workflow.comsol.models.gec_ccp.java import (
    APPLY_CLASS,
    BASELINE_EXPORT_CLASS,
    BASELINE_RUN_CLASS,
    EXTERNAL_EXPORT_CLASS,
    EXTERNAL_RUN_CLASS,
    NATIVE_EEDF_AUDIT_CLASS,
    generate_apply_java,
    generate_export_java,
    generate_run_java,
)
from swarm_workflow.comsol.models.gec_ccp.planning.evidence import ClosureEvidence


@dataclass(frozen=True, slots=True)
class PreparedGecCcpLayout:
    """Filesystem layout and immutable plan assembled before source emission."""

    plan: GecCcpPlan
    paths: dict[str, Path]
    external_directory: Path
    baseline_directory: Path | None


def prepare_layout(
    mapping: GecCcpMapping,
    contract: MphContract,
    *,
    write_java: bool,
) -> PreparedGecCcpLayout:
    output = mapping.results.output_directory
    output.mkdir(parents=True, exist_ok=True)
    if mapping.model.baseline_output_mph is not None:
        mapping.model.baseline_output_mph.parent.mkdir(parents=True, exist_ok=True)
    mapping.model.output_mph.parent.mkdir(parents=True, exist_ok=True)
    external_dir = output / "swarm_tables"
    external_dir.mkdir(parents=True, exist_ok=True)
    baseline_dir: Path | None = None
    paths: dict[str, Path] = {
        "apply": output / f"{APPLY_CLASS}.java",
        "external_run": output / f"{EXTERNAL_RUN_CLASS}.java",
        "external_export": output / f"{EXTERNAL_EXPORT_CLASS}.java",
    }
    if mapping.run.include_builtin_reference:
        baseline_dir = output / "builtin_druyvesteyn"
        baseline_dir.mkdir(parents=True, exist_ok=True)
        paths.update(
            {
                "baseline_run": output / f"{BASELINE_RUN_CLASS}.java",
                "baseline_export": output / f"{BASELINE_EXPORT_CLASS}.java",
            }
        )
    else:
        for class_name in (BASELINE_RUN_CLASS, BASELINE_EXPORT_CLASS):
            (output / f"{class_name}.java").unlink(missing_ok=True)
    native_eedf_audit = _prepare_native_eedf_audit(
        mapping,
        output=output,
        write_java=write_java,
    )
    if native_eedf_audit is not None:
        paths["native_eedf_audit"] = native_eedf_audit.java_path
    result_directories = [external_dir]
    if baseline_dir is not None:
        result_directories.append(baseline_dir)
    expected = tuple(
        directory / name
        for directory in result_directories
        for name in (
            "domain_period_average.csv",
            "axis_period_average.csv",
            "radial_period_average.csv",
            "axis_phase_resolved.csv",
            "domain_phase_closure.csv",
            "electrode_waveform.csv",
            "closure_phase.csv",
            "closure_phase_radial.csv",
            "conservation_volume.csv",
            "conservation_wall.csv",
            "conservation_terminal_power.csv",
        )
    )
    plan_json = output / "gec_ccp_plan.json"
    plan = GecCcpPlan(
        mapping=mapping,
        contract=contract,
        output_directory=output,
        plan_json=plan_json,
        baseline_run_java=paths.get("baseline_run"),
        baseline_export_java=paths.get("baseline_export"),
        apply_java=paths["apply"],
        external_run_java=paths["external_run"],
        external_export_java=paths["external_export"],
        native_eedf_audit=native_eedf_audit,
        expected_result_files=expected,
    )
    return PreparedGecCcpLayout(
        plan=plan,
        paths=paths,
        external_directory=external_dir,
        baseline_directory=baseline_dir,
    )


def _prepare_native_eedf_audit(
    mapping: GecCcpMapping,
    *,
    output: Path,
    write_java: bool,
) -> ComsolEedfAuditPlan | None:
    if not _uses_function_eedf(mapping.closure):
        return None
    function_spec = mapping.closure.function_eedf
    if function_spec is None:
        raise GecCcpWorkflowError("Function-EEDF closure lacks its import spec")
    return prepare_comsol_eedf_audit(
        model_path=mapping.model.output_mph,
        table_path=mapping.bundle.path / function_spec.table,
        output_directory=output / "native_function_eedf_audit",
        component=mapping.model.component,
        physics=mapping.model.physics,
        function_tag=function_spec.function_tag,
        class_name=NATIVE_EEDF_AUDIT_CLASS,
        moment_mean_count=9,
        require_model_exists=False,
        write_java=write_java,
    )


def materialize_java(
    layout: PreparedGecCcpLayout,
    *,
    preintegrated_rate_rows: list[dict[str, Any]] | None,
    write_java: bool,
) -> None:
    if not write_java:
        return
    plan = layout.plan
    mapping = plan.mapping
    if mapping.run.include_builtin_reference:
        _materialize_baseline_java(layout)
    plan.apply_java.write_text(
        generate_apply_java(
            mapping,
            input_mph=mapping.model.input_mph,
            preintegrated_rate_rows=preintegrated_rate_rows,
        ),
        encoding="utf-8",
    )
    plan.external_run_java.write_text(
        generate_run_java(
            mapping,
            run_role="external",
            class_name=EXTERNAL_RUN_CLASS,
            input_mph=mapping.model.output_mph,
            output_mph=mapping.model.output_mph,
            study=mapping.model.external_study,
            solution=mapping.model.external_solution,
            conversion_source_study=mapping.model.external_study,
            time_periodic_feature=mapping.model.external_time_periodic_feature,
            convert_periodic_solution=(
                mapping.run.external_waveform_dataset != mapping.run.phase_dataset
            ),
        ),
        encoding="utf-8",
    )
    plan.external_export_java.write_text(
        generate_export_java(
            mapping,
            class_name=EXTERNAL_EXPORT_CLASS,
            input_mph=mapping.model.output_mph,
            output_dir=layout.external_directory,
            period_dataset=mapping.run.external_period_dataset,
            waveform_dataset=mapping.run.external_waveform_dataset,
            include_external_closure_audit=True,
        ),
        encoding="utf-8",
    )


def _materialize_baseline_java(layout: PreparedGecCcpLayout) -> None:
    plan = layout.plan
    mapping = plan.mapping
    if (
        plan.baseline_run_java is None
        or plan.baseline_export_java is None
        or mapping.model.baseline_output_mph is None
        or mapping.run.period_dataset is None
        or mapping.run.baseline_waveform_dataset is None
        or layout.baseline_directory is None
    ):
        raise GecCcpWorkflowError("built-in reference configuration is incomplete")
    plan.baseline_run_java.write_text(
        generate_run_java(
            mapping,
            run_role="baseline",
            class_name=BASELINE_RUN_CLASS,
            input_mph=mapping.model.input_mph,
            output_mph=mapping.model.baseline_output_mph,
            time_periodic_feature=mapping.model.time_periodic_feature,
        ),
        encoding="utf-8",
    )
    plan.baseline_export_java.write_text(
        generate_export_java(
            mapping,
            class_name=BASELINE_EXPORT_CLASS,
            input_mph=mapping.model.baseline_output_mph,
            output_dir=layout.baseline_directory,
            period_dataset=mapping.run.period_dataset,
            waveform_dataset=mapping.run.baseline_waveform_dataset,
        ),
        encoding="utf-8",
    )


def write_plan_manifest(
    layout: PreparedGecCcpLayout,
    evidence: ClosureEvidence,
    *,
    write_java: bool,
) -> None:
    write_json(layout.plan.plan_json, _plan_payload(layout, evidence, write_java))


def _plan_payload(
    layout: PreparedGecCcpLayout,
    evidence: ClosureEvidence,
    write_java: bool,
) -> dict[str, Any]:
    plan = layout.plan
    mapping = plan.mapping
    bundle_summary = evidence.bundle
    nonlinear = _nonlinear_plan_metadata(mapping)
    return {
        "stage": "prepare-gec-ccp",
        "status": "ready",
        "result_role": mapping.results.role,
        "physical_target_eligible": mapping.results.role == "physical_target",
        "mapping": {
            "path": str(mapping.path),
            "sha256": hashlib.sha256(mapping.path.read_bytes()).hexdigest(),
        },
        "bundle": bundle_summary,
        "model": {
            "input_mph": {
                "path": str(mapping.model.input_mph),
                "sha256": hashlib.sha256(
                    mapping.model.input_mph.read_bytes()
                ).hexdigest(),
                "size_bytes": mapping.model.input_mph.stat().st_size,
            },
            "contract": asdict(plan.contract),
        },
        "closure": _closure_payload(mapping, evidence),
        "solve": _solve_payload(mapping, nonlinear),
        "bindings": _bindings_payload(mapping),
        "generated_java": {
            key: {
                "path": str(value),
                "materialized": bool(write_java and value.is_file()),
                "sha256": (
                    hashlib.sha256(value.read_bytes()).hexdigest()
                    if write_java and value.is_file()
                    else None
                ),
            }
            for key, value in layout.paths.items()
        },
        "native_function_eedf_audit": _native_eedf_audit_payload(
            plan.native_eedf_audit
        ),
        "expected_results": [str(path) for path in plan.expected_result_files],
    }


def _closure_payload(
    mapping: GecCcpMapping,
    evidence: ClosureEvidence,
) -> dict[str, Any]:
    bundle_summary = evidence.bundle
    return {
        "electron_transport": mapping.closure.electron_transport,
        "zero_field_isotropization_Td": mapping.closure.zero_field_isotropization_Td,
        "transport_tensor": _transport_tensor_metadata(
            mapping,
            source=str(bundle_summary.get("source")),
        ),
        "transport_inputs": _transport_input_contract(
            mapping.closure, source=str(bundle_summary.get("source"))
        ),
        "gradient_response_policy": mapping.closure.gradient_response_policy,
        "gradient_response": evidence.restricted_gradient_response,
        "thermal_diffusion_model": mapping.closure.thermal_diffusion_model,
        "reaction_model": mapping.closure.reaction_model,
        "external_rate_processes": list(mapping.closure.external_rate_processes),
        "elastic_energy_loss_model": mapping.closure.elastic_energy_loss_model,
        "elastic_energy_loss": bundle_summary.get("elastic_energy_loss"),
        "source_field": mapping.closure.source_field,
        "lookup_jacobian": mapping.closure.lookup_jacobian,
        "coefficient_representation": mapping.closure.interpolation,
        "data_processing": {
            "swarm_source_postprocess": bundle_summary["source_policy"]["postprocess"],
            "transport_coupling": (
                "positive_log_value_shape_preserving_piecewise_cubic_"
                "Hermite_C1_on_log_mean_energy"
                if _uses_transport(mapping.closure)
                else "not_applicable"
            ),
            "transport_extrapolation": "smooth_constant_asymptote",
            "function_eedf_projection": (
                "moment_preserving_nonnegative_native_2d_projection"
                if _uses_function_eedf(mapping.closure)
                else "not_applicable"
            ),
            "function_eedf_interpolation": (
                "COMSOL_native_2d_linear_C0"
                if _uses_function_eedf(mapping.closure)
                else "not_applicable"
            ),
            "function_eedf_extrapolation": (
                "constant" if _uses_function_eedf(mapping.closure) else "not_applicable"
            ),
        },
        "function_eedf": bundle_summary.get("function_eedf"),
        "upstream_eedf_evidence": bundle_summary.get("upstream_eedf_evidence"),
        "eedf_artifacts": bundle_summary.get("eedf_artifacts"),
        "two_term_joint_consistency": evidence.two_term_joint_consistency,
        "upstream_eedf_rate_consistency": evidence.upstream_eedf_rate_consistency,
        "reaction_handling": _reaction_handling_metadata(mapping),
        "preintegrated_rate_closure": evidence.preintegrated_rate_summary,
        "active_rate_support": evidence.active_rate_support,
        "mc_rate_censoring": evidence.mc_rate_censoring,
        "mc_rate_interpolation": evidence.mc_rate_interpolation,
        "active_closure_support": evidence.active_closure_support,
        "support_policy": mapping.run.support_policy,
        "low_energy_guard": evidence.low_energy_guard,
        "validation_mean_energy_floor_eV": (
            mapping.run.validation_mean_energy_floor_eV
        ),
        "wall_closure": {
            "owner": "input_COMSOL_GEC_model",
            "model": "WallDriftDiffusion_Te_thermal_velocity",
            "external_eedf_half_range_moments_consumed": False,
            "angular_half_range_response_identified": False,
            "scope": "shared_boundary_model_not_external_swarm_closure",
        },
    }


def _nonlinear_plan_metadata(mapping: GecCcpMapping) -> dict[str, Any]:
    if mapping.run.nonlinear_globalization == "double_dogleg":
        baseline_method = "COMSOL_double_dogleg_trust_region"
        external_method = "COMSOL_double_dogleg_trust_region"
    elif mapping.run.nonlinear_globalization == "automatic_newton_no_recovery":
        baseline_method = "model_defined_automatic_Newton_without_minimum_step_recovery"
        external_method = "model_defined_automatic_Newton_without_minimum_step_recovery"
    else:
        baseline_method = "model_defined"
        external_method = "model_defined"
    solve_roles = ["external"]
    if mapping.run.include_builtin_reference:
        solve_roles.append("baseline")
    globalization: dict[str, Any] = {
        "method": mapping.run.nonlinear_globalization,
        "scope": solve_roles,
        "residual_scaling": (
            "fieldwise"
            if mapping.run.nonlinear_globalization == "double_dogleg"
            else "model_defined"
        ),
        "physical_equations_changed": False,
    }
    if mapping.run.nonlinear_globalization == "automatic_newton_no_recovery":
        globalization.update(
            {
                "evaluation_role": "diagnostic",
                "equation_preserving": True,
                "minimum_step_recovery": "off",
                "all_other_nonlinear_settings": "model_defined",
            }
        )
    return {
        "baseline_method": baseline_method,
        "external_method": external_method,
        "solve_roles": solve_roles,
        "globalization": globalization,
    }


def _solve_payload(
    mapping: GecCcpMapping,
    nonlinear: dict[str, Any],
) -> dict[str, Any]:
    return {
        "power_W": mapping.run.power_W,
        "include_builtin_reference": mapping.run.include_builtin_reference,
        "native_physics_initial_values": True,
        "saved_solution_dependency": False,
        "power_sweep": False,
        "coefficient_sweep": False,
        **(
            {"baseline_nonlinear_method": nonlinear["baseline_method"]}
            if mapping.run.include_builtin_reference
            else {}
        ),
        "external_nonlinear_method": nonlinear["external_method"],
        "nonlinear_solver_overrides": (
            mapping.run.nonlinear_globalization != "model_defined"
        ),
        "nonlinear_globalization": nonlinear["globalization"],
        "external_model_is_separate_copy": True,
        "external_study_policy": (
            "reuse_native_study"
            if mapping.model.external_study == mapping.model.study
            else "dedicated_external_study"
        ),
        "periodic_solver_exact_identities": {
            "scope": nonlinear["solve_roles"],
            "equation_view_locks": {
                "ptp.Mn": "0.04[kg/mol]",
                "ptp.ebar": "exp(En_per-Ne_per)*1[V]",
                "ptp.Te": "2*ptp.ebar/3",
            },
            "reason": (
                "remove_exact_algebraic_cancellation_and_exp_ratio_"
                "underflow_during_the_periodic_nonlinear_solve"
            ),
            "applied_during_periodic_solve": True,
            "removed_before_conversion": True,
            "removed_before_final_save": True,
            "expected_saved_lock_count": 0,
            "physical_equations_changed": False,
        },
        "native_heavy_species_ownership": {
            "owner": "comsol_plasma_interface",
            "species": ["Ar", "Ar_1p"],
            "mass_fraction_mapping_override": False,
            "domain_weak_override": False,
            "boundary_weak_override": False,
            "custom_initialization": False,
            "surface_reaction": {
                "tag": "sr1",
                "formula": "Ar+=>Ar",
                "owner": "comsol_plasma_interface",
            },
            "time_periodic_to_time_dependent": "native",
        },
        "source_stabilization": mapping.run.source_stabilization,
        "reaction_source_stabilization": (mapping.run.reaction_source_stabilization),
    }


def _bindings_payload(mapping: GecCcpMapping) -> dict[str, Any]:
    return {
        "component": mapping.model.component,
        "physics": mapping.model.physics,
        "plasma_feature": mapping.model.plasma_feature,
        **(
            {"baseline_study": mapping.model.study}
            if mapping.run.include_builtin_reference
            else {}
        ),
        "external_study": mapping.model.external_study,
        "external_solution": mapping.model.external_solution,
        "conversion_study": mapping.model.conversion_study,
        "datasets": {
            **(
                {
                    "period_baseline": mapping.run.period_dataset,
                    "waveform_baseline": mapping.run.baseline_waveform_dataset,
                }
                if mapping.run.include_builtin_reference
                else {}
            ),
            "period_external": mapping.run.external_period_dataset,
            "phase": mapping.run.phase_dataset,
            "waveform_external": mapping.run.external_waveform_dataset,
        },
    }


def _native_eedf_audit_payload(
    audit: ComsolEedfAuditPlan | None,
) -> dict[str, Any] | None:
    if audit is None:
        return None
    return {
        "active_table": {
            "path": str(audit.table_path),
            "sha256": audit.table_sha256,
            "import_contract": audit.import_contract.import_settings(),
        },
        "query_csv": {
            "path": str(audit.query_path),
            "sha256": audit.query_sha256,
            "points": audit.point_count,
            "moment_means": audit.moment_mean_count,
        },
        "expected_runtime_outputs": [
            str(audit.values_path),
            str(audit.contract_path),
            str(audit.values_path.parent / "comsol_eedf_audit.json"),
        ],
    }
