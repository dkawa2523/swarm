from __future__ import annotations

import ast
from dataclasses import asdict, replace
from pathlib import Path
from types import ModuleType

import pytest

import swarm_workflow.comsol.models.gec_ccp as facade
import swarm_workflow.comsol.models.gec_ccp.validation.bundle as bundle_validation
import swarm_workflow.comsol.models.gec_ccp.validation.bundle_context as bundle_context
import swarm_workflow.comsol.models.gec_ccp.validation.bundle_guards as bundle_guards
import swarm_workflow.comsol.models.gec_ccp.closure as closure_contract
import swarm_workflow.comsol.models.gec_ccp.closure_arguments as closure_arguments
import swarm_workflow.comsol.models.gec_ccp.config_values as config_values
import swarm_workflow.comsol.models.gec_ccp.audits.conservation as conservation_audit
import swarm_workflow.comsol.models.gec_ccp.contracts as contracts
import swarm_workflow.comsol.models.gec_ccp.data as data_helpers
import swarm_workflow.comsol.models.gec_ccp.validation.eedf as eedf_validation
import swarm_workflow.comsol.models.gec_ccp.validation.eedf_artifacts as eedf_artifacts
import swarm_workflow.comsol.models.gec_ccp.validation.elastic_energy as elastic_energy
import swarm_workflow.comsol.models.gec_ccp.execution as execution
import swarm_workflow.comsol.models.gec_ccp.execution.context as execution_context
import swarm_workflow.comsol.models.gec_ccp.execution.inputs as execution_inputs
import swarm_workflow.comsol.models.gec_ccp.execution.pipeline as execution_pipeline
import swarm_workflow.comsol.models.gec_ccp.execution.postsolve as execution_postsolve
import swarm_workflow.comsol.models.gec_ccp.execution.preflight as execution_preflight
import swarm_workflow.comsol.models.gec_ccp.execution.runtime as execution_runtime
import swarm_workflow.comsol.models.gec_ccp.execution.solver as execution_solver
import swarm_workflow.comsol.models.gec_ccp.execution.status as execution_status
import swarm_workflow.comsol.models.gec_ccp.audits.function_eedf as function_run_audit
import swarm_workflow.comsol.models.gec_ccp.java as java_generation
import swarm_workflow.comsol.models.gec_ccp.java_apply_sections as java_apply_sections
import swarm_workflow.comsol.models.gec_ccp.java_export_sections as java_export_sections
import swarm_workflow.comsol.models.gec_ccp.java_support as java_support
import swarm_workflow.comsol.models.gec_ccp.validation.joint_consistency as joint_consistency
import swarm_workflow.comsol.models.gec_ccp.mapping as mapping_parser
import swarm_workflow.comsol.models.gec_ccp.mapping_sections as mapping_sections
import swarm_workflow.comsol.models.gec_ccp.validation.mc_bundle as mc_bundle_validation
import swarm_workflow.comsol.models.gec_ccp.validation.mc_quality as mc_quality
import swarm_workflow.comsol.models.gec_ccp.mph as mph_contract
import swarm_workflow.comsol.models.gec_ccp.planning.evidence as plan_evidence
import swarm_workflow.comsol.models.gec_ccp.prepare as preparation
import swarm_workflow.comsol.models.gec_ccp.prepare_stages as prepare_stages
import swarm_workflow.comsol.models.gec_ccp.validation.provenance as bundle_provenance
import swarm_workflow.comsol.models.gec_ccp.validation.quality as bundle_quality
import swarm_workflow.comsol.models.gec_ccp.audits.saved_model as saved_model_audit
import swarm_workflow.comsol.models.gec_ccp.validation.table_inputs as table_inputs
import swarm_workflow.comsol.models.gec_ccp.audits.transport_run as transport_run_audit
import swarm_workflow.comsol.models.gec_ccp.audits.transport_support as transport_support
import swarm_workflow.comsol.models.gec_ccp.audits.transport_values as transport_values


PRODUCTION_MAPS = (
    "argon_gec_ccp_two_term_function_eedf.yaml",
    "argon_gec_ccp_monte_carlo_function_eedf_restricted_lmea.yaml",
    "argon_gec_ccp_propagator_function_eedf.yaml",
)
PUBLIC_FACADE_NAMES = {
    "GecCcpPlan",
    "GecCcpQualityError",
    "GecCcpRunSummary",
    "GecCcpWorkflowError",
    "execute_gec_ccp_run",
    "format_gec_ccp_plan",
    "format_gec_ccp_summary",
    "prepare_gec_ccp_run",
}
MOVED_OWNER_FUNCTIONS = {
    "validate_manifest": bundle_provenance,
    "validate_active_tables": table_inputs,
    "validate_elastic_energy_loss": elastic_energy,
    "validate_upstream_eedf_evidence": eedf_artifacts,
    "validate_active_function_eedf": eedf_artifacts,
    "validate_bundle_quality": bundle_quality,
    "build_closure_evidence": plan_evidence,
    "prepare_gec_ccp_run": preparation,
    "execute_gec_ccp_run": execution_pipeline,
    "prepare_execution": execution_preflight,
    "run_solver_steps": execution_solver,
    "validate_execution_outputs": execution_runtime,
    "run_postsolve_audits": execution_postsolve,
    "evaluate_run_quality": execution_postsolve,
    "write_final_status": execution_status,
    "enforce_quality_acceptance": execution_status,
    "status_contract_metadata": execution_status,
    "canonical_result_artifacts": execution_status,
    "validate_gec_plan_inputs": execution_inputs,
    "_comsol_major_minor": execution_runtime,
    "_validated_comsol_runtime": execution_runtime,
    "_validated_saved_mph_versions": execution_runtime,
    "_comsol_blocker_code": execution_solver,
    "audit_gec_ccp_conservation_run": conservation_audit,
    "_audit_domain_phase_state": conservation_audit,
    "_conservation_balance": conservation_audit,
    "_read_comsol_numeric_csv": conservation_audit,
    "_read_fixed_comsol_table": conservation_audit,
    "audit_gec_ccp_function_eedf_run": function_run_audit,
    "_function_eedf_rate_range_quality": function_run_audit,
    "_audit_function_eedf_rates": function_run_audit,
    "audit_gec_ccp_closure_support_run": transport_run_audit,
    "audit_gec_ccp_transport_run": transport_run_audit,
    "_phase_expression_columns": transport_values,
    "_independent_elastic_energy_loss_value_audit": transport_values,
    "_independent_rate_coefficient_values": transport_values,
    "_independent_external_rate_value_audit": transport_values,
    "_interpolate_domain_state_to_radial_cut": transport_values,
    "_independent_isotropic_transport_tensor": transport_values,
    "_independent_field_aligned_transport_tensor": transport_values,
    "_independent_transport_tensors": transport_values,
    "_independent_transport_value_audit": transport_values,
    "_independent_domain_transport_value_audit": transport_values,
    "_axisymmetric_node_volume_weights": transport_support,
    "_guard_log_pchip_values": transport_support,
    "_low_energy_guard_influence_audit": transport_support,
    "_operating_mean_energy_support_audit": transport_support,
    "_saved_sw_logeps_audit": saved_model_audit,
    "_audit_transport_binding": saved_model_audit,
    "_audit_saved_elastic_energy_loss_binding": saved_model_audit,
    "_audit_saved_reaction_handling": saved_model_audit,
    "generate_apply_java": java_generation,
    "build_apply_java_context": java_apply_sections,
    "render_closure_arguments": java_apply_sections,
    "render_function_eedf": java_apply_sections,
    "closure_argument_range": closure_arguments,
    "smooth_log_energy_argument": closure_arguments,
    "inline_interpolation_lines": java_support,
    "java_header": java_support,
    "data_export_lines": java_support,
    "numerical_table_export_lines": java_support,
    "build_export_java_context": java_export_sections,
    "render_field_exports": java_export_sections,
    "render_conservation_exports": java_export_sections,
    "prepare_layout": prepare_stages,
    "materialize_java": prepare_stages,
    "write_plan_manifest": prepare_stages,
    "generate_run_java": java_generation,
    "generate_export_java": java_generation,
    "_mc_rate_censoring_audit": mc_bundle_validation,
    "_mc_rate_interpolation_audit": mc_bundle_validation,
    "independent_bundle_quality_audit": mc_quality,
    "_reaction_cross_sections_from_model_xml": eedf_validation,
    "_audit_gec_argon_cross_section_identity": eedf_validation,
    "_rate_significance": eedf_validation,
    "_rate_significance_metadata": eedf_validation,
    "_audit_function_eedf_source_rates": eedf_validation,
    "_integrate_native_function_eedf_rate": eedf_validation,
    "_read_function_eedf": eedf_validation,
    "_audit_upstream_eedf_artifact": eedf_validation,
    "_native_function_eedf_moment_audit": eedf_validation,
    "_native_function_eedf_row": eedf_validation,
    "_build_dense_function_eedf_rate_closure": joint_consistency,
    "_build_upstream_eedf_rate_consistency_audit": joint_consistency,
    "_build_two_term_joint_consistency_audit": joint_consistency,
    "_independent_log_piecewise_cubic_transport": joint_consistency,
    "_active_closure_mean_energy_support": joint_consistency,
    "_validate_two_term_temporal_growth_transport_contract": bundle_guards,
    "_validate_low_energy_guard_bundle": bundle_guards,
}
MOVED_BUNDLE_CONSTANTS = {
    "GEC_PROVENANCE_HASH_KEYS": bundle_context,
    "GEC_OPTIONAL_PROVENANCE_HASH_KEYS": bundle_context,
    "REQUIRED_TABLE_COLUMNS": table_inputs,
    "REQUIRED_RATE_EVIDENCE_COLUMNS": table_inputs,
    "GEC_LOW_ENERGY_GUARD_RELATIVE_INFLUENCE_LIMIT": bundle_guards,
    "GEC_ELASTIC_ENERGY_LOSS_TABLE": bundle_guards,
    "GEC_ELASTIC_ENERGY_LOSS_COLUMN": bundle_guards,
    "REQUIRED_MC_QUALITY_COLUMNS": mc_bundle_validation,
    "FUNCTION_EEDF_RATE_SIGNIFICANCE_FRACTION": eedf_validation,
    "FUNCTION_EEDF_RATE_RELATIVE_ERROR_LIMIT": eedf_validation,
    "FUNCTION_EEDF_RATE_P95_RELATIVE_ERROR_LIMIT": eedf_validation,
    "FUNCTION_EEDF_RATE_NORMALIZED_RMSE_LIMIT": eedf_validation,
    "UPSTREAM_EEDF_RATE_AUDIT_TABLE": eedf_validation,
    "GEC_TWO_TERM_JOINT_TRANSPORT_RELATIVE_ERROR_LIMIT": joint_consistency,
    "GEC_TWO_TERM_TRANSPORT_PROCESS_TYPES": joint_consistency,
    "GEC_CLOSURE_SUPPORT_MARGIN_LOG_FRACTION": joint_consistency,
    "GEC_CLOSURE_SUPPORT_MARGIN_LOG_ABSOLUTE": joint_consistency,
}
MOVED_RUN_AUDIT_CONSTANTS = {
    "GEC_INDEPENDENT_TRANSPORT_RELATIVE_TOLERANCE": transport_values,
    "GEC_REACTION_PHYSICS_CONTRACT": saved_model_audit,
}


def _top_level_definitions(module: ModuleType) -> set[str]:
    tree = ast.parse(Path(module.__file__).read_text(encoding="utf-8"))
    return {
        node.name
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
    }


def _imported_modules(module: ModuleType) -> set[str]:
    tree = ast.parse(Path(module.__file__).read_text(encoding="utf-8"))
    current_package = (
        module.__name__
        if Path(module.__file__).name == "__init__.py"
        else module.__name__.rsplit(".", 1)[0]
    )
    names: set[str] = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.ImportFrom):
            continue
        if node.level == 0:
            if node.module is not None:
                names.add(node.module)
            continue
        package_parts = current_package.split(".")
        base = package_parts[: len(package_parts) - node.level + 1]
        if node.module is not None:
            names.add(".".join((*base, node.module)))
        else:
            names.update(".".join((*base, alias.name)) for alias in node.names)
    names.update(
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    )
    return names


def _function_line_count(function: object) -> int:
    module = __import__(function.__module__, fromlist=["*"])
    tree = ast.parse(Path(module.__file__).read_text(encoding="utf-8"))
    definition = next(
        node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == function.__name__
    )
    return definition.end_lineno - definition.lineno + 1


def test_gec_facade_has_only_high_level_public_composition() -> None:
    assert set(facade.__all__) == PUBLIC_FACADE_NAMES
    assert facade.prepare_gec_ccp_run is preparation.prepare_gec_ccp_run
    assert facade.execute_gec_ccp_run is execution.execute_gec_ccp_run
    assert facade.GecCcpPlan is contracts.GecCcpPlan
    assert facade.GecCcpRunSummary is contracts.GecCcpRunSummary
    assert facade.GecCcpWorkflowError is contracts.GecCcpWorkflowError
    assert facade.GecCcpQualityError is contracts.GecCcpQualityError
    assert _top_level_definitions(facade) == {
        "format_gec_ccp_plan",
        "format_gec_ccp_summary",
    }
    for obsolete_export in (
        "_validate_gec_plan_inputs",
        "_audit_transport_binding",
        "_read_comsol_numeric_csv",
        "generate_apply_java",
        "load_gec_ccp_mapping",
        "inspect_gec_ccp_mph",
        "GEC_PROVENANCE_HASH_KEYS",
    ):
        assert obsolete_export not in facade.__dict__
    assert set(execution.__all__) == {"execute_gec_ccp_run"}
    assert _top_level_definitions(execution) == set()


def test_moved_functions_have_one_concrete_owner() -> None:
    modules = (
        facade,
        bundle_guards,
        bundle_provenance,
        preparation,
        execution,
        execution_context,
        execution_inputs,
        execution_pipeline,
        execution_postsolve,
        execution_preflight,
        execution_runtime,
        execution_solver,
        execution_status,
        eedf_validation,
        eedf_artifacts,
        elastic_energy,
        conservation_audit,
        function_run_audit,
        java_generation,
        java_apply_sections,
        java_export_sections,
        java_support,
        closure_arguments,
        joint_consistency,
        mc_bundle_validation,
        mc_quality,
        plan_evidence,
        prepare_stages,
        bundle_quality,
        saved_model_audit,
        table_inputs,
        transport_run_audit,
        transport_support,
        transport_values,
    )
    definitions = {module: _top_level_definitions(module) for module in modules}
    for name, owner in MOVED_OWNER_FUNCTIONS.items():
        assert name in definitions[owner]
        assert all(
            name not in names
            for module, names in definitions.items()
            if module is not owner
        )
        if owner is not bundle_validation:
            assert name not in bundle_validation.__dict__
    for name, owner in MOVED_BUNDLE_CONSTANTS.items():
        assert name in owner.__dict__
        assert name not in bundle_validation.__dict__
    for name, owner in MOVED_RUN_AUDIT_CONSTANTS.items():
        assert name in owner.__dict__


def test_obsolete_flat_gec_modules_are_removed() -> None:
    root = Path(__file__).resolve().parents[1] / "swarm_workflow"
    assert not (root / "comsol_gec_ccp.py").exists()
    assert not any(root.glob("gec_ccp_*.py"))


def test_orchestrators_use_direct_owner_objects() -> None:
    assert execution.execute_gec_ccp_run is execution_pipeline.execute_gec_ccp_run
    assert (
        execution_pipeline.execute_gec_ccp_run.__globals__["prepare_execution"]
        is execution_preflight.prepare_execution
    )
    assert execution_preflight.prepare_execution.__globals__["prepare_gec_ccp_run"] is (
        preparation.prepare_gec_ccp_run
    )
    assert (
        execution_pipeline.execute_gec_ccp_run.__globals__["run_postsolve_audits"]
        is execution_postsolve.run_postsolve_audits
    )
    assert (
        execution_postsolve.run_postsolve_audits.__globals__[
            "audit_gec_ccp_transport_run"
        ]
        is transport_run_audit.audit_gec_ccp_transport_run
    )
    assert preparation.prepare_gec_ccp_run.__globals__["prepare_stages"] is (
        prepare_stages
    )
    assert prepare_stages.materialize_java.__globals__["generate_apply_java"] is (
        java_generation.generate_apply_java
    )
    assert java_generation.generate_export_java.__globals__[
        "java_export_sections"
    ] is java_export_sections
    assert (
        preparation.prepare_gec_ccp_run.__globals__["build_closure_evidence"]
        is plan_evidence.build_closure_evidence
    )
    assert (
        plan_evidence.build_closure_evidence.__globals__["validate_gec_ccp_bundle"]
        is bundle_validation.validate_gec_ccp_bundle
    )
    assert (
        plan_evidence.build_closure_evidence.__globals__[
            "_validate_low_energy_guard_bundle"
        ]
        is bundle_guards._validate_low_energy_guard_bundle
    )
    assert (
        execution_inputs.validate_gec_plan_inputs.__globals__[
            "_validate_low_energy_guard_bundle"
        ]
        is bundle_guards._validate_low_energy_guard_bundle
    )
    assert bundle_validation.validate_gec_ccp_bundle.__globals__["_provenance"] is (
        bundle_provenance
    )
    assert bundle_validation.validate_gec_ccp_bundle.__globals__["_table_inputs"] is (
        table_inputs
    )
    assert table_inputs.validate_active_tables.__globals__["_bundle_guards"] is (
        bundle_guards
    )
    assert (
        joint_consistency._build_two_term_joint_consistency_audit.__globals__[
            "_validate_two_term_temporal_growth_transport_contract"
        ]
        is bundle_guards._validate_two_term_temporal_growth_transport_contract
    )
    assert (
        plan_evidence.build_closure_evidence.__globals__["_mc_rate_censoring_audit"]
        is mc_bundle_validation._mc_rate_censoring_audit
    )
    assert (
        bundle_quality.validate_bundle_quality.__globals__[
            "independent_bundle_quality_audit"
        ]
        is mc_quality.independent_bundle_quality_audit
    )
    assert (
        plan_evidence.build_closure_evidence.__globals__[
            "_build_two_term_joint_consistency_audit"
        ]
        is joint_consistency._build_two_term_joint_consistency_audit
    )
    assert preparation.prepare_gec_ccp_run.__globals__["_validate_contract"] is (
        mph_contract._validate_contract
    )


def test_gec_layers_follow_the_dependency_direction() -> None:
    package = "swarm_workflow.comsol.models.gec_ccp"

    def imports_layer(module: ModuleType, *layers: str) -> bool:
        prefixes = tuple(f"{package}.{layer}" for layer in layers)
        return any(
            name == prefix or name.startswith(prefix + ".")
            for name in _imported_modules(module)
            for prefix in prefixes
        )

    domain_modules = (
        contracts,
        data_helpers,
        config_values,
        mapping_parser,
        mapping_sections,
        closure_contract,
        closure_arguments,
    )
    for module in domain_modules:
        assert not imports_layer(
            module,
            "validation",
            "audits",
            "prepare",
            "execution",
            "java",
            "plots",
        )

    validation_modules = (
        bundle_guards,
        bundle_validation,
        bundle_context,
        bundle_provenance,
        eedf_validation,
        eedf_artifacts,
        elastic_energy,
        joint_consistency,
        mc_bundle_validation,
        mc_quality,
        bundle_quality,
        table_inputs,
    )
    for module in validation_modules:
        assert not imports_layer(
            module,
            "audits",
            "prepare",
            "execution",
            "java",
            "plots",
        )

    for module in (
        java_generation,
        java_apply_sections,
        java_export_sections,
        java_support,
        mph_contract,
    ):
        assert not imports_layer(
            module,
            "audits",
            "prepare",
            "execution",
            "plots",
        )

    assert not imports_layer(preparation, "audits", "execution", "plots")
    assert not imports_layer(prepare_stages, "audits", "execution", "plots")
    assert not imports_layer(plan_evidence, "audits", "execution", "plots")
    for module in (
        conservation_audit,
        function_run_audit,
        saved_model_audit,
        transport_run_audit,
        transport_support,
        transport_values,
    ):
        assert not imports_layer(module, "prepare", "execution", "plots")
    for module in (
        execution,
        execution_context,
        execution_inputs,
        execution_pipeline,
        execution_postsolve,
        execution_preflight,
        execution_runtime,
        execution_solver,
        execution_status,
    ):
        assert not imports_layer(module, "plots")


def test_execution_pipeline_dependencies_are_acyclic() -> None:
    prefix = "swarm_workflow.comsol.models.gec_ccp.execution"

    def imports_execution(module: ModuleType, *names: str) -> bool:
        targets = {f"{prefix}.{name}" for name in names}
        return bool(_imported_modules(module).intersection(targets))

    assert not imports_execution(
        execution_context,
        "inputs",
        "postsolve",
        "preflight",
        "runtime",
        "solver",
        "status",
        "pipeline",
    )
    assert not imports_execution(
        execution_inputs,
        "postsolve",
        "preflight",
        "runtime",
        "solver",
        "status",
        "pipeline",
    )
    assert not imports_execution(
        execution_status,
        "inputs",
        "postsolve",
        "preflight",
        "runtime",
        "solver",
        "pipeline",
    )
    for module in (
        execution_preflight,
        execution_solver,
        execution_runtime,
        execution_postsolve,
    ):
        assert not imports_execution(module, "pipeline")


def test_validation_and_execution_orchestrators_remain_bounded() -> None:
    assert _function_line_count(bundle_validation.validate_gec_ccp_bundle) <= 150
    assert _function_line_count(execution_pipeline.execute_gec_ccp_run) <= 150
    assert _function_line_count(preparation.prepare_gec_ccp_run) <= 150
    assert _function_line_count(java_generation.generate_apply_java) <= 150
    assert _function_line_count(java_generation.generate_export_java) <= 150
    for stage in (
        execution_preflight.prepare_execution,
        execution_solver.run_solver_steps,
        execution_runtime.validate_execution_outputs,
        execution_postsolve.run_postsolve_audits,
        execution_postsolve.evaluate_run_quality,
        execution_status.write_final_status,
        execution_status.enforce_quality_acceptance,
        bundle_provenance.validate_manifest,
        table_inputs.validate_active_tables,
        elastic_energy.validate_elastic_energy_loss,
        eedf_artifacts.validate_upstream_eedf_evidence,
        eedf_artifacts.validate_active_function_eedf,
        bundle_quality.validate_bundle_quality,
        mc_quality.independent_bundle_quality_audit,
        saved_model_audit._audit_transport_binding,
        prepare_stages.prepare_layout,
        prepare_stages.materialize_java,
        prepare_stages.write_plan_manifest,
        java_apply_sections.build_apply_java_context,
        java_apply_sections.render_closure_arguments,
        java_apply_sections.render_function_eedf,
        java_export_sections.build_export_java_context,
        java_export_sections.render_field_exports,
        java_export_sections.render_conservation_exports,
    ):
        assert _function_line_count(stage) <= 250, stage.__name__


def test_only_package_root_composes_the_public_workflow_api() -> None:
    implementation_modules = (
        bundle_validation,
        bundle_context,
        bundle_guards,
        bundle_provenance,
        contracts,
        data_helpers,
        closure_contract,
        execution,
        execution_context,
        execution_inputs,
        execution_pipeline,
        execution_postsolve,
        execution_preflight,
        execution_runtime,
        execution_solver,
        execution_status,
        eedf_validation,
        eedf_artifacts,
        elastic_energy,
        conservation_audit,
        function_run_audit,
        java_generation,
        java_apply_sections,
        java_export_sections,
        java_support,
        closure_arguments,
        joint_consistency,
        mapping_parser,
        mapping_sections,
        mc_bundle_validation,
        mc_quality,
        mph_contract,
        plan_evidence,
        preparation,
        prepare_stages,
        bundle_quality,
        saved_model_audit,
        table_inputs,
        transport_run_audit,
        transport_support,
        transport_values,
    )
    for module in implementation_modules:
        assert "swarm_workflow.comsol.models.gec_ccp" not in _imported_modules(module)


@pytest.mark.parametrize("name", PRODUCTION_MAPS)
def test_production_mapping_contract_remains_stable(name: str) -> None:
    root = Path(__file__).resolve().parents[1]
    path = root / "comsol_modes" / "maps" / name

    mapping = mapping_parser.load_gec_ccp_mapping(path)
    payload = asdict(mapping)
    contract = mph_contract.inspect_gec_ccp_mph(mapping.model.input_mph)
    mph_contract._validate_contract(mapping, contract)
    transport_contract = closure_contract._transport_input_contract(
        mapping.closure,
        source=mapping.bundle.expected_source,
    )

    assert type(mapping) is contracts.GecCcpMapping
    assert type(contract) is contracts.MphContract
    assert transport_contract["muN"]["source"] == "external_swarm_table"
    assert set(payload) == {
        "path",
        "root",
        "model",
        "bundle",
        "reactions",
        "closure",
        "run",
        "results",
        "logs",
    }


def test_mapping_error_keeps_contract_exception_identity(tmp_path: Path) -> None:
    path = tmp_path / "mapping.yaml"
    path.write_text("schema_version: 1\n", encoding="utf-8")

    with pytest.raises(
        contracts.GecCcpWorkflowError,
        match=(
            "GEC CCP mapping requires schema_version: 2; migrate "
            "run.powers_W to the single direct-solve setting run.power_W"
        ),
    ):
        mapping_parser.load_gec_ccp_mapping(path)


def test_invalid_mph_contract_keeps_validation_error() -> None:
    root = Path(__file__).resolve().parents[1]
    mapping = mapping_parser.load_gec_ccp_mapping(
        root / "comsol_modes" / "maps" / "argon_gec_ccp_two_term_function_eedf.yaml"
    )
    contract = mph_contract.inspect_gec_ccp_mph(mapping.model.input_mph)

    with pytest.raises(
        contracts.GecCcpWorkflowError,
        match=("local MPH does not match the GEC CCP mapping: physics tag=invalid"),
    ):
        mph_contract._validate_contract(
            mapping,
            replace(contract, physics_tag="invalid"),
        )
