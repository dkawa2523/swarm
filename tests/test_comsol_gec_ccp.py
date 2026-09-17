from __future__ import annotations

import csv
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from types import SimpleNamespace
from zipfile import ZipFile

import numpy as np
import pytest

import electron_swarm.solvers.monte_carlo.evidence as mc_evidence
from electron_swarm.solvers.monte_carlo.source_identity import (
    monte_carlo_source_sha256,
)
import swarm_workflow.comsol.models.gec_ccp.closure as closure_contract
import swarm_workflow.comsol.models.gec_ccp.audits.conservation as conservation_audit
import swarm_workflow.comsol.models.gec_ccp.contracts as contracts
import swarm_workflow.comsol.models.gec_ccp.data as data_helpers
import swarm_workflow.comsol.models.gec_ccp.validation.eedf as eedf_validation
import swarm_workflow.comsol.models.gec_ccp.execution as execution
import swarm_workflow.comsol.models.gec_ccp.execution.postsolve as execution_postsolve
import swarm_workflow.comsol.models.gec_ccp.execution.preflight as execution_preflight
import swarm_workflow.comsol.models.gec_ccp.execution.runtime as execution_runtime
import swarm_workflow.comsol.models.gec_ccp.execution.solver as execution_solver
import swarm_workflow.comsol.models.gec_ccp.audits.function_eedf as function_run_audit
import swarm_workflow.comsol.models.gec_ccp.closure_arguments as closure_arguments
import swarm_workflow.comsol.models.gec_ccp.validation.joint_consistency as joint_consistency
import swarm_workflow.comsol.models.gec_ccp.validation.mc_bundle as mc_bundle_validation
import swarm_workflow.comsol.models.gec_ccp.validation.mc_quality as mc_quality
import swarm_workflow.comsol.models.gec_ccp.validation.provenance as provenance_validation
import swarm_workflow.comsol.models.gec_ccp.audits.saved_model as saved_model_audit
import swarm_workflow.comsol.models.gec_ccp.audits.transport_support as transport_support
import swarm_workflow.comsol.models.gec_ccp.audits.transport_values as transport_value_audits
import swarm_workflow.tables.contracts as table_contracts
from swarm_workflow.comsol.runtime import ComsolAdapterError
from swarm_workflow.quality.policy import QualityThresholds, quality_thresholds_payload
from swarm_workflow.comsol.models.gec_ccp.closure import (
    _transport_audit_components,
    _transport_property_tensors,
)
from swarm_workflow.comsol.models.gec_ccp.contracts import (
    GecCcpQualityError,
    GecCcpWorkflowError,
)
from swarm_workflow.comsol.models.gec_ccp.execution.inputs import (
    validate_gec_plan_inputs,
)
from swarm_workflow.comsol.models.gec_ccp.java import (
    generate_apply_java,
    generate_run_java,
)
from swarm_workflow.comsol.models.gec_ccp.mapping import load_gec_ccp_mapping
from swarm_workflow.comsol.models.gec_ccp.mph import inspect_gec_ccp_mph
from swarm_workflow.comsol.models.gec_ccp.prepare import prepare_gec_ccp_run
from swarm_workflow.comsol.models.gec_ccp.audits.saved_model import (
    _audit_transport_binding,
)
from swarm_workflow.comsol.input.function_eedf import (
    piecewise_linear_weighted_moments,
)
from swarm_workflow.quality.solver import PropagatorTargetRequirement


@pytest.mark.parametrize(
    ("name", "source", "estimator_schema", "reaction_model"),
    (
        (
            "argon_gec_ccp_two_term_function_eedf.yaml",
            "two_term",
            None,
            "function_eedf",
        ),
        (
            "argon_gec_ccp_monte_carlo_function_eedf_restricted_lmea.yaml",
            "monte_carlo",
            mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION,
            "function_eedf_preintegrated_inelastic",
        ),
    ),
)
def test_production_function_eedf_maps_share_restricted_lmea_contract(
    name: str,
    source: str,
    estimator_schema: str | None,
    reaction_model: str,
) -> None:
    root = Path(__file__).resolve().parents[1]
    mapping = load_gec_ccp_mapping(root / "comsol_modes" / "maps" / name)

    assert mapping.bundle.expected_source == source
    assert mapping.bundle.expected_field_type == "dc"
    assert "flux" in mapping.bundle.expected_transport_definition
    assert mapping.closure.electron_transport == "swarm_mobility_einstein"
    assert mapping.closure.thermal_diffusion_model == ("off_restricted_diagonal")
    assert mapping.closure.reaction_model == reaction_model
    assert mapping.closure.elastic_energy_loss_model == ("external_solver_native")
    assert mapping.closure.lookup_jacobian == "exact"
    assert mapping.bundle.expected_mc_transport_estimator_schema == (estimator_schema)
    assert mapping.bundle.expected_mc_tail_estimator_schema == (
        mc_evidence.WEIGHTED_MC_TAIL_ESTIMATOR_SCHEMA_VERSION
        if source == "monte_carlo"
        else None
    )
    assert mapping.bundle.expected_mc_qualification_profile == (
        "function_eedf_restricted_lmea" if source == "monte_carlo" else None
    )
    assert mapping.results.role == "physical_target"
    assert mapping.run.include_builtin_reference is False
    assert mapping.model.baseline_output_mph is None
    assert mapping.run.period_dataset is None
    assert mapping.run.baseline_waveform_dataset is None


def test_propagator_map_uses_common_adaptive_function_spreadsheet() -> None:
    root = Path(__file__).resolve().parents[1]
    mapping = load_gec_ccp_mapping(
        root / "comsol_modes" / "maps" / "argon_gec_ccp_propagator_function_eedf.yaml"
    )

    assert mapping.bundle.expected_source == "propagator"
    assert mapping.closure.reaction_model == "function_eedf"
    assert mapping.closure.function_eedf is not None
    assert mapping.closure.function_eedf.table == "eedf_f0_comsol_2d.csv"
    assert mapping.closure.function_eedf.interpolation == (
        "structured_spreadsheet_linear_projection"
    )
    assert mapping.results.role == "physical_target"


def test_gec_ccp_owns_an_explicit_propagator_target_requirement() -> None:
    assert provenance_validation._GEC_CCP_PROPAGATOR_TARGET_REQUIREMENT == (
        PropagatorTargetRequirement(
            target="argon_gec_ccp_restricted_local_mean_energy",
            fields_Td=(3000.0, 3500.0, 4000.0),
            operating_range_bracket_Td=(3000.0, 3500.0),
            table_support_cap_Td=4000.0,
            mixture_id=0,
            medium_grid=(300, 48),
            fine_grid=(600, 72),
            max_scalar_refinement_limit=0.01,
            max_eedf_weighted_L1_limit=0.02,
        )
    )


def test_default_gec_plan_omits_builtin_reference_work(tmp_path: Path) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    obsolete_reference_fields = (
        "baseline_output_mph:",
        "include_builtin_reference:",
        "period_dataset:",
        "baseline_waveform_dataset:",
    )
    source = mapping_path.read_text(encoding="utf-8")
    mapping_path.write_text(
        "\n".join(
            line
            for line in source.splitlines()
            if not line.strip().startswith(obsolete_reference_fields)
        )
        + "\n",
        encoding="utf-8",
    )

    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    payload = json.loads(plan.plan_json.read_text(encoding="utf-8"))

    assert plan.baseline_run_java is None
    assert plan.baseline_export_java is None
    assert len(plan.expected_result_files) == 11
    assert {path.parent.name for path in plan.expected_result_files} == {"swarm_tables"}
    assert payload["solve"]["include_builtin_reference"] is False
    assert payload["solve"]["nonlinear_globalization"]["scope"] == ["external"]
    assert "baseline_run" not in payload["generated_java"]
    assert "baseline_export" not in payload["generated_java"]


def test_builtin_reference_fields_require_explicit_opt_in(
    tmp_path: Path,
) -> None:
    mapping_path, _ = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "  include_builtin_reference: true\n", ""
        ),
        encoding="utf-8",
    )

    with pytest.raises(
        GecCcpWorkflowError,
        match="baseline_output_mph is only valid",
    ):
        load_gec_ccp_mapping(mapping_path)


def test_axisymmetric_node_volume_weights_integrate_linear_2pi_r() -> None:
    radial = np.asarray([0.0, 1.0, 0.0, 1.0])
    axial = np.asarray([0.0, 0.0, 1.0, 1.0])

    weights, metadata = transport_support._axisymmetric_node_volume_weights(
        radial, axial
    )

    assert np.all(weights >= 0.0)
    assert float(np.sum(weights)) == pytest.approx(math.pi)
    assert metadata["node_count"] == 4
    assert metadata["triangle_count"] == 2


def test_axisymmetric_node_volume_weights_mask_gec_electrode_voids() -> None:
    coordinates = np.asarray(
        [
            [0.0, 0.0],
            [0.0, 0.0254],
            [0.0538, 0.0],
            [0.0538, 0.0254],
            [0.0538, -0.0381],
            [0.0538, 0.0635],
            [0.1016, -0.0381],
            [0.1016, 0.0635],
        ]
    )

    weights, metadata = transport_support._axisymmetric_node_volume_weights(
        coordinates[:, 0], coordinates[:, 1]
    )
    expected = math.pi * (0.0538**2 * 0.0254 + (0.1016**2 - 0.0538**2) * 0.1016)

    assert float(np.sum(weights)) == pytest.approx(expected)
    assert metadata["geometry"] == "stepped_gec_gas_domain"
    assert metadata["reconstructed_volume_m3"] == pytest.approx(
        metadata["inferred_exact_volume_m3"]
    )
    assert metadata["triangle_count"] < metadata["candidate_triangle_count"]


def test_mc_restricted_lmea_bundle_audit_retains_inactive_transport_failures() -> None:
    profile = "function_eedf_restricted_lmea"
    fields = [10.0, 20.0, 30.0, 40.0]
    inactive_reasons = [
        "reduced_diffusion_L_rse_unavailable_or_exceeds_threshold",
        "mc_weighted_growth_lag_not_converged:diffusion_L_m2_s",
    ]
    rows = [
        {
            "E_over_N_Td": field,
            "passed": 0,
            "qualification_profile": profile,
            "active_closure_quality_passed": 1,
            "active_closure_failure_reasons_json": "[]",
            "failure_reasons_json": json.dumps(inactive_reasons),
            "eedf_normalization_error": 1.0e-12,
            "valid_replicates": 4,
            "solver_transport_replicates": 4,
            "solver_diagnostics_available": 1,
            "mobility_rse": 0.02,
            "max_major_rate_rse": 0.03,
            "solver_mean_energy_stationarity_relative_ci95_bound": 0.02,
            "solver_mobility_stationarity_absolute_log_drift": 0.02,
            "solver_transport_mean_energy_max_relative_ci95_bound": 0.02,
            "solver_transport_mean_energy_relative_tolerance": 0.1,
            "solver_population_growth_relative_tolerance": 0.1,
            "solver_population_growth_gate_mode": "positive_paired_log",
            "solver_population_growth_max_relative_ci95_bound": 0.02,
            "quality_source": "monte_carlo_aggregate_quality",
        }
        for field in fields
    ]
    source_policy = {
        "qualification_profile": profile,
        "qualification_outputs": [
            "function_eedf",
            "reduced_mobility",
            "direct_elastic_energy_loss",
        ],
        "additional_mc_evidence": [
            "particle_diffusion_L_T",
            "energy_mobility",
            "restricted_energy_diffusion_L_T",
            "direct_reaction_rates",
        ],
        "transport_eligible_anchor_points": 4,
        "transport_eligible_E_over_N_Td": fields,
        "transport_excluded_E_over_N_Td": {},
        "full_transport_eligible_anchor_points": 0,
        "full_transport_eligible_E_over_N_Td": [],
        "full_transport_excluded_E_over_N_Td": {
            str(int(field)): inactive_reasons for field in fields
        },
        "active_mean_energy_anchors_eV": [1.0, 2.0, 3.0, 4.0],
        "maximum_adjacent_mean_energy_ratio": 2.0,
        "maximum_adjacent_mean_energy_gap_eV": 1.0,
    }
    thresholds = quality_thresholds_payload(
        QualityThresholds(
            mobility_rse=0.08,
            diffusion_rse=0.12,
            major_rate_rse=0.12,
        )
    )

    audit = mc_quality.independent_bundle_quality_audit(
        rows,
        source="monte_carlo",
        thresholds=thresholds,
        mc_sampling_plan=[{"e_over_n_Td": field} for field in fields],
        mc_estimator_schema=(
            mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
        ),
        source_policy=source_policy,
        transport_rows=[
            {"E_over_N_Td": field, "mean_energy_eV": field / 10.0} for field in fields
        ],
        valid_ranges={"E_over_N_Td": [10.0, 40.0]},
    )

    assert audit["passed"] is True
    assert all(item["failure_reasons"] == [] for item in audit["rows"])
    assert audit["transport_anchor_coverage"]["policy"].startswith(
        "qualified_active_closure"
    )


def test_weighted_growth_v6_quality_audit_accepts_positive_and_sparse_rows() -> None:
    thresholds = quality_thresholds_payload(QualityThresholds())
    row = {
        "E_over_N_Td": "1000",
        "passed": "1",
        "eedf_normalization_error": "0",
        "failure_reasons_json": "[]",
        "aggregate_failure_reasons_json": "[]",
        "aggregate_quality_passed": "1",
        "uncertainty_available": "1",
        "solver_diagnostics_available": "1",
        "solver_transport_qualified": "1",
        "valid_replicates": "4",
        "solver_transport_replicates": "4",
        "mobility_rse": "0.01",
        "energy_mobility_rse": "0.01",
        "diffusion_L_rse": "0.01",
        "diffusion_T_rse": "0.01",
        "energy_diffusion_L_rse": "0.01",
        "energy_diffusion_T_rse": "0.01",
        "max_major_rate_rse": "0.01",
        "solver_origin_stationarity_limiting_field": "mobility_m2_V_s",
        "solver_origin_stationarity_limiting_relative_ci95_bound": "0.05",
        "solver_origin_stationarity_limiting_relative_tolerance": "0.1",
        "solver_lag_convergence_max_relative_ci95_bound": "0.05",
        "solver_lag_convergence_relative_tolerance": "0.15",
        "solver_transport_mean_energy_max_relative_ci95_bound": "0.01",
        "solver_transport_mean_energy_relative_tolerance": "0.1",
        "solver_population_growth_gate_mode": "positive_paired_log",
        "solver_population_growth_max_relative_ci95_bound": "0.05",
        "solver_population_growth_poisson_interval_max_ratio": "1.2",
        "solver_population_growth_sparse_max_metric": "0.5",
        "solver_population_growth_relative_tolerance": "0.1",
        "quality_source": "monte_carlo_aggregate_quality",
    }
    sampling_plan = [
        {
            "e_over_n_Td": 1000.0,
            "particles": 512,
            "warmup_collisions": 8192,
            "max_collisions": 32768,
            "tail_max_collisions": 16384,
            "replicas": 4,
            "transport_correlation_lag_barriers": 64,
            "transport_estimator": "single_field",
        }
    ]

    audit = mc_quality.independent_bundle_quality_audit(
        [row],
        source="monte_carlo",
        thresholds=thresholds,
        mc_sampling_plan=sampling_plan,
        mc_estimator_schema=(
            mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
        ),
    )

    assert audit["passed"] is True

    legacy_weighted_lag_tolerance = dict(row)
    legacy_weighted_lag_tolerance["solver_lag_convergence_relative_tolerance"] = "0.25"
    audit = mc_quality.independent_bundle_quality_audit(
        [legacy_weighted_lag_tolerance],
        source="monte_carlo",
        thresholds=thresholds,
        mc_sampling_plan=sampling_plan,
        mc_estimator_schema=(
            mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
        ),
    )
    assert audit["passed"] is False
    assert (
        "solver_lag_convergence_max_relative_ci95_bound"
        in (audit["rows"][0]["failure_reasons"])
    )

    sparse = dict(row)
    sparse["solver_population_growth_gate_mode"] = "sparse_negligible_transport"
    sparse["solver_population_growth_max_relative_ci95_bound"] = "0.4"
    sparse["solver_population_growth_poisson_interval_max_ratio"] = "0.8"
    sparse["solver_population_growth_sparse_max_metric"] = "0.02"
    audit = mc_quality.independent_bundle_quality_audit(
        [sparse],
        source="monte_carlo",
        thresholds=thresholds,
        mc_sampling_plan=sampling_plan,
        mc_estimator_schema=(
            mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
        ),
    )
    assert audit["passed"] is True

    diffusion_limited = dict(row)
    diffusion_limited["solver_origin_stationarity_limiting_field"] = "diffusion_L_m2_s"
    diffusion_limited["solver_origin_stationarity_limiting_relative_ci95_bound"] = (
        "0.115"
    )
    diffusion_limited["solver_origin_stationarity_limiting_relative_tolerance"] = "0.25"
    audit = mc_quality.independent_bundle_quality_audit(
        [diffusion_limited],
        source="monte_carlo",
        thresholds=thresholds,
        mc_sampling_plan=sampling_plan,
        mc_estimator_schema=(
            mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
        ),
    )
    assert audit["passed"] is True

    invalid_field = dict(diffusion_limited)
    invalid_field["solver_origin_stationarity_limiting_field"] = "unknown"
    audit = mc_quality.independent_bundle_quality_audit(
        [invalid_field],
        source="monte_carlo",
        thresholds=thresholds,
        mc_sampling_plan=sampling_plan,
        mc_estimator_schema=(
            mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
        ),
    )
    assert audit["passed"] is False
    assert (
        "solver_origin_stationarity_limiting_field"
        in (audit["rows"][0]["failure_reasons"])
    )

    mismatched_tolerance = dict(row)
    mismatched_tolerance["solver_origin_stationarity_limiting_relative_tolerance"] = (
        "0.25"
    )
    audit = mc_quality.independent_bundle_quality_audit(
        [mismatched_tolerance],
        source="monte_carlo",
        thresholds=thresholds,
        mc_sampling_plan=sampling_plan,
        mc_estimator_schema=(
            mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
        ),
    )
    assert audit["passed"] is False
    assert (
        "solver_origin_stationarity_limiting_relative_ci95_bound"
        in (audit["rows"][0]["failure_reasons"])
    )

    for bound_name, rejected_value in (
        ("solver_lag_convergence_max_relative_ci95_bound", "0.3"),
        ("solver_population_growth_max_relative_ci95_bound", "0.2"),
    ):
        rejected = dict(row)
        rejected[bound_name] = rejected_value
        audit = mc_quality.independent_bundle_quality_audit(
            [rejected],
            source="monte_carlo",
            thresholds=thresholds,
            mc_sampling_plan=sampling_plan,
            mc_estimator_schema=(
                mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
            ),
        )
        assert audit["passed"] is False
        assert bound_name in audit["rows"][0]["failure_reasons"]

    for required_name in (
        "solver_lag_convergence_max_relative_ci95_bound",
        "solver_lag_convergence_relative_tolerance",
        "solver_population_growth_max_relative_ci95_bound",
        "solver_population_growth_relative_tolerance",
    ):
        missing = dict(row)
        missing[required_name] = ""
        with pytest.raises(GecCcpWorkflowError, match=required_name):
            mc_quality.independent_bundle_quality_audit(
                [missing],
                source="monte_carlo",
                thresholds=thresholds,
                mc_sampling_plan=sampling_plan,
                mc_estimator_schema=(
                    mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
                ),
            )

    inflated = dict(row)
    inflated["solver_population_growth_max_relative_ci95_bound"] = "0.2"
    inflated["solver_population_growth_relative_tolerance"] = "1.0"
    audit = mc_quality.independent_bundle_quality_audit(
        [inflated],
        source="monte_carlo",
        thresholds=thresholds,
        mc_sampling_plan=sampling_plan,
        mc_estimator_schema=(
            mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
        ),
    )
    assert audit["passed"] is False
    assert (
        "solver_population_growth_max_relative_ci95_bound"
        in audit["rows"][0]["failure_reasons"]
    )

    planned_fields = (500.0, 1000.0, 1500.0, 2000.0)
    planned = [
        {
            "e_over_n_Td": field,
            "particles": 512,
            "warmup_collisions": 8192,
            "max_collisions": 32768,
            "tail_max_collisions": 16384,
            "replicas": 4,
            "transport_correlation_lag_barriers": 64,
        }
        for field in planned_fields
    ]

    def quality_row(field: float) -> dict[str, str]:
        result = dict(row)
        result["E_over_N_Td"] = str(field)
        return result

    endpoint_trimmed = mc_quality.independent_bundle_quality_audit(
        [quality_row(1000.0), quality_row(1500.0)],
        source="monte_carlo",
        thresholds=thresholds,
        mc_sampling_plan=planned,
        mc_estimator_schema=(
            mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
        ),
    )
    coverage = endpoint_trimmed["transport_anchor_coverage"]
    assert coverage["qualified_E_over_N_Td"] == [1000.0, 1500.0]
    assert coverage["failed_outside_exported_support_E_over_N_Td"] == [
        500.0,
        2000.0,
    ]

    with pytest.raises(
        GecCcpWorkflowError,
        match="crosses failed internal planned E/N anchors: 1000",
    ):
        mc_quality.independent_bundle_quality_audit(
            [
                quality_row(500.0),
                quality_row(1500.0),
                quality_row(2000.0),
            ],
            source="monte_carlo",
            thresholds=thresholds,
            mc_sampling_plan=planned,
            mc_estimator_schema=(
                mc_evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
            ),
        )


def test_external_elastic_loss_supports_function_eedf_mobility_target(
    tmp_path: Path,
) -> None:
    mapping_path, _ = _write_fixture(tmp_path)
    assert (
        load_gec_ccp_mapping(mapping_path).closure.elastic_energy_loss_model
        == "comsol_mratio"
    )

    target = (
        mapping_path.read_text(encoding="utf-8")
        .replace(
            "  electron_transport: comsol_specify_all_restricted",
            "  electron_transport: swarm_mobility_einstein",
        )
        .replace(
            "  reaction_model: external_rates",
            "  reaction_model: function_eedf\n"
            "  elastic_energy_loss_model: external_solver_native\n"
            "  function_eedf:\n"
            "    table: eedf_f0_comsol_2d.csv\n"
            "    function_tag: sw_eedf_test\n"
            "    interpolation: structured_spreadsheet_linear_projection\n"
            "    extrapolation: constant",
        )
    )
    mapping_path.write_text(target, encoding="utf-8")
    closure = load_gec_ccp_mapping(mapping_path).closure
    assert closure.reaction_model == "function_eedf"
    assert closure.electron_transport == "swarm_mobility_einstein"
    assert closure.elastic_energy_loss_model == "external_solver_native"

    incompatible = target.replace(
        "  electron_transport: swarm_mobility_einstein",
        "  electron_transport: comsol",
    )
    mapping_path.write_text(incompatible, encoding="utf-8")
    with pytest.raises(
        GecCcpWorkflowError,
        match="external_solver_native elastic energy loss is incompatible",
    ):
        load_gec_ccp_mapping(mapping_path)


def test_gec_ccp_requires_explicit_result_role(tmp_path: Path) -> None:
    mapping_path, _ = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "  role: physical_target\n",
            "",
        ),
        encoding="utf-8",
    )

    with pytest.raises(
        GecCcpWorkflowError,
        match="role must be a non-empty string",
    ):
        load_gec_ccp_mapping(mapping_path)


def test_lookup_rejects_duplicate_anchor_instead_of_overwriting() -> None:
    rows = [
        {"mean_energy_eV": "1", "coefficient": "2"},
        {"mean_energy_eV": "1.0", "coefficient": "3"},
    ]

    with pytest.raises(
        GecCcpWorkflowError,
        match="contains duplicate mean_energy_eV",
    ):
        data_helpers._lookup(rows, "mean_energy_eV", "coefficient")


def test_two_term_hybrid_binds_einstein_particle_and_external_energy_transport(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8")
        .replace(
            "electron_transport: comsol_specify_all_restricted",
            "electron_transport: swarm_hybrid_einstein_de",
        )
        .replace("  zero_field_isotropization_Td: 0.1\n", ""),
        encoding="utf-8",
    )

    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    source = generate_apply_java(plan.mapping)

    assert '"SpecifyElectronDensityAndEnergy", "SpecifyAll"' in source
    assert "sw_log_muN_e" in source
    assert "sw_log_muenN_e" in source
    assert "sw_log_DenN_e" in source
    assert "sw_log_DeN_L_e" not in source
    assert "sw_log_DeN_T_e" not in source
    assert "ptp.Te" in source
    plan_data = json.loads(plan.plan_json.read_text(encoding="utf-8"))
    closure = plan_data["closure"]
    assert closure["electron_transport"] == "swarm_hybrid_einstein_de"
    assert closure["transport_inputs"]["DeN"] == {
        "source": "comsol_einstein_from_external_mobility",
        "external_table_column": None,
        "formula": "DeN=muN*ptp.Te",
    }
    assert closure["transport_inputs"]["muenN"]["source"] == ("external_swarm_table")
    assert closure["wall_closure"] == {
        "owner": "input_COMSOL_GEC_model",
        "model": "WallDriftDiffusion_Te_thermal_velocity",
        "external_eedf_half_range_moments_consumed": False,
        "angular_half_range_response_identified": False,
        "scope": "shared_boundary_model_not_external_swarm_closure",
    }
    assert closure["transport_inputs"]["DenN"]["source"] == ("external_swarm_table")
    assert closure["transport_tensor"]["external_particle_diffusion_consumed"] is False
    assert closure["data_processing"]["swarm_source_postprocess"] == "none"


def test_external_two_term_elastic_loss_generation_export_and_support(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    _enable_external_elastic_loss(mapping_path, bundle)

    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    java = plan.apply_java.read_text(encoding="utf-8")
    external_export = plan.external_export_java.read_text(encoding="utf-8")
    baseline_export = plan.baseline_export_java.read_text(encoding="utf-8")
    payload = json.loads(plan.plan_json.read_text(encoding="utf-8"))

    assert java.count('feature("eir1").active(false)') == 1
    assert java.count('create("swElLoss", "GeneralPowerDeposition", 2)') == 1
    assert 'feature("swElLoss").selection().set(new int[]{1})' in java
    assert "-ptp.ne*ptp.n_wAr*exp(sw_logKel(sw_logeps))" in java
    assert "1[eV*m^3/s]" in java
    assert '.feature("eir1").set("kf"' not in java
    assert java.count('"SpecifyReactionUsing", "RateConstant"') == 2
    assert 'model.func("sw_logKel").set("interp", "piecewisecubic")' in java
    assert "ptp.kf_1" not in external_export
    assert "ptp.kf_2" in external_export
    assert "ptp.kf_3" in external_export
    assert "ptp.n_wAr" in external_export
    assert "ptp.Qgen" not in external_export
    assert "-ptp.ne*ptp.n_wAr*exp(sw_logKel(sw_logeps))" in external_export
    assert "ptp.xdintop_ptp3" in external_export
    assert "ptp.kf_1" in baseline_export
    assert "0[W/m^3]" in baseline_export
    support = payload["closure"]["active_closure_support"]
    assert support["supports_mean_energy_eV"]["elastic_energy_loss"] == [
        1.0,
        2.0,
    ]
    assert (
        payload["closure"]["elastic_energy_loss"]["energy_equation_owner"]
        == "GeneralPowerDeposition:swElLoss"
    )


def test_monte_carlo_hybrid_uses_energy_LT_but_einstein_particle_diffusion(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    _enable_external_elastic_loss(mapping_path, bundle, source="monte_carlo")

    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    java = plan.apply_java.read_text(encoding="utf-8")
    payload = json.loads(plan.plan_json.read_text(encoding="utf-8"))
    transport = payload["closure"]["transport_inputs"]

    assert "sw_log_DenN_L_e" in java
    assert "sw_log_DenN_T_e" in java
    assert "sw_log_DeN_L_e" not in java
    assert "sw_log_DeN_T_e" not in java
    assert "sw_logeps_el" in java
    assert transport["DeN"]["formula"] == "DeN=muN*ptp.Te"
    assert transport["DenN"]["external_table_column"] == [
        "reduced_electron_energy_diffusion_L_m2_s_m3",
        "reduced_electron_energy_diffusion_T_m2_s_m3",
    ]
    tensor = payload["closure"]["transport_tensor"]
    assert tensor["external_particle_diffusion_consumed"] is False
    assert tensor["full_gradient_response_identified"] is False
    assert tensor["zero_field_isotropization_Td"] == 0.1
    statistics = payload["closure"]["elastic_energy_loss"]["statistical_evidence"]
    assert statistics["minimum_valid_replicates"] == 3
    assert statistics["maximum_relative_standard_error"] == 0.05
    assert statistics["relative_standard_error_qualified"] is True


def test_monte_carlo_elastic_loss_rejects_unqualified_RSE(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    _enable_external_elastic_loss(mapping_path, bundle, source="monte_carlo")
    table = bundle / "elastic_energy_loss_vs_mean_energy.csv"
    rows = list(csv.DictReader(table.open(encoding="utf-8")))
    for row in rows:
        coefficient = float(row["elastic_energy_loss_rate_coefficient_eV_m3_s"])
        standard_error = 0.5 * coefficient
        critical = float(row["ci95_critical_value"])
        row["elastic_energy_loss_standard_error_eV_m3_s"] = standard_error
        row["elastic_energy_loss_relative_standard_error"] = 0.5
        row["elastic_energy_loss_ci95_low_eV_m3_s"] = (
            coefficient - critical * standard_error
        )
        row["elastic_energy_loss_ci95_high_eV_m3_s"] = (
            coefficient + critical * standard_error
        )
    with table.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["tables"][table.name]["sha256"] = sha256(table.read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(GecCcpWorkflowError, match="CI/RSE contract"):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_gec_ccp_rejects_obsolete_full_transport_id_with_migration(
    tmp_path: Path,
) -> None:
    mapping_path, _ = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "electron_transport: comsol_specify_all_restricted",
            "electron_transport: swarm_field_aligned_full",
        ),
        encoding="utf-8",
    )

    with pytest.raises(
        GecCcpWorkflowError,
        match="swarm_field_aligned_full is obsolete; migrate to "
        "comsol_specify_all_restricted",
    ):
        load_gec_ccp_mapping(mapping_path)


def test_gec_ccp_rejects_obsolete_expected_postprocess_key(
    tmp_path: Path,
) -> None:
    mapping_path, _ = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "  expected_transport_definition: flux\n",
            "  expected_transport_definition: flux\n"
            "  expected_postprocess: regularized\n",
        ),
        encoding="utf-8",
    )

    with pytest.raises(
        GecCcpWorkflowError,
        match=r"unknown bundle key\(s\): expected_postprocess",
    ):
        load_gec_ccp_mapping(mapping_path)


def test_gec_ccp_requires_explicit_gradient_response_policy(
    tmp_path: Path,
) -> None:
    mapping_path, _ = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "  gradient_response_policy: standard_local_energy\n", ""
        ),
        encoding="utf-8",
    )

    with pytest.raises(
        GecCcpWorkflowError,
        match="gradient_response_policy must be a non-empty string",
    ):
        load_gec_ccp_mapping(mapping_path)


@pytest.mark.parametrize(
    "electron_transport",
    (
        "comsol",
        "swarm_mobility_einstein",
        "swarm_hybrid_einstein_de",
        "comsol_specify_all_restricted",
    ),
)
def test_gec_ccp_require_full_fails_closed_for_every_transport_mode(
    tmp_path: Path,
    electron_transport: str,
) -> None:
    mapping_path, _ = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8")
        .replace(
            "electron_transport: comsol_specify_all_restricted",
            f"electron_transport: {electron_transport}",
        )
        .replace(
            "gradient_response_policy: standard_local_energy",
            "gradient_response_policy: require_full",
        ),
        encoding="utf-8",
    )

    with pytest.raises(
        GecCcpWorkflowError,
        match="require_full is unavailable",
    ):
        load_gec_ccp_mapping(mapping_path)


def test_two_term_rejects_mc_only_zero_field_tensor_scale(
    tmp_path: Path,
) -> None:
    mapping_path, _ = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "electron_transport: comsol_specify_all_restricted",
            "electron_transport: comsol_specify_all_restricted\n"
            "  zero_field_isotropization_Td: 0.1",
        ),
        encoding="utf-8",
    )

    with pytest.raises(
        GecCcpWorkflowError,
        match="only valid for a monte_carlo",
    ):
        load_gec_ccp_mapping(mapping_path)


def test_monte_carlo_physical_target_requires_estimator_schema(
    tmp_path: Path,
) -> None:
    mapping_path, _ = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8")
        .replace("expected_source: two_term", "expected_source: monte_carlo")
        .replace(
            "electron_transport: comsol_specify_all_restricted",
            "electron_transport: comsol_specify_all_restricted\n"
            "  zero_field_isotropization_Td: 0.1",
        ),
        encoding="utf-8",
    )

    with pytest.raises(
        GecCcpWorkflowError,
        match="must declare bundle.expected_mc_transport_estimator_schema",
    ):
        load_gec_ccp_mapping(mapping_path)


def test_two_term_required_transport_kernel_rejects_legacy_bundle(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "  expected_transport_definition: flux",
            "  expected_transport_definition: flux\n"
            "  expected_two_term_transport_kernel_schema: "
            "two_term_temporal_growth_transport.v1",
        ),
        encoding="utf-8",
    )

    with pytest.raises(
        GecCcpWorkflowError,
        match="lacks the required temporal-growth",
    ):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_two_term_pt_transport_kernel_is_reintegrated_with_shared_model(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    _enable_two_term_temporal_growth_fixture(
        mapping_path,
        bundle,
        growth_frequency_s_inv=0.0,
    )

    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    payload = json.loads(plan.plan_json.read_text(encoding="utf-8"))
    bundle_kernel = payload["bundle"]["two_term_transport_kernel"]
    assert bundle_kernel["schema"] == ("two_term_temporal_growth_transport.v1")
    joint = payload["closure"]["two_term_joint_consistency"]
    assert joint["passed"] is True
    assert joint["temporal_growth_transport"] == bundle_kernel
    assert joint["transport_kernel_contract"]["effective_momentum_model"] == (
        "nu_m_eff(epsilon)=nu_m(epsilon)+growth_frequency"
    )
    assert joint["transport_kernel_contract"]["aggregate_floor_m2"] == (
        pytest.approx(1.0e-24)
    )


def test_two_term_pt_transport_audit_rejects_unshifted_coefficients(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    _enable_two_term_temporal_growth_fixture(
        mapping_path,
        bundle,
        growth_frequency_s_inv=1.0e9,
    )

    with pytest.raises(
        GecCcpWorkflowError,
        match="transport joint consistency audit failed",
    ):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_monte_carlo_field_aligned_transport_requires_zero_field_scale(
    tmp_path: Path,
) -> None:
    mapping_path, _ = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "expected_source: two_term",
            "expected_source: monte_carlo\n"
            "  expected_mc_transport_estimator_schema: direct_mc_transport.v4",
        ),
        encoding="utf-8",
    )

    with pytest.raises(
        GecCcpWorkflowError,
        match="must be a positive number",
    ):
        load_gec_ccp_mapping(mapping_path)


def test_gec_ccp_rejects_postprocessed_solver_bundle(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["source_policy"]["postprocess"] = "regularized"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(
        GecCcpWorkflowError,
        match=r"requires source_policy\.postprocess=none",
    ):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_gec_ccp_rejects_monte_carlo_bundle_without_raw_rate_evidence(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    _enable_monte_carlo_fixture(mapping_path, bundle)
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    del manifest["tables"]["rate_evidence.csv"]
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(
        GecCcpWorkflowError,
        match="pure Monte Carlo GEC input requires canonical raw trajectory",
    ):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_gec_ccp_rejects_non_power_of_two_mc_correlation_lag(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    _enable_monte_carlo_fixture(mapping_path, bundle)
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    sampling_plan = json.loads(manifest["hashes"]["mc_sampling_plan_json"])
    sampling_plan[0]["transport_correlation_lag_barriers"] = 6
    sampling_json = json.dumps(
        sampling_plan,
        sort_keys=True,
        separators=(",", ":"),
    )
    manifest["hashes"]["mc_sampling_plan_json"] = sampling_json
    manifest["hashes"]["mc_sampling_plan_sha256"] = sha256(
        sampling_json.encode("utf-8")
    ).hexdigest()
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(
        GecCcpWorkflowError,
        match="invalid mc_sampling_plan_json",
    ):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_two_term_grad_diffusivity_generates_and_records_restricted_matrix(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "thermal_diffusion_model: off_restricted_diagonal",
            "thermal_diffusion_model: comsol_grad_diffusivity",
        ),
        encoding="utf-8",
    )

    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    java = generate_apply_java(plan.mapping)
    manifest = json.loads(plan.plan_json.read_text(encoding="utf-8"))
    response = manifest["closure"]["gradient_response"]

    assert '.prop("ElectronProperties").set("IncludeThermalDiffusion", true)' in java
    assert response["response_matrix_assumptions"] == {
        "A_nn": "De",
        "A_nE": "dDe/dln(mean_energy)",
        "A_un": "Den",
        "A_uE": "Den",
    }
    assert response["full_gradient_response_identified"] is False
    assert response["full_physical_eligible"] is False
    assert response["scalar_diffusion_audit"]["passed"] is True


def test_two_term_grad_diffusivity_rejects_anisotropic_diffusion(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "thermal_diffusion_model: off_restricted_diagonal",
            "thermal_diffusion_model: comsol_grad_diffusivity",
        ),
        encoding="utf-8",
    )
    table = bundle / "transport_vs_mean_energy.csv"
    rows = list(csv.DictReader(table.open(encoding="utf-8")))
    rows[0]["reduced_diffusion_T_m2_s_m3"] = str(
        1.01 * float(rows[0]["reduced_diffusion_L_m2_s_m3"])
    )
    _write_csv(table, list(rows[0]), [[row[name] for name in rows[0]] for row in rows])
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["tables"]["transport_vs_mean_energy.csv"]["sha256"] = sha256(
        table.read_bytes()
    ).hexdigest()
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(
        GecCcpWorkflowError,
        match="requires longitudinal and transverse particle/energy diffusion",
    ):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_two_term_scalar_binding_rejects_anisotropy_with_thermal_off(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    table = bundle / "transport_vs_mean_energy.csv"
    rows = list(csv.DictReader(table.open(encoding="utf-8")))
    rows[0]["reduced_electron_energy_diffusion_T_m2_s_m3"] = str(
        1.01 * float(rows[0]["reduced_electron_energy_diffusion_L_m2_s_m3"])
    )
    _write_csv(
        table,
        list(rows[0]),
        [[row[name] for name in rows[0]] for row in rows],
    )
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["tables"]["transport_vs_mean_energy.csv"]["sha256"] = sha256(
        table.read_bytes()
    ).hexdigest()
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(
        GecCcpWorkflowError,
        match="requires longitudinal and transverse particle/energy diffusion",
    ):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_external_rates_bind_positive_suffix_without_artificial_floor(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    rates_path = bundle / "rates_vs_mean_energy.csv"
    rows = list(csv.DictReader(rates_path.open(encoding="utf-8")))
    rows.extend(
        [
            {
                "mean_energy_eV": "0.5",
                "process_type": "elastic",
                "rate_coefficient_m3_s": "5e-15",
            },
            {
                "mean_energy_eV": "0.5",
                "process_type": "excitation",
                "rate_coefficient_m3_s": "0",
            },
            {
                "mean_energy_eV": "0.5",
                "process_type": "ionization",
                "rate_coefficient_m3_s": "0",
            },
        ]
    )
    _write_csv(
        rates_path,
        list(rows[0]),
        [[row[name] for name in rows[0]] for row in rows],
    )
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["tables"]["rates_vs_mean_energy.csv"]["sha256"] = sha256(
        rates_path.read_bytes()
    ).hexdigest()
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    mapping = load_gec_ccp_mapping(mapping_path, bundle_path=bundle)
    support = closure_contract._active_rate_support_metadata(mapping, rows)
    assert support["artificial_floor"] is False
    assert support["processes"]["excitation"]["zero_rows_observed"] == 1
    assert support["processes"]["ionization"]["zero_rows_observed"] == 1
    assert support["processes"]["excitation"]["positive_support_mean_energy_eV"] == [
        1.0,
        2.0,
    ]
    assert joint_consistency._active_closure_mean_energy_support(
        mapping, preintegrated_rate_rows=rows
    )["common_intersection_mean_energy_eV"] == [1.0, 2.0]
    plan_data = json.loads(plan.plan_json.read_text(encoding="utf-8"))
    evidence = plan_data["closure"]["upstream_eedf_evidence"]
    assert evidence["active_in_comsol"] is False
    assert evidence["binding_status"] == "inactive_evidence_only"
    assert evidence["role"] == "upstream_projected_eedf_evidence"
    assert evidence["bundle_artifact_role"] == ("canonical_comsol_function_eedf_input")
    assert evidence["bundle_comsol_import_capable"] is True
    assert "canonical_comsol_input" not in evidence
    assert (
        evidence["table_sha256"]
        == sha256((bundle / "eedf_f0_comsol_2d.csv").read_bytes()).hexdigest()
    )
    joint = plan_data["closure"]["two_term_joint_consistency"]
    assert "active_function_eedf" not in joint
    assert (
        joint["upstream_projected_eedf_evidence"]["sha256"]
        == (evidence["table_sha256"])
    )
    assert joint["cross_section_identity"]["passed"] is True
    rate_consistency = joint["external_rate_consistency"]
    assert rate_consistency["passed"] is True
    assert rate_consistency["eedf_active_in_comsol"] is False
    assert rate_consistency["eedf_role"] == ("inactive_upstream_consistency_evidence")
    assert rate_consistency["active_inelastic_rate_owner"] == (
        "two_term_integrated_rates:rates_vs_mean_energy.csv"
    )
    assert rate_consistency["active_rate_values_modified"] is False
    assert Path(rate_consistency["audit_csv"]).name == ("upstream_eedf_rate_audit.csv")
    assert set(rate_consistency["processes"]) == {
        "excitation",
        "ionization",
    }
    assert all(
        process["p95_relative_error_limit"] == pytest.approx(0.05)
        and process["maximum_relative_error_limit"] == pytest.approx(0.10)
        and process["normalized_rmse_limit"] == pytest.approx(0.02)
        for process in rate_consistency["processes"].values()
    )
    assert joint["closure_ownership"] == {
        "projected_eedf": "inactive_upstream_consistency_evidence",
        "inelastic_rate_owner": ("two_term_integrated_rates:rates_vs_mean_energy.csv"),
        "reintegration_modifies_active_rates": False,
    }
    assert joint["failure_reasons"] == []
    upstream_rates = plan_data["closure"]["upstream_eedf_rate_consistency"]
    assert upstream_rates["passed"] is True
    assert upstream_rates["active_inelastic_rate_owner"] == (
        "two_term_integrated_rates:rates_vs_mean_energy.csv"
    )


def test_monte_carlo_uses_standard_restricted_energy_density_closure(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    _enable_monte_carlo_fixture(mapping_path, bundle)

    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    java = generate_apply_java(plan.mapping)
    response = json.loads(plan.plan_json.read_text(encoding="utf-8"))["closure"][
        "gradient_response"
    ]
    assert response["comsol_executable_eligible"] is True
    assert response["restricted_closure_eligible"] is True
    assert response["full_physical_eligible"] is False
    assert response["independent_2x2_gradient_response_identified"] is False
    assert response["closure_scope"] == (
        "standard_COMSOL_electron_density_and_energy_density_flux"
    )
    assert response["response_matrix_provenance"]["A_un"] == (
        "Monte_Carlo_density_packet_energy_diffusion_for_COMSOL_energy_density_flux"
    )
    assert response["response_matrix_provenance"]["A_uE"] == (
        "same_MC_Den_applied_to_grad_log_energy_density_by_standard_COMSOL_closure"
    )
    assert "sw_logeps_transport" in java
    assert "sw_logeps_rate_ionization" in java
    assert "sw_knorm_ionization_e" in java
    assert "sw_logk_ionization_e" not in java
    closure = json.loads(plan.plan_json.read_text(encoding="utf-8"))["closure"]
    upstream_rates = closure["upstream_eedf_rate_consistency"]
    assert upstream_rates["passed"] is True
    assert upstream_rates["source"] == "monte_carlo"
    assert upstream_rates["projected_eedf"]["active_in_comsol"] is False
    assert upstream_rates["projected_eedf"]["role"] == (
        "inactive_upstream_consistency_evidence"
    )
    assert upstream_rates["active_inelastic_rate_owner"] == (
        "monte_carlo_direct_trajectory_rates:rates_vs_mean_energy.csv"
    )
    assert upstream_rates["reintegration_modifies_active_rates"] is False
    assert upstream_rates["cross_section_identity"]["passed"] is True
    assert upstream_rates["rate_consistency"]["artificial_rate_floor"] is False
    assert (
        Path(upstream_rates["rate_consistency"]["audit_csv"]).name
        == "upstream_eedf_rate_audit.csv"
    )
    assert upstream_rates["magnitude_policy"] == {
        "significance_fraction_of_process_peak": 1.0e-5,
        "below_significance_floor": ("record_only_for_pointwise_relative_error"),
        "artificial_rate_floor": False,
        "censored_zero_handling": (
            "preserve_zero_and_separately_retain_MC_upper_bound_evidence"
        ),
    }
    assert closure["mc_rate_interpolation"]["passed"] is True
    assert closure["mc_rate_censoring"]["physics_qualification"]["passed"] is True
    assert response["source_interpretation"] == (
        "Monte_Carlo_direct_particle_and_energy_density_packet_"
        "coefficients_under_standard_COMSOL_restricted_closure;_"
        "not_an_independent_2x2_gradient_response"
    )
    handling = json.loads(plan.plan_json.read_text(encoding="utf-8"))["closure"][
        "reaction_handling"
    ]
    assert all(
        item["rate_source"] == "monte_carlo_trajectory_time_average_sigma_v"
        for item in handling
    )
    assert all(
        item["eedf_source"]
        == "independent_monte_carlo_eedf_consistency_evidence_not_rate_owner"
        for item in handling
    )


@pytest.mark.parametrize(
    ("corruption", "message"),
    (
        ("source_policy", "eligible anchors disagree with quality.csv"),
        ("transport", "transport table anchors disagree with quality.csv"),
        ("valid_range", "valid E/N range disagrees"),
    ),
)
def test_monte_carlo_transport_anchor_contract_crosschecks_bundle_artifacts(
    tmp_path: Path,
    corruption: str,
    message: str,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    _enable_monte_carlo_fixture(mapping_path, bundle)
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if corruption == "source_policy":
        manifest["source_policy"]["transport_eligible_E_over_N_Td"] = [10.0]
    elif corruption == "valid_range":
        manifest["valid_ranges"]["E_over_N_Td"] = [10.0, 30.0]
    else:
        transport_path = bundle / "transport_vs_mean_energy.csv"
        rows = list(csv.DictReader(transport_path.open(encoding="utf-8")))
        rows[0]["E_over_N_Td"] = "15"
        _write_csv(
            transport_path,
            list(rows[0]),
            [[row[name] for name in rows[0]] for row in rows],
        )
        manifest["tables"]["transport_vs_mean_energy.csv"]["sha256"] = sha256(
            transport_path.read_bytes()
        ).hexdigest()
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(GecCcpWorkflowError, match=message):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_monte_carlo_binds_censored_zero_rates_without_floor(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    _enable_monte_carlo_fixture(mapping_path, bundle)
    rates_path = bundle / "rates_vs_mean_energy.csv"
    rates = list(csv.DictReader(rates_path.open(encoding="utf-8")))
    for row in rates:
        if float(row["mean_energy_eV"]) == 1.0 and row["process_type"] == "ionization":
            row["rate_coefficient_m3_s"] = "0"
    ionization_peak = max(
        float(row["rate_coefficient_m3_s"])
        for row in rates
        if row["process_type"] == "ionization"
    )
    relevance_fraction = 1.0e-5
    upper = 0.5 * relevance_fraction * ionization_peak
    exposure = -math.log(0.05) / upper
    _write_csv(
        rates_path,
        list(rates[0]),
        [[row[name] for name in rates[0]] for row in rates],
    )
    evidence_path = bundle / "rate_evidence.csv"
    evidence = list(csv.DictReader(evidence_path.open(encoding="utf-8")))
    for row in evidence:
        if float(row["mean_energy_eV"]) == 1.0 and row["process_type"] == "ionization":
            row.update(
                {
                    "rate_coefficient_mean_m3_s": "0",
                    "rate_coefficient_ci95_low_m3_s": "",
                    "rate_coefficient_ci95_high_m3_s": "",
                    "estimate_status": "censored_all_zero",
                    "uncertainty_available": "False",
                    "pooled_event_count": "0",
                    "pooled_target_exposure_s_m3": str(exposure),
                    "pooled_all_zero_upper_95_m3_s": str(upper),
                    "pooled_zero_event_status": "all_zero_upper_95",
                }
            )
    _write_csv(
        evidence_path,
        list(evidence[0]),
        [[row[name] for name in evidence[0]] for row in evidence],
    )
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["tables"]["rates_vs_mean_energy.csv"]["sha256"] = sha256(
        rates_path.read_bytes()
    ).hexdigest()
    manifest["tables"]["rate_evidence.csv"]["sha256"] = sha256(
        evidence_path.read_bytes()
    ).hexdigest()
    manifest["quality_thresholds"]["required_rate_min_process_peak_fraction"] = (
        relevance_fraction
    )
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    closure = json.loads(plan.plan_json.read_text(encoding="utf-8"))["closure"]
    java = generate_apply_java(plan.mapping)

    assert closure["active_closure_support"]["common_intersection_mean_energy_eV"] == [
        1.0,
        2.0,
    ]
    assert (
        closure["active_rate_support"]["processes"]["ionization"][
            "rows_omitted_before_positive_suffix"
        ]
        == 0
    )
    assert (
        closure["active_rate_support"]["nonnegative_zero_preserving_active_support"]
        is True
    )
    interpolation = closure["mc_rate_interpolation"]
    assert interpolation["passed"] is True
    assert interpolation["processes"]["ionization"]["minimum_normalized_rate"] >= 0.0
    assert (
        interpolation["processes"]["ionization"]["maximum_zero_anchor_absolute_error"]
        == 0.0
    )
    censoring = closure["mc_rate_censoring"]
    assert censoring["nominal_numerical_binding"]["passed"] is True
    assert censoring["physics_qualification"]["passed"] is True
    assert censoring["censored_anchor_count"] == 1
    censored = censoring["physics_qualification"]["censored_anchors"][0]
    assert censored["passed"] is True
    assert censored["pooled_all_zero_upper_95_m3_s"] == pytest.approx(upper)
    assert censored["required_upper_limit_m3_s"] == pytest.approx(
        relevance_fraction * ionization_peak
    )
    upstream = closure["upstream_eedf_rate_consistency"]
    assert upstream["passed"] is True
    ionization = upstream["rate_consistency"]["processes"]["ionization"]
    assert ionization["zero_rate_anchors"] == 1
    assert ionization["relevant_zero_rate_anchors"] == 0
    assert upstream["rate_consistency"]["artificial_rate_floor"] is False
    with Path(upstream["rate_consistency"]["audit_csv"]).open(
        encoding="utf-8", newline=""
    ) as handle:
        audit_rows = list(csv.DictReader(handle))
    censored_row = next(
        row
        for row in audit_rows
        if row["process_type"] == "ionization" and float(row["mean_energy_eV"]) == 1.0
    )
    assert censored_row["solver_rate_is_zero"] == "True"
    assert censored_row["relevant"] == "False"
    assert censored_row["relative_error_if_relevant"] == ""
    assert censored_row["passed"] == "True"
    assert "sw_knorm_ionization_e" in java
    assert '"0.00000000000000000e+00"' in java
    assert "sw_logeps_rate_ionization" in java
    assert "sw_logk_ionization_e" not in java


def test_monte_carlo_function_eedf_keeps_direct_rate_censoring_as_evidence(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    _enable_monte_carlo_fixture(mapping_path, bundle)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "  reaction_model: external_rates",
            "  reaction_model: function_eedf\n"
            "  function_eedf:\n"
            "    table: eedf_f0_comsol_2d.csv\n"
            "    function_tag: sw_eedf_monte_carlo\n"
            "    interpolation: structured_spreadsheet_linear_projection\n"
            "    extrapolation: constant",
        ),
        encoding="utf-8",
    )
    _set_fixture_censored_ionization(
        bundle,
        relevance_fraction=1.0e-5,
        upper_fraction_of_peak=0.5e-5,
    )

    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    closure = json.loads(plan.plan_json.read_text(encoding="utf-8"))["closure"]
    censoring = closure["mc_rate_censoring"]

    assert censoring["physics_qualification"]["passed"] is True
    assert censoring["censored_anchor_count"] == 1
    assert censoring["nominal_numerical_binding"] == {
        "passed": True,
        "estimator": "monte_carlo_trajectory_time_average_sigma_v",
        "role": "independent_function_eedf_tail_evidence",
        "active_comsol_rate_input": False,
        "zero_estimates_bound_exactly": True,
        "artificial_floor": False,
        "two_term_substitution": False,
    }


def _set_fixture_censored_ionization(
    bundle: Path,
    *,
    relevance_fraction: float,
    upper_fraction_of_peak: float,
    pooled_event_count: int = 0,
    exposure_scale: float = 1.0,
) -> None:
    rates_path = bundle / "rates_vs_mean_energy.csv"
    rates = list(csv.DictReader(rates_path.open(encoding="utf-8")))
    for row in rates:
        if float(row["mean_energy_eV"]) == 1.0 and row["process_type"] == "ionization":
            row["rate_coefficient_m3_s"] = "0"
    peak = max(
        float(row["rate_coefficient_m3_s"])
        for row in rates
        if row["process_type"] == "ionization"
    )
    upper = upper_fraction_of_peak * peak
    exposure = exposure_scale * -math.log(0.05) / upper
    _write_csv(
        rates_path,
        list(rates[0]),
        [[row[name] for name in rates[0]] for row in rates],
    )

    evidence_path = bundle / "rate_evidence.csv"
    evidence = list(csv.DictReader(evidence_path.open(encoding="utf-8")))
    for row in evidence:
        if float(row["mean_energy_eV"]) == 1.0 and row["process_type"] == "ionization":
            row.update(
                {
                    "rate_coefficient_mean_m3_s": "0",
                    "rate_coefficient_ci95_low_m3_s": "",
                    "rate_coefficient_ci95_high_m3_s": "",
                    "estimate_status": "censored_all_zero",
                    "uncertainty_available": "False",
                    "pooled_event_count": str(pooled_event_count),
                    "pooled_target_exposure_s_m3": str(exposure),
                    "pooled_all_zero_upper_95_m3_s": str(upper),
                    "pooled_zero_event_status": "all_zero_upper_95",
                }
            )
    _write_csv(
        evidence_path,
        list(evidence[0]),
        [[row[name] for name in evidence[0]] for row in evidence],
    )
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["tables"]["rates_vs_mean_energy.csv"]["sha256"] = sha256(
        rates_path.read_bytes()
    ).hexdigest()
    manifest["tables"]["rate_evidence.csv"]["sha256"] = sha256(
        evidence_path.read_bytes()
    ).hexdigest()
    manifest["quality_thresholds"]["required_rate_min_process_peak_fraction"] = (
        relevance_fraction
    )
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")


@pytest.mark.parametrize(
    ("relevance_fraction", "upper_fraction_of_peak"),
    ((1.0e-5, 1.0e-5), (1.0e-5, 2.0e-5), (0.0, 1.0e-8)),
)
def test_monte_carlo_censored_rate_equal_above_or_zero_fraction_fails_closed(
    tmp_path: Path,
    relevance_fraction: float,
    upper_fraction_of_peak: float,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    _enable_monte_carlo_fixture(mapping_path, bundle)
    _set_fixture_censored_ionization(
        bundle,
        relevance_fraction=relevance_fraction,
        upper_fraction_of_peak=upper_fraction_of_peak,
    )

    with pytest.raises(
        GecCcpWorkflowError,
        match="censored-rate qualification failed before COMSOL execution",
    ):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


@pytest.mark.parametrize(
    ("pooled_event_count", "exposure_scale"),
    ((1, 1.0), (0, 2.0)),
)
def test_monte_carlo_censored_rate_count_or_exposure_mismatch_fails_closed(
    tmp_path: Path,
    pooled_event_count: int,
    exposure_scale: float,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    _enable_monte_carlo_fixture(mapping_path, bundle)
    _set_fixture_censored_ionization(
        bundle,
        relevance_fraction=1.0e-5,
        upper_fraction_of_peak=0.5e-5,
        pooled_event_count=pooled_event_count,
        exposure_scale=exposure_scale,
    )

    with pytest.raises(
        GecCcpWorkflowError,
        match="lacks a self-consistent pooled 95% zero-event upper bound",
    ):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_monte_carlo_upstream_eedf_rate_audit_fails_on_trajectory_rate_drift(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    _enable_monte_carlo_fixture(mapping_path, bundle)
    rates_path = bundle / "rates_vs_mean_energy.csv"
    rates = list(csv.DictReader(rates_path.open(encoding="utf-8")))
    for row in rates:
        if row["process_type"] == "excitation":
            row["rate_coefficient_m3_s"] = str(
                1.2 * float(row["rate_coefficient_m3_s"])
            )
    _write_csv(
        rates_path,
        list(rates[0]),
        [[row[name] for name in rates[0]] for row in rates],
    )
    evidence_path = bundle / "rate_evidence.csv"
    evidence = list(csv.DictReader(evidence_path.open(encoding="utf-8")))
    for row in evidence:
        if row["process_type"] != "excitation":
            continue
        for column in (
            "rate_coefficient_mean_m3_s",
            "rate_coefficient_ci95_low_m3_s",
            "rate_coefficient_ci95_high_m3_s",
        ):
            row[column] = str(1.2 * float(row[column]))
    _write_csv(
        evidence_path,
        list(evidence[0]),
        [[row[name] for name in evidence[0]] for row in evidence],
    )
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["tables"]["rates_vs_mean_energy.csv"]["sha256"] = sha256(
        rates_path.read_bytes()
    ).hexdigest()
    manifest["tables"]["rate_evidence.csv"]["sha256"] = sha256(
        evidence_path.read_bytes()
    ).hexdigest()
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(
        GecCcpWorkflowError,
        match=("projected-EEDF/rate consistency audit failed: rate_excitation="),
    ):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_monte_carlo_upstream_eedf_rate_audit_fails_on_mph_xs_drift(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    _enable_monte_carlo_fixture(mapping_path, bundle)
    model_path = mapping_path.parent.parent / "argon_gec_ccp.mph"
    with ZipFile(model_path) as archive:
        model_xml = archive.read("dmodel.xml").decode("utf-8")
    excitation = re.search(
        r'(<PhysicsFeature op="ElectronImpactReaction" tag="eir2">.*?'
        r'<param param="ydata" value=")([^"]+)("/>)',
        model_xml,
        flags=re.DOTALL,
    )
    assert excitation is not None
    changed_values = excitation.group(2).replace("'0'", "'1e-99'", 1)
    model_xml = (
        model_xml[: excitation.start(2)]
        + changed_values
        + model_xml[excitation.end(2) :]
    )
    with ZipFile(model_path, "w") as archive:
        archive.writestr("dmodel.xml", model_xml)

    with pytest.raises(
        GecCcpWorkflowError,
        match="cross_section_identity=excitation",
    ):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_gec_ccp_contract_and_java_preserve_solved_mean_energy(
    tmp_path: Path,
) -> None:
    mapping, bundle = _write_fixture(tmp_path)

    plan = prepare_gec_ccp_run(mapping, bundle_path=bundle)
    source = generate_apply_java(plan.mapping)

    assert plan.contract.original_eedf == "Druyvesteyn"
    assert plan.contract.physics_operation == "ColdPlasmaTimePeriodic"
    assert set(plan.contract.reaction_features) == {"eir1", "eir2", "eir3"}
    assert '"SpecifyElectronDensityAndEnergy", "SpecifyAll"' in source
    assert "sw_logeps" in source
    assert "En_per-Ne_per" in source
    assert "sw_ne_safe" not in source
    assert "ptp.en/sw_ne_safe" not in source
    run_source = generate_run_java(
        plan.mapping,
        run_role="external",
        class_name="ExternalRun",
        input_mph=plan.mapping.model.output_mph,
        output_mph=plan.mapping.model.output_mph,
        study=plan.mapping.model.external_study,
        solution=plan.mapping.model.external_solution,
        time_periodic_feature=(plan.mapping.model.external_time_periodic_feature),
    )
    assert '.set("ptp.ebar",' in run_source
    assert "exp(En_per-Ne_per)*1[V]" in run_source
    assert '.set("ptp.Te",' in run_source
    assert "2*ptp.ebar/3" in run_source
    assert 'set("rstepabs"' not in run_source
    assert '.set("ptp.wAr_1p",' not in run_source
    assert "Ar_1p.weak$1" not in run_source
    assert "Ar_1p.weak$2" not in run_source
    assert "model.init().create(" not in run_source
    assert '.set("WAr_1p_per", ' not in run_source
    assert "native heavy-species equations" in run_source
    assert "nojac" not in run_source
    assert "/ptp.wAr_1p" not in run_source
    assert "/(1-ptp.wAr_1p)" not in run_source
    assert "sw_log_muN_e" in source
    assert "sw_log_DeN_e" in source
    assert "sw_log_DeN_L_e" not in source
    assert "sw_log_DeN_T_e" not in source
    assert "sw_log_muenN_e" in source
    assert "sw_log_DenN_e" in source
    assert "sw_log_DenN_L_e" not in source
    assert "sw_log_DenN_T_e" not in source
    assert "ptp.Er" not in source
    assert "ptp.Ez" not in source
    assert '.prop("ElectronProperties").set("ReducedProps", true)' in source
    assert '.prop("ElectronProperties").set("TensorElectronProps", false)' in source
    assert '.prop("ElectronProperties").set("IncludeThermalDiffusion", false)' in source
    assert '.set("muN",' in source
    assert '.set("DeN",' in source
    assert '.set("muenN",' in source
    assert '.set("DenN",' in source
    assert "/max(ptp.Nn,1[1/m^3])" not in source
    assert "nojac" not in source
    assert '"SourceStabilization", false' in source
    assert '"ReactionSourceStabilization", false' in source
    assert "SpecifyMeanElectronEnergy" not in source
    assert source.count('"RateConstantForm", "UseRate"') == 3
    assert source.count('"SpecifyReactionUsing", "RateConstant"') == 3
    assert source.count("6.02214076e23[1/mol]*exp(sw_logk_") == 3
    assert source.count('*1[m^3/s]"') == 3
    assert source.count('.set("interp", "piecewisecubic")') == 7
    assert "remapPeriodicLogEnergySolution" not in source
    assert "promotePeriodicSolution" not in source
    assert '"SpecifyReactionUsing", "UseLookupTable"' not in source
    assert '"SpecifyReactionUsing", "UseCrossSectionData"' not in source
    assert '"xratedata"' not in source
    assert '"yratedata"' not in source
    assert "xtownratedata" not in source
    assert "sw_mobility_blend" not in source
    run_source = generate_run_java(
        plan.mapping,
        run_role="baseline",
        class_name="Run",
        input_mph=plan.mapping.model.output_mph,
        output_mph=plan.mapping.model.output_mph,
    )
    assert 'model.param().set("P0", "1[W]")' in run_source
    assert run_source.count('model.study("std1").run()') == 1
    assert run_source.count('model.study("std2").run()') == 1
    assert "for (" not in run_source
    assert "blend" not in run_source
    assert '"dtech"' not in run_source
    assert "Ar_1p.weak" not in run_source
    assert "WAr_1p_per" not in run_source
    assert "swExtSeg" not in run_source
    assert "LowerLimit" not in run_source
    assert "UpperLimit" not in run_source
    assert '.feature("tper").set("useinitsol", false)' in run_source
    assert '.feature("tper").set("initmethod", "init")' in run_source
    assert '.feature("v1").set("initsol", "zero")' in run_source
    assert 'model.sol("sol1").clearSolutionData()' in run_source
    assert '"SourceStabilization", false' in run_source
    assert '"ReactionSourceStabilization", false' in run_source
    external_source = plan.external_run_java.read_text(encoding="utf-8")
    assert 'model.study("std3").run()' in external_source
    assert '.feature("tptd").set("notstudy", "std3")' in external_source
    assert '.feature("tptd").set("notstudystep", "tper1")' in external_source
    assert '.feature("fc1").set("dtech"' not in external_source
    assert '.feature("fc1").set("initsteph"' not in external_source
    assert '.feature("fc1").set("minsteph"' not in external_source
    assert '.feature("fc1").set("rstepabs"' not in external_source
    assert '.featureInfo("info").set("ptp.Mn",' in external_source
    assert 'new String[]{"0.04[kg/mol]"}' in external_source
    assert "log(ptp.ebar/1[V])" in external_source
    physical_argument = external_source.index(
        "the same closure remains defined on std2 physical-time data"
    )
    external_restore = external_source.index(
        "Restore COMSOL's native definitions before std2"
    )
    external_conversion = external_source.index('model.study("std2").run();')
    refreshed_definitions = external_source.index(
        "model.sol(solutionTag).updateSolution();"
    )
    assert (
        external_restore
        < physical_argument
        < refreshed_definitions
        < external_conversion
    )
    assert external_source.count("model.sol(solutionTag).updateSolution();") == 1
    for generated, periodic_study in (
        (run_source, "std1"),
        (external_source, "std3"),
    ):
        periodic_run = generated.index(f'model.study("{periodic_study}").run();')
        restore = generated.index("Restore COMSOL's native definitions before std2")
        conversion_run = generated.index('model.study("std2").run();')
        final_save = generated.rindex("model.save(")
        assert periodic_run < restore < conversion_run < final_save
        assert generated.count('.removeLock("ptp.Mn")') == 1
        assert generated.count('.removeLock("ptp.ebar")') == 1
        assert generated.count('.removeLock("ptp.Te")') == 1
        assert "Ar_1p.weak" not in generated
    external_export = plan.external_export_java.read_text(encoding="utf-8")
    assert '.set("data", "dset4")' in external_export
    assert '"data", "dset5"' in external_export
    assert 'export().create("swDomainAvg", "Data")' in external_export
    assert '"WAr_1p_per"' not in external_export
    assert "domain_period_average.csv" in external_export
    assert any(
        path.name == "domain_period_average.csv" for path in plan.expected_result_files
    )
    assert any(
        path.name == "conservation_volume.csv" for path in plan.expected_result_files
    )
    assert 'numerical().create("swConservationVolume", "IntSurface")' in (
        external_export
    )
    assert '"intvolume", true' in external_export
    assert 'numerical().create("swConservationWall", "IntLine")' in (external_export)
    assert '"intsurface", true' in external_export
    assert "2*pi*r" not in external_export
    assert "ptp.ne_bnd+ptp.tne_bnd" in external_export
    assert "ptp.en_bnd+ptp.ten_bnd" in external_export
    assert "ptp.seflux" in external_export
    assert "down(ptp.gflux_ne_avr)" in external_export
    assert "down(ptp.gflux_enr)" in external_export
    assert "ptp.R_wAr_1p_av" in external_export
    baseline_export = plan.baseline_export_java.read_text(encoding="utf-8")

    def domain_phase_expressions(java: str) -> list[str]:
        line = next(
            item
            for item in java.splitlines()
            if 'export("swDomainPhaseClosure").set("expr"' in item
        )
        payload = line.split("new String[]{", maxsplit=1)[1].rsplit("});", maxsplit=1)[
            0
        ]
        return json.loads(f"[{payload}]")

    external_domain = domain_phase_expressions(external_export)
    baseline_domain = domain_phase_expressions(baseline_export)
    active_transport = {
        "ptp.muerr",
        "ptp.muezr",
        "ptp.muerz",
        "ptp.muezz",
        "ptp.Derr",
        "ptp.Dezr",
        "ptp.Derz",
        "ptp.Dezz",
        "ptp.muenrr",
        "ptp.muenzr",
        "ptp.muenrz",
        "ptp.muenzz",
        "ptp.Denrr",
        "ptp.Denzr",
        "ptp.Denrz",
        "ptp.Denzz",
    }
    assert active_transport <= set(external_domain)
    assert active_transport.isdisjoint(baseline_domain)
    assert set(external_domain) - set(baseline_domain) == active_transport
    assert plan.plan_json.exists()
    plan_data = json.loads(plan.plan_json.read_text(encoding="utf-8"))
    assert plan_data["closure"]["electron_transport"] == (
        "comsol_specify_all_restricted"
    )
    assert plan_data["closure"]["zero_field_isotropization_Td"] is None
    tensor = plan_data["closure"]["transport_tensor"]
    assert tensor["thermal_diffusion"] is False
    assert tensor["gradient_closure"] == "comsol_specify_all_restricted"
    assert tensor["full_gradient_response_identified"] is False
    assert "standard_local_energy_restricted" in tensor["scope"]
    response = plan_data["closure"]["gradient_response"]
    assert response["comsol_executable_eligible"] is True
    assert response["restricted_closure_eligible"] is True
    assert response["full_physical_eligible"] is False
    assert response["full_gradient_response_identified"] is False
    assert response["response_matrix_assumptions"] == {
        "A_nn": "De",
        "A_nE": "0",
        "A_un": "Den",
        "A_uE": "Den",
    }
    assert plan_data["closure"]["reaction_model"] == "external_rates"
    assert plan_data["closure"]["source_field"] == "steady_dc"
    assert plan_data["closure"]["lookup_jacobian"] == "exact"
    solve = plan_data["solve"]
    assert solve["include_builtin_reference"] is True
    assert solve["native_physics_initial_values"] is True
    assert solve["saved_solution_dependency"] is False
    assert solve["power_sweep"] is False
    assert solve["coefficient_sweep"] is False
    assert solve["baseline_nonlinear_method"] == "model_defined"
    assert solve["external_nonlinear_method"] == "model_defined"
    assert solve["nonlinear_solver_overrides"] is False
    assert solve["nonlinear_globalization"]["scope"] == [
        "external",
        "baseline",
    ]
    assert solve["periodic_solver_exact_identities"]["scope"] == [
        "external",
        "baseline",
    ]
    assert plan_data["bindings"]["external_solution"] == "sol3"
    assert plan_data["bindings"]["datasets"] == {
        "period_baseline": "dset1",
        "period_external": "dset4",
        "phase": "dset3",
        "waveform_baseline": "dset2",
        "waveform_external": "dset5",
    }
    assert len(plan_data["mapping"]["sha256"]) == 64
    assert len(plan_data["model"]["input_mph"]["sha256"]) == 64
    assert plan_data["model"]["input_mph"]["size_bytes"] > 0
    assert all(
        len(metadata["sha256"]) == 64
        for metadata in plan_data["generated_java"].values()
    )
    assert all(
        metadata["materialized"] is True
        for metadata in plan_data["generated_java"].values()
    )
    plan_inputs = validate_gec_plan_inputs(plan)
    assert plan_inputs["bundle_artifacts_verified"] is True
    assert len(plan_inputs["generated_java_sha256"]) == 5
    assert len(plan_data["bundle"]["manifest_sha256"]) == 64
    frozen_artifacts = plan_data["bundle"]["artifact_verification"]["sha256"]
    bundle_manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    assert set(frozen_artifacts) == set(bundle_manifest["tables"])
    assert plan_inputs["bundle_artifact_sha256"] == frozen_artifacts
    assert set(plan_data["bundle"]["hashes"]) == {
        "workflow_config_sha256",
        "base_config_sha256",
        "cross_sections_sha256",
    }
    assert plan_data["bundle"]["quality_thresholds"] == (
        quality_thresholds_payload(QualityThresholds())
    )


def test_gec_ccp_nonwriting_prepare_does_not_claim_stale_java(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    materialized = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    materialized.external_run_java.write_text(
        "// stale generated source\n", encoding="utf-8"
    )

    plan = prepare_gec_ccp_run(
        mapping_path,
        bundle_path=bundle,
        write_java=False,
    )
    plan_data = json.loads(plan.plan_json.read_text(encoding="utf-8"))

    assert plan.external_run_java.read_text(encoding="utf-8") == (
        "// stale generated source\n"
    )
    assert all(
        metadata["materialized"] is False
        for metadata in plan_data["generated_java"].values()
    )
    assert all(
        metadata["sha256"] is None for metadata in plan_data["generated_java"].values()
    )
    with pytest.raises(
        GecCcpWorkflowError,
        match="Java artifact is not materialized",
    ):
        validate_gec_plan_inputs(plan)


def test_gec_ccp_preflight_mismatch_replaces_stale_completed_status(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    status_path = plan.output_directory / "gec_ccp_run_status.json"
    status_path.write_text('{"status":"completed"}', encoding="utf-8")
    plan.apply_java.write_text(
        plan.apply_java.read_text(encoding="utf-8") + "// tampered\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        execution_preflight,
        "prepare_gec_ccp_run",
        lambda *args, **kwargs: plan,
    )

    with pytest.raises(
        GecCcpWorkflowError,
        match="generated Java changed",
    ):
        execution.execute_gec_ccp_run(mapping_path, bundle_path=bundle)

    status = json.loads(status_path.read_text(encoding="utf-8"))
    assert status["status"] == "blocked"
    assert status["blocker_code"] == "plan_input_provenance_mismatch"
    assert status["solve_status"] == "not_started"
    assert status["gradient_response"]["scalar_diffusion_audit"]["passed"] is True


def test_gec_ccp_plan_rejects_bundle_artifact_changed_after_prepare(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    rates_path = bundle / "rates_vs_mean_energy.csv"
    rates_path.write_text(
        rates_path.read_text(encoding="utf-8") + "\n",
        encoding="utf-8",
    )

    with pytest.raises(
        GecCcpWorkflowError,
        match=(
            "bundle artifact changed after GEC plan generation: "
            "rates_vs_mean_energy.csv"
        ),
    ):
        validate_gec_plan_inputs(plan)


def test_gec_ccp_prepare_failure_replaces_stale_completed_status(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping = load_gec_ccp_mapping(mapping_path, bundle_path=bundle)
    status_path = mapping.results.output_directory / "gec_ccp_run_status.json"
    status_path.parent.mkdir(parents=True, exist_ok=True)
    status_path.write_text(
        '{"status":"completed","quality_accepted":true}',
        encoding="utf-8",
    )
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["status"] = "failed"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(
        GecCcpWorkflowError,
        match="bundle manifest status is not ok",
    ):
        execution.execute_gec_ccp_run(mapping_path, bundle_path=bundle)

    status = json.loads(status_path.read_text(encoding="utf-8"))
    assert status["status"] == "blocked"
    assert status["blocker_code"] == "preflight_validation_failed"
    assert status["solve_status"] == "not_started"
    assert status["quality_status"] == "not_evaluated"
    assert status["quality_accepted_for_declared_closure"] is None
    assert status["physical_target_accepted"] is False


def test_gec_ccp_revalidates_bundle_after_apply_before_external_solve(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    monkeypatch.setattr(
        execution_preflight,
        "prepare_gec_ccp_run",
        lambda *args, **kwargs: plan,
    )
    operations: list[str] = []

    def mutate_after_apply(
        *args: object,
        operation: str,
        **kwargs: object,
    ) -> object:
        operations.append(operation)
        if operation == "apply":
            rates_path = bundle / "rates_vs_mean_energy.csv"
            rates_path.write_text(
                rates_path.read_text(encoding="utf-8") + "\n",
                encoding="utf-8",
            )
        return SimpleNamespace(operation=operation)

    monkeypatch.setattr(
        execution_solver,
        "execute_generated_comsol_java",
        mutate_after_apply,
    )

    with pytest.raises(
        GecCcpWorkflowError,
        match=(
            "bundle artifact changed after GEC plan generation: "
            "rates_vs_mean_energy.csv"
        ),
    ):
        execution.execute_gec_ccp_run(mapping_path, bundle_path=bundle)

    assert operations == ["apply"]
    status = json.loads(
        (plan.output_directory / "gec_ccp_run_status.json").read_text(encoding="utf-8")
    )
    assert status["status"] == "blocked"
    assert status["blocker_code"] == "plan_input_provenance_mismatch"
    assert status["solve_status"] == "external_not_started"
    assert status["completed_operations"] == operations


def test_gec_ccp_rejects_mc_provenance_on_two_term_bundle(tmp_path: Path) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["hashes"].update(
        {
            "mc_sampling_plan_json": "[]",
            "mc_sampling_plan_sha256": "d" * 64,
            "mc_transport_estimator_schema_version": "direct_mc_transport.v3",
        }
    )
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(
        GecCcpWorkflowError,
        match="must not contain MC sampling provenance",
    ):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_gec_ccp_completed_solve_fails_closed_when_quality_is_rejected(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    monkeypatch.setattr(
        execution_preflight,
        "prepare_gec_ccp_run",
        lambda *args, **kwargs: plan,
    )

    def fake_execute(*args: object, operation: str, **kwargs: object) -> object:
        for result in plan.expected_result_files:
            result.write_text("% synthetic test result\n0,0\n", encoding="utf-8")
        plan.mapping.model.output_mph.write_bytes(b"saved external mph")
        return SimpleNamespace(operation=operation)

    monkeypatch.setattr(execution_solver, "execute_generated_comsol_java", fake_execute)
    monkeypatch.setattr(
        execution_runtime,
        "_validated_comsol_runtime",
        lambda *args, **kwargs: {
            "version": "6.4",
            "build": "429",
            "executable": "comsolbatch.exe",
            "input_model_upgrade_required": False,
        },
    )
    monkeypatch.setattr(
        execution_runtime,
        "_validated_saved_mph_versions",
        lambda *args, **kwargs: {"baseline": "6.4", "external": "6.4"},
    )

    def failed_conservation_audit(*args: object, **kwargs: object) -> Path:
        audit = plan.output_directory / "conservation_audit.json"
        audit.write_text(json.dumps({"status": "failed"}), encoding="utf-8")
        return audit

    monkeypatch.setattr(
        execution_postsolve,
        "audit_gec_ccp_conservation_run",
        failed_conservation_audit,
    )

    def passed_transport_audit(*args: object, **kwargs: object) -> Path:
        audit = plan.output_directory / "transport_audit.json"
        audit.write_text(json.dumps({"status": "passed"}), encoding="utf-8")
        return audit

    monkeypatch.setattr(
        execution_postsolve,
        "audit_gec_ccp_transport_run",
        passed_transport_audit,
    )

    with pytest.raises(GecCcpQualityError, match="physics-quality acceptance"):
        execution.execute_gec_ccp_run(mapping_path, bundle_path=bundle)

    status = json.loads(
        (plan.output_directory / "gec_ccp_run_status.json").read_text(encoding="utf-8")
    )
    assert status["status"] == "rejected"
    assert status["solve_status"] == "completed"
    assert status["quality_status"] == "failed"
    assert status["quality_accepted"] is False
    assert status["quality_acceptance_scope"] == ("declared_restricted_closure")
    assert status["quality_accepted_for_declared_closure"] is False
    assert status["full_gradient_response_identified"] is False
    assert status["full_physical_eligible"] is False
    assert (
        status["gradient_response"]["independent_2x2_gradient_response_identified"]
        is False
    )
    assert status["rejection_code"] == "physics_quality_failed"
    assert [item["path"] for item in status["results"]] == sorted(
        path.resolve().relative_to(plan.output_directory.resolve()).as_posix()
        for path in plan.expected_result_files
    )
    assert all(item["size_bytes"] > 0 for item in status["results"])
    assert all(len(item["sha256"]) == 64 for item in status["results"])


def test_gec_ccp_diagnostic_acceptance_is_not_physical_target_acceptance(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "role: physical_target",
            "role: diagnostic_control",
        ),
        encoding="utf-8",
    )
    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    plan_payload = json.loads(plan.plan_json.read_text(encoding="utf-8"))
    assert plan_payload["result_role"] == "diagnostic_control"
    assert plan_payload["physical_target_eligible"] is False
    monkeypatch.setattr(
        execution_preflight,
        "prepare_gec_ccp_run",
        lambda *args, **kwargs: plan,
    )

    def fake_execute(*args: object, operation: str, **kwargs: object) -> object:
        for result in plan.expected_result_files:
            result.write_text(
                "% synthetic test result\n0,0\n",
                encoding="utf-8",
            )
        plan.mapping.model.output_mph.write_bytes(b"saved external mph")
        return SimpleNamespace(operation=operation)

    monkeypatch.setattr(
        execution_solver,
        "execute_generated_comsol_java",
        fake_execute,
    )
    monkeypatch.setattr(
        execution_runtime,
        "_validated_comsol_runtime",
        lambda *args, **kwargs: {
            "version": "6.4",
            "build": "429",
            "executable": "comsolbatch.exe",
            "input_model_upgrade_required": False,
        },
    )
    monkeypatch.setattr(
        execution_runtime,
        "_validated_saved_mph_versions",
        lambda *args, **kwargs: {"baseline": "6.4", "external": "6.4"},
    )

    def passed_audit(name: str) -> Path:
        audit = plan.output_directory / name
        audit.write_text(json.dumps({"status": "passed"}), encoding="utf-8")
        return audit

    monkeypatch.setattr(
        execution_postsolve,
        "audit_gec_ccp_conservation_run",
        lambda *args: passed_audit("conservation_audit.json"),
    )
    monkeypatch.setattr(
        execution_postsolve,
        "audit_gec_ccp_transport_run",
        lambda *args: passed_audit("transport_audit.json"),
    )

    summary = execution.execute_gec_ccp_run(
        mapping_path,
        bundle_path=bundle,
    )

    status = json.loads(summary.status_json.read_text(encoding="utf-8"))
    assert status["status"] == "completed"
    assert status["quality_accepted"] is True
    assert status["quality_accepted_for_declared_closure"] is True
    assert status["result_role"] == "diagnostic_control"
    assert status["physical_target_accepted"] is False


def test_gec_ccp_external_rates_without_transport_run_support_audit(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "electron_transport: comsol_specify_all_restricted",
            "electron_transport: comsol",
        ),
        encoding="utf-8",
    )
    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    monkeypatch.setattr(
        execution_preflight,
        "prepare_gec_ccp_run",
        lambda *args, **kwargs: plan,
    )

    def fake_execute(*args: object, operation: str, **kwargs: object) -> object:
        for result in plan.expected_result_files:
            result.write_text("% synthetic test result\n0,0\n", encoding="utf-8")
        plan.mapping.model.output_mph.write_bytes(b"saved external mph")
        return SimpleNamespace(operation=operation)

    monkeypatch.setattr(execution_solver, "execute_generated_comsol_java", fake_execute)
    monkeypatch.setattr(
        execution_runtime,
        "_validated_comsol_runtime",
        lambda *args, **kwargs: {
            "version": "6.4",
            "build": "429",
            "executable": "comsolbatch.exe",
            "input_model_upgrade_required": False,
        },
    )
    monkeypatch.setattr(
        execution_runtime,
        "_validated_saved_mph_versions",
        lambda *args, **kwargs: {"baseline": "6.4", "external": "6.4"},
    )

    def passed_audit(name: str, payload: dict[str, object]) -> Path:
        audit = plan.output_directory / name
        audit.write_text(json.dumps(payload), encoding="utf-8")
        return audit

    monkeypatch.setattr(
        execution_postsolve,
        "audit_gec_ccp_conservation_run",
        lambda *args: passed_audit("conservation_audit.json", {"status": "passed"}),
    )
    support_calls: list[object] = []

    def passed_support_audit(candidate: object) -> Path:
        support_calls.append(candidate)
        return passed_audit(
            "closure_support_audit.json",
            {"status": "passed", "support": {"passed": True}},
        )

    monkeypatch.setattr(
        execution_postsolve,
        "audit_gec_ccp_closure_support_run",
        passed_support_audit,
    )

    summary = execution.execute_gec_ccp_run(mapping_path, bundle_path=bundle)

    status = json.loads(summary.status_json.read_text(encoding="utf-8"))
    assert support_calls == [plan]
    assert status["status"] == "completed"
    assert status["closure_support_quality_status"] == "passed"
    assert status["external_closure_support_qualification"] == {"passed": True}
    assert status["transport_quality_status"] == "not_applicable"


def test_gec_ccp_function_eedf_requires_native_saved_function_audit(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "reaction_model: external_rates",
            "reaction_model: function_eedf\n"
            "  interpolation: differentiable_log_piecewise_cubic\n"
            "  function_eedf:\n"
            "    table: eedf_f0_comsol_2d.csv\n"
            "    function_tag: sw_eedf_test\n"
            "    interpolation: structured_spreadsheet_linear_projection\n"
            "    extrapolation: constant",
        ),
        encoding="utf-8",
    )
    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    assert plan.native_eedf_audit is not None
    monkeypatch.setattr(
        execution_preflight,
        "prepare_gec_ccp_run",
        lambda *args, **kwargs: plan,
    )

    def fake_execute(*args: object, operation: str, **kwargs: object) -> object:
        for result in plan.expected_result_files:
            result.write_text("% synthetic test result\n0,0\n", encoding="utf-8")
        plan.mapping.model.output_mph.write_bytes(b"saved external mph")
        return SimpleNamespace(
            operation=operation,
            stdout_paths=(tmp_path / f"{operation}.txt",),
        )

    monkeypatch.setattr(execution_solver, "execute_generated_comsol_java", fake_execute)
    monkeypatch.setattr(
        execution_solver, "extract_comsol_eedf_audit_log", lambda *args: None
    )
    monkeypatch.setattr(
        execution_solver,
        "analyze_comsol_eedf_audit",
        lambda *args, **kwargs: {"passed": False, "reason": "binding mismatch"},
    )
    monkeypatch.setattr(
        execution_runtime,
        "_validated_comsol_runtime",
        lambda *args, **kwargs: {
            "version": "6.4",
            "build": "429",
            "executable": "comsolbatch.exe",
            "input_model_upgrade_required": False,
        },
    )
    monkeypatch.setattr(
        execution_runtime,
        "_validated_saved_mph_versions",
        lambda *args, **kwargs: {"baseline": "6.4", "external": "6.4"},
    )

    def audit_file(name: str) -> Path:
        path = plan.output_directory / name
        path.write_text(json.dumps({"status": "passed"}), encoding="utf-8")
        return path

    monkeypatch.setattr(
        execution_postsolve,
        "audit_gec_ccp_conservation_run",
        lambda *args: audit_file("conservation_audit.json"),
    )
    monkeypatch.setattr(
        execution_postsolve,
        "audit_gec_ccp_function_eedf_run",
        lambda *args: audit_file("function_eedf_audit.json"),
    )
    monkeypatch.setattr(
        execution_postsolve,
        "audit_gec_ccp_transport_run",
        lambda *args: audit_file("transport_audit.json"),
    )

    with pytest.raises(GecCcpQualityError, match="native_function_eedf"):
        execution.execute_gec_ccp_run(mapping_path, bundle_path=bundle)

    status = json.loads(
        (plan.output_directory / "gec_ccp_run_status.json").read_text(encoding="utf-8")
    )
    assert status["status"] == "rejected"
    assert status["native_function_eedf_quality_status"] == "failed"
    assert "gec_native_eedf_audit" in status["operations"]


def test_gec_ccp_keeps_native_eedf_audit_when_external_compile_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "reaction_model: external_rates",
            "reaction_model: function_eedf\n"
            "  interpolation: differentiable_log_piecewise_cubic\n"
            "  function_eedf:\n"
            "    table: eedf_f0_comsol_2d.csv\n"
            "    function_tag: sw_eedf_test\n"
            "    interpolation: structured_spreadsheet_linear_projection\n"
            "    extrapolation: constant",
        ),
        encoding="utf-8",
    )
    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    assert plan.native_eedf_audit is not None
    monkeypatch.setattr(
        execution_preflight,
        "prepare_gec_ccp_run",
        lambda *args, **kwargs: plan,
    )

    def fake_execute(*args: object, operation: str, **kwargs: object) -> object:
        if operation == "gec_external_run":
            failure_log = tmp_path / "external_failure.txt"
            failure_log.write_text(
                "Unknown property.\n - Property: nargs\n", encoding="utf-8"
            )
            raise ComsolAdapterError(
                "external compile failed",
                step_result={"stdout": str(failure_log)},
            )
        stdout = tmp_path / f"{operation}.txt"
        stdout.write_text("synthetic\n", encoding="utf-8")
        return SimpleNamespace(operation=operation, stdout_paths=(stdout,))

    monkeypatch.setattr(execution_solver, "execute_generated_comsol_java", fake_execute)
    monkeypatch.setattr(
        execution_solver, "extract_comsol_eedf_audit_log", lambda *args: None
    )
    monkeypatch.setattr(
        execution_solver,
        "analyze_comsol_eedf_audit",
        lambda *args, **kwargs: {"passed": True},
    )

    with pytest.raises(ComsolAdapterError, match="compile failed"):
        execution.execute_gec_ccp_run(mapping_path, bundle_path=bundle)

    audit_path = plan.native_eedf_audit.values_path.parent / "comsol_eedf_audit.json"
    assert json.loads(audit_path.read_text(encoding="utf-8"))["status"] == "passed"
    status = json.loads(
        (plan.output_directory / "gec_ccp_run_status.json").read_text(encoding="utf-8")
    )
    assert status["blocker_code"] == "comsol_execution_failed"
    assert status["native_function_eedf_quality_status"] == "passed"
    assert status["native_function_eedf_audit"] == str(audit_path)


@pytest.mark.parametrize(
    ("small_rate", "expected_pass", "significant_samples"),
    ((4.04e-6, True, 1), (2.0e-5, False, 2)),
)
def test_function_eedf_reintegration_uses_magnitude_aware_rate_quality(
    small_rate: float,
    expected_pass: bool,
    significant_samples: int,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source_rates = {1.0: 1.0, 2.0: small_rate}
    comsol_rates = np.asarray([1.0, source_rates[2.0] * 1.152])
    monkeypatch.setattr(
        function_run_audit,
        "_integrate_native_function_eedf_rate",
        lambda grid, mean, energy, sigma: source_rates[float(mean)],
    )

    result = function_run_audit._audit_function_eedf_rates(
        object(),
        {"excitation": (np.asarray([0.0]), np.asarray([0.0]))},
        (contracts.GecReactionSpec("excitation", "eir2", "excitation"),),
        np.asarray([[1.0, 2.0]]),
        (comsol_rates * 6.02214076e23).reshape(2, 1),
        {"eir2": [0]},
    )

    process = result["processes"]["excitation"]
    assert process["passed"] is expected_pass
    assert process["normalized_rmse"] < 0.02
    assert process["significance_fraction_of_process_peak"] == 1.0e-5
    assert process["significant_samples"] == significant_samples
    assert process["record_only_samples"] == 2 - significant_samples
    assert process["below_significance_floor_policy"] == (
        "record_only_for_pointwise_relative_error"
    )


@pytest.mark.parametrize("expected_source", ["two_term", "monte_carlo"])
def test_function_eedf_source_consistency_uses_raw_solver_rates(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    expected_source: str,
) -> None:
    source_rates = {1.0: 1.0, 2.0: 4.04e-6}
    derived_rates = {1.0: 1.0, 2.0: source_rates[2.0] * 1.152}
    _write_csv(
        tmp_path / "rates_vs_mean_energy.csv",
        [
            "E_over_N_Td",
            "mean_energy_eV",
            "process_type",
            "rate_coefficient_m3_s",
        ],
        [
            [10.0, 1.0, "excitation", source_rates[1.0]],
            [20.0, 2.0, "excitation", source_rates[2.0]],
        ],
    )
    monkeypatch.setattr(
        eedf_validation,
        "_integrate_native_function_eedf_rate",
        lambda grid, mean, energy, sigma: derived_rates[float(mean)],
    )
    mapping = SimpleNamespace(
        bundle=SimpleNamespace(path=tmp_path, expected_source=expected_source),
        results=SimpleNamespace(output_directory=tmp_path),
    )

    result = eedf_validation._audit_function_eedf_source_rates(
        mapping,
        object(),
        {"excitation": (np.asarray([0.0]), np.asarray([0.0]))},
        (contracts.GecReactionSpec("excitation", "eir2", "excitation"),),
    )

    process = result["processes"]["excitation"]
    assert process["passed"] is True
    assert result["source"] == f"{expected_source}_rates"
    assert process["normalized_rmse"] < 0.02
    assert process["p95_relative_error_significant"] == pytest.approx(0.0)
    assert process["p95_relative_error_limit"] == pytest.approx(0.05)
    assert process["significance_fraction_of_process_peak"] == 1.0e-5
    assert process["significant_samples"] == 1
    with (tmp_path / "function_eedf_rate_audit.csv").open(
        encoding="utf-8", newline=""
    ) as handle:
        rows = list(csv.DictReader(handle))
    assert rows[1]["relevant"] == "False"
    assert rows[1]["criterion"] == "below_1e-5_process_peak_record_only"
    assert rows[1]["passed"] == "True"


def test_function_eedf_source_consistency_fails_on_relevant_p95(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    means = np.arange(1.0, 21.0)
    _write_csv(
        tmp_path / "rates_vs_mean_energy.csv",
        [
            "E_over_N_Td",
            "mean_energy_eV",
            "process_type",
            "rate_coefficient_m3_s",
        ],
        [[mean, mean, "excitation", 1.0] for mean in means],
    )
    monkeypatch.setattr(
        eedf_validation,
        "_integrate_native_function_eedf_rate",
        lambda _grid, mean, _energy, _sigma: 1.06 if mean <= 2.0 else 1.0,
    )
    mapping = SimpleNamespace(
        bundle=SimpleNamespace(path=tmp_path, expected_source="two_term"),
        results=SimpleNamespace(output_directory=tmp_path),
    )

    result = eedf_validation._audit_function_eedf_source_rates(
        mapping,
        object(),
        {"excitation": (np.asarray([0.0]), np.asarray([0.0]))},
        (contracts.GecReactionSpec("excitation", "eir2", "excitation"),),
        audit_csv_name="upstream_eedf_rate_audit.csv",
    )

    process = result["processes"]["excitation"]
    assert process["maximum_relative_error_significant"] < 0.10
    assert process["normalized_rmse"] < 0.02
    assert process["p95_relative_error_significant"] > 0.05
    assert process["passed"] is False
    assert result["passed"] is False


def test_rate_significance_keeps_the_1e_minus_5_boundary_in_quality_gate() -> None:
    _, floor, significant = eedf_validation._rate_significance(
        np.asarray([1.0, 0.999e-5, 1.0e-5]),
        np.asarray([1.0, 0.999e-5, 1.0e-5]),
    )

    assert floor == pytest.approx(1.0e-5)
    assert significant.tolist() == [True, False, True]


def test_domain_phase_state_audit_requires_finite_positive_states(
    tmp_path: Path,
) -> None:
    path = tmp_path / "domain_phase_closure.csv"
    path.write_text(
        "% R,Z,ptp.ne (1/m^3) @ t=0,ptp.ebar (V) @ t=0,"
        "ptp.wAr_1p (1) @ t=0\n"
        "0,0,1e12,4.5,1e-8\n"
        "0.1,0.2,2e12,6.0,2e-8\n",
        encoding="utf-8",
    )

    accepted = conservation_audit._audit_domain_phase_state(path)

    assert accepted["passed"] is True
    assert accepted["quantities"]["electron_density_m3"]["spatial_rows"] == 2

    path.write_text(
        "% R,Z,ptp.ne (1/m^3) @ t=0,ptp.ebar (V) @ t=0,"
        "ptp.wAr_1p (1) @ t=0\n"
        "0,0,0,4.5,1e-8\n",
        encoding="utf-8",
    )
    rejected = conservation_audit._audit_domain_phase_state(path)

    assert rejected["passed"] is False
    assert rejected["quantities"]["electron_density_m3"]["passed"] is False


def test_repository_gec_ccp_model_reports_druyvesteyn() -> None:
    model = Path(__file__).parents[1] / "comsol_modes" / "argon_gec_ccp.mph"
    contract = inspect_gec_ccp_mph(model)

    assert contract.original_eedf == "Druyvesteyn"
    assert contract.physics_tag == "ptp"
    assert contract.plasma_feature == "pes1"
    assert contract.common_heavy_species_molar_mass_kg_mol == 0.04
    assert dict(contract.heavy_species_molar_masses_kg_mol) == {
        "Ar": 0.04,
        "Ars": 0.04,
        "Ar_1p": 0.04,
    }
    assert {species.tag for species in contract.heavy_species if species.enabled} == {
        "Ar",
        "Ar_1p",
    }
    assert contract.heavy_species_selection == "BaseGeometry"
    assert contract.axisymmetric is True
    assert contract.heavy_species_formulation == "FEMLogLinear"
    assert contract.heavy_species_diffusion_model == "MixtureAveraged"
    assert contract.heavy_species_migration is True
    assert contract.heavy_species_convection is False
    assert contract.mixture_diffusion_correction is False
    assert contract.ion_tensor_properties is False
    assert contract.ion_electric_field_time_model == "Instantaneous"
    assert [
        (reaction.tag, reaction.formula)
        for reaction in contract.surface_reactions
        if reaction.enabled
    ] == [("sr1", "Ar+=>Ar")]


def test_gec_ccp_rejects_nonbinary_active_heavy_species(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    model = tmp_path / "argon_gec_ccp.mph"
    with ZipFile(model) as archive:
        xml = archive.read("dmodel.xml").decode("utf-8")
    xml = xml.replace("    <entityFlags>DISABLED</entityFlags>\n", "", 1)
    with ZipFile(model, "w") as archive:
        archive.writestr("dmodel.xml", xml)

    with pytest.raises(
        GecCcpWorkflowError,
        match="active heavy species must be exactly binary Ar and Ar_1p",
    ):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


@pytest.mark.parametrize(
    ("old", "new", "message"),
    (
        ("FEMLogLinear", "FEMLinear", "heavy-species formulation"),
        ("MixtureAveraged", "Fick", "heavy-species diffusion model"),
        (
            'param="Migration" value="1|1,\'1\'"',
            'param="Migration" value="1|1,\'0\'"',
            "heavy-species migration",
        ),
        (
            'param="Convection" value="1|1,\'0\'"',
            'param="Convection" value="1|1,\'1\'"',
            "heavy-species convection",
        ),
    ),
)
def test_gec_ccp_rejects_incompatible_heavy_transport_contract(
    tmp_path: Path, old: str, new: str, message: str
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    model = tmp_path / "argon_gec_ccp.mph"
    with ZipFile(model) as archive:
        xml = archive.read("dmodel.xml").decode("utf-8")
    assert old in xml
    with ZipFile(model, "w") as archive:
        archive.writestr("dmodel.xml", xml.replace(old, new, 1))

    with pytest.raises(GecCcpWorkflowError, match=message):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_gec_ccp_rejects_molar_mass_identity_lock_for_mixed_masses(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    model = tmp_path / "argon_gec_ccp.mph"
    with ZipFile(model) as archive:
        xml = archive.read("dmodel.xml").decode("utf-8")
    xml = xml.replace(
        '<PhysicsFeature op="Species" tag="Ar_1p">\n'
        '    <param param="sType" value="1|1,\'ion\'"/>\n'
        '    <param param="M" value="1|1,\'0.04[kg/mol]\'"/>',
        '<PhysicsFeature op="Species" tag="Ar_1p">\n'
        '    <param param="sType" value="1|1,\'ion\'"/>\n'
        '    <param param="M" value="1|1,\'0.05[kg/mol]\'"/>',
    )
    with ZipFile(model, "w") as archive:
        archive.writestr("dmodel.xml", xml)

    with pytest.raises(
        GecCcpWorkflowError,
        match="heavy-species molar masses are not uniformly",
    ):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


@pytest.mark.parametrize(
    ("thermal_model", "thermal_flag"),
    (
        ("off_restricted_diagonal", "0"),
        ("comsol_grad_diffusivity", "1"),
    ),
)
def test_function_eedf_transport_audit_checks_saved_binding_and_values(
    tmp_path: Path,
    thermal_model: str,
    thermal_flag: str,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "thermal_diffusion_model: off_restricted_diagonal",
            f"thermal_diffusion_model: {thermal_model}",
        ),
        encoding="utf-8",
    )
    mapping = load_gec_ccp_mapping(mapping_path, bundle_path=bundle)
    feature_params = [
        '<param param="SpecifyElectronDensityAndEnergy" '
        "value=\"1|1,'SpecifyAll'\"></param>"
    ]
    headers = [
        "R",
        "Z",
        "ptp.ebar (V) @ t=0",
        "ptp.Nn (1/m^3) @ t=0",
    ]
    row = ["0", "0", "1.5", "1e20"]
    table_rows = list(
        csv.DictReader((bundle / "transport_vs_mean_energy.csv").open(encoding="utf-8"))
    )
    reduced = {
        name: joint_consistency._independent_log_piecewise_cubic_transport(
            table_rows,
            column,
            np.asarray([1.5]),
            support_minimum_eV=1.0,
            support_maximum_eV=2.0,
        )[0]
        for _, name, column, _ in closure_contract._active_transport_functions(
            mapping.closure, source=mapping.bundle.expected_source
        )
    }
    quantity_value = {
        "muN": reduced["muN"] / 1.0e20,
        "DeN": reduced["DeN"] / 1.0e20,
        "muenN": reduced["muenN"] / 1.0e20,
        "DenN": reduced["DenN"] / 1.0e20,
    }
    for property_name, _, tensor, _ in _transport_property_tensors(
        mapping.closure, source=mapping.bundle.expected_source
    ):
        feature_params.append(
            f'<param param="{property_name}" '
            f"value=\"1|1,'{','.join(tensor)}'\"></param>"
        )
    for (
        quantity,
        component,
        actual_expression,
        expected_expression,
        _,
    ) in _transport_audit_components(
        mapping.closure, source=mapping.bundle.expected_source
    ):
        headers.extend(
            [
                f"{actual_expression} (unit) @ t=0",
                f"{expected_expression} (unit) @ t=0",
            ]
        )
        value = quantity_value[quantity] if component in {"rr", "phiphi", "zz"} else 0.0
        row.extend([f"{value:.17e}", f"{value:.17e}"])
        log_argument = closure_arguments.smooth_log_energy_argument(
            "log(ptp.ebar/1[V])", minimum_eV=1.0, maximum_eV=2.0
        )
    model_xml = (
        '<PhysicsFeature tag="pes1">'
        + "".join(feature_params)
        + "<FeatureInfo></FeatureInfo>"
        + "</PhysicsFeature>"
        + '<PhysicsFeature op="Species" tag="Ar_1p">'
        + "<FeatureInfo></FeatureInfo></PhysicsFeature>"
        + '<PhysicsFeature op="SurfaceReaction" tag="sr1">'
        + '<param param="formula" value="1|1,\'Ar+=&gt;Ar\'"></param>'
        + "</PhysicsFeature>"
        + '<param param="ReducedProps" value="1|1,\'1\'"></param>'
        + '<param param="TensorElectronProps" value="1|1,\'0\'"></param>'
        + '<param param="IncludeThermalDiffusion" value="1|1,\''
        + thermal_flag
        + "'\"></param>"
        + '<Expr tag="swClosureVars"><expressions name="sw_logeps" expr="'
        + log_argument
        + '"></expressions></Expr>'
    )
    closure_path = (
        mapping.results.output_directory / "swarm_tables" / "closure_phase.csv"
    )
    radial_closure_path = (
        mapping.results.output_directory / "swarm_tables" / "closure_phase_radial.csv"
    )
    closure_path.parent.mkdir(parents=True, exist_ok=True)
    domain_path = (
        mapping.results.output_directory / "swarm_tables" / "domain_phase_closure.csv"
    )

    def write_domain(mean_energy_eV: float, *, corrupt_transport: bool = False) -> None:
        energy = np.asarray([[mean_energy_eV]], dtype=float)
        transport_reduced = {
            name: float(
                joint_consistency._independent_log_piecewise_cubic_transport(
                    table_rows,
                    column,
                    np.asarray([mean_energy_eV]),
                    support_minimum_eV=1.0,
                    support_maximum_eV=2.0,
                )[0]
            )
            for _, name, column, _ in closure_contract._active_transport_functions(
                mapping.closure, source=mapping.bundle.expected_source
            )
        }
        transport_values: list[float] = []
        transport_headers: list[str] = []
        for quantity, component, expression, _, _ in _transport_audit_components(
            mapping.closure, source=mapping.bundle.expected_source
        ):
            transport_headers.append(f"{expression} (unit) @ t=0")
            value = (
                transport_reduced[quantity] / 1.0e20
                if component in {"rr", "zz"}
                else 0.0
            )
            if corrupt_transport and expression == "ptp.Denzz":
                value *= 2.0
            transport_values.append(value)
        elastic = 6.02214076e23 * float(
            transport_value_audits._independent_rate_coefficient_values(
                mapping,
                "elastic",
                energy,
                [1.0, 2.0],
            )[0, 0]
        )
        excitation = 6.02214076e23 * float(
            transport_value_audits._independent_rate_coefficient_values(
                mapping,
                "excitation",
                energy,
                [1.0, 2.0],
            )[0, 0]
        )
        ionization = 6.02214076e23 * float(
            transport_value_audits._independent_rate_coefficient_values(
                mapping,
                "ionization",
                energy,
                [1.0, 2.0],
            )[0, 0]
        )
        with domain_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                [
                    "% R",
                    "Z",
                    "ptp.Nn (1/m^3) @ t=0",
                    "ptp.Er (V/m) @ t=0",
                    "ptp.Ez (V/m) @ t=0",
                    "ptp.ebar (V) @ t=0",
                    *transport_headers,
                    "ptp.kf_1 (m^3/(mol*s)) @ t=0",
                    "ptp.kf_2 (m^3/(mol*s)) @ t=0",
                    "ptp.kf_3 (m^3/(mol*s)) @ t=0",
                ]
            )
            writer.writerow(
                [
                    "0",
                    "0",
                    "1e20",
                    "1",
                    "0",
                    mean_energy_eV,
                    *transport_values,
                    elastic,
                    excitation,
                    ionization,
                ]
            )

    write_domain(1.5)

    def write_closure() -> None:
        for path in (closure_path, radial_closure_path):
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow([f"% {headers[0]}", *headers[1:]])
                writer.writerow(row)

    write_closure()

    audit = _audit_transport_binding(mapping, model_xml)

    assert audit["passed"] is True, json.dumps(audit, indent=2)
    assert audit["saved_sw_logeps"]["passed"] is True
    assert audit["native_heavy_species_ownership"]["passed"] is True
    assert audit["native_heavy_species_ownership"]["custom_species_lock_count"] == 0
    assert audit["periodic_solver_equation_view_restoration"]["passed"] is True
    assert audit["independent_bundle_value_audit"]["passed"] is True
    assert audit["orientation_regularization"] == {
        "status": "not_applicable",
        "diagnostic_only": True,
    }
    custom_weak = model_xml.replace(
        '<PhysicsFeature op="Species" tag="Ar_1p">'
        "<FeatureInfo></FeatureInfo></PhysicsFeature>",
        '<PhysicsFeature op="Species" tag="Ar_1p"><FeatureInfo>'
        '<lock param="root.comp1.ptp.Ar_1p.weak$1" '
        "value=\"1|1,'0'\"></lock></FeatureInfo></PhysicsFeature>",
    )
    assert (
        _audit_transport_binding(mapping, custom_weak)[
            "native_heavy_species_ownership"
        ]["passed"]
        is False
    )
    saved_periodic_lock = model_xml.replace(
        '<PhysicsFeature tag="pes1">',
        '<PhysicsFeature tag="pes1"><FeatureInfo>'
        '<lock param="ptp.Mn" value="1|1,\'0.04[kg/mol]\'"></lock>'
        "</FeatureInfo>",
    )
    assert (
        _audit_transport_binding(mapping, saved_periodic_lock)[
            "periodic_solver_equation_view_restoration"
        ]["passed"]
        is False
    )
    custom_initialization = model_xml + (
        '<CommonFeature op="Initialization" tag="customIonInit">'
        '<expressions name="WAr_1p_per" expr="0"></expressions>'
        "</CommonFeature>"
    )
    assert (
        _audit_transport_binding(mapping, custom_initialization)[
            "native_heavy_species_ownership"
        ]["passed"]
        is False
    )
    write_domain(1.000001)
    endpoint_audit = _audit_transport_binding(mapping, model_xml)
    assert (
        endpoint_audit["operating_mean_energy_range"]["strictly_inside_common_support"]
        is True
    )
    assert (
        endpoint_audit["operating_mean_energy_range"]["endpoint_margin_passed"] is False
    )
    write_domain(1.5)
    write_domain(1.5, corrupt_transport=True)
    assert _audit_transport_binding(mapping, model_xml)["passed"] is False


def test_transport_audit_requires_one_exact_saved_log_argument(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping = load_gec_ccp_mapping(mapping_path, bundle_path=bundle)
    expression = closure_arguments.smooth_log_energy_argument(
        "log(ptp.ebar/1[V])", minimum_eV=1.0, maximum_eV=2.0
    )
    variable = (
        '<Expr tag="swClosureVars"><expressions name="sw_logeps" expr="'
        + expression
        + '"></expressions></Expr>'
    )

    assert (
        saved_model_audit._saved_sw_logeps_audit(
            variable, minimum_eV=1.0, maximum_eV=2.0
        )["passed"]
        is True
    )
    assert (
        saved_model_audit._saved_sw_logeps_audit(
            variable + variable, minimum_eV=1.0, maximum_eV=2.0
        )["passed"]
        is False
    )
    assert (
        saved_model_audit._saved_sw_logeps_audit(
            variable.replace("log(ptp.ebar/1[V])", "log(ptp.en/ptp.ne)"),
            minimum_eV=1.0,
            maximum_eV=2.0,
        )["passed"]
        is False
    )

    support = joint_consistency._active_closure_mean_energy_support(mapping)
    assert support["passed"] is True
    assert support["common_intersection_mean_energy_eV"] == [1.0, 2.0]
    assert support["accepted_interior_mean_energy_eV"][0] > 1.0
    assert support["accepted_interior_mean_energy_eV"][1] < 2.0


def test_transport_support_intersects_transport_rates_and_function_eedf(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "  reaction_model: external_rates",
            "  reaction_model: function_eedf_preintegrated_inelastic\n"
            "  function_eedf:\n"
            "    table: eedf_f0_comsol_2d.csv\n"
            "    function_tag: sw_eedf_two_term\n"
            "    interpolation: structured_spreadsheet_linear_projection\n"
            "    extrapolation: constant",
        ),
        encoding="utf-8",
    )
    transport_path = bundle / "transport_vs_mean_energy.csv"
    rows = list(csv.reader(transport_path.open(encoding="utf-8")))
    rows[1][0] = "1.1"
    rows[2][0] = "1.9"
    _write_csv(transport_path, rows[0], rows[1:])
    mapping = load_gec_ccp_mapping(mapping_path, bundle_path=bundle)

    support = joint_consistency._active_closure_mean_energy_support(mapping)

    assert support["passed"] is True
    assert support["common_intersection_mean_energy_eV"] == [1.1, 1.9]
    assert support["supports_mean_energy_eV"]["function_eedf"] == [1.0, 2.0]
    assert support["supports_mean_energy_eV"]["rate:excitation"] == [1.0, 2.0]
    assert support["supports_mean_energy_eV"]["rate:ionization"] == [1.0, 2.0]


def test_mc_independent_transport_audit_checks_mixed_field_tensor(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8")
        .replace(
            "expected_source: two_term",
            "expected_source: monte_carlo\n"
            "  expected_mc_transport_estimator_schema: direct_mc_transport.v4",
        )
        .replace(
            "electron_transport: comsol_specify_all_restricted",
            "electron_transport: comsol_specify_all_restricted\n"
            "  zero_field_isotropization_Td: 0.1",
        ),
        encoding="utf-8",
    )
    transport_path = bundle / "transport_vs_mean_energy.csv"
    rows = list(csv.reader(transport_path.open(encoding="utf-8")))
    rows[1][3] = "10"
    rows[2][3] = "11"
    rows[1][7] = "20"
    rows[2][7] = "21"
    _write_csv(transport_path, rows[0], rows[1:])
    mapping = load_gec_ccp_mapping(mapping_path, bundle_path=bundle)
    table_rows = list(csv.DictReader(transport_path.open(encoding="utf-8")))
    reduced = {
        name: joint_consistency._independent_log_piecewise_cubic_transport(
            table_rows,
            column,
            np.asarray([1.5]),
            support_minimum_eV=1.0,
            support_maximum_eV=2.0,
        )[0]
        for _, name, column, _ in closure_contract._active_transport_functions(
            mapping.closure, source=mapping.bundle.expected_source
        )
    }
    neutral = 1.0e20
    er = 3.0e4
    ez = 4.0e4
    zero = 0.0

    def isotropic(value: float) -> list[float]:
        return [value, zero, zero, zero, value, zero, zero, zero, value]

    def field_tensor(longitudinal: float, transverse: float) -> list[float]:
        field_squared = er**2 + ez**2
        denominator = field_squared + (0.1e-21 * neutral) ** 2
        trace_part = field_squared / (3.0 * denominator)
        delta = longitudinal - transverse
        mean = (longitudinal + 2.0 * transverse) / 3.0
        rr = mean + delta * (er**2 / denominator - trace_part)
        pp = mean - delta * trace_part
        zz = mean + delta * (ez**2 / denominator - trace_part)
        rz = delta * er * ez / denominator
        return [rr, zero, rz, zero, pp, zero, rz, zero, zz]

    tensors = {
        "muN": isotropic(reduced["muN"]),
        "DeN": field_tensor(reduced["DeN_L"], reduced["DeN_T"]),
        "muenN": isotropic(reduced["muenN"]),
        "DenN": field_tensor(reduced["DenN_L"], reduced["DenN_T"]),
    }
    radial_headers = [
        "R",
        "Z",
        "ptp.ebar (V) @ t=0",
        "ptp.Nn (1/m^3) @ t=0",
    ]
    radial_row = [0.0, 0.0, 1.5, neutral]
    component_index = {
        name: index for index, name in enumerate(closure_contract.GEC_TENSOR_COMPONENTS)
    }
    for quantity, component, actual, _, _ in _transport_audit_components(
        mapping.closure, source=mapping.bundle.expected_source
    ):
        radial_headers.append(f"{actual} (unit) @ t=0")
        radial_row.append(tensors[quantity][component_index[component]] / neutral)
    domain_headers = [
        "R",
        "Z",
        "ptp.Nn (1/m^3) @ t=0",
        "ptp.Er (V/m) @ t=0",
        "ptp.Ez (V/m) @ t=0",
        "ptp.ebar (V) @ t=0",
    ]
    domain_values = np.asarray([[0.0, 0.0, neutral, er, ez, 1.5]])
    result = transport_value_audits._independent_transport_value_audit(
        mapping,
        radial_headers,
        np.asarray([radial_row]),
        domain_headers,
        domain_values,
        joint_consistency._active_closure_mean_energy_support(mapping),
    )

    assert result["passed"] is True, result
    assert result["mixed_Er_Ez_sample_required"] is True
    assert all(
        item["tensor_components_checked"]
        == list(closure_contract.GEC_AXISYMMETRIC_ACTIVE_TRANSPORT_COMPONENTS)
        for item in result["quantities"].values()
    )
    radial_row[radial_headers.index("ptp.Dezr (unit) @ t=0")] *= 1.5
    rejected = transport_value_audits._independent_transport_value_audit(
        mapping,
        radial_headers,
        np.asarray([radial_row]),
        domain_headers,
        domain_values,
        joint_consistency._active_closure_mean_energy_support(mapping),
    )
    assert rejected["passed"] is False


def test_physical_target_transport_audit_checks_every_domain_phase(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    _enable_monte_carlo_fixture(mapping_path, bundle)
    transport_path = bundle / "transport_vs_mean_energy.csv"
    rows = list(csv.reader(transport_path.open(encoding="utf-8")))
    rows[1][3], rows[2][3] = "10", "11"
    rows[1][6], rows[2][6] = "40", "41"
    rows[1][7], rows[2][7] = "20", "21"
    _write_csv(transport_path, rows[0], rows[1:])
    mapping = load_gec_ccp_mapping(mapping_path, bundle_path=bundle)
    support = joint_consistency._active_closure_mean_energy_support(mapping)

    phase_labels = ["t=0", "t=T/4", "t=T/2"]
    energy = np.asarray([[1.2, 1.5, 1.8], [1.3, 1.6, 1.7]])
    neutral = np.asarray([[1.0e20, 1.1e20, 1.2e20], [2.0e20, 2.1e20, 2.2e20]])
    radial_field = np.asarray([[3.0e4, 0.0, 0.0], [5.0e4, -3.0e4, 2.0e4]])
    axial_field = np.asarray([[4.0e4, 5.0e4, 0.0], [0.0, 4.0e4, -2.0e4]])
    flat_energy = energy.reshape(-1)
    flat_neutral = neutral.reshape(-1)
    flat_er = radial_field.reshape(-1)
    flat_ez = axial_field.reshape(-1)
    table_rows = list(csv.DictReader(transport_path.open(encoding="utf-8")))
    reduced = {
        name: joint_consistency._independent_log_piecewise_cubic_transport(
            table_rows,
            column,
            flat_energy,
            support_minimum_eV=1.0,
            support_maximum_eV=2.0,
        )
        for _, name, column, _ in closure_contract._active_transport_functions(
            mapping.closure, source=mapping.bundle.expected_source
        )
    }

    def isotropic(values: np.ndarray) -> np.ndarray:
        zero = np.zeros_like(values)
        return (
            np.stack(
                [values, zero, zero, zero, values, zero, zero, zero, values],
                axis=1,
            )
            / flat_neutral[:, None]
        )

    def field_aligned(longitudinal: np.ndarray, transverse: np.ndarray) -> np.ndarray:
        field_squared = flat_er**2 + flat_ez**2
        denominator = field_squared + (0.1e-21 * flat_neutral) ** 2
        trace_part = field_squared / (3.0 * denominator)
        delta = longitudinal - transverse
        mean = (longitudinal + 2.0 * transverse) / 3.0
        rr = mean + delta * (flat_er**2 / denominator - trace_part)
        pp = mean - delta * trace_part
        zz = mean + delta * (flat_ez**2 / denominator - trace_part)
        rz = delta * flat_er * flat_ez / denominator
        zero = np.zeros_like(longitudinal)
        return (
            np.stack([rr, zero, rz, zero, pp, zero, rz, zero, zz], axis=1)
            / flat_neutral[:, None]
        )

    tensors = {
        "muN": isotropic(reduced["muN"]),
        "DeN": field_aligned(reduced["DeN_L"], reduced["DeN_T"]),
        "muenN": isotropic(reduced["muenN"]),
        "DenN": field_aligned(reduced["DenN_L"], reduced["DenN_T"]),
    }
    audit_components = _transport_audit_components(
        mapping.closure, source=mapping.bundle.expected_source
    )
    component_index = {
        name: index for index, name in enumerate(closure_contract.GEC_TENSOR_COMPONENTS)
    }
    headers = ["R", "Z"]
    values = [[0.01, 0.02], [0.09, 0.05]]
    for phase, label in enumerate(phase_labels):
        for expression in ("ptp.Nn", "ptp.Er", "ptp.Ez", "ptp.ebar"):
            headers.append(f"{expression} (unit) @ {label}")
        for row in range(2):
            values[row].extend(
                [
                    neutral[row, phase],
                    radial_field[row, phase],
                    axial_field[row, phase],
                    energy[row, phase],
                ]
            )
        for quantity, component, expression, _, _ in audit_components:
            headers.append(f"{expression} (unit) @ {label}")
            for row in range(2):
                sample = row * len(phase_labels) + phase
                values[row].append(
                    tensors[quantity][sample, component_index[component]]
                )
    domain_values = np.asarray(values, dtype=float)

    result = transport_value_audits._independent_domain_transport_value_audit(
        mapping, headers, domain_values, support
    )

    assert result["passed"] is True, result
    assert result["acceptance_scope"] == "full_domain_all_exported_rf_phases"
    assert result["exported_equation_active_component_count"] == 16
    assert result["coverage"] == {
        "spatial_node_count": 2,
        "rf_phase_count": 3,
        "required_sample_count": 6,
        "finite_positive_state_count": 6,
        "within_bundle_support_count": 6,
        "qualified_low_endpoint_continuation_count": 0,
        "above_bundle_support_count": 0,
        "evaluated_sample_count": 6,
        "fraction": 1.0,
        "complete": True,
        "qualified_low_endpoint_continuation": False,
    }
    assert (
        sum(item["required_value_count"] for item in result["quantities"].values())
        == 96
    )
    assert all(
        item["finite_value_count"] == 24
        and item["maximum_normalized_error"] <= 1.0
        and item["p95_normalized_error"] is not None
        and item["worst_sample"] is not None
        for item in result["quantities"].values()
    )

    corrupted = domain_values.copy()
    corrupt_header = "ptp.Derz (unit) @ t=T/2"
    corrupted[0, headers.index(corrupt_header)] = 1.0
    rejected = transport_value_audits._independent_domain_transport_value_audit(
        mapping, headers, corrupted, support
    )
    assert rejected["passed"] is False
    assert rejected["worst_sample"]["R_m"] == pytest.approx(0.01)
    assert rejected["worst_sample"]["Z_m"] == pytest.approx(0.02)
    assert rejected["worst_sample"]["phase"] == "t=T/2"
    assert rejected["worst_sample"]["quantity"] == "DeN"
    assert rejected["worst_sample"]["component"] == "rz"

    incomplete = domain_values.copy()
    incomplete[1, headers.index("ptp.ebar (unit) @ t=T/4")] = np.nan
    incomplete_result = (
        transport_value_audits._independent_domain_transport_value_audit(
            mapping, headers, incomplete, support
        )
    )
    assert incomplete_result["passed"] is False
    assert incomplete_result["coverage"]["evaluated_sample_count"] == 5
    assert incomplete_result["coverage"]["fraction"] == pytest.approx(5.0 / 6.0)
    assert incomplete_result["coverage"]["complete"] is False


def test_mobility_target_transport_audit_checks_only_active_mu(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "electron_transport: comsol_specify_all_restricted",
            "electron_transport: swarm_mobility_einstein",
        ),
        encoding="utf-8",
    )
    mapping = load_gec_ccp_mapping(mapping_path, bundle_path=bundle)
    support = joint_consistency._active_closure_mean_energy_support(mapping)
    table_rows = list(
        csv.DictReader((bundle / "transport_vs_mean_energy.csv").open(encoding="utf-8"))
    )
    reduced_mu = float(
        joint_consistency._independent_log_piecewise_cubic_transport(
            table_rows,
            "reduced_mobility_m2_V_s_m3",
            np.asarray([1.5]),
            support_minimum_eV=1.0,
            support_maximum_eV=2.0,
        )[0]
    )
    local_mu = reduced_mu / 1.0e20
    headers = [
        "R",
        "Z",
        "ptp.Nn (1/m^3) @ t=0",
        "ptp.Er (V/m) @ t=0",
        "ptp.Ez (V/m) @ t=0",
        "ptp.ebar (V) @ t=0",
        "ptp.muerr (m^2/(V*s)) @ t=0",
        "ptp.muezr (m^2/(V*s)) @ t=0",
        "ptp.muerz (m^2/(V*s)) @ t=0",
        "ptp.muezz (m^2/(V*s)) @ t=0",
    ]
    values = np.asarray(
        [[0.01, 0.02, 1.0e20, 3.0e4, 4.0e4, 1.5, local_mu, 0.0, 0.0, local_mu]],
        dtype=float,
    )

    result = transport_value_audits._independent_domain_transport_value_audit(
        mapping, headers, values, support
    )

    assert result["passed"] is True, result
    assert result["expected_equation_active_component_count"] == 4
    assert result["exported_equation_active_component_count"] == 4
    assert result["active_external_quantities"] == ["muN"]
    assert result["active_transport_contract_complete"] is True

    high_values = values.copy()
    high_values[0, headers.index("ptp.ebar (V) @ t=0")] = 2.2
    high_result = transport_value_audits._independent_domain_transport_value_audit(
        mapping,
        headers,
        high_values,
        support,
        qualified_low_endpoint_continuation=True,
    )
    assert high_result["passed"] is False
    assert high_result["coverage"]["above_bundle_support_count"] == 1


def test_function_eedf_rate_audit_excludes_external_elastic_owner() -> None:
    root = Path(__file__).resolve().parents[1]
    mapping = load_gec_ccp_mapping(
        root / "comsol_modes" / "maps" / "argon_gec_ccp_two_term_function_eedf.yaml"
    )

    active = closure_contract._function_eedf_rate_reactions(mapping)

    assert [(index, reaction.process_type) for index, reaction in active] == [
        (2, "excitation"),
        (3, "ionization"),
    ]


def test_function_eedf_rate_range_rejects_only_significant_negative_rate() -> None:
    negligible = function_run_audit._function_eedf_rate_range_quality(
        np.asarray([-1.0e-6, 1.0])
    )
    significant = function_run_audit._function_eedf_rate_range_quality(
        np.asarray([-2.0e-5, 1.0])
    )

    assert negligible["passed"] is True
    assert significant["passed"] is False


def test_gec_ccp_closure_ablation_changes_only_requested_inputs(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    original = mapping_path.read_text(encoding="utf-8")

    mapping_path.write_text(
        original.replace(
            "electron_transport: comsol_specify_all_restricted",
            "electron_transport: comsol",
        ),
        encoding="utf-8",
    )
    rates_plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    rates_source = generate_apply_java(rates_plan.mapping)
    assert '"SpecifyElectronDensityAndEnergy", "SpecifyAll"' not in rates_source
    assert "sw_muN_e" not in rates_source
    assert rates_source.count('"RateConstantForm", "UseRate"') == 3
    assert "remapPeriodicLogEnergySolution" not in rates_source
    rates_manifest = json.loads(rates_plan.plan_json.read_text(encoding="utf-8"))
    assert rates_manifest["closure"]["electron_transport"] == "comsol"
    assert rates_manifest["closure"]["reaction_model"] == "external_rates"
    assert rates_manifest["closure"]["function_eedf"] is None

    mapping_path.write_text(
        original.replace(
            "reaction_model: external_rates",
            "reaction_model: comsol_eedf",
        ),
        encoding="utf-8",
    )
    transport_plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    transport_source = generate_apply_java(transport_plan.mapping)
    assert '"SpecifyElectronDensityAndEnergy", "SpecifyAll"' in transport_source
    assert "sw_log_muN_e" in transport_source
    assert '"RateConstantForm", "UseRate"' not in transport_source
    assert "remapPeriodicLogEnergySolution(model, new double[]{" not in (
        transport_source
    )
    transport_manifest = json.loads(
        transport_plan.plan_json.read_text(encoding="utf-8")
    )
    assert transport_manifest["closure"]["electron_transport"] == (
        "comsol_specify_all_restricted"
    )
    assert transport_manifest["closure"]["reaction_model"] == "comsol_eedf"
    assert transport_manifest["closure"]["function_eedf"] is None
    assert transport_manifest["solve"]["native_physics_initial_values"] is True

    mapping_path.write_text(
        original.replace(
            "electron_transport: comsol_specify_all_restricted",
            "electron_transport: swarm_mobility_einstein",
        ),
        encoding="utf-8",
    )
    mobility_plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    mobility_source = generate_apply_java(mobility_plan.mapping)
    assert '"SpecifyElectronDensityAndEnergy", "SpecifyMueOnly"' in mobility_source
    assert "sw_log_muN_e" in mobility_source
    assert "sw_log_DeN_L_e" not in mobility_source
    assert "sw_log_DeN_T_e" not in mobility_source
    assert "sw_log_muenN_e" not in mobility_source
    assert "sw_log_DenN_L_e" not in mobility_source
    assert "sw_log_DenN_T_e" not in mobility_source
    assert mobility_source.count('"RateConstantForm", "UseRate"') == 3
    mobility_manifest = json.loads(mobility_plan.plan_json.read_text(encoding="utf-8"))
    assert mobility_manifest["closure"]["electron_transport"] == (
        "swarm_mobility_einstein"
    )
    assert mobility_manifest["closure"]["function_eedf"] is None


@pytest.mark.parametrize(
    ("electron_transport", "expected_mode", "all_coefficients"),
    (
        ("comsol_specify_all_restricted", "SpecifyAll", True),
        ("swarm_mobility_einstein", "SpecifyMueOnly", False),
    ),
)
def test_gec_ccp_two_term_function_eedf_uses_structured_spreadsheet(
    tmp_path: Path,
    electron_transport: str,
    expected_mode: str,
    all_coefficients: bool,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    source_text = mapping_path.read_text(encoding="utf-8")
    source_text = source_text.replace(
        "electron_transport: comsol_specify_all_restricted",
        f"electron_transport: {electron_transport}",
    ).replace(
        "reaction_model: external_rates",
        "reaction_model: function_eedf\n"
        "  interpolation: differentiable_log_piecewise_cubic\n"
        "  function_eedf:\n"
        "    table: eedf_f0_comsol_2d.csv\n"
        "    function_tag: sw_eedf_test\n"
        "    interpolation: structured_spreadsheet_linear_projection\n"
        "    extrapolation: constant",
    )
    mapping_path.write_text(source_text, encoding="utf-8")
    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    java = generate_apply_java(plan.mapping)
    run_java = plan.external_run_java.read_text(encoding="utf-8")

    assert f'"SpecifyElectronDensityAndEnergy", "{expected_mode}"' in java
    assert "sw_log_muN_e" in java
    assert ("sw_log_DeN_e" in java) is all_coefficients
    assert "sw_log_DeN_L_e" not in java
    assert "sw_log_DeN_T_e" not in java
    assert ("sw_log_muenN_e" in java) is all_coefficients
    assert ("sw_log_DenN_e" in java) is all_coefficients
    assert "sw_log_DenN_L_e" not in java
    assert "sw_log_DenN_T_e" not in java
    assert 'func().create("sw_eedf_test", "Interpolation")' in java
    assert "sw_eedf_test_shapes" not in java
    assert '"Analytic"' not in java
    assert '.set("struct", "grid")' not in java
    assert '.set("scaledata", "auto")' in java
    assert java.count('.set("nargs", 2)') == 1
    assert "eedf_f0_binding_seed.csv" not in java
    assert '.set("struct", "spreadsheet")' in java
    assert ".discardData();" not in java
    assert (
        '.set("funcnametable", new String[][]{new String[]{'
        '"sw_eedf_test", "1"}})' in java
    )
    assert '.set("argunit", "eV,eV")' in java
    assert '.set("fununit", "1")' in java
    assert '.set("interp", "linear")' in java
    assert '.set("extrap", "const")' in java
    assert ".importData()" in java
    assert '.prop("EEDFSettings").set("eedf", "sw_eedf_test")' in java
    assert java.count('"SpecifyReactionUsing", "UseCrossSectionData"') == 3
    assert java.count('"eedf", "FromPhysicsInterfaceProperty"') == 3
    assert '"SpecifyReactionUsing", "RateConstant"' not in java
    assert '"RateConstantForm", "UseRate"' not in java
    assert "remapPeriodicLogEnergySolution(model," not in java
    assert "promotePeriodicSolution(model," not in java
    assert "SpecifyMeanElectronEnergy" not in java
    assert 'model.param().set("P0", "1[W]")' in run_java
    assert run_java.count('model.study("std3").run()') == 1
    assert "for (double" not in run_java
    assert run_java.count("for (String solutionTag : model.sol().tags())") == 1

    manifest = json.loads(plan.plan_json.read_text(encoding="utf-8"))
    assert manifest["solve"]["coefficient_sweep"] is False
    assert manifest["closure"]["electron_transport"] == electron_transport
    assert manifest["closure"]["reaction_model"] == "function_eedf"
    assert manifest["closure"]["function_eedf"]["representation"] == (
        "structured_spreadsheet_linear_projection"
    )
    assert manifest["closure"]["function_eedf"]["continuity"] == "C0"
    assert manifest["closure"]["function_eedf"]["data_structure"] == ("spreadsheet")
    assert manifest["closure"]["function_eedf"]["grid_shape"] == [2, 201]
    assert manifest["bundle"]["function_eedf"]["mean_energy_range_eV"] == [
        1.0,
        2.0,
    ]
    joint = manifest["closure"]["two_term_joint_consistency"]
    if all_coefficients:
        assert "active_function_eedf" not in joint
        assert "upstream_projected_eedf_evidence" in joint
        assert joint["status"] == "passed"
        assert joint["passed"] is True
        assert joint["relative_error_limit"] == pytest.approx(0.04)
        assert set(joint["coefficients"]) == {
            "muN",
            "DeN_L",
            "DeN_T",
            "muenN",
            "DenN_L",
            "DenN_T",
        }
        assert all(
            item["maximum_relative_error"] < 1.0e-6
            for item in joint["coefficients"].values()
        )
        assert joint["cross_section_source"]["aggregate_floor_m2"] == (
            pytest.approx(1.0e-24)
        )
        kernel = joint["transport_kernel_contract"]
        assert kernel["schema"] == "two_term_active_function_transport.v1"
        assert kernel["process_selection"] == [
            "elastic",
            "excitation",
            "ionization",
        ]
        assert kernel["quadrature"] == (
            "8_point_Gauss_Legendre_per_active_energy_interval"
        )
        assert kernel["transport_interpolation"] == (
            "log_PCHIP_piecewise_cubic_Hermite_C1"
        )
        assert len(kernel["sha256"]) == 64
        assert len(joint["evidence_rows_sha256"]) == 64
        assert joint["thermal_diffusion_enabled"] is False
        assert "two-gradient response matrix" in joint["scope_limit"]
    else:
        assert joint is None
    assert manifest["solve"]["native_physics_initial_values"] is True
    assert plan.native_eedf_audit is not None
    assert plan.native_eedf_audit.java_path.exists()
    assert plan.native_eedf_audit.query_path.exists()
    native = manifest["native_function_eedf_audit"]
    assert native["active_table"]["sha256"] == (plan.native_eedf_audit.table_sha256)
    assert native["query_csv"]["points"] == plan.native_eedf_audit.point_count
    assert "native_eedf_audit" in manifest["generated_java"]
    plan_inputs = validate_gec_plan_inputs(plan)
    assert plan_inputs["active_function_eedf_table_sha256"] == (
        plan.native_eedf_audit.table_sha256
    )


def test_gec_ccp_partial_external_rates_bind_only_selected_process(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    _enable_monte_carlo_fixture(mapping_path, bundle)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8")
        .replace(
            "  reaction_model: external_rates",
            "  reaction_model: external_rates\n  external_rate_processes: [ionization]",
        )
        .replace("role: physical_target", "role: ablation"),
        encoding="utf-8",
    )

    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    java = generate_apply_java(plan.mapping)
    closure = json.loads(plan.plan_json.read_text(encoding="utf-8"))["closure"]
    handling = closure["reaction_handling"]

    assert closure["external_rate_processes"] == ["ionization"]
    assert java.count('"SpecifyReactionUsing", "RateConstant"') == 1
    assert 'feature("eir1").set("SpecifyReactionUsing"' not in java
    assert 'feature("eir2").set("SpecifyReactionUsing"' not in java
    assert 'feature("eir3").set("SpecifyReactionUsing", "RateConstant")' in java
    assert "sw_knorm_ionization_e" in java
    assert "sw_knorm_excitation_e" not in java
    assert [item["binding"] for item in handling] == [
        "UseCrossSectionData",
        "UseCrossSectionData",
        "RateConstant",
    ]
    assert [item["rate_source"] for item in handling] == [
        "comsol_cross_section_integral",
        "comsol_cross_section_integral",
        "monte_carlo_trajectory_time_average_sigma_v",
    ]
    upstream = closure["upstream_eedf_rate_consistency"]
    assert set(upstream["rate_consistency"]["processes"]) == {"ionization"}


def test_gec_ccp_hybrid_function_eedf_binds_only_inelastic_rates(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "reaction_model: external_rates",
            "reaction_model: function_eedf_preintegrated_inelastic\n"
            "  interpolation: differentiable_log_piecewise_cubic\n"
            "  function_eedf:\n"
            "    table: eedf_f0_comsol_2d.csv\n"
            "    function_tag: sw_eedf_test\n"
            "    interpolation: structured_spreadsheet_linear_projection\n"
            "    extrapolation: constant",
        ),
        encoding="utf-8",
    )

    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    java = generate_apply_java(plan.mapping)
    closure_metadata = json.loads(plan.plan_json.read_text(encoding="utf-8"))["closure"]
    metadata = closure_metadata["reaction_handling"]

    assert java.count('"SpecifyReactionUsing", "UseCrossSectionData"') == 1
    assert java.count('"eedf", "FromPhysicsInterfaceProperty"') == 1
    assert java.count('"SpecifyReactionUsing", "RateConstant"') == 2
    assert java.count('"RateConstantForm", "UseRate"') == 2
    assert "sw_logk_elastic_e" not in java
    assert "sw_logk_excitation_e" in java
    assert "sw_logk_ionization_e" in java
    assert '.feature("eir1").set("SpecifyReactionUsing", ' in java
    assert '.feature("eir2").set("SpecifyReactionUsing", "RateConstant")' in java
    assert '.feature("eir3").set("SpecifyReactionUsing", "RateConstant")' in java
    assert '.remove("eir' not in java
    assert '.set("de"' not in java
    assert [item["binding"] for item in metadata] == [
        "UseCrossSectionData",
        "RateConstant",
        "RateConstant",
    ]
    assert all(
        item["energy_loss"] == "preserved_existing_ElectronImpactReaction_de"
        for item in metadata
    )
    assert [item["rate_source"] for item in metadata] == [
        "comsol_cross_section_integral",
        "active_function_eedf_dense_reintegration",
        "active_function_eedf_dense_reintegration",
    ]
    dense = closure_metadata["preintegrated_rate_closure"]
    assert dense["source"] == "active_function_eedf_dense_reintegration"
    assert dense["mean_energy_grid_points"] == 2
    assert dense["strictly_positive_without_floor"] is True
    assert dense["raw_bundle_rates_role"] == (
        "independent_source_and_uncertainty_evidence"
    )
    assert set(dense["processes"]) == {"excitation", "ionization"}
    assert all(item["points"] == 2 for item in dense["processes"].values())
    assert len(dense["closure_sha256"]) == 64
    assert "1.00000000000000000e-20" not in java
    assert "2.00000000000000000e-20" not in java
    blocks = []
    for reaction, binding, de in zip(
        plan.mapping.reactions,
        ("UseCrossSectionData", "RateConstant", "RateConstant"),
        ("0", "11.50", "15.80"),
        strict=True,
    ):
        extra = (
            '<param param="eedf" value="1|1,\'FromPhysicsInterfaceProperty\'"/>'
            if binding == "UseCrossSectionData"
            else '<param param="RateConstantForm" value="1|1,\'UseRate\'"/>'
            f'<param param="kf" value="1|1,\'{contracts.AVOGADRO_PER_MOL}'
            f"*exp(sw_logk_{reaction.process_type}_e(sw_logeps))"
            "*1[m^3/s]'\"/>"
        )
        formula = saved_model_audit.GEC_REACTION_PHYSICS_CONTRACT[
            reaction.process_type
        ]["formula"].replace(">", "&gt;")
        blocks.append(
            '<PhysicsFeature op="ElectronImpactReaction" '
            f'tag="{reaction.feature}">'
            f'<param param="formula" value="1|1,\'{formula}\'"/>'
            f'<param param="de" value="1|1,\'{de}\'"/>'
            '<param param="SpecifyReactionUsing" '
            f"value=\"1|1,'{binding}'\"/>{extra}</PhysicsFeature>"
        )
    model_xml = "".join(blocks)
    audit = saved_model_audit._audit_saved_reaction_handling(plan.mapping, model_xml)
    assert audit["passed"] is True
    assert audit["actual_binding_counts"] == {
        "UseCrossSectionData": 1,
        "RateConstant": 2,
    }
    assert all(item["energy_loss_preserved"] for item in audit["reactions"])
    assert (
        saved_model_audit._audit_saved_reaction_handling(
            plan.mapping,
            model_xml.replace("UseCrossSectionData", "RateConstant", 1),
        )["passed"]
        is False
    )


def test_dense_function_eedf_rate_omits_only_leading_zero_support(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "  reaction_model: external_rates",
            "  reaction_model: function_eedf_preintegrated_inelastic\n"
            "  function_eedf:\n"
            "    table: eedf_f0_comsol_2d.csv\n"
            "    function_tag: sw_eedf_test\n"
            "    interpolation: structured_spreadsheet_linear_projection\n"
            "    extrapolation: constant",
        ),
        encoding="utf-8",
    )
    mapping = load_gec_ccp_mapping(mapping_path, bundle_path=bundle)
    monkeypatch.setattr(
        joint_consistency,
        "_read_function_eedf",
        lambda _path: SimpleNamespace(
            mean_energies_eV=np.asarray([1.0, 2.0, 3.0, 4.0, 5.0])
        ),
    )
    monkeypatch.setattr(
        joint_consistency,
        "_integrate_native_function_eedf_rate",
        lambda _grid, mean, *_args: 0.0 if mean in {1.0, 3.0} else mean * 1.0e-20,
    )

    rows, summary = joint_consistency._build_dense_function_eedf_rate_closure(
        mapping,
        input_mph=mapping.model.input_mph,
    )

    assert len(rows) == 4
    assert {float(row["mean_energy_eV"]) for row in rows} == {4.0, 5.0}
    assert all(float(row["rate_coefficient_m3_s"]) > 0.0 for row in rows)
    assert all(
        process["rows_omitted_before_contiguous_positive_support"] == 3
        and process["zero_rows_observed"] == 2
        and process["positive_rows_discarded_before_support"] == 1
        and process["positive_support_mean_energy_eV"] == [4.0, 5.0]
        for process in summary["processes"].values()
    )
    assert summary["strictly_positive_without_floor"] is True
    assert "accepted_solution" in summary["leading_zero_rate_policy"]


def test_two_term_joint_consistency_fails_closed_on_transport_drift(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "reaction_model: external_rates",
            "reaction_model: function_eedf\n"
            "  interpolation: differentiable_log_piecewise_cubic\n"
            "  function_eedf:\n"
            "    table: eedf_f0_comsol_2d.csv\n"
            "    function_tag: sw_eedf_test\n"
            "    interpolation: structured_spreadsheet_linear_projection\n"
            "    extrapolation: constant",
        ),
        encoding="utf-8",
    )
    transport_path = bundle / "transport_vs_mean_energy.csv"
    rows = list(csv.reader(transport_path.open(encoding="utf-8")))
    mobility_index = rows[0].index("reduced_mobility_m2_V_s_m3")
    for row in rows[1:]:
        row[mobility_index] = str(float(row[mobility_index]) * 1.10)
    _write_csv(transport_path, rows[0], rows[1:])
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["tables"]["transport_vs_mean_energy.csv"]["sha256"] = sha256(
        transport_path.read_bytes()
    ).hexdigest()
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(
        GecCcpWorkflowError,
        match="joint consistency audit failed: muN",
    ):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_two_term_joint_consistency_fails_closed_on_external_rate_drift(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    rates_path = bundle / "rates_vs_mean_energy.csv"
    rows = list(csv.DictReader(rates_path.open(encoding="utf-8")))
    for row in rows:
        if row["process_type"] == "excitation":
            row["rate_coefficient_m3_s"] = str(
                1.2 * float(row["rate_coefficient_m3_s"])
            )
    _write_csv(
        rates_path,
        list(rows[0]),
        [[row[name] for name in rows[0]] for row in rows],
    )
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["tables"]["rates_vs_mean_energy.csv"]["sha256"] = sha256(
        rates_path.read_bytes()
    ).hexdigest()
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(
        GecCcpWorkflowError,
        match="projected-EEDF/rate consistency audit failed: rate_excitation=",
    ):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_two_term_joint_consistency_fails_closed_on_mph_cross_section_drift(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    model_path = mapping_path.parent.parent / "argon_gec_ccp.mph"
    with ZipFile(model_path) as archive:
        model_xml = archive.read("dmodel.xml").decode("utf-8")
    excitation = re.search(
        r'(<PhysicsFeature op="ElectronImpactReaction" tag="eir2">.*?'
        r'<param param="ydata" value=")([^"]+)("/>)',
        model_xml,
        flags=re.DOTALL,
    )
    assert excitation is not None
    changed_values = excitation.group(2).replace("'0'", "'1e-99'", 1)
    model_xml = (
        model_xml[: excitation.start(2)]
        + changed_values
        + model_xml[excitation.end(2) :]
    )
    with ZipFile(model_path, "w") as archive:
        archive.writestr("dmodel.xml", model_xml)

    with pytest.raises(
        GecCcpWorkflowError,
        match="cross_section_identity=excitation",
    ):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_two_term_joint_consistency_formula_rejects_monte_carlo(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8")
        .replace(
            "expected_source: two_term",
            "expected_source: monte_carlo\n"
            "  expected_mc_transport_estimator_schema: direct_mc_transport.v4",
        )
        .replace(
            "electron_transport: comsol_specify_all_restricted",
            "electron_transport: comsol_specify_all_restricted\n"
            "  zero_field_isotropization_Td: 0.1",
        )
        .replace(
            "reaction_model: external_rates",
            "reaction_model: function_eedf\n"
            "  interpolation: differentiable_log_piecewise_cubic\n"
            "  function_eedf:\n"
            "    table: eedf_f0_comsol_2d.csv\n"
            "    function_tag: sw_eedf_test\n"
            "    interpolation: structured_spreadsheet_linear_projection\n"
            "    extrapolation: constant",
        ),
        encoding="utf-8",
    )
    mapping = load_gec_ccp_mapping(mapping_path, bundle_path=bundle)

    with pytest.raises(GecCcpWorkflowError, match="only valid for a two_term"):
        joint_consistency._build_two_term_joint_consistency_audit(
            mapping,
            input_mph=mapping.model.input_mph,
        )


def test_gec_ccp_function_eedf_rejects_obsolete_wide_analytic_contract(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "reaction_model: external_rates",
            "reaction_model: function_eedf\n"
            "  function_eedf:\n"
            "    table: eedf_f0_vs_mean_energy.csv\n"
            "    function_tag: sw_eedf_test\n"
            "    interpolation: dimensionless_c1\n"
            "    extrapolation: constant",
        ),
        encoding="utf-8",
    )

    with pytest.raises(GecCcpWorkflowError, match="canonical artifact"):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_gec_ccp_shape_preserving_log_piecewise_cubic_is_positive_and_exact_jacobian(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "lookup_jacobian: exact",
            "lookup_jacobian: exact\n"
            "  interpolation: differentiable_log_piecewise_cubic",
        ),
        encoding="utf-8",
    )

    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    source = generate_apply_java(plan.mapping)

    assert source.count('.set("interp", "piecewisecubic")') == 7
    assert "sw_log_muN_e" in source
    assert "sw_log_DeN_e" in source
    assert "sw_log_DeN_L_e" not in source
    assert "sw_log_DeN_T_e" not in source
    assert "exp(sw_log_muN_e(" in source
    assert "En_per-Ne_per" in source
    assert "ptp.en/sw_ne_safe" not in source
    assert '"(En_per-Ne_per)"' in source
    assert "sqrt((((En_per-Ne_per))" not in source
    assert "sw_Nn_safe" not in source
    assert "nojac" not in source
    manifest = json.loads(plan.plan_json.read_text(encoding="utf-8"))
    assert manifest["closure"]["coefficient_representation"] == (
        "differentiable_log_piecewise_cubic"
    )


def test_gec_ccp_rejects_obsolete_piecewise_transport_lookup(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "interpolation: differentiable_log_piecewise_cubic",
            "interpolation: piecewise_cubic",
        ),
        encoding="utf-8",
    )

    with pytest.raises(GecCcpWorkflowError, match="piecewise_cubic is obsolete"):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_gec_ccp_rejects_lagged_jacobian_for_log_spline(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "lookup_jacobian: exact",
            "lookup_jacobian: lagged\n"
            "  interpolation: differentiable_log_piecewise_cubic",
        ),
        encoding="utf-8",
    )

    with pytest.raises(
        GecCcpWorkflowError,
        match="lookup_jacobian must be exact",
    ):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_gec_ccp_mobility_einstein_allows_unreported_energy_transport(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "electron_transport: comsol_specify_all_restricted",
            "electron_transport: swarm_mobility_einstein",
        ),
        encoding="utf-8",
    )
    _write_csv(
        bundle / "transport_vs_mean_energy.csv",
        [
            "mean_energy_eV",
            "reduced_mobility_m2_V_s_m3",
            "reduced_diffusion_L_m2_s_m3",
            "reduced_diffusion_T_m2_s_m3",
            "reduced_electron_energy_mobility_m2_V_s_m3",
            "reduced_electron_energy_diffusion_m2_s_m3",
            "reduced_electron_energy_diffusion_L_m2_s_m3",
            "reduced_electron_energy_diffusion_T_m2_s_m3",
        ],
        [
            [1, 10, "", "", "", "", "", ""],
            [2, 11, "", "", "", "", "", ""],
        ],
    )
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["tables"]["transport_vs_mean_energy.csv"]["sha256"] = sha256(
        (bundle / "transport_vs_mean_energy.csv").read_bytes()
    ).hexdigest()
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)

    assert plan.mapping.closure.electron_transport == "swarm_mobility_einstein"


def test_gec_ccp_lagged_lookup_is_not_a_public_mode(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "lookup_jacobian: exact",
            "lookup_jacobian: lagged",
        ),
        encoding="utf-8",
    )

    with pytest.raises(GecCcpWorkflowError, match="must be exact"):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_gec_ccp_rejects_removed_nonlinear_globalization_override(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "  power_W: 1.0\n",
            "  power_W: 1.0\n  external_nonlinear_globalization: magic\n",
        ),
        encoding="utf-8",
    )

    with pytest.raises(GecCcpWorkflowError, match="external_nonlinear_globalization"):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_gec_ccp_rejects_obsolete_initialization_and_solver_paths(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "run:\n",
            "initialization:\n  mode: promote_external_solution\nrun:\n",
        ),
        encoding="utf-8",
    )
    with pytest.raises(GecCcpWorkflowError, match="initialization is obsolete"):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_gec_ccp_rejects_bundle_source_mismatch(tmp_path: Path) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["source"] = "monte_carlo_smoothed"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(GecCcpWorkflowError, match="bundle source mismatch"):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_gec_ccp_rejects_modified_manifest_artifact(tmp_path: Path) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    quality_path = bundle / "quality.csv"
    quality_path.write_text(
        quality_path.read_text(encoding="utf-8") + "\n",
        encoding="utf-8",
    )

    with pytest.raises(GecCcpWorkflowError, match="artifact SHA-256 mismatch"):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_gec_ccp_rejects_unknown_schema_key(tmp_path: Path) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "  expected_source: two_term\n",
            "  expected_source: two_term\n  expected_soruce: typo\n",
        ),
        encoding="utf-8",
    )

    with pytest.raises(GecCcpWorkflowError, match="unknown bundle key"):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_gec_ccp_rejects_obsolete_closure_mode(tmp_path: Path) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "closure:\n",
            "closure:\n  mode: transport_and_eedf\n",
        ),
        encoding="utf-8",
    )
    with pytest.raises(GecCcpWorkflowError, match="closure.mode is obsolete"):
        prepare_gec_ccp_run(mapping_path, bundle_path=bundle)


def test_external_elastic_saved_binding_and_numeric_value_audits(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    _enable_external_elastic_loss(mapping_path, bundle)
    mapping = load_gec_ccp_mapping(mapping_path, bundle_path=bundle)
    qgen = "-ptp.ne*ptp.n_wAr*exp(sw_logKel(sw_logeps))*1[eV*m^3/s]"
    model_xml = (
        """
<Model>
  <PhysicsFeature op="ElectronImpactReaction" tag="eir1">
    <entityFlags>DISABLED</entityFlags>
    <param param="formula" value="1|1,'e+Ar=&gt;e+Ar'"/>
    <param param="SpecifyReactionUsing" value="1|1,'UseCrossSectionData'"/>
    <param param="de" value="1|1,'0'"/>
  </PhysicsFeature>
  <PhysicsFeature op="ElectronImpactReaction" tag="eir2">
    <param param="formula" value="1|1,'e+Ar=&gt;e+Ars'"/>
    <param param="SpecifyReactionUsing" value="1|1,'RateConstant'"/>
    <param param="RateConstantForm" value="1|1,'UseRate'"/>
    <param param="kf" value="1|1,'6.02214076e23[1/mol]*exp(sw_logk_excitation_e(sw_logeps))*1[m^3/s]'"/>
    <param param="de" value="1|1,'11.5'"/>
  </PhysicsFeature>
  <PhysicsFeature op="ElectronImpactReaction" tag="eir3">
    <param param="formula" value="1|1,'e+Ar=&gt;2e+Ar+'"/>
    <param param="SpecifyReactionUsing" value="1|1,'RateConstant'"/>
    <param param="RateConstantForm" value="1|1,'UseRate'"/>
    <param param="kf" value="1|1,'6.02214076e23[1/mol]*exp(sw_logk_ionization_e(sw_logeps))*1[m^3/s]'"/>
    <param param="de" value="1|1,'15.8'"/>
  </PhysicsFeature>
  <PhysicsFeature op="GeneralPowerDeposition" tag="swElLoss">
    <selection selType="GEOMDIM">
      <explicit dim="2" hDim="2" geom="/geom/geom1" entities="2,1"/>
    </selection>
    <param param="Qgen" value="1|1,'"""
        + qgen
        + """'"/>
  </PhysicsFeature>
</Model>
"""
    )
    ownership = saved_model_audit._audit_saved_reaction_handling(mapping, model_xml)
    assert ownership["passed"] is True, ownership
    assert ownership["actual_binding_counts"] == {
        "UseCrossSectionData": 0,
        "RateConstant": 2,
    }
    assert ownership["reactions"][0]["active"] is False
    assert ownership["elastic_energy_loss"]["passed"] is True
    assert (
        saved_model_audit._audit_saved_reaction_handling(
            mapping, model_xml.replace("<entityFlags>DISABLED</entityFlags>", "")
        )["passed"]
        is False
    )

    headers = [
        "R",
        "Z",
        "ptp.ebar (V) @ t=0",
        "ptp.ne (1/m^3) @ t=0",
        "ptp.n_wAr (1/m^3) @ t=0",
        ("-ptp.ne*ptp.n_wAr*exp(sw_logKel(sw_logeps))*1[eV*m^3/s] (W/m^3) @ t=0"),
    ]
    expected = -2.0e15 * 1.0e20 * 1.0e-18 * 1.602176634e-19
    values = np.asarray([[0.0, 0.0, 1.5, 2.0e15, 1.0e20, expected]])
    support = joint_consistency._active_closure_mean_energy_support(mapping)
    numeric = transport_value_audits._independent_elastic_energy_loss_value_audit(
        mapping, headers, values, support
    )
    assert numeric["passed"] is True, numeric
    values[0, -1] *= -1.0
    assert (
        transport_value_audits._independent_elastic_energy_loss_value_audit(
            mapping, headers, values, support
        )["passed"]
        is False
    )


def test_conservation_audit_counts_external_elastic_Qgen(tmp_path: Path) -> None:
    geometry = 0.00260191657233603
    for name, volume in (("swarm_tables", [2, 2, 10, -3, -2, geometry]),):
        directory = tmp_path / name
        directory.mkdir()
        _write_csv(
            directory / "conservation_volume.csv",
            [f"v{index}" for index in range(6)],
            [volume],
        )
        _write_csv(
            directory / "conservation_wall.csv",
            [f"w{index}" for index in range(8)],
            [[-2, -2, -5, 5, 0, 0, 2, 5]],
        )
        _write_csv(
            directory / "conservation_terminal_power.csv",
            ["terminal", "prescribed", "residual", "current"],
            [[10, 10, 0, 0]],
        )
        _write_csv(
            directory / "domain_phase_closure.csv",
            [
                "% R",
                "Z",
                "ptp.ne (1/m^3) @ t=0",
                "ptp.ebar (V) @ t=0",
                "ptp.wAr_1p (1) @ t=0",
            ],
            [[0, 0, 1.0e15, 2.0, 1.0e-8]],
        )
    output = conservation_audit.audit_gec_ccp_conservation_run(
        SimpleNamespace(
            output_directory=tmp_path,
            mapping=SimpleNamespace(
                run=SimpleNamespace(include_builtin_reference=False)
            ),
        )
    )
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["status"] == "passed"
    channels = payload["cases"]["external"]["electron_energy_gross_channels"]
    assert channels["external_elastic_signed_power_W"] == -2
    assert (
        payload["cases"]["external"]["balances"]["electron_energy"]["volume_component"]
        == 5
    )


def test_model_defined_globalization_emits_no_solver_override(
    tmp_path: Path,
) -> None:
    mapping_path, _ = _write_fixture(tmp_path)
    mapping = load_gec_ccp_mapping(mapping_path)

    baseline = generate_run_java(
        mapping,
        run_role="baseline",
        class_name="BaselineRun",
        input_mph=mapping.model.input_mph,
        output_mph=mapping.model.baseline_output_mph,
        study=mapping.model.study,
        solution="sol1",
        time_periodic_feature=mapping.model.time_periodic_feature,
    )
    external = generate_run_java(
        mapping,
        run_role="external",
        class_name="ExternalRun",
        input_mph=mapping.model.output_mph,
        output_mph=mapping.model.output_mph,
        study=mapping.model.external_study,
        solution=mapping.model.external_solution,
        time_periodic_feature=(mapping.model.external_time_periodic_feature),
    )

    assert mapping.run.nonlinear_globalization == "model_defined"
    for source in (baseline, external):
        assert '.feature("fc1").set("dtech"' not in source
        assert '.feature("fc1").set("resscale"' not in source


def test_double_dogleg_is_symmetric_cold_direct_solve_contract(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "run:\n",
            "run:\n  nonlinear_globalization: double_dogleg\n",
            1,
        ),
        encoding="utf-8",
    )

    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    baseline = plan.baseline_run_java.read_text(encoding="utf-8")
    external = plan.external_run_java.read_text(encoding="utf-8")

    for source, solution in ((baseline, "sol1"), (external, "sol3")):
        prefix = f'model.sol("{solution}").feature("s1").feature("fc1")'
        assert source.count(f'{prefix}.set("dtech", "ddog");') == 1
        assert source.count(f'{prefix}.set("resscale", "scalefieldwise");') == 1
        assert source.count('.set("dtech", "ddog");') == 1
        assert source.count('.set("resscale", "scalefieldwise");') == 1

    solve = json.loads(plan.plan_json.read_text(encoding="utf-8"))["solve"]
    assert solve["baseline_nonlinear_method"] == ("COMSOL_double_dogleg_trust_region")
    assert solve["external_nonlinear_method"] == ("COMSOL_double_dogleg_trust_region")
    assert solve["nonlinear_solver_overrides"] is True
    assert solve["nonlinear_globalization"] == {
        "method": "double_dogleg",
        "scope": ["external", "baseline"],
        "residual_scaling": "fieldwise",
        "physical_equations_changed": False,
    }
    assert solve["native_physics_initial_values"] is True
    assert solve["saved_solution_dependency"] is False
    assert solve["power_sweep"] is False
    assert solve["coefficient_sweep"] is False
    assert not any("seed" in str(key).lower() for key in solve)


def test_automatic_newton_no_recovery_is_symmetric_diagnostic_contract(
    tmp_path: Path,
) -> None:
    mapping_path, bundle = _write_fixture(tmp_path)
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "run:\n",
            "run:\n  nonlinear_globalization: automatic_newton_no_recovery\n",
            1,
        ),
        encoding="utf-8",
    )

    plan = prepare_gec_ccp_run(mapping_path, bundle_path=bundle)
    assert plan.mapping.run.nonlinear_globalization == "automatic_newton_no_recovery"
    baseline = plan.baseline_run_java.read_text(encoding="utf-8")
    external = plan.external_run_java.read_text(encoding="utf-8")

    for source, solution in ((baseline, "sol1"), (external, "sol3")):
        prefix = f'model.sol("{solution}").feature("s1").feature("fc1")'
        assert source.count(f'{prefix}.set("useminsteprecovery", "off");') == 1
        assert source.count('set("useminsteprecovery", "off");') == 1
        assert '.feature("fc1").set("dtech"' not in source
        assert '.feature("fc1").set("resscale"' not in source
        assert '.feature("fc1").set("initstep"' not in source
        assert '.feature("fc1").set("initsteph"' not in source
        assert '.feature("fc1").set("minstep",' not in source
        assert '.feature("fc1").set("minsteph",' not in source
        assert '.feature("fc1").set("maxiter"' not in source
        assert "Ar_1p.weak" not in source
        assert '.set("WAr_1p_per",' not in source

    solve = json.loads(plan.plan_json.read_text(encoding="utf-8"))["solve"]
    assert solve["baseline_nonlinear_method"] == (
        "model_defined_automatic_Newton_without_minimum_step_recovery"
    )
    assert solve["external_nonlinear_method"] == (
        "model_defined_automatic_Newton_without_minimum_step_recovery"
    )
    assert solve["nonlinear_solver_overrides"] is True
    assert solve["nonlinear_globalization"] == {
        "method": "automatic_newton_no_recovery",
        "scope": ["external", "baseline"],
        "residual_scaling": "model_defined",
        "physical_equations_changed": False,
        "evaluation_role": "diagnostic",
        "equation_preserving": True,
        "minimum_step_recovery": "off",
        "all_other_nonlinear_settings": "model_defined",
    }
    assert solve["native_physics_initial_values"] is True
    assert solve["saved_solution_dependency"] is False
    assert solve["power_sweep"] is False
    assert solve["coefficient_sweep"] is False


def _write_fixture(root: Path) -> tuple[Path, Path]:
    maps = root / "maps"
    maps.mkdir()
    model = root / "argon_gec_ccp.mph"
    with (
        Path(__file__).resolve().parents[1]
        / "examples"
        / "cross_sections"
        / "argon_application_library.csv"
    ).open(encoding="utf-8", newline="") as handle:
        application_rows = list(csv.DictReader(handle))
    fixture_cross_sections = {
        process_type: (
            np.asarray(
                [
                    float(row["energy_eV"])
                    for row in application_rows
                    if row["type"] == process_type
                ]
            ),
            np.asarray(
                [
                    float(row["cross_section_m2"])
                    for row in application_rows
                    if row["type"] == process_type
                ]
            ),
        )
        for process_type in ("elastic", "excitation", "ionization")
    }

    def xml_array(values: np.ndarray) -> str:
        encoded = ",".join(f"'{float(value):.17g}'" for value in values)
        return f"1|{len(values)},{encoded}"

    elastic_energy, elastic_sigma = fixture_cross_sections["elastic"]
    excitation_energy, excitation_sigma = fixture_cross_sections["excitation"]
    ionization_energy, ionization_sigma = fixture_cross_sections["ionization"]
    xml = f"""
<Model>
  <GeomSequence tag="geom1"><axisymmetric>true</axisymmetric></GeomSequence>
  <Physics op="ColdPlasmaTimePeriodic" tag="ptp"/>
  <PhysicsFeature op="ElectronImpactReaction" tag="eir1">
    <param param="eedf" value="1|1,'FromPhysicsInterfaceProperty'"/>
    <param param="xdata" value="{xml_array(elastic_energy)}"/>
    <param param="ydata" value="{xml_array(elastic_sigma)}"/>
  </PhysicsFeature>
  <PhysicsFeature op="ElectronImpactReaction" tag="eir2">
    <param param="xdata" value="{xml_array(excitation_energy)}"/>
    <param param="ydata" value="{xml_array(excitation_sigma)}"/>
  </PhysicsFeature>
  <PhysicsFeature op="ElectronImpactReaction" tag="eir3">
    <param param="xdata" value="{xml_array(ionization_energy)}"/>
    <param param="ydata" value="{xml_array(ionization_sigma)}"/>
  </PhysicsFeature>
  <PhysicsFeature op="Species" tag="Ar">
    <param param="sType" value="1|1,'neutral'"/>
    <param param="M" value="1|1,'0.04[kg/mol]'"/>
    <param param="FromMassConstraint" value="1|1,'1'"/>
    <param param="DT" value="1|1,'0'"/>
    <param param="z" value="1|1,'0'"/>
  </PhysicsFeature>
  <PhysicsFeature op="Species" tag="Ars">
    <entityFlags>DISABLED</entityFlags>
    <param param="sType" value="1|1,'neutral'"/>
    <param param="M" value="1|1,'0.04[kg/mol]'"/>
    <param param="FromMassConstraint" value="1|1,'0'"/>
    <param param="DT" value="1|1,'0'"/>
    <param param="z" value="1|1,'0'"/>
  </PhysicsFeature>
  <PhysicsFeature op="Species" tag="Ar_1p">
    <param param="sType" value="1|1,'ion'"/>
    <param param="M" value="1|1,'0.04[kg/mol]'"/>
    <param param="FromMassConstraint" value="1|1,'0'"/>
    <param param="DT" value="1|1,'0'"/>
    <param param="z" value="1|1,'1'"/>
  </PhysicsFeature>
  <PhysicsFeature op="SurfaceReaction" tag="sr1">
    <param param="formula" value="1|1,'Ar+=&gt;Ar'"/>
  </PhysicsFeature>
  <PhysicsFeature op="PlasmaEsModel" tag="pes1"/>
  <PhysicsProp>
    <param param="eedf" value="1|1,'Druyvesteyn'"/>
    <param param="HeavySpeciesSelection" value="1|1,'BaseGeometry'"/>
    <param param="Formulation" value="1|1,'FEMLogLinear'"/>
    <param param="DiffusionModel" value="1|1,'MixtureAveraged'"/>
    <param param="Migration" value="1|1,'1'"/>
    <param param="Convection" value="1|1,'0'"/>
    <param param="MixtureDiffusionCorrection" value="1|1,'0'"/>
    <param param="IonTensorProps" value="1|1,'0'"/>
    <param param="ElectricFieldAppliedToIons" value="1|1,'Instantaneous'"/>
  </PhysicsProp>
  <Study tag="std1"/>
  <Study tag="std3"/>
  <Study tag="std2"/>
  <DatasetFeature tag="dset1"/>
  <DatasetFeature tag="dset2"/>
  <DatasetFeature tag="dset3"/>
  <DatasetFeature tag="dset4"/>
  <DatasetFeature tag="dset5"/>
  <DatasetFeature tag="cln1"/>
  <DatasetFeature tag="cln2"/>
</Model>
""".strip()
    with ZipFile(model, "w") as archive:
        archive.writestr("dmodel.xml", xml)

    bundle = root / "bundle"
    bundle.mkdir()
    transport_columns = [
        "mean_energy_eV",
        "E_over_N_Td",
        "reduced_mobility_m2_V_s_m3",
        "reduced_diffusion_L_m2_s_m3",
        "reduced_diffusion_T_m2_s_m3",
        "reduced_electron_energy_mobility_m2_V_s_m3",
        "reduced_electron_energy_diffusion_m2_s_m3",
        "reduced_electron_energy_diffusion_L_m2_s_m3",
        "reduced_electron_energy_diffusion_T_m2_s_m3",
    ]
    _write_csv(
        bundle / "transport_vs_mean_energy.csv",
        transport_columns,
        [
            [1, 10, 10, 20, 20, 30, 40, 40, 40],
            [2, 20, 11, 21, 21, 31, 41, 41, 41],
        ],
    )
    rate_columns = [
        "mean_energy_eV",
        "process_type",
        "rate_coefficient_m3_s",
    ]
    _write_csv(
        bundle / "rates_vs_mean_energy.csv",
        rate_columns,
        [
            [1, "elastic", 1e-14],
            [2, "elastic", 2e-14],
            [1, "excitation", 1e-20],
            [2, "excitation", 2e-20],
            [1, "ionization", 1e-22],
            [2, "ionization", 2e-22],
        ],
    )
    rate_evidence_columns = [
        "E_over_N_Td",
        "mean_energy_eV",
        "species",
        "process",
        "process_type",
        "rate_coefficient_mean_m3_s",
        "rate_coefficient_ci95_low_m3_s",
        "rate_coefficient_ci95_high_m3_s",
        "estimate_status",
        "uncertainty_available",
        "pooled_event_count",
        "pooled_target_exposure_s_m3",
        "pooled_all_zero_upper_95_m3_s",
        "pooled_zero_event_status",
    ]
    _write_csv(
        bundle / "rate_evidence.csv",
        rate_evidence_columns,
        [
            [
                10,
                1,
                "Ar",
                "elastic",
                "elastic",
                1e-14,
                0.9e-14,
                1.1e-14,
                "estimate",
                True,
                10,
                1e20,
                0.0,
                "observed_events",
            ],
            [
                20,
                2,
                "Ar",
                "elastic",
                "elastic",
                2e-14,
                1.8e-14,
                2.2e-14,
                "estimate",
                True,
                10,
                1e20,
                0.0,
                "observed_events",
            ],
            [
                10,
                1,
                "Ar",
                "excitation",
                "excitation",
                1e-20,
                0.9e-20,
                1.1e-20,
                "estimate",
                True,
                10,
                1e20,
                0.0,
                "observed_events",
            ],
            [
                20,
                2,
                "Ar",
                "excitation",
                "excitation",
                2e-20,
                1.8e-20,
                2.2e-20,
                "estimate",
                True,
                10,
                1e20,
                0.0,
                "observed_events",
            ],
            [
                10,
                1,
                "Ar",
                "ionization",
                "ionization",
                1e-22,
                0.9e-22,
                1.1e-22,
                "estimate",
                True,
                10,
                1e20,
                0.0,
                "observed_events",
            ],
            [
                20,
                2,
                "Ar",
                "ionization",
                "ionization",
                2e-22,
                1.8e-22,
                2.2e-22,
                "estimate",
                True,
                10,
                1e20,
                0.0,
                "observed_events",
            ],
        ],
    )
    function_eedf_columns = [
        "dimensionless_energy",
        "shape_0000",
        "shape_0001",
    ]
    _write_csv(
        bundle / "eedf_f0_vs_mean_energy.csv",
        function_eedf_columns,
        [
            [0.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
            [2.0, 0.0, 0.0],
        ],
    )
    active_eedf_columns = [
        "electron_energy_eV",
        "mean_energy_eV",
        "eepf_eV_m32",
    ]
    active_energy = np.linspace(0.0, 100.0, 201)
    active_means = np.asarray([1.0, 2.0])

    def discrete_exponential(requested_mean: float) -> np.ndarray:
        lower_decay = 1.0e-3
        upper_decay = 100.0
        for _ in range(100):
            decay = 0.5 * (lower_decay + upper_decay)
            trial = np.exp(-decay * active_energy)
            normalization, first_moment = piecewise_linear_weighted_moments(
                active_energy, trial
            )
            if first_moment / normalization > requested_mean:
                lower_decay = decay
            else:
                upper_decay = decay
        result = np.exp(-0.5 * (lower_decay + upper_decay) * active_energy)
        normalization, _ = piecewise_linear_weighted_moments(active_energy, result)
        return result / normalization

    active_values_array = np.asarray(
        [discrete_exponential(float(mean)) for mean in active_means]
    )
    rate_nodes, rate_weights = np.polynomial.legendre.leggauss(8)
    rate_gamma = np.sqrt(2.0 * 1.602176634e-19 / 9.1093837015e-31)

    def projected_rate(
        f0: np.ndarray,
        cross_energy: np.ndarray,
        cross_sigma: np.ndarray,
    ) -> float:
        knots = np.unique(
            np.concatenate(
                (
                    active_energy,
                    cross_energy[
                        (cross_energy > active_energy[0])
                        & (cross_energy < active_energy[-1])
                    ],
                )
            )
        )
        rate_left = knots[:-1, None]
        rate_right = knots[1:, None]
        rate_half_width = 0.5 * (rate_right - rate_left)
        rate_energy = 0.5 * (rate_right + rate_left) + rate_half_width * rate_nodes
        sampled_f0 = np.interp(rate_energy, active_energy, f0)
        sampled_sigma = np.interp(
            rate_energy,
            cross_energy,
            cross_sigma,
            left=0.0,
            right=float(cross_sigma[-1]),
        )
        return float(
            np.sum(
                rate_half_width
                * rate_weights
                * rate_gamma
                * sampled_sigma
                * rate_energy
                * sampled_f0
            )
        )

    consistent_rate_rows: list[list[object]] = []
    consistent_rate_evidence_rows: list[list[object]] = []
    for process_type in ("elastic", "excitation", "ionization"):
        cross_energy, cross_sigma = fixture_cross_sections[process_type]
        for mean, field, f0 in zip(
            active_means,
            (10.0, 20.0),
            active_values_array,
            strict=True,
        ):
            rate = projected_rate(f0, cross_energy, cross_sigma)
            consistent_rate_rows.append([mean, process_type, rate])
            consistent_rate_evidence_rows.append(
                [
                    field,
                    mean,
                    "Ar",
                    process_type,
                    process_type,
                    rate,
                    0.9 * rate,
                    1.1 * rate,
                    "estimate",
                    True,
                    10,
                    1e20,
                    0.0,
                    "observed_events",
                ]
            )
    _write_csv(
        bundle / "rates_vs_mean_energy.csv",
        rate_columns,
        consistent_rate_rows,
    )
    _write_csv(
        bundle / "rate_evidence.csv",
        rate_evidence_columns,
        consistent_rate_evidence_rows,
    )
    # Keep the fixture's projected EEDF, rates, and full two-term transport
    # physically co-consistent with the application-library cross sections
    # embedded above. The monotone exponential tail resolves both inelastic
    # channels while keeping the requested 1 and 2 eV moments exact.
    nodes, weights = np.polynomial.legendre.leggauss(8)
    left = active_energy[:-1, None]
    right = active_energy[1:, None]
    half_width = 0.5 * (right - left)
    sample_energy = 0.5 * (right + left) + half_width * nodes
    fraction = (sample_energy - left) / (right - left)
    gamma = np.sqrt(2.0 * 1.602176634e-19 / 9.1093837015e-31)
    sigma_m = np.zeros_like(sample_energy)
    for cross_energy, cross_sigma in fixture_cross_sections.values():
        sigma_m += np.interp(sample_energy, cross_energy, cross_sigma)
    sigma_m = np.maximum(
        sigma_m,
        table_contracts.TWO_TERM_MOMENTUM_CROSS_SECTION_FLOOR_M2,
    )
    consistent_transport_rows = []
    for mean, field, f0 in zip(
        active_means, (10.0, 20.0), active_values_array, strict=True
    ):
        f0_left = f0[:-1, None]
        f0_right = f0[1:, None]
        f0_sample = f0_left + (f0_right - f0_left) * fraction
        derivative = (f0_right - f0_left) / (right - left)

        def integrate(values: np.ndarray) -> float:
            return float(np.sum(half_width * weights * values))

        mobility = -gamma / 3.0 * integrate(sample_energy / sigma_m * derivative)
        diffusion = gamma / 3.0 * integrate(sample_energy / sigma_m * f0_sample)
        energy_mobility = (
            -gamma
            / (3.0 * mean)
            * integrate(sample_energy * sample_energy / sigma_m * derivative)
        )
        energy_diffusion = (
            gamma
            / (3.0 * mean)
            * integrate(sample_energy * sample_energy / sigma_m * f0_sample)
        )
        consistent_transport_rows.append(
            [
                mean,
                field,
                mobility,
                diffusion,
                diffusion,
                energy_mobility,
                energy_diffusion,
                energy_diffusion,
                energy_diffusion,
            ]
        )
    _write_csv(
        bundle / "transport_vs_mean_energy.csv",
        transport_columns,
        consistent_transport_rows,
    )
    projected_norm_error = 0.0
    projected_mean_error = 0.0
    for requested_mean, row_values in zip(
        active_means, active_values_array, strict=True
    ):
        normalization, first_moment = piecewise_linear_weighted_moments(
            active_energy, row_values
        )
        projected_norm_error = max(projected_norm_error, abs(normalization - 1.0))
        projected_mean_error = max(
            projected_mean_error,
            abs(first_moment / normalization - requested_mean) / requested_mean,
        )
    _write_csv(
        bundle / "eedf_f0_comsol_2d.csv",
        active_eedf_columns,
        [
            [energy, mean, active_values_array[mean_index, energy_index]]
            for mean_index, mean in enumerate(active_means)
            for energy_index, energy in enumerate(active_energy)
        ],
    )
    quality_columns = [
        "E_over_N_Td",
        "passed",
        "eedf_normalization_error",
    ]
    _write_csv(
        bundle / "quality.csv",
        quality_columns,
        [[10, 1, 0], [20, 1, 1e-16]],
    )
    (bundle / "manifest.json").write_text(
        json.dumps(
            {
                "status": "ok",
                "source": "two_term",
                "hashes": {
                    "workflow_config_sha256": "c" * 64,
                    "base_config_sha256": "a" * 64,
                    "cross_sections_sha256": "b" * 64,
                },
                "quality_thresholds": quality_thresholds_payload(QualityThresholds()),
                "source_policy": {
                    "field_type": "dc",
                    "transport_definition": "flux",
                    "rf_frequency_Hz": None,
                    "postprocess": "none",
                },
                "valid_ranges": {"mean_energy_eV": [1, 2]},
                "monotonicity": {"mean_energy_strictly_monotonic": True},
                "tables": {
                    "transport_vs_mean_energy.csv": {
                        "columns": transport_columns,
                        "artifact_role": "canonical_comsol_coefficient_input",
                        "sha256": sha256(
                            (bundle / "transport_vs_mean_energy.csv").read_bytes()
                        ).hexdigest(),
                    },
                    "rates_vs_mean_energy.csv": {
                        "columns": rate_columns,
                        "artifact_role": "canonical_comsol_coefficient_input",
                        "sha256": sha256(
                            (bundle / "rates_vs_mean_energy.csv").read_bytes()
                        ).hexdigest(),
                    },
                    "rate_evidence.csv": {
                        "columns": rate_evidence_columns,
                        "artifact_role": "raw_swarm_rate_evidence",
                        "statistics": {
                            "replicate_interval": "two_sided_student_t_95",
                            "zero_event_confidence": 0.95,
                            "zero_event_upper_bound": (
                                "poisson_zero_count_over_pooled_target_exposure"
                            ),
                        },
                        "sha256": sha256(
                            (bundle / "rate_evidence.csv").read_bytes()
                        ).hexdigest(),
                    },
                    "eedf_f0_vs_mean_energy.csv": {
                        "columns": function_eedf_columns,
                        "artifact_role": ("canonicalized_function_eedf_evidence"),
                        "canonical_comsol_input": False,
                        "representation": ("dimensionless_shape_preserving_c1_convex"),
                        "anchor_mean_energy_eV": [1.0, 2.0],
                        "shape_columns": ["shape_0000", "shape_0001"],
                        "normalization_error_max": 1.0e-12,
                        "mean_energy_relative_error_max": 1.0e-12,
                        "nonnegative_minimum": 0.0,
                        "mean_axis_derivative_jump_max": 0.0,
                        "sha256": sha256(
                            (bundle / "eedf_f0_vs_mean_energy.csv").read_bytes()
                        ).hexdigest(),
                    },
                    "eedf_f0_comsol_2d.csv": {
                        "argument": "electron_energy_eV,mean_energy_eV",
                        "argument_order": [
                            "electron_energy_eV",
                            "mean_energy_eV",
                        ],
                        "artifact_role": ("canonical_comsol_function_eedf_input"),
                        "canonical_comsol_input": True,
                        "columns": active_eedf_columns,
                        "format": "csv",
                        "representation": (
                            "physical_2d_adaptive_moment_rate_projected"
                        ),
                        "grid_axis_order": [
                            "electron_energy_eV",
                            "mean_energy_eV",
                        ],
                        "row_order": ("mean_energy_major_then_electron_energy"),
                        "energy_grid_points": len(active_energy),
                        "mean_energy_grid_points": len(active_means),
                        "grid_shape": [len(active_means), len(active_energy)],
                        "electron_energy_range_eV": [0.0, 100.0],
                        "mean_energy_range_eV": [1.0, 2.0],
                        "projected_normalization_error_max": (projected_norm_error),
                        "projected_mean_energy_relative_error_max": (
                            projected_mean_error
                        ),
                        "projected_nonnegative_minimum": float(
                            np.min(active_values_array)
                        ),
                        "comsol_import": {
                            "source": "file",
                            "struct": "spreadsheet",
                            "nargs": 2,
                            "argunit": "eV,eV",
                            "fununit": "1",
                            "interp": "linear",
                            "extrap": "const",
                            "funcnametable_position": "1",
                            "scaledata": "auto",
                        },
                        "sha256": sha256(
                            (bundle / "eedf_f0_comsol_2d.csv").read_bytes()
                        ).hexdigest(),
                    },
                    "quality.csv": {
                        "columns": quality_columns,
                        "artifact_role": "swarm_quality_audit",
                        "sha256": sha256(
                            (bundle / "quality.csv").read_bytes()
                        ).hexdigest(),
                    },
                },
            }
        ),
        encoding="utf-8",
    )
    cross_sections = root / "argon_cross_sections.csv"
    _write_csv(
        cross_sections,
        [
            "species",
            "process",
            "type",
            "threshold_eV",
            "mass_amu",
            "energy_eV",
            "cross_section_m2",
        ],
        [
            ["Ar", "e+Ar=>2e+Ar+", "ionization", 15.8, 39.948, 0, 0],
            ["Ar", "e+Ar=>2e+Ar+", "ionization", 15.8, 39.948, 15.8, 0],
            ["Ar", "e+Ar=>2e+Ar+", "ionization", 15.8, 39.948, 20, 1e-20],
            ["Ar", "e+Ar=>2e+Ar+", "ionization", 15.8, 39.948, 100, 1e-20],
        ],
    )
    mapping = maps / "gec.yaml"
    mapping.write_text(
        """
schema_version: 2
model:
  input_mph: ../argon_gec_ccp.mph
  baseline_output_mph: ../work/baseline.mph
  external_output_mph: ../work/external.mph
  component: comp1
  physics: ptp
  plasma_feature: pes1
  time_periodic_study: std1
  time_periodic_feature: tper
  external_time_periodic_study: std3
  external_time_periodic_feature: tper1
  external_solution: sol3
  conversion_study: std2
  expected_original_eedf: Druyvesteyn
bundle:
  path: ../bundle
  expected_source: two_term
  expected_field_type: dc
  expected_transport_definition: flux
reactions:
  - {name: elastic, feature: eir1, process_type: elastic}
  - {name: excitation, feature: eir2, process_type: excitation}
  - {name: ionization, feature: eir3, process_type: ionization}
closure:
  lookup_jacobian: exact
  interpolation: differentiable_log_piecewise_cubic
  electron_transport: comsol_specify_all_restricted
  gradient_response_policy: standard_local_energy
  thermal_diffusion_model: off_restricted_diagonal
  reaction_model: external_rates
run:
  power_parameter: P0
  power_W: 1.0
  include_builtin_reference: true
  source_stabilization: false
  reaction_source_stabilization: false
  axis_dataset: cln1
  radial_dataset: cln2
  period_dataset: dset1
  external_period_dataset: dset4
  phase_dataset: dset3
  baseline_waveform_dataset: dset2
  external_waveform_dataset: dset5
results:
  role: physical_target
  output_directory: ../results
logs:
  path: ../logs
""".lstrip(),
        encoding="utf-8",
    )
    return mapping, bundle


def _write_csv(path: Path, columns: list[str], rows: list[list[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        writer.writerows(rows)


def _enable_monte_carlo_fixture(mapping_path: Path, bundle: Path) -> None:
    """Give the compact fixture the canonical direct-MC provenance contract."""

    transport_definition = "mc_flux_particle_tracking_fixed_population"
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8")
        .replace(
            "expected_source: two_term",
            "expected_source: monte_carlo\n"
            "  expected_mc_transport_estimator_schema: direct_mc_transport.v4",
        )
        .replace(
            "expected_transport_definition: flux",
            f"expected_transport_definition: {transport_definition}",
        )
        .replace(
            "electron_transport: comsol_specify_all_restricted",
            "electron_transport: comsol_specify_all_restricted\n"
            "  zero_field_isotropization_Td: 0.1",
        ),
        encoding="utf-8",
    )

    quality_columns = [
        "E_over_N_Td",
        "passed",
        "eedf_normalization_error",
        *mc_bundle_validation.REQUIRED_MC_QUALITY_COLUMNS,
    ]
    quality_values = {
        "aggregate_quality_passed": 1,
        "failure_reasons_json": "[]",
        "aggregate_failure_reasons_json": "[]",
        "mobility_rse": 0.01,
        "diffusion_L_rse": 0.01,
        "diffusion_T_rse": 0.01,
        "energy_mobility_rse": 0.01,
        "energy_diffusion_L_rse": 0.01,
        "energy_diffusion_T_rse": 0.01,
        "max_major_rate_rse": 0.01,
        "valid_replicates": 3,
        "uncertainty_available": 1,
        "solver_diagnostics_available": 1,
        "solver_transport_qualified": 1,
        "solver_transport_replicates": 3,
        "solver_origin_stationarity_limiting_field": "mean_energy_eV",
        "solver_origin_stationarity_limiting_relative_ci95_bound": 0.05,
        "solver_mean_energy_stationarity_relative_ci95_bound": 0.05,
        "solver_mobility_stationarity_absolute_log_drift": 0.02,
        "solver_lag_convergence_max_relative_ci95_bound": 0.05,
        "solver_transport_mean_energy_max_relative_ci95_bound": 0.05,
        "solver_population_growth_max_relative_ci95_bound": "",
        "solver_population_growth_gate_mode": "",
        "solver_population_growth_poisson_interval_max_ratio": "",
        "solver_population_growth_sparse_max_metric": "",
        "solver_origin_stationarity_limiting_relative_tolerance": 0.10,
        "solver_lag_convergence_relative_tolerance": 0.25,
        "solver_transport_mean_energy_relative_tolerance": 0.10,
        "solver_population_growth_relative_tolerance": "",
        "quality_source": "monte_carlo_aggregate_quality",
    }
    quality_path = bundle / "quality.csv"
    _write_csv(
        quality_path,
        quality_columns,
        [
            [
                field,
                1,
                0.0,
                *(
                    quality_values[name]
                    for name in mc_bundle_validation.REQUIRED_MC_QUALITY_COLUMNS
                ),
            ]
            for field in (10.0, 20.0)
        ],
    )
    sampling_plan = [
        {
            "e_over_n_Td": field,
            "particles": 64,
            "warmup_collisions": 100,
            "max_collisions": 200,
            "tail_max_collisions": 200,
            "replicas": 3,
            "transport_correlation_lag_barriers": 64,
            "transport_estimator": "single_field",
        }
        for field in (10.0, 20.0)
    ]
    sampling_json = json.dumps(
        sampling_plan,
        sort_keys=True,
        separators=(",", ":"),
    )

    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["source"] = "monte_carlo"
    manifest["source_policy"].update(
        {
            "transport_definition": transport_definition,
            "transport_eligible_anchor_points": 2,
            "transport_eligible_E_over_N_Td": [10.0, 20.0],
            "transport_excluded_E_over_N_Td": {},
        }
    )
    manifest["valid_ranges"]["E_over_N_Td"] = [10.0, 20.0]
    manifest["tables"]["quality.csv"].update(
        {
            "columns": quality_columns,
            "sha256": sha256(quality_path.read_bytes()).hexdigest(),
        }
    )
    manifest["hashes"].update(
        {
            "mc_sampling_plan_json": sampling_json,
            "mc_sampling_plan_sha256": sha256(
                sampling_json.encode("utf-8")
            ).hexdigest(),
            "mc_transport_estimator_schema_version": (
                mc_evidence.DIRECT_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
            ),
            "mc_eedf_estimator_schema_version": (
                mc_evidence.MC_EEDF_ESTIMATOR_SCHEMA_VERSION
            ),
            "mc_solver_source_sha256": monte_carlo_source_sha256(),
        }
    )
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")


def _enable_two_term_temporal_growth_fixture(
    mapping_path: Path,
    bundle: Path,
    *,
    growth_frequency_s_inv: float,
) -> None:
    mapping_path.write_text(
        mapping_path.read_text(encoding="utf-8").replace(
            "  expected_transport_definition: flux",
            "  expected_transport_definition: flux\n"
            "  expected_two_term_transport_kernel_schema: "
            "two_term_temporal_growth_transport.v1",
        ),
        encoding="utf-8",
    )
    table = bundle / "transport_vs_mean_energy.csv"
    with table.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        columns = list(reader.fieldnames or ())
    evidence_columns = list(table_contracts.TWO_TERM_TEMPORAL_GROWTH_TRANSPORT_COLUMNS)
    density = 1.0e22
    for row in rows:
        row.update(
            {
                "gas_number_density_m3": density,
                "temporal_growth_frequency_s_inv": growth_frequency_s_inv,
                "reduced_temporal_growth_frequency_m3_s": (
                    growth_frequency_s_inv / density
                ),
            }
        )
    _write_csv(
        table,
        [*columns, *evidence_columns],
        [
            [row.get(column, "") for column in [*columns, *evidence_columns]]
            for row in rows
        ],
    )
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["source_policy"]["two_term_transport_kernel"] = {
        "schema": "two_term_temporal_growth_transport.v1",
        "correction_applied_every_anchor": True,
        "effective_momentum_model": (
            "nu_m_eff(epsilon)=nu_m(epsilon)+growth_frequency"
        ),
        "growth_frequency_source": "converged_temporal_growth_eigenvalue",
        "momentum_cross_section_floor_policy": (
            "sum_raw_process_cross_sections_then_single_floor"
        ),
        "minimum_momentum_cross_section_m2": 1.0e-24,
        "evidence_columns": evidence_columns,
        "anchor_points": len(rows),
        "gas_number_density_m3": [density, density],
        "growth_frequency_s_inv": [
            growth_frequency_s_inv,
            growth_frequency_s_inv,
        ],
    }
    manifest["tables"]["transport_vs_mean_energy.csv"].update(
        {
            "columns": [*columns, *evidence_columns],
            "sha256": sha256(table.read_bytes()).hexdigest(),
        }
    )
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")


def _enable_external_elastic_loss(
    mapping_path: Path,
    bundle: Path,
    *,
    source: str = "two_term",
) -> None:
    columns = [
        "mean_energy_eV",
        "E_over_N_Td",
        "E_over_N_V_m2",
        "elastic_energy_loss_rate_coefficient_eV_m3_s",
        "elastic_energy_loss_standard_error_eV_m3_s",
        "elastic_energy_loss_relative_standard_error",
        "elastic_energy_loss_ci95_low_eV_m3_s",
        "elastic_energy_loss_ci95_high_eV_m3_s",
        "ci95_critical_value",
        "estimate_status",
        "valid_replicates",
        "uncertainty_available",
    ]
    if source == "monte_carlo":
        critical = 4.302652729911275
        rows = [
            [
                mean,
                field,
                field * 1.0e-21,
                coefficient,
                0.05 * coefficient,
                0.05,
                coefficient * (1.0 - 0.05 * critical),
                coefficient * (1.0 + 0.05 * critical),
                critical,
                "replicated_estimate",
                3,
                1,
            ]
            for mean, field, coefficient in (
                (1.0, 10.0, 1.0e-18),
                (2.0, 20.0, 1.0e-18),
            )
        ]
        physics_contract = {
            "schema": "swarm.elastic_energy_loss.v1",
            "symbol": "K_epsilon_el",
            "estimator": (
                "trajectory_event_energy_change_per_target_density_residence_time"
            ),
            "aggregation": "sum_target_fraction_times_process_coefficient",
            "source_eedf": ("same_sampled_trajectory_event_and_residence_measure"),
            "event_models": ["maxwellian_target_exact_binary_collision_isotropic"],
            "neutral_thermal_motion_model": (
                "maxwellian_relative_speed_exact_binary_collision"
            ),
            "gas_temperature_terms_included": True,
            "gas_temperature_K": 300.0,
            "sign_convention": "positive_is_net_electron_energy_loss",
            "uncertainty": "independent_replica_student_t_95",
        }
    else:
        rows = [
            [
                mean,
                field,
                field * 1.0e-21,
                1.0e-18,
                "",
                "",
                "",
                "",
                "",
                "deterministic_operator_moment",
                1,
                0,
            ]
            for mean, field in ((1.0, 10.0), (2.0, 20.0))
        ]
        physics_contract = {
            "schema": "swarm.elastic_energy_loss.v1",
            "symbol": "K_epsilon_el",
            "estimator": ("same_discrete_elastic_collision_operator_energy_moment"),
            "operator": (
                "native_finite_volume_scharfetter_gummel_elastic_A_D_zero_field"
            ),
            "operator_isolation": (
                "zero_field_reassembly_from_same_elastic_A_D_coefficients"
            ),
            "source_eedf": "same_solved_eedf",
            "neutral_thermal_motion_model": ("finite_temperature_fokker_planck"),
            "gas_temperature_terms_included": True,
            "gas_temperature_K": 300.0,
            "sign_convention": "positive_is_net_electron_energy_loss",
            "uncertainty": "deterministic_kinetic_solver",
        }
    table = bundle / "elastic_energy_loss_vs_mean_energy.csv"
    _write_csv(table, columns, rows)
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["source"] = source
    manifest["source_policy"]["transport_definition"] = (
        "mc_flux_particle_tracking_fixed_population"
        if source == "monte_carlo"
        else "flux"
    )
    if source == "monte_carlo":
        manifest["source_policy"].update(
            {
                "transport_eligible_anchor_points": 2,
                "transport_eligible_E_over_N_Td": [10.0, 20.0],
                "transport_excluded_E_over_N_Td": {},
            }
        )
        manifest["valid_ranges"]["E_over_N_Td"] = [10.0, 20.0]
        quality_columns = [
            "E_over_N_Td",
            "passed",
            "eedf_normalization_error",
            *mc_bundle_validation.REQUIRED_MC_QUALITY_COLUMNS,
        ]
        quality_values = {
            "aggregate_quality_passed": 1,
            "failure_reasons_json": "[]",
            "aggregate_failure_reasons_json": "[]",
            "mobility_rse": 0.01,
            "diffusion_L_rse": 0.01,
            "diffusion_T_rse": 0.01,
            "energy_mobility_rse": 0.01,
            "energy_diffusion_L_rse": 0.01,
            "energy_diffusion_T_rse": 0.01,
            "max_major_rate_rse": 0.01,
            "valid_replicates": 3,
            "uncertainty_available": 1,
            "solver_diagnostics_available": 1,
            "solver_transport_qualified": 1,
            "solver_transport_replicates": 3,
            "solver_origin_stationarity_limiting_field": "mean_energy_eV",
            "solver_origin_stationarity_limiting_relative_ci95_bound": 0.05,
            "solver_mean_energy_stationarity_relative_ci95_bound": 0.05,
            "solver_mobility_stationarity_absolute_log_drift": 0.02,
            "solver_lag_convergence_max_relative_ci95_bound": 0.05,
            "solver_transport_mean_energy_max_relative_ci95_bound": 0.05,
            "solver_population_growth_max_relative_ci95_bound": "",
            "solver_population_growth_gate_mode": "",
            "solver_population_growth_poisson_interval_max_ratio": "",
            "solver_population_growth_sparse_max_metric": "",
            "solver_origin_stationarity_limiting_relative_tolerance": 0.10,
            "solver_lag_convergence_relative_tolerance": 0.25,
            "solver_transport_mean_energy_relative_tolerance": 0.10,
            "solver_population_growth_relative_tolerance": "",
            "quality_source": "monte_carlo_aggregate_quality",
        }
        quality_table = bundle / "quality.csv"
        _write_csv(
            quality_table,
            quality_columns,
            [
                [
                    field,
                    1,
                    0.0,
                    *(
                        quality_values[name]
                        for name in mc_bundle_validation.REQUIRED_MC_QUALITY_COLUMNS
                    ),
                ]
                for field in (10.0, 20.0)
            ],
        )
        manifest["tables"]["quality.csv"].update(
            {
                "columns": quality_columns,
                "sha256": sha256(quality_table.read_bytes()).hexdigest(),
            }
        )
        sampling_plan = [
            {
                "e_over_n_Td": field,
                "particles": 64,
                "warmup_collisions": 100,
                "max_collisions": 200,
                "tail_max_collisions": 200,
                "replicas": 3,
                "transport_correlation_lag_barriers": 64,
                "transport_estimator": "single_field",
            }
            for field in (10.0, 20.0)
        ]
        sampling_json = json.dumps(
            sampling_plan,
            sort_keys=True,
            separators=(",", ":"),
        )
        manifest["hashes"].update(
            {
                "mc_sampling_plan_json": sampling_json,
                "mc_sampling_plan_sha256": sha256(
                    sampling_json.encode("utf-8")
                ).hexdigest(),
                "mc_transport_estimator_schema_version": (
                    mc_evidence.DIRECT_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION
                ),
                "mc_eedf_estimator_schema_version": (
                    mc_evidence.MC_EEDF_ESTIMATOR_SCHEMA_VERSION
                ),
                "mc_solver_source_sha256": monte_carlo_source_sha256(),
            }
        )
    manifest["tables"][table.name] = {
        "argument": "mean_energy_eV",
        "columns": columns,
        "artifact_role": "canonical_comsol_elastic_energy_loss_input",
        "sha256": sha256(table.read_bytes()).hexdigest(),
        "physics_contract": physics_contract,
    }
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    text = mapping_path.read_text(encoding="utf-8")
    text = (
        text.replace(
            "  expected_source: two_term",
            f"  expected_source: {source}"
            + (
                "\n  expected_mc_transport_estimator_schema: direct_mc_transport.v4"
                if source == "monte_carlo"
                else ""
            ),
        )
        .replace(
            "  expected_transport_definition: flux",
            "  expected_transport_definition: "
            + manifest["source_policy"]["transport_definition"],
        )
        .replace(
            "  electron_transport: comsol_specify_all_restricted",
            "  electron_transport: swarm_hybrid_einstein_de"
            + (
                "\n  zero_field_isotropization_Td: 0.1"
                if source == "monte_carlo"
                else ""
            ),
        )
        .replace(
            "  reaction_model: external_rates",
            "  reaction_model: external_rates\n"
            "  elastic_energy_loss_model: external_solver_native",
        )
    )
    mapping_path.write_text(text, encoding="utf-8")
