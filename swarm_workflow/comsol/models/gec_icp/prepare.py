"""Prepare a provenance-bound GEC-ICP COMSOL run."""

from __future__ import annotations

import csv
from dataclasses import asdict
import hashlib
from pathlib import Path
from typing import Any

from swarm_workflow._io import write_json
from swarm_workflow.comsol.input.eedf_audit import (
    prepare_comsol_eedf_audit,
    render_comsol_eedf_audit_java,
)
from swarm_workflow.comsol.java import compose_java_mains

from .bundle import validate_gec_icp_bundle
from .closure_readback import (
    READBACK_CLASS,
    generate_closure_readback_java,
    prepare_gec_icp_closure_readback,
)
from .cross_sections import (
    CANONICAL_ARGON_CROSS_SECTIONS,
    COMSOL_ARGON_IMPORT_REFERENCE,
)
from .java import (
    APPLY_CLASS,
    CLOSURE_CLASS,
    NATIVE_EEDF_AUDIT_CLASS,
    SOLVE_CLASS,
    generate_apply_java,
    generate_solve_java,
)
from .mph import validate_gec_icp_mph
from .run_contracts import GecIcpPlan, GecIcpWorkflowError
from .run_mapping import load_gec_icp_run_mapping


_RESULT_NAMES = (
    "convergence_volume.csv",
    "convergence_mean_energy_min.csv",
    "convergence_mean_energy_max.csv",
    "convergence_coil_power.csv",
    "convergence.json",
    "closure_readback.json",
    "comsol_eedf_audit.json",
)


def prepare_gec_icp_run(
    mapping_path: str | Path,
    *,
    bundle_path: str | Path | None = None,
    write_java: bool = True,
) -> GecIcpPlan:
    """Validate model and solver input, then materialize an execution plan."""

    mapping = load_gec_icp_run_mapping(mapping_path, bundle_path=bundle_path)
    mph_contract = validate_gec_icp_mph(mapping.model_mapping)
    bundle_evidence = validate_gec_icp_bundle(
        mapping.bundle.path,
        expected_source=mapping.bundle.expected_source,
        expected_pressure_Pa=mapping.run.pressure_Pa,
        expected_temperature_K=mapping.run.gas_temperature_K,
        expected_transport_definition=mapping.bundle.expected_transport_definition,
        expected_mc_qualification_profile=(
            mapping.bundle.expected_mc_qualification_profile
        ),
    )
    output = mapping.output_directory
    output.mkdir(parents=True, exist_ok=True)
    mapping.run.output_mph.parent.mkdir(parents=True, exist_ok=True)
    mapping.log_path.mkdir(parents=True, exist_ok=True)
    closure_java = output / f"{CLOSURE_CLASS}.java"
    apply_java = output / f"{APPLY_CLASS}.java"
    solve_java = output / f"{SOLVE_CLASS}.java"
    closure_readback = prepare_gec_icp_closure_readback(
        mapping,
        output_directory=output,
        write_java=False,
    )
    eedf_table = mapping.bundle.path / mapping.closure.function_eedf.table
    native_eedf_audit = prepare_comsol_eedf_audit(
        model_path=mapping.run.output_mph,
        table_path=eedf_table,
        output_directory=output / "native_eedf_audit",
        component=mapping.model_mapping.model.component,
        physics=mapping.model_mapping.model.plasma_physics,
        function_tag=mapping.closure.function_eedf.function_tag,
        class_name=NATIVE_EEDF_AUDIT_CLASS,
        require_model_exists=False,
        write_java=False,
    )
    closure_support_java = (
        apply_java,
        closure_readback.java_path,
        native_eedf_audit.java_path,
    )
    if write_java:
        with native_eedf_audit.query_path.open(encoding="utf-8", newline="") as stream:
            audit_rows = list(csv.DictReader(stream))
        native_source = render_comsol_eedf_audit_java(
            class_name=NATIVE_EEDF_AUDIT_CLASS,
            model_path=mapping.run.output_mph,
            component=mapping.model_mapping.model.component,
            physics=mapping.model_mapping.model.plasma_physics,
            function_tag=mapping.closure.function_eedf.function_tag,
            audit_rows=audit_rows,
            table_sha256=native_eedf_audit.table_sha256,
            query_sha256=native_eedf_audit.query_sha256,
        )
        stage_sources = (
            (APPLY_CLASS, apply_java, generate_apply_java(mapping)),
            (
                READBACK_CLASS,
                closure_readback.java_path,
                generate_closure_readback_java(closure_readback),
            ),
            (NATIVE_EEDF_AUDIT_CLASS, native_eedf_audit.java_path, native_source),
        )
        for _, path, source in stage_sources:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(source, encoding="utf-8")
        closure_java.write_text(
            compose_java_mains(
                CLOSURE_CLASS,
                tuple((class_name, source) for class_name, _, source in stage_sources),
            ),
            encoding="utf-8",
        )
        solve_java.write_text(generate_solve_java(mapping), encoding="utf-8")
    expected = tuple(output / name for name in _RESULT_NAMES)
    plan_json = output / "gec_icp_plan.json"
    input_mph = mapping.model_mapping.model.input_mph
    if not COMSOL_ARGON_IMPORT_REFERENCE.is_file():
        raise GecIcpWorkflowError(
            "repository COMSOL cross-section import reference is missing: "
            f"{COMSOL_ARGON_IMPORT_REFERENCE}"
        )
    plan_payload = {
        "schema": "swarm.gec_icp_comsol_plan.v1",
        "status": "ready",
        "model_adapter": "gec_icp",
        "model_mapping": str(mapping.model_mapping_path),
        "mapping": str(mapping.path),
        "source": mapping.bundle.expected_source,
        "input_mph": _artifact(input_mph),
        "cross_section_evidence": {
            "canonical_ground_state_csv": _artifact(CANONICAL_ARGON_CROSS_SECTIONS),
            "repository_comsol_import_reference": _artifact(
                COMSOL_ARGON_IMPORT_REFERENCE
            ),
            "offline_mph_contract": (
                "eir1/eir2/eir4 species, process, type, threshold, energy, and "
                "cross-section arrays match canonical_ground_state_csv; "
                "COMSOL-owned eir3/eir5 arrays are structurally validated"
            ),
        },
        "mph_contract": _jsonable(asdict(mph_contract)),
        "bundle": bundle_evidence,
        "closure_ownership": {
            "external_swarm": [
                "function_eedf_f0_of_energy_and_mean_energy",
                "reduced_electron_mobility",
            ],
            "comsol": [
                "particle_diffusion_from_SpecifyMueOnly",
                "electron_energy_transport",
                "all_electron_impact_rate_integrals",
                "elastic_energy_exchange",
                "Ar_and_Ars_chemistry",
                "induction_current_RF_field",
                "walls_and_heavy_species_transport",
            ],
            "approximation": (
                "homogeneous_DC_Swarm_tables_used_as_a_local_mean_energy_"
                "closure_for_the_13.56_MHz_frequency_transient_model"
            ),
        },
        "run": _jsonable(asdict(mapping.run)),
        "generated_java": {
            "closure_entry_point": (
                _artifact(closure_java) if write_java else str(closure_java)
            ),
            "closure_support_stages": [
                _artifact(path) if write_java else str(path)
                for path in closure_support_java
            ],
            "solve_and_export": (
                _artifact(solve_java) if write_java else str(solve_java)
            ),
        },
        "closure_readback": {
            "model": str(mapping.run.output_mph),
            "transport": _artifact(closure_readback.transport_path),
            "anchor_count": closure_readback.anchor_count,
            "report": str(closure_readback.report_path),
            "stage_order": (
                "apply_then_saved_model_readback_then_native_function_eedf_"
                "within_one_comsol_process"
            ),
        },
        "coefficient_sweep": False,
        "cold_start": True,
        "expected_results": [str(path) for path in expected],
    }
    write_json(plan_json, plan_payload)
    return GecIcpPlan(
        mapping=mapping,
        mph_contract=mph_contract,
        bundle_evidence=bundle_evidence,
        output_directory=output,
        plan_json=plan_json,
        closure_java=closure_java,
        closure_support_java=closure_support_java,
        closure_readback=closure_readback,
        native_eedf_audit=native_eedf_audit,
        solve_java=solve_java,
        expected_result_files=expected,
    )


def _artifact(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise GecIcpWorkflowError(f"planned input does not exist: {path}")
    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def _jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


__all__ = ["prepare_gec_icp_run"]
