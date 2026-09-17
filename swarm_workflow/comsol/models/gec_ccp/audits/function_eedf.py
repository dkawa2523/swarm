"""Audit saved GEC-CCP Function-EEDF runs and reconstructed rates."""

from __future__ import annotations

import hashlib
from pathlib import Path
import re
from typing import Any
from zipfile import ZipFile

import numpy as np

from swarm_workflow._io import write_json
from swarm_workflow.comsol.input.function_eedf import ComsolFunctionEedfGrid
from swarm_workflow.comsol.input.mean_energy import MeanEnergyArgument
from swarm_workflow.comsol.models.gec_ccp.closure import (
    _function_eedf_rate_reactions,
    _uses_function_eedf,
)
from swarm_workflow.comsol.models.gec_ccp.audits.conservation import (
    _read_comsol_numeric_csv,
)
from swarm_workflow.comsol.models.gec_ccp.contracts import (
    FUNCTION_EEDF_PREINTEGRATED_INELASTIC,
    GecCcpPlan,
    GecCcpWorkflowError,
    GecReactionSpec,
)
from swarm_workflow.comsol.models.gec_ccp.validation.eedf import (
    FUNCTION_EEDF_RATE_NORMALIZED_RMSE_LIMIT,
    FUNCTION_EEDF_RATE_RELATIVE_ERROR_LIMIT,
    FUNCTION_EEDF_RATE_SIGNIFICANCE_FRACTION,
    _audit_function_eedf_source_rates,
    _audit_gec_argon_cross_section_identity,
    _integrate_native_function_eedf_rate,
    _native_function_eedf_moment_audit,
    _rate_significance,
    _rate_significance_metadata,
    _reaction_cross_sections_from_model_xml,
    _read_function_eedf,
)
from swarm_workflow.comsol.models.gec_ccp.closure_arguments import (
    closure_argument_range,
)
from swarm_workflow.comsol.models.gec_ccp.validation.joint_consistency import (
    _active_closure_mean_energy_support,
)
from swarm_workflow.comsol.models.gec_ccp.mph import _property_values
from swarm_workflow.comsol.models.gec_ccp.audits.saved_model import (
    _audit_saved_reaction_handling,
)


def _function_eedf_rate_range_quality(values: np.ndarray) -> dict[str, Any]:
    """Accept finite rates with negative noise below the significance floor."""

    finite = bool(np.all(np.isfinite(values)))
    minimum = float(np.min(values))
    maximum = float(np.max(values))
    negative_tolerance = max(maximum, 0.0) * FUNCTION_EEDF_RATE_SIGNIFICANCE_FRACTION
    passed = bool(finite and maximum > 0.0 and minimum >= -negative_tolerance)
    return {
        "passed": passed,
        "finite": finite,
        "minimum_m3_per_mol_s": minimum,
        "maximum_m3_per_mol_s": maximum,
        "negative_tolerance_m3_per_mol_s": negative_tolerance,
        "negative_tolerance_fraction_of_process_peak": (
            FUNCTION_EEDF_RATE_SIGNIFICANCE_FRACTION
        ),
    }


def audit_gec_ccp_function_eedf_run(plan: GecCcpPlan) -> Path:
    """Verify the saved EEDF binding and its converged reaction rates.

    The audit is deliberately independent of solver success alone.  It checks
    the saved MPH binding, native or preintegrated reaction mode, embedded data,
    finite phase-local rates, and that every exported phase-local mean energy
    lies strictly inside the tabulated Function-EEDF range.
    """

    mapping = plan.mapping
    spec = mapping.closure.function_eedf
    if not _uses_function_eedf(mapping.closure) or spec is None:
        raise GecCcpWorkflowError(
            "Function-EEDF audit requires closure.reaction_model=function_eedf"
        )
    phase_path = (
        mapping.results.output_directory / "swarm_tables" / "domain_phase_closure.csv"
    )
    headers, values = _read_comsol_numeric_csv(phase_path)
    energy_indices = [
        index for index, header in enumerate(headers) if header.startswith("ptp.ebar ")
    ]
    function_reactions = _function_eedf_rate_reactions(mapping)
    if not function_reactions:
        raise GecCcpWorkflowError(
            "Function-EEDF closure has no active EEDF-derived reaction"
        )
    rate_indices = {
        reaction.feature: [
            index
            for index, header in enumerate(headers)
            if header.startswith(f"ptp.kf_{reaction_index} ")
        ]
        for reaction_index, reaction in function_reactions
    }
    if not energy_indices or any(
        len(indices) != len(energy_indices) for indices in rate_indices.values()
    ):
        raise GecCcpWorkflowError(
            "domain phase closure export lacks mean energy or reaction rates"
        )
    closure_columns = [
        *energy_indices,
        *(index for indices in rate_indices.values() for index in indices),
    ]
    if not np.all(np.isfinite(values[:, closure_columns])):
        raise GecCcpWorkflowError(
            "Function-EEDF phase closure export contains a nonfinite value"
        )

    function_table_path = mapping.bundle.path / spec.table
    function_grid = _read_function_eedf(function_table_path)
    tabulated_means = function_grid.mean_energies_eV
    table_min = float(np.min(tabulated_means))
    table_max = float(np.max(tabulated_means))
    operating_energy = values[:, energy_indices]
    operating_min = float(np.min(operating_energy))
    operating_max = float(np.max(operating_energy))
    inside_range = table_min <= operating_min and operating_max <= table_max

    rate_ranges: dict[str, dict[str, Any]] = {}
    rates_valid = True
    active_reactions = tuple(reaction for _, reaction in function_reactions)
    for reaction in active_reactions:
        reaction_values = values[:, rate_indices[reaction.feature]]
        quality = _function_eedf_rate_range_quality(reaction_values)
        rates_valid = rates_valid and bool(quality["passed"])
        rate_ranges[reaction.process_type] = quality

    try:
        with ZipFile(mapping.model.output_mph) as archive:
            model_xml = archive.read("dmodel.xml").decode("utf-8", errors="replace")
            embedded_resources = sorted(
                name for name in archive.namelist() if name.startswith("resources/")
            )
            embedded_resource_hashes = {
                name: hashlib.sha256(archive.read(name)).hexdigest()
                for name in embedded_resources
            }
    except (OSError, KeyError) as exc:
        raise GecCcpWorkflowError(
            f"cannot audit saved Function-EEDF MPH: {mapping.model.output_mph}"
        ) from exc
    eedf_values = _property_values(model_xml, "eedf")
    active_eedf = next(
        (
            value
            for value in reversed(eedf_values)
            if value != "FromPhysicsInterfaceProperty"
        ),
        None,
    )
    saved_reaction_handling = _audit_saved_reaction_handling(mapping, model_xml)
    binding_counts = saved_reaction_handling["actual_binding_counts"]
    cross_section_binding_count = binding_counts["UseCrossSectionData"]
    rate_constant_binding_count = binding_counts["RateConstant"]
    reaction_eedf_binding_count = saved_reaction_handling[
        "active_reaction_eedf_bindings"
    ]
    function_block_match = re.search(
        rf'<FunctionFeature\b[^>]*op="Interpolation"[^>]*tag="'
        rf'{re.escape(spec.function_tag)}"[^>]*>(.*?)</FunctionFeature>',
        model_xml,
        flags=re.DOTALL,
    )
    function_block = function_block_match.group(1) if function_block_match else ""
    function_table_sha256 = hashlib.sha256(function_table_path.read_bytes()).hexdigest()
    matching_embedded_resources = sorted(
        name
        for name, digest in embedded_resource_hashes.items()
        if digest == function_table_sha256
    )
    column_type_match = re.search(
        r'valueMatrix="([^"]*)" name="p:columnType"', function_block
    )
    column_types = column_type_match.group(1) if column_type_match else ""
    two_dimensional = (
        column_types.count("'arg'") == 2 and column_types.count("'value'") >= 1
    )
    function_resource_valid = bool(
        function_block
        and f'value="{Path(spec.table).name}" name="p:importedname"' in function_block
        and 'value="file" name="p:source"' in function_block
        and two_dimensional
        and spec.function_tag in function_block
        and matching_embedded_resources
    )
    cross_section_tables = _reaction_cross_sections_from_model_xml(
        model_xml, mapping.reactions
    )
    cross_section_identity = _audit_gec_argon_cross_section_identity(
        cross_section_tables, mapping.reactions
    )
    rate_arguments: dict[str, MeanEnergyArgument] = {}
    if mapping.closure.reaction_model == FUNCTION_EEDF_PREINTEGRATED_INELASTIC:
        supports = _active_closure_mean_energy_support(mapping)[
            "rate_supports_mean_energy_eV"
        ]
        rate_arguments = {
            name: MeanEnergyArgument(*closure_argument_range(mapping, support))
            for name, support in supports.items()
        }
    rate_reintegration = _audit_function_eedf_rates(
        function_grid,
        cross_section_tables,
        active_reactions,
        operating_energy,
        values,
        rate_indices,
        rate_arguments=rate_arguments,
    )
    source_rate_consistency = _audit_function_eedf_source_rates(
        mapping,
        function_grid,
        cross_section_tables,
        active_reactions,
    )
    model_binding_valid = bool(
        active_eedf == spec.function_tag
        and saved_reaction_handling["passed"]
        and function_resource_valid
    )

    status = (
        "passed"
        if (
            inside_range
            and rates_valid
            and model_binding_valid
            and cross_section_identity["passed"]
            and rate_reintegration["passed"]
            and source_rate_consistency["passed"]
        )
        else "failed"
    )
    audit_path = plan.output_directory / "function_eedf_audit.json"
    write_json(
        audit_path,
        {
            "stage": "audit-gec-ccp-function-eedf",
            "status": status,
            "model": {
                "path": str(mapping.model.output_mph),
                "sha256": hashlib.sha256(
                    mapping.model.output_mph.read_bytes()
                ).hexdigest(),
                "active_eedf_token": active_eedf,
                "expected_function_tag": spec.function_tag,
                "cross_section_reaction_bindings": cross_section_binding_count,
                "rate_constant_reaction_bindings": rate_constant_binding_count,
                "reaction_eedf_interface_bindings": reaction_eedf_binding_count,
                "reaction_handling": saved_reaction_handling,
                "function_table_sha256": function_table_sha256,
                "function_resource_exact_match": function_resource_valid,
                "function_operation": ("Interpolation" if function_block else None),
                "function_arguments": 2 if two_dimensional else None,
                "matching_embedded_resources": matching_embedded_resources,
                "binding_valid": model_binding_valid,
            },
            "function_table": {
                "path": str(function_table_path),
                "sha256": function_table_sha256,
                "mean_energy_range_eV": [table_min, table_max],
                "representation": spec.interpolation,
                "continuity": "C0",
                "projected_moment_audit": (
                    _native_function_eedf_moment_audit(function_grid)
                ),
            },
            "converged_operating_range": {
                "phase_local_mean_energy_eV": [operating_min, operating_max],
                "strictly_inside_function_range": inside_range,
                "exported_domain_nodes": int(values.shape[0]),
                "exported_phase_samples_per_node": len(energy_indices),
            },
            "rates": {
                "closure_rate_source": mapping.closure.reaction_model,
                "active_function_reactions": [
                    reaction.process_type for reaction in active_reactions
                ],
                "finite_nonnegative_with_significance_tolerance": rates_valid,
                "ranges": rate_ranges,
                "cross_section_identity": cross_section_identity,
                "independent_reintegration": rate_reintegration,
                "swarm_source_consistency": source_rate_consistency,
            },
        },
    )
    if status != "passed":
        raise GecCcpWorkflowError(
            "converged COMSOL result failed the Function-EEDF usage audit; "
            f"see {audit_path}"
        )
    return audit_path


def _audit_function_eedf_rates(
    function_grid: ComsolFunctionEedfGrid,
    cross_sections: dict[str, tuple[np.ndarray, np.ndarray]],
    reactions: tuple[GecReactionSpec, ...],
    operating_energy: np.ndarray,
    exported_values: np.ndarray,
    rate_indices: dict[str, list[int]],
    *,
    rate_arguments: dict[str, MeanEnergyArgument] | None = None,
) -> dict[str, Any]:
    flat_energy = operating_energy.reshape(-1)
    order = np.argsort(flat_energy, kind="stable")
    sample_count = min(257, order.size)
    selected = order[
        np.unique(np.rint(np.linspace(0, order.size - 1, sample_count)).astype(int))
    ]
    sampled_means = flat_energy[selected]
    avogadro = 6.02214076e23
    process_results: dict[str, Any] = {}
    for reaction in reactions:
        cross_energy, cross_sigma = cross_sections[reaction.process_type]
        argument = (rate_arguments or {}).get(reaction.process_type)
        reference_means = argument.values(sampled_means) if argument else sampled_means
        expected = np.asarray(
            [
                _integrate_native_function_eedf_rate(
                    function_grid,
                    float(mean),
                    cross_energy,
                    cross_sigma,
                )
                for mean in reference_means
            ],
            dtype=float,
        )
        actual = (
            exported_values[:, rate_indices[reaction.feature]].reshape(-1)[selected]
            / avogadro
        )
        difference = actual - expected
        process_peak, significance_floor, significant = _rate_significance(
            expected, actual
        )
        relative = np.abs(difference[significant]) / np.maximum(
            np.abs(expected[significant]), 1.0e-300
        )
        normalized_rmse = float(
            np.sqrt(np.mean(np.square(difference)))
            / max(np.sqrt(np.mean(np.square(expected))), 1.0e-300)
        )
        maximum_relative = float(np.max(relative)) if relative.size else 0.0
        passed = (
            normalized_rmse <= FUNCTION_EEDF_RATE_NORMALIZED_RMSE_LIMIT
            and maximum_relative <= FUNCTION_EEDF_RATE_RELATIVE_ERROR_LIMIT
        )
        process_results[reaction.process_type] = {
            "reference_mean_energy_range_eV": [
                float(np.min(reference_means)),
                float(np.max(reference_means)),
            ],
            "coefficient_argument_range_eV": (
                [argument.minimum_eV, argument.maximum_eV] if argument else None
            ),
            "rate_unit": "m^3/s per target particle",
            "samples": int(expected.size),
            "reference_rate_range_m3_s": [
                float(np.min(expected)),
                float(np.max(expected)),
            ],
            "comsol_rate_range_m3_s": [float(np.min(actual)), float(np.max(actual))],
            "normalized_rmse": normalized_rmse,
            "normalized_rmse_limit": FUNCTION_EEDF_RATE_NORMALIZED_RMSE_LIMIT,
            "normalized_rmse_scope": "all_samples",
            "maximum_relative_error_significant": maximum_relative,
            "maximum_relative_error_limit": (FUNCTION_EEDF_RATE_RELATIVE_ERROR_LIMIT),
            **_rate_significance_metadata(
                process_peak=process_peak,
                significance_floor=significance_floor,
                significant=significant,
            ),
            "passed": passed,
        }
    return {
        "passed": all(item["passed"] for item in process_results.values()),
        "method": (
            "independent sigma(E)*v(E)*sqrt(E)*f0(E,meanE) quadrature "
            "of the active physical 2D linear projection; COMSOL molar kf "
            "divided by Avogadro constant"
        ),
        "processes": process_results,
    }
