"""Materialize dimensionless C1 and active COMSOL Function-EEDF tables."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from swarm_workflow._io import write_csv as _write_csv
from swarm_workflow.comsol.input.function_eedf import (
    MEAN_AXIS_RATE_IMPORTANCE_FRACTION,
    SOURCE_MEAN_RELATIVE_ERROR_LIMIT,
    RATE_IMPORTANCE_FRACTION,
    RATE_SCALED_ERROR_TOLERANCE,
    SHAPE_TOTAL_VARIATION_TOLERANCE,
    C1FunctionEedf,
    ComsolEedfImportContract,
    ComsolFunctionEedfGrid,
    FunctionEedfError,
    build_c1_function_eedf,
    project_c1_function_eedf_to_comsol_grid,
    read_collision_rate_kernels,
)

from .contracts import (
    COMSOL_FUNCTION_EEDF_TABLE,
    FUNCTION_EEDF_TABLE,
    ComsolExportError,
)
from .manifest import _sha256_file


def _write_c1_function_eedf_tables(
    output_dir: Path,
    *,
    source: Path,
    grouped: dict[float, list[tuple[float, float, float]]],
    source_support: dict[str, Any],
    rate_kernel_path: Path | None = None,
) -> dict[str, dict[str, Any]]:
    try:
        representation = build_c1_function_eedf(grouped)
    except FunctionEedfError as exc:
        raise ComsolExportError(f"invalid Function EEDF: {exc}") from exc

    if (
        representation.normalization_error_max > 1.0e-8
        or representation.unit_mean_error_max > 1.0e-8
        or representation.minimum_value < -1.0e-14
    ):
        raise ComsolExportError(
            "C1 Function EEDF failed the input contract: "
            f"normalization error="
            f"{representation.normalization_error_max:.3e}, "
            f"unit-mean error={representation.unit_mean_error_max:.3e}, "
            f"minimum={representation.minimum_value:.3e}"
        )

    shape_columns = representation.shape_columns

    rows = [
        {
            "dimensionless_energy": float(x_value),
            **{
                column: float(shape[point_index])
                for column, shape in zip(
                    shape_columns,
                    representation.shapes,
                    strict=True,
                )
            },
        }
        for point_index, x_value in enumerate(representation.dimensionless_energy)
    ]

    name = FUNCTION_EEDF_TABLE

    columns = ("dimensionless_energy", *shape_columns)

    _write_csv(output_dir / name, columns, rows)

    try:
        rate_kernels = (
            read_collision_rate_kernels(rate_kernel_path)
            if rate_kernel_path is not None
            else ()
        )
        projected = project_c1_function_eedf_to_comsol_grid(
            representation,
            rate_kernels=rate_kernels,
        )
    except FunctionEedfError as exc:
        raise ComsolExportError(
            f"invalid projected COMSOL Function EEDF: {exc}"
        ) from exc

    # Plasma Function-EEDF consumes one two-argument spreadsheet contract.
    # Solver identity is provenance, never a serialization switch.
    import_contract = ComsolEedfImportContract("spreadsheet")
    active_name = COMSOL_FUNCTION_EEDF_TABLE

    active_columns = (
        "electron_energy_eV",
        "mean_energy_eV",
        "eepf_eV_m32",
    )

    _write_csv(
        output_dir / active_name,
        active_columns,
        (
            {
                "electron_energy_eV": float(energy),
                "mean_energy_eV": float(mean_energy),
                "eepf_eV_m32": float(projected.values_eV_m32[mean_index, energy_index]),
            }
            for mean_index, mean_energy in enumerate(projected.mean_energies_eV)
            for energy_index, energy in enumerate(projected.electron_energies_eV)
        ),
    )

    common_metadata = _c1_evidence_metadata(
        representation=representation,
        shape_columns=shape_columns,
        source=source,
        source_support=source_support,
    )
    active_metadata = _active_c1_metadata(
        active_columns=active_columns,
        active_name=active_name,
        import_contract=import_contract,
        name=name,
        output_dir=output_dir,
        projected=projected,
        representation=representation,
        source=source,
        source_support=source_support,
    )

    artifacts = {
        name: {
            **common_metadata,
            "columns": list(columns),
            "format": "csv",
            "canonical_comsol_input": False,
            "artifact_role": "canonicalized_function_eedf_evidence",
            "sha256": _sha256_file(output_dir / name),
        },
        active_name: active_metadata,
    }

    return artifacts


def _c1_evidence_metadata(
    *,
    representation: C1FunctionEedf,
    shape_columns: tuple[str, ...],
    source: Path,
    source_support: dict[str, Any],
) -> dict[str, Any]:
    return {
        "argument": "electron_energy_eV,mean_energy_eV",
        "argument_order": ["electron_energy_eV", "mean_energy_eV"],
        "derived_for_comsol_interpolation": True,
        "definition": "x=electron_energy_eV/mean_energy_eV; h_i(x)=mean_energy_i^(3/2)*f0_i(mean_energy_i*x); each h_i is a moment-preserving anchor shape",
        "source_reconstruction": "conservative_bin_probability_mass_to_finite_piecewise_f0",
        "source_artifact": source.name,
        "source_sha256": _sha256_file(source),
        "representation": "dimensionless_shape_preserving_c1_convex",
        "grid": "common_dimensionless_energy_with_wide_shape_columns",
        "dimensionless_energy_grid_points": len(representation.dimensionless_energy),
        "mean_energy_cases": len(representation.mean_energies_eV),
        "anchor_mean_energy_eV": representation.mean_energies_eV.tolist(),
        "shape_columns": list(shape_columns),
        "normalization_error_max": representation.normalization_error_max,
        "mean_energy_relative_error_max": representation.unit_mean_error_max,
        "nonnegative_minimum": representation.minimum_value,
        "source_mean_energy_relative_error_max": representation.source_mean_relative_error_max,
        "source_mean_energy_relative_error_limit": SOURCE_MEAN_RELATIVE_ERROR_LIMIT,
        "source_energy_max_eV": representation.source_energy_max_eV,
        "energy_axis_interpolation": "shape_preserving_piecewise_cubic_Hermite_C1",
        "mean_energy_axis_interpolation": "offline_C1_shape_family; active_table_uses_error_controlled_materialized_mean_nodes",
        "constraint_invariant": "each h_i has integral(sqrt(x)*h_i,dx)=1 and integral(x^(3/2)*h_i,dx)=1",
        "source_support": source_support,
        "tail_policy": "retain_all_source_bins_then_zero_beyond_source_support",
        "units": {
            "dimensionless_energy": "1",
            **{column: "1" for column in shape_columns},
        },
    }


def _active_c1_metadata(
    *,
    active_columns: tuple[str, str, str],
    active_name: str,
    import_contract: ComsolEedfImportContract,
    name: str,
    output_dir: Path,
    projected: ComsolFunctionEedfGrid,
    representation: C1FunctionEedf,
    source: Path,
    source_support: dict[str, Any],
) -> dict[str, Any]:
    return {
        "argument": "electron_energy_eV,mean_energy_eV",
        "argument_order": ["electron_energy_eV", "mean_energy_eV"],
        "artifact_role": "canonical_comsol_function_eedf_input",
        "canonical_comsol_input": True,
        "columns": list(active_columns),
        "definition": "physical f0(E,m) on independently adaptive energy and mean-energy axes; every solver anchor is retained and nonnegative row projection preserves normalization and mean energy",
        "format": import_contract.format,
        "representation": import_contract.representation,
        "source_artifact": source.name,
        "source_sha256": _sha256_file(source),
        "projection_evidence_artifact": name,
        "projection_evidence_sha256": _sha256_file(output_dir / name),
        "anchor_mean_energy_eV": representation.mean_energies_eV.tolist(),
        "energy_grid_points": len(projected.electron_energies_eV),
        "mean_energy_grid_points": len(projected.mean_energies_eV),
        "grid_shape": [
            len(projected.mean_energies_eV),
            len(projected.electron_energies_eV),
        ],
        "grid_axis_order": ["electron_energy_eV", "mean_energy_eV"],
        "row_order": import_contract.row_order,
        "electron_energy_range_eV": [
            float(projected.electron_energies_eV[0]),
            float(projected.electron_energies_eV[-1]),
        ],
        "energy_support_policy": "full_source_support_no_solver_specific_tail_cutoff",
        "source_energy_max_eV": representation.source_energy_max_eV,
        "support_expansion_beyond_source_eV": max(
            0.0,
            float(projected.electron_energies_eV[-1])
            - representation.source_energy_max_eV,
        ),
        "mean_energy_range_eV": [
            float(projected.mean_energies_eV[0]),
            float(projected.mean_energies_eV[-1]),
        ],
        "projected_normalization_error_max": projected.normalization_error_max,
        "projected_mean_energy_relative_error_max": projected.mean_energy_relative_error_max,
        "projected_nonnegative_minimum": projected.minimum_value,
        "shape_total_variation_error_max": (projected.shape_total_variation_error_max),
        "shape_total_variation_tolerance": SHAPE_TOTAL_VARIATION_TOLERANCE,
        "rate_scaled_error_max": projected.rate_scaled_error_max,
        "rate_scaled_error_tolerance": (
            RATE_SCALED_ERROR_TOLERANCE if projected.rate_kernel_count else None
        ),
        "rate_importance_fraction": (
            RATE_IMPORTANCE_FRACTION if projected.rate_kernel_count else None
        ),
        "rate_error_norm": (
            "abs(projected-reference)/max(reference,global_kernel_peak*importance_fraction); values_unmodified"
            if projected.rate_kernel_count
            else None
        ),
        "rate_kernel_count": projected.rate_kernel_count,
        "adaptive_energy_axis": True,
        "adaptive_mean_energy_axis": True,
        "mean_energy_axis": "error_controlled_c1_materialization_with_solver_anchors",
        "solver_anchor_mean_energy_count": len(representation.mean_energies_eV),
        "materialized_mean_nodes_include_all_solver_anchors": True,
        "mean_axis_shape_total_variation_error_max": (
            projected.mean_axis_shape_total_variation_error_max
        ),
        "mean_axis_shape_total_variation_tolerance": (SHAPE_TOTAL_VARIATION_TOLERANCE),
        "mean_axis_rate_scaled_error_max": (projected.mean_axis_rate_scaled_error_max),
        "mean_axis_rate_scaled_error_tolerance": (
            RATE_SCALED_ERROR_TOLERANCE if projected.rate_kernel_count else None
        ),
        "mean_axis_rate_importance_fraction": (
            MEAN_AXIS_RATE_IMPORTANCE_FRACTION if projected.rate_kernel_count else None
        ),
        "interpolation_contract": "tensor_product_linear on error-controlled axes; every materialized mean row has unit normalization and first moment equal to its mean, so linear mean interpolation preserves both moments",
        "source_support": source_support,
        "mean_energy_support_policy": (
            "no_clamp_no_synthetic_extension; values outside the tabulated "
            "range require new solver anchors"
        ),
        "units": {
            "electron_energy_eV": "eV",
            "mean_energy_eV": "eV",
            "eepf_eV_m32": "1/eV^(3/2)",
        },
        "sha256": _sha256_file(output_dir / active_name),
        "comsol_import": import_contract.import_settings(),
    }
