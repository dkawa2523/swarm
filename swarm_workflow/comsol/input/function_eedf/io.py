"""Read and validate canonical COMSOL Function-EEDF artifacts."""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .contracts import (
    COMSOL_EEDF_COLUMNS,
    SOURCE_MEAN_RELATIVE_ERROR_LIMIT,
    C1FunctionEedf,
    ComsolEedfImportContract,
    ComsolFunctionEedfGrid,
    FunctionEedfError,
)
from .moments import (
    _pchip_minimum,
    pchip_weighted_moments,
    piecewise_linear_weighted_moments,
)


def read_c1_function_eedf(
    table_path: str | Path,
    metadata: dict[str, Any],
) -> C1FunctionEedf:
    """Read and validate the canonical wide shape table and its manifest entry."""

    path = Path(table_path)
    try:
        means = np.asarray(metadata["anchor_mean_energy_eV"], dtype=float)
        shape_columns = tuple(str(value) for value in metadata["shape_columns"])
        source_energy_max = float(metadata["source_energy_max_eV"])
        source_error = float(metadata["source_mean_energy_relative_error_max"])
        source_error_limit = float(metadata["source_mean_energy_relative_error_limit"])
        norm_error = float(metadata["normalization_error_max"])
        mean_error = float(metadata["mean_energy_relative_error_max"])
        minimum = float(metadata["nonnegative_minimum"])
    except (KeyError, TypeError, ValueError) as exc:
        raise FunctionEedfError("invalid canonical Function-EEDF metadata") from exc
    if metadata.get("representation") != "dimensionless_shape_preserving_c1_convex":
        raise FunctionEedfError("unsupported Function-EEDF representation")
    if len(means) < 2 or len(shape_columns) != len(means):
        raise FunctionEedfError("Function-EEDF anchor/shape counts differ")
    if (
        np.any(~np.isfinite(means))
        or np.any(means <= 0.0)
        or np.any(np.diff(means) <= 0.0)
    ):
        raise FunctionEedfError("Function-EEDF anchors are invalid")

    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        expected = ("dimensionless_energy", *shape_columns)
        if tuple(reader.fieldnames or ()) != expected:
            raise FunctionEedfError(
                "Function-EEDF wide-table columns do not match metadata"
            )
        rows = list(reader)
    try:
        x_values = np.asarray(
            [float(row["dimensionless_energy"]) for row in rows], dtype=float
        )
        shapes = tuple(
            np.asarray([float(row[column]) for row in rows], dtype=float)
            for column in shape_columns
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise FunctionEedfError(
            "Function-EEDF table contains a nonnumeric value"
        ) from exc
    if len(x_values) < 3 or x_values[0] != 0.0 or np.any(np.diff(x_values) <= 0.0):
        raise FunctionEedfError("Function-EEDF dimensionless grid is invalid")
    if any(len(shape) != len(x_values) for shape in shapes):
        raise FunctionEedfError("Function-EEDF shape length differs from its grid")
    if any(np.any(~np.isfinite(shape)) or np.any(shape < -1.0e-14) for shape in shapes):
        raise FunctionEedfError("Function-EEDF shape is nonfinite or negative")
    if (
        not all(
            math.isfinite(value)
            for value in (
                source_energy_max,
                source_error,
                source_error_limit,
                norm_error,
                mean_error,
                minimum,
            )
        )
        or source_energy_max <= 0.0
        or source_error_limit != SOURCE_MEAN_RELATIVE_ERROR_LIMIT
        or source_error > source_error_limit
        or norm_error > 1.0e-8
        or mean_error > 1.0e-8
        or minimum < -1.0e-14
    ):
        raise FunctionEedfError("Function-EEDF manifest audit failed")
    for shape in shapes:
        normalization, first_moment = pchip_weighted_moments(x_values, shape)
        if (
            abs(normalization - 1.0) > 1.0e-8
            or abs(first_moment / normalization - 1.0) > 1.0e-8
            or _pchip_minimum(x_values, shape) < -1.0e-14
            or shape[-2] != 0.0
            or shape[-1] != 0.0
        ):
            raise FunctionEedfError("Function-EEDF numerical audit failed")
    return C1FunctionEedf(
        mean_energies_eV=means,
        dimensionless_energy=x_values,
        shapes=shapes,
        source_energy_max_eV=source_energy_max,
        source_mean_relative_error_max=source_error,
        normalization_error_max=norm_error,
        unit_mean_error_max=mean_error,
        minimum_value=minimum,
    )


def read_comsol_function_eedf_grid(
    table_path: str | Path,
    metadata: Mapping[str, Any] | None = None,
) -> ComsolFunctionEedfGrid:
    """Read a rectangular physical EEDF from canonical spreadsheet serialization.

    The format does not permit scattered points, duplicate coordinates,
    reordered rows, or a different moment/positivity contract.
    """

    path = Path(table_path)
    contract = ComsolEedfImportContract.from_path(path)
    energy_grid, mean_grid, value_rows = _read_eedf_spreadsheet_axes(path)
    grid = _physical_eedf_grid(energy_grid, mean_grid, value_rows)
    if metadata is not None:
        _validate_comsol_function_eedf_metadata(grid, metadata, contract)
    return grid


def _read_eedf_spreadsheet_axes(
    path: Path,
) -> tuple[np.ndarray, np.ndarray, list[np.ndarray]]:
    energies: list[float] = []
    means: list[float] = []
    rows: list[np.ndarray] = []
    current_energy: list[float] = []
    current_values: list[float] = []
    current_mean: float | None = None

    def finish_row() -> None:
        if current_mean is None:
            return
        if not energies:
            energies.extend(current_energy)
        elif energies != current_energy:
            raise FunctionEedfError("Function-EEDF spreadsheet is not rectangular")
        means.append(current_mean)
        rows.append(np.asarray(current_values, dtype=float))
        current_energy.clear()
        current_values.clear()

    try:
        with path.open("r", encoding="utf-8-sig", newline="") as stream:
            reader = csv.reader(line for line in stream if line.strip())
            if tuple(next(reader)) != COMSOL_EEDF_COLUMNS:
                raise FunctionEedfError(
                    "Function-EEDF spreadsheet columns are not canonical"
                )
            for row in reader:
                if len(row) != 3:
                    raise FunctionEedfError(
                        "Function-EEDF spreadsheet requires exactly three columns"
                    )
                energy, mean, value = (float(item) for item in row)
                if current_mean is not None and mean != current_mean:
                    finish_row()
                current_mean = mean
                current_energy.append(energy)
                current_values.append(value)
            finish_row()
    except (OSError, ValueError, StopIteration) as exc:
        raise FunctionEedfError(f"invalid Function-EEDF spreadsheet: {exc}") from exc
    return np.asarray(energies), np.asarray(means), rows


def _physical_eedf_grid(
    energy_grid: np.ndarray,
    mean_grid: np.ndarray,
    value_rows: list[np.ndarray],
) -> ComsolFunctionEedfGrid:
    if len(mean_grid) < 2 or not value_rows:
        raise FunctionEedfError(
            "COMSOL Function-EEDF requires at least two mean-energy rows"
        )
    if (
        len(energy_grid) < 3
        or energy_grid[0] != 0.0
        or np.any(np.diff(energy_grid) <= 0.0)
        or np.any(np.diff(mean_grid) <= 0.0)
    ):
        raise FunctionEedfError(
            "COMSOL Function-EEDF axes must be strictly increasing and the "
            "energy axis must start at zero with at least three points"
        )
    if len(value_rows) != len(mean_grid) or any(
        len(row) != len(energy_grid) for row in value_rows
    ):
        raise FunctionEedfError(
            "COMSOL Function-EEDF structured data block is not rectangular"
        )
    values = np.vstack(value_rows)
    if (
        np.any(~np.isfinite(energy_grid))
        or np.any(~np.isfinite(mean_grid))
        or np.any(~np.isfinite(values))
        or np.any(energy_grid < 0.0)
        or np.any(mean_grid <= 0.0)
        or np.any(values < 0.0)
    ):
        raise FunctionEedfError(
            "COMSOL Function-EEDF values must have E>=0, meanE>0, and f0>=0"
        )

    normalization_error = 0.0
    mean_error = 0.0
    for requested_mean, row in zip(mean_grid, values, strict=True):
        normalization, first_moment = piecewise_linear_weighted_moments(
            energy_grid, row
        )
        if normalization <= 0.0:
            raise FunctionEedfError("COMSOL Function-EEDF row has zero normalization")
        normalization_error = max(normalization_error, abs(normalization - 1.0))
        mean_error = max(
            mean_error,
            abs(first_moment / normalization - requested_mean) / requested_mean,
        )
    minimum = float(np.min(values))
    return ComsolFunctionEedfGrid(
        electron_energies_eV=energy_grid,
        mean_energies_eV=mean_grid,
        values_eV_m32=values,
        normalization_error_max=normalization_error,
        mean_energy_relative_error_max=mean_error,
        minimum_value=minimum,
    )


def _validate_comsol_function_eedf_metadata(
    grid: ComsolFunctionEedfGrid,
    metadata: Mapping[str, Any],
    contract: ComsolEedfImportContract,
) -> None:
    expected_columns = list(COMSOL_EEDF_COLUMNS)
    expected_shape = [
        len(grid.mean_energies_eV),
        len(grid.electron_energies_eV),
    ]
    comsol_import = metadata.get("comsol_import")
    if (
        metadata.get("representation") != contract.representation
        or metadata.get("artifact_role") != "canonical_comsol_function_eedf_input"
        or metadata.get("canonical_comsol_input") is not True
        or metadata.get("format") != contract.format
        or metadata.get("columns") != expected_columns
        or metadata.get("grid_shape") != expected_shape
        or metadata.get("mean_energy_grid_points") != expected_shape[0]
        or metadata.get("energy_grid_points") != expected_shape[1]
        or metadata.get("grid_axis_order") != ["electron_energy_eV", "mean_energy_eV"]
        or metadata.get("row_order") != contract.row_order
        or not isinstance(comsol_import, Mapping)
        or any(
            comsol_import.get(key) != expected
            for key, expected in contract.import_settings().items()
        )
    ):
        raise FunctionEedfError(
            "invalid canonical COMSOL Function-EEDF physical-table metadata"
        )
    reported = (
        ("projected_normalization_error_max", grid.normalization_error_max),
        (
            "projected_mean_energy_relative_error_max",
            grid.mean_energy_relative_error_max,
        ),
        ("projected_nonnegative_minimum", grid.minimum_value),
    )
    try:
        matches = all(
            math.isclose(float(metadata[name]), actual, rel_tol=1.0e-9, abs_tol=1.0e-13)
            for name, actual in reported
        )
    except (KeyError, TypeError, ValueError):
        matches = False
    if not matches:
        raise FunctionEedfError(
            "COMSOL Function-EEDF metadata audit disagrees with the physical table"
        )
