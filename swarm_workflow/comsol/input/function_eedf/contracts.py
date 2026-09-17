"""Typed contracts for canonical COMSOL Function-EEDF representations."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np


# The source-bin quadrature may differ slightly from a solver's independently
# integrated mean energy. Larger discrepancies would make the exponential
# moment match conceal a materially inconsistent source distribution.
SOURCE_MEAN_RELATIVE_ERROR_LIMIT = 1.0e-3

COMSOL_EEDF_COLUMNS = ("electron_energy_eV", "mean_energy_eV", "eepf_eV_m32")


class FunctionEedfError(ValueError):
    """Raised when an EEDF family cannot satisfy its kinetic constraints."""


@dataclass(frozen=True, slots=True)
class ComsolEedfImportContract:
    """The canonical solver-independent COMSOL Function-EEDF import."""

    structure: str = "spreadsheet"

    def __post_init__(self) -> None:
        if self.structure != "spreadsheet":
            raise FunctionEedfError(
                "canonical Function-EEDF input requires spreadsheet serialization"
            )

    @classmethod
    def from_path(cls, path: str | Path) -> ComsolEedfImportContract:
        try:
            with Path(path).open("r", encoding="utf-8-sig", newline="") as stream:
                first = next((line.strip() for line in stream if line.strip()), "")
        except OSError as exc:
            raise FunctionEedfError(
                f"cannot read COMSOL Function-EEDF: {path}"
            ) from exc
        if tuple(next(csv.reader([first]))) == COMSOL_EEDF_COLUMNS:
            return cls("spreadsheet")
        raise FunctionEedfError(
            "COMSOL Function-EEDF requires canonical spreadsheet columns"
        )

    @property
    def format(self) -> str:
        return "csv"

    @property
    def representation(self) -> str:
        return "physical_2d_adaptive_moment_rate_projected"

    @property
    def row_order(self) -> str:
        return "mean_energy_major_then_electron_energy"

    def import_settings(self) -> dict[str, Any]:
        return {
            "source": "file",
            "struct": self.structure,
            "nargs": 2,
            "argunit": "eV,eV",
            "fununit": "1",
            "interp": "linear",
            "extrap": "const",
            "funcnametable_position": "1",
            "scaledata": "auto",
        }


@dataclass(frozen=True, slots=True)
class C1FunctionEedf:
    """Dimensionless EEDF shapes with invariant zeroth and first moments."""

    mean_energies_eV: np.ndarray
    dimensionless_energy: np.ndarray
    shapes: tuple[np.ndarray, ...]
    source_energy_max_eV: float
    source_mean_relative_error_max: float
    normalization_error_max: float
    unit_mean_error_max: float
    minimum_value: float

    @property
    def shape_columns(self) -> tuple[str, ...]:
        return tuple(f"shape_{index:04d}" for index in range(len(self.shapes)))


@dataclass(frozen=True, slots=True)
class ComsolFunctionEedfGrid:
    """Physical 2D EEDF grid that preserves moments under linear interpolation."""

    electron_energies_eV: np.ndarray
    mean_energies_eV: np.ndarray
    values_eV_m32: np.ndarray
    normalization_error_max: float
    mean_energy_relative_error_max: float
    minimum_value: float
    shape_total_variation_error_max: float = 0.0
    rate_scaled_error_max: float | None = None
    rate_kernel_count: int = 0
    mean_axis_shape_total_variation_error_max: float = 0.0
    mean_axis_rate_scaled_error_max: float | None = None


@dataclass(frozen=True, slots=True)
class CollisionRateKernel:
    """One physical collision kernel used to control EEDF projection error."""

    species: str
    process: str
    process_type: str
    threshold_eV: float
    electron_energies_eV: np.ndarray
    cross_sections_m2: np.ndarray
    high_energy_extrapolation: Literal["zero", "hold", "error"] = "zero"

    @property
    def identity(self) -> str:
        return f"{self.species}:{self.process}:{self.process_type}"
