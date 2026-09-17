"""Gas-condition and cross-section input parsers."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal, cast

from electron_swarm.core.config import (
    ConditionsConfig,
    CrossSectionFileConfig,
    CrossSectionsConfig,
    GasComponent,
)
from electron_swarm.core.config_sections.common import (
    as_mapping_section,
    as_path,
    reject_unknown_fields,
    strict_float,
    string_value,
    validate_literal,
)


def parse_conditions(raw: dict[str, Any]) -> ConditionsConfig:
    cond_raw = as_mapping_section(raw, "conditions", "conditions")
    reject_unknown_fields(
        cond_raw,
        {
            "gas_temperature_K",
            "pressure_Pa",
            "gas_number_density_m3",
            "length_scale_m",
            "gas_mixture",
            "species",
            "mass_amu",
        },
        "conditions",
    )
    gas_mixture_raw = cond_raw.get("gas_mixture", [])
    if gas_mixture_raw is None:
        gas_mixture_raw = []
    if not isinstance(gas_mixture_raw, list):
        raise ValueError("conditions.gas_mixture must be a list")
    for index, item in enumerate(gas_mixture_raw):
        if not isinstance(item, dict):
            raise ValueError("conditions.gas_mixture entries must be mappings")
        reject_unknown_fields(
            item,
            {"species", "fraction", "mass_amu"},
            f"conditions.gas_mixture[{index}]",
        )
        if "species" not in item:
            raise ValueError(f"conditions.gas_mixture[{index}].species is required")

    if "gas_mixture" not in cond_raw:
        gas_mixture = [
            GasComponent(
                species=string_value(
                    cond_raw.get("species", "Ar"),
                    "conditions.species",
                ),
                fraction=1.0,
                mass_amu=strict_float(
                    cond_raw.get("mass_amu", 39.948),
                    "conditions.mass_amu",
                ),
            )
        ]
    else:
        gas_mixture = [
            GasComponent(
                species=string_value(
                    item["species"],
                    f"conditions.gas_mixture[{index}].species",
                ),
                fraction=strict_float(
                    item.get("fraction", 1.0),
                    f"conditions.gas_mixture[{index}].fraction",
                ),
                mass_amu=strict_float(
                    item.get("mass_amu", 39.948),
                    f"conditions.gas_mixture[{index}].mass_amu",
                ),
            )
            for index, item in enumerate(gas_mixture_raw)
        ]

    return ConditionsConfig(
        gas_temperature_K=strict_float(
            cond_raw.get("gas_temperature_K", 300.0),
            "conditions.gas_temperature_K",
        ),
        pressure_Pa=(
            strict_float(cond_raw["pressure_Pa"], "conditions.pressure_Pa")
            if cond_raw.get("pressure_Pa") is not None
            else None
        ),
        gas_number_density_m3=(
            strict_float(
                cond_raw["gas_number_density_m3"],
                "conditions.gas_number_density_m3",
            )
            if cond_raw.get("gas_number_density_m3") is not None
            else None
        ),
        length_scale_m=(
            strict_float(cond_raw["length_scale_m"], "conditions.length_scale_m")
            if cond_raw.get("length_scale_m") is not None
            else None
        ),
        gas_mixture=gas_mixture,
    )


def parse_cross_sections(raw: dict[str, Any], base: Path) -> CrossSectionsConfig:
    xs_raw = as_mapping_section(raw, "cross_sections", "cross_sections")
    reject_unknown_fields(
        xs_raw,
        {"format", "files", "high_energy_extrapolation"},
        "cross_sections",
    )
    xs_format = string_value(xs_raw.get("format", "csv"), "cross_sections.format")
    files_raw = xs_raw.get("files", [])
    if not isinstance(files_raw, list):
        raise ValueError("cross_sections.files must be a list")

    xs_files: list[CrossSectionFileConfig] = []
    for index, item in enumerate(files_raw):
        if not isinstance(item, dict):
            raise ValueError("cross_sections.files entries must be mappings")
        reject_unknown_fields(
            item,
            {"path", "species", "format"},
            f"cross_sections.files[{index}]",
        )
        xs_files.append(
            CrossSectionFileConfig(
                path=as_path(item["path"], base) or Path(item["path"]),
                species=item.get("species"),
                format=item.get("format", xs_format),
            )
        )
    return CrossSectionsConfig(
        format=xs_format,
        files=xs_files,
        high_energy_extrapolation=cast(
            Literal["zero", "hold", "error"],
            validate_literal(
                string_value(
                    xs_raw.get("high_energy_extrapolation", "zero"),
                    "cross_sections.high_energy_extrapolation",
                ),
                {"zero", "hold", "error"},
                "cross_sections.high_energy_extrapolation",
            ),
        ),
    )
