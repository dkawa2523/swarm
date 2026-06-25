"""MCIG reference output readers."""

from __future__ import annotations

import pandas as pd

from .common import (
    ExternalReferenceConfig,
    ReferenceCaseResult,
    load_canonical_reference_csv,
    reference_cases_from_frame,
)


def _load_mcig_csv(config: ExternalReferenceConfig) -> list[ReferenceCaseResult]:
    if not config.path.exists():
        raise FileNotFoundError(f"MCIG reference output not found: {config.path}")
    frame = pd.read_csv(config.path)
    rename = {
        "E_eV": "energy_eV",
        "energy": "energy_eV",
        "E/N_Td": "E_over_N_Td",
        "eedf": "eedf_eV_inv",
        "pdf_eV_inv": "eedf_eV_inv",
        "eepf": "eepf_eV_m32",
    }
    frame = frame.rename(columns={key: value for key, value in rename.items() if key in frame})
    if "case_id" not in frame:
        frame["case_id"] = "external_reference"
    return reference_cases_from_frame(
        frame,
        reference_id="mcig",
        convention=config.eedf_convention,
        metadata={
            "reference_format": "mcig_csv",
            "reference_role": "monte_carlo_reference",
            "uncertainty_unavailable": config.uncertainty == "unavailable",
            "angular_model": config.angular_model,
        },
    )


def load_mcig_reference(config: ExternalReferenceConfig) -> list[ReferenceCaseResult]:
    metadata = {
        "reference_format": config.format,
        "reference_role": "monte_carlo_reference",
        "uncertainty_unavailable": config.uncertainty == "unavailable",
        "angular_model": config.angular_model,
    }
    if config.format == "electron_swarm_reference_csv":
        return load_canonical_reference_csv(
            config.path,
            reference_id="mcig",
            convention=config.eedf_convention,
            metadata=metadata,
        )
    if config.format == "mcig_csv":
        return _load_mcig_csv(config)
    raise ValueError(f"unsupported MCIG reference format: {config.format}")
