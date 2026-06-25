"""BOLSIG+ reference output readers."""

from __future__ import annotations

from io import StringIO
from pathlib import Path

import pandas as pd

from .common import (
    ExternalReferenceConfig,
    ReferenceCaseResult,
    load_canonical_reference_csv,
    reference_cases_from_frame,
)


def _load_bolsig_text(path: Path, *, convention: str) -> list[ReferenceCaseResult]:
    if not path.exists():
        raise FileNotFoundError(f"BOLSIG+ reference output not found: {path}")
    metadata: dict[str, str] = {}
    table_lines: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("#"):
            line = line[1:].strip()
        if ":" in line and not table_lines:
            key, value = line.split(":", 1)
            metadata[key.strip()] = value.strip()
            continue
        table_lines.append(line)
    if not table_lines:
        raise ValueError("BOLSIG+ text output did not contain an EEDF table")

    header = table_lines[0].replace(",", " ").split()
    data_lines = [" ".join(header), *table_lines[1:]]
    frame = pd.read_csv(StringIO("\n".join(data_lines)), sep=r"\s+", engine="python")
    if "E_over_N_Td" not in frame:
        frame["E_over_N_Td"] = float(metadata.get("E_over_N_Td", metadata.get("E/N_Td", 0.0)))
    if "case_id" not in frame:
        frame["case_id"] = metadata.get("case_id", "external_reference")
    for key in ("mean_energy_eV", "drift_velocity_m_s", "mobility_m2_V_s"):
        if key in metadata and key not in frame:
            frame[key] = float(metadata[key])

    return reference_cases_from_frame(
        frame,
        reference_id="bolsig_plus",
        convention=convention,
        metadata={
            "reference_format": "bolsig_text",
            "reference_role": "two_term_reference",
            "angular_model": "unknown",
        },
    )


def load_bolsig_reference(config: ExternalReferenceConfig) -> list[ReferenceCaseResult]:
    metadata = {
        "reference_format": config.format,
        "reference_role": "two_term_reference",
        "angular_model": config.angular_model,
    }
    if config.format == "electron_swarm_reference_csv":
        return load_canonical_reference_csv(
            config.path,
            reference_id="bolsig_plus",
            convention=config.eedf_convention,
            metadata=metadata,
        )
    if config.format == "bolsig_text":
        return _load_bolsig_text(config.path, convention=config.eedf_convention)
    raise ValueError(f"unsupported BOLSIG+ reference format: {config.format}")
