"""Long-form transport coefficient writer.

The summary CSV keeps legacy flux-compatible columns.  This writer emits a
separate long table so flux, bulk, and source-gradient coefficients can be
inspected without widening the summary schema.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from electron_swarm.core.config import OutputConfig
from electron_swarm.core.results import SwarmCaseResult, SwarmRunResult


def _append_value(
    rows: list[dict[str, object]],
    case: SwarmCaseResult,
    *,
    kind: str,
    quantity: str,
    value: float | None,
    unit: str,
    method: str,
    validity: str,
) -> None:
    if value is None:
        return
    rows.append(
        {
            "solver": case.solver,
            "case_id": case.case_id,
            "E_over_N_Td": case.e_over_n_Td,
            "kind": kind,
            "quantity": quantity,
            "value": float(value),
            "unit": unit,
            "method": method,
            "validity": validity,
        }
    )


def _transport_method_and_validity(case: SwarmCaseResult) -> tuple[str, str]:
    if case.transport is None or case.transport.metadata is None:
        return "", ""
    meta = case.transport.metadata
    return meta.coefficient_definition, meta.swarm_condition


def _append_drift_diffusion(
    rows: list[dict[str, object]],
    case: SwarmCaseResult,
    *,
    kind: str,
    transport: object,
    method: str,
    validity: str,
) -> None:
    _append_value(rows, case, kind=kind, quantity="drift_velocity", value=getattr(transport, "drift_velocity_m_s"), unit="m/s", method=method, validity=validity)
    _append_value(rows, case, kind=kind, quantity="mobility", value=getattr(transport, "mobility_m2_V_s"), unit="m2/(V s)", method=method, validity=validity)
    _append_value(rows, case, kind=kind, quantity="diffusion_longitudinal", value=getattr(transport, "diffusion_longitudinal_m2_s"), unit="m2/s", method=method, validity=validity)
    _append_value(rows, case, kind=kind, quantity="diffusion_transverse", value=getattr(transport, "diffusion_transverse_m2_s"), unit="m2/s", method=method, validity=validity)
    _append_value(rows, case, kind=kind, quantity="characteristic_energy_longitudinal", value=getattr(transport, "characteristic_energy_longitudinal_eV"), unit="eV", method=method, validity=validity)
    _append_value(rows, case, kind=kind, quantity="characteristic_energy_transverse", value=getattr(transport, "characteristic_energy_transverse_eV"), unit="eV", method=method, validity=validity)


def _transport_rows(case: SwarmCaseResult) -> list[dict[str, object]]:
    if case.transport is None:
        return []
    rows: list[dict[str, object]] = []
    method, validity = _transport_method_and_validity(case)
    _append_drift_diffusion(rows, case, kind="flux", transport=case.transport.flux, method=method, validity=validity)
    if case.transport.bulk is not None:
        _append_drift_diffusion(rows, case, kind="bulk", transport=case.transport.bulk, method=method, validity=validity)
    source = case.transport.source
    _append_value(rows, case, kind="source", quantity="ionization_frequency", value=source.ionization_frequency_s_inv, unit="1/s", method=method, validity=validity)
    _append_value(rows, case, kind="source", quantity="attachment_frequency", value=source.attachment_frequency_s_inv, unit="1/s", method=method, validity=validity)
    _append_value(rows, case, kind="source", quantity="effective_growth_frequency", value=source.effective_growth_frequency_s_inv, unit="1/s", method=method, validity=validity)
    _append_value(rows, case, kind="source", quantity="gradient_velocity", value=source.gradient_velocity_m_s, unit="m/s", method=method, validity=validity)
    _append_value(rows, case, kind="source", quantity="curvature_diffusion_longitudinal", value=source.curvature_diffusion_longitudinal_m2_s, unit="m2/s", method=method, validity=validity)
    _append_value(rows, case, kind="source", quantity="curvature_diffusion_transverse", value=source.curvature_diffusion_transverse_m2_s, unit="m2/s", method=method, validity=validity)
    return rows


def transport_frame(result: SwarmRunResult) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for case in result.cases:
        rows.extend(_transport_rows(case))
    return pd.DataFrame(rows)


def write_transport_outputs(result: SwarmRunResult, output: OutputConfig) -> Path | None:
    frame = transport_frame(result)
    if frame.empty:
        return None
    path = output.directory / f"{output.base_name}_transport.csv"
    frame.to_csv(path, index=False, float_format=output.float_format)
    return path
