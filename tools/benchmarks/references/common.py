"""Shared external-reference ingestion helpers.

External references are benchmark inputs, not product solver modes.  Readers
normalize all distributions to EEDF F(E) in 1/eV before comparison.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from hashlib import sha256
from io import StringIO
from typing import Any, Literal, cast

import numpy as np
import pandas as pd

from electron_swarm.core.results import RateResult, SwarmCaseResult
from electron_swarm.grids.energy import cell_edges_from_centers, reported_cell_edges
from electron_swarm.core.yaml_input import read_yaml_mapping
from tools.benchmarks.eedf_compare import (
    cell_widths_from_centers,
    EedfComparisonCase,
    compare_eedf_cases,
    normalize_eedf,
)


SCALAR_COLUMNS = {
    "mean_energy_eV",
    "drift_velocity_m_s",
    "mobility_m2_V_s",
    "diffusion_L_m2_s",
    "diffusion_T_m2_s",
    "net_ionization_frequency_s",
}
UNCERTAINTY_COLUMNS = {
    "mean_energy_eV": "mean_energy_eV_ci95",
    "drift_velocity_m_s": "drift_velocity_m_s_ci95",
    "mobility_m2_V_s": "mobility_m2_V_s_ci95",
    "diffusion_L_m2_s": "diffusion_L_m2_s_ci95",
    "diffusion_T_m2_s": "diffusion_T_m2_s_ci95",
    "net_ionization_frequency_s": "net_ionization_frequency_s_ci95",
}

ReferenceId = Literal["bolsig_plus", "mcig"]
ReferenceFormat = Literal[
    "electron_swarm_reference_csv",
    "bolsig_text",
    "mcig_csv",
]
ReferenceEedfConvention = Literal["eedf", "eepf"]
ReferenceUncertainty = Literal["unavailable", "reported"]
ReferenceAngularModel = Literal["isotropic", "mcig_default", "unknown"]


@dataclass(slots=True)
class ExternalReferenceConfig:
    id: ReferenceId
    path: Path
    format: ReferenceFormat = "electron_swarm_reference_csv"
    eedf_convention: ReferenceEedfConvention = "eedf"
    uncertainty: ReferenceUncertainty = "unavailable"
    angular_model: ReferenceAngularModel = "unknown"


@dataclass(slots=True)
class ReferenceCaseResult:
    reference_id: str
    case_id: str
    e_over_n_Td: float
    energy_eV: np.ndarray = field(repr=False)
    eedf_eV_inv: np.ndarray = field(repr=False)
    energy_widths_eV: np.ndarray | None = field(default=None, repr=False)
    rates: dict[str, float] = field(default_factory=dict)
    scalars: dict[str, float] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)


def _validate_literal(value: str, allowed: set[str], field_name: str) -> str:
    if value not in allowed:
        raise ValueError(f"{field_name} must be one of {sorted(allowed)}")
    return value


def parse_external_reference_configs(
    raw: dict[str, Any],
    *,
    base: Path,
) -> list[ExternalReferenceConfig]:
    refs_raw = raw.get("references", {}) or {}
    if not isinstance(refs_raw, dict):
        raise ValueError("references must be a mapping")
    if set(refs_raw) - {"external"}:
        raise ValueError("references accepts only the external list")
    external_raw = refs_raw.get("external", []) or []
    if not isinstance(external_raw, list):
        raise ValueError("references.external must be a list")
    configs: list[ExternalReferenceConfig] = []
    for index, item in enumerate(external_raw):
        if not isinstance(item, dict):
            raise ValueError(f"references.external[{index}] must be a mapping")
        unknown = set(item) - {
            "id",
            "path",
            "format",
            "eedf_convention",
            "uncertainty",
            "angular_model",
        }
        if unknown:
            raise ValueError(
                f"Unsupported references.external[{index}] fields: {sorted(unknown)}"
            )
        if "id" not in item:
            raise ValueError(f"references.external[{index}].id is required")
        if "path" not in item:
            raise ValueError(f"references.external[{index}].path is required")
        path = Path(str(item["path"]))
        if not path.is_absolute():
            path = (base / path).resolve()
        configs.append(
            ExternalReferenceConfig(
                id=cast(
                    ReferenceId,
                    _validate_literal(
                        str(item["id"]),
                        {"bolsig_plus", "mcig"},
                        f"references.external[{index}].id",
                    ),
                ),
                path=path,
                format=cast(
                    ReferenceFormat,
                    _validate_literal(
                        str(item.get("format", "electron_swarm_reference_csv")),
                        {"electron_swarm_reference_csv", "bolsig_text", "mcig_csv"},
                        f"references.external[{index}].format",
                    ),
                ),
                eedf_convention=cast(
                    ReferenceEedfConvention,
                    _validate_literal(
                        str(item.get("eedf_convention", "eedf")),
                        {"eedf", "eepf"},
                        f"references.external[{index}].eedf_convention",
                    ),
                ),
                uncertainty=cast(
                    ReferenceUncertainty,
                    _validate_literal(
                        str(item.get("uncertainty", "unavailable")),
                        {"unavailable", "reported"},
                        f"references.external[{index}].uncertainty",
                    ),
                ),
                angular_model=cast(
                    ReferenceAngularModel,
                    _validate_literal(
                        str(item.get("angular_model", "unknown")),
                        {"isotropic", "mcig_default", "unknown"},
                        f"references.external[{index}].angular_model",
                    ),
                ),
            )
        )
    return configs


def load_external_reference_configs(path: Path) -> list[ExternalReferenceConfig]:
    data = read_yaml_mapping(path)
    return parse_external_reference_configs(data, base=Path(path).resolve().parent)


def eepf_to_eedf(energy_eV: np.ndarray, eepf_eV_m32: np.ndarray) -> np.ndarray:
    energy = np.asarray(energy_eV, dtype=float)
    eepf = np.asarray(eepf_eV_m32, dtype=float)
    if np.any(energy < 0.0):
        raise ValueError("reference energy_eV must be nonnegative")
    return eepf * np.sqrt(np.maximum(energy, 0.0))


def _validate_energy(energy_eV: np.ndarray) -> None:
    if (
        energy_eV.ndim != 1
        or len(energy_eV) < 2
        or not np.all(np.isfinite(energy_eV))
        or np.any(energy_eV < 0.0)
        or np.any(np.diff(energy_eV) <= 0.0)
    ):
        raise ValueError(
            "reference energy_eV must be finite, nonnegative, and strictly increasing"
        )


def _normal_reference_eedf(
    energy_eV: np.ndarray,
    values: np.ndarray,
    *,
    convention: str,
    widths: np.ndarray | None = None,
) -> tuple[np.ndarray, float]:
    _validate_energy(energy_eV)
    raw = eepf_to_eedf(energy_eV, values) if convention == "eepf" else values
    if np.any(raw < 0):
        raise ValueError("reference distribution must be nonnegative")
    return normalize_eedf(
        energy_eV,
        raw,
        widths if widths is not None else cell_widths_from_centers(energy_eV),
    )


def _first_float(group: pd.DataFrame, column: str) -> float | None:
    if column not in group:
        return None
    supplied = group[column].replace("", np.nan).dropna()
    values = pd.to_numeric(supplied, errors="raise")
    if values.empty:
        return None
    numbers = values.to_numpy(dtype=float)
    if not np.isfinite(numbers).all() or not np.allclose(
        numbers, numbers[0], rtol=1e-12, atol=0
    ):
        raise ValueError(
            f"reference scalar {column} must be finite and consistent within a case"
        )
    return float(numbers[0])


def _result_from_group(
    reference_id: str,
    group: pd.DataFrame,
    *,
    convention: str,
    metadata: dict[str, Any],
) -> ReferenceCaseResult:
    if "energy_eV" not in group:
        raise ValueError("reference CSV requires energy_eV")
    eedf_columns = [name for name in ("eedf_eV_inv", "eepf_eV_m32") if name in group]
    if len(eedf_columns) != 1:
        raise ValueError(
            "reference CSV requires exactly one of eedf_eV_inv or eepf_eV_m32"
        )
    eedf_column = eedf_columns[0]
    effective_convention = "eepf" if eedf_column == "eepf_eV_m32" else "eedf"
    if convention != effective_convention:
        raise ValueError(
            f"declared {convention} convention contradicts reference column {eedf_column}"
        )

    ordered = group.sort_values("energy_eV")
    energy = pd.to_numeric(ordered["energy_eV"], errors="raise").to_numpy(dtype=float)
    values = pd.to_numeric(ordered[eedf_column], errors="raise").to_numpy(dtype=float)
    widths = None
    if "energy_width_eV" in ordered:
        widths = pd.to_numeric(ordered["energy_width_eV"], errors="raise").to_numpy(
            dtype=float
        )
        reported_cell_edges(energy, widths)
    eedf, original_norm = _normal_reference_eedf(
        energy,
        values,
        convention=effective_convention,
        widths=widths,
    )
    integration_widths = cell_widths_from_centers(energy) if widths is None else widths
    derived_mean = float(np.sum(energy * eedf * integration_widths))

    scalars: dict[str, float] = {}
    for column in SCALAR_COLUMNS:
        value = _first_float(ordered, column)
        if value is not None:
            scalars[column] = value
    mean_reported = "mean_energy_eV" in scalars
    scalars.setdefault("mean_energy_eV", derived_mean)

    rates: dict[str, float] = {}
    for column in ordered.columns:
        if column.startswith("rate__"):
            value = _first_float(ordered, column)
            if value is not None:
                if value < 0:
                    raise ValueError(f"reference rate {column} must be nonnegative")
                rates[column.removeprefix("rate__")] = value

    reported_mean = scalars.get("mean_energy_eV", derived_mean)
    result_metadata = {
        **metadata,
        "eedf_convention_input": eedf_column,
        "distribution_geometry": "reported_cells"
        if widths is not None
        else "inferred_midpoint_cells",
        "mean_energy_source": "reported" if mean_reported else "derived_eedf",
        "eedf_normalization_before": original_norm,
        "eedf_mean_energy_eV": derived_mean,
        "reported_mean_energy_difference_eV": abs(reported_mean - derived_mean),
    }
    scalar_ci95: dict[str, float] = {}
    for scalar, column in UNCERTAINTY_COLUMNS.items():
        value = _first_float(ordered, column)
        if value is not None and value >= 0.0:
            scalar_ci95[scalar] = value
        elif value is not None:
            raise ValueError(f"reference uncertainty {column} must be nonnegative")
    if scalar_ci95:
        result_metadata["scalar_ci95"] = scalar_ci95
    return ReferenceCaseResult(
        reference_id=reference_id,
        case_id=str(ordered["case_id"].iloc[0]),
        e_over_n_Td=float(ordered["E_over_N_Td"].iloc[0]),
        energy_eV=energy,
        eedf_eV_inv=eedf,
        energy_widths_eV=widths,
        rates=rates,
        scalars=scalars,
        metadata=result_metadata,
    )


def load_canonical_reference_csv(
    path: Path,
    *,
    reference_id: str,
    convention: str,
    metadata: dict[str, Any] | None = None,
) -> list[ReferenceCaseResult]:
    if not path.exists():
        raise FileNotFoundError(f"external reference output not found: {path}")
    content = path.read_bytes()
    frame = pd.read_csv(
        StringIO(content.decode("utf-8-sig")),
        dtype={"case_id": str},
        keep_default_na=False,
    )
    return reference_cases_from_frame(
        frame,
        reference_id=reference_id,
        convention=convention,
        metadata={
            **(metadata or {}),
            "source_path": str(path.resolve()),
            "source_sha256": sha256(content).hexdigest(),
        },
    )


def reference_cases_from_frame(
    frame: pd.DataFrame,
    *,
    reference_id: str,
    convention: str,
    metadata: dict[str, Any] | None = None,
) -> list[ReferenceCaseResult]:
    missing = {"case_id", "E_over_N_Td", "energy_eV"} - set(frame.columns)
    if missing:
        raise ValueError(f"reference table missing required columns: {sorted(missing)}")
    frame = frame.copy()
    if (
        frame.empty
        or frame["case_id"].isna().any()
        or frame["case_id"].astype(str).str.strip().eq("").any()
    ):
        raise ValueError("reference table requires nonempty case identifiers and rows")
    frame["E_over_N_Td"] = pd.to_numeric(frame["E_over_N_Td"], errors="raise")
    fields = frame["E_over_N_Td"].to_numpy(dtype=float)
    if not np.isfinite(fields).all() or np.any(fields <= 0):
        raise ValueError("reference E/N must be finite and positive")
    return [
        _result_from_group(
            reference_id,
            group,
            convention=convention,
            metadata=metadata or {},
        )
        for _, group in frame.groupby(["case_id", "E_over_N_Td"], sort=False)
    ]


def reference_comparison_case(reference: ReferenceCaseResult) -> EedfComparisonCase:
    energy = reference.energy_eV
    eedf = reference.eedf_eV_inv
    scalar = reference.scalars
    if reference.energy_widths_eV is None:
        edges, widths = cell_edges_from_centers(energy)
    else:
        widths = reference.energy_widths_eV
        edges = reported_cell_edges(energy, widths)
    mean_energy = scalar.get("mean_energy_eV")
    if mean_energy is None:
        mean_energy = reference.metadata.get("eedf_mean_energy_eV")
    if mean_energy is None:
        mean_energy = float(
            np.sum(energy * eedf * widths) / max(float(np.sum(eedf * widths)), 1.0e-300)
        )
    rates = [
        RateResult(
            solver=f"reference:{reference.reference_id}",
            case_id=reference.case_id,
            e_over_n_Td=reference.e_over_n_Td,
            species="reference",
            process=name,
            process_type=name
            if name in {"ionization", "excitation", "attachment"}
            else "reference",
            threshold_eV=None,
            rate_coefficient_m3_s=value,
        )
        for name, value in reference.rates.items()
    ]
    return EedfComparisonCase(
        solver=f"reference:{reference.reference_id}",
        case_id=reference.case_id,
        e_over_n_Td=reference.e_over_n_Td,
        mean_energy_eV=float(mean_energy),
        drift_velocity_m_s=scalar.get("drift_velocity_m_s"),
        mobility_m2_V_s=scalar.get("mobility_m2_V_s"),
        diffusion_L_m2_s=scalar.get("diffusion_L_m2_s"),
        diffusion_T_m2_s=scalar.get("diffusion_T_m2_s"),
        energy_eV=energy,
        eedf=eedf,
        energy_widths_eV=widths,
        energy_edges_eV=edges,
        rates=rates,
        metadata={
            **reference.metadata,
            "solver_method": reference.reference_id,
            "physics_level": "external_reference",
            "reference_id": reference.reference_id,
            "angular_model": reference.metadata.get("angular_model", "unknown"),
        },
    )


def reference_comparison_metrics(
    reference: ReferenceCaseResult,
    candidate: SwarmCaseResult | EedfComparisonCase,
) -> dict[str, float]:
    reference_case = reference_comparison_case(reference)
    metrics = dict(compare_eedf_cases(reference_case, candidate).metrics)
    stale = candidate.metadata.get("electron_electron_transport_stale") is True
    if stale:
        metrics.pop("drift_velocity_relative_difference", None)
    for attr, name in (
        ("mobility_m2_V_s", "mobility_relative_difference"),
        ("diffusion_L_m2_s", "diffusion_L_relative_difference"),
        ("diffusion_T_m2_s", "diffusion_T_relative_difference"),
    ):
        reference_value = getattr(reference_case, attr)
        candidate_value = getattr(candidate, attr)
        if stale or reference_value is None or candidate_value is None:
            continue
        if abs(float(reference_value)) > 0.0:
            metrics[name] = abs(float(candidate_value) - float(reference_value)) / max(
                abs(float(reference_value)),
                1.0e-300,
            )
        elif candidate_value == 0:
            metrics[name] = 0.0
    # Reference process names and explicit type aggregates own rate identity.
    # Never compare an unmatched rate against a manufactured zero.
    metrics.pop("major_rate_relative_difference", None)
    metrics.pop("rate_weighted_eedf_error", None)
    differences = {}
    for name, value in reference.rates.items():
        aggregate = name in {"ionization", "excitation", "attachment"}
        matches = [
            rate
            for rate in candidate.rates
            if (rate.process_type == name if aggregate else rate.process == name)
        ]
        if not matches:
            continue
        if not aggregate and len(matches) != 1:
            raise ValueError(f"ambiguous external reference rate identity: {name}")
        actual = sum(rate.mixture_weighted_rate_m3_s for rate in matches)
        if value == 0 and actual != 0:
            continue
        difference = abs(actual - value) / abs(value) if value else 0.0
        metrics[f"{name}_rate_relative_difference"] = difference
        differences[name] = difference
    scale = max(reference.rates.values(), default=0.0)
    important = [
        name for name, value in reference.rates.items() if value >= 0.01 * scale
    ]
    if important and all(name in differences for name in important):
        metrics["major_rate_relative_difference"] = max(
            differences[name] for name in important
        )
        metrics["rate_weighted_eedf_error"] = metrics["major_rate_relative_difference"]
    return metrics


def reference_comparison_status(
    metrics, thresholds, *, confidence_status="unknown", angular_model_status="match"
) -> dict:
    """Coverage and evidence gates, with no missing-as-zero or CI override."""
    missing = [
        name
        for name in thresholds
        if name not in metrics or not np.isfinite(metrics[name])
    ]
    failed = [
        name
        for name, bound in thresholds.items()
        if name not in missing and metrics[name] > bound
    ]
    if (
        not missing
        and confidence_status == "pass_within_mc_uncertainty"
        and set(failed) <= {"mean_energy_relative_difference"}
    ):
        status = "PASS_WITHIN_MC_UNCERTAINTY"
    else:
        status = "FAIL" if failed else "PARTIAL" if missing else "PASS"
    if angular_model_status != "match" and status != "FAIL":
        status = "DEGRADED"
    return {
        "status": status,
        "unavailable_metrics": ";".join(missing),
        "failed_metrics": ";".join(failed),
    }


def combined_reference_status(statuses) -> str:
    values = set(statuses)
    if not values or values <= {"missing", "UNAVAILABLE"}:
        return "UNAVAILABLE"
    for status in ("FAIL", "DEGRADED"):
        if status in values:
            return status
    if values & {"PARTIAL", "missing", "UNAVAILABLE"}:
        return "PARTIAL"
    return (
        "PASS_WITHIN_MC_UNCERTAINTY"
        if "PASS_WITHIN_MC_UNCERTAINTY" in values
        else "PASS"
    )


def load_reference_cases(config: ExternalReferenceConfig) -> list[ReferenceCaseResult]:
    if config.id == "bolsig_plus":
        from .bolsig import load_bolsig_reference

        return load_bolsig_reference(config)
    if config.id == "mcig":
        from .mcig import load_mcig_reference

        return load_mcig_reference(config)
    raise ValueError(f"unsupported reference id: {config.id}")
