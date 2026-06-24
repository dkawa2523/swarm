"""Triage internal-MC EEDF differences without fitting external references.

This development tool compares raw EEDF probability mass, quantiles, and
selected internal-MC physics metadata.  It does not treat two-term, MCIG, or any
other MC output as a truth source and it does not smooth or retune product
results.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from electron_swarm import load_config, run
from electron_swarm.core.config import SwarmConfig
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.core.numerics import widths_from_centers
from electron_swarm.diagnostics.eedf_compare import compare_eedf_cases, normalize_eedf
from electron_swarm.references.common import (
    reference_cases_from_frame,
    reference_to_swarm_case,
)
from tools.benchmark_common import (
    case_audit_metadata,
    metadata_value,
    variant_config,
    write_csv,
)


DEFAULT_BANDS = (
    (0.0, 5.0),
    (5.0, 10.0),
    (10.0, 20.0),
    (20.0, 40.0),
    (40.0, 80.0),
    (80.0, float("inf")),
)
TAIL_THRESHOLDS = (20.0, 40.0, 80.0)

METRIC_FIELDS = [
    "source",
    "case_id",
    "E_over_N_Td",
    "mean_energy_eV",
    "E50_eV",
    "E90_eV",
    "E99_eV",
    "tail_survival_gt_20_eV",
    "tail_survival_gt_40_eV",
    "tail_survival_gt_80_eV",
    "transport_definition",
    "mc_population_model",
    "swarm_population_treatment",
    "eedf_estimator",
    "nonconservative_growth_treatment",
    "ionization_event_fraction",
    "mc_physical_branching_gap_eV",
    "mc_secondary_electron_count",
    "mc_branching_resample_count",
    "mc_population_total_weight_final",
    "mc_population_log_growth_estimate_s_inv",
    "mc_population_weight_cv",
    "mc_population_resampling_energy_adjustment_eV",
    "mc_energy_balance_status",
    "mc_tail_uncertainty_status",
    "mc_energy_samples_above_xs_max_fraction",
]

BAND_FIELDS = [
    "source",
    "case_id",
    "E_over_N_Td",
    "band",
    "energy_left_eV",
    "energy_right_eV",
    "probability_mass",
    "mean_energy_contribution_eV",
    "survival_probability_at_left",
    "sample_count",
    "effective_sample_count",
    "relative_standard_error",
]

PAIR_FIELDS = [
    "E_over_N_Td",
    "reference_source",
    "candidate_source",
    "eedf_relative_l1",
    "log_tail_error",
    "tail_probability_difference",
    "mean_energy_reference_eV",
    "mean_energy_candidate_eV",
    "mean_energy_signed_difference_eV",
    "mean_energy_signed_relative_difference",
    "E50_difference_eV",
    "E90_difference_eV",
    "E99_difference_eV",
]

ASSESSMENT_FIELDS = [
    "E_over_N_Td",
    "reference_source",
    "candidate_source",
    "classification",
    "dominant_band",
    "dominant_band_mean_energy_difference_eV",
    "evidence",
]


@dataclass(slots=True)
class Distribution:
    source: str
    case_id: str
    e_over_n_Td: float
    energy_eV: np.ndarray = field(repr=False)
    eedf_eV_inv: np.ndarray = field(repr=False)
    widths_eV: np.ndarray = field(repr=False)
    metadata: dict[str, Any] = field(default_factory=dict)
    sample_count: np.ndarray | None = field(default=None, repr=False)
    effective_sample_count: np.ndarray | None = field(default=None, repr=False)

    @property
    def mean_energy_eV(self) -> float:
        return float(np.sum(self.energy_eV * self.eedf_eV_inv * self.widths_eV))


def _finite_widths(energy: np.ndarray, widths: np.ndarray | None) -> np.ndarray:
    if widths is not None:
        arr = np.asarray(widths, dtype=float)
        if arr.shape == energy.shape and np.all(np.isfinite(arr)) and np.all(arr > 0.0):
            return arr
    return widths_from_centers(energy)


def _normalized_distribution(
    source: str,
    case_id: str,
    e_over_n_Td: float,
    energy_eV: np.ndarray,
    eedf_eV_inv: np.ndarray,
    widths_eV: np.ndarray | None,
    *,
    metadata: dict[str, Any] | None = None,
    sample_count: np.ndarray | None = None,
    effective_sample_count: np.ndarray | None = None,
) -> Distribution:
    energy = np.asarray(energy_eV, dtype=float)
    widths = _finite_widths(energy, widths_eV)
    eedf, _ = normalize_eedf(energy, np.asarray(eedf_eV_inv, dtype=float), widths)
    return Distribution(
        source=source,
        case_id=case_id,
        e_over_n_Td=float(e_over_n_Td),
        energy_eV=energy,
        eedf_eV_inv=eedf,
        widths_eV=widths,
        metadata=dict(metadata or {}),
        sample_count=None if sample_count is None else np.asarray(sample_count, dtype=float),
        effective_sample_count=(
            None
            if effective_sample_count is None
            else np.asarray(effective_sample_count, dtype=float)
        ),
    )


def _source_for_case(case: SwarmCaseResult) -> str:
    metadata = case_audit_metadata(case)
    if (
        case.solver == "monte_carlo"
        and metadata_value(metadata, "monte_carlo_backend") == "internal"
    ):
        return "internal_mc"
    return case.solver.removeprefix("reference:")


def _case_distribution_metadata(case: SwarmCaseResult) -> dict[str, Any]:
    return case_audit_metadata(case)


def _distribution_from_case(
    case: SwarmCaseResult,
    *,
    source: str | None = None,
) -> Distribution:
    return _normalized_distribution(
        source or _source_for_case(case),
        case.case_id,
        case.e_over_n_Td,
        case.energy_eV,
        case.eedf,
        case.energy_widths_eV,
        metadata=_case_distribution_metadata(case),
        sample_count=case.eedf_counts,
        effective_sample_count=case.eedf_effective_counts,
    )


def _float_or_blank(value: Any) -> float | str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return ""
    return number if np.isfinite(number) else ""


def _float_metadata(metadata: dict[str, Any], key: str) -> float:
    value = _float_or_blank(metadata_value(metadata, key))
    return float(value) if value != "" else 0.0


def _summary_metadata(path: Path | None) -> dict[tuple[str, float], dict[str, Any]]:
    if path is None:
        return {}
    frame = pd.read_csv(path)
    rows: dict[tuple[str, float], dict[str, Any]] = {}
    for record in frame.to_dict("records"):
        rows[(str(record["case_id"]), float(record["E_over_N_Td"]))] = dict(record)
    return rows


def _distributions_from_internal_tail(
    tail_path: Path,
    *,
    summary_path: Path | None = None,
) -> list[Distribution]:
    frame = pd.read_csv(tail_path)
    eedf_column = "eedf" if "eedf" in frame.columns else "eedf_eV_inv"
    if eedf_column not in frame.columns:
        raise ValueError("internal tail CSV requires eedf or eedf_eV_inv")
    summary = _summary_metadata(summary_path)
    distributions: list[Distribution] = []
    for (case_id, eover), group in frame.groupby(["case_id", "E_over_N_Td"], sort=False):
        ordered = group.sort_values("energy_eV")
        metadata = summary.get((str(case_id), float(eover)), {})
        widths = (
            pd.to_numeric(ordered["energy_width_eV"], errors="coerce").to_numpy(float)
            if "energy_width_eV" in ordered
            else None
        )
        distributions.append(
            _normalized_distribution(
                "internal_mc",
                str(case_id),
                float(eover),
                pd.to_numeric(ordered["energy_eV"], errors="raise").to_numpy(float),
                pd.to_numeric(ordered[eedf_column], errors="raise").to_numpy(float),
                widths,
                metadata=metadata,
                sample_count=(
                    pd.to_numeric(ordered["sample_count"], errors="coerce").to_numpy(float)
                    if "sample_count" in ordered
                    else None
                ),
                effective_sample_count=(
                    pd.to_numeric(
                        ordered["effective_sample_count"], errors="coerce"
                    ).to_numpy(float)
                    if "effective_sample_count" in ordered
                    else None
                ),
            )
        )
    return distributions


def _reference_distributions(path: Path, reference_id: str) -> list[Distribution]:
    cases = reference_cases_from_frame(
        pd.read_csv(path),
        reference_id=reference_id,
        convention="eedf",
    )
    return [
        _distribution_from_case(reference_to_swarm_case(case), source=reference_id)
        for case in cases
    ]


def _variant_config(cfg: SwarmConfig, solver: str) -> SwarmConfig:
    return variant_config(cfg, solver)


def _run_config_distributions(
    config_path: Path,
    *,
    include_internal: bool,
    include_two_term: bool,
) -> list[Distribution]:
    cfg = load_config(config_path)
    distributions: list[Distribution] = []
    for solver, include in (
        ("two_term", include_two_term),
        ("monte_carlo", include_internal),
    ):
        if not include:
            continue
        result = run(
            _variant_config(cfg, solver),
            write=False,
            collect_diagnostics=True,
        )
        distributions.extend(_distribution_from_case(case) for case in result.cases)
    return distributions


def _bin_edges(dist: Distribution) -> tuple[np.ndarray, np.ndarray]:
    left = dist.energy_eV - 0.5 * dist.widths_eV
    right = dist.energy_eV + 0.5 * dist.widths_eV
    if len(left):
        left[0] = max(0.0, float(left[0]))
    return left, right


def _band_integrals(
    dist: Distribution,
    left_eV: float,
    right_eV: float,
) -> tuple[float, float]:
    left_edges, right_edges = _bin_edges(dist)
    overlap_left = np.maximum(left_edges, left_eV)
    overlap_right = np.minimum(right_edges, right_eV)
    overlap = np.clip(overlap_right - overlap_left, 0.0, None)
    mass = float(np.sum(dist.eedf_eV_inv * overlap))
    mean_contribution = float(
        np.sum(0.5 * dist.eedf_eV_inv * (overlap_right**2 - overlap_left**2) * (overlap > 0.0))
    )
    return mass, mean_contribution


def _tail_probability(dist: Distribution, threshold_eV: float) -> float:
    mass, _ = _band_integrals(dist, threshold_eV, float("inf"))
    return mass


def _energy_quantile(dist: Distribution, quantile: float) -> float:
    left, right = _bin_edges(dist)
    mass = dist.eedf_eV_inv * np.maximum(right - left, 0.0)
    total = max(float(np.sum(mass)), 1.0e-300)
    target = quantile * total
    cumulative = 0.0
    for bin_left, bin_right, density, bin_mass in zip(
        left, right, dist.eedf_eV_inv, mass, strict=False
    ):
        if cumulative + bin_mass >= target:
            remaining = target - cumulative
            if density > 0.0:
                return float(bin_left + remaining / density)
            return float(0.5 * (bin_left + bin_right))
        cumulative += float(bin_mass)
    return float(dist.energy_eV[-1])


def _ionization_event_fraction(metadata: dict[str, Any]) -> float | str:
    ion = _float_metadata(metadata, "mc_collision_count_ionization")
    accepted = sum(
        _float_metadata(metadata, key)
        for key in (
            "mc_collision_count_elastic",
            "mc_collision_count_excitation",
            "mc_collision_count_ionization",
            "mc_collision_count_attachment",
            "mc_collision_count_superelastic",
        )
    )
    if accepted <= 0.0:
        return ""
    return ion / accepted


def _quantile_rows(distributions: list[Distribution]) -> list[dict[str, object]]:
    return [
        {
            "source": dist.source,
            "case_id": dist.case_id,
            "E_over_N_Td": dist.e_over_n_Td,
            "E50_eV": _energy_quantile(dist, 0.50),
            "E90_eV": _energy_quantile(dist, 0.90),
            "E99_eV": _energy_quantile(dist, 0.99),
        }
        for dist in distributions
    ]


def _metrics_rows(distributions: list[Distribution]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for dist in distributions:
        metadata = dist.metadata
        row: dict[str, object] = {
            "source": dist.source,
            "case_id": dist.case_id,
            "E_over_N_Td": dist.e_over_n_Td,
            "mean_energy_eV": dist.mean_energy_eV,
            "E50_eV": _energy_quantile(dist, 0.50),
            "E90_eV": _energy_quantile(dist, 0.90),
            "E99_eV": _energy_quantile(dist, 0.99),
            "transport_definition": metadata_value(metadata, "transport_definition"),
            "mc_population_model": metadata_value(metadata, "mc_population_model"),
            "swarm_population_treatment": metadata_value(
                metadata, "swarm_population_treatment"
            ),
            "eedf_estimator": metadata_value(metadata, "eedf_estimator"),
            "nonconservative_growth_treatment": metadata_value(
                metadata, "nonconservative_growth_treatment"
            ),
            "ionization_event_fraction": _ionization_event_fraction(metadata),
            "mc_physical_branching_gap_eV": metadata_value(
                metadata, "mc_physical_branching_gap_eV"
            ),
            "mc_secondary_electron_count": metadata_value(
                metadata, "mc_secondary_electron_count"
            ),
            "mc_branching_resample_count": metadata_value(
                metadata, "mc_branching_resample_count"
            ),
            "mc_population_total_weight_final": metadata_value(
                metadata, "mc_population_total_weight_final"
            ),
            "mc_population_log_growth_estimate_s_inv": metadata_value(
                metadata, "mc_population_log_growth_estimate_s_inv"
            ),
            "mc_population_weight_cv": metadata_value(
                metadata, "mc_population_weight_cv"
            ),
            "mc_population_resampling_energy_adjustment_eV": metadata_value(
                metadata, "mc_population_resampling_energy_adjustment_eV"
            ),
            "mc_energy_balance_status": metadata_value(
                metadata, "mc_energy_balance_status"
            ),
            "mc_tail_uncertainty_status": metadata_value(
                metadata, "mc_tail_uncertainty_status"
            ),
            "mc_energy_samples_above_xs_max_fraction": metadata_value(
                metadata, "mc_energy_samples_above_xs_max_fraction"
            ),
        }
        for threshold in TAIL_THRESHOLDS:
            row[f"tail_survival_gt_{int(threshold)}_eV"] = _tail_probability(
                dist, threshold
            )
        rows.append(row)
    return rows


def _band_rows(distributions: list[Distribution]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for dist in distributions:
        for left, right in DEFAULT_BANDS:
            mass, mean_contribution = _band_integrals(dist, left, right)
            mask = (dist.energy_eV >= left) & (dist.energy_eV < right)
            sample_count = ""
            effective_count = ""
            relative_error = ""
            if dist.sample_count is not None and len(dist.sample_count) == len(mask):
                sample_count = float(np.nansum(dist.sample_count[mask]))
            if (
                dist.effective_sample_count is not None
                and len(dist.effective_sample_count) == len(mask)
            ):
                effective_value = float(np.nansum(dist.effective_sample_count[mask]))
                effective_count = effective_value
                if effective_value > 0.0:
                    relative_error = 1.0 / effective_value**0.5
            right_label = "inf" if not np.isfinite(right) else f"{right:g}"
            rows.append(
                {
                    "source": dist.source,
                    "case_id": dist.case_id,
                    "E_over_N_Td": dist.e_over_n_Td,
                    "band": f"{left:g}-{right_label}",
                    "energy_left_eV": left,
                    "energy_right_eV": right if np.isfinite(right) else "",
                    "probability_mass": mass,
                    "mean_energy_contribution_eV": mean_contribution,
                    "survival_probability_at_left": _tail_probability(dist, left),
                    "sample_count": sample_count,
                    "effective_sample_count": effective_count,
                    "relative_standard_error": relative_error,
                }
            )
    return rows


def _case_from_distribution(dist: Distribution) -> SwarmCaseResult:
    return SwarmCaseResult(
        solver=dist.source,
        case_id=dist.case_id,
        e_over_n_Td=dist.e_over_n_Td,
        mean_energy_eV=dist.mean_energy_eV,
        drift_velocity_m_s=0.0,
        mobility_m2_V_s=0.0,
        reduced_mobility_m2_V_s_m3=0.0,
        diffusion_L_m2_s=0.0,
        diffusion_T_m2_s=0.0,
        reduced_diffusion_L_m2_s_m3=0.0,
        reduced_diffusion_T_m2_s_m3=0.0,
        net_ionization_frequency_s=0.0,
        effective_townsend_m2=0.0,
        energy_eV=dist.energy_eV,
        eedf=dist.eedf_eV_inv,
        eepf=dist.eedf_eV_inv / np.sqrt(np.maximum(dist.energy_eV, 1.0e-300)),
        energy_widths_eV=dist.widths_eV,
        eedf_counts=None
        if dist.sample_count is None
        else dist.sample_count.astype(int),
        eedf_effective_counts=dist.effective_sample_count,
        metadata=dist.metadata,
    )


def _group_by_eover(distributions: list[Distribution]) -> dict[float, list[Distribution]]:
    grouped: dict[float, list[Distribution]] = {}
    for dist in distributions:
        grouped.setdefault(float(dist.e_over_n_Td), []).append(dist)
    return grouped


def _comparison_rows(distributions: list[Distribution]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for eover, group in _group_by_eover(distributions).items():
        internal = next((dist for dist in group if dist.source == "internal_mc"), None)
        if internal is None:
            continue
        internal_case = _case_from_distribution(internal)
        for reference in group:
            if reference is internal:
                continue
            reference_case = _case_from_distribution(reference)
            metrics = compare_eedf_cases(reference_case, internal_case).metrics
            signed = internal.mean_energy_eV - reference.mean_energy_eV
            rows.append(
                {
                    "E_over_N_Td": eover,
                    "reference_source": reference.source,
                    "candidate_source": internal.source,
                    "eedf_relative_l1": metrics.get("eedf_relative_l1", ""),
                    "log_tail_error": metrics.get("log_tail_error", ""),
                    "tail_probability_difference": metrics.get(
                        "tail_probability_difference", ""
                    ),
                    "mean_energy_reference_eV": reference.mean_energy_eV,
                    "mean_energy_candidate_eV": internal.mean_energy_eV,
                    "mean_energy_signed_difference_eV": signed,
                    "mean_energy_signed_relative_difference": signed
                    / max(abs(reference.mean_energy_eV), 1.0e-300),
                    "E50_difference_eV": metrics.get("E50_difference_eV", ""),
                    "E90_difference_eV": metrics.get("E90_difference_eV", ""),
                    "E99_difference_eV": metrics.get("E99_difference_eV", ""),
                }
            )
    return rows


def _dominant_band(
    bands: list[dict[str, object]],
    eover: float,
    reference_source: str,
    candidate_source: str,
) -> tuple[str, float]:
    ref = {
        str(row["band"]): float(row["mean_energy_contribution_eV"])
        for row in bands
        if row["source"] == reference_source and float(row["E_over_N_Td"]) == eover
    }
    cand = {
        str(row["band"]): float(row["mean_energy_contribution_eV"])
        for row in bands
        if row["source"] == candidate_source and float(row["E_over_N_Td"]) == eover
    }
    best_band = ""
    best_diff = 0.0
    for band, ref_value in ref.items():
        diff = cand.get(band, 0.0) - ref_value
        if abs(diff) > abs(best_diff):
            best_band = band
            best_diff = diff
    return best_band, best_diff


def _assessment_rows(
    distributions: list[Distribution],
    comparisons: list[dict[str, object]],
    bands: list[dict[str, object]],
) -> list[dict[str, object]]:
    by_key = {
        (dist.source, float(dist.e_over_n_Td)): dist for dist in distributions
    }
    rows: list[dict[str, object]] = []
    for row in comparisons:
        eover = float(row["E_over_N_Td"])
        candidate_source = str(row["candidate_source"])
        reference_source = str(row["reference_source"])
        candidate = by_key.get((candidate_source, eover))
        metadata = candidate.metadata if candidate is not None else {}
        xs_fraction = _float_metadata(metadata, "mc_energy_samples_above_xs_max_fraction")
        tail_status = str(metadata_value(metadata, "mc_tail_uncertainty_status"))
        tail_probability = _tail_probability(candidate, 40.0) if candidate is not None else 0.0
        branching_gap = _float_metadata(metadata, "mc_physical_branching_gap_eV")
        ion_fraction = _ionization_event_fraction(metadata)
        band, band_diff = _dominant_band(
            bands, eover, reference_source, candidate_source
        )
        if xs_fraction > 0.01 or (
            tail_status not in {"", "ok"} and tail_probability > 1.0e-8
        ):
            classification = "numerical_tail_or_xs_range_limited"
        elif branching_gap > 0.0 and ion_fraction != "" and float(ion_fraction) > 0.01:
            classification = "nonconservative_population_or_branching_difference"
        else:
            classification = "model_definition_difference_or_reference_definition"
        evidence = (
            f"signed_mean_diff={float(row['mean_energy_signed_difference_eV']):.6g} eV; "
            f"dominant_band={band}; branching_gap={branching_gap:.6g} eV; "
            f"ionization_event_fraction={ion_fraction}; xs_above_fraction={xs_fraction:.3g}"
        )
        rows.append(
            {
                "E_over_N_Td": eover,
                "reference_source": reference_source,
                "candidate_source": candidate_source,
                "classification": classification,
                "dominant_band": band,
                "dominant_band_mean_energy_difference_eV": band_diff,
                "evidence": evidence,
            }
        )
    return rows


def _write_report(
    path: Path,
    *,
    metrics: list[dict[str, object]],
    comparisons: list[dict[str, object]],
    assessments: list[dict[str, object]],
) -> None:
    lines = [
        "# Internal MC Physics Triage",
        "",
        "This development report compares raw EEDF probability mass and model "
        "metadata. It does not declare two-term, MCIG, or another MC result as "
        "the truth, and it does not fit or smooth internal MC output.",
        "",
        "## Distributions",
        "",
    ]
    for row in metrics:
        lines.append(
            "- {source} E/N={E_over_N_Td:g} Td: mean={mean_energy_eV:.6g} eV, "
            "E90={E90_eV:.6g} eV, tail>40={tail_survival_gt_40_eV:.6g}".format(
                **row
            )
        )
    if comparisons:
        lines.extend(["", "## Internal MC Comparisons", ""])
        for row in comparisons:
            lines.append(
                "- {reference_source} -> {candidate_source} E/N={E_over_N_Td:g} Td: "
                "signed mean diff={mean_energy_signed_difference_eV:.6g} eV, "
                "L1={eedf_relative_l1:.6g}".format(**row)
            )
    if assessments:
        lines.extend(["", "## Difference Classification", ""])
        for row in assessments:
            lines.append(
                "- E/N={E_over_N_Td:g} Td vs {reference_source}: "
                "{classification}; {evidence}".format(**row)
            )
    lines.extend(
        [
            "",
            "## Interpretation Guardrails",
            "",
            "- Nonconservative ionization can separate flux and bulk swarm quantities.",
            "- Fixed-particle single-daughter MC is not a full growth-population or "
            "branching MC model.",
            "- Differences should be read through band mass, quantiles, survival "
            "probability, and the EEDF sampling definition before assigning a "
            "code defect.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def _parse_reference_arg(value: str) -> tuple[str, Path]:
    if "=" in value:
        reference_id, path = value.split("=", 1)
        if not reference_id:
            raise ValueError("reference id must not be empty")
        return reference_id, Path(path)
    path = Path(value)
    return path.stem, path


def run_physics_triage(
    *,
    output_dir: Path,
    config_path: Path | None = None,
    internal_tail: Path | None = None,
    internal_summary: Path | None = None,
    references: list[tuple[str, Path]] | None = None,
    include_internal: bool = True,
    include_two_term: bool = True,
) -> tuple[Path, ...]:
    distributions: list[Distribution] = []
    if config_path is not None:
        distributions.extend(
            _run_config_distributions(
                config_path,
                include_internal=include_internal,
                include_two_term=include_two_term,
            )
        )
    if internal_tail is not None:
        distributions.extend(
            _distributions_from_internal_tail(
                internal_tail,
                summary_path=internal_summary,
            )
        )
    for reference_id, path in references or []:
        distributions.extend(_reference_distributions(path, reference_id))
    if not distributions:
        raise ValueError("no distributions to triage")

    output_dir.mkdir(parents=True, exist_ok=True)
    metrics = _metrics_rows(distributions)
    quantiles = _quantile_rows(distributions)
    bands = _band_rows(distributions)
    comparisons = _comparison_rows(distributions)
    assessments = _assessment_rows(distributions, comparisons, bands)

    metric_path = output_dir / "internal_mc_physics_triage_metrics.csv"
    quantile_path = output_dir / "internal_mc_physics_triage_quantiles.csv"
    band_path = output_dir / "internal_mc_physics_triage_band_mass.csv"
    comparison_path = output_dir / "internal_mc_physics_triage_pair_metrics.csv"
    assessment_path = output_dir / "internal_mc_physics_triage_assessment.csv"
    report_path = output_dir / "internal_mc_physics_triage_report.md"
    write_csv(metric_path, metrics, METRIC_FIELDS)
    write_csv(
        quantile_path,
        quantiles,
        ["source", "case_id", "E_over_N_Td", "E50_eV", "E90_eV", "E99_eV"],
    )
    write_csv(band_path, bands, BAND_FIELDS)
    write_csv(comparison_path, comparisons, PAIR_FIELDS)
    write_csv(assessment_path, assessments, ASSESSMENT_FIELDS)
    _write_report(
        report_path,
        metrics=metrics,
        comparisons=comparisons,
        assessments=assessments,
    )
    return (
        metric_path,
        quantile_path,
        band_path,
        comparison_path,
        assessment_path,
        report_path,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--internal-tail", "--internal", type=Path)
    parser.add_argument("--internal-summary", type=Path)
    parser.add_argument(
        "--reference",
        "--references",
        action="append",
        default=[],
        help="External reference as ID=PATH or PATH. CSV uses electron_swarm_reference_csv columns.",
    )
    parser.add_argument("--mcig-reference", type=Path)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/benchmarks/internal_mc_physics_triage"),
    )
    parser.add_argument("--skip-run-internal", action="store_true")
    parser.add_argument("--skip-run-two-term", action="store_true")
    args = parser.parse_args()

    references = [_parse_reference_arg(value) for value in args.reference]
    if args.mcig_reference is not None:
        references.append(("mcig", args.mcig_reference))
    paths = run_physics_triage(
        output_dir=args.output_dir,
        config_path=args.config,
        internal_tail=args.internal_tail,
        internal_summary=args.internal_summary,
        references=references,
        include_internal=not args.skip_run_internal,
        include_two_term=not args.skip_run_two_term,
    )
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
