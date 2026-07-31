"""Reproducible analysis for the external Swarm--COMSOL benchmark.

Raw inputs are read-only.  All derived tables, QA records, source hashes,
figure contracts, and figures are written below the report directory selected
in ``analysis_inputs.json``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sqlite3
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch
import numpy as np
import pandas as pd


REPORT_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CONFIG = Path(__file__).with_name("analysis_inputs.json")

BLUE = "#1769AA"
BLUE_LIGHT = "#8DB8D8"
BLUE_OPEN = "#DDECF6"
CHARCOAL = "#263238"
GREY = "#78838C"
LIGHT_GREY = "#D7DDE1"
GOLD = "#C79318"
GOLD_LIGHT = "#F1D895"
PINK = "#B85C7A"
WHITE = "#FFFFFF"
BOLTZMANN_CONSTANT_J_K = 1.380649e-23
ELEMENTARY_CHARGE_C = 1.602176634e-19
ELECTRON_MASS_KG = 9.1093837139e-31
NET_FLUX_SPEED_RATIO_GUIDE = 0.1
CURRENT_RETAINED_FRACTIONS: tuple[float, ...] = (
    0.50,
    0.60,
    0.70,
    0.80,
    0.85,
    0.90,
    0.95,
    1.00,
)
TYPICAL_DRIFT_DIFFUSION_GUIDANCE_TD = 500.0

QUANTITIES: tuple[tuple[str, str, str], ...] = (
    ("electron_density", "Electron density", r"m$^{-3}$"),
    ("mean_electron_energy", "Mean electron energy", "eV"),
    ("electric_potential", "Electric potential", "V"),
    ("E_over_N", "E/N", "Td"),
    ("electron_current_density", "Electron current density", r"A m$^{-2}$"),
    (
        "total_current_density",
        r"Electron + Ar$^+$ conductive current",
        r"A m$^{-2}$",
    ),
    (
        "excitation_source",
        "Direct excitation source (eir2)",
        r"m$^{-3}$ s$^{-1}$",
    ),
    (
        "ionization_source",
        "Direct ionization source (eir4)",
        r"m$^{-3}$ s$^{-1}$",
    ),
)

REQUIRED_PROFILE_COLUMNS = {
    "x",
    "electron_density",
    "mean_electron_energy",
    "electric_potential",
    "electron_current_density",
    "ion_current_density",
    "total_current_density",
    "excitation_source",
    "ionization_source",
    "E_over_N",
    "applied_voltage",
    "gas_pressure",
}


@dataclass(frozen=True)
class ProfileSource:
    path: Path
    label: str
    mesh_elements: int
    path_id: str


def _json_default(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, Path):
        return value.as_posix()
    raise TypeError(f"Cannot serialize {type(value)!r}")


def write_json(path: Path, payload: Mapping[str, Any] | Sequence[Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, default=_json_default) + "\n",
        encoding="utf-8",
    )


def resolve_path(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else PROJECT_ROOT / path


def relative_path(path: Path) -> str:
    try:
        return path.resolve().relative_to(PROJECT_ROOT.resolve()).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    resolved = path.resolve()
    config = json.loads(resolved.read_text(encoding="utf-8"))
    config["_config_path"] = resolved
    return config


def profile_source(
    item: Mapping[str, Any], *, path_id: str, default_label: str
) -> ProfileSource:
    return ProfileSource(
        path=resolve_path(item["path"]),
        label=str(item.get("label", default_label)),
        mesh_elements=int(item.get("mesh_elements", 0)),
        path_id=str(item.get("path_id", path_id)),
    )


def read_profile(source: ProfileSource) -> pd.DataFrame:
    frame = pd.read_csv(source.path)
    missing = sorted(REQUIRED_PROFILE_COLUMNS - set(frame.columns))
    if missing:
        raise ValueError(f"{source.path} is missing profile columns: {missing}")
    numeric = frame.loc[:, sorted(REQUIRED_PROFILE_COLUMNS)].apply(
        pd.to_numeric, errors="coerce"
    )
    frame.loc[:, numeric.columns] = numeric
    frame = frame.sort_values("x", kind="stable").reset_index(drop=True)
    return frame


def convert_e_over_n_to_td(values: np.ndarray, unit: str) -> np.ndarray:
    normalized = unit.lower().replace(" ", "").replace("^", "")
    if normalized in {"td", "townsend"}:
        return values.astype(float)
    if normalized in {"vm2", "v*m2", "v·m2"}:
        return values.astype(float) / 1.0e-21
    raise ValueError(f"Unsupported E/N profile unit: {unit!r}")


def _deduplicated_xy(frame: pd.DataFrame, field: str) -> tuple[np.ndarray, np.ndarray]:
    selected = frame.loc[:, ["x", field]].dropna().sort_values("x")
    if selected["x"].duplicated().any():
        selected = selected.groupby("x", as_index=False, sort=True)[field].mean()
    return selected["x"].to_numpy(float), selected[field].to_numpy(float)


def common_grid(
    left: pd.DataFrame, right: pd.DataFrame, field: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lx, ly = _deduplicated_xy(left, field)
    rx, ry = _deduplicated_xy(right, field)
    lower = max(lx.min(), rx.min())
    upper = min(lx.max(), rx.max())
    if not upper > lower:
        raise ValueError(f"No spatial overlap for {field}")
    grid = np.unique(
        np.concatenate(
            (
                lx[(lx >= lower) & (lx <= upper)],
                rx[(rx >= lower) & (rx <= upper)],
                np.array([lower, upper]),
            )
        )
    )
    return grid, np.interp(grid, lx, ly), np.interp(grid, rx, ry)


def trapezoid_weights(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or len(x) < 2 or np.any(np.diff(x) <= 0):
        raise ValueError("x must be a strictly increasing one-dimensional grid")
    weights = np.empty_like(x)
    weights[0] = 0.5 * (x[1] - x[0])
    weights[-1] = 0.5 * (x[-1] - x[-2])
    weights[1:-1] = 0.5 * (x[2:] - x[:-2])
    return weights


def _validated_piecewise_linear_arrays(
    x: np.ndarray, *values: np.ndarray
) -> tuple[np.ndarray, ...]:
    grid = np.asarray(x, dtype=float)
    arrays = tuple(np.asarray(value, dtype=float) for value in values)
    if (
        grid.ndim != 1
        or len(grid) < 2
        or np.any(np.diff(grid) <= 0)
        or any(value.shape != grid.shape for value in arrays)
    ):
        raise ValueError(
            "piecewise-linear arrays require one-dimensional values on a "
            "strictly increasing grid"
        )
    if not np.all(np.isfinite(grid)) or any(
        not np.all(np.isfinite(value)) for value in arrays
    ):
        raise ValueError("piecewise-linear integration requires finite values")
    return (grid, *arrays)


def piecewise_linear_integral(x: np.ndarray, values: np.ndarray) -> float:
    """Integrate the continuous piecewise-linear reconstruction exactly."""

    grid, y = _validated_piecewise_linear_arrays(x, values)
    return float(np.sum(0.5 * np.diff(grid) * (y[:-1] + y[1:])))


def piecewise_linear_product_integral(
    x: np.ndarray, left: np.ndarray, right: np.ndarray
) -> float:
    """Integrate the product of two piecewise-linear reconstructions exactly."""

    grid, a, b = _validated_piecewise_linear_arrays(x, left, right)
    cell_integrals = np.diff(grid) * (
        2.0 * a[:-1] * b[:-1]
        + a[:-1] * b[1:]
        + a[1:] * b[:-1]
        + 2.0 * a[1:] * b[1:]
    ) / 6.0
    return float(np.sum(cell_integrals))


def piecewise_linear_squared_integral(
    x: np.ndarray, values: np.ndarray
) -> float:
    return piecewise_linear_product_integral(x, values, values)


def trapezoidal_nodal_correlation(
    x: np.ndarray, left: np.ndarray, right: np.ndarray
) -> float:
    """Legacy sensitivity calculation using trapezoidal nodal weights."""

    weights = trapezoid_weights(x)
    total_weight = weights.sum()
    left_mean = float(np.dot(weights, left) / total_weight)
    right_mean = float(np.dot(weights, right) / total_weight)
    left_centered = left - left_mean
    right_centered = right - right_mean
    denominator = math.sqrt(
        float(np.dot(weights, left_centered**2))
        * float(np.dot(weights, right_centered**2))
    )
    if denominator == 0.0:
        return float("nan")
    return float(np.dot(weights, left_centered * right_centered) / denominator)


def weighted_correlation(x: np.ndarray, left: np.ndarray, right: np.ndarray) -> float:
    """Pearson correlation of continuous piecewise-linear reconstructions."""

    grid, left_values, right_values = _validated_piecewise_linear_arrays(
        x, left, right
    )
    length = float(grid[-1] - grid[0])
    left_mean = piecewise_linear_integral(grid, left_values) / length
    right_mean = piecewise_linear_integral(grid, right_values) / length
    left_centered = left_values - left_mean
    right_centered = right_values - right_mean
    denominator = math.sqrt(
        piecewise_linear_squared_integral(grid, left_centered)
        * piecewise_linear_squared_integral(grid, right_centered)
    )
    if denominator == 0.0:
        return float("nan")
    return (
        piecewise_linear_product_integral(
            grid, left_centered, right_centered
        )
        / denominator
    )


def comparison_metrics(
    external: pd.DataFrame, reference: pd.DataFrame, e_over_n_unit: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    metric_rows: list[dict[str, Any]] = []
    aligned = pd.DataFrame()
    for field, label, unit in QUANTITIES:
        x, ext, ref = common_grid(external, reference, field)
        if field == "E_over_N":
            ext = convert_e_over_n_to_td(ext, e_over_n_unit)
            ref = convert_e_over_n_to_td(ref, e_over_n_unit)
        difference = ext - ref
        numerator = piecewise_linear_squared_integral(x, difference)
        denominator = piecewise_linear_squared_integral(x, ref)
        relative_l2 = math.sqrt(numerator / denominator) if denominator > 0 else np.nan
        ext_integral = piecewise_linear_integral(x, ext)
        ref_integral = piecewise_linear_integral(x, ref)
        integral_ratio = ext_integral / ref_integral if ref_integral != 0 else np.nan
        trapezoidal_numerator = float(np.trapezoid(difference**2, x))
        trapezoidal_denominator = float(np.trapezoid(ref**2, x))
        trapezoidal_relative_l2 = (
            math.sqrt(trapezoidal_numerator / trapezoidal_denominator)
            if trapezoidal_denominator > 0
            else np.nan
        )
        legacy_denominator = float(np.linalg.norm(ref))
        legacy_l2 = (
            float(np.linalg.norm(ext - ref) / legacy_denominator)
            if legacy_denominator > 0
            else np.nan
        )
        legacy_corr = (
            float(np.corrcoef(ext, ref)[0, 1])
            if np.std(ext) > 0 and np.std(ref) > 0
            else np.nan
        )
        metric_rows.append(
            {
                "quantity": field,
                "quantity_label": label,
                "unit": unit,
                "relative_L2_spatial_weighted": relative_l2,
                "relative_L2_trapezoidal_of_squared_samples_sensitivity": (
                    trapezoidal_relative_l2
                ),
                "integral_external": ext_integral,
                "integral_reference": ref_integral,
                "integral_ratio_external_over_reference": integral_ratio,
                "signed_integral_external": ext_integral,
                "signed_integral_reference": ref_integral,
                "signed_integral_ratio_external_over_reference": integral_ratio,
                "shape_correlation_spatial_weighted": weighted_correlation(x, ext, ref),
                "shape_correlation_trapezoidal_nodal_sensitivity": (
                    trapezoidal_nodal_correlation(x, ext, ref)
                ),
                "legacy_L2_unweighted_discrete": legacy_l2,
                "legacy_correlation_unweighted_discrete": legacy_corr,
                "common_grid_rows": len(x),
                "x_min_m": float(x.min()),
                "x_max_m": float(x.max()),
            }
        )
        quantity_aligned = pd.DataFrame(
            {
                "quantity": field,
                "x_m": x,
                "external_value": ext,
                "reference_value": ref,
                "unit": unit,
            }
        )
        aligned = pd.concat((aligned, quantity_aligned), ignore_index=True)
    return pd.DataFrame(metric_rows), aligned


def clip_piecewise_linear_interval(
    x: np.ndarray,
    values: Sequence[np.ndarray],
    lower: float,
    upper: float,
) -> tuple[np.ndarray, list[np.ndarray]]:
    grid, *series = _validated_piecewise_linear_arrays(x, *values)
    clipped_lower = max(float(lower), float(grid.min()))
    clipped_upper = min(float(upper), float(grid.max()))
    if not clipped_upper > clipped_lower:
        raise ValueError("requested interval has no positive-length overlap")
    inner = grid[(grid > clipped_lower) & (grid < clipped_upper)]
    clipped_x = np.unique(
        np.concatenate(([clipped_lower], inner, [clipped_upper]))
    )
    return clipped_x, [
        np.interp(clipped_x, grid, value) for value in series
    ]


def regional_error_attribution(
    external: pd.DataFrame,
    reference: pd.DataFrame,
    e_over_n_unit: str,
) -> pd.DataFrame:
    """Partition comparison error into three geometric axial windows."""

    rows: list[dict[str, Any]] = []
    for field, label, unit in QUANTITIES:
        x, ext, ref = common_grid(external, reference, field)
        if field == "E_over_N":
            ext = convert_e_over_n_to_td(ext, e_over_n_unit)
            ref = convert_e_over_n_to_td(ref, e_over_n_unit)
        lower = float(x.min())
        upper = float(x.max())
        length = upper - lower
        regions = (
            ("cathode_side_10_percent", lower, lower + 0.1 * length),
            (
                "central_80_percent",
                lower + 0.1 * length,
                upper - 0.1 * length,
            ),
            ("anode_side_10_percent", upper - 0.1 * length, upper),
        )
        total_error_norm = piecewise_linear_squared_integral(x, ext - ref)
        total_reference_norm = piecewise_linear_squared_integral(x, ref)
        for region, region_lower, region_upper in regions:
            region_x, (region_ext, region_ref) = clip_piecewise_linear_interval(
                x, (ext, ref), region_lower, region_upper
            )
            region_error_norm = piecewise_linear_squared_integral(
                region_x, region_ext - region_ref
            )
            region_reference_norm = piecewise_linear_squared_integral(
                region_x, region_ref
            )
            external_integral = piecewise_linear_integral(region_x, region_ext)
            reference_integral = piecewise_linear_integral(region_x, region_ref)
            rows.append(
                {
                    "quantity": field,
                    "quantity_label": label,
                    "unit": unit,
                    "region": region,
                    "x_min_m": region_lower,
                    "x_max_m": region_upper,
                    "region_length_m": region_upper - region_lower,
                    "squared_error_integral": region_error_norm,
                    "squared_error_share": (
                        region_error_norm / total_error_norm
                        if total_error_norm > 0
                        else np.nan
                    ),
                    "region_normalized_relative_L2": (
                        math.sqrt(region_error_norm / region_reference_norm)
                        if region_reference_norm > 0
                        else np.nan
                    ),
                    "reference_squared_norm_integral": region_reference_norm,
                    "reference_norm_share": (
                        region_reference_norm / total_reference_norm
                        if total_reference_norm > 0
                        else np.nan
                    ),
                    "signed_integral_external": external_integral,
                    "signed_integral_reference": reference_integral,
                    "signed_integral_ratio_external_over_reference": (
                        external_integral / reference_integral
                        if reference_integral != 0
                        else np.nan
                    ),
                    "integration_definition": (
                        "exact_product_integral_of_union_grid_piecewise_linear_"
                        "reconstruction"
                    ),
                }
            )
    return pd.DataFrame(rows)


def clip_profile_interval(
    frame: pd.DataFrame, field: str, lower: float, upper: float
) -> tuple[np.ndarray, np.ndarray]:
    x, y = _deduplicated_xy(frame, field)
    lower = max(lower, float(x.min()))
    upper = min(upper, float(x.max()))
    inner = x[(x > lower) & (x < upper)]
    clipped_x = np.unique(np.concatenate(([lower], inner, [upper])))
    return clipped_x, np.interp(clipped_x, x, y)


def spatial_rsd(x: np.ndarray, values: np.ndarray) -> tuple[float, float, float]:
    """Population RSD of a continuous piecewise-linear reconstruction."""

    grid, y = _validated_piecewise_linear_arrays(x, values)
    length = float(grid[-1] - grid[0])
    mean = piecewise_linear_integral(grid, y) / length
    variance = piecewise_linear_squared_integral(grid, y - mean) / length
    standard_deviation = math.sqrt(max(variance, 0.0))
    rsd = standard_deviation / abs(mean) if mean != 0 else np.nan
    return mean, standard_deviation, rsd


def spatial_rsd_trapezoidal_sensitivity(
    x: np.ndarray, values: np.ndarray
) -> tuple[float, float, float]:
    """Legacy trapezoidal-of-squared-samples RSD sensitivity calculation."""

    weights = trapezoid_weights(x)
    mean = float(np.dot(weights, values) / weights.sum())
    variance = float(np.dot(weights, (values - mean) ** 2) / weights.sum())
    standard_deviation = math.sqrt(max(variance, 0.0))
    rsd = standard_deviation / abs(mean) if mean != 0 else np.nan
    return mean, standard_deviation, rsd


def current_rsd_sensitivity(
    profiles: Sequence[tuple[ProfileSource, pd.DataFrame]],
    retained_fractions: Sequence[float] = CURRENT_RETAINED_FRACTIONS,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for source, frame in profiles:
        domain_lower = float(frame["x"].min())
        domain_upper = float(frame["x"].max())
        length = domain_upper - domain_lower
        for retained_fraction in retained_fractions:
            retained = float(retained_fraction)
            if not 0.0 < retained <= 1.0:
                raise ValueError("retained fractions must lie in (0, 1]")
            excluded_each_side = 0.5 * (1.0 - retained)
            lower = domain_lower + excluded_each_side * length
            upper = domain_upper - excluded_each_side * length
            x, current = clip_profile_interval(
                frame, "total_current_density", lower, upper
            )
            mean, standard_deviation, rsd = spatial_rsd(x, current)
            rows.append(
                {
                    "path": source.path_id,
                    "path_label": source.label,
                    "mesh_elements": source.mesh_elements,
                    "central_retained_fraction": retained,
                    "central_retained_percent": 100.0 * retained,
                    "excluded_fraction_each_boundary": excluded_each_side,
                    "x_min_m": lower,
                    "x_max_m": upper,
                    "piecewise_linear_exact_mean_A_m2": mean,
                    "piecewise_linear_exact_standard_deviation_A_m2": (
                        standard_deviation
                    ),
                    "piecewise_linear_exact_relative_standard_deviation": rsd,
                    "metric_definition": (
                        "piecewise_linear_exact_spatial_population_RSD"
                    ),
                }
            )
    return pd.DataFrame(rows)


def current_uniformity(
    profiles: Sequence[tuple[ProfileSource, pd.DataFrame]],
    fallback_path: Path | None,
) -> tuple[pd.DataFrame, bool]:
    rows: list[dict[str, Any]] = []
    covered: set[tuple[str, int, str]] = set()
    for source, frame in profiles:
        x = frame["x"].to_numpy(float)
        domain_lower = float(np.nanmin(x))
        domain_upper = float(np.nanmax(x))
        length = domain_upper - domain_lower
        intervals = {
            "full_domain": (domain_lower, domain_upper),
            "central_80_percent": (
                domain_lower + 0.1 * length,
                domain_upper - 0.1 * length,
            ),
        }
        for region, (lower, upper) in intervals.items():
            clipped_x, current = clip_profile_interval(
                frame, "total_current_density", lower, upper
            )
            mean, standard_deviation, rsd = spatial_rsd(clipped_x, current)
            (
                trapezoidal_mean,
                trapezoidal_standard_deviation,
                trapezoidal_rsd,
            ) = spatial_rsd_trapezoidal_sensitivity(clipped_x, current)
            rows.append(
                {
                    "path": source.path_id,
                    "path_label": source.label,
                    "mesh_elements": source.mesh_elements,
                    "domain": region,
                    "x_min_m": lower,
                    "x_max_m": upper,
                    "weighted_mean_A_m2": mean,
                    "weighted_standard_deviation_A_m2": standard_deviation,
                    "relative_standard_deviation": rsd,
                    "trapezoidal_nodal_mean_A_m2_sensitivity": trapezoidal_mean,
                    "trapezoidal_of_squared_samples_standard_deviation_A_m2_sensitivity": (
                        trapezoidal_standard_deviation
                    ),
                    "trapezoidal_of_squared_samples_RSD_sensitivity": (
                        trapezoidal_rsd
                    ),
                    "metric_source": "raw_profile",
                    "metric_definition": (
                        "piecewise_linear_exact_spatial_population_RSD"
                    ),
                    "comparable_to_spatial_weighted": True,
                }
            )
            covered.add((source.path_id, source.mesh_elements, region))
    fallback_used = False
    if fallback_path is not None and fallback_path.exists():
        fallback = pd.read_csv(fallback_path)
        for record in fallback.to_dict("records"):
            key = (
                str(record.get("path", "external_swarm")),
                int(record["mesh_elements"]),
                str(record["domain"]),
            )
            fallback_comparable = str(
                record.get("comparable_to_spatial_weighted", "false")
            ).strip().lower() in {"true", "1", "yes"}
            if key in covered and fallback_comparable:
                continue
            fallback_used = True
            rows.append(
                {
                    "path": key[0],
                    "path_label": "External Swarm tables (archived LEA)",
                    "mesh_elements": key[1],
                    "domain": key[2],
                    "x_min_m": np.nan,
                    "x_max_m": np.nan,
                    "weighted_mean_A_m2": np.nan,
                    "weighted_standard_deviation_A_m2": np.nan,
                    "relative_standard_deviation": float(
                        record["relative_standard_deviation"]
                    ),
                    "trapezoidal_nodal_mean_A_m2_sensitivity": np.nan,
                    "trapezoidal_of_squared_samples_standard_deviation_A_m2_sensitivity": (
                        np.nan
                    ),
                    "trapezoidal_of_squared_samples_RSD_sensitivity": np.nan,
                    "metric_source": str(
                        record.get("evidence_status", "summary_fallback")
                    ),
                    "metric_definition": str(
                        record.get(
                            "metric_definition",
                            "historical_summary_definition_unspecified",
                        )
                    ),
                    "comparable_to_spatial_weighted": fallback_comparable,
                }
            )
    result = pd.DataFrame(rows).sort_values(
        ["path", "mesh_elements", "domain"], kind="stable"
    )
    return result.reset_index(drop=True), fallback_used


def normalize_runtime(runtime: pd.DataFrame) -> pd.DataFrame:
    required = {"scope", "path", "seconds", "repetition"}
    missing = required - set(runtime.columns)
    if missing:
        raise ValueError(f"Runtime CSV is missing columns: {sorted(missing)}")
    normalized = runtime.copy()
    normalized["seconds"] = pd.to_numeric(normalized["seconds"], errors="coerce")
    normalized["repetition"] = pd.to_numeric(
        normalized["repetition"], errors="coerce"
    ).astype("Int64")
    scopes = set(normalized["scope"].astype(str))
    external = normalized[normalized["path"] == "external_swarm"]
    additions: list[pd.DataFrame] = []
    if "solve_only" not in scopes and not external[external["scope"] == "solve"].empty:
        solve = external[external["scope"] == "solve"].copy()
        solve["scope"] = "solve_only"
        additions.append(solve)
    if "end_to_end" not in scopes:
        stage_names = {"apply", "verify", "solve", "export"}
        stages = external[external["scope"].isin(stage_names)]
        if not stages.empty:
            pivot = stages.pivot_table(
                index=["path", "repetition"],
                columns="scope",
                values="seconds",
                aggfunc="first",
            )
            complete = pivot.dropna(subset=sorted(stage_names), how="any")
            if not complete.empty:
                derived = complete.loc[:, sorted(stage_names)].sum(axis=1).reset_index()
                derived = derived.rename(columns={0: "seconds"})
                derived["scope"] = "end_to_end"
                derived["evidence_status"] = "derived_from_stage_sum"
                additions.append(derived)
    if additions:
        normalized = pd.concat((normalized, *additions), ignore_index=True)
    return normalized


def runtime_metrics(
    runtime: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    normalized = normalize_runtime(runtime)
    summary = (
        normalized.groupby(["scope", "path"], sort=True)["seconds"]
        .agg(["count", "median", "min", "max"])
        .reset_index()
        .rename(
            columns={
                "count": "repetitions",
                "median": "median_seconds",
                "min": "min_seconds",
                "max": "max_seconds",
            }
        )
    )
    scope_interpretations = {
        "solve_only": (
            "historical COMSOL run-stage timing; not isolated nonlinear-solver time"
        ),
        "end_to_end": "historical complete-path observation",
        "apply": "external table apply stage",
        "verify": "external table verification stage",
        "solve": "external COMSOL run stage",
        "export": "external result export stage",
    }
    summary["scope_interpretation"] = summary["scope"].map(
        scope_interpretations
    )
    speed_rows: list[dict[str, Any]] = []
    for scope in ("solve_only", "end_to_end"):
        subset = summary[summary["scope"] == scope].set_index("path")
        if {"external_swarm", "builtin_reference"} <= set(subset.index):
            external = float(subset.loc["external_swarm", "median_seconds"])
            builtin = float(subset.loc["builtin_reference", "median_seconds"])
            speed_rows.append(
                {
                    "scope": scope,
                    "external_median_seconds": external,
                    "builtin_median_seconds": builtin,
                    "speedup_builtin_over_external": builtin / external,
                    "external_repetitions": int(
                        subset.loc["external_swarm", "repetitions"]
                    ),
                    "builtin_repetitions": int(
                        subset.loc["builtin_reference", "repetitions"]
                    ),
                }
            )
    return summary, pd.DataFrame(speed_rows)


def _piecewise_condition_measure(
    x: np.ndarray, control: np.ndarray, predicate: str, threshold: float
) -> float:
    if predicate not in {"below", "above"}:
        raise ValueError(predicate)
    measure = 0.0
    for x0, x1, y0, y1 in zip(x[:-1], x[1:], control[:-1], control[1:]):
        interval = x1 - x0
        left = y0 < threshold if predicate == "below" else y0 > threshold
        right = y1 < threshold if predicate == "below" else y1 > threshold
        if left and right:
            measure += interval
        elif left != right and y1 != y0:
            crossing_fraction = (threshold - y0) / (y1 - y0)
            crossing_fraction = float(np.clip(crossing_fraction, 0.0, 1.0))
            measure += (
                interval * crossing_fraction
                if left
                else interval * (1.0 - crossing_fraction)
            )
    return measure


def _condition_integral_fraction(
    x: np.ndarray,
    control: np.ndarray,
    values: np.ndarray,
    predicate: str,
    threshold: float,
) -> float:
    crossings: list[float] = []
    for crossing_values, crossing_threshold in (
        (control, threshold),
        (values, 0.0),
    ):
        for x0, x1, y0, y1 in zip(
            x[:-1], x[1:], crossing_values[:-1], crossing_values[1:]
        ):
            if y1 == y0:
                continue
            fraction = (crossing_threshold - y0) / (y1 - y0)
            if 0.0 < fraction < 1.0:
                crossings.append(float(x0 + fraction * (x1 - x0)))
    augmented_x = np.unique(np.concatenate((x, np.asarray(crossings))))
    augmented_control = np.interp(augmented_x, x, control)
    augmented_values = np.abs(np.interp(augmented_x, x, values))
    selected_integral = 0.0
    total_integral = float(np.trapezoid(augmented_values, augmented_x))
    for index in range(len(augmented_x) - 1):
        midpoint_control = 0.5 * (
            augmented_control[index] + augmented_control[index + 1]
        )
        selected = (
            midpoint_control < threshold
            if predicate == "below"
            else midpoint_control > threshold
        )
        if selected:
            selected_integral += 0.5 * (
                augmented_values[index] + augmented_values[index + 1]
            ) * (augmented_x[index + 1] - augmented_x[index])
    return selected_integral / total_integral if total_integral > 0 else np.nan


def clamp_metrics(
    external: pd.DataFrame, bundle: Path, e_over_n_unit: str
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, float]]:
    mean_table = pd.read_csv(bundle / "mean_energy_vs_en.csv")
    computed_transport_path = bundle / "transport_vs_en.csv"
    computed_transport = (
        pd.read_csv(computed_transport_path)
        if computed_transport_path.exists()
        else mean_table
    )
    lookup_transport_path = bundle / "transport_vs_mean_energy.csv"
    lookup_transport = (
        pd.read_csv(lookup_transport_path)
        if lookup_transport_path.exists()
        else computed_transport
    )
    en_min = float(mean_table["E_over_N_Td"].min())
    en_max = float(mean_table["E_over_N_Td"].max())
    computed_energy_min = float(computed_transport["mean_energy_eV"].min())
    computed_energy_max = float(computed_transport["mean_energy_eV"].max())
    lookup_energy_min = float(lookup_transport["mean_energy_eV"].min())
    lookup_energy_max = float(lookup_transport["mean_energy_eV"].max())
    boundary_policy_rows = lookup_transport.loc[
        np.isclose(
            pd.to_numeric(lookup_transport["E_over_N_Td"], errors="coerce"),
            0.0,
        )
    ]
    boundary_policy_constant_max = (
        float(boundary_policy_rows["mean_energy_eV"].max())
        if not boundary_policy_rows.empty
        else lookup_energy_min
    )
    x, en_values_raw = _deduplicated_xy(external, "E_over_N")
    en_values = convert_e_over_n_to_td(en_values_raw, e_over_n_unit)
    length = float(x.max() - x.min())
    energy_x, energy = _deduplicated_xy(external, "mean_electron_energy")
    metric_specs = (
        ("E_over_N_below_table", x, en_values, "below", en_min, "Td"),
        ("E_over_N_above_table", x, en_values, "above", en_max, "Td"),
        (
            "mean_energy_below_computed_domain",
            energy_x,
            energy,
            "below",
            computed_energy_min,
            "eV",
        ),
        (
            "mean_energy_above_computed_domain",
            energy_x,
            energy,
            "above",
            computed_energy_max,
            "eV",
        ),
        (
            "mean_energy_in_constant_boundary_policy_region",
            energy_x,
            energy,
            "below",
            boundary_policy_constant_max,
            "eV",
        ),
        (
            "mean_energy_below_lookup_argument_domain",
            energy_x,
            energy,
            "below",
            lookup_energy_min,
            "eV",
        ),
        (
            "mean_energy_above_lookup_argument_domain",
            energy_x,
            energy,
            "above",
            lookup_energy_max,
            "eV",
        ),
    )
    rows: list[dict[str, Any]] = []
    for metric, axis, control, predicate, threshold, unit in metric_specs:
        rows.append(
            {
                "metric": metric,
                "predicate": predicate,
                "threshold": threshold,
                "threshold_unit": unit,
                "spatial_length_m": _piecewise_condition_measure(
                    axis, control, predicate, threshold
                ),
                "spatial_fraction": _piecewise_condition_measure(
                    axis, control, predicate, threshold
                )
                / float(axis.max() - axis.min()),
                "node_fraction": float(
                    np.mean(control < threshold)
                    if predicate == "below"
                    else np.mean(control > threshold)
                ),
            }
        )
    impact_rows: list[dict[str, Any]] = []
    impact_specs = (
        ("E_over_N_below_table", x, en_values, "below", en_min),
        (
            "mean_energy_below_computed_domain",
            energy_x,
            energy,
            "below",
            computed_energy_min,
        ),
        (
            "mean_energy_in_constant_boundary_policy_region",
            energy_x,
            energy,
            "below",
            boundary_policy_constant_max,
        ),
    )
    for region, control_x, control, predicate, threshold in impact_specs:
        for field, label, unit in QUANTITIES:
            value_x, values = _deduplicated_xy(external, field)
            if field == "E_over_N":
                values = convert_e_over_n_to_td(values, e_over_n_unit)
            if not np.array_equal(value_x, control_x):
                values = np.interp(control_x, value_x, values)
            impact_rows.append(
                {
                    "quantity": field,
                    "quantity_label": label,
                    "quantity_unit": unit,
                    "region": region,
                    "absolute_integral_share": _condition_integral_fraction(
                        control_x,
                        control,
                        values,
                        predicate,
                        threshold,
                    ),
                }
            )
    context = {
        "table_E_over_N_min_Td": en_min,
        "table_E_over_N_max_Td": en_max,
        "profile_E_over_N_min_Td": float(en_values.min()),
        "profile_E_over_N_max_Td": float(en_values.max()),
        "computed_mean_energy_min_eV": computed_energy_min,
        "computed_mean_energy_max_eV": computed_energy_max,
        "lookup_mean_energy_min_eV": lookup_energy_min,
        "lookup_mean_energy_max_eV": lookup_energy_max,
        "boundary_policy_constant_max_mean_energy_eV": (
            boundary_policy_constant_max
        ),
        "profile_mean_energy_min_eV": float(energy.min()),
        "profile_mean_energy_max_eV": float(energy.max()),
        "domain_length_m": length,
    }
    return pd.DataFrame(rows), pd.DataFrame(impact_rows), context


def _unit_interval_real_roots(coefficients: np.ndarray) -> np.ndarray:
    """Return numerically real roots strictly inside the unit interval."""

    coefficients = np.asarray(coefficients, dtype=float)
    if not np.all(np.isfinite(coefficients)):
        raise ValueError("polynomial coefficients must be finite")
    scale = float(np.max(np.abs(coefficients), initial=0.0))
    if scale == 0.0:
        return np.array([], dtype=float)
    normalized = coefficients / scale
    while len(normalized) > 1 and abs(normalized[-1]) < 1.0e-13:
        normalized = normalized[:-1]
    if len(normalized) <= 1:
        return np.array([], dtype=float)
    roots = np.polynomial.polynomial.polyroots(normalized)
    real = [
        float(root.real)
        for root in roots
        if abs(root.imag) <= 1.0e-9 * max(1.0, abs(root.real))
        and 1.0e-12 < root.real < 1.0 - 1.0e-12
    ]
    return np.asarray(sorted(set(round(root, 14) for root in real)), dtype=float)


def _net_flux_speed_ratio_values(
    electron_density_m3: np.ndarray,
    electron_current_density_A_m2: np.ndarray,
    mean_energy_eV: np.ndarray,
) -> np.ndarray:
    density = np.asarray(electron_density_m3, dtype=float)
    current = np.asarray(electron_current_density_A_m2, dtype=float)
    energy = np.asarray(mean_energy_eV, dtype=float)
    if np.any(density <= 0.0) or np.any(energy <= 0.0):
        raise ValueError(
            "net-flux-speed diagnostic requires positive density and mean energy"
        )
    net_flux_speed = np.abs(current) / (ELEMENTARY_CHARGE_C * density)
    energy_equivalent_rms_speed = np.sqrt(
        2.0 * ELEMENTARY_CHARGE_C * energy / ELECTRON_MASS_KG
    )
    return net_flux_speed / energy_equivalent_rms_speed


def _net_flux_ratio_measure_above(
    x: np.ndarray,
    electron_density_m3: np.ndarray,
    electron_current_density_A_m2: np.ndarray,
    mean_energy_eV: np.ndarray,
    threshold: float,
) -> float:
    """Exact length measure for a threshold on piecewise-linear inputs.

    Squaring the positive diagnostic ratio produces a cubic polynomial on
    each interval. Its real roots partition the interval exactly.
    """

    grid, density, current, energy = _validated_piecewise_linear_arrays(
        x,
        electron_density_m3,
        electron_current_density_A_m2,
        mean_energy_eV,
    )
    if threshold < 0.0 or np.any(density <= 0.0) or np.any(energy <= 0.0):
        raise ValueError("threshold must be nonnegative and state values positive")
    selected_length = 0.0
    physical_constant = (
        threshold**2
        * ELEMENTARY_CHARGE_C**2
        * (2.0 * ELEMENTARY_CHARGE_C / ELECTRON_MASS_KG)
    )
    for index, interval_length in enumerate(np.diff(grid)):
        j_end = current[index : index + 2]
        n_end = density[index : index + 2]
        e_end = energy[index : index + 2]
        j_scale = max(float(np.max(np.abs(j_end))), np.finfo(float).tiny)
        n_scale = max(float(np.max(np.abs(n_end))), np.finfo(float).tiny)
        e_scale = max(float(np.max(np.abs(e_end))), np.finfo(float).tiny)
        j_poly = np.array([j_end[0], j_end[1] - j_end[0]]) / j_scale
        n_poly = np.array([n_end[0], n_end[1] - n_end[0]]) / n_scale
        e_poly = np.array([e_end[0], e_end[1] - e_end[0]]) / e_scale
        scaled_constant = (
            physical_constant * n_scale**2 * e_scale / j_scale**2
        )
        polynomial = np.polynomial.polynomial.polyadd(
            np.polynomial.polynomial.polymul(j_poly, j_poly),
            -scaled_constant
            * np.polynomial.polynomial.polymul(
                np.polynomial.polynomial.polymul(n_poly, n_poly),
                e_poly,
            ),
        )
        breakpoints = np.unique(
            np.concatenate(([0.0], _unit_interval_real_roots(polynomial), [1.0]))
        )
        for left, right in zip(breakpoints[:-1], breakpoints[1:]):
            midpoint = 0.5 * (left + right)
            if np.polynomial.polynomial.polyval(midpoint, polynomial) > 0.0:
                selected_length += interval_length * (right - left)
    return float(selected_length)


def _net_flux_ratio_piecewise_max(
    x: np.ndarray,
    electron_density_m3: np.ndarray,
    electron_current_density_A_m2: np.ndarray,
    mean_energy_eV: np.ndarray,
) -> float:
    """Maximum diagnostic ratio over piecewise-linear input profiles."""

    grid, density, current, energy = _validated_piecewise_linear_arrays(
        x,
        electron_density_m3,
        electron_current_density_A_m2,
        mean_energy_eV,
    )
    if np.any(density <= 0.0) or np.any(energy <= 0.0):
        raise ValueError("state values must be positive")
    maximum = 0.0
    for index in range(len(grid) - 1):
        j_end = current[index : index + 2]
        n_end = density[index : index + 2]
        e_end = energy[index : index + 2]
        j_scale = max(float(np.max(np.abs(j_end))), np.finfo(float).tiny)
        n_scale = max(float(np.max(np.abs(n_end))), np.finfo(float).tiny)
        e_scale = max(float(np.max(np.abs(e_end))), np.finfo(float).tiny)
        j_poly = np.array([j_end[0], j_end[1] - j_end[0]]) / j_scale
        n_poly = np.array([n_end[0], n_end[1] - n_end[0]]) / n_scale
        e_poly = np.array([e_end[0], e_end[1] - e_end[0]]) / e_scale
        j_derivative = np.polynomial.polynomial.polyder(j_poly)
        n_derivative = np.polynomial.polynomial.polyder(n_poly)
        e_derivative = np.polynomial.polynomial.polyder(e_poly)
        critical_polynomial = np.polynomial.polynomial.polyadd(
            2.0
            * np.polynomial.polynomial.polymul(
                np.polynomial.polynomial.polymul(j_derivative, n_poly),
                e_poly,
            ),
            -2.0
            * np.polynomial.polynomial.polymul(
                np.polynomial.polynomial.polymul(n_derivative, j_poly),
                e_poly,
            ),
        )
        critical_polynomial = np.polynomial.polynomial.polyadd(
            critical_polynomial,
            -np.polynomial.polynomial.polymul(
                np.polynomial.polynomial.polymul(e_derivative, j_poly),
                n_poly,
            ),
        )
        candidates = np.unique(
            np.concatenate(
                ([0.0, 1.0], _unit_interval_real_roots(critical_polynomial))
            )
        )
        candidate_density = np.polynomial.polynomial.polyval(
            candidates, n_poly
        ) * n_scale
        candidate_current = np.polynomial.polynomial.polyval(
            candidates, j_poly
        ) * j_scale
        candidate_energy = np.polynomial.polynomial.polyval(
            candidates, e_poly
        ) * e_scale
        maximum = max(
            maximum,
            float(
                np.max(
                    _net_flux_speed_ratio_values(
                        candidate_density,
                        candidate_current,
                        candidate_energy,
                    )
                )
            ),
        )
    return maximum


def fluid_applicability_diagnostics(
    profiles: Sequence[tuple[ProfileSource, pd.DataFrame]],
    *,
    gas_temperature_K: float,
    pressure_Pa: float,
    operational_ratio_guide: float = NET_FLUX_SPEED_RATIO_GUIDE,
) -> pd.DataFrame:
    """Quantify two checkable premises of the drift-diffusion model.

    ``n_e/N`` is only an electron-fraction proxy for weak ionization. The
    flux-speed ratio contains drift, diffusion, and boundary flux and is
    therefore a heuristic scale-separation diagnostic, not a validity test.
    """

    neutral_density = pressure_Pa / (
        BOLTZMANN_CONSTANT_J_K * gas_temperature_K
    )
    rows: list[dict[str, Any]] = []
    for source, frame in profiles:
        x, density = _deduplicated_xy(frame, "electron_density")
        current_x, current = _deduplicated_xy(
            frame, "electron_current_density"
        )
        energy_x, energy = _deduplicated_xy(frame, "mean_electron_energy")
        if not np.array_equal(current_x, x):
            current = np.interp(x, current_x, current)
        if not np.array_equal(energy_x, x):
            energy = np.interp(x, energy_x, energy)
        domain_length = float(x[-1] - x[0])
        above_guide = _net_flux_ratio_measure_above(
            x,
            density,
            current,
            energy,
            operational_ratio_guide,
        )
        above_unity = _net_flux_ratio_measure_above(
            x,
            density,
            current,
            energy,
            1.0,
        )
        rows.append(
            {
                "path": source.path_id,
                "path_label": source.label,
                "mesh_elements": source.mesh_elements,
                "gas_temperature_K": gas_temperature_K,
                "pressure_Pa": pressure_Pa,
                "pressure_mTorr": pressure_Pa / 0.133322,
                "neutral_number_density_m3": neutral_density,
                "max_electron_density_m3": float(np.max(density)),
                "max_electron_to_neutral_ratio": float(
                    np.max(density) / neutral_density
                ),
                "net_flux_speed_ratio_definition": (
                    "|J_e|/(e n_e sqrt(2 e mean_energy/m_e))"
                ),
                "net_flux_speed_ratio_operational_guide": (
                    operational_ratio_guide
                ),
                "max_net_flux_speed_ratio": _net_flux_ratio_piecewise_max(
                    x, density, current, energy
                ),
                "spatial_fraction_net_flux_speed_ratio_above_guide": (
                    above_guide / domain_length
                ),
                "spatial_fraction_net_flux_speed_ratio_above_unity": (
                    above_unity / domain_length
                ),
                "interpretation": (
                    "weak-ionization proxy plus heuristic boundary/fluid "
                    "scale-separation flag; not a hard validity criterion"
                ),
            }
        )
    return pd.DataFrame(rows)


def operating_point_diagnostics(
    profiles: Sequence[tuple[ProfileSource, pd.DataFrame]],
    e_over_n_unit: str,
    *,
    gas_temperature_K: float,
    pressure_Pa: float,
    ballast_resistance_ohm: float = 10_000.0,
    high_field_threshold_Td: float = TYPICAL_DRIFT_DIFFUSION_GUIDANCE_TD,
) -> pd.DataFrame:
    """Extract circuit and signed-potential diagnostics from archived profiles."""

    neutral_density = pressure_Pa / (
        BOLTZMANN_CONSTANT_J_K * gas_temperature_K
    )
    rows: list[dict[str, Any]] = []
    for source, frame in profiles:
        x, potential = _deduplicated_xy(frame, "electric_potential")
        en_x, en_raw = _deduplicated_xy(frame, "E_over_N")
        en_Td = convert_e_over_n_to_td(en_raw, e_over_n_unit)
        if not np.array_equal(en_x, x):
            en_Td = np.interp(x, en_x, en_Td)
        source_voltage = float(
            pd.to_numeric(frame["applied_voltage"], errors="raise").iloc[0]
        )
        gap_voltage = float(potential[-1] - potential[0])
        maximum_index = int(np.argmax(potential))
        total_variation_from_field = (
            piecewise_linear_integral(x, en_Td)
            * neutral_density
            * 1.0e-21
        )
        high_field_share = _condition_integral_fraction(
            x,
            en_Td,
            en_Td,
            "above",
            high_field_threshold_Td,
        )
        ballast_current = (
            source_voltage - gap_voltage
        ) / ballast_resistance_ohm
        rows.append(
            {
                "path": source.path_id,
                "path_label": source.label,
                "mesh_elements": source.mesh_elements,
                "source_voltage_V": source_voltage,
                "cathode_potential_V": float(potential[0]),
                "anode_terminal_potential_V": float(potential[-1]),
                "gap_endpoint_voltage_V": gap_voltage,
                "ballast_drop_V": source_voltage - gap_voltage,
                "ballast_resistance_ohm": ballast_resistance_ohm,
                "ballast_current_A": ballast_current,
                "internal_max_potential_V": float(potential[maximum_index]),
                "internal_max_potential_x_m": float(x[maximum_index]),
                "anode_side_potential_overshoot_V": float(
                    potential[maximum_index] - potential[-1]
                ),
                "axial_absolute_potential_variation_V": (
                    total_variation_from_field
                ),
                "high_field_threshold_Td": high_field_threshold_Td,
                "high_field_absolute_potential_variation_share": (
                    high_field_share
                ),
                "field_reversal_interpretation": (
                    "internal potential maximum is a reversal proxy; signed "
                    "electric field was not exported"
                ),
            }
        )
    return pd.DataFrame(rows)


def townsend_source_representation_audit(
    archived_export: pd.DataFrame,
    intended_bundle: Path,
    *,
    gas_temperature_K: float,
    pressure_Pa: float,
) -> pd.DataFrame:
    """Compare archived Townsend-flux sources with a rate-form counterfactual."""

    rates = pd.read_csv(intended_bundle / "rates_vs_mean_energy.csv")
    neutral_density = pressure_Pa / (
        BOLTZMANN_CONSTANT_J_K * gas_temperature_K
    )
    x, mean_energy = _deduplicated_xy(
        archived_export, "mean_electron_energy"
    )
    _, density = _deduplicated_xy(archived_export, "electron_density")
    _, electron_current = _deduplicated_xy(
        archived_export, "electron_current_density"
    )
    specs = (
        (
            "excitation",
            "eir2 direct excitation",
            "excitation_source",
            "excitation_townsend",
        ),
        (
            "ionization",
            "eir4 direct ionization",
            "ionization_source",
            "ionization_townsend",
        ),
    )
    rows: list[dict[str, Any]] = []
    for process_type, label, source_column, townsend_column in specs:
        process = (
            rates.loc[rates["process_type"] == process_type]
            .groupby("mean_energy_eV", as_index=False)[
                [
                    "mixture_weighted_rate_m3_s",
                    "mixture_weighted_reduced_townsend_m2",
                ]
            ]
            .sum()
            .sort_values("mean_energy_eV")
        )
        rate = np.interp(
            mean_energy,
            process["mean_energy_eV"],
            process["mixture_weighted_rate_m3_s"],
        )
        exported_townsend = archived_export[townsend_column].to_numpy(float)
        observed_source = archived_export[source_column].to_numpy(float)
        reconstructed_townsend_source = (
            exported_townsend
            * neutral_density
            * np.abs(electron_current)
            / ELEMENTARY_CHARGE_C
        )
        counterfactual_rate_source = (
            neutral_density * density * rate
        )
        observed_integral = piecewise_linear_integral(x, observed_source)
        reconstructed_integral = piecewise_linear_integral(
            x, reconstructed_townsend_source
        )
        counterfactual_integral = piecewise_linear_integral(
            x, counterfactual_rate_source
        )
        rows.append(
            {
                "process_type": process_type,
                "channel_label": label,
                "observed_townsend_source_integral_m2_s": observed_integral,
                "reconstructed_townsend_source_integral_m2_s": (
                    reconstructed_integral
                ),
                "townsend_reconstruction_relative_difference": (
                    abs(reconstructed_integral - observed_integral)
                    / abs(observed_integral)
                    if observed_integral != 0.0
                    else np.nan
                ),
                "counterfactual_rate_source_integral_m2_s": (
                    counterfactual_integral
                ),
                "counterfactual_rate_over_observed_townsend": (
                    counterfactual_integral / observed_integral
                    if observed_integral != 0.0
                    else np.nan
                ),
                "townsend_source_definition": (
                    "(alpha/N) N |J_e|/e"
                ),
                "counterfactual_rate_source_definition": "N n_e k",
                "interpretation": (
                    "postprocess representation sensitivity; not a matched "
                    "COMSOL rerun and not evidence that rate form is superior"
                ),
            }
        )
    return pd.DataFrame(rows)


def partial_electron_energy_audit(
    profiles: Sequence[tuple[ProfileSource, pd.DataFrame]],
) -> pd.DataFrame:
    """Compute a deliberately incomplete profile-derived energy audit."""

    rows: list[dict[str, Any]] = []
    for source, frame in profiles:
        x, potential = _deduplicated_xy(frame, "electric_potential")
        current_x, electron_current = _deduplicated_xy(
            frame, "electron_current_density"
        )
        excitation_x, excitation = _deduplicated_xy(
            frame, "excitation_source"
        )
        ionization_x, ionization = _deduplicated_xy(
            frame, "ionization_source"
        )
        if not np.array_equal(current_x, x):
            electron_current = np.interp(x, current_x, electron_current)
        if not np.array_equal(excitation_x, x):
            excitation = np.interp(x, excitation_x, excitation)
        if not np.array_equal(ionization_x, x):
            ionization = np.interp(x, ionization_x, ionization)
        field_power = float(
            np.sum(
                -0.5
                * (electron_current[:-1] + electron_current[1:])
                * np.diff(potential)
            )
        )
        excitation_loss = (
            piecewise_linear_integral(x, excitation)
            * 11.5
            * ELEMENTARY_CHARGE_C
        )
        ionization_loss = (
            piecewise_linear_integral(x, ionization)
            * 15.8
            * ELEMENTARY_CHARGE_C
        )
        direct_loss = excitation_loss + ionization_loss
        rows.append(
            {
                "path": source.path_id,
                "path_label": source.label,
                "electron_field_power_signed_W_m2": field_power,
                "eir2_threshold_weighted_loss_W_m2": excitation_loss,
                "eir4_threshold_weighted_loss_W_m2": ionization_loss,
                "direct_eir2_eir4_loss_sum_W_m2": direct_loss,
                "direct_loss_to_field_power_ratio": (
                    direct_loss / field_power if field_power != 0.0 else np.nan
                ),
                "audit_scope": (
                    "partial only: excludes energy flux, elastic exchange, "
                    "superelastic, stepwise/Penning, secondary emission, "
                    "boundary terms, and time derivative"
                ),
            }
        )
    return pd.DataFrame(rows)


def tail_convergence_audit(
    database: Path,
    *,
    tail_probability_target: float,
    edge_to_peak_target: float,
    max_energy_limit_eV: float,
) -> pd.DataFrame:
    """Audit the solver diagnostics used to qualify the formal Swarm table.

    ``quality.csv`` verifies exported table normalization but does not by
    itself establish convergence of the energy-space solve.  The workflow
    database retains the per-case diagnostics needed for that separate gate.
    """

    query = """
        SELECT solver, e_over_n_Td, replicate, diagnostics_json
        FROM cases
        ORDER BY solver, e_over_n_Td, replicate
    """
    rows: list[dict[str, Any]] = []
    with sqlite3.connect(database) as connection:
        records = connection.execute(query).fetchall()
    if not records:
        raise ValueError(f"No solver cases found in {database}")
    for solver, e_over_n_Td, replicate, diagnostics_json in records:
        payload = json.loads(diagnostics_json)
        diagnostic = payload.get(str(solver), {})
        if not isinstance(diagnostic, Mapping):
            raise ValueError(
                f"Missing {solver!r} diagnostics at E/N={e_over_n_Td} Td"
            )
        converged = bool(diagnostic.get("converged", False))
        tail_probability = float(diagnostic.get("tail_probability", np.nan))
        edge_to_peak = float(diagnostic.get("edge_to_peak", np.nan))
        grid_max_eV = float(diagnostic.get("grid_max_eV", np.nan))
        tail_pass = bool(
            np.isfinite(tail_probability)
            and tail_probability <= tail_probability_target
        )
        edge_pass = bool(
            np.isfinite(edge_to_peak) and edge_to_peak <= edge_to_peak_target
        )
        grid_limit_hit = bool(
            np.isfinite(grid_max_eV)
            and grid_max_eV >= max_energy_limit_eV * (1.0 - 1.0e-9)
        )
        rows.append(
            {
                "solver": str(solver),
                "E_over_N_Td": float(e_over_n_Td),
                "replicate": int(replicate),
                "solver_converged": converged,
                "iterations": diagnostic.get("iterations"),
                "residual_L1": diagnostic.get("residual_L1"),
                "adaptive_cycles": diagnostic.get("adaptive_cycles"),
                "grid_max_eV": grid_max_eV,
                "max_energy_limit_eV": max_energy_limit_eV,
                "grid_limit_hit": grid_limit_hit,
                "tail_probability": tail_probability,
                "tail_probability_target": tail_probability_target,
                "tail_pass": tail_pass,
                "edge_to_peak": edge_to_peak,
                "edge_to_peak_target": edge_to_peak_target,
                "edge_pass": edge_pass,
                "formal_tail_grid_gate_pass": bool(
                    tail_pass and edge_pass and not grid_limit_hit
                ),
                "formal_solver_convergence_gate_pass": converged,
                "formal_swarm_case_gate_pass": bool(
                    converged and tail_pass and edge_pass and not grid_limit_hit
                ),
            }
        )
    return pd.DataFrame(rows)


def drift_diffusion_regime_metrics(
    profiles: Sequence[tuple[ProfileSource, pd.DataFrame]],
    e_over_n_unit: str,
    gas_temperature_K: float,
    pressure_Pa: float,
    threshold_Td: float = TYPICAL_DRIFT_DIFFUSION_GUIDANCE_TD,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, float]]:
    """Describe profile support above COMSOL's typical 500 Td guidance.

    The threshold is contextual guidance for drift-diffusion use, not a hard
    validity boundary or an activation criterion.
    """

    number_density_m3 = pressure_Pa / (
        BOLTZMANN_CONSTANT_J_K * gas_temperature_K
    )
    electric_field_threshold_V_m = number_density_m3 * threshold_Td * 1.0e-21
    metric_rows: list[dict[str, Any]] = []
    impact_rows: list[dict[str, Any]] = []
    direct_source_labels = {
        "excitation_source": "Direct excitation source (eir2)",
        "ionization_source": "Direct ionization source (eir4)",
    }
    for source, frame in profiles:
        x, en_raw = _deduplicated_xy(frame, "E_over_N")
        en_Td = convert_e_over_n_to_td(en_raw, e_over_n_unit)
        domain_length = float(x.max() - x.min())
        above_length = _piecewise_condition_measure(
            x, en_Td, "above", threshold_Td
        )
        metric_rows.append(
            {
                "path": source.path_id,
                "path_label": source.label,
                "mesh_elements": source.mesh_elements,
                "regime": "E_over_N_above_typical_drift_diffusion_guidance",
                "threshold_Td": threshold_Td,
                "threshold_interpretation": (
                    "COMSOL typical drift-diffusion guidance; contextual, "
                    "not a hard validity cutoff"
                ),
                "gas_temperature_K": gas_temperature_K,
                "pressure_Pa": pressure_Pa,
                "neutral_number_density_m3": number_density_m3,
                "equivalent_electric_field_threshold_V_m": (
                    electric_field_threshold_V_m
                ),
                "domain_length_m": domain_length,
                "spatial_length_m": above_length,
                "spatial_fraction": above_length / domain_length,
                "node_count": int(np.sum(en_Td > threshold_Td)),
                "node_total": int(len(en_Td)),
                "node_fraction": float(np.mean(en_Td > threshold_Td)),
            }
        )
        for field, label, unit in QUANTITIES:
            value_x, values = _deduplicated_xy(frame, field)
            if not np.array_equal(value_x, x):
                values = np.interp(x, value_x, values)
            if field == "E_over_N":
                values = convert_e_over_n_to_td(values, e_over_n_unit)
            impact_rows.append(
                {
                    "path": source.path_id,
                    "path_label": source.label,
                    "quantity": field,
                    "quantity_label": direct_source_labels.get(field, label),
                    "quantity_unit": unit,
                    "regime": (
                        "E_over_N_above_typical_drift_diffusion_guidance"
                    ),
                    "threshold_Td": threshold_Td,
                    "absolute_integral_share": _condition_integral_fraction(
                        x, en_Td, values, "above", threshold_Td
                    ),
                    "integration_definition": (
                        "absolute_integral_of_piecewise_linear_reconstruction_"
                        "split_at_field_threshold_and_value_zero_crossings"
                    ),
                }
            )
    context = {
        "threshold_Td": threshold_Td,
        "gas_temperature_K": gas_temperature_K,
        "pressure_Pa": pressure_Pa,
        "neutral_number_density_m3": number_density_m3,
        "equivalent_electric_field_threshold_V_m": (
            electric_field_threshold_V_m
        ),
    }
    return pd.DataFrame(metric_rows), pd.DataFrame(impact_rows), context


def historical_activation_audit(
    archived_export: pd.DataFrame,
    intended_bundle: Path,
    *,
    relative_error_threshold: float = 0.01,
    formulation: str = "LocalEnergyApproximationE",
) -> pd.DataFrame:
    """Compare intended historical tables with archived exported values.

    This is deliberately not called an activation verification: the archived
    execution has unresolved Java/class provenance.  In the archived LEA
    formulation the mean-energy(E/N) lookup is inactive, while transport and
    Townsend values can be checked for numerical consistency.
    """

    mean_energy_table = pd.read_csv(intended_bundle / "mean_energy_vs_en.csv")
    transport_table = pd.read_csv(
        intended_bundle / "transport_vs_mean_energy.csv"
    )
    rates_table = pd.read_csv(intended_bundle / "rates_vs_mean_energy.csv")

    def prepare_table(
        frame: pd.DataFrame,
        x_column: str,
        y_column: str,
        process_type: str | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        selected = frame
        if process_type is not None:
            selected = selected[selected["process_type"] == process_type]
            selected = selected.groupby(x_column, as_index=False)[y_column].sum()
        selected = (
            selected.loc[:, [x_column, y_column]]
            .dropna()
            .sort_values(x_column, kind="stable")
        )
        if selected[x_column].duplicated().any():
            selected = selected.groupby(x_column, as_index=False)[y_column].mean()
        return (
            selected[x_column].to_numpy(float),
            selected[y_column].to_numpy(float),
        )

    specs = (
        {
            "quantity": "mean_energy_from_E_over_N",
            "quantity_label": "Mean energy from E/N",
            "table": mean_energy_table,
            "table_file": "mean_energy_vs_en.csv",
            "x_column": "E_over_N_V_m2",
            "y_column": "mean_energy_eV",
            "query_column": "E_over_N",
            "export_column": "mean_electron_energy",
            "process_type": None,
            "formulation_activity": "inactive_in_archived_LEA",
        },
        {
            "quantity": "reduced_mobility_from_mean_energy",
            "quantity_label": "Reduced mobility from mean energy",
            "table": transport_table,
            "table_file": "transport_vs_mean_energy.csv",
            "x_column": "mean_energy_eV",
            "y_column": "reduced_mobility_m2_V_s_m3",
            "query_column": "mean_electron_energy",
            "export_column": "reduced_mobility",
            "process_type": None,
            "formulation_activity": (
                "active_transport_in_archived_LEA_path_provenance_unresolved"
            ),
        },
        {
            "quantity": "direct_excitation_townsend_from_mean_energy",
            "quantity_label": "Direct excitation Townsend (eir2) from mean energy",
            "table": rates_table,
            "table_file": "rates_vs_mean_energy.csv",
            "x_column": "mean_energy_eV",
            "y_column": "mixture_weighted_reduced_townsend_m2",
            "query_column": "mean_electron_energy",
            "export_column": "excitation_townsend",
            "process_type": "excitation",
            "formulation_activity": (
                "active_Townsend_in_archived_LEA_path_provenance_unresolved"
            ),
        },
        {
            "quantity": "direct_ionization_townsend_from_mean_energy",
            "quantity_label": "Direct ionization Townsend (eir4) from mean energy",
            "table": rates_table,
            "table_file": "rates_vs_mean_energy.csv",
            "x_column": "mean_energy_eV",
            "y_column": "mixture_weighted_reduced_townsend_m2",
            "query_column": "mean_electron_energy",
            "export_column": "ionization_townsend",
            "process_type": "ionization",
            "formulation_activity": (
                "active_Townsend_in_archived_LEA_path_provenance_unresolved"
            ),
        },
    )
    rows: list[dict[str, Any]] = []
    for spec in specs:
        for column in (spec["query_column"], spec["export_column"]):
            if column not in archived_export:
                raise ValueError(
                    f"Archived export lacks activation-audit column {column!r}"
                )
        table_x, table_y = prepare_table(
            spec["table"],
            spec["x_column"],
            spec["y_column"],
            spec["process_type"],
        )
        query = archived_export[spec["query_column"]].to_numpy(float)
        observed = archived_export[spec["export_column"]].to_numpy(float)
        intended = np.interp(query, table_x, table_y)
        denominator_floor = max(
            float(np.max(np.abs(intended))) * 1.0e-12,
            np.finfo(float).tiny,
        )
        relative_error = np.abs(observed - intended) / np.maximum(
            np.abs(intended), denominator_floor
        )
        above_threshold = relative_error > relative_error_threshold
        numerically_matching = not bool(np.any(above_threshold))
        if spec["formulation_activity"] == "inactive_in_archived_LEA":
            status = (
                "inactive_in_archived_LEA_expected_mismatch_"
                "activation_unproven_stale_class"
                if not numerically_matching
                else "inactive_in_archived_LEA_coincidental_numeric_match_"
                "activation_unproven_stale_class"
            )
        elif numerically_matching:
            status = (
                "numerically_consistent_but_activation_unproven_stale_class"
            )
        else:
            status = "mismatch_and_activation_unproven_stale_class"
        rows.append(
            {
                "quantity": spec["quantity"],
                "quantity_label": spec["quantity_label"],
                "intended_table": (
                    f"{relative_path(intended_bundle)}/{spec['table_file']}"
                ),
                "table_argument": spec["x_column"],
                "archived_query_column": spec["query_column"],
                "archived_export_column": spec["export_column"],
                "process_type": spec["process_type"] or "",
                "comsol_formulation": formulation,
                "formulation_activity": spec["formulation_activity"],
                "rows_evaluated": len(relative_error),
                "relative_error_denominator_floor": denominator_floor,
                "relative_error_median": float(np.median(relative_error)),
                "relative_error_p95": float(np.quantile(relative_error, 0.95)),
                "relative_error_max": float(np.max(relative_error)),
                "relative_error_threshold": relative_error_threshold,
                "points_above_1_percent": int(np.sum(above_threshold)),
                "fraction_above_1_percent": float(np.mean(above_threshold)),
                "numeric_match_within_1_percent": numerically_matching,
                "status": status,
                "audit_scope": (
                    "intended_historical_bundle_vs_archived_export; "
                    "not_activation_verification"
                ),
                "provenance_limitation": (
                    "stale-class risk prevents binding the archived export to "
                    "the intended bundle/class"
                ),
            }
        )
    return pd.DataFrame(rows)


def qa_profile(
    source: ProfileSource, frame: pd.DataFrame, expected_rows: int | None
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    def add(
        check: str,
        passed: bool,
        observed: Any,
        expected: str,
        severity: str,
        detail: str,
    ) -> None:
        rows.append(
            {
                "source_id": source.path_id,
                "source_path": relative_path(source.path),
                "check": check,
                "passed": bool(passed),
                "observed": observed,
                "expected": expected,
                "severity_if_failed": severity,
                "detail": detail,
            }
        )

    add(
        "row_count",
        expected_rows is None or len(frame) == expected_rows,
        len(frame),
        str(expected_rows) if expected_rows is not None else "not constrained",
        "medium",
        "Unexpected profile length can indicate a changed export or mesh.",
    )
    numeric = frame.select_dtypes(include=[np.number])
    nonfinite = int((~np.isfinite(numeric.to_numpy(float))).sum())
    add(
        "all_numeric_values_finite",
        nonfinite == 0,
        nonfinite,
        "0 non-finite cells",
        "critical",
        "Non-finite values invalidate spatial integrals and comparisons.",
    )
    duplicated_x = int(frame["x"].duplicated().sum())
    add(
        "x_unique",
        duplicated_x == 0,
        duplicated_x,
        "0 duplicate x values",
        "high",
        "Duplicate coordinates make integration weights ambiguous.",
    )
    nonincreasing = int(np.sum(np.diff(frame["x"].to_numpy(float)) <= 0))
    add(
        "x_strictly_increasing",
        nonincreasing == 0,
        nonincreasing,
        "0 non-increasing intervals",
        "high",
        "A strictly increasing spatial coordinate is required.",
    )
    for field in ("electron_density", "excitation_source", "ionization_source"):
        negative = int((frame[field] < 0).sum())
        add(
            f"{field}_nonnegative",
            negative == 0,
            negative,
            "0 negative rows",
            "high",
            "Negative density or source values are physically invalid for this report.",
        )
    for field in ("applied_voltage", "gas_pressure"):
        unique = frame[field].dropna().nunique()
        add(
            f"{field}_constant",
            unique == 1,
            unique,
            "1 distinct value",
            "high",
            "Benchmark conditions must be constant across the exported profile.",
        )
    return rows


def qa_sources(
    config: Mapping[str, Any],
    external_source: ProfileSource,
    external: pd.DataFrame,
    reference_source: ProfileSource,
    reference: pd.DataFrame,
    mesh_profiles: Sequence[tuple[ProfileSource, pd.DataFrame]],
    runtime: pd.DataFrame,
    runtime_summary: pd.DataFrame,
    current_summary: pd.DataFrame,
    activation_audit: pd.DataFrame,
    tail_audit: pd.DataFrame,
    fallback_used: bool,
    source_paths: Sequence[tuple[str, Path]],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    expected_rows = config.get("expected", {}).get("profile_rows")
    qa_rows = qa_profile(external_source, external, expected_rows)
    qa_rows.extend(qa_profile(reference_source, reference, expected_rows))
    voltage_match = math.isclose(
        float(external["applied_voltage"].iloc[0]),
        float(reference["applied_voltage"].iloc[0]),
        rel_tol=0.0,
        abs_tol=1.0e-12,
    )
    pressure_match = math.isclose(
        float(external["gas_pressure"].iloc[0]),
        float(reference["gas_pressure"].iloc[0]),
        rel_tol=0.0,
        abs_tol=1.0e-9,
    )
    for check, passed, ext_value, ref_value, unit in (
        (
            "comparison_voltage_match",
            voltage_match,
            float(external["applied_voltage"].iloc[0]),
            float(reference["applied_voltage"].iloc[0]),
            "V",
        ),
        (
            "comparison_pressure_match",
            pressure_match,
            float(external["gas_pressure"].iloc[0]),
            float(reference["gas_pressure"].iloc[0]),
            "Pa",
        ),
    ):
        qa_rows.append(
            {
                "source_id": "comparison",
                "source_path": (
                    f"{relative_path(external_source.path)} | "
                    f"{relative_path(reference_source.path)}"
                ),
                "check": check,
                "passed": passed,
                "observed": f"external={ext_value:g}; reference={ref_value:g} {unit}",
                "expected": "equal",
                "severity_if_failed": "critical",
                "detail": "Unmatched conditions confound the closure comparison.",
            }
        )
    condition_audit = config.get("condition_audit", {})
    if isinstance(condition_audit, Mapping):
        historical_match = bool(
            condition_audit.get("historical_closure_gas_conditions_match", False)
        )
        historical_swarm = condition_audit.get("historical_external_swarm", {})
        historical_comsol = condition_audit.get("historical_comsol_profiles", {})
        qa_rows.append(
            {
                "source_id": "condition_audit",
                "source_path": relative_path(config["_config_path"]),
                "check": "historical_closure_gas_conditions_match",
                "passed": historical_match,
                "observed": (
                    f"Swarm={historical_swarm}; COMSOL={historical_comsol}"
                ),
                "expected": "identical gas_temperature_K and pressure_Pa",
                "severity_if_failed": "high",
                "detail": (
                    "The 2026-07-14 profile difference cannot be attributed to "
                    "closure alone when the historical Swarm and COMSOL gas "
                    "conditions differ."
                ),
            }
        )
    activation_config = config.get("historical_activation_audit", {})
    if isinstance(activation_config, Mapping) and activation_config:
        activation_provenance_verified = bool(
            activation_config.get("provenance_verified", False)
        )
        mismatch_rows = int(
            (~activation_audit["numeric_match_within_1_percent"]).sum()
        )
        qa_rows.append(
            {
                "source_id": "historical_activation_audit",
                "source_path": str(activation_config.get("intended_bundle", "")),
                "check": "historical_activation_provenance_verified",
                "passed": activation_provenance_verified,
                "observed": (
                    f"verified={activation_provenance_verified}; "
                    f"numeric_mismatch_rows={mismatch_rows}; "
                    "archived formulation=LocalEnergyApproximationE"
                ),
                "expected": (
                    "clean class bound by hash to intended bundle and archived export"
                ),
                "severity_if_failed": "high",
                "detail": (
                    "This is an intended-bundle-versus-export audit only. "
                    "Stale-class risk prevents activation proof; the archived "
                    "LEA formulation also leaves the mean-energy(E/N) table inactive."
                ),
            }
        )
        provenance_checks = (
            (
                "historical_swarm_config_sha256",
                historical_swarm.get("evidence_path")
                if isinstance(historical_swarm, Mapping)
                else None,
                historical_swarm.get("source_config_sha256")
                if isinstance(historical_swarm, Mapping)
                else None,
            ),
            (
                "formal_executed_config_sha256",
                condition_audit.get("formal_swarm", {}).get("executed_config_path")
                if isinstance(condition_audit.get("formal_swarm"), Mapping)
                else None,
                condition_audit.get("formal_swarm", {}).get(
                    "executed_config_sha256"
                )
                if isinstance(condition_audit.get("formal_swarm"), Mapping)
                else None,
            ),
            (
                "portable_repro_config_sha256",
                condition_audit.get("formal_swarm", {}).get(
                    "portable_repro_config_path"
                )
                if isinstance(condition_audit.get("formal_swarm"), Mapping)
                else None,
                condition_audit.get("formal_swarm", {}).get(
                    "portable_repro_config_sha256"
                )
                if isinstance(condition_audit.get("formal_swarm"), Mapping)
                else None,
            ),
            (
                "formal_swarm_database_sha256",
                condition_audit.get("formal_swarm", {}).get("database_path")
                if isinstance(condition_audit.get("formal_swarm"), Mapping)
                else None,
                condition_audit.get("formal_swarm", {}).get("database_sha256")
                if isinstance(condition_audit.get("formal_swarm"), Mapping)
                else None,
            ),
        )
        for check_name, path_value, expected_hash in provenance_checks:
            evidence_path = resolve_path(str(path_value)) if path_value else None
            actual_hash = (
                sha256_file(evidence_path)
                if evidence_path is not None and evidence_path.is_file()
                else None
            )
            qa_rows.append(
                {
                    "source_id": "condition_provenance",
                    "source_path": (
                        relative_path(evidence_path)
                        if evidence_path is not None
                        else "missing"
                    ),
                    "check": check_name,
                    "passed": bool(
                        expected_hash
                        and actual_hash
                        and actual_hash == str(expected_hash).lower()
                    ),
                    "observed": actual_hash,
                    "expected": expected_hash,
                    "severity_if_failed": "critical",
                    "detail": (
                        "The packaged condition evidence must remain bound to "
                        "the exact reviewed source by SHA-256."
                    ),
                }
            )
        formal_verified = bool(
            condition_audit.get("formal_run_conditions_verified", False)
        )
        qa_rows.append(
            {
                "source_id": "condition_audit",
                "source_path": relative_path(config["_config_path"]),
                "check": "formal_run_conditions_verified",
                "passed": formal_verified,
                "observed": formal_verified,
                "expected": True,
                "severity_if_failed": "high",
                "detail": (
                    "Configured target equality is not execution evidence; the "
                    "clean COMSOL run must export and verify the requested conditions."
                ),
            }
        )
    expected_repetitions = int(
        config.get("expected", {}).get("runtime_repetitions_per_path", 3)
    )
    for scope in ("solve_only", "end_to_end"):
        for path_id in ("external_swarm", "builtin_reference"):
            selected = runtime_summary[
                (runtime_summary["scope"] == scope)
                & (runtime_summary["path"] == path_id)
            ]
            repetitions = int(selected["repetitions"].iloc[0]) if len(selected) else 0
            qa_rows.append(
                {
                    "source_id": "runtime",
                    "source_path": str(config["runtime"]),
                    "check": f"{scope}_{path_id}_repetitions",
                    "passed": repetitions >= expected_repetitions,
                    "observed": repetitions,
                    "expected": f">={expected_repetitions}",
                    "severity_if_failed": "high",
                    "detail": "Timing claims require independent process repetitions.",
                }
            )
    expected_meshes = {
        int(value) for value in config.get("expected", {}).get("mesh_elements", [])
    }
    raw_meshes = {
        source.mesh_elements for source, _ in mesh_profiles if source.mesh_elements
    }
    qa_rows.append(
        {
            "source_id": "mesh_profiles",
            "source_path": "; ".join(relative_path(item.path) for item, _ in mesh_profiles),
            "check": "raw_mesh_profiles_available",
            "passed": expected_meshes <= raw_meshes and not fallback_used,
            "observed": f"raw={sorted(raw_meshes)}; summary_fallback={fallback_used}",
            "expected": f"raw profiles for {sorted(expected_meshes)}",
            "severity_if_failed": "high",
            "detail": "Summary-only mesh statistics cannot independently reproduce edge-region sensitivity.",
        }
    )
    runtime_nonfinite = int((~np.isfinite(runtime["seconds"].to_numpy(float))).sum())
    qa_rows.append(
        {
            "source_id": "runtime",
            "source_path": str(config["runtime"]),
            "check": "runtime_finite_and_positive",
            "passed": runtime_nonfinite == 0 and bool((runtime["seconds"] > 0).all()),
            "observed": f"nonfinite={runtime_nonfinite}; nonpositive={int((runtime['seconds'] <= 0).sum())}",
            "expected": "all finite and >0",
            "severity_if_failed": "critical",
            "detail": "Invalid duration values invalidate runtime summaries.",
        }
    )
    tail_grid_failures = int(
        (~tail_audit["formal_tail_grid_gate_pass"]).sum()
    )
    solver_failures = int(
        (~tail_audit["formal_solver_convergence_gate_pass"]).sum()
    )
    qa_rows.append(
        {
            "source_id": "formal_swarm_database",
            "source_path": str(config.get("formal_swarm_database", "")),
            "check": "formal_swarm_tail_convergence",
            "passed": (
                tail_grid_failures == 0
                and solver_failures == 0
                and len(tail_audit) > 0
            ),
            "observed": (
                f"cases={len(tail_audit)}; tail/grid failures="
                f"{tail_grid_failures}; solver-convergence failures="
                f"{solver_failures}; "
                f"max E/N={tail_audit['E_over_N_Td'].max():g} Td"
            ),
            "expected": (
                "every case converged, tail<=target, edge/peak<=target, "
                "and energy-grid limit not reached"
            ),
            "severity_if_failed": "critical",
            "detail": (
                "Export normalization alone does not qualify an EEDF table; "
                "per-case solver and tail diagnostics must pass."
            ),
        }
    )
    qa = pd.DataFrame(qa_rows)
    hash_rows: list[dict[str, Any]] = []
    seen: set[Path] = set()
    for source_id, path in source_paths:
        resolved = path.resolve()
        if resolved in seen or not resolved.is_file():
            continue
        seen.add(resolved)
        hash_rows.append(
            {
                "source_id": source_id,
                "path": relative_path(resolved),
                "bytes": resolved.stat().st_size,
                "sha256": sha256_file(resolved),
            }
        )
    hashes = pd.DataFrame(hash_rows)
    critical_or_high = qa["severity_if_failed"].isin(["critical", "high"])
    blocking_failures = qa[critical_or_high & ~qa["passed"]]
    status = {
        "data_quality_status": (
            "pass" if blocking_failures.empty else "share_with_caveats"
        ),
        "checks": int(len(qa)),
        "passed_checks": int(qa["passed"].sum()),
        "failed_checks": int((~qa["passed"]).sum()),
        "blocking_failures": blocking_failures["check"].tolist(),
        "raw_profile_files": len({source.path for source, _ in mesh_profiles}),
        "runtime_rows": len(runtime),
        "current_summary_rows": len(current_summary),
        "historical_activation_audit_rows": len(activation_audit),
    }
    return qa, hashes, status


def set_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.0,
            "axes.titlesize": 10.5,
            "axes.labelsize": 9.0,
            "axes.edgecolor": CHARCOAL,
            "axes.labelcolor": CHARCOAL,
            "xtick.color": CHARCOAL,
            "ytick.color": CHARCOAL,
            "text.color": CHARCOAL,
            "axes.facecolor": WHITE,
            "figure.facecolor": WHITE,
            "grid.color": LIGHT_GREY,
            "grid.linewidth": 0.6,
            "axes.linewidth": 0.8,
            "legend.frameon": False,
            "savefig.facecolor": WHITE,
        }
    )


def panel_label(axis: plt.Axes, label: str) -> None:
    axis.text(
        -0.11,
        1.07,
        label,
        transform=axis.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
    )


def finish_figure(
    figure: plt.Figure,
    output_base: Path,
    source_note: str,
    *,
    tight_rect: tuple[float, float, float, float] = (0.0, 0.035, 1.0, 0.97),
) -> list[Path]:
    figure.text(0.01, 0.008, source_note, ha="left", va="bottom", fontsize=6.8, color=GREY)
    figure.tight_layout(rect=tight_rect)
    paths = [output_base.with_suffix(".png"), output_base.with_suffix(".pdf")]
    figure.savefig(paths[0], dpi=320, bbox_inches="tight")
    figure.savefig(paths[1], bbox_inches="tight")
    plt.close(figure)
    return paths


def figure_workflow(
    figure_dir: Path, evidence_status: str
) -> tuple[list[Path], dict[str, Any]]:
    figure, axis = plt.subplots(figsize=(11.2, 6.2))
    axis.set_xlim(0, 12)
    axis.set_ylim(0, 8)
    axis.axis("off")
    nodes = [
        (0.4, 5.6, 2.0, 1.2, "schema v2 YAML\nsolver: two_term", BLUE_OPEN, BLUE),
        (3.0, 5.6, 2.0, 1.2, "SQLite\nSwarm results", BLUE_OPEN, BLUE),
        (5.6, 5.6, 2.0, 1.2, "Canonical tables\n+ manifest hashes", BLUE_OPEN, BLUE),
        (
            8.2,
            5.6,
            2.0,
            1.2,
            "COMSOL bundle\n7 lookup quantities",
            BLUE_OPEN,
            BLUE,
        ),
        (8.2, 2.6, 2.0, 1.2, "MPH apply\nmode + property verify", "#F0F2F4", CHARCOAL),
        (5.6, 2.6, 2.0, 1.2, "COMSOL spatial PDEs\nwalls, species, field", "#F0F2F4", CHARCOAL),
        (3.0, 2.6, 2.0, 1.2, "Voltage continuation\nsolve", "#F0F2F4", CHARCOAL),
        (0.4, 2.6, 2.0, 1.2, "Profiles + timing\ncomparison + QA", "#F0F2F4", CHARCOAL),
    ]
    for x, y, width, height, text_value, fill, edge in nodes:
        axis.add_patch(
            FancyBboxPatch(
                (x, y),
                width,
                height,
                boxstyle="round,pad=0.08,rounding_size=0.08",
                facecolor=fill,
                edgecolor=edge,
                linewidth=1.4,
            )
        )
        axis.text(x + width / 2, y + height / 2, text_value, ha="center", va="center")
    arrow_pairs = [
        ((2.4, 6.2), (3.0, 6.2), "data"),
        ((5.0, 6.2), (5.6, 6.2), "data"),
        ((7.6, 6.2), (8.2, 6.2), "data"),
        ((9.2, 5.6), (9.2, 3.8), "apply"),
        ((8.2, 3.2), (7.6, 3.2), "tables"),
        ((5.6, 3.2), (5.0, 3.2), "solve"),
        ((3.0, 3.2), (2.4, 3.2), "export"),
    ]
    for start, end, label in arrow_pairs:
        axis.add_patch(
            FancyArrowPatch(
                start,
                end,
                arrowstyle="-|>",
                mutation_scale=12,
                linewidth=1.2,
                color=CHARCOAL,
            )
        )
        axis.text(
            (start[0] + end[0]) / 2,
            (start[1] + end[1]) / 2 + 0.18,
            label,
            ha="center",
            fontsize=7.5,
            color=GREY,
        )
    seven = [
        "mean energy",
        "reduced mobility",
        "longitudinal reduced diffusion",
        "energy mobility",
        "energy diffusion",
        "excitation Townsend",
        "ionization Townsend",
    ]
    axis.text(
        10.55,
        6.2,
        "Written lookup quantities",
        color=BLUE,
        fontweight="bold",
        va="top",
    )
    axis.text(
        10.55,
        5.85,
        "\n".join(f"{index}. {name}" for index, name in enumerate(seven, 1)),
        fontsize=7.7,
        va="top",
        linespacing=1.35,
    )
    axis.text(
        6.0,
        1.28,
        "Selected-channel hybrid boundary: Swarm supplies seven lookup quantities "
        "(seven X-Y pairs; 14 property arrays); COMSOL retains the spatial model.",
        ha="center",
        fontsize=9,
        fontweight="bold",
    )
    axis.text(
        6.0,
        0.82,
        "Archived MPH uses LEA: six transport/Townsend quantities are candidates; "
        "mean-energy(E/N) is inactive. Residual reaction channels remain under "
        "COMSOL/Maxwellian treatment; historical activation is not proven.",
        ha="center",
        fontsize=7.5,
        color=GREY,
    )
    axis.set_title(
        "Data flow and responsibility boundary", loc="left", fontsize=13, pad=12
    )
    source_note = (
        "Source: schema v2 configuration, bundle manifest, COMSOL mapping and "
        f"Java generator. Evidence status: {evidence_status}."
    )
    paths = finish_figure(
        figure,
        figure_dir / "fig01_workflow_boundary",
        source_note,
        tight_rect=(0.0, 0.04, 1.0, 0.94),
    )
    contract = {
        "analytical_question": "Which physics is supplied by Swarm and which remains in COMSOL?",
        "takeaway": (
            "Seven lookup quantities (seven X-Y pairs, 14 COMSOL arrays) cross "
            "a selected-channel hybrid interface; the active subset is "
            "formulation-dependent and the archived path is LEA."
        ),
        "family": "process flow",
        "variant": "responsibility-boundary diagram",
        "fields": seven,
        "palette": "blue for Swarm, charcoal for COMSOL; labels retain meaning in grayscale",
    }
    return paths, contract


def _shade_lookup_range(
    axis: plt.Axes, en_min: float, en_max: float, profile_min: float, profile_max: float
) -> None:
    overlap_min = max(en_min, profile_min)
    overlap_max = min(en_max, profile_max)
    if overlap_max > overlap_min:
        axis.axvspan(overlap_min, overlap_max, color=BLUE, alpha=0.07, lw=0)
    if profile_min < en_min:
        axis.axvspan(profile_min, en_min, color=GOLD, alpha=0.22, lw=0)
    if profile_max > en_max:
        axis.axvspan(en_max, profile_max, color=PINK, alpha=0.12, lw=0)
    axis.axvline(en_min, color=GOLD, linestyle=":", linewidth=1)
    axis.axvline(en_max, color=GREY, linestyle=":", linewidth=1)


def figure_lookup(
    figure_dir: Path,
    bundle: Path,
    clamp: pd.DataFrame,
    clamp_context: Mapping[str, float],
    evidence_status: str,
    profile_evidence_status: str,
    lookup_evidence_status: str,
) -> tuple[list[Path], dict[str, Any]]:
    transport_en_path = bundle / "transport_vs_en.csv"
    transport_mean_path = bundle / "transport_vs_mean_energy.csv"
    rates_mean_path = bundle / "rates_vs_mean_energy.csv"
    transport_en = pd.read_csv(transport_en_path)
    transport_mean = pd.read_csv(transport_mean_path)
    rates_mean = pd.read_csv(rates_mean_path)
    transport_en = (
        transport_en[transport_en["E_over_N_Td"] > 0]
        .sort_values("E_over_N_Td")
        .drop_duplicates("E_over_N_Td", keep="last")
    )
    positive_transport_mean = transport_mean[
        transport_mean["mean_energy_eV"] > 0
    ].copy()
    positive_rates_mean = rates_mean[rates_mean["mean_energy_eV"] > 0].copy()
    excitation = (
        positive_rates_mean[
            positive_rates_mean["process_type"] == "excitation"
        ]
        .groupby("mean_energy_eV", as_index=False)[
            "mixture_weighted_reduced_townsend_m2"
        ]
        .sum()
    )
    ionization = (
        positive_rates_mean[
            positive_rates_mean["process_type"] == "ionization"
        ]
        .groupby("mean_energy_eV", as_index=False)[
            "mixture_weighted_reduced_townsend_m2"
        ]
        .sum()
    )
    figure, axes = plt.subplots(2, 3, figsize=(12.0, 7.3), sharex=False)
    plots = [
        (
            "reduced_mobility_m2_V_s_m3",
            "Reduced mobility (LEA-active candidate)",
            r"m$^{-1}$ V$^{-1}$ s$^{-1}$",
        ),
        (
            "reduced_diffusion_L_m2_s_m3",
            "Longitudinal diffusion (LEA-active candidate)",
            r"m$^{-1}$ s$^{-1}$",
        ),
        (
            "reduced_electron_energy_mobility_m2_V_s_m3",
            "Energy mobility (LEA-active candidate)",
            r"m$^{-1}$ V$^{-1}$ s$^{-1}$",
        ),
        (
            "reduced_electron_energy_diffusion_m2_s_m3",
            "Energy diffusion (LEA-active candidate)",
            r"m$^{-1}$ s$^{-1}$",
        ),
    ]
    en_min = clamp_context["table_E_over_N_min_Td"]
    en_max = clamp_context["table_E_over_N_max_Td"]
    profile_min = clamp_context["profile_E_over_N_min_Td"]
    profile_max = clamp_context["profile_E_over_N_max_Td"]
    computed_energy_min = clamp_context["computed_mean_energy_min_eV"]
    computed_energy_max = clamp_context["computed_mean_energy_max_eV"]
    lookup_energy_min = clamp_context["lookup_mean_energy_min_eV"]
    profile_energy_min = clamp_context["profile_mean_energy_min_eV"]
    profile_energy_max = clamp_context["profile_mean_energy_max_eV"]
    boundary_policy_max = clamp_context[
        "boundary_policy_constant_max_mean_energy_eV"
    ]

    axis = axes.flat[0]
    axis.plot(
        transport_en["E_over_N_Td"],
        transport_en["mean_energy_eV"],
        color=BLUE,
        linewidth=1.7,
    )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_title("Mean energy lookup (inactive in archived LEA)", loc="left")
    axis.set_ylabel("eV")
    axis.set_xlabel("E/N (Td)")
    axis.grid(True, which="major", alpha=0.65)
    _shade_lookup_range(axis, en_min, en_max, profile_min, profile_max)
    panel_label(axis, "A")

    def shade_mean_energy_use(axis: plt.Axes) -> None:
        overlap_min = max(lookup_energy_min, profile_energy_min)
        overlap_max = min(computed_energy_max, profile_energy_max)
        if overlap_max > overlap_min:
            axis.axvspan(
                max(overlap_min, computed_energy_min),
                overlap_max,
                color=BLUE,
                alpha=0.07,
                lw=0,
            )
        policy_upper = min(profile_energy_max, computed_energy_min)
        if policy_upper > profile_energy_min:
            axis.axvspan(
                profile_energy_min,
                min(boundary_policy_max, policy_upper),
                color=GOLD,
                alpha=0.27,
                lw=0,
            )
            transition_lower = max(profile_energy_min, boundary_policy_max)
            if policy_upper > transition_lower:
                axis.axvspan(
                    transition_lower,
                    policy_upper,
                    color=GOLD_LIGHT,
                    alpha=0.35,
                    lw=0,
                )
        axis.axvline(
            boundary_policy_max,
            color=GOLD,
            linestyle=":",
            linewidth=0.9,
        )
        axis.axvline(
            computed_energy_min,
            color=GOLD,
            linestyle="--",
            linewidth=0.9,
        )

    for index, (field, title, unit) in enumerate(plots, start=1):
        axis = axes.flat[index]
        axis.plot(
            positive_transport_mean["mean_energy_eV"],
            positive_transport_mean[field],
            color=BLUE,
            linewidth=1.7,
        )
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_title(title, loc="left")
        axis.set_ylabel(unit)
        axis.set_xlabel("Mean electron energy (eV)")
        axis.grid(True, which="major", alpha=0.65)
        shade_mean_energy_use(axis)
        if index == 1:
            axis.legend(
                handles=[
                    Patch(
                        facecolor=BLUE,
                        alpha=0.07,
                        edgecolor="none",
                        label="Profile use within computed support",
                    ),
                    Patch(
                        facecolor=GOLD,
                        alpha=0.27,
                        edgecolor="none",
                        label="Constant/zero boundary policy",
                    ),
                    Patch(
                        facecolor=GOLD_LIGHT,
                        alpha=0.35,
                        edgecolor="none",
                        label="Below first computed Swarm point",
                    ),
                    Line2D(
                        [0],
                        [0],
                        color=GOLD,
                        linestyle=":",
                        linewidth=0.9,
                        label=f"Policy boundary {boundary_policy_max:g} eV",
                    ),
                    Line2D(
                        [0],
                        [0],
                        color=GOLD,
                        linestyle="--",
                        linewidth=0.9,
                        label=f"First computed point {computed_energy_min:.3g} eV",
                    ),
                ],
                loc="best",
                fontsize=5.7,
            )
        panel_label(axis, chr(ord("A") + index))
    axis = axes.flat[5]
    axis.plot(
        excitation["mean_energy_eV"],
        excitation["mixture_weighted_reduced_townsend_m2"],
        color=BLUE,
        linewidth=1.7,
        label="Direct excitation (eir2)",
    )
    axis.plot(
        ionization["mean_energy_eV"],
        ionization["mixture_weighted_reduced_townsend_m2"],
        color=CHARCOAL,
        linewidth=1.5,
        linestyle="--",
        label="Direct ionization (eir4)",
    )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_title("Reduced Townsend (LEA-active candidates)", loc="left")
    axis.set_ylabel(r"m$^2$")
    axis.set_xlabel("Mean electron energy (eV)")
    axis.legend(loc="best")
    axis.grid(True, which="major", alpha=0.65)
    shade_mean_energy_use(axis)
    panel_label(axis, "F")
    low_field_fraction = float(
        clamp.loc[
            clamp["metric"] == "E_over_N_below_table", "spatial_fraction"
        ].iloc[0]
    )
    below_computed_fraction = float(
        clamp.loc[
            clamp["metric"] == "mean_energy_below_computed_domain",
            "spatial_fraction",
        ].iloc[0]
    )
    constant_policy_fraction = float(
        clamp.loc[
            clamp["metric"]
            == "mean_energy_in_constant_boundary_policy_region",
            "spatial_fraction",
        ].iloc[0]
    )
    figure.suptitle(
        "Formulation-aware lookup support: LFA field relation versus LEA-active tables\n"
        f"A: inactive-in-LEA E/N coverage ({low_field_fraction:.3%} below "
        f"{en_min:.4g} Td); B–F: {below_computed_fraction:.3%} below first "
        f"computed Swarm energy {computed_energy_min:.4g} eV, including "
        f"{constant_policy_fraction:.3%} in the ≤{boundary_policy_max:.3g} eV "
        "constant/zero boundary-policy region",
        x=0.055,
        ha="left",
        fontsize=11.5,
    )
    source_note = (
        f"Lookup source: {relative_path(transport_en_path)}, "
        f"{relative_path(transport_mean_path)}, and "
        f"{relative_path(rates_mean_path)} ({lookup_evidence_status}); "
        f"used-range source: external COMSOL profile ({profile_evidence_status}). "
        "Zero Townsend floor rows are omitted on the log y-axis. "
        f"Overall report status: {evidence_status}."
    )
    paths = finish_figure(
        figure,
        figure_dir / "fig02_lookup_domain",
        source_note,
        tight_rect=(0.0, 0.045, 1.0, 0.90),
    )
    contract = {
        "analytical_question": (
            "Do the formulation-specific lookup arguments used by COMSOL remain "
            "inside computed Swarm support or rely on an explicit boundary policy?"
        ),
        "takeaway": (
            f"The inactive LEA E/N relation has {low_field_fraction:.3%} "
            "below its computed field floor, whereas the active LEA tables "
            f"query mean energies below computed support over "
            f"{below_computed_fraction:.3%} of domain length."
        ),
        "family": "ordered relationship",
        "variant": (
            "one LFA-axis panel and five LEA-axis log-log lookup panels with "
            "computed-support and boundary-policy bands"
        ),
        "fields": ["mean_energy_eV"]
        + [field for field, _, _ in plots]
        + ["excitation reduced Townsend", "ionization reduced Townsend"],
        "palette": (
            "single blue root, charcoal comparison, gold constant-policy "
            "region, pale-gold transition"
        ),
        "data_vintages": {
            "lookup_tables": lookup_evidence_status,
            "COMSOL_used_range": profile_evidence_status,
        },
    }
    return paths, contract


def figure_profiles(
    figure_dir: Path,
    external: pd.DataFrame,
    reference: pd.DataFrame,
    external_label: str,
    reference_label: str,
    e_over_n_unit: str,
    evidence_status: str,
) -> tuple[list[Path], dict[str, Any]]:
    figure, axes = plt.subplots(2, 2, figsize=(11.4, 7.4), sharex=True)
    fields = [
        ("electron_density", "Electron density", r"m$^{-3}$", True),
        ("mean_electron_energy", "Mean electron energy", "eV", False),
        ("electric_potential", "Electric potential", "V", False),
        ("E_over_N", "Reduced electric field", "Td", True),
    ]
    for index, (field, title, unit, log_scale) in enumerate(fields):
        axis = axes.flat[index]
        x, ext, ref = common_grid(external, reference, field)
        if field == "E_over_N":
            ext = convert_e_over_n_to_td(ext, e_over_n_unit)
            ref = convert_e_over_n_to_td(ref, e_over_n_unit)
        axis.plot(x, ext, color=BLUE, linewidth=1.7, label=external_label)
        axis.plot(
            x,
            ref,
            color=CHARCOAL,
            linewidth=1.45,
            linestyle="--",
            label=reference_label,
        )
        if log_scale and np.all(ext > 0) and np.all(ref > 0):
            axis.set_yscale("log")
        axis.set_title(title, loc="left")
        axis.set_ylabel(unit)
        axis.grid(True, alpha=0.65)
        axis.axvline(float(x.min()), color=LIGHT_GREY, linewidth=0.8)
        axis.axvline(float(x.max()), color=LIGHT_GREY, linewidth=0.8)
        panel_label(axis, chr(ord("A") + index))
    axes.flat[0].text(
        0.01,
        0.04,
        "Cathode",
        transform=axes.flat[0].transAxes,
        fontsize=7,
        color=GREY,
        ha="left",
    )
    axes.flat[0].text(
        0.99,
        0.04,
        "Anode",
        transform=axes.flat[0].transAxes,
        fontsize=7,
        color=GREY,
        ha="right",
    )
    external_terminal = float(
        external.sort_values("x")["electric_potential"].iloc[-1]
    )
    reference_terminal = float(
        reference.sort_values("x")["electric_potential"].iloc[-1]
    )
    axes.flat[2].text(
        0.02,
        0.04,
        (
            f"Anode endpoint: external {external_terminal:.3f} V; "
            f"reference {reference_terminal:.3f} V"
        ),
        transform=axes.flat[2].transAxes,
        fontsize=6.8,
        color=GREY,
        va="bottom",
    )
    for axis in axes[-1, :]:
        axis.set_xlabel("Axial position x (m)")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.98, 0.985))
    voltage = float(external["applied_voltage"].iloc[0])
    pressure = float(external["gas_pressure"].iloc[0])
    figure.suptitle(
        "Final-time spatial state and field profiles\n"
        f"Source setting {voltage:g} V, {pressure:g} Pa; terminal gap voltages differ",
        x=0.065,
        ha="left",
        fontsize=12,
    )
    source_note = (
        "Source: external and built-in profile CSVs; linear interpolation on the "
        f"union of overlapping x grids. Evidence status: {evidence_status}."
    )
    paths = finish_figure(
        figure,
        figure_dir / "fig03_profile_state",
        source_note,
        tight_rect=(0.0, 0.045, 1.0, 0.91),
    )
    contract = {
        "analytical_question": "Where do state and field profiles agree or differ?",
        "takeaway": (
            "Profile shape and integrated magnitude must be read separately; "
            "the same 200 V source setting does not produce the same terminal "
            "gap voltage because ballast feedback differs."
        ),
        "family": "ordered relationship",
        "variant": "four aligned two-series line panels",
        "fields": [item[0] for item in fields],
        "palette": "blue solid external; charcoal dashed reference",
    }
    return paths, contract


def figure_sources_currents(
    figure_dir: Path,
    external: pd.DataFrame,
    reference: pd.DataFrame,
    external_label: str,
    reference_label: str,
    current_summary: pd.DataFrame,
    evidence_status: str,
) -> tuple[list[Path], dict[str, Any]]:
    figure, axes = plt.subplots(2, 3, figsize=(13.0, 7.5))
    line_fields = [
        (
            "excitation_source",
            "Direct excitation source (eir2)",
            r"m$^{-3}$ s$^{-1}$",
        ),
        (
            "ionization_source",
            "Direct ionization source (eir4)",
            r"m$^{-3}$ s$^{-1}$",
        ),
        ("electron_current_density", "Electron current density", r"A m$^{-2}$"),
        (
            "total_current_density",
            r"Electron + Ar$^+$ conductive current density",
            r"A m$^{-2}$",
        ),
    ]
    for index, (field, title, unit) in enumerate(line_fields):
        axis = axes.flat[index]
        x, ext, ref = common_grid(external, reference, field)
        axis.plot(x, ext, color=BLUE, linewidth=1.7, label=external_label)
        axis.plot(
            x,
            ref,
            color=CHARCOAL,
            linewidth=1.4,
            linestyle="--",
            label=reference_label,
        )
        if "current" in field:
            axis.axhline(0, color=GREY, linewidth=0.7)
            finite_nonzero = np.abs(
                np.concatenate((ext[np.isfinite(ext)], ref[np.isfinite(ref)]))
            )
            finite_nonzero = finite_nonzero[finite_nonzero > 0]
            if len(finite_nonzero):
                linthresh = max(
                    float(np.percentile(finite_nonzero, 25)),
                    float(np.max(finite_nonzero)) * 1.0e-8,
                )
                axis.set_yscale("symlog", linthresh=linthresh)
                axis.text(
                    0.02,
                    0.04,
                    "symlog: retains raw boundary nodes",
                    transform=axis.transAxes,
                    fontsize=6.2,
                    color=GREY,
                )
        axis.set_title(title, loc="left")
        axis.set_xlabel("Axial position x (m)")
        axis.set_ylabel(unit)
        axis.grid(True, alpha=0.65)
        panel_label(axis, chr(ord("A") + index))
    weighted_axis = axes.flat[4]
    weighted = current_summary[
        current_summary["metric_definition"]
        == "piecewise_linear_exact_spatial_population_RSD"
    ]
    weighted_domains = ["full_domain", "central_80_percent"]
    domain_positions = np.arange(len(weighted_domains), dtype=float)
    path_styles = (
        ("external_swarm", external_label, BLUE, "o", "-"),
        ("builtin_reference", reference_label, CHARCOAL, "s", "--"),
    )
    for path_index, (path_id, label, color, marker, line_style) in enumerate(
        path_styles
    ):
        selected = weighted[weighted["path"] == path_id]
        values = np.array(
            [
                float(
                    selected.loc[
                        selected["domain"] == domain,
                        "relative_standard_deviation",
                    ].iloc[0]
                )
                if not selected[selected["domain"] == domain].empty
                else np.nan
                for domain in weighted_domains
            ]
        )
        positions = domain_positions + (-0.10 if path_index == 0 else 0.10)
        weighted_axis.plot(
            positions,
            values,
            color=color,
            marker=marker,
            linestyle=line_style,
            linewidth=1.2,
            markersize=5.5,
            markerfacecolor=WHITE if path_index else color,
            label=label,
        )
        for position, value in zip(positions, values):
            if np.isfinite(value) and value > 0:
                weighted_axis.annotate(
                    f"{100.0 * value:.3g}%",
                    (position, value),
                    xytext=(0, 5),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=6.8,
                    color=color,
                )
    weighted_axis.set_yscale("log")
    weighted_axis.set_xticks(
        domain_positions, ["Full domain", "Central 80%\n(geometric window)"]
    )
    weighted_axis.set_xlabel("200 elements · raw profiles")
    weighted_axis.set_ylabel("Exact piecewise-linear current RSD")
    weighted_axis.set_title("Conductive-current spatial constancy", loc="left")
    weighted_axis.grid(True, axis="y", which="both", alpha=0.65)
    weighted_axis.legend(loc="best", fontsize=6.7)
    weighted_axis.yaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda value, _: f"{100.0 * value:g}%")
    )
    panel_label(weighted_axis, "E")

    legacy_axis = axes.flat[5]
    mesh_plot = current_summary[
        current_summary["path"] == "external_swarm"
    ].copy()
    legacy = mesh_plot[
        mesh_plot["metric_definition"] == "unweighted_discrete_population_RSD"
    ]
    legacy_meshes = sorted(legacy["mesh_elements"].unique())
    width = 0.34
    positions = np.arange(len(legacy_meshes))
    for offset, domain, label, color, hatch in (
        (-width / 2, "full_domain", "Full domain", GOLD, ""),
        (width / 2, "central_80_percent", "Central 80%", BLUE, "//"),
    ):
        values = [
            float(
                legacy.loc[
                    (legacy["mesh_elements"] == mesh)
                    & (legacy["domain"] == domain),
                    "relative_standard_deviation",
                ].iloc[0]
            )
            if not legacy[
                (legacy["mesh_elements"] == mesh)
                & (legacy["domain"] == domain)
            ].empty
            else np.nan
            for mesh in legacy_meshes
        ]
        legacy_axis.bar(
            positions + offset,
            values,
            width,
            color=color,
            edgecolor=CHARCOAL,
            linewidth=0.7,
            hatch=hatch,
            label=label,
        )
        for position, value in zip(positions + offset, values):
            if np.isfinite(value) and value > 0:
                legacy_axis.annotate(
                    f"{100.0 * value:.3g}%",
                    (position, value),
                    xytext=(0, 4),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=6.5,
                    color=CHARCOAL,
                )
    legacy_axis.set_xticks(positions, [str(mesh) for mesh in legacy_meshes])
    legacy_axis.set_xlabel("Mesh elements · summary only")
    legacy_axis.set_ylabel("Unweighted discrete current RSD")
    legacy_axis.set_title("Legacy mesh summary — not comparable to E", loc="left")
    legacy_axis.legend(loc="best")
    legacy_axis.set_yscale("log")
    legacy_axis.grid(True, axis="y", which="both", alpha=0.65)
    legacy_axis.yaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda value, _: f"{100.0 * value:g}%")
    )
    panel_label(legacy_axis, "F")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.985, 0.985))
    figure.suptitle(
        "Electron-impact sources, current profiles, and mesh sensitivity",
        x=0.055,
        ha="left",
        fontsize=12,
    )
    source_note = (
        "Source: external/reference profile CSVs and current_uniformity.csv. "
        f"Evidence status: {evidence_status}. Panel E uses exact piecewise-linear "
        "raw-profile "
        "RSD; panel F is an unweighted legacy summary and is not numerically "
        "comparable with panel E. Current panels retain raw endpoint nodes on "
        "a symmetric-log axis; no nodal-recovery study has classified the "
        "endpoint spikes as physical or numerical."
    )
    paths = finish_figure(
        figure,
        figure_dir / "fig04_sources_currents_mesh",
        source_note,
        tight_rect=(0.0, 0.045, 1.0, 0.93),
    )
    contract = {
        "analytical_question": (
            "How do direct-source and conductive-current profiles differ, and "
            "how sensitive is current constancy to path, geometric window, and mesh?"
        ),
        "takeaway": (
            "External and reference exact piecewise-linear RSD are compared in "
            "full and central-80% geometric windows; "
            "the saved evidence cannot support a weighted 200/400 comparison, "
            "and electron plus Ar+ conduction is not a complete terminal-current "
            "or displacement-current audit."
        ),
        "family": "ordered relationship and grouped comparison",
        "variant": (
            "four line panels, one raw exact-piecewise-linear RSD panel, and one "
            "separate legacy unweighted mesh-summary panel"
        ),
        "fields": [item[0] for item in line_fields]
        + [
            "mesh_elements",
            "domain",
            "relative_standard_deviation",
            "metric_definition",
        ],
        "palette": "blue external, charcoal reference, gold/blue legacy domains",
    }
    return paths, contract


def figure_metric_summary(
    figure_dir: Path, metrics: pd.DataFrame, evidence_status: str
) -> tuple[list[Path], dict[str, Any]]:
    plot = metrics.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(plot))
    figure, axes = plt.subplots(1, 2, figsize=(12.0, 5.9), sharey=True)
    axes[0].barh(
        y,
        plot["relative_L2_spatial_weighted"],
        color=BLUE,
        edgecolor=CHARCOAL,
        linewidth=0.7,
    )
    axes[0].set_yticks(y, plot["quantity_label"])
    axes[0].set_xlabel(r"$[\int(q_{ext}-q_{ref})^2dx/\int q_{ref}^2dx]^{1/2}$")
    axes[0].set_title("Exact piecewise-linear relative L2", loc="left")
    axes[0].grid(True, axis="x", alpha=0.65)
    panel_label(axes[0], "A")
    ratios = plot[
        "signed_integral_ratio_external_over_reference"
    ].to_numpy(float)
    axes[1].axvline(1.0, color=CHARCOAL, linewidth=1.1, linestyle="--", label="Unity")
    axes[1].hlines(y, 1.0, ratios, color=BLUE_LIGHT, linewidth=2.0)
    axes[1].scatter(
        ratios,
        y,
        s=42,
        facecolor=WHITE,
        edgecolor=BLUE,
        linewidth=1.5,
        zorder=3,
    )
    axes[1].set_xlabel(r"$\int q_{ext}dx / \int q_{ref}dx$")
    axes[1].set_title("Signed spatial integral ratio", loc="left")
    axes[1].grid(True, axis="x", alpha=0.65)
    panel_label(axes[1], "B")
    figure.suptitle(
        "Spatial comparison metrics by physical quantity\n"
        "Reference denominator is the COMSOL built-in Boltzmann result; it is not ground truth",
        x=0.08,
        ha="left",
        fontsize=12,
    )
    source_note = (
        "Source: comparison_metrics.csv; exact integration of union-grid "
        "piecewise-linear reconstructions. Current ratios are signed. "
        f"Evidence status: {evidence_status}."
    )
    paths = finish_figure(
        figure,
        figure_dir / "fig05_metric_summary",
        source_note,
        tight_rect=(0.0, 0.05, 1.0, 0.88),
    )
    contract = {
        "analytical_question": "Which quantities differ in pointwise shape/magnitude and spatial integral?",
        "takeaway": (
            "Exact piecewise-linear L2 and signed integral ratio answer "
            "different questions and are shown in separate panels."
        ),
        "family": "comparison and benchmark",
        "variant": "horizontal bars plus dot-to-unity plot",
        "fields": [
            "relative_L2_spatial_weighted",
            "signed_integral_ratio_external_over_reference",
        ],
        "palette": "single blue root with charcoal unity reference",
    }
    return paths, contract


def figure_spatial_robustness(
    figure_dir: Path,
    regional: pd.DataFrame,
    current_sensitivity: pd.DataFrame,
    evidence_status: str,
) -> tuple[list[Path], dict[str, Any]]:
    figure, axes = plt.subplots(1, 2, figsize=(13.0, 6.2))

    quantity_order = [
        field for field, _, _ in QUANTITIES
    ][::-1]
    quantity_labels = {
        field: label for field, label, _ in QUANTITIES
    }
    region_styles = (
        ("cathode_side_10_percent", "Cathode-side 10%", GOLD, ""),
        ("central_80_percent", "Central 80% (geometric)", BLUE, "//"),
        ("anode_side_10_percent", "Anode-side 10%", PINK, ".."),
    )
    y = np.arange(len(quantity_order))
    left = np.zeros(len(quantity_order), dtype=float)
    for region, label, color, hatch in region_styles:
        values = np.array(
            [
                float(
                    regional.loc[
                        (regional["quantity"] == quantity)
                        & (regional["region"] == region),
                        "squared_error_share",
                    ].iloc[0]
                )
                for quantity in quantity_order
            ]
        )
        axes[0].barh(
            y,
            values,
            left=left,
            color=color,
            edgecolor=CHARCOAL,
            linewidth=0.45,
            hatch=hatch,
            label=label,
        )
        left += values
    axes[0].set_yticks(
        y, [quantity_labels[quantity] for quantity in quantity_order]
    )
    axes[0].set_xlim(0.0, 1.0)
    axes[0].xaxis.set_major_formatter(
        matplotlib.ticker.PercentFormatter(xmax=1.0, decimals=0)
    )
    axes[0].set_xlabel("Share of full-domain squared-error integral")
    axes[0].set_title("Spatial attribution of comparison error", loc="left")
    axes[0].legend(loc="lower center", bbox_to_anchor=(0.5, 1.01), ncol=1)
    axes[0].grid(True, axis="x", alpha=0.55)
    panel_label(axes[0], "A")

    path_styles = (
        ("external_swarm", BLUE, "o", "-"),
        ("builtin_reference", CHARCOAL, "s", "--"),
    )
    for path_id, color, marker, line_style in path_styles:
        selected = current_sensitivity[
            current_sensitivity["path"] == path_id
        ].sort_values("central_retained_percent")
        if selected.empty:
            continue
        label = str(selected["path_label"].iloc[0])
        axes[1].plot(
            selected["central_retained_percent"],
            selected["piecewise_linear_exact_relative_standard_deviation"],
            color=color,
            marker=marker,
            linestyle=line_style,
            linewidth=1.45,
            markersize=4.5,
            markerfacecolor=WHITE if path_id == "builtin_reference" else color,
            label=label,
        )
        for retained in (80.0, 100.0):
            point = selected[
                np.isclose(selected["central_retained_percent"], retained)
            ]
            if point.empty:
                continue
            value = float(
                point[
                    "piecewise_linear_exact_relative_standard_deviation"
                ].iloc[0]
            )
            axes[1].annotate(
                f"{100.0 * value:.3g}%",
                (retained, value),
                xytext=(-3 if retained == 100.0 else 3, 6),
                textcoords="offset points",
                ha="right" if retained == 100.0 else "left",
                fontsize=6.8,
                color=color,
            )
    axes[1].axvline(80.0, color=GREY, linestyle=":", linewidth=0.9)
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Central geometric length retained (%)")
    axes[1].set_ylabel("Exact piecewise-linear conductive-current RSD")
    axes[1].set_title("Current constancy versus retained window", loc="left")
    axes[1].yaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda value, _: f"{100.0 * value:g}%")
    )
    axes[1].grid(True, which="both", alpha=0.55)
    axes[1].legend(loc="best")
    panel_label(axes[1], "B")

    figure.suptitle(
        "Spatial robustness diagnostics\n"
        "Error attribution and geometric-window sensitivity are distinct diagnostics",
        x=0.06,
        ha="left",
        fontsize=12,
    )
    source_note = (
        "Source: regional_error_attribution.csv and current_rsd_sensitivity.csv; "
        "exact products of union-grid piecewise-linear reconstructions. "
        f"Evidence status: {evidence_status}."
    )
    paths = finish_figure(
        figure,
        figure_dir / "fig07_spatial_robustness",
        source_note,
        tight_rect=(0.0, 0.05, 1.0, 0.88),
    )
    contract = {
        "analytical_question": (
            "Where does profile error accumulate, and is conductive-current "
            "RSD robust to the geometric central-window choice?"
        ),
        "takeaway": (
            "Different quantities concentrate error in different regions; "
            "current RSD changes sharply only when edge regions enter the window."
        ),
        "family": "composition and sensitivity",
        "variant": "100% stacked horizontal bars plus highlighted two-series line",
        "fields": [
            "quantity",
            "region",
            "squared_error_share",
            "central_retained_percent",
            "piecewise_linear_exact_relative_standard_deviation",
        ],
        "palette": (
            "gold cathode-side, blue central geometric window, pink anode-side; "
            "blue external and charcoal reference"
        ),
    }
    return paths, contract


def figure_activation_regime(
    figure_dir: Path,
    activation_audit: pd.DataFrame,
    regime_impact: pd.DataFrame,
    regime_context: Mapping[str, float],
    evidence_status: str,
) -> tuple[list[Path], dict[str, Any]]:
    figure, axes = plt.subplots(1, 2, figsize=(13.5, 6.4))

    activation = activation_audit.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(activation))
    median = activation["relative_error_median"].to_numpy(float)
    p95 = activation["relative_error_p95"].to_numpy(float)
    maximum = activation["relative_error_max"].to_numpy(float)
    positive_values = np.concatenate((median, p95, maximum))
    positive_values = positive_values[
        np.isfinite(positive_values) & (positive_values > 0)
    ]
    plot_floor = (
        max(float(positive_values.min()) * 0.5, 1.0e-18)
        if len(positive_values)
        else 1.0e-18
    )
    axes[0].hlines(
        y,
        np.maximum(median, plot_floor),
        np.maximum(maximum, plot_floor),
        color=BLUE_LIGHT,
        linewidth=2.0,
    )
    axes[0].scatter(
        np.maximum(median, plot_floor),
        y,
        s=32,
        color=BLUE,
        marker="o",
        label="Median",
        zorder=3,
    )
    axes[0].scatter(
        np.maximum(p95, plot_floor),
        y,
        s=36,
        facecolor=WHITE,
        edgecolor=BLUE,
        linewidth=1.2,
        marker="o",
        label="P95",
        zorder=3,
    )
    axes[0].scatter(
        np.maximum(maximum, plot_floor),
        y,
        s=38,
        color=CHARCOAL,
        marker="|",
        linewidth=1.5,
        label="Maximum",
        zorder=3,
    )
    axes[0].axvline(
        0.01,
        color=GOLD,
        linestyle="--",
        linewidth=1.2,
        label="1% threshold",
    )
    axes[0].set_xscale("log")
    axes[0].set_yticks(y, activation["quantity_label"])
    axes[0].set_xlabel("Relative error")
    axes[0].set_title(
        "Table/export consistency audit\n(not activation verification)",
        loc="left",
    )
    axes[0].grid(True, which="both", axis="x", alpha=0.55)
    axes[0].legend(loc="lower right", fontsize=7)
    for index, row in activation.iterrows():
        axes[0].annotate(
            f"{int(row['points_above_1_percent'])}/{int(row['rows_evaluated'])} >1%",
            (max(float(row["relative_error_max"]), plot_floor), index),
            xytext=(4, 0),
            textcoords="offset points",
            ha="left",
            va="center",
            fontsize=6.5,
            color=GREY,
        )
    panel_label(axes[0], "A")

    short_labels = {
        "electron_density": "Electron density",
        "mean_electron_energy": "Mean energy",
        "electric_potential": "Potential",
        "E_over_N": "E/N",
        "electron_current_density": "Electron current",
        "total_current_density": "Conductive current",
        "excitation_source": "Direct excitation (eir2)",
        "ionization_source": "Direct ionization (eir4)",
    }
    quantity_order = [field for field, _, _ in QUANTITIES][::-1]
    y = np.arange(len(quantity_order), dtype=float)
    for path_index, (path_id, color, marker) in enumerate(
        (
            ("external_swarm", BLUE, "o"),
            ("builtin_reference", CHARCOAL, "s"),
        )
    ):
        selected = regime_impact[regime_impact["path"] == path_id]
        values = np.array(
            [
                float(
                    selected.loc[
                        selected["quantity"] == quantity,
                        "absolute_integral_share",
                    ].iloc[0]
                )
                for quantity in quantity_order
            ]
        )
        positions = y + (-0.12 if path_index == 0 else 0.12)
        axes[1].scatter(
            values,
            positions,
            s=34,
            color=color if path_index == 0 else WHITE,
            edgecolor=color,
            marker=marker,
            linewidth=1.2,
            label=str(selected["path_label"].iloc[0]),
            zorder=3,
        )
        for value, position in zip(values, positions):
            axes[1].annotate(
                f"{100.0 * value:.3g}%",
                (value, position),
                xytext=(4, 0),
                textcoords="offset points",
                ha="left",
                va="center",
                fontsize=6.2,
                color=color,
            )
    axes[1].set_xscale("log")
    axes[1].set_yticks(
        y, [short_labels[quantity] for quantity in quantity_order]
    )
    axes[1].set_xlabel("Share of absolute full-domain spatial integral")
    axes[1].set_title(
        f"Contribution where E/N > {regime_context['threshold_Td']:g} Td",
        loc="left",
    )
    axes[1].xaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda value, _: f"{100.0 * value:g}%")
    )
    axes[1].grid(True, which="both", axis="x", alpha=0.55)
    axes[1].legend(loc="lower left", fontsize=7)
    panel_label(axes[1], "B")

    figure.suptitle(
        "Historical table/export consistency and high-field regime\n"
        "Archived MPH uses LEA: mean-energy(E/N) is inactive; stale class "
        "prevents activation proof",
        x=0.055,
        ha="left",
        fontsize=12,
    )
    source_note = (
        "Panel A: intended historical bundle versus archived export, not an "
        "activation verification. Panel B: 500 Td is typical COMSOL "
        "drift-diffusion guidance, not a hard validity cutoff. "
        f"Evidence status: {evidence_status}."
    )
    paths = finish_figure(
        figure,
        figure_dir / "fig08_activation_regime",
        source_note,
        tight_rect=(0.0, 0.055, 1.0, 0.87),
    )
    contract = {
        "analytical_question": (
            "Do archived exports numerically match the intended historical "
            "tables, and which spatial integrals are concentrated above 500 Td?"
        ),
        "takeaway": (
            "Some archived quantities are numerically consistent with the "
            "intended tables, but only mobility and eir2/eir4 have active-channel "
            "profile evidence; stale-class provenance prevents activation proof."
        ),
        "family": "uncertainty and regime comparison",
        "variant": (
            "log-scale median/P95/max error ranges plus paired log-scale "
            "absolute-integral-share dots"
        ),
        "fields": [
            "relative_error_median",
            "relative_error_p95",
            "relative_error_max",
            "points_above_1_percent",
            "absolute_integral_share",
        ],
        "palette": "blue intended/export audit and external; charcoal reference; gold threshold",
    }
    return paths, contract


def figure_fluid_applicability(
    figure_dir: Path,
    profiles: Sequence[tuple[ProfileSource, pd.DataFrame]],
    *,
    gas_temperature_K: float,
    pressure_Pa: float,
    evidence_status: str,
) -> tuple[list[Path], dict[str, Any]]:
    """Plot weak-ionization and net-flux scale-separation diagnostics."""

    neutral_density = pressure_Pa / (
        BOLTZMANN_CONSTANT_J_K * gas_temperature_K
    )
    figure, axes = plt.subplots(1, 2, figsize=(11.4, 4.6), sharex=True)
    styles = (
        (BLUE, "-", 1.8),
        (CHARCOAL, "--", 1.6),
    )
    maxima: list[str] = []
    for (source, frame), (color, linestyle, linewidth) in zip(
        profiles, styles, strict=False
    ):
        x, density = _deduplicated_xy(frame, "electron_density")
        current_x, current = _deduplicated_xy(
            frame, "electron_current_density"
        )
        energy_x, energy = _deduplicated_xy(frame, "mean_electron_energy")
        display_x = np.linspace(float(x.min()), float(x.max()), 5001)
        display_density = np.interp(display_x, x, density)
        display_current = np.interp(display_x, current_x, current)
        display_energy = np.interp(display_x, energy_x, energy)
        ratio = _net_flux_speed_ratio_values(
            display_density, display_current, display_energy
        )
        distance_mm = (display_x - float(x.min())) * 1.0e3
        label = (
            "External Swarm tables"
            if source.path_id == "external_swarm"
            else "Built-in Boltzmann reference"
        )
        axes[0].plot(
            distance_mm,
            display_density / neutral_density,
            color=color,
            linestyle=linestyle,
            linewidth=linewidth,
            label=label,
        )
        axes[1].plot(
            distance_mm,
            ratio,
            color=color,
            linestyle=linestyle,
            linewidth=linewidth,
            label=label,
        )
        maxima.append(
            f"{source.path_id}: max(nₑ/N)={np.max(density / neutral_density):.3g}"
        )
    axes[0].set_yscale("log")
    axes[0].set_title("Weak-ionization proxy", loc="left")
    axes[0].set_ylabel(r"Electron-to-neutral ratio $n_e/N$")
    axes[0].legend(loc="best")
    axes[0].grid(True, which="major", alpha=0.65)
    panel_label(axes[0], "A")
    axes[1].set_yscale("log")
    axes[1].set_title("Net electron-flux speed / energy speed", loc="left")
    axes[1].set_ylabel(
        r"$|J_e|/[e n_e\sqrt{2e\bar{\varepsilon}/m_e}]$"
    )
    axes[1].axhline(
        NET_FLUX_SPEED_RATIO_GUIDE,
        color=GOLD,
        linestyle=":",
        linewidth=1.2,
        label="0.1 operational guide",
    )
    axes[1].axhline(
        1.0,
        color=GREY,
        linestyle=":",
        linewidth=1.0,
        label="unity",
    )
    axes[1].legend(loc="best")
    axes[1].grid(True, which="major", alpha=0.65)
    panel_label(axes[1], "B")
    for axis in axes:
        axis.set_xlabel("Distance from cathode (mm)")
    figure.suptitle(
        "Checkable fluid-model premises and boundary-local scale separation",
        x=0.055,
        ha="left",
        fontsize=12,
    )
    figure.text(
        0.5,
        0.905,
        (
            f"p={pressure_Pa / 0.133322:.3g} mTorr, "
            f"Tg={gas_temperature_K:.5g} K, "
            f"N={neutral_density:.4g} m⁻³; "
            "0.1 is an audit guide, not a COMSOL validity cutoff"
        ),
        ha="center",
        fontsize=7.7,
        color=GREY,
    )
    source_note = (
        "Source: archived final-time electron density, mean-energy, and "
        "conductive electron-current profiles. The speed ratio includes drift, "
        "diffusion, and boundary flux; it is not a direct drift-velocity measure. "
        f"Evidence status: {evidence_status}."
    )
    paths = finish_figure(
        figure,
        figure_dir / "fig09_fluid_applicability",
        source_note,
        tight_rect=(0.0, 0.06, 1.0, 0.88),
    )
    contract = {
        "analytical_question": (
            "Which drift-diffusion premises can be checked directly from the "
            "archived profiles, and where does net flux cease to be small "
            "relative to the energy-equivalent electron speed?"
        ),
        "takeaway": (
            "Electron fractions are below 6e-6, supporting weak ionization, "
            "while the heuristic net-flux ratio becomes order unity only in "
            "short boundary regions for both paths."
        ),
        "family": "physical applicability diagnostic",
        "variant": "two-panel log-scale axial profiles",
        "fields": [
            "electron_density",
            "electron_current_density",
            "mean_electron_energy",
            "gas_temperature_K",
            "pressure_Pa",
        ],
        "palette": "blue external, charcoal reference, gold operational guide",
        "annotations": maxima,
    }
    return paths, contract


def figure_runtime(
    figure_dir: Path,
    runtime: pd.DataFrame,
    summary: pd.DataFrame,
    expected_repetitions: int,
    evidence_status: str,
) -> tuple[list[Path], dict[str, Any]]:
    normalized = normalize_runtime(runtime)
    figure, axes = plt.subplots(1, 2, figsize=(12.0, 5.5))
    scope_order = ["solve_only", "end_to_end"]
    path_order = ["external_swarm", "builtin_reference"]
    primary_repetitions = summary[
        summary["scope"].isin(scope_order)
        & summary["path"].isin(path_order)
    ]["repetitions"].to_numpy(int)
    single_observation = bool(
        len(primary_repetitions) and np.max(primary_repetitions) == 1
    )
    labels = {
        "solve_only": "COMSOL run stage",
        "end_to_end": "End to end",
        "external_swarm": "External",
        "builtin_reference": "Built-in ref.",
    }
    width = 0.34
    x = np.arange(len(scope_order))
    for path_index, path_id in enumerate(path_order):
        medians: list[float] = []
        errors_low: list[float] = []
        errors_high: list[float] = []
        for scope in scope_order:
            selected = summary[
                (summary["scope"] == scope) & (summary["path"] == path_id)
            ]
            if selected.empty:
                medians.append(np.nan)
                errors_low.append(0.0)
                errors_high.append(0.0)
            else:
                record = selected.iloc[0]
                medians.append(float(record["median_seconds"]))
                errors_low.append(
                    float(record["median_seconds"] - record["min_seconds"])
                )
                errors_high.append(
                    float(record["max_seconds"] - record["median_seconds"])
                )
        positions = x + (path_index - 0.5) * width
        color = BLUE if path_id == "external_swarm" else CHARCOAL
        hatch = "" if path_id == "external_swarm" else "//"
        if single_observation:
            axes[0].scatter(
                positions,
                medians,
                s=58,
                marker="o" if path_id == "external_swarm" else "s",
                facecolor=color if path_id == "external_swarm" else WHITE,
                edgecolor=color,
                linewidth=1.3,
                label=labels[path_id],
                zorder=4,
            )
        else:
            axes[0].bar(
                positions,
                medians,
                width,
                color=color,
                alpha=0.9,
                edgecolor=CHARCOAL,
                linewidth=0.7,
                hatch=hatch,
                yerr=np.vstack((errors_low, errors_high)),
                capsize=3,
                label=labels[path_id],
            )
            for scope_index, scope in enumerate(scope_order):
                observations = normalized[
                    (normalized["scope"] == scope)
                    & (normalized["path"] == path_id)
                ]["seconds"].to_numpy(float)
                offsets = np.linspace(-0.055, 0.055, max(len(observations), 1))
                axes[0].scatter(
                    positions[scope_index] + offsets[: len(observations)],
                    observations,
                    s=19,
                    facecolor=WHITE,
                    edgecolor=color,
                    linewidth=1.0,
                    zorder=4,
                )
    axes[0].set_xticks(x, [labels[item] for item in scope_order])
    axes[0].set_ylabel("Elapsed time (s)")
    axes[0].set_title(
        "Path runtime: run-stage and end-to-end observations", loc="left"
    )
    axes[0].grid(True, axis="y", alpha=0.65)
    axes[0].legend(loc="best")
    panel_label(axes[0], "A")
    stages = ["apply", "verify", "solve", "export"]
    stage_summary = summary[
        (summary["path"] == "external_swarm") & summary["scope"].isin(stages)
    ].set_index("scope")
    stage_median = [
        float(stage_summary.loc[stage, "median_seconds"])
        if stage in stage_summary.index
        else np.nan
        for stage in stages
    ]
    stage_min = [
        float(stage_summary.loc[stage, "min_seconds"])
        if stage in stage_summary.index
        else np.nan
        for stage in stages
    ]
    stage_max = [
        float(stage_summary.loc[stage, "max_seconds"])
        if stage in stage_summary.index
        else np.nan
        for stage in stages
    ]
    stage_x = np.arange(len(stages))
    if single_observation:
        axes[1].scatter(
            stage_x,
            stage_median,
            s=58,
            facecolor=[BLUE_OPEN, BLUE_LIGHT, BLUE, "#28547A"],
            edgecolor=CHARCOAL,
            linewidth=1.1,
            zorder=4,
        )
    else:
        axes[1].bar(
            stage_x,
            stage_median,
            color=[BLUE_OPEN, BLUE_LIGHT, BLUE, "#28547A"],
            edgecolor=CHARCOAL,
            linewidth=0.7,
            yerr=np.vstack(
                (
                    np.asarray(stage_median) - np.asarray(stage_min),
                    np.asarray(stage_max) - np.asarray(stage_median),
                )
            ),
            capsize=3,
        )
        for stage_index, stage in enumerate(stages):
            observations = normalized[
                (normalized["scope"] == stage)
                & (normalized["path"] == "external_swarm")
            ]["seconds"].to_numpy(float)
            offsets = np.linspace(-0.055, 0.055, max(len(observations), 1))
            axes[1].scatter(
                stage_x[stage_index] + offsets[: len(observations)],
                observations,
                s=19,
                facecolor=WHITE,
                edgecolor=CHARCOAL,
                linewidth=0.9,
                zorder=4,
            )
    axes[1].set_xticks(stage_x, [stage.title() for stage in stages])
    axes[1].set_ylabel("Elapsed time (s)")
    axes[1].set_title("External end-to-end stage breakdown", loc="left")
    axes[1].grid(True, axis="y", alpha=0.65)
    panel_label(axes[1], "B")
    count_values = summary[
        summary["scope"].isin(scope_order)
    ]["repetitions"].to_numpy(int)
    min_count = int(count_values.min()) if len(count_values) else 0
    qualifier = (
        "formal repetition gate met"
        if min_count >= expected_repetitions
        else f"provisional: minimum n={min_count}, required n={expected_repetitions}"
    )
    figure.suptitle(
        (
            "Runtime observations and external workflow breakdown\n"
            f"{qualifier}; solver architectures, DOF, and continuation are unmatched"
        ),
        x=0.055,
        ha="left",
        fontsize=12,
    )
    if single_observation:
        source_note = (
            "Source: runtime observations CSV. Markers are single observations "
            "(n=1); no distribution or uncertainty interval is shown. Historical "
            "42.914 s is a COMSOL run-stage time, not isolated nonlinear-solver "
            "time; the CSV does not bind the 370 s entry to an identified raw log. "
            f"Evidence status: {evidence_status}."
        )
    else:
        source_note = (
            "Source: runtime observations CSV. Bars are medians; whiskers span "
            "min–max; open circles are repetitions. Historical 42.914 s is a "
            "COMSOL run-stage time, not isolated nonlinear-solver time; the CSV "
            "does not bind the 370 s entry to an identified raw log. "
            f"Evidence status: {evidence_status}."
        )
    paths = finish_figure(
        figure,
        figure_dir / "fig06_runtime_breakdown",
        source_note,
        tight_rect=(0.0, 0.05, 1.0, 0.90),
    )
    contract = {
        "analytical_question": "Is the historical timing difference present in run-stage and complete-workflow scopes?",
        "takeaway": (
            f"Runtime distributions require n>={expected_repetitions}; this "
            f"dataset has minimum n={min_count}, and the two paths do not match "
            "DOF, energy-space dimension, or continuation sequence."
        ),
        "family": "uncertainty and benchmark",
        "variant": (
            "single-observation dots and stage dots"
            if single_observation
            else "grouped median/range bars with repetition points and stage breakdown"
        ),
        "fields": ["scope", "path", "seconds", "repetition"],
        "palette": "blue external, charcoal reference, ordered blue stage shades",
    }
    return paths, contract


def generate_figures(
    figure_dir: Path,
    evidence_status: str,
    profile_evidence_status: str,
    lookup_evidence_status: str,
    bundle: Path,
    external_source: ProfileSource,
    external: pd.DataFrame,
    reference_source: ProfileSource,
    reference: pd.DataFrame,
    e_over_n_unit: str,
    metrics: pd.DataFrame,
    clamp: pd.DataFrame,
    clamp_context: Mapping[str, float],
    current_summary: pd.DataFrame,
    regional: pd.DataFrame,
    current_sensitivity: pd.DataFrame,
    activation_audit: pd.DataFrame,
    regime_impact: pd.DataFrame,
    regime_context: Mapping[str, float],
    fluid_diagnostics: pd.DataFrame,
    runtime: pd.DataFrame,
    runtime_summary: pd.DataFrame,
    expected_repetitions: int,
    input_hashes: pd.DataFrame,
) -> list[dict[str, Any]]:
    figure_dir.mkdir(parents=True, exist_ok=True)
    set_plot_style()
    builders = [
        (
            "fig01_workflow_boundary",
            lambda: figure_workflow(
                figure_dir, "implementation_schema_v2_2026-07-31"
            ),
            ["config", "bundle_manifest", "analysis_code"],
        ),
        (
            "fig02_lookup_domain",
            lambda: figure_lookup(
                figure_dir,
                bundle,
                clamp,
                clamp_context,
                evidence_status,
                profile_evidence_status,
                lookup_evidence_status,
            ),
            ["bundle_transport", "bundle_rates", "external_profile"],
        ),
        (
            "fig03_profile_state",
            lambda: figure_profiles(
                figure_dir,
                external,
                reference,
                external_source.label,
                reference_source.label,
                e_over_n_unit,
                evidence_status,
            ),
            ["external_profile", "reference_profile"],
        ),
        (
            "fig04_sources_currents_mesh",
            lambda: figure_sources_currents(
                figure_dir,
                external,
                reference,
                external_source.label,
                reference_source.label,
                current_summary,
                evidence_status,
            ),
            ["external_profile", "reference_profile", "mesh_profiles"],
        ),
        (
            "fig05_metric_summary",
            lambda: figure_metric_summary(figure_dir, metrics, evidence_status),
            ["external_profile", "reference_profile"],
        ),
        (
            "fig06_runtime_breakdown",
            lambda: figure_runtime(
                figure_dir,
                runtime,
                runtime_summary,
                expected_repetitions,
                evidence_status,
            ),
            ["runtime"],
        ),
        (
            "fig07_spatial_robustness",
            lambda: figure_spatial_robustness(
                figure_dir,
                regional,
                current_sensitivity,
                evidence_status,
            ),
            ["external_profile", "reference_profile"],
        ),
        (
            "fig08_activation_regime",
            lambda: figure_activation_regime(
                figure_dir,
                activation_audit,
                regime_impact,
                regime_context,
                evidence_status,
            ),
            [
                "external_profile",
                "reference_profile",
                "historical_activation_mean_energy",
                "historical_activation_transport_mean_energy",
                "historical_activation_rates_mean_energy",
            ],
        ),
        (
            "fig09_fluid_applicability",
            lambda: figure_fluid_applicability(
                figure_dir,
                [
                    (external_source, external),
                    (reference_source, reference),
                ],
                gas_temperature_K=float(regime_context["gas_temperature_K"]),
                pressure_Pa=float(regime_context["pressure_Pa"]),
                evidence_status=evidence_status,
            ),
            ["external_profile", "reference_profile"],
        ),
    ]
    metadata: list[dict[str, Any]] = []
    hash_by_id = {
        row["source_id"]: {"path": row["path"], "sha256": row["sha256"]}
        for row in input_hashes.to_dict("records")
    }
    for figure_id, builder, source_ids in builders:
        paths, contract = builder()
        figure_evidence_status = (
            "implementation_schema_v2_2026-07-31"
            if figure_id == "fig01_workflow_boundary"
            else evidence_status
        )
        contract["output_footprint"] = "320-dpi PNG and vector PDF"
        contract["data_sufficiency"] = (
            "Rendered from every available row; missing formal repetitions or raw "
            "mesh profiles remain explicitly labelled."
        )
        metadata.append(
            {
                "figure_id": figure_id,
                "outputs": [
                    {
                        "path": relative_path(path),
                        "bytes": path.stat().st_size,
                        "sha256": sha256_file(path),
                    }
                    for path in paths
                ],
                "chart_contract": contract,
                "source_metadata": [
                    hash_by_id[source_id]
                    for source_id in source_ids
                    if source_id in hash_by_id
                ],
                "evidence_status": figure_evidence_status,
                "data_vintages": (
                    {
                        "lookup_tables": lookup_evidence_status,
                        "COMSOL_used_range": profile_evidence_status,
                    }
                    if figure_id == "fig02_lookup_domain"
                    else {
                        "implementation": "implementation_schema_v2_2026-07-31"
                    }
                    if figure_id == "fig01_workflow_boundary"
                    else {
                        "profile_and_runtime_evidence": profile_evidence_status
                    }
                ),
            }
        )
    return metadata


def run_analysis(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    evidence_status = str(config["evidence_status"])
    profile_evidence_status = str(
        config.get("profile_evidence_status", evidence_status)
    )
    lookup_evidence_status = str(
        config.get("lookup_evidence_status", evidence_status)
    )
    output_data = resolve_path(config["outputs"]["data_dir"])
    figure_dir = resolve_path(config["outputs"]["figure_dir"])
    output_data.mkdir(parents=True, exist_ok=True)
    external_source = profile_source(
        config["comparison"]["external"],
        path_id="external_swarm",
        default_label="External Swarm tables (archived LEA)",
    )
    reference_source = profile_source(
        config["comparison"]["reference"],
        path_id="builtin_reference",
        default_label="COMSOL built-in Boltzmann reference",
    )
    external = read_profile(external_source)
    reference = read_profile(reference_source)
    mesh_sources: list[tuple[ProfileSource, pd.DataFrame]] = []
    for index, item in enumerate(config.get("mesh_profiles", [])):
        source = profile_source(
            item,
            path_id=f"mesh_profile_{index}",
            default_label=f"Mesh profile {index + 1}",
        )
        mesh_sources.append((source, read_profile(source)))
    fallback_path = (
        resolve_path(config["mesh_summary_fallback"])
        if config.get("mesh_summary_fallback")
        else None
    )
    runtime_path = resolve_path(config["runtime"])
    runtime = pd.read_csv(runtime_path)
    runtime["seconds"] = pd.to_numeric(runtime["seconds"], errors="coerce")
    bundle = resolve_path(config["bundle"])
    e_over_n_unit = str(config.get("profile_units", {}).get("E_over_N", "V m^2"))
    metrics, aligned = comparison_metrics(external, reference, e_over_n_unit)
    metrics.insert(0, "evidence_status", evidence_status)
    regional = regional_error_attribution(
        external, reference, e_over_n_unit
    )
    regional.insert(0, "evidence_status", evidence_status)
    current_summary, fallback_used = current_uniformity(mesh_sources, fallback_path)
    current_summary.insert(0, "evidence_status", evidence_status)
    current_sensitivity = current_rsd_sensitivity(mesh_sources)
    current_sensitivity.insert(0, "evidence_status", evidence_status)
    runtime_summary, runtime_speedup = runtime_metrics(runtime)
    runtime_summary.insert(0, "evidence_status", evidence_status)
    if not runtime_speedup.empty:
        runtime_speedup.insert(0, "evidence_status", evidence_status)
    clamp, clamp_impact, clamp_context = clamp_metrics(
        external, bundle, e_over_n_unit
    )
    clamp.insert(0, "evidence_status", evidence_status)
    clamp_impact.insert(0, "evidence_status", evidence_status)
    historical_comsol_conditions = config.get("condition_audit", {}).get(
        "historical_comsol_profiles", {}
    )
    regime, regime_impact, regime_context = drift_diffusion_regime_metrics(
        mesh_sources,
        e_over_n_unit,
        gas_temperature_K=float(
            historical_comsol_conditions.get("gas_temperature_K", 293.15)
        ),
        pressure_Pa=float(
            historical_comsol_conditions.get(
                "pressure_Pa", external["gas_pressure"].iloc[0]
            )
        ),
    )
    regime.insert(0, "evidence_status", evidence_status)
    regime_impact.insert(0, "evidence_status", evidence_status)
    fluid_diagnostics = fluid_applicability_diagnostics(
        mesh_sources,
        gas_temperature_K=float(
            historical_comsol_conditions.get("gas_temperature_K", 293.15)
        ),
        pressure_Pa=float(
            historical_comsol_conditions.get(
                "pressure_Pa", external["gas_pressure"].iloc[0]
            )
        ),
    )
    fluid_diagnostics.insert(0, "evidence_status", evidence_status)
    operating_point = operating_point_diagnostics(
        mesh_sources,
        e_over_n_unit,
        gas_temperature_K=float(
            historical_comsol_conditions.get("gas_temperature_K", 293.15)
        ),
        pressure_Pa=float(
            historical_comsol_conditions.get(
                "pressure_Pa", external["gas_pressure"].iloc[0]
            )
        ),
    )
    operating_point.insert(0, "evidence_status", evidence_status)
    partial_energy = partial_electron_energy_audit(mesh_sources)
    partial_energy.insert(0, "evidence_status", evidence_status)
    tail_config = config.get("tail_convergence", {})
    formal_swarm_database = resolve_path(config["formal_swarm_database"])
    tail_audit = tail_convergence_audit(
        formal_swarm_database,
        tail_probability_target=float(
            tail_config.get("tail_probability_target", 1.0e-9)
        ),
        edge_to_peak_target=float(
            tail_config.get("edge_to_peak_target", 1.0e-10)
        ),
        max_energy_limit_eV=float(
            tail_config.get("max_energy_limit_eV", 20000.0)
        ),
    )
    tail_audit.insert(0, "evidence_status", lookup_evidence_status)
    activation_config = config.get("historical_activation_audit", {})
    intended_historical_bundle = resolve_path(
        activation_config["intended_bundle"]
    )
    activation_export_path = resolve_path(
        activation_config.get("archived_export", external_source.path)
    )
    activation_export_source = ProfileSource(
        path=activation_export_path,
        label="Archived external COMSOL export",
        mesh_elements=external_source.mesh_elements,
        path_id="historical_activation_export",
    )
    activation_export = read_profile(activation_export_source)
    activation_audit = historical_activation_audit(
        activation_export,
        intended_historical_bundle,
        relative_error_threshold=float(
            activation_config.get("relative_error_threshold", 0.01)
        ),
        formulation=str(
            activation_config.get(
                "comsol_formulation", "LocalEnergyApproximationE"
            )
        ),
    )
    activation_audit.insert(0, "evidence_status", evidence_status)
    townsend_representation = townsend_source_representation_audit(
        activation_export,
        intended_historical_bundle,
        gas_temperature_K=float(
            historical_comsol_conditions.get("gas_temperature_K", 293.15)
        ),
        pressure_Pa=float(
            historical_comsol_conditions.get(
                "pressure_Pa", external["gas_pressure"].iloc[0]
            )
        ),
    )
    townsend_representation.insert(0, "evidence_status", evidence_status)
    bundle_transport_path = (
        bundle / "transport_vs_en.csv"
        if (bundle / "transport_vs_en.csv").exists()
        else bundle / "transport_vs_mean_energy.csv"
    )
    bundle_rates_path = (
        bundle / "rates_vs_en.csv"
        if (bundle / "rates_vs_en.csv").exists()
        else bundle / "rates_vs_mean_energy.csv"
    )
    source_paths: list[tuple[str, Path]] = [
        ("config", config["_config_path"]),
        ("analysis_code", Path(__file__)),
        ("external_profile", external_source.path),
        ("reference_profile", reference_source.path),
        ("runtime", runtime_path),
        ("bundle_manifest", bundle / "manifest.json"),
        ("bundle_transport", bundle_transport_path),
        ("bundle_rates", bundle_rates_path),
        ("bundle_mean_energy", bundle / "mean_energy_vs_en.csv"),
        ("formal_swarm_database", formal_swarm_database),
        (
            "historical_activation_manifest",
            intended_historical_bundle / "manifest.json",
        ),
        (
            "historical_activation_mean_energy",
            intended_historical_bundle / "mean_energy_vs_en.csv",
        ),
        (
            "historical_activation_transport_mean_energy",
            intended_historical_bundle / "transport_vs_mean_energy.csv",
        ),
        (
            "historical_activation_rates_mean_energy",
            intended_historical_bundle / "rates_vs_mean_energy.csv",
        ),
        ("historical_activation_export", activation_export_path),
    ]
    condition_audit = config.get("condition_audit", {})
    if isinstance(condition_audit, Mapping):
        historical_swarm = condition_audit.get("historical_external_swarm", {})
        formal_swarm = condition_audit.get("formal_swarm", {})
        for source_id, entry, key in (
            (
                "historical_swarm_config",
                historical_swarm,
                "evidence_path",
            ),
            (
                "formal_executed_config",
                formal_swarm,
                "executed_config_path",
            ),
            (
                "portable_repro_config",
                formal_swarm,
                "portable_repro_config_path",
            ),
        ):
            if isinstance(entry, Mapping) and entry.get(key):
                source_paths.append((source_id, resolve_path(str(entry[key]))))
    if fallback_path:
        source_paths.append(("mesh_summary_fallback", fallback_path))
    for source, _ in mesh_sources:
        source_paths.append(("mesh_profiles", source.path))
    qa, source_hashes, qa_status = qa_sources(
        config,
        external_source,
        external,
        reference_source,
        reference,
        mesh_sources,
        runtime,
        runtime_summary,
        current_summary,
        activation_audit,
        tail_audit,
        fallback_used,
        source_paths,
    )
    output_tables = {
        "comparison_metrics": metrics,
        "aligned_profiles_long": aligned,
        "regional_error_attribution": regional,
        "current_uniformity": current_summary,
        "current_rsd_sensitivity": current_sensitivity,
        "runtime_summary": runtime_summary,
        "runtime_speedup": runtime_speedup,
        "clamp_metrics": clamp,
        "clamp_impact": clamp_impact,
        "regime_metrics": regime,
        "regime_impact": regime_impact,
        "fluid_applicability_diagnostics": fluid_diagnostics,
        "operating_point_diagnostics": operating_point,
        "townsend_source_representation_audit": townsend_representation,
        "partial_electron_energy_audit": partial_energy,
        "tail_convergence_audit": tail_audit,
        "historical_activation_audit": activation_audit,
        "qa_summary": qa,
        "source_hashes": source_hashes,
    }
    for name, frame in output_tables.items():
        frame.to_csv(output_data / f"{name}.csv", index=False)
    expected_repetitions = int(
        config.get("expected", {}).get("runtime_repetitions_per_path", 3)
    )
    figure_metadata = generate_figures(
        figure_dir,
        evidence_status,
        profile_evidence_status,
        lookup_evidence_status,
        bundle,
        external_source,
        external,
        reference_source,
        reference,
        e_over_n_unit,
        metrics,
        clamp,
        clamp_context,
        current_summary,
        regional,
        current_sensitivity,
        activation_audit,
        regime_impact,
        regime_context,
        fluid_diagnostics,
        runtime,
        runtime_summary,
        expected_repetitions,
        source_hashes,
    )
    write_json(output_data / "figure_metadata.json", figure_metadata)
    generated_at = datetime.now(timezone.utc).astimezone().isoformat()
    low_clamp_fraction = float(
        clamp.loc[
            clamp["metric"] == "E_over_N_below_table", "spatial_fraction"
        ].iloc[0]
    )
    formal_runtime = runtime_summary[
        runtime_summary["scope"].isin(["solve_only", "end_to_end"])
        & runtime_summary["path"].isin(["external_swarm", "builtin_reference"])
    ]
    summary = {
        "analysis_version": 2,
        "generated_at": generated_at,
        "evidence_status": evidence_status,
        "data_vintages": {
            "profile_comparison": profile_evidence_status,
            "runtime_and_mesh": profile_evidence_status,
            "lookup_tables": lookup_evidence_status,
            "lookup_range_overlay": (
                f"lookup curves: {lookup_evidence_status}; COMSOL-used field range: "
                f"{profile_evidence_status}"
            ),
            "historical_activation_audit": (
                "intended historical bundle versus archived export; stale-class "
                "provenance unresolved"
            ),
        },
        "config": relative_path(config["_config_path"]),
        "method_definitions": {
            "relative_L2_spatial_weighted": (
                "sqrt(integral((external-reference)^2 dx) / "
                "integral(reference^2 dx)); exact integral of the product of "
                "union-grid piecewise-linear reconstructions"
            ),
            "relative_L2_trapezoidal_of_squared_samples_sensitivity": (
                "legacy sensitivity calculation that applies the trapezoidal "
                "rule to squared nodal samples"
            ),
            "signed_integral_ratio_external_over_reference": (
                "signed integral(external dx) / signed integral(reference dx); "
                "piecewise-linear integrals on the same grid"
            ),
            "shape_correlation_spatial_weighted": (
                "Pearson correlation from exact means, variances, and covariance "
                "of union-grid piecewise-linear reconstructions"
            ),
            "current_RSD": (
                "exact piecewise-linear spatial population standard deviation "
                "divided by absolute piecewise-linear spatial mean"
            ),
            "central_80_percent": (
                "geometric window after excluding 10% of physical length at "
                "each boundary; it is not asserted to be a bulk/sheath boundary"
            ),
            "regional_error_attribution": (
                "exact squared-error share, region-normalized L2, reference-norm "
                "share, and signed integral ratio for cathode-side 10%, central "
                "80%, and anode-side 10% geometric regions"
            ),
            "high_field_regime": (
                "absolute-integral shares where E/N exceeds 500 Td; 500 Td is "
                "typical COMSOL drift-diffusion guidance, not a hard validity cutoff"
            ),
            "fluid_applicability_diagnostics": (
                "max(n_e/N) plus exact spatial support of "
                "|J_e|/[e n_e sqrt(2 e mean_energy/m_e)] above operational "
                "guides; the latter includes drift, diffusion, and boundary flux "
                "and is not a hard validity test"
            ),
            "operating_point_diagnostics": (
                "source voltage, endpoint gap voltage, ballast drop/current, "
                "internal potential maximum, and absolute axial potential "
                "variation derived from final-time profiles"
            ),
            "townsend_source_representation_audit": (
                "archived Townsend-flux source compared with N n_e k(mean_energy) "
                "as a postprocess counterfactual; not a matched COMSOL rerun"
            ),
            "partial_electron_energy_audit": (
                "signed integral J_e E dx and threshold-weighted eir2/eir4 "
                "losses only; all other volume, boundary, flux, and transient "
                "energy terms are omitted"
            ),
            "tail_convergence_audit": (
                "per-case workflow-database diagnostics; every formal point "
                "must report solver convergence, tail probability <=1e-9, "
                "edge/peak <=1e-10, and no 20 keV grid-limit hit"
            ),
            "historical_activation_audit": (
                "intended historical bundle interpolated at archived export "
                "arguments; this is not activation verification because stale-class "
                "provenance is unresolved and mean-energy(E/N) is inactive in LEA"
            ),
        },
        "comparison": {
            "external_profile": relative_path(external_source.path),
            "reference_profile": relative_path(reference_source.path),
            "quantities": len(metrics),
        },
        "condition_audit": config.get("condition_audit", {}),
        "historical_activation_audit": config.get(
            "historical_activation_audit", {}
        ),
        "clamp_context": clamp_context,
        "regime_context": regime_context,
        "low_E_over_N_clamp_spatial_fraction": low_clamp_fraction,
        "qa": qa_status,
        "formal_gates": {
            "runtime_repetition_gate_met": bool(
                len(formal_runtime) == 4
                and (formal_runtime["repetitions"] >= expected_repetitions).all()
            ),
            "formal_swarm_tail_convergence_gate_met": bool(
                len(tail_audit) > 0
                and tail_audit["formal_tail_grid_gate_pass"].all()
            ),
            "formal_swarm_solver_convergence_gate_met": bool(
                len(tail_audit) > 0
                and tail_audit[
                    "formal_solver_convergence_gate_pass"
                ].all()
            ),
            "formal_swarm_all_case_gates_met": bool(
                len(tail_audit) > 0
                and tail_audit["formal_swarm_case_gate_pass"].all()
            ),
            "raw_mesh_200_400_gate_met": bool(
                {
                    int(value)
                    for value in config.get("expected", {}).get("mesh_elements", [])
                }
                <= {source.mesh_elements for source, _ in mesh_sources}
                and not fallback_used
            ),
            "profile_voltage_pressure_match": bool(
                qa.loc[
                    qa["check"].isin(
                        ["comparison_voltage_match", "comparison_pressure_match"]
                    ),
                    "passed",
                ].all()
            ),
            "historical_closure_gas_conditions_match": bool(
                qa.loc[
                    qa["check"] == "historical_closure_gas_conditions_match",
                    "passed",
                ].all()
            ),
            "formal_target_conditions_configured_to_match": bool(
                config.get("condition_audit", {}).get(
                    "formal_target_conditions_configured_to_match", False
                )
            ),
            "formal_run_conditions_verified": bool(
                config.get("condition_audit", {}).get(
                    "formal_run_conditions_verified", False
                )
            ),
            "historical_activation_provenance_verified": bool(
                qa.loc[
                    qa["check"]
                    == "historical_activation_provenance_verified",
                    "passed",
                ].all()
                and not qa.loc[
                    qa["check"]
                    == "historical_activation_provenance_verified"
                ].empty
            ),
            "profile_numeric_QA_passed": bool(
                qa.loc[
                    qa["check"].isin(
                        [
                            "all_numeric_values_finite",
                            "electron_density_nonnegative",
                            "excitation_source_nonnegative",
                            "ionization_source_nonnegative",
                        ]
                    ),
                    "passed",
                ].all()
            ),
        },
        "derived_tables": {
            name: {
                "path": relative_path(output_data / f"{name}.csv"),
                "rows": len(frame),
                "sha256": sha256_file(output_data / f"{name}.csv"),
            }
            for name, frame in output_tables.items()
        },
        "figures": {
            item["figure_id"]: item["outputs"] for item in figure_metadata
        },
    }
    write_json(output_data / "analysis_summary.json", summary)
    return summary


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="Input JSON; paths inside it are relative to the repository root.",
    )
    args = parser.parse_args(argv)
    summary = run_analysis(args.config)
    print(
        json.dumps(
            {
                "status": summary["qa"]["data_quality_status"],
                "evidence_status": summary["evidence_status"],
                "derived_tables": len(summary["derived_tables"]),
                "figures": len(summary["figures"]),
                "low_E_over_N_clamp_spatial_fraction": summary[
                    "low_E_over_N_clamp_spatial_fraction"
                ],
                "formal_gates": summary["formal_gates"],
            },
            indent=2,
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
