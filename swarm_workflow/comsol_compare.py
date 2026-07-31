"""Small, reproducible comparison of two exported COMSOL profiles."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import math
from pathlib import Path
from typing import Any

import numpy as np

from ._io import write_csv, write_json


PROFILE_COLUMNS = (
    "electron_density",
    "mean_electron_energy",
    "electric_potential",
    "E_over_N",
    "electron_current_density",
    "total_current_density",
    "excitation_source",
    "ionization_source",
)
CONDITION_COLUMNS = ("applied_voltage", "gas_pressure")


class ComsolCompareError(ValueError):
    """Raised when exported COMSOL profiles cannot be compared."""


@dataclass(frozen=True, slots=True)
class ComsolComparisonSummary:
    output_dir: Path
    metrics_csv: Path
    summary_json: Path
    quantities: tuple[str, ...]
    speedup: float | None
    profile_plot: Path | None


def compare_comsol_profiles(
    external_csv: str | Path,
    reference_csv: str | Path,
    *,
    output_dir: str | Path,
    external_runtime_s: float | None = None,
    reference_runtime_s: float | None = None,
    create_plot: bool = False,
) -> ComsolComparisonSummary:
    """Compare canonical profiles on their common spatial interval."""
    external_path = Path(external_csv).resolve()
    reference_path = Path(reference_csv).resolve()
    output = Path(output_dir).resolve()
    external = _read_profile(external_path)
    reference = _read_profile(reference_path)
    conditions = _matching_conditions(external, reference)
    quantities = tuple(
        name for name in PROFILE_COLUMNS if name in external and name in reference
    )
    if not quantities:
        raise ComsolCompareError("profiles have no common supported result columns")

    x_reference = reference["x"]
    lower = max(float(external["x"][0]), float(x_reference[0]))
    upper = min(float(external["x"][-1]), float(x_reference[-1]))
    mask = (x_reference >= lower) & (x_reference <= upper)
    if int(np.count_nonzero(mask)) < 2:
        raise ComsolCompareError("profiles do not share a usable spatial interval")
    x = x_reference[mask]

    metric_rows: list[dict[str, Any]] = []
    for quantity in quantities:
        candidate = np.interp(x, external["x"], external[quantity])
        expected = reference[quantity][mask]
        difference = candidate - expected
        expected_norm = float(np.linalg.norm(expected))
        l2_error = (
            float(np.linalg.norm(difference) / expected_norm)
            if expected_norm > 0.0
            else None
        )
        expected_integral = float(np.trapezoid(expected, x))
        candidate_integral = float(np.trapezoid(candidate, x))
        integral_ratio = (
            candidate_integral / expected_integral
            if abs(expected_integral) > 0.0
            else None
        )
        correlation = _correlation(candidate, expected)
        metric_rows.append(
            {
                "quantity": quantity,
                "l2_relative_error": l2_error,
                "shape_correlation": correlation,
                "integral_ratio_external_over_reference": integral_ratio,
                "external_min": float(np.min(candidate)),
                "external_max": float(np.max(candidate)),
                "reference_min": float(np.min(expected)),
                "reference_max": float(np.max(expected)),
            }
        )

    speedup = _speedup(external_runtime_s, reference_runtime_s)
    current_conservation = None
    if "total_current_density" in external and "total_current_density" in reference:
        current_conservation = {
            "external": _current_conservation(external),
            "reference": _current_conservation(reference),
        }
    metrics_csv = output / "comparison_metrics.csv"
    summary_json = output / "comparison_summary.json"
    profile_plot = output / "spatial_profile_comparison.png" if create_plot else None
    write_csv(
        metrics_csv,
        (
            "quantity",
            "l2_relative_error",
            "shape_correlation",
            "integral_ratio_external_over_reference",
            "external_min",
            "external_max",
            "reference_min",
            "reference_max",
        ),
        metric_rows,
    )
    if profile_plot is not None:
        _plot_profiles(
            profile_plot,
            x,
            external,
            reference,
            mask,
            quantities,
            metric_rows,
            conditions,
        )
    write_json(
        summary_json,
        {
            "external_csv": str(external_path),
            "reference_csv": str(reference_path),
            "common_x_range_m": [float(x[0]), float(x[-1])],
            "operating_conditions": conditions,
            "quantities": list(quantities),
            "runtime": {
                "external_s": external_runtime_s,
                "reference_s": reference_runtime_s,
                "speedup_reference_over_external": speedup,
            },
            "current_conservation": current_conservation,
            "profile_plot": str(profile_plot) if profile_plot is not None else None,
        },
    )
    return ComsolComparisonSummary(
        output_dir=output,
        metrics_csv=metrics_csv,
        summary_json=summary_json,
        quantities=quantities,
        speedup=speedup,
        profile_plot=profile_plot,
    )


def _matching_conditions(
    external: dict[str, np.ndarray],
    reference: dict[str, np.ndarray],
) -> dict[str, float]:
    conditions: dict[str, float] = {}
    for name in CONDITION_COLUMNS:
        if name not in external or name not in reference:
            raise ComsolCompareError(
                f"both profiles must contain operating-condition column {name!r}"
            )
        external_value = float(np.mean(external[name]))
        reference_value = float(np.mean(reference[name]))
        if not math.isclose(
            external_value,
            reference_value,
            rel_tol=1.0e-9,
            abs_tol=1.0e-12,
        ):
            raise ComsolCompareError(
                f"operating condition {name!r} differs: "
                f"external={external_value:.12g}, reference={reference_value:.12g}"
            )
        conditions[name] = external_value
    return conditions


def format_comsol_comparison(summary: ComsolComparisonSummary) -> str:
    speedup = "not recorded" if summary.speedup is None else f"{summary.speedup:.3g}x"
    lines = [
        "COMSOL comparison complete",
        f"quantities: {len(summary.quantities)}",
        f"speedup: {speedup}",
        f"metrics: {summary.metrics_csv}",
        f"summary: {summary.summary_json}",
    ]
    if summary.profile_plot is not None:
        lines.append(f"plot: {summary.profile_plot}")
    return "\n".join(lines)


def _plot_profiles(
    path: Path,
    x: np.ndarray,
    external: dict[str, np.ndarray],
    reference: dict[str, np.ndarray],
    reference_mask: np.ndarray,
    quantities: tuple[str, ...],
    metric_rows: list[dict[str, Any]],
    conditions: dict[str, float],
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import ScalarFormatter
    except ImportError as exc:
        raise ComsolCompareError(
            "profile plotting requires the 'plot' extra: pip install -e .[plot]"
        ) from exc

    labels = {
        "electron_density": ("Electron density", r"m$^{-3}$", 1.0, True),
        "mean_electron_energy": ("Mean electron energy", "eV", 1.0, False),
        "electric_potential": ("Electric potential", "V", 1.0, False),
        "E_over_N": ("Reduced electric field", "Td", 1.0e21, True),
        "electron_current_density": (
            "Electron current density",
            r"A m$^{-2}$",
            1.0,
            False,
        ),
        "total_current_density": (
            "Total current density",
            r"A m$^{-2}$",
            1.0,
            False,
        ),
        "excitation_source": ("Excitation source", r"m$^{-3}$ s$^{-1}$", 1.0, True),
        "ionization_source": ("Ionization source", r"m$^{-3}$ s$^{-1}$", 1.0, True),
    }
    ncols = 2
    nrows = math.ceil(len(quantities) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(13, 3.6 * nrows), sharex=True)
    axes_array = np.atleast_1d(axes).ravel()
    x_cm = x * 100.0
    metrics = {row["quantity"]: row for row in metric_rows}

    for axis, quantity in zip(axes_array, quantities, strict=False):
        title, unit, scale, use_log = labels[quantity]
        candidate = np.interp(x, external["x"], external[quantity]) * scale
        expected = reference[quantity][reference_mask] * scale
        axis.plot(x_cm, candidate, color="#0072B2", linewidth=1.8, label="External Swarm")
        axis.plot(
            x_cm,
            expected,
            color="#D55E00",
            linewidth=1.6,
            linestyle="--",
            label="COMSOL built-in",
        )
        if use_log and np.all(candidate > 0.0) and np.all(expected > 0.0):
            axis.set_yscale("log")
        elif not use_log:
            formatter = ScalarFormatter(useMathText=True)
            formatter.set_powerlimits((-3, 4))
            axis.yaxis.set_major_formatter(formatter)
        row = metrics[quantity]
        correlation = row["shape_correlation"]
        ratio = row["integral_ratio_external_over_reference"]
        correlation_text = "n/a" if correlation is None else f"{correlation:.3f}"
        ratio_text = "n/a" if ratio is None else f"{ratio:.3f}"
        axis.text(
            0.02,
            0.04,
            f"corr={correlation_text}   L2={row['l2_relative_error']:.3f}   integral={ratio_text}",
            transform=axis.transAxes,
            fontsize=8.5,
            bbox={"facecolor": "white", "alpha": 0.82, "edgecolor": "0.8"},
        )
        axis.set_title(title)
        axis.set_ylabel(unit)
        axis.grid(True, alpha=0.25)
        if quantity == "total_current_density":
            lower = x_cm[0] + 0.1 * (x_cm[-1] - x_cm[0])
            upper = x_cm[-1] - 0.1 * (x_cm[-1] - x_cm[0])
            bulk = (x_cm >= lower) & (x_cm <= upper)
            inset = axis.inset_axes((0.54, 0.54, 0.43, 0.39))
            inset.plot(x_cm[bulk], candidate[bulk], color="#0072B2", linewidth=1.2)
            inset.plot(
                x_cm[bulk],
                expected[bulk],
                color="#D55E00",
                linewidth=1.1,
                linestyle="--",
            )
            inset.set_title("Central 80%", fontsize=8)
            inset.tick_params(labelsize=7)
            inset.grid(True, alpha=0.2)

    for axis in axes_array[len(quantities) :]:
        axis.set_visible(False)
    for axis in axes_array[-ncols:]:
        if axis.get_visible():
            axis.set_xlabel("Axial position (cm)")
    axes_array[0].legend(loc="best", frameon=True)
    fig.suptitle(
        "COMSOL spatial profiles at "
        f"{conditions['applied_voltage']:.6g} V and "
        f"{conditions['gas_pressure']:.6g} Pa",
        fontsize=14,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.98))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def _read_profile(path: Path) -> dict[str, np.ndarray]:
    if not path.is_file():
        raise ComsolCompareError(f"result CSV does not exist: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as fp:
        rows = list(csv.DictReader(fp))
    if len(rows) < 2 or "x" not in (rows[0] if rows else {}):
        raise ComsolCompareError(f"result CSV requires x and at least two rows: {path}")
    columns: dict[str, np.ndarray] = {}
    for name in ("x", *PROFILE_COLUMNS, *CONDITION_COLUMNS):
        if name not in rows[0]:
            continue
        try:
            values = np.asarray([float(row[name]) for row in rows], dtype=float)
        except (TypeError, ValueError, KeyError) as exc:
            raise ComsolCompareError(f"{path}: column {name} is not numeric") from exc
        if not np.all(np.isfinite(values)):
            raise ComsolCompareError(f"{path}: column {name} contains non-finite values")
        columns[name] = values
    order = np.argsort(columns["x"])
    x = columns["x"][order]
    if np.any(np.diff(x) <= 0.0):
        raise ComsolCompareError(f"{path}: x values must be unique")
    return {name: values[order] for name, values in columns.items()}


def _correlation(candidate: np.ndarray, reference: np.ndarray) -> float | None:
    if float(np.std(candidate)) == 0.0 or float(np.std(reference)) == 0.0:
        return None
    value = float(np.corrcoef(candidate, reference)[0, 1])
    return value if math.isfinite(value) else None


def _speedup(external_s: float | None, reference_s: float | None) -> float | None:
    for value, name in ((external_s, "external"), (reference_s, "reference")):
        if value is not None and (not math.isfinite(value) or value <= 0.0):
            raise ComsolCompareError(f"{name} runtime must be positive and finite")
    if external_s is None or reference_s is None:
        return None
    return reference_s / external_s


def _current_conservation(profile: dict[str, np.ndarray]) -> dict[str, Any]:
    """Summarize total current globally and away from the two sheath regions."""
    x = profile["x"]
    current = profile["total_current_density"]
    span = float(x[-1] - x[0])
    lower = float(x[0] + 0.1 * span)
    upper = float(x[-1] - 0.1 * span)
    bulk = current[(x >= lower) & (x <= upper)]

    def stats(values: np.ndarray) -> dict[str, float | None]:
        mean = float(np.mean(values))
        standard_deviation = float(np.std(values))
        return {
            "mean_A_per_m2": mean,
            "standard_deviation_A_per_m2": standard_deviation,
            "relative_standard_deviation": (
                standard_deviation / abs(mean) if mean != 0.0 else None
            ),
            "min_A_per_m2": float(np.min(values)),
            "max_A_per_m2": float(np.max(values)),
        }

    return {
        "full_domain": stats(current),
        "bulk_central_80_percent": {
            "x_range_m": [lower, upper],
            **stats(bulk),
        },
    }
