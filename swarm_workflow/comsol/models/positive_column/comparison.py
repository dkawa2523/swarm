"""Unit-aware, spatially weighted comparison of saved positive-column results."""

from dataclasses import dataclass
import math
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from ...._io import write_csv, write_json
from .contracts import PositiveColumnWorkflowError
from .results import read_positive_column_results
from .result_units import QUANTITIES

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
    """The saved results do not define a compatible spatial comparison."""


@dataclass(frozen=True, slots=True)
class ComsolComparisonSummary:
    output_dir: Path
    metrics_csv: Path
    summary_json: Path
    quantities: tuple[str, ...]
    speedup: float | None
    profile_plot: Path | None
    profile_svg: Path | None
    profiles_csv: Path


def compare_comsol_profiles(
    external_results: str | Path,
    reference_results: str | Path,
    *,
    output_dir: str | Path,
    external_solution: int | None = None,
    reference_solution: int | None = None,
    external_runtime_s: float | None = None,
    reference_runtime_s: float | None = None,
    create_plot: bool = False,
) -> ComsolComparisonSummary:
    """Compare one selected solution on the complete common spatial interval.

    Inputs are result summaries with CSV hashes and units. The output directory
    must be new; the complete comparison is published by one directory rename.
    """
    output = Path(output_dir).resolve()
    if output.exists():
        raise ComsolCompareError(
            "comparison output already exists; choose a new output directory"
        )
    try:
        artifacts = [
            read_positive_column_results(path)
            for path in (external_results, reference_results)
        ]
        external, external_id = _profile(artifacts[0], external_solution)
        reference, reference_id = _profile(artifacts[1], reference_solution)
        _matching_physics(artifacts)
    except (PositiveColumnWorkflowError, OSError, ValueError, KeyError) as error:
        raise ComsolCompareError(str(error)) from error
    conditions = _matching_conditions(external, reference)
    quantities = tuple(
        name for name in PROFILE_COLUMNS if name in external and name in reference
    )
    if not quantities:
        raise ComsolCompareError("profiles have no common supported result columns")
    lower = max(external["x"][0], reference["x"][0])
    upper = min(external["x"][-1], reference["x"][-1])
    if lower >= upper:
        raise ComsolCompareError("profiles do not share a usable spatial interval")
    x = np.unique(
        np.r_[
            lower,
            upper,
            external["x"][(external["x"] > lower) & (external["x"] < upper)],
            reference["x"][(reference["x"] > lower) & (reference["x"] < upper)],
        ]
    )
    metric_rows, profile_rows, curves = [], [], {}
    for name in quantities:
        candidate = np.interp(x, external["x"], external[name])
        expected = np.interp(x, reference["x"], reference[name])
        difference = candidate - expected
        if not np.all(np.isfinite(difference)):
            raise ComsolCompareError(
                "profile differences exceed the floating-point range"
            )
        l2, ratio = _relative_metrics(x, candidate, expected)
        metric_rows.append(
            {
                "quantity": name,
                "unit": QUANTITIES[name].unit,
                "l2_relative_error": l2,
                "shape_correlation": _correlation(x, candidate, expected),
                "integral_ratio_external_over_reference": ratio,
                "external_min": float(candidate.min()),
                "external_max": float(candidate.max()),
                "reference_min": float(expected.min()),
                "reference_max": float(expected.max()),
            }
        )
        profile_rows.extend(
            dict(
                quantity=name,
                x_m=float(position),
                external=float(a),
                reference=float(b),
                difference=float(d),
                unit=QUANTITIES[name].unit,
            )
            for position, a, b, d in zip(
                x, candidate, expected, difference, strict=True
            )
        )
        curves[name] = (candidate, expected)
    runtimes = [external_runtime_s, reference_runtime_s]
    if (runtimes[0] is None) != (runtimes[1] is None):
        raise ComsolCompareError("supply both solve runtimes or neither")
    scope = "user_supplied_solve_times"
    if runtimes[0] is None:
        runtimes = [
            artifact.metadata.get("timing_s", {}).get("run") for artifact in artifacts
        ]
        scope = "saved_run_stage"
    speedup = _speedup(*runtimes)
    names = (
        "comparison_metrics.csv",
        "comparison_profiles.csv",
        "comparison_summary.json",
    )
    plot_names = (
        ("spatial_profile_comparison.png", "spatial_profile_comparison.svg")
        if create_plot
        else ()
    )
    metadata = {
        "artifact_schema": "swarm.positive_column_comparison.v1",
        "evidence_scope": "saved_result_comparison; no new solve or physical acceptance",
        "sources": {
            key: {
                "summary": artifact.metadata["summary_path"],
                "summary_sha256": artifact.metadata["summary_sha256"],
                "csv_sha256": artifact.metadata["result_csv_sha256"],
                "selection": selection,
                "source_units": artifact.metadata["source_units"],
            }
            for key, artifact, selection in zip(
                ("external", "reference"),
                artifacts,
                (external_id, reference_id),
                strict=True,
            )
        },
        "physical_context": artifacts[0].metadata["physical_context"],
        "common_x_range_m": [float(lower), float(upper)],
        "operating_conditions": conditions,
        "quantities": list(quantities),
        "units": {name: QUANTITIES[name].unit for name in quantities},
        "interpolation": "piecewise linear on the union of native nodes and overlap endpoints; no extrapolation",
        "l2_measure": "sqrt(integral((external-reference)^2 dx) / integral(reference^2 dx))",
        "runtime": {
            "external_s": runtimes[0],
            "reference_s": runtimes[1],
            "scope": scope,
            "speedup_reference_over_external": speedup,
        },
        "conductive_current_uniformity": {
            "external": _current_uniformity(external),
            "reference": _current_uniformity(reference),
        }
        if all("total_current_density" in profile for profile in (external, reference))
        else None,
        "formulations": {
            key: artifact.metadata.get("mean_energy_formulation", {}).get("mode")
            for key, artifact in zip(("external", "reference"), artifacts, strict=True)
        },
        "files": list(names[:2] + plot_names),
    }
    # Errors leave prior comparisons and source artifacts untouched.
    output.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(
        prefix=".positive-column-comparison-", dir=output.parent
    ) as temporary:
        staging = Path(temporary).resolve()
        if staging.parent != output.parent:
            raise ComsolCompareError("comparison staging escaped its destination")
        write_csv(staging / names[0], tuple(metric_rows[0]), metric_rows)
        write_csv(
            staging / names[1],
            ("quantity", "x_m", "external", "reference", "difference", "unit"),
            profile_rows,
        )
        if create_plot:
            _plot_profiles(staging, x, curves, metric_rows, conditions)
        write_json(staging / names[2], metadata)
        staging.rename(output)
    return ComsolComparisonSummary(
        output,
        output / names[0],
        output / names[2],
        quantities,
        speedup,
        output / plot_names[0] if create_plot else None,
        output / plot_names[1] if create_plot else None,
        output / names[1],
    )


def _profile(artifact, selection):
    available = sorted({row["solution_index"] for row in artifact.rows})
    if selection is None:
        if len(available) != 1:
            raise ComsolCompareError(
                f"profile comparison requires one saved solution; choose a solution index from {available}"
            )
        selection = available[0]
    if (
        isinstance(selection, bool)
        or not isinstance(selection, int)
        or selection not in available
    ):
        raise ComsolCompareError(
            f"saved solution {selection!r} is unavailable; choose from {available}"
        )
    rows = [row for row in artifact.rows if row["solution_index"] == selection]
    identity = {
        name: rows[0][name]
        for name in (
            "dataset",
            "solution_index",
            "time_s",
            "parameter_name",
            "parameter_unit",
            "parameter_value",
        )
    }
    if identity["parameter_name"] not in {
        None,
        "t",
        artifact.metadata.get("voltage_parameter"),
    }:
        raise ComsolCompareError(
            "profile comparison supports saved time or the configured voltage parameter; other parameters are unsupported"
        )
    if len(rows) < 2 or "x" not in artifact.units:
        raise ComsolCompareError(
            "spatial comparison requires x and at least two points"
        )
    columns = {
        name: np.array([row[name] for row in rows], dtype=float)
        for name in ("x", *PROFILE_COLUMNS, *CONDITION_COLUMNS)
        if name in rows[0]
    }
    if identity["time_s"] is not None:
        columns["time_s"] = np.array([row["time_s"] for row in rows])
    order = np.argsort(columns["x"])
    if np.any(np.diff(columns["x"][order]) <= 0):
        raise ComsolCompareError("spatial x values must be unique")
    return {name: values[order] for name, values in columns.items()}, identity


def _matching_physics(artifacts):
    contexts = []
    for artifact in artifacts:
        context = artifact.metadata.get("physical_context")
        if (
            not isinstance(context, dict)
            or "gas_temperature_K" not in context
            or "gas_mixture" not in context
        ):
            raise ComsolCompareError("result is missing gas physical_context")
        contexts.append(
            {
                **context,
                "gas_mixture": sorted(
                    context["gas_mixture"], key=lambda item: item["species"]
                ),
            }
        )
    if contexts[0] != contexts[1]:
        raise ComsolCompareError("gas physical_context differs between results")


def _matching_conditions(external, reference):
    conditions = {}
    names = CONDITION_COLUMNS + (
        ("time_s",) if "time_s" in external or "time_s" in reference else ()
    )
    for name in names:
        if name not in external or name not in reference:
            raise ComsolCompareError(
                f"both profiles must contain operating-condition column {name!r}"
            )
        for profile in (external, reference):
            if not np.allclose(profile[name], profile[name][0], rtol=1e-9, atol=1e-12):
                raise ComsolCompareError(
                    f"operating condition {name!r} must be constant within a profile"
                )
        a, b = float(external[name][0]), float(reference[name][0])
        if not math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12):
            raise ComsolCompareError(
                f"operating condition {name!r} differs: external={a:g}, reference={b:g}"
            )
        conditions[name] = a
    return conditions


def _square_integral(x, y):
    return float(np.sum(np.diff(x) * (y[:-1] ** 2 + y[:-1] * y[1:] + y[1:] ** 2) / 3))


def _relative_metrics(x, a, b):
    # Normalize units/magnitudes before squaring: the same length factor cancels
    # in both ratios. Native curve values stay unchanged in the exported CSV.
    x = (x - x[0]) / (x[-1] - x[0])
    scale = float(np.max(np.abs(b)))
    if scale == 0:
        return None, None
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            a, b = a / scale, b / scale
            delta = a - b
            magnitude = float(np.max(np.abs(delta)))
            l2 = (
                magnitude
                * math.sqrt(
                    _square_integral(x, delta / magnitude) / _square_integral(x, b)
                )
                if magnitude
                else 0.0
            )
            integral = float(np.trapezoid(b, x))
            ratio = float(np.trapezoid(a, x)) / integral if integral != 0 else None
    except (FloatingPointError, OverflowError) as error:
        raise ComsolCompareError(
            "profile metrics exceed the floating-point range"
        ) from error
    if not math.isfinite(l2) or (ratio is not None and not math.isfinite(ratio)):
        raise ComsolCompareError("profile metrics exceed the floating-point range")
    return l2, ratio


def _correlation(x, a, b):
    if np.all(a == a[0]) or np.all(b == b[0]):
        return None
    x = (x - x[0]) / (x[-1] - x[0])
    a, b = a / np.max(np.abs(a)), b / np.max(np.abs(b))
    a = a - np.trapezoid(a, x)
    b = b - np.trapezoid(b, x)
    aa, bb = _square_integral(x, a), _square_integral(x, b)
    if aa == 0 or bb == 0:
        return None
    ab = float(
        np.sum(
            np.diff(x)
            * (
                2 * a[:-1] * b[:-1]
                + a[:-1] * b[1:]
                + a[1:] * b[:-1]
                + 2 * a[1:] * b[1:]
            )
            / 6
        )
    )
    return float(np.clip(ab / math.sqrt(aa * bb), -1.0, 1.0))


def _speedup(external, reference):
    for value in (external, reference):
        if value is not None and (
            isinstance(value, bool) or not math.isfinite(value) or value <= 0
        ):
            raise ComsolCompareError("solve runtimes must be positive and finite")
    if external is None or reference is None:
        return None
    speedup = reference / external
    if not math.isfinite(speedup):
        raise ComsolCompareError("runtime ratio exceeds the floating-point range")
    return speedup


def _current_uniformity(profile):
    native = profile["x"]
    current = profile["total_current_density"]

    def stats(lower, upper):
        x = np.unique(np.r_[lower, native[(native > lower) & (native < upper)], upper])
        values = np.interp(x, native, current)
        fraction = (x - lower) / (upper - lower)
        scale = float(np.max(np.abs(values)))
        normalized = values / scale if scale else values
        normalized_mean = float(np.trapezoid(normalized, fraction))
        mean = normalized_mean * scale
        std = (
            math.sqrt(_square_integral(fraction, normalized - normalized_mean)) * scale
        )
        return {
            "mean_A_per_m2": mean,
            "standard_deviation_A_per_m2": std,
            "relative_standard_deviation": std / abs(mean) if mean != 0 else None,
            "min_A_per_m2": float(values.min()),
            "max_A_per_m2": float(values.max()),
        }

    lower, upper = float(native[0]), float(native[-1])
    a, b = lower + 0.1 * (upper - lower), upper - 0.1 * (upper - lower)
    return {
        "measure": "piecewise-linear spatial integrals",
        "full_domain": stats(lower, upper),
        "central_80_percent": {"x_range_m": [a, b], **stats(a, b)},
    }


def _plot_profiles(directory, x, curves, metrics, conditions):
    try:
        from matplotlib.figure import Figure
    except ImportError as error:
        raise ComsolCompareError(
            "profile plotting requires electron-swarm[plot]"
        ) from error
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib import rc_context

    with rc_context({"svg.fonttype": "none", "text.usetex": False}):
        fig = Figure(figsize=(12, 3.4 * math.ceil(len(curves) / 2) + 0.7))
        FigureCanvasAgg(fig)
        axes = np.atleast_1d(fig.subplots(math.ceil(len(curves) / 2), 2)).ravel()
        for axis, row in zip(axes, metrics):
            quantity = QUANTITIES[row["quantity"]]
            a, b = curves[row["quantity"]]
            axis.plot(x, a * quantity.display_scale, marker=".", label="External")
            axis.plot(
                x,
                b * quantity.display_scale,
                linestyle="--",
                marker=".",
                label="Reference",
            )
            if quantity.logarithmic and np.all(a > 0) and np.all(b > 0):
                axis.set_yscale("log")
            axis.set_title(quantity.label)
            axis.set_ylabel(quantity.display_unit or quantity.unit)
            axis.set_xlabel("Axial position [m]")
            axis.grid(True, alpha=0.2)
            axis.legend()
            value = row["l2_relative_error"]
            axis.text(
                0.02,
                0.02,
                f"Spatial L2: {value:.4g}"
                if value is not None
                else "Spatial L2: undefined (zero reference)",
                transform=axis.transAxes,
                fontsize=8,
            )
        for axis in axes[len(curves) :]:
            axis.set_visible(False)
        fig.suptitle(
            f"Positive column | {conditions['applied_voltage']:g} V; {conditions['gas_pressure']:g} Pa"
            + (f"; t={conditions['time_s']:g} s" if "time_s" in conditions else "")
        )
        fig.tight_layout(rect=(0, 0.025, 1, 0.95))
        fig.text(
            0.02,
            0.008,
            "Saved solutions; piecewise-linear comparison on common support. No independent physical validation is inferred.",
            fontsize=8,
        )
        for suffix in ("png", "svg"):
            fig.savefig(directory / f"spatial_profile_comparison.{suffix}", dpi=160)
        fig.clear()


def format_comsol_comparison(summary):
    return "\n".join(
        [
            "Positive-column comparison written",
            f"quantities: {len(summary.quantities)}",
            "speedup: "
            + (
                f"{summary.speedup:.3g}x"
                if summary.speedup is not None
                else "not recorded"
            ),
            f"metrics: {summary.metrics_csv}",
            f"profiles: {summary.profiles_csv}",
            f"summary: {summary.summary_json}",
            *(
                [f"plot: {summary.profile_plot}", f"svg: {summary.profile_svg}"]
                if summary.profile_plot
                else []
            ),
        ]
    )
