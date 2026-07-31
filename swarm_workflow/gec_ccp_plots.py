"""Reproducible static plots for the argon GEC CCP workflow."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any


# Chart contract:
# - Coefficient plots answer which DC-swarm closure COMSOL receives.
# - Spatial plots compare built-in Druyvesteyn and Swarm-table closures.
# - Waveform plots compare electrode voltage/current over an RF period.
# - COMSOL plots are skipped, never fabricated, when result CSVs are absent.


class GecCcpPlotError(RuntimeError):
    """Raised when GEC CCP plots cannot be produced from supplied data."""


@dataclass(frozen=True, slots=True)
class GecCcpPlotSummary:
    output_directory: Path
    figures: tuple[Path, ...]
    manifest: Path
    skipped: tuple[str, ...]


def plot_gec_ccp_results(
    bundle_path: str | Path,
    *,
    output_dir: str | Path,
    comsol_results_dir: str | Path | None = None,
) -> GecCcpPlotSummary:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError as exc:
        raise GecCcpPlotError(
            "plot-gec-ccp requires the optional matplotlib dependency"
        ) from exc

    bundle = Path(bundle_path).resolve()
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    figures: list[Path] = []
    skipped: list[str] = []

    transport = _numeric_rows(bundle / "transport_vs_mean_energy.csv")
    rates = _rows(bundle / "rates_vs_mean_energy.csv")
    eedf = _numeric_rows(bundle / "eedf.csv")
    if not transport or not rates or not eedf:
        raise GecCcpPlotError("bundle lacks transport, rate, or EEDF data")

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    mean_energy = np.asarray([row["mean_energy_eV"] for row in transport])
    axes[0, 0].plot(
        mean_energy,
        [row["reduced_mobility_m2_V_s_m3"] for row in transport],
        color="#1769aa",
        linewidth=2,
    )
    axes[0, 0].set_ylabel(r"$\mu_e N$ [m$^{-1}$ V$^{-1}$ s$^{-1}$]")
    axes[0, 0].set_title("Reduced electron mobility")

    axes[0, 1].plot(
        mean_energy,
        [row["reduced_diffusion_L_m2_s_m3"] for row in transport],
        color="#1769aa",
        linewidth=2,
        label="particle",
    )
    axes[0, 1].plot(
        mean_energy,
        [
            row["reduced_electron_energy_diffusion_m2_s_m3"]
            for row in transport
        ],
        color="#c04b35",
        linewidth=2,
        linestyle="--",
        label="energy",
    )
    axes[0, 1].set_ylabel(r"Reduced diffusion [m$^{-1}$ s$^{-1}$]")
    axes[0, 1].set_title("Reduced diffusion closure")
    axes[0, 1].legend(frameon=False)

    rate_colors = {
        "elastic": "#555555",
        "excitation": "#e69500",
        "ionization": "#8f3985",
    }
    rate_styles = {"elastic": "-", "excitation": "--", "ionization": "-."}
    for process in ("elastic", "excitation", "ionization"):
        selected = [
            row
            for row in rates
            if row.get("process_type") == process
            and _positive(row.get("rate_coefficient_m3_s"))
        ]
        axes[1, 0].plot(
            [float(row["mean_energy_eV"]) for row in selected],
            [float(row["rate_coefficient_m3_s"]) for row in selected],
            color=rate_colors[process],
            linestyle=rate_styles[process],
            linewidth=2,
            label=process,
        )
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_ylabel(r"Rate coefficient [m$^3$ s$^{-1}$]")
    axes[1, 0].set_title("Imported reaction rates")
    axes[1, 0].legend(frameon=False)

    available_en = sorted({row["E_over_N_Td"] for row in eedf})
    requested = (10.0, 100.0, 1000.0, 10000.0)
    selected_en = sorted({_nearest(available_en, value) for value in requested})
    palette = ("#555555", "#1769aa", "#e69500", "#8f3985")
    for color, en_value in zip(palette, selected_en):
        selected = [
            row
            for row in eedf
            if math.isclose(row["E_over_N_Td"], en_value, rel_tol=1e-10)
            and row["eedf"] > 0.0
        ]
        axes[1, 1].plot(
            [row["electron_energy_eV"] for row in selected],
            [row["eedf"] for row in selected],
            color=color,
            linewidth=1.8,
            label=f"{en_value:g} Td",
        )
    axes[1, 1].set_xscale("log")
    axes[1, 1].set_yscale("log")
    axes[1, 1].set_ylabel(r"EEDF [eV$^{-1}$]")
    axes[1, 1].set_title("Two-term EEDF")
    axes[1, 1].legend(frameon=False)

    for axis in axes.flat:
        axis.set_xlabel("Mean energy [eV]" if axis is not axes[1, 1] else "Energy [eV]")
        axis.grid(True, which="both", alpha=0.2)
    for axis in (axes[0, 0], axes[0, 1], axes[1, 0]):
        axis.set_xscale("log")
    fig.suptitle(
        "Argon GEC CCP input closure\n"
        "Steady DC two-term Swarm tables supplied to the RF mean-energy model",
        fontsize=14,
    )
    figures.extend(_save(fig, output / "swarm_closure_overview", plt))

    if comsol_results_dir is None:
        skipped.append("COMSOL comparison: no result directory supplied")
    else:
        results = Path(comsol_results_dir).resolve()
        baseline = results / "builtin_druyvesteyn"
        external = results / "swarm_tables"
        axis_files = (
            baseline / "axis_period_average.csv",
            external / "axis_period_average.csv",
        )
        if all(path.exists() for path in axis_files):
            figures.extend(
                _plot_spatial_comparison(
                    axis_files[0],
                    axis_files[1],
                    output / "axis_period_average_comparison",
                    plt,
                )
            )
        else:
            skipped.append("axis COMSOL comparison: result CSVs are absent")
        waveform_files = (
            baseline / "electrode_waveform.csv",
            external / "electrode_waveform.csv",
        )
        if all(path.exists() for path in waveform_files):
            figures.extend(
                _plot_waveform_comparison(
                    waveform_files[0],
                    waveform_files[1],
                    output / "electrode_waveform_comparison",
                    plt,
                )
            )
        else:
            skipped.append("waveform COMSOL comparison: result CSVs are absent")

    manifest = output / "plot_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "stage": "plot-gec-ccp",
                "status": "ok",
                "bundle": str(bundle),
                "figures": [str(path) for path in figures],
                "skipped": skipped,
                "comparison_labels": {
                    "baseline": "COMSOL built-in Druyvesteyn",
                    "external": "Swarm mean-energy lookup tables",
                },
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return GecCcpPlotSummary(output, tuple(figures), manifest, tuple(skipped))


def _plot_spatial_comparison(
    baseline_path: Path,
    external_path: Path,
    output_stem: Path,
    plt: Any,
) -> list[Path]:
    baseline = _comsol_numeric_matrix(baseline_path)
    external = _comsol_numeric_matrix(external_path)
    n_fields = min(baseline.shape[1], external.shape[1]) - 1
    if n_fields < 1:
        raise GecCcpPlotError("COMSOL spatial CSV has no result fields")
    labels = (
        "Electron density",
        "Mean electron energy",
        "Potential",
        "Ionization source",
        "Absorbed power",
    )
    units = (r"m$^{-3}$", "eV", "V", r"m$^{-3}$ s$^{-1}$", r"W m$^{-3}$")
    count = min(n_fields, len(labels))
    fig, axes = plt.subplots(count, 1, figsize=(8, 2.5 * count), constrained_layout=True)
    if count == 1:
        axes = [axes]
    for index, axis in enumerate(axes):
        column = index + 1
        axis.plot(
            baseline[:, 0] * 100.0,
            baseline[:, column],
            color="#555555",
            linestyle="--",
            linewidth=2,
            label="built-in Druyvesteyn",
        )
        axis.plot(
            external[:, 0] * 100.0,
            external[:, column],
            color="#1769aa",
            linewidth=2,
            label="Swarm tables",
        )
        axis.set_ylabel(f"{labels[index]}\n[{units[index]}]")
        axis.grid(True, alpha=0.2)
    axes[-1].set_xlabel("Axis position [cm]")
    axes[0].legend(frameon=False)
    fig.suptitle("Argon GEC CCP period-averaged axial profiles")
    return _save(fig, output_stem, plt)


def _plot_waveform_comparison(
    baseline_path: Path,
    external_path: Path,
    output_stem: Path,
    plt: Any,
) -> list[Path]:
    baseline = _comsol_numeric_matrix(baseline_path)
    external = _comsol_numeric_matrix(external_path)
    if min(baseline.shape[1], external.shape[1]) < 3:
        raise GecCcpPlotError("COMSOL waveform CSV needs coordinate, voltage, current")
    fig, axes = plt.subplots(2, 1, figsize=(9, 6), sharex=True, constrained_layout=True)
    for data, color, style, label in (
        (baseline, "#555555", "--", "built-in Druyvesteyn"),
        (external, "#1769aa", "-", "Swarm tables"),
    ):
        phase = _normalized_coordinate(data[:, 0])
        axes[0].plot(phase, data[:, 1], color=color, linestyle=style, label=label)
        axes[1].plot(phase, data[:, 2], color=color, linestyle=style, label=label)
    axes[0].set_ylabel("Electrode voltage [V]")
    axes[1].set_ylabel("Electrode current [A]")
    axes[1].set_xlabel("RF period fraction")
    axes[0].legend(frameon=False)
    for axis in axes:
        axis.grid(True, alpha=0.2)
    fig.suptitle("Argon GEC CCP powered-electrode waveform")
    return _save(fig, output_stem, plt)


def _save(fig: Any, stem: Path, plt: Any) -> list[Path]:
    paths = [stem.with_suffix(".png"), stem.with_suffix(".svg")]
    fig.savefig(paths[0], dpi=180, bbox_inches="tight")
    fig.savefig(paths[1], bbox_inches="tight")
    plt.close(fig)
    return paths


def _rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise GecCcpPlotError(f"missing plot input: {path}")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _numeric_rows(path: Path) -> list[dict[str, float]]:
    rows = _rows(path)
    parsed: list[dict[str, float]] = []
    for row in rows:
        parsed_row: dict[str, float] = {}
        for key, value in row.items():
            try:
                number = float(value)
            except (TypeError, ValueError):
                continue
            if math.isfinite(number):
                parsed_row[key] = number
        parsed.append(parsed_row)
    return parsed


def _comsol_numeric_matrix(path: Path) -> Any:
    try:
        import numpy as np
    except ImportError as exc:
        raise GecCcpPlotError("numpy is required") from exc
    lines = [
        line
        for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
        if line.strip() and not line.lstrip().startswith("%")
    ]
    if not lines:
        raise GecCcpPlotError(f"COMSOL CSV is empty: {path}")
    delimiter = "," if "," in lines[0] else None
    data = np.genfromtxt(lines, delimiter=delimiter, invalid_raise=False)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    data = data[~np.isnan(data).all(axis=1)]
    return data


def _nearest(values: list[float], target: float) -> float:
    return min(values, key=lambda value: abs(math.log(value / target)))


def _positive(value: object) -> bool:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(number) and number > 0.0


def _normalized_coordinate(values: Any) -> Any:
    low = values.min()
    span = values.max() - low
    return (values - low) / span if span > 0.0 else values * 0.0
