"""Plot helpers for the unified output schema."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from electron_swarm.core.config import OutputConfig
from electron_swarm.core.results import SwarmRunResult
from electron_swarm.core.solver_registry import solver_plot_color, solver_plot_label


def _positive_curve(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mask = np.isfinite(x) & np.isfinite(y) & (y > 0.0)
    if not np.any(mask):
        return np.array([], dtype=float), np.array([], dtype=float)
    return np.asarray(x[mask], dtype=float), np.asarray(y[mask], dtype=float)


def _truncate_mc_curve(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mask = np.isfinite(x) & np.isfinite(y)
    if not np.any(mask):
        return np.array([], dtype=float), np.array([], dtype=float)
    x = np.asarray(x[mask], dtype=float)
    y = np.asarray(y[mask], dtype=float)
    positive_idx = np.flatnonzero(y > 0.0)
    if positive_idx.size == 0:
        return np.array([], dtype=float), np.array([], dtype=float)
    start = int(positive_idx[0])
    nonpositive_after_start = np.flatnonzero(y[start:] <= 0.0)
    stop = start + int(nonpositive_after_start[0]) if nonpositive_after_start.size else len(y)
    return x[start:stop], y[start:stop]


def write_plots(result: SwarmRunResult, output: OutputConfig) -> dict[str, Path]:
    if not output.write_plots:
        return {}
    output.directory.mkdir(parents=True, exist_ok=True)
    summary = pd.DataFrame([case.summary_dict() for case in result.cases])
    paths: dict[str, Path] = {}
    if summary.empty:
        return paths

    def line_plot(
        y: str, ylabel: str, filename: str, yscale: str | None = None
    ) -> None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        for solver, frame in summary.groupby("solver"):
            frame = frame.sort_values("E_over_N_Td")
            ax.plot(frame["E_over_N_Td"], frame[y], marker="o", label=solver)
        ax.set_xlabel("E/N [Td]")
        ax.set_ylabel(ylabel)
        if yscale:
            ax.set_yscale(yscale)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend()
        fig.tight_layout()
        path = output.directory / filename
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths[filename] = path

    line_plot(
        "mean_energy_eV",
        "Mean electron energy [eV]",
        f"{output.base_name}_mean_energy.png",
    )
    line_plot(
        "drift_velocity_m_s",
        "Drift velocity [m/s]",
        f"{output.base_name}_drift_velocity.png",
    )
    line_plot(
        "reduced_mobility_m2_V_s_m3",
        "Reduced mobility mu*N [m^-1 V^-1 s^-1]",
        f"{output.base_name}_reduced_mobility.png",
    )
    line_plot(
        "reduced_diffusion_L_m2_s_m3",
        "Reduced longitudinal diffusion D_L*N [m^-1 s^-1]",
        f"{output.base_name}_reduced_diffusion.png",
    )

    eedf_cases = sorted(result.cases, key=lambda case: (case.e_over_n_Td, case.solver))
    unique_en = sorted({float(case.e_over_n_Td) for case in eedf_cases})
    ncols = min(3, max(1, len(unique_en)))
    nrows = int(math.ceil(len(unique_en) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(5.2 * ncols, 3.8 * nrows),
        sharex=False,
        sharey=True,
    )
    axes_arr = np.atleast_1d(axes).ravel()

    ymax = max(
        (
            float(np.max(case.eedf[np.isfinite(case.eedf) & (case.eedf > 0.0)]))
            for case in eedf_cases
            if np.any(np.isfinite(case.eedf) & (case.eedf > 0.0))
        ),
        default=1.0e-5,
    )
    ymin = 1.0e-6
    ymax = max(ymax * 1.15, ymin * 10.0)

    for ax, e_over_n_Td in zip(axes_arr, unique_en, strict=False):
        cases_at_en = [
            case for case in eedf_cases if float(case.e_over_n_Td) == float(e_over_n_Td)
        ]
        for case in sorted(cases_at_en, key=lambda item: item.solver):
            color = solver_plot_color(case.solver)
            label = solver_plot_label(case.solver)
            if case.solver == "monte_carlo":
                x, y = _truncate_mc_curve(case.energy_eV, case.eedf)
                if x.size == 0:
                    continue
                ax.plot(
                    x,
                    y,
                    color=color,
                    linewidth=1.9,
                    linestyle="--",
                    marker="s",
                    markersize=3.5,
                    markevery=max(1, len(x) // 20),
                    label=label,
                )
            else:
                x, y = _positive_curve(case.energy_eV, case.eedf)
                if x.size == 0:
                    continue
                ax.plot(
                    x,
                    y,
                    color=color,
                    linewidth=2.2,
                    label=label,
                )
        ax.set_title(f"{e_over_n_Td:g} Td")
        ax.set_xlabel("Electron energy [eV]")
        ax.set_ylabel("EEDF [eV^-1]")
        ax.set_yscale("log")
        ax.set_ylim(ymin, ymax)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize="small")

    for ax in axes_arr[len(unique_en) :]:
        ax.set_visible(False)

    fig.suptitle("EEDF by E/N")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.97))
    path = output.directory / f"{output.base_name}_eedf.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths["eedf_png"] = path
    return paths
