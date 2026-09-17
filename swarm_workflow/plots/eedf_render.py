"""Static rendering for EEDF distributions and comparison metrics."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from swarm_workflow.plots.eedf_contracts import COLORS, EedfCase, EedfDataset
from swarm_workflow.plots.eedf_data import case_at
from swarm_workflow.plots.eedf_metrics import half_density_on, mass_on, union_edges


def _stairs(
    ax: object,
    case: EedfCase,
    edges: np.ndarray,
    **kwargs: object,
) -> None:
    density = mass_on(case, edges) / np.diff(edges)
    ax.stairs(density, edges, **kwargs)


def _uncertainty_band(
    ax: object,
    case: EedfCase,
    edges: np.ndarray,
    *,
    floor: float,
    color: str,
    alpha: float,
) -> None:
    half = half_density_on(case, edges)
    if half is None:
        return
    density = mass_on(case, edges) / np.diff(edges)
    lower = np.maximum(density - half, floor)
    upper = np.maximum(density + half, floor)
    ax.fill_between(
        edges[:-1],
        lower,
        upper,
        step="post",
        color=color,
        alpha=alpha,
        linewidth=0.0,
    )


def _quantile(case: EedfCase, probability: float) -> float:
    mass = case.density_eV_inv * case.widths_eV / case.source_normalization
    cumulative = np.cumsum(mass)
    index = min(
        int(np.searchsorted(cumulative, probability, side="left")),
        len(mass) - 1,
    )
    prior = 0.0 if index == 0 else float(cumulative[index - 1])
    fraction = (probability - prior) / max(float(mass[index]), 1.0e-300)
    return float(
        case.edges_eV[index] + np.clip(fraction, 0.0, 1.0) * case.widths_eV[index]
    )


def plot_distributions(
    output: Path,
    representative_fields: Sequence[float],
    statuses: Mapping[float, bool],
    mc: EedfDataset,
    two_term: EedfDataset,
    propagator: EedfDataset,
    tail_thresholds_eV: Sequence[float],
) -> tuple[Path, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    panel_cases = [
        case_at(dataset, field)
        for field in representative_fields
        for dataset in (mc, two_term, propagator)
    ]
    y_max = (
        max(
            float(np.max(case.density_eV_inv / case.source_normalization))
            for case in panel_cases
        )
        * 1.6
    )
    y_floor = max(y_max * 1.0e-9, 1.0e-12)
    columns = 2
    rows = int(math.ceil(len(representative_fields) / columns))
    figure_height = 4.25 * rows + 1.35
    fig, axes = plt.subplots(
        rows,
        columns,
        figsize=(12.8, figure_height),
        squeeze=False,
    )
    for ax, field in zip(axes.ravel(), representative_fields, strict=False):
        tt_case = case_at(two_term, field)
        prop_case = case_at(propagator, field)
        mc_case = case_at(mc, field)
        edges = union_edges((tt_case, prop_case, mc_case))
        panel_x_max = max(
            max(_quantile(case, 0.999999) for case in (tt_case, prop_case, mc_case))
            * 1.25,
            max(tail_thresholds_eV, default=0.0) * 1.25,
        )
        _stairs(
            ax,
            tt_case,
            edges,
            color=COLORS["two_term"],
            linewidth=2.2,
            linestyle=(0, (2, 1.5)),
            label="two-term",
        )
        _stairs(
            ax,
            prop_case,
            edges,
            color=COLORS["propagator"],
            linewidth=1.9,
            linestyle=(0, (5, 2)),
            label="propagator P1",
        )
        _uncertainty_band(
            ax,
            mc_case,
            edges,
            floor=y_floor,
            color=COLORS["monte_carlo"],
            alpha=0.16,
        )
        _stairs(
            ax,
            mc_case,
            edges,
            color=COLORS["monte_carlo"],
            linewidth=2.0,
            label="Monte Carlo (raw aggregate)",
        )
        source_note = (
            "MC active-closure gate: "
            f"{'passed' if statuses[field] else 'not passed'}\n"
            "raw MC retained in comparison"
        )
        for threshold in tail_thresholds_eV:
            ax.axvline(
                float(threshold),
                color="#9CA3AF",
                linewidth=0.7,
                alpha=0.65,
            )
        ax.set_yscale("log")
        ax.set_xlim(0.0, panel_x_max)
        ax.set_ylim(y_floor, y_max)
        ax.grid(
            axis="both",
            which="major",
            color="#D1D5DB",
            linewidth=0.6,
            alpha=0.65,
        )
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(
            f"E/N = {field:g} Td",
            loc="left",
            fontsize=11,
            fontweight="bold",
        )
        ax.text(
            0.98,
            0.96,
            source_note,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8.3,
            color="#374151",
        )
        ax.set_xlabel("electron energy, E (eV)")
        ax.set_ylabel("probability-density EEDF, F(E) (eV$^{-1}$)")
    for ax in axes.ravel()[len(representative_fields) :]:
        ax.axis("off")
    legend = [
        Line2D(
            [],
            [],
            color=COLORS["monte_carlo"],
            linewidth=2.0,
            label="Monte Carlo (raw aggregate)",
        ),
        Patch(
            facecolor=COLORS["monte_carlo"],
            alpha=0.16,
            label="MC pointwise 95% t interval",
        ),
        Line2D(
            [],
            [],
            color=COLORS["two_term"],
            linewidth=2.2,
            linestyle=(0, (2, 1.5)),
            label="two-term",
        ),
        Line2D(
            [],
            [],
            color=COLORS["propagator"],
            linewidth=1.9,
            linestyle=(0, (5, 2)),
            label="propagator P1",
        ),
    ]
    fig.legend(
        handles=legend,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.0 - 0.72 / figure_height),
        ncol=4,
        frameon=False,
        fontsize=9,
    )
    fig.suptitle(
        "EEDF comparison across independent raw solver outputs",
        x=0.07,
        y=1.0 - 0.08 / figure_height,
        ha="left",
        va="top",
        fontsize=16,
        fontweight="bold",
    )
    fig.text(
        0.07,
        1.0 - 0.44 / figure_height,
        "All curves are normalized F(E); each panel uses a common conservative "
        "grid and a field-specific energy range.",
        ha="left",
        va="top",
        fontsize=9.5,
        color="#4B5563",
    )
    fig.subplots_adjust(
        left=0.08,
        right=0.98,
        bottom=0.58 / figure_height,
        top=1.0 - 1.35 / figure_height,
        hspace=0.31,
        wspace=0.18,
    )
    png = output / "eedf_physical_comparison.png"
    svg = output / "eedf_physical_comparison.svg"
    fig.savefig(png, dpi=220, facecolor="white")
    fig.savefig(svg, facecolor="white")
    plt.close(fig)
    return png, svg


def plot_metrics(
    output: Path,
    rows: Sequence[Mapping[str, object]],
) -> tuple[Path, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pairs = {
        "raw_mc_vs_two_term": (
            "raw MC vs two-term",
            COLORS["monte_carlo"],
            "-",
            "o",
        ),
        "propagator_vs_two_term": (
            "propagator P1 vs two-term",
            COLORS["propagator"],
            (0, (5, 2)),
            "s",
        ),
        "raw_mc_vs_propagator": (
            "raw MC vs propagator P1",
            COLORS["mc_vs_propagator"],
            (0, (2, 1.5)),
            "^",
        ),
    }
    metrics = (
        ("total_variation", "physical-energy total variation"),
        ("shape_total_variation", "mean-scaled shape total variation"),
        (
            "reconstructed_mean_energy_relative_difference",
            "relative mean-energy difference",
        ),
    )
    fig, axes = plt.subplots(3, 1, figsize=(11.5, 10.2), sharex=True)
    for ax, (metric, label) in zip(axes, metrics, strict=True):
        for pair, (pair_label, color, linestyle, marker) in pairs.items():
            selected_rows = sorted(
                (row for row in rows if row["comparison"] == pair),
                key=lambda row: float(row["E_over_N_Td"]),
            )
            x = np.asarray([float(row["E_over_N_Td"]) for row in selected_rows])
            y = np.asarray([float(row[metric]) for row in selected_rows])
            ax.plot(
                x,
                y,
                color=color,
                linewidth=1.8,
                linestyle=linestyle,
                marker=marker,
                markersize=4,
                label=pair_label,
            )
        ax.set_ylabel(label)
        ax.set_ylim(bottom=0.0)
        ax.grid(axis="y", color="#D1D5DB", linewidth=0.7, alpha=0.7)
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xscale("log")
    axes[-1].set_xlabel("E/N (Td)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.925),
        ncol=2,
        frameon=False,
        fontsize=9,
    )
    fig.suptitle(
        "EEDF differences across the common reduced-field anchors",
        x=0.08,
        y=0.985,
        ha="left",
        fontsize=16,
        fontweight="bold",
    )
    fig.text(
        0.08,
        0.95,
        "Every metric compares independent raw solver outputs. Total variation "
        "is bounded in [0, 1].",
        ha="left",
        fontsize=9.5,
        color="#4B5563",
    )
    fig.subplots_adjust(
        left=0.12,
        right=0.98,
        bottom=0.08,
        top=0.86,
        hspace=0.18,
    )
    png = output / "eedf_difference_metrics.png"
    svg = output / "eedf_difference_metrics.svg"
    fig.savefig(png, dpi=220, facecolor="white")
    fig.savefig(svg, facecolor="white")
    plt.close(fig)
    return png, svg
