"""Static figures for validated GEC-ICP saved-solution comparisons."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

from .comparison_analysis_contracts import (
    GecIcpComparisonCaseData,
    METRIC_KEYS,
    METRIC_LABELS,
)


_PALETTE = ("#1F2937", "#2563EB", "#D97706", "#4D7C0F")


def render_time_history(
    output_directory: Path,
    cases: Sequence[GecIcpComparisonCaseData],
    *,
    common_time_s: float,
) -> tuple[Path, Path]:
    """Plot all absolute histories without pretending terminal times align."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 2, figsize=(12.4, 11.0), squeeze=False)
    inventory_metrics = {
        "electron_inventory",
        "metastable_inventory",
        "ion_inventory",
    }
    for ax, key in zip(axes.ravel(), METRIC_KEYS, strict=True):
        for index, case in enumerate(cases):
            x_ms = [row.time_s * 1.0e3 for row in case.time_series]
            y = [float(getattr(row, key)) for row in case.time_series]
            ax.plot(
                x_ms,
                y,
                color=_PALETTE[index],
                linewidth=2.0 if index == 0 else 1.7,
                marker="o",
                markersize=2.7,
                label=case.label,
            )
            ax.scatter(
                [case.terminal.time_s * 1.0e3],
                [float(getattr(case.terminal, key))],
                color=_PALETTE[index],
                marker="s",
                s=26,
                zorder=4,
            )
        if key in inventory_metrics:
            ax.set_yscale("log")
        ax.set_xscale("symlog", linthresh=1.0e-4, linscale=0.5)
        ax.set_xlim(left=0.0)
        ax.axvline(
            common_time_s * 1.0e3,
            color="#9CA3AF",
            linewidth=1.0,
            linestyle=(0, (4, 3)),
        )
        ax.set_title(METRIC_LABELS[key], loc="left", fontsize=11, fontweight="bold")
        ax.set_xlabel("physical time (ms)")
        ax.grid(which="major", color="#D1D5DB", linewidth=0.65, alpha=0.7)
        ax.spines[["top", "right"]].set_visible(False)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.928),
        ncol=min(len(cases), 4),
        frameon=False,
        fontsize=9.5,
    )
    fig.suptitle(
        "GEC-ICP saved-solution time histories",
        x=0.07,
        y=0.99,
        ha="left",
        fontsize=16,
        fontweight="bold",
    )
    fig.text(
        0.07,
        0.955,
        f"Dashed line: exact {common_time_s * 1.0e3:g} ms common comparison "
        "time; square: each model's own terminal state.",
        ha="left",
        fontsize=9.5,
        color="#4B5563",
    )
    fig.subplots_adjust(
        left=0.09,
        right=0.98,
        bottom=0.07,
        top=0.87,
        hspace=0.32,
        wspace=0.22,
    )
    png = output_directory / "icp_time_history.png"
    svg = output_directory / "icp_time_history.svg"
    fig.savefig(
        png, dpi=220, facecolor="white", metadata={"Software": "swarm_workflow"}
    )
    fig.savefig(svg, facecolor="white", metadata={"Creator": "swarm_workflow"})
    plt.close(fig)
    return png, svg


def render_common_time_relative_change(
    output_directory: Path,
    cases: Sequence[GecIcpComparisonCaseData],
    *,
    baseline_case_id: str,
    common_time_s: float,
) -> tuple[Path, Path]:
    """Plot signed changes at the exact common time only."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    baseline = next(case for case in cases if case.case_id == baseline_case_id)
    candidates = [case for case in cases if case.case_id != baseline_case_id]
    x = np.arange(len(METRIC_KEYS), dtype=float)
    group_width = 0.78
    width = group_width / len(candidates)
    fig, ax = plt.subplots(figsize=(12.3, 6.4))
    for index, case in enumerate(candidates):
        changes = [
            100.0
            * (float(getattr(case.common, key)) - float(getattr(baseline.common, key)))
            / abs(float(getattr(baseline.common, key)))
            for key in METRIC_KEYS
        ]
        offset = (index - 0.5 * (len(candidates) - 1)) * width
        ax.bar(
            x + offset,
            changes,
            width=width * 0.9,
            color=_PALETTE[index + 1],
            label=case.label,
        )
    ax.axhline(0.0, color="#111827", linewidth=0.9)
    ax.set_xticks(x, [METRIC_LABELS[key] for key in METRIC_KEYS])
    ax.tick_params(axis="x", labelrotation=20)
    for tick in ax.get_xticklabels():
        tick.set_horizontalalignment("right")
    ax.set_ylabel(f"signed change from {baseline.label} (%)")
    ax.grid(axis="y", color="#D1D5DB", linewidth=0.7, alpha=0.75)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper right", frameon=False, ncol=min(len(candidates), 3))
    fig.suptitle(
        f"Like-for-like comparison at {common_time_s * 1.0e3:g} ms",
        x=0.075,
        y=0.98,
        ha="left",
        fontsize=16,
        fontweight="bold",
    )
    fig.text(
        0.075,
        0.93,
        "All bars use the same physical time. Positive values exceed the saved "
        f"{baseline.label} baseline.",
        ha="left",
        fontsize=9.5,
        color="#4B5563",
    )
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.25, top=0.86)
    png = output_directory / "icp_common_time_relative_change.png"
    svg = output_directory / "icp_common_time_relative_change.svg"
    fig.savefig(
        png, dpi=220, facecolor="white", metadata={"Software": "swarm_workflow"}
    )
    fig.savefig(svg, facecolor="white", metadata={"Creator": "swarm_workflow"})
    plt.close(fig)
    return png, svg


__all__ = ("render_common_time_relative_change", "render_time_history")
