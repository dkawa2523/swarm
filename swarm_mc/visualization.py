# Copyright (c) 2020-2021 ETH Zurich
"""
Convenience plotting helpers that operate directly on the CSV outputs.
"""

from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plot_time_series(csv_path: str, show: bool = True, block: bool = True,
                     save_path: str = None) -> plt.Figure:
    """
    Plot key quantities from ``*_temporal_evolution.csv``.
    """
    df = pd.read_csv(csv_path)

    fig, axes = plt.subplots(2, 2, figsize=(10, 6))
    df.plot(x="time", y="num_electrons", ax=axes[0, 0], label="electrons")
    df.plot(x="time", y="num_cations", ax=axes[0, 0], label="cations")
    df.plot(x="time", y="num_anions", ax=axes[0, 0], label="anions")
    axes[0, 0].set_ylabel("count")

    df.plot(x="time", y="mean_energy", ax=axes[0, 1], color="C1")
    axes[0, 1].set_ylabel("mean energy (eV)")

    df.plot(x="time", y="num_collisions", ax=axes[1, 0], color="C2")
    axes[1, 0].set_ylabel("collisions")

    df.plot(x="time", y="mean_position_z", ax=axes[1, 1], color="C3",
            label="z-position")
    axes[1, 1].set_ylabel("mean position z (m)")

    for ax in axes.ravel():
        ax.grid(True)

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    if not show:
        plt.close(fig)

    if show:
        plt.show(block=block)
    return fig


def _safe_positive(values, floor: float = 1e-30):
    arr = values.to_numpy()
    return np.clip(arr, floor, None)


def _log_limits(values, floor_ratio: float = 1e-2, padding: float = 1.5):
    arr = np.asarray(values)
    pos = arr[np.isfinite(arr) & (arr > 0)]
    if pos.size == 0:
        return (1e-30, 1)
    lo = pos.min()
    hi = pos.max()
    lo = max(lo, hi * floor_ratio)
    return (lo / padding, hi * padding)


def plot_energy_distribution(csv_path: str, show: bool = True,
                             block: bool = True, save_path: str = None) -> plt.Figure:
    """
    Plot eedf and eepf from ``*_energy_distribution.csv`` with log-scaled y axes.
    """
    df = pd.read_csv(csv_path)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    eedf = _safe_positive(df["eedf (eV-1)"])
    axes[0].plot(df["energy (eV)"], eedf)
    axes[0].set_xlabel("energy (eV)")
    axes[0].set_ylabel("eedf (eV-1)")
    axes[0].set_yscale("log")
    axes[0].set_ylim(_log_limits(df["eedf (eV-1)"]))
    axes[0].grid(True, which="both")

    eepf = _safe_positive(df["eepf (eV-3/2)"])
    axes[1].plot(df["energy (eV)"], eepf)
    axes[1].set_xlabel("energy (eV)")
    axes[1].set_ylabel("eepf (eV-3/2)")
    axes[1].set_yscale("log")
    axes[1].set_ylim(_log_limits(df["eepf (eV-3/2)"]))
    axes[1].grid(True, which="both")

    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    if not show:
        plt.close(fig)
    if show:
        plt.show(block=block)
    return fig


def plot_energy_distribution_combined(csv_path: str, show: bool = True,
                                      block: bool = True, save_path: str = None
                                      ) -> plt.Figure:
    """
    Plot EEDF and EEPF together using twin y-axes (both log scaled).
    """
    df = pd.read_csv(csv_path)
    energy = df["energy (eV)"]
    eedf = _safe_positive(df["eedf (eV-1)"])
    eepf = _safe_positive(df["eepf (eV-3/2)"])

    fig, ax1 = plt.subplots(figsize=(7, 4))
    l1, = ax1.plot(energy, eedf, color="C0", label="EEDF")
    ax1.set_xlabel("energy (eV)")
    ax1.set_ylabel("eedf (eV-1)", color="C0")
    ax1.set_yscale("log")
    ax1.set_ylim(_log_limits(df["eedf (eV-1)"]))
    ax1.tick_params(axis="y", labelcolor="C0")
    ax1.grid(True, which="both")

    ax2 = ax1.twinx()
    l2, = ax2.plot(energy, eepf, color="C1", label="EEPF")
    ax2.set_ylabel("eepf (eV-3/2)", color="C1")
    ax2.set_yscale("log")
    ax2.set_ylim(_log_limits(df["eepf (eV-3/2)"]))
    ax2.tick_params(axis="y", labelcolor="C1")

    fig.tight_layout()
    fig.legend(handles=[l1, l2], loc="upper right")
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    if not show:
        plt.close(fig)
    if show:
        plt.show(block=block)
    return fig


def plot_summary(csv_path: str, x: str = "E/N (Td)",
                 y: str = "flux drift velocity (m.s-1)",
                 show: bool = True, block: bool = True,
                 save_path: str = None) -> Optional[plt.Figure]:
    """
    Plot a sweep summary (e.g. drift velocity vs reduced electric field).
    """
    path = Path(csv_path)
    if not path.exists():
        return None

    df = pd.read_csv(path)
    fig, ax = plt.subplots(figsize=(7, 4))
    df.plot(x=x, y=y, style="o-", ax=ax)
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.grid(True)
    fig.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
    if not show:
        plt.close(fig)
    if show:
        plt.show(block=block)
    return fig


def plot_summary_metrics(csv_path: str, metrics, show: bool = False,
                         block: bool = True, save_dir: Optional[Path] = None):
    """
    Plot multiple summary metrics vs E/N. Metrics is iterable of
    (column_name, filename_suffix, ylabel).
    Returns list of generated figures.
    """
    path = Path(csv_path)
    if not path.exists():
        return []
    df = pd.read_csv(path)
    figures = []
    for col, fname, ylabel in metrics:
        if col not in df.columns:
            continue
        fig, ax = plt.subplots(figsize=(7, 4))
        df.plot(x="E/N (Td)", y=col, style="o-", ax=ax)
        ax.set_xlabel("E/N (Td)")
        ax.set_ylabel(ylabel)
        ax.grid(True)
        fig.tight_layout()
        if save_dir:
            save_dir.mkdir(parents=True, exist_ok=True)
            fig.savefig(Path(save_dir) / fname, bbox_inches="tight")
        if not show:
            plt.close(fig)
        if show:
            plt.show(block=block)
        figures.append(fig)
    return figures
