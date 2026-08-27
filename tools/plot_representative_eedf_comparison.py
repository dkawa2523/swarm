"""Plot representative MC, two-term, and COMSOL analytic EEDFs.

The swarm bundle stores the probability-density convention ``F(E)`` in
``1/eV``. COMSOL documents its assumed Maxwellian and Druyvesteyn functions
as ``f0(E)`` in ``eV^(-3/2)``. This script converts the analytic functions via
``F(E) = sqrt(E) * f0(E)`` before plotting them beside the swarm results.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


REPRESENTATIVE_EN_TD = (20.0, 100.0, 300.0, 1000.0, 4000.0)
COMSOL_EEDF_REFERENCE = (
    "https://doc.comsol.com/6.4/doc/com.comsol.help.plasma/"
    "plasma_ug_boltzmann.06.09.html"
)

COLORS = {
    "monte_carlo": "#2563EB",
    "two_term": "#E07A1F",
    "maxwellian": "#4B5563",
    "druyvesteyn": "#7C6A0A",
}


@dataclass(frozen=True)
class EedfCase:
    solver: str
    e_over_n_td: float
    energy_eV: np.ndarray
    widths_eV: np.ndarray
    probability_eV_inv: np.ndarray
    reported_mean_energy_eV: float
    source_normalization: float
    calculated_mean_energy_eV: float


def comsol_f0_eV_m32(
    energy_eV: np.ndarray,
    mean_energy_eV: float,
    power_g: float,
) -> np.ndarray:
    """Return COMSOL's generalized analytic ``f0(E)`` convention."""

    energy = np.asarray(energy_eV, dtype=float)
    if mean_energy_eV <= 0.0:
        raise ValueError("mean energy must be positive")
    if not 1.0 <= power_g <= 2.0:
        raise ValueError("COMSOL analytic EEDF power g must be in [1, 2]")
    gamma_5 = math.gamma(5.0 / (2.0 * power_g))
    gamma_3 = math.gamma(3.0 / (2.0 * power_g))
    beta_1 = gamma_5**1.5 * gamma_3**-2.5
    beta_2 = gamma_5 / gamma_3
    scaled = np.maximum(energy, 0.0) * beta_2 / mean_energy_eV
    return (
        power_g
        * mean_energy_eV**-1.5
        * beta_1
        * np.exp(-(scaled**power_g))
    )


def comsol_probability_eedf_eV_inv(
    energy_eV: np.ndarray,
    mean_energy_eV: float,
    power_g: float,
) -> np.ndarray:
    """Return COMSOL's analytic EEDF as normalized probability density F(E)."""

    energy = np.asarray(energy_eV, dtype=float)
    return np.sqrt(np.maximum(energy, 0.0)) * comsol_f0_eV_m32(
        energy,
        mean_energy_eV,
        power_g,
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _format_en(e_over_n_td: float) -> str:
    if float(e_over_n_td).is_integer():
        return f"{int(e_over_n_td):04d}"
    return str(e_over_n_td).replace(".", "p")


def _read_cases(
    bundle: Path,
    solver: str,
    representative_en_td: Iterable[float],
) -> dict[float, EedfCase]:
    eedf_path = bundle / "eedf.csv"
    mean_path = bundle / "mean_energy_vs_en.csv"
    eedf = pd.read_csv(eedf_path)
    means = pd.read_csv(mean_path)
    cases: dict[float, EedfCase] = {}
    for requested in representative_en_td:
        matching = means[np.isclose(means["E_over_N_Td"], requested)]
        if len(matching) != 1:
            raise ValueError(
                f"{solver} bundle does not contain exactly one {requested:g} Td case"
            )
        reported_mean = float(matching.iloc[0]["mean_energy_eV"])
        rows = eedf[np.isclose(eedf["E_over_N_Td"], requested)].copy()
        if rows.empty:
            raise ValueError(f"{solver} EEDF lacks {requested:g} Td")
        rows = rows.sort_values("electron_energy_eV")
        energy = rows["electron_energy_eV"].to_numpy(dtype=float)
        widths = rows["energy_width_eV"].to_numpy(dtype=float)
        probability = rows["eedf"].to_numpy(dtype=float)
        if (
            np.any(~np.isfinite(energy))
            or np.any(~np.isfinite(widths))
            or np.any(~np.isfinite(probability))
            or np.any(np.diff(energy) <= 0.0)
            or np.any(widths <= 0.0)
            or np.any(probability < 0.0)
        ):
            raise ValueError(f"{solver} {requested:g} Td EEDF is invalid")
        source_norm = float(np.sum(probability * widths))
        if source_norm <= 0.0:
            raise ValueError(f"{solver} {requested:g} Td EEDF has zero mass")
        normalized = probability / source_norm
        calculated_mean = float(np.sum(energy * normalized * widths))
        cases[float(requested)] = EedfCase(
            solver=solver,
            e_over_n_td=float(requested),
            energy_eV=energy,
            widths_eV=widths,
            probability_eV_inv=normalized,
            reported_mean_energy_eV=reported_mean,
            source_normalization=source_norm,
            calculated_mean_energy_eV=calculated_mean,
        )
    return cases


def _plot_positive(ax, x: np.ndarray, y: np.ndarray, *args, **kwargs) -> None:
    mask = np.isfinite(x) & np.isfinite(y) & (y > 0.0)
    ax.plot(x[mask], y[mask], *args, **kwargs)


def _configure_axes(ax, *, xlabel: str, ylabel: str, x_max: float) -> None:
    ax.set_yscale("log")
    ax.set_xlim(0.0, x_max)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", color="#D1D5DB", linewidth=0.7, alpha=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#6B7280")
    ax.spines["bottom"].set_color("#6B7280")
    ax.tick_params(colors="#374151", labelsize=9)
    ax.tick_params(axis="x", labelbottom=True)


def _plot_detail_panel(
    ax,
    mc: EedfCase,
    two_term: EedfCase,
    reference_mean_eV: float,
    reference_label: str,
    x_max: float,
) -> None:
    energy = np.linspace(0.0, x_max, 1600)
    maxwell = comsol_probability_eedf_eV_inv(energy, reference_mean_eV, 1.0)
    druyvesteyn = comsol_probability_eedf_eV_inv(energy, reference_mean_eV, 2.0)
    _plot_positive(
        ax,
        mc.energy_eV,
        mc.probability_eV_inv,
        color=COLORS["monte_carlo"],
        linewidth=2.0,
        label="Monte Carlo",
    )
    _plot_positive(
        ax,
        two_term.energy_eV,
        two_term.probability_eV_inv,
        color=COLORS["two_term"],
        linewidth=2.0,
        label="two-term",
    )
    _plot_positive(
        ax,
        energy,
        maxwell,
        color=COLORS["maxwellian"],
        linewidth=1.8,
        linestyle="--",
        label="COMSOL Maxwellian",
    )
    _plot_positive(
        ax,
        energy,
        druyvesteyn,
        color=COLORS["druyvesteyn"],
        linewidth=1.8,
        linestyle=(0, (1.5, 1.5)),
        label="COMSOL Druyvesteyn",
    )
    _configure_axes(
        ax,
        xlabel="Electron energy, ε (eV)",
        ylabel="Probability-density EEDF, F(ε) (eV⁻¹)",
        x_max=x_max,
    )
    peak = max(
        float(np.max(mc.probability_eV_inv)),
        float(np.max(two_term.probability_eV_inv)),
        float(np.max(maxwell)),
        float(np.max(druyvesteyn)),
    )
    ax.set_ylim(max(peak * 1.0e-7, 1.0e-10), peak * 1.6)
    ax.set_title(
        f"COMSOL analytic mean = {reference_label}: {reference_mean_eV:.4g} eV\n"
        f"Actual EEDF means: MC {mc.reported_mean_energy_eV:.4g} eV | "
        f"2T {two_term.reported_mean_energy_eV:.4g} eV",
        loc="left",
        color="#111827",
        fontsize=11,
        fontweight="bold",
    )


def _plot_detail(
    output_dir: Path,
    e_over_n_td: float,
    mc: EedfCase,
    two_term: EedfCase,
) -> tuple[Path, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x_max = 12.0 * max(
        mc.reported_mean_energy_eV,
        two_term.reported_mean_energy_eV,
    )
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 6.1), sharex=True, sharey=True)
    _plot_detail_panel(
        axes[0],
        mc,
        two_term,
        mc.reported_mean_energy_eV,
        "MC mean",
        x_max,
    )
    _plot_detail_panel(
        axes[1],
        mc,
        two_term,
        two_term.reported_mean_energy_eV,
        "two-term mean",
        x_max,
    )
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.855),
        ncol=4,
        frameon=False,
        fontsize=9.5,
    )
    fig.suptitle(
        f"EEDF comparison at E/N = {e_over_n_td:g} Td",
        x=0.065,
        y=0.98,
        ha="left",
        fontsize=15,
        fontweight="bold",
        color="#111827",
    )
    fig.text(
        0.065,
        0.91,
        "COMSOL analytic f₀ converted to F(ε)=√ε f₀; both panels use the same physical EEDF data.",
        ha="left",
        fontsize=9.5,
        color="#4B5563",
    )
    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.105, top=0.69, wspace=0.12)
    stem = f"eedf_comparison_{_format_en(e_over_n_td)}_Td"
    png = output_dir / f"{stem}.png"
    svg = output_dir / f"{stem}.svg"
    fig.savefig(png, dpi=220, facecolor="white")
    fig.savefig(svg, facecolor="white")
    plt.close(fig)
    return png, svg


def _plot_overview(
    output_dir: Path,
    representative_en_td: tuple[float, ...],
    mc_cases: dict[float, EedfCase],
    two_term_cases: dict[float, EedfCase],
) -> tuple[Path, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 2, figsize=(12.8, 13.8), sharex=True, sharey=True)
    flat_axes = axes.ravel()
    scaled_energy = np.linspace(0.0, 12.0, 1400)
    maxwell = comsol_probability_eedf_eV_inv(scaled_energy, 1.0, 1.0)
    druyvesteyn = comsol_probability_eedf_eV_inv(scaled_energy, 1.0, 2.0)
    for ax, e_over_n_td in zip(flat_axes, representative_en_td, strict=False):
        mc = mc_cases[e_over_n_td]
        two_term = two_term_cases[e_over_n_td]
        _plot_positive(
            ax,
            mc.energy_eV / mc.reported_mean_energy_eV,
            mc.reported_mean_energy_eV * mc.probability_eV_inv,
            color=COLORS["monte_carlo"],
            linewidth=2.0,
            label="Monte Carlo",
        )
        _plot_positive(
            ax,
            two_term.energy_eV / two_term.reported_mean_energy_eV,
            two_term.reported_mean_energy_eV * two_term.probability_eV_inv,
            color=COLORS["two_term"],
            linewidth=2.0,
            label="two-term",
        )
        _plot_positive(
            ax,
            scaled_energy,
            maxwell,
            color=COLORS["maxwellian"],
            linewidth=1.8,
            linestyle="--",
            label="COMSOL Maxwellian",
        )
        _plot_positive(
            ax,
            scaled_energy,
            druyvesteyn,
            color=COLORS["druyvesteyn"],
            linewidth=1.8,
            linestyle=(0, (1.5, 1.5)),
            label="COMSOL Druyvesteyn",
        )
        _configure_axes(
            ax,
            xlabel="Scaled energy, ε/⟨ε⟩",
            ylabel="Scaled EEDF, ⟨ε⟩F(ε)",
            x_max=12.0,
        )
        ax.set_ylim(1.0e-7, 1.2)
        ax.set_title(
            f"E/N = {e_over_n_td:g} Td",
            loc="left",
            fontsize=12,
            fontweight="bold",
            color="#111827",
        )
        ax.text(
            0.98,
            0.95,
            f"⟨ε⟩ MC {mc.reported_mean_energy_eV:.3g} eV\n"
            f"⟨ε⟩ 2T {two_term.reported_mean_energy_eV:.3g} eV",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8.5,
            color="#4B5563",
        )
    for ax in flat_axes[len(representative_en_td) :]:
        ax.axis("off")
        ax.text(
            0.04,
            0.92,
            "How to read the overview",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=13,
            fontweight="bold",
            color="#111827",
        )
        ax.text(
            0.04,
            0.78,
            "• Each MC and two-term curve uses its own ⟨ε⟩.\n\n"
            "• Maxwellian and Druyvesteyn are COMSOL's analytic\n"
            "  shapes with the same zeroth and first moments.\n\n"
            "• Separation after scaling is a shape difference, not\n"
            "  a mean-energy difference.\n\n"
            "• The individual files show physical energy in eV.",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=10.5,
            linespacing=1.25,
            color="#4B5563",
        )
    handles, labels = flat_axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.925),
        ncol=4,
        frameon=False,
        fontsize=10,
    )
    fig.suptitle(
        "EEDF shape comparison at representative reduced fields",
        x=0.075,
        y=0.982,
        ha="left",
        fontsize=17,
        fontweight="bold",
        color="#111827",
    )
    fig.text(
        0.075,
        0.948,
        "Each physical EEDF is scaled by its own mean energy; analytic COMSOL shapes are therefore uniquely moment matched.",
        ha="left",
        fontsize=10,
        color="#4B5563",
    )
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.07, top=0.88, hspace=0.28, wspace=0.16)
    png = output_dir / "eedf_comparison_representative_en.png"
    svg = output_dir / "eedf_comparison_representative_en.svg"
    fig.savefig(png, dpi=220, facecolor="white")
    fig.savefig(svg, facecolor="white")
    plt.close(fig)
    return png, svg


def _analytic_audit(mean_energy_eV: float, power_g: float) -> dict[str, float]:
    upper = 100.0 * mean_energy_eV
    energy = np.concatenate(
        (
            np.asarray([0.0]),
            np.geomspace(mean_energy_eV * 1.0e-12, upper, 200_000),
        )
    )
    probability = comsol_probability_eedf_eV_inv(
        energy,
        mean_energy_eV,
        power_g,
    )
    normalization = float(np.trapezoid(probability, energy))
    calculated_mean = float(np.trapezoid(energy * probability, energy))
    return {
        "normalization": normalization,
        "normalization_error": abs(normalization - 1.0),
        "calculated_mean_energy_eV": calculated_mean,
        "mean_energy_relative_error": abs(calculated_mean - mean_energy_eV)
        / mean_energy_eV,
    }


def generate(
    mc_bundle: Path,
    two_term_bundle: Path,
    output_dir: Path,
    representative_en_td: tuple[float, ...] = REPRESENTATIVE_EN_TD,
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    mc_cases = _read_cases(mc_bundle, "monte_carlo", representative_en_td)
    two_term_cases = _read_cases(two_term_bundle, "two_term", representative_en_td)
    outputs: list[Path] = []
    overview = _plot_overview(
        output_dir,
        representative_en_td,
        mc_cases,
        two_term_cases,
    )
    outputs.extend(overview)
    summary_rows: list[dict[str, float | str]] = []
    detail_files: dict[str, list[str]] = {}
    for e_over_n_td in representative_en_td:
        mc = mc_cases[e_over_n_td]
        two_term = two_term_cases[e_over_n_td]
        detail = _plot_detail(output_dir, e_over_n_td, mc, two_term)
        outputs.extend(detail)
        detail_files[f"{e_over_n_td:g}"] = [path.name for path in detail]
        summary_rows.append(
            {
                "E_over_N_Td": e_over_n_td,
                "mc_reported_mean_energy_eV": mc.reported_mean_energy_eV,
                "mc_calculated_mean_energy_eV": mc.calculated_mean_energy_eV,
                "mc_source_normalization": mc.source_normalization,
                "two_term_reported_mean_energy_eV": two_term.reported_mean_energy_eV,
                "two_term_calculated_mean_energy_eV": two_term.calculated_mean_energy_eV,
                "two_term_source_normalization": two_term.source_normalization,
                "solver_mean_energy_relative_difference": abs(
                    mc.reported_mean_energy_eV
                    - two_term.reported_mean_energy_eV
                )
                / two_term.reported_mean_energy_eV,
                "detail_png": detail[0].name,
                "detail_svg": detail[1].name,
            }
        )
    summary_path = output_dir / "eedf_comparison_summary.csv"
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    outputs.append(summary_path)
    analytic_audits = {
        "maxwellian_g1": _analytic_audit(1.0, 1.0),
        "druyvesteyn_g2": _analytic_audit(1.0, 2.0),
    }
    maximum_source_norm_error = max(
        abs(case.source_normalization - 1.0)
        for case in (*mc_cases.values(), *two_term_cases.values())
    )
    maximum_source_mean_error = max(
        abs(case.calculated_mean_energy_eV - case.reported_mean_energy_eV)
        / case.reported_mean_energy_eV
        for case in (*mc_cases.values(), *two_term_cases.values())
    )
    manifest = {
        "title": "Representative EEDF comparison for current COMSOL/two-term/MC spatial case",
        "representative_E_over_N_Td": list(representative_en_td),
        "selection_basis": (
            "common MC and two-term anchors spanning low, intermediate, inelastic, "
            "and high reduced-field regimes"
        ),
        "convention": {
            "plotted_quantity": "F(electron_energy) probability-density EEDF",
            "plotted_units": "eV^-1",
            "normalization": "integral F(E) dE = 1",
            "mean_energy": "integral E F(E) dE",
            "comsol_conversion": "F(E) = sqrt(E) * f0(E)",
            "comsol_native_f0_units": "eV^-3/2",
        },
        "comsol_analytic_reference": {
            "documentation": COMSOL_EEDF_REFERENCE,
            "generalized_formula": (
                "f0(E)=g*phi^(-3/2)*beta1*exp(-(E*beta2/phi)^g); "
                "beta1=Gamma(5/(2g))^(3/2)*Gamma(3/(2g))^(-5/2); "
                "beta2=Gamma(5/(2g))/Gamma(3/(2g))"
            ),
            "maxwellian_power_g": 1.0,
            "druyvesteyn_power_g": 2.0,
            "phi_definition": "mean electron energy in eV",
            "detail_plot_moment_basis": (
                "left panel matches the MC mean; right panel matches the two-term mean"
            ),
            "overview_moment_basis": (
                "each physical curve is nondimensionalized by its own mean energy"
            ),
        },
        "sources": {
            "monte_carlo": {
                "path": str((mc_bundle / "eedf.csv").resolve()),
                "sha256": _sha256(mc_bundle / "eedf.csv"),
                "mean_energy_path": str(
                    (mc_bundle / "mean_energy_vs_en.csv").resolve()
                ),
                "mean_energy_sha256": _sha256(mc_bundle / "mean_energy_vs_en.csv"),
            },
            "two_term": {
                "path": str((two_term_bundle / "eedf.csv").resolve()),
                "sha256": _sha256(two_term_bundle / "eedf.csv"),
                "mean_energy_path": str(
                    (two_term_bundle / "mean_energy_vs_en.csv").resolve()
                ),
                "mean_energy_sha256": _sha256(
                    two_term_bundle / "mean_energy_vs_en.csv"
                ),
            },
        },
        "qa": {
            "maximum_source_normalization_error": maximum_source_norm_error,
            "maximum_source_mean_energy_relative_error": maximum_source_mean_error,
            "analytic_distribution_audits": analytic_audits,
        },
        "outputs": {
            "overview": [path.name for path in overview],
            "details_by_E_over_N_Td": detail_files,
            "summary": summary_path.name,
        },
    }
    manifest_path = output_dir / "eedf_comparison_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    outputs.append(manifest_path)
    return outputs


def _parser(repo_root: Path) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    base = repo_root / "outputs" / "argon_gec_ccp"
    parser.add_argument(
        "--mc-bundle",
        type=Path,
        default=base / "monte_carlo_function_eedf_restricted_lmea_bundle" / "mixture_0000",
    )
    parser.add_argument(
        "--two-term-bundle",
        type=Path,
        default=base / "two_term_function_eedf_restricted_lmea_bundle" / "mixture_0000",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=base / "plots_current_comsol_two_term_mc_spatial",
    )
    parser.add_argument(
        "--representative-en",
        type=float,
        nargs="+",
        default=list(REPRESENTATIVE_EN_TD),
        metavar="TD",
    )
    return parser


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    args = _parser(repo_root).parse_args()
    outputs = generate(
        args.mc_bundle,
        args.two_term_bundle,
        args.output_dir,
        tuple(float(value) for value in args.representative_en),
    )
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
