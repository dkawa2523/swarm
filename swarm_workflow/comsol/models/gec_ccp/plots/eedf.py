"""Operating-EEDF and tail-distribution rendering for GEC CCP results."""

from __future__ import annotations

from functools import partial
import json
import math
from pathlib import Path
from typing import Any

from swarm_workflow.comsol.models.gec_ccp.plots import contracts as _contracts
from swarm_workflow.comsol.models.gec_ccp.plots import data as _data
from swarm_workflow.comsol.input.function_eedf import (
    ComsolFunctionEedfGrid,
    evaluate_comsol_function_eedf_grid,
    FunctionEedfError,
    read_comsol_function_eedf_grid,
)


def plot_closure_overview(
    bundle: Path,
    output: Path,
    plt: Any,
    *,
    closure_source: str,
    reaction_model: str | None,
) -> list[Path]:
    """Render the bundle-level transport, rate, and EEDF overview."""

    import numpy as np

    transport = _data._numeric_rows(bundle / "transport_vs_mean_energy.csv")
    rates = _data._rows(bundle / "rates_vs_mean_energy.csv")
    eedf = _data._numeric_rows(bundle / "eedf.csv")
    if not transport or not rates or not eedf:
        raise _contracts.GecCcpPlotError("bundle lacks transport, rate, or EEDF data")
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    mean_energy = np.asarray([row["mean_energy_eV"] for row in transport])
    axes[0, 0].plot(
        mean_energy,
        [row["reduced_mobility_m2_V_s_m3"] for row in transport],
        color="#1769aa",
        linewidth=2,
        label=r"particle $\mu_eN$",
    )
    energy_mobility_column = "reduced_electron_energy_mobility_m2_V_s_m3"
    if not all(energy_mobility_column in row for row in transport):
        raise _contracts.GecCcpPlotError(
            "full Swarm closure lacks reduced electron-energy mobility"
        )
    axes[0, 0].plot(
        mean_energy,
        [row[energy_mobility_column] for row in transport],
        color="#c04b35",
        linewidth=2,
        linestyle="--",
        label=r"energy $\mu_\varepsilon N$",
    )
    axes[0, 0].set_ylabel(r"Reduced mobility [m$^{-1}$ V$^{-1}$ s$^{-1}$]")
    axes[0, 0].set_title("Particle and electron-energy mobility")
    axes[0, 0].legend(frameon=False)
    diffusion_columns = {
        "particle_L": "reduced_diffusion_L_m2_s_m3",
        "particle_T": "reduced_diffusion_T_m2_s_m3",
        "energy_L": "reduced_electron_energy_diffusion_L_m2_s_m3",
        "energy_T": "reduced_electron_energy_diffusion_T_m2_s_m3",
    }
    if not all(
        column in row for row in transport for column in diffusion_columns.values()
    ):
        raise _contracts.GecCcpPlotError(
            "full Swarm closure lacks particle/energy L,T diffusion"
        )
    diffusion = {
        name: np.asarray([row[column] for row in transport], dtype=float)
        for name, column in diffusion_columns.items()
    }
    if closure_source == "two_term":
        _plot_scalar_diffusion(axes[0, 1], mean_energy, diffusion, np)
    else:
        _plot_tensor_diffusion(axes[0, 1], mean_energy, diffusion)
    axes[0, 1].set_ylabel(r"Reduced diffusion [m$^{-1}$ s$^{-1}$]")
    axes[0, 1].set_title("Particle and electron-energy diffusion")
    axes[0, 1].legend(frameon=False, fontsize=8)
    _plot_rate_overview(axes[1, 0], rates, reaction_model)
    _plot_eedf_overview(axes[1, 1], eedf, closure_source)
    for axis in axes.flat:
        axis.set_xlabel("Mean energy [eV]" if axis is not axes[1, 1] else "Energy [eV]")
        axis.grid(True, which="both", alpha=0.2)
    for axis in (axes[0, 0], axes[0, 1], axes[1, 0]):
        axis.set_xscale("log")
    fig.suptitle(
        "Argon GEC CCP input closure\n"
        f"Steady DC {closure_source} tables supplied to the RF mean-energy model",
        fontsize=14,
    )
    return _data._save(fig, output / "swarm_closure_overview", plt)


def _plot_scalar_diffusion(
    axis: Any, mean_energy: Any, diffusion: dict[str, Any], np: Any
) -> None:
    if not (
        np.allclose(
            diffusion["particle_L"],
            diffusion["particle_T"],
            rtol=1.0e-12,
            atol=0.0,
        )
        and np.allclose(
            diffusion["energy_L"],
            diffusion["energy_T"],
            rtol=1.0e-12,
            atol=0.0,
        )
    ):
        raise _contracts.GecCcpPlotError("two_term scalar transport requires D_L = D_T")
    axis.plot(
        mean_energy,
        diffusion["particle_L"],
        color="#1769aa",
        linewidth=2,
        label=r"particle $D_LN=D_TN$ (scalar)",
    )
    axis.plot(
        mean_energy,
        diffusion["energy_L"],
        color="#c04b35",
        linewidth=2,
        linestyle="--",
        label=r"energy $D_{\varepsilon L}N=D_{\varepsilon T}N$ (scalar)",
    )


def _plot_tensor_diffusion(
    axis: Any,
    mean_energy: Any,
    diffusion: dict[str, Any],
) -> None:
    for name, color, linestyle, label in (
        ("particle_L", "#1769aa", "-", r"particle $D_LN$"),
        ("particle_T", "#1769aa", "--", r"particle $D_TN$"),
        ("energy_L", "#c04b35", "-", r"energy $D_{\varepsilon L}N$"),
        ("energy_T", "#c04b35", "--", r"energy $D_{\varepsilon T}N$"),
    ):
        axis.plot(
            mean_energy,
            diffusion[name],
            color=color,
            linewidth=2,
            linestyle=linestyle,
            label=label,
        )


def _plot_rate_overview(
    axis: Any, rates: list[dict[str, str]], reaction_model: str | None
) -> None:
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
            and _data._positive(row.get("rate_coefficient_m3_s"))
        ]
        axis.plot(
            [float(row["mean_energy_eV"]) for row in selected],
            [float(row["rate_coefficient_m3_s"]) for row in selected],
            color=rate_colors[process],
            linestyle=rate_styles[process],
            linewidth=2,
            label=process,
        )
    axis.set_yscale("log")
    axis.set_ylabel(r"Rate coefficient [m$^3$ s$^{-1}$]")
    if reaction_model == "function_eedf":
        axis.set_title("COMSOL cross-section rates (Swarm rate table: audit)")
    elif reaction_model == _contracts.HYBRID_FUNCTION_EEDF_REACTION_MODEL:
        axis.set_title("Elastic Function-EEDF integral; preintegrated inelastic rates")
    elif reaction_model == "external_rates":
        axis.set_title("Imported Swarm reaction rates")
    else:
        axis.set_title("Swarm reaction-rate table audit")
    axis.legend(frameon=False)


def _plot_eedf_overview(
    axis: Any,
    eedf: list[dict[str, float]],
    closure_source: str,
) -> None:
    available_en = sorted({row["E_over_N_Td"] for row in eedf})
    selected_en = sorted(
        {
            _data._nearest(available_en, value)
            for value in (10.0, 100.0, 1000.0, 10000.0)
        }
    )
    palette = ("#555555", "#1769aa", "#e69500", "#8f3985")
    for color, en_value in zip(palette, selected_en):
        selected = [
            row
            for row in eedf
            if math.isclose(row["E_over_N_Td"], en_value, rel_tol=1e-10)
            and row["eedf"] > 0.0
        ]
        axis.plot(
            [row["electron_energy_eV"] for row in selected],
            [row["eedf"] for row in selected],
            color=color,
            linewidth=1.8,
            label=f"{en_value:g} Td",
        )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_ylabel(r"EEDF [eV$^{-1}$]")
    axis.set_title(f"{closure_source} EEDF audit")
    axis.legend(frameon=False)


def plot_operating_eedf_if_available(
    bundle: Path,
    external_results: Path,
    output: Path,
    plt: Any,
    *,
    reaction_model: str | None,
    function_eedf_table: str | None,
) -> tuple[dict[str, Any], list[Path], list[str]]:
    """Render the operating EEDF when both required COMSOL tables exist."""

    closure_files = (
        external_results / "closure_phase.csv",
        external_results / "axis_phase_resolved.csv",
    )
    if not all(path.exists() for path in closure_files):
        return (
            {},
            [],
            ["operating EEDF audit: external closure/phase CSVs are absent"],
        )
    summary, figures = _plot_operating_eedf(
        bundle,
        closure_files[0],
        closure_files[1],
        output / "operating_eedf",
        plt,
        reaction_model=reaction_model,
        function_eedf_table=function_eedf_table,
    )
    return {"operating_eedf": summary}, figures, []


def _plot_operating_eedf(
    dc_bundle: Path,
    closure_phase_path: Path,
    axis_phase_path: Path,
    output_stem: Path,
    plt: Any,
    *,
    reaction_model: str | None,
    function_eedf_table: str | None,
) -> tuple[dict[str, Any], list[Path]]:
    """Plot EEDFs at electron-density-weighted operating quantiles.

    Function-EEDF runs use the canonical f0(E, meanE) table that COMSOL
    actually imports, converted to F(E)=sqrt(E)f0.  Other reaction modes use
    raw ``eedf.csv`` only as a Swarm diagnostic.
    """
    import numpy as np

    mean_energy = _data._comsol_field_matrix(closure_phase_path, "ptp.ebar")
    electron_density = _data._comsol_field_matrix(axis_phase_path, "ptp.ne")
    if mean_energy.shape != electron_density.shape:
        raise _contracts.GecCcpPlotError(
            "COMSOL ebar and ne phase samples have different shapes"
        )
    values = mean_energy.reshape(-1)
    weights = electron_density.reshape(-1)
    valid = (
        np.isfinite(values) & np.isfinite(weights) & (values > 0.0) & (weights > 0.0)
    )
    values = values[valid]
    weights = weights[valid]
    if values.size == 0 or float(np.sum(weights)) <= 0.0:
        raise _contracts.GecCcpPlotError(
            "COMSOL phase results contain no valid ebar/ne samples"
        )

    active_function_eedf = bool(
        reaction_model in _contracts.FUNCTION_EEDF_REACTION_MODELS
        and function_eedf_table
    )
    if active_function_eedf:
        active_table_path = dc_bundle / str(function_eedf_table)
        bundle_manifest = json.loads(
            (dc_bundle / "manifest.json").read_text(encoding="utf-8")
        )
        table_metadata = bundle_manifest.get("tables", {}).get(str(function_eedf_table))
        if not isinstance(table_metadata, dict):
            raise _contracts.GecCcpPlotError(
                "active Function-EEDF metadata is missing from bundle manifest"
            )
        try:
            active_eedf = read_comsol_function_eedf_grid(
                active_table_path,
                table_metadata,
            )
        except FunctionEedfError as exc:
            raise _contracts.GecCcpPlotError(str(exc)) from exc
        table_mean_energy = active_eedf.mean_energies_eV
        dc_slice_at = partial(_function_eedf_slice, active_eedf)
        dc_tail_curve_at = partial(_function_eedf_tail_curve, active_eedf)
        if reaction_model == _contracts.HYBRID_FUNCTION_EEDF_REACTION_MODEL:
            eedf_role = (
                "active COMSOL Function-EEDF input for elastic cross-section "
                "integration; excitation/ionization use preintegrated Swarm "
                "rates; plotted as F(E)=sqrt(E)*f0(E,meanE)"
            )
        else:
            eedf_role = (
                "active COMSOL Function-EEDF input; plotted as F(E)=sqrt(E)*f0(E,meanE)"
            )
        dc_legend_role = "active Function-EEDF"
    else:
        raw_rows = _data._numeric_rows(dc_bundle / "eedf.csv")
        if not raw_rows:
            raise _contracts.GecCcpPlotError("bundle raw eedf.csv is empty")
        table_mean_energy = np.asarray(
            sorted({row["mean_energy_eV"] for row in raw_rows}), dtype=float
        )
        dc_slice_at = partial(_nearest_eedf_slice, dc_bundle)
        dc_tail_curve_at = partial(_eedf_tail_curve, dc_bundle)
        eedf_role = (
            "raw Swarm diagnostic selected by COMSOL mean energy; not imposed "
            "when the reaction model does not use Function-EEDF"
        )
        dc_legend_role = "raw diagnostic"
    table_low = float(np.min(table_mean_energy))
    table_high = float(np.max(table_mean_energy))
    in_range = (values >= table_low) & (values <= table_high)
    quantile_probabilities = (0.01, 0.10, 0.50, 0.90, 0.99)
    quantile_values = _weighted_quantiles(values, weights, quantile_probabilities)
    quantiles = {
        f"q{int(probability * 100):02d}": float(value)
        for probability, value in zip(quantile_probabilities, quantile_values)
    }

    representative = (
        ("q10", quantiles["q10"]),
        ("q50", quantiles["q50"]),
        ("q90", quantiles["q90"]),
    )
    colors = ("#1769aa", "#e69500", "#8f3985")
    dc_tail_curves = {
        threshold: dc_tail_curve_at(threshold) for threshold in (11.5, 15.8)
    }
    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True, constrained_layout=True)
    slice_metrics: list[dict[str, Any]] = []
    for (label, target), color in zip(representative, colors):
        dc_slice = dc_slice_at(target)
        dc_positive = (dc_slice["energy_eV"] > 0.0) & (dc_slice["eedf_eV_inv"] > 0.0)
        axes[0].plot(
            dc_slice["energy_eV"][dc_positive],
            dc_slice["eedf_eV_inv"][dc_positive],
            color=color,
            linewidth=2,
            label=(
                f"{label} {dc_legend_role}: mean={dc_slice['mean_energy_eV']:.2f} eV"
            ),
        )
        axes[1].plot(
            dc_slice["energy_eV"],
            _tail_survival(dc_slice),
            color=color,
            linewidth=2,
        )
        entry: dict[str, Any] = {
            "quantile": label,
            "comsol_target_mean_energy_eV": target,
            "dc": _eedf_slice_metrics(dc_slice),
            "tail_probability_interpolated_at_target": {
                "dc_above_11_5_eV": _interpolate_tail_curve(
                    dc_tail_curves[11.5], target
                ),
                "dc_above_15_8_eV": _interpolate_tail_curve(
                    dc_tail_curves[15.8], target
                ),
            },
        }
        slice_metrics.append(entry)

    for axis in axes:
        axis.axvline(
            11.5,
            color="#777777",
            linestyle=":",
            linewidth=1.2,
            label="11.5 eV excitation threshold" if axis is axes[0] else None,
        )
        axis.axvline(
            15.8,
            color="#222222",
            linestyle="-.",
            linewidth=1.2,
            label="15.8 eV ionization threshold" if axis is axes[0] else None,
        )
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.grid(True, which="both", alpha=0.2)
    axes[0].set_ylabel(r"EEDF [eV$^{-1}$]")
    axes[1].set_ylabel("Probability above energy")
    axes[1].set_xlabel("Electron energy [eV]")
    axes[0].set_title(
        "Active Function-EEDF at COMSOL operating quantiles"
        if active_function_eedf
        else "Raw Swarm EEDF diagnostic at COMSOL operating quantiles"
    )
    axes[1].set_title("High-energy tail probability")
    axes[0].legend(frameon=False, ncol=2, fontsize=8)
    axes[0].set_xlim(5.0e-2, 60.0)
    axes[0].set_ylim(1.0e-12, 1.0)
    axes[1].set_ylim(1.0e-12, 1.1)
    fig.suptitle(
        "Axial x RF-phase samples, weighted by electron density\n"
        + (
            "solid: canonical Function-EEDF used by COMSOL"
            if active_function_eedf
            else "solid: raw steady-DC Swarm diagnostic"
        )
    )

    raw_fraction = float(np.mean(in_range))
    weighted_fraction = float(np.sum(weights[in_range]) / np.sum(weights))
    summary = {
        "sample_scope": "COMSOL center-axis nodes x RF phase",
        "nodes": int(mean_energy.shape[0]),
        "phase_samples": int(mean_energy.shape[1]),
        "valid_samples": int(values.size),
        "weighting": "electron density",
        "raw_min_mean_energy_eV": float(np.min(values)),
        "raw_max_mean_energy_eV": float(np.max(values)),
        "density_weighted_quantiles_eV": quantiles,
        "swarm_table_mean_energy_range_eV": [table_low, table_high],
        "raw_in_range_fraction": raw_fraction,
        "density_weighted_in_range_fraction": weighted_fraction,
        "coverage_status": (
            "pass_no_endpoint_hold_in_sampled_solution"
            if raw_fraction == 1.0
            else "fail_extend_swarm_table_before_acceptance"
        ),
        "additional_swarm_points_required": raw_fraction < 1.0,
        "eedf_role": eedf_role,
        "active_function_eedf_table": (
            str(dc_bundle / str(function_eedf_table)) if active_function_eedf else None
        ),
        "representative_slices": slice_metrics,
    }
    return summary, _data._save(fig, output_stem, plt)


def _weighted_quantiles(values: Any, weights: Any, probabilities: Any) -> Any:
    import numpy as np

    order = np.argsort(values, kind="stable")
    sorted_values = values[order]
    sorted_weights = weights[order]
    cumulative = np.cumsum(sorted_weights)
    cumulative = (cumulative - 0.5 * sorted_weights) / cumulative[-1]
    return np.interp(probabilities, cumulative, sorted_values)


def _nearest_eedf_slice(bundle: Path, target_mean_energy_eV: float) -> dict[str, Any]:
    import numpy as np

    rows = _data._numeric_rows(bundle / "eedf.csv")
    available = sorted({row["mean_energy_eV"] for row in rows})
    selected_mean = min(
        available,
        key=lambda value: abs(math.log(value / target_mean_energy_eV)),
    )
    selected = [
        row
        for row in rows
        if math.isclose(row["mean_energy_eV"], selected_mean, rel_tol=1.0e-10)
    ]
    selected.sort(key=lambda row: row["electron_energy_eV"])
    return {
        "mean_energy_eV": float(selected_mean),
        "E_over_N_Td": float(selected[0]["E_over_N_Td"]),
        "energy_eV": np.asarray(
            [row["electron_energy_eV"] for row in selected], dtype=float
        ),
        "width_eV": np.asarray(
            [row["energy_width_eV"] for row in selected], dtype=float
        ),
        "eedf_eV_inv": np.asarray([row["eedf"] for row in selected], dtype=float),
    }


def _function_eedf_slice(
    representation: ComsolFunctionEedfGrid,
    target_mean_energy_eV: float,
) -> dict[str, Any]:
    """Evaluate the exact physical Function-EEDF table imported by COMSOL."""
    import numpy as np

    means = representation.mean_energies_eV
    target = float(np.clip(target_mean_energy_eV, means[0], means[-1]))
    energies = representation.electron_energies_eV
    f0 = evaluate_comsol_function_eedf_grid(
        representation,
        energies,
        target,
    )
    distribution = np.sqrt(np.maximum(energies, 0.0)) * np.maximum(f0, 0.0)
    quadrature_width = np.empty_like(energies)
    quadrature_width[0] = 0.5 * (energies[1] - energies[0])
    quadrature_width[-1] = 0.5 * (energies[-1] - energies[-2])
    quadrature_width[1:-1] = 0.5 * (energies[2:] - energies[:-2])
    return {
        "mean_energy_eV": target,
        "E_over_N_Td": None,
        "energy_eV": energies,
        "width_eV": quadrature_width,
        "eedf_eV_inv": distribution,
        "source": "active_comsol_function_eedf_physical_grid",
    }


def _function_eedf_tail_curve(
    representation: ComsolFunctionEedfGrid,
    threshold_eV: float,
) -> tuple[Any, Any]:
    import numpy as np

    means = representation.mean_energies_eV.tolist()
    probabilities = [
        _tail_probability(_function_eedf_slice(representation, mean), threshold_eV)
        for mean in means
    ]
    return (
        np.asarray(means, dtype=float),
        np.maximum(np.asarray(probabilities, dtype=float), 1.0e-300),
    )


def _tail_survival(eedf_slice: dict[str, Any]) -> Any:
    import numpy as np

    mass = np.maximum(eedf_slice["eedf_eV_inv"] * eedf_slice["width_eV"], 0.0)
    total = float(np.sum(mass))
    return np.maximum(np.cumsum(mass[::-1])[::-1] / max(total, 1.0e-300), 1.0e-300)


def _eedf_slice_metrics(eedf_slice: dict[str, Any]) -> dict[str, float]:
    result: dict[str, Any] = {
        "table_mean_energy_eV": float(eedf_slice["mean_energy_eV"]),
        "E_over_N_Td": (
            None
            if eedf_slice.get("E_over_N_Td") is None
            else float(eedf_slice["E_over_N_Td"])
        ),
        "tail_probability_above_11_5_eV": _tail_probability(eedf_slice, 11.5),
        "tail_probability_above_15_8_eV": _tail_probability(eedf_slice, 15.8),
    }
    if "source" in eedf_slice:
        result["source"] = eedf_slice["source"]
    return result


def _tail_probability(eedf_slice: dict[str, Any], threshold_eV: float) -> float:
    import numpy as np

    mass = np.maximum(eedf_slice["eedf_eV_inv"] * eedf_slice["width_eV"], 0.0)
    total = float(np.sum(mass))
    selected = eedf_slice["energy_eV"] >= threshold_eV
    return float(np.sum(mass[selected]) / max(total, 1.0e-300))


def _eedf_tail_curve(
    bundle: Path,
    threshold_eV: float,
) -> tuple[Any, Any]:
    import numpy as np

    rows = _data._numeric_rows(bundle / "eedf.csv")
    mean_energies = sorted({row["mean_energy_eV"] for row in rows})
    probabilities: list[float] = []
    for mean_energy in mean_energies:
        selected = [
            row
            for row in rows
            if math.isclose(row["mean_energy_eV"], mean_energy, rel_tol=1.0e-10)
        ]
        mass = np.asarray(
            [max(row["eedf"], 0.0) * row["energy_width_eV"] for row in selected]
        )
        energy = np.asarray([row["electron_energy_eV"] for row in selected])
        total = float(np.sum(mass))
        probabilities.append(
            float(np.sum(mass[energy >= threshold_eV]) / max(total, 1.0e-300))
        )
    return (
        np.asarray(mean_energies, dtype=float),
        np.maximum(np.asarray(probabilities, dtype=float), 1.0e-300),
    )


def _interpolate_tail_curve(curve: tuple[Any, Any], target: float) -> float:
    import numpy as np

    x, y = curve
    return float(_positive_log_interpolate(np.asarray([target]), x, y)[0])


def _positive_log_interpolate(x: Any, known_x: Any, known_y: Any) -> Any:
    import numpy as np

    return np.exp(np.interp(np.log(x), np.log(known_x), np.log(known_y)))
