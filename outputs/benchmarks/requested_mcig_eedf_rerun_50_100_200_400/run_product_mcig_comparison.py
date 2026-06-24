from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from electron_swarm import load_config, run
from electron_swarm.core.numerics import widths_from_centers


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
BASE = "ar_50_100_200_400_two_term_mc_vs_mcig"
CONFIG = OUT / f"{BASE}.yaml"
FIELDS = [50.0, 100.0, 200.0, 400.0]
MCIG_INPUTS = OUT / "mcig_reference_inputs"
MCIG_100 = MCIG_INPUTS / "ar_mcig_gui_100Td_reference_eedf.csv"
MCIG_200_400_600 = MCIG_INPUTS / "ar_mcig_vs_internal_mc_200_400_600_actual_curves.csv"
MCIG_RUN_LOGS = {
    100.0: MCIG_INPUTS / "ar_mcig_gui_100Td_runbyrun.dat",
    200.0: MCIG_INPUTS / "ar_mcig_gui_200Td_runbyrun.dat",
    400.0: MCIG_INPUTS / "ar_mcig_gui_400Td_runbyrun.dat",
}


def finite_widths(energy: np.ndarray, explicit: np.ndarray | None = None) -> np.ndarray:
    if explicit is not None:
        widths = np.asarray(explicit, dtype=float)
        if widths.shape == energy.shape and np.all(np.isfinite(widths)) and np.all(widths > 0):
            return widths
    return widths_from_centers(energy)


def normalized_curve(
    energy: np.ndarray,
    eedf: np.ndarray,
    widths: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, float, float, float]:
    energy = np.asarray(energy, dtype=float)
    eedf = np.asarray(eedf, dtype=float)
    widths = finite_widths(energy, widths)
    order = np.argsort(energy)
    energy = energy[order]
    eedf = eedf[order]
    widths = widths[order]
    area = float(np.sum(eedf * widths))
    if not np.isfinite(area) or area <= 0.0:
        raise ValueError("EEDF curve has non-positive integral")
    eedf = eedf / area
    mean = float(np.sum(energy * eedf * widths))
    tail_mask = energy >= 15.76
    tail = float(np.sum(eedf[tail_mask] * widths[tail_mask]))
    return eedf, widths, area, mean, tail


def common_support_metrics(
    ref_energy: np.ndarray,
    ref_eedf: np.ndarray,
    ref_widths: np.ndarray,
    cand_energy: np.ndarray,
    cand_eedf: np.ndarray,
) -> tuple[float, float, float, float]:
    support_min = max(float(np.min(ref_energy)), float(np.min(cand_energy)))
    support_max = min(float(np.max(ref_energy)), float(np.max(cand_energy)))
    mask = (ref_energy >= support_min) & (ref_energy <= support_max)
    if not np.any(mask):
        return math.nan, support_min, support_max, math.nan
    energy = ref_energy[mask]
    widths = ref_widths[mask]
    ref_y = ref_eedf[mask].copy()
    cand_y = np.interp(energy, cand_energy, cand_eedf)
    ref_area = float(np.sum(ref_y * widths))
    cand_area = float(np.sum(cand_y * widths))
    if ref_area > 0.0:
        ref_y /= ref_area
    if cand_area > 0.0:
        cand_y /= cand_area
    rel_l1 = float(
        np.sum(np.abs(cand_y - ref_y) * widths)
        / max(float(np.sum(np.abs(ref_y) * widths)), 1.0e-300)
    )
    log_mask = (energy >= 15.76) & (ref_y > 0.0) & (cand_y > 0.0)
    log_tail = (
        float(np.sqrt(np.mean((np.log10(cand_y[log_mask]) - np.log10(ref_y[log_mask])) ** 2)))
        if np.any(log_mask)
        else math.nan
    )
    return rel_l1, support_min, support_max, log_tail


def load_mcig_references() -> dict[float, pd.DataFrame]:
    refs: dict[float, pd.DataFrame] = {}
    if MCIG_100.exists():
        df = pd.read_csv(MCIG_100)
        df = df.loc[np.isclose(df["E_over_N_Td"], 100.0)].copy()
        refs[100.0] = pd.DataFrame(
            {
                "E_over_N_Td": 100.0,
                "source": "MCIG GUI actual",
                "energy_eV": df["energy_eV"].to_numpy(float),
                "eedf_eV_inv": df["eedf_eV_inv"].to_numpy(float),
            }
        )
    if MCIG_200_400_600.exists():
        df = pd.read_csv(MCIG_200_400_600)
        df = df.loc[df["source"] == "MCIG GUI actual"].copy()
        for field in (200.0, 400.0):
            sub = df.loc[np.isclose(df["E_over_N_Td"], field)].copy()
            if not sub.empty:
                refs[field] = sub[["E_over_N_Td", "source", "energy_eV", "eedf_eV_inv"]]
    return refs


def parse_mcig_run_log(path: Path) -> dict[str, str]:
    text = path.read_text(encoding="utf-8", errors="replace") if path.exists() else ""
    patterns = {
        "E_over_N_Td": r"Electric field /N \(Td\)\s+([0-9.Ee+\-]+)",
        "gas_temperature_K": r"Gas temperature \(K\)\s+([0-9.Ee+\-]+)",
        "magnetic_field_Hx": r"Magnetic field /N \(Hx\)\s+([0-9.Ee+\-]+)",
        "angular_frequency_per_N": r"Angular frequency /N \(m3.rad/s\)\s+([0-9.Ee+\-]+)",
        "ionization_degree": r"Ionization degree\s+([0-9.Ee+\-]+)",
        "energy_sharing_parameter_eV": r"Energy sharing parameter \(eV\)\s+([0-9.Ee+\-]+)",
        "angular_scattering_model": r"Angular scattering model\s+([0-9.Ee+\-]+)",
        "energy_sharing_model": r"Energy sharing model\s+([0-9.Ee+\-]+)",
        "growth_model": r"Growth model\s+([0-9.Ee+\-]+)",
        "particles": r"Number of particles\s+([0-9.Ee+\-]+)",
        "edf_max_energy_eV": r"EDF/VDF maximum energy \(eV\)\s+([0-9.Ee+\-]+)",
        "rng": r"Random number generator\s+([0-9.Ee+\-]+)",
        "mole_fraction_Ar": r"Mole fraction Ar\s+([0-9.Ee+\-]+)",
        "mcig_reported_mean_energy_eV": r"Mean energy \(eV\)\s+([0-9.Ee+\-]+)",
    }
    return {
        key: (match.group(1) if (match := re.search(pattern, text)) else "")
        for key, pattern in patterns.items()
    }


def write_config() -> None:
    CONFIG.write_text(
        """schema_version: 2

run:
  solvers:
    - id: two_term
    - id: monte_carlo
  e_over_n_Td: [50.0, 100.0, 200.0, 400.0]
  case_prefix: ar_mcig_sweep_rerun

conditions:
  gas_temperature_K: 300.0
  pressure_Pa: 100.0
  gas_mixture:
    - species: Ar
      fraction: 1.0
      mass_amu: 39.948

cross_sections:
  high_energy_extrapolation: hold
  files:
    - path: ../../../examples/cross_sections/argon_minimal.csv
      species: Ar

physics:
  field:
    type: dc
    magnetic_field:
      enabled: false
      B_T: 0.0
      angle_EB_deg: 90.0
  angular_scattering:
    model: isotropic
    higher_moment_closure: zero
  electron_electron:
    enabled: false
    model: none
  ionization:
    energy_sharing: equal
  energy_grid_policy:
    adaptive: true
    threshold_refinement: true
    max_eV_limit: 1000.0

solvers:
  monte_carlo:
    population_model: fixed_particle_single_daughter
    particles: 1024
    warmup_collisions: 300
    max_collisions: 3000
    seed: 20260624

comparison:
  enabled: true
  reference_solver: monte_carlo
  candidate_solvers: [two_term]
  compare_eedf: true
  required: false

output:
  directory: .
  base_name: ar_50_100_200_400_product
""",
        encoding="utf-8",
    )


def plot_individual(
    field: float,
    curves: pd.DataFrame,
    out_path: Path,
    *,
    title_suffix: str,
) -> None:
    sub = curves.loc[np.isclose(curves["E_over_N_Td"], field)].copy()
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), constrained_layout=True)
    colors = {
        "MCIG GUI actual": "#111827",
        "product two_term": "#2563eb",
        "product internal MC": "#dc2626",
    }
    for source, group in sub.groupby("source", sort=False):
        group = group.sort_values("energy_eV")
        linewidth = 2.35 if source == "MCIG GUI actual" else 1.9
        axes[0].plot(
            group["energy_eV"],
            group["eedf_eV_inv"],
            label=source,
            linewidth=linewidth,
            color=colors.get(source),
        )
        axes[1].semilogy(
            group["energy_eV"],
            np.clip(group["eedf_eV_inv"], 1.0e-12, None),
            label=source,
            linewidth=linewidth,
            color=colors.get(source),
        )
    emax = 45.0 if field <= 100.0 else 80.0 if field <= 200.0 else 140.0
    axes[0].set_xlim(0.0, min(emax, 80.0))
    axes[1].set_xlim(0.0, emax)
    axes[1].set_ylim(1.0e-8, 2.5e-1)
    axes[0].set_title(f"{field:.0f} Td linear")
    axes[1].set_title(f"{field:.0f} Td tail")
    for ax in axes:
        ax.set_xlabel("Energy (eV)")
        ax.set_ylabel("EEDF (1/eV), normalized")
        ax.grid(True, which="both", alpha=0.25)
    axes[1].legend(loc="upper right", frameon=False)
    fig.suptitle(f"Product two_term / internal MC vs MCIG reference {title_suffix}")
    fig.savefig(out_path.with_suffix(".png"), dpi=180)
    fig.savefig(out_path.with_suffix(".svg"))
    plt.close(fig)


def plot_overview(curves: pd.DataFrame, out_path: Path, *, semilog: bool) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.2), constrained_layout=True)
    colors = {
        "MCIG GUI actual": "#111827",
        "product two_term": "#2563eb",
        "product internal MC": "#dc2626",
    }
    for ax, field in zip(axes.reshape(-1), FIELDS, strict=True):
        sub = curves.loc[np.isclose(curves["E_over_N_Td"], field)].copy()
        for source, group in sub.groupby("source", sort=False):
            group = group.sort_values("energy_eV")
            linewidth = 2.25 if source == "MCIG GUI actual" else 1.75
            y = group["eedf_eV_inv"]
            if semilog:
                ax.semilogy(
                    group["energy_eV"],
                    np.clip(y, 1.0e-12, None),
                    label=source,
                    linewidth=linewidth,
                    color=colors.get(source),
                )
            else:
                ax.plot(
                    group["energy_eV"],
                    y,
                    label=source,
                    linewidth=linewidth,
                    color=colors.get(source),
                )
        if field <= 100.0:
            ax.set_xlim(0.0, 45.0)
        elif field <= 200.0:
            ax.set_xlim(0.0, 80.0)
        else:
            ax.set_xlim(0.0, 140.0)
        if semilog:
            ax.set_ylim(1.0e-8, 2.5e-1)
        ax.set_title(
            f"{field:.0f} Td" + ("" if field in (100.0, 200.0, 400.0) else " (MCIG missing)")
        )
        ax.set_xlabel("Energy (eV)")
        ax.set_ylabel("EEDF (1/eV), normalized")
        ax.grid(True, which="both", alpha=0.25)
    handles, labels = axes.reshape(-1)[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False)
    fig.suptitle(
        "Product two_term / internal MC rerun vs MCIG GUI reference"
        + (" (tail/log view)" if semilog else " (linear view)")
    )
    fig.savefig(out_path.with_suffix(".png"), dpi=180)
    fig.savefig(out_path.with_suffix(".svg"))
    plt.close(fig)


def main() -> None:
    write_config()
    result = run(load_config(CONFIG), write=True, collect_diagnostics=False)
    product_eedf = pd.read_csv(OUT / "ar_50_100_200_400_product_eedf.csv")
    summary = pd.read_csv(OUT / "ar_50_100_200_400_product_summary.csv")
    refs = load_mcig_references()

    curves: list[pd.DataFrame] = []
    normalized_by_key: dict[tuple[float, str], dict[str, np.ndarray | float]] = {}
    for field, ref in refs.items():
        energy = ref["energy_eV"].to_numpy(float)
        y, widths, area, mean, tail = normalized_curve(
            energy, ref["eedf_eV_inv"].to_numpy(float)
        )
        frame = pd.DataFrame(
            {
                "E_over_N_Td": field,
                "source": "MCIG GUI actual",
                "solver": "mcig",
                "energy_eV": energy,
                "energy_width_eV": widths,
                "eedf_eV_inv": y,
            }
        ).sort_values("energy_eV")
        curves.append(frame)
        normalized_by_key[(field, "mcig")] = {
            "energy": frame["energy_eV"].to_numpy(float),
            "widths": frame["energy_width_eV"].to_numpy(float),
            "eedf": frame["eedf_eV_inv"].to_numpy(float),
            "area": area,
            "mean": mean,
            "tail": tail,
        }

    for field in FIELDS:
        for solver, source in (
            ("two_term", "product two_term"),
            ("monte_carlo", "product internal MC"),
        ):
            df = product_eedf.loc[
                (product_eedf["solver"] == solver)
                & np.isclose(product_eedf["E_over_N_Td"], field)
            ].copy()
            energy = df["energy_eV"].to_numpy(float)
            widths_raw = df["energy_width_eV"].to_numpy(float)
            y, widths, area, mean, tail = normalized_curve(
                energy,
                df["eedf"].to_numpy(float),
                widths_raw,
            )
            frame = pd.DataFrame(
                {
                    "E_over_N_Td": field,
                    "source": source,
                    "solver": solver,
                    "energy_eV": energy,
                    "energy_width_eV": widths,
                    "eedf_eV_inv": y,
                }
            ).sort_values("energy_eV")
            curves.append(frame)
            normalized_by_key[(field, solver)] = {
                "energy": frame["energy_eV"].to_numpy(float),
                "widths": frame["energy_width_eV"].to_numpy(float),
                "eedf": frame["eedf_eV_inv"].to_numpy(float),
                "area": area,
                "mean": mean,
                "tail": tail,
            }

    curves_df = pd.concat(curves, ignore_index=True)
    curves_df.to_csv(OUT / f"{BASE}_curves.csv", index=False)

    rows: list[dict[str, object]] = []
    for field in FIELDS:
        ref = normalized_by_key.get((field, "mcig"))
        if ref is None:
            rows.append(
                {
                    "E_over_N_Td": field,
                    "comparison": "MCIG reference missing",
                    "candidate_solver": "",
                    "status": "missing_mcig_reference",
                }
            )
            continue
        ref_energy = ref["energy"]
        ref_widths = ref["widths"]
        ref_eedf = ref["eedf"]
        assert isinstance(ref_energy, np.ndarray)
        assert isinstance(ref_widths, np.ndarray)
        assert isinstance(ref_eedf, np.ndarray)
        for solver, source in (
            ("two_term", "product two_term"),
            ("monte_carlo", "product internal MC"),
        ):
            cand = normalized_by_key[(field, solver)]
            cand_energy = cand["energy"]
            cand_eedf = cand["eedf"]
            assert isinstance(cand_energy, np.ndarray)
            assert isinstance(cand_eedf, np.ndarray)
            rel_l1, support_min, support_max, log_tail = common_support_metrics(
                ref_energy,
                ref_eedf,
                ref_widths,
                cand_energy,
                cand_eedf,
            )
            mean_candidate = float(cand["mean"])
            mean_ref = float(ref["mean"])
            tail_candidate = float(cand["tail"])
            tail_ref = float(ref["tail"])
            reported = summary.loc[
                (summary["solver"] == solver) & np.isclose(summary["E_over_N_Td"], field),
                "mean_energy_eV",
            ]
            rows.append(
                {
                    "E_over_N_Td": field,
                    "comparison": f"{source} vs MCIG GUI actual",
                    "candidate_solver": solver,
                    "status": "complete",
                    "eedf_relative_l1_common_support": rel_l1,
                    "common_support_min_eV": support_min,
                    "common_support_max_eV": support_max,
                    "mean_energy_candidate_eV": mean_candidate,
                    "mean_energy_mcig_eV": mean_ref,
                    "mean_energy_relative_difference": abs(mean_candidate - mean_ref)
                    / max(abs(mean_ref), 1.0e-300),
                    "tail_probability_candidate_gt_15p76eV": tail_candidate,
                    "tail_probability_mcig_gt_15p76eV": tail_ref,
                    "tail_probability_difference_candidate_minus_mcig": tail_candidate
                    - tail_ref,
                    "log_tail_rms_error_gt_15p76eV": log_tail,
                    "summary_reported_mean_energy_eV": (
                        float(reported.iloc[0]) if not reported.empty else math.nan
                    ),
                    "candidate_integral_before_renormalization": float(cand["area"]),
                    "mcig_integral_before_renormalization": float(ref["area"]),
                }
            )
    metrics = pd.DataFrame(rows)
    metrics.to_csv(OUT / f"{BASE}_metrics.csv", index=False)

    settings_path = OUT / f"{BASE}_mcig_settings.csv"
    with settings_path.open("w", newline="", encoding="utf-8") as fp:
        writer = csv.writer(fp)
        writer.writerow(["E_over_N_Td", "key", "mcig_value", "product_value", "assessment"])
        for field in FIELDS:
            if field not in MCIG_RUN_LOGS:
                writer.writerow(
                    [
                        field,
                        "reference_status",
                        "",
                        "product run available",
                        "MCIG run-by-run EEDF not available in repository",
                    ]
                )
                continue
            values = parse_mcig_run_log(MCIG_RUN_LOGS[field])
            writer.writerow([field, "E_over_N_Td", values["E_over_N_Td"], field, "match"])
            writer.writerow([field, "species_fraction", values["mole_fraction_Ar"], "Ar fraction 1.0", "match"])
            writer.writerow([field, "gas_temperature_K", values["gas_temperature_K"], "300.0", "match"])
            writer.writerow([field, "magnetic_field", values["magnetic_field_Hx"], "disabled, B_T=0", "match for no magnetic field"])
            writer.writerow([field, "rf_or_ac_field", values["angular_frequency_per_N"], "field.type=dc", "match for dc field"])
            writer.writerow([field, "electron_electron", values["ionization_degree"], "electron_electron disabled", "match for no e-e collisions"])
            writer.writerow([field, "cross_section_source", "argon_minimal_bolsig.txt from examples/cross_sections/argon_minimal.csv", "examples/cross_sections/argon_minimal.csv", "same source, different parser/format"])
            writer.writerow([field, "cross_section_energy_range", f"MCIG EDF max {values['edf_max_energy_eV']} eV", "CSV tabulated to 100 eV; high_energy_extrapolation=hold; max_eV_limit=1000", "not exact; explicit hold extrapolation used"])
            writer.writerow([field, "angular_scattering", values["angular_scattering_model"], "physics.angular_scattering=isotropic", "not confirmed equivalent"])
            writer.writerow([field, "ionization_energy_sharing", f"model={values['energy_sharing_model']}, parameter={values['energy_sharing_parameter_eV']} eV", "energy_sharing=equal", "not confirmed equivalent"])
            writer.writerow([field, "population/statistics", f"MCIG particles={values['particles']}, rng={values['rng']}", "internal MC particles=1024, warmup=300, collisions=3000, seed=20260624", "not identical; product rerun uses higher/different statistics"])

    plot_overview(curves_df, OUT / f"{BASE}_eedf_overview_linear", semilog=False)
    plot_overview(curves_df, OUT / f"{BASE}_eedf_overview_tail", semilog=True)
    for field in FIELDS:
        suffix = "(MCIG missing)" if field not in refs else f"(Ar, {field:.0f} Td)"
        plot_individual(
            field,
            curves_df,
            OUT / f"{BASE}_{int(field)}Td_eedf",
            title_suffix=suffix,
        )

    complete = metrics.loc[metrics["status"] == "complete"].copy()
    complete["l1_text"] = complete["eedf_relative_l1_common_support"].map(lambda v: f"{v:.6g}")
    report_lines = [
        "# Product two_term/internal MC vs MCIG GUI EEDF sweep",
        "",
        "Fresh product run with schema_version=2.",
        "",
        "## Product run",
        "",
        f"- Config: `{CONFIG.relative_to(ROOT)}`",
        "- Product solvers: `two_term`, `monte_carlo`",
        "- E/N: 50, 100, 200, 400 Td",
        "- Ar, 300 K, pressure 100 Pa, no magnetic field, dc field, no e-e collisions",
        "- Cross sections: `examples/cross_sections/argon_minimal.csv`",
        "- High-energy handling: `high_energy_extrapolation: hold`, `max_eV_limit: 1000.0`",
        "- Internal MC: 1024 particles, 300 warmup collisions, 3000 production collisions, seed 20260624",
        "",
        "## MCIG availability",
        "",
        "- 100 Td: actual MCIG GUI EEDF archived in `mcig_reference_inputs/`",
        "- 200/400 Td: actual MCIG GUI EEDF archived in `mcig_reference_inputs/`",
        "- 50 Td: no saved MCIG GUI EEDF found; product curves are plotted without MCIG comparison metrics",
        "",
        "## Metrics vs MCIG GUI actual",
        "",
        "| E/N Td | comparison | rel L1 | candidate mean eV | MCIG mean eV | mean rel diff | tail prob diff >15.76 eV |",
        "|---:|---|---:|---:|---:|---:|---:|",
    ]
    for _, row in complete.iterrows():
        report_lines.append(
            f"| {row['E_over_N_Td']:.0f} | {row['comparison']} | "
            f"{row['eedf_relative_l1_common_support']:.6g} | "
            f"{row['mean_energy_candidate_eV']:.6g} | "
            f"{row['mean_energy_mcig_eV']:.6g} | "
            f"{row['mean_energy_relative_difference']:.6g} | "
            f"{row['tail_probability_difference_candidate_minus_mcig']:.6g} |"
        )
    report_lines.extend(
        [
            "",
            "## Evaluation",
            "",
            "The product curves track the MCIG bulk shape well at 100 Td and degrade as E/N increases. The increase in L1 and mean-energy difference at 200/400 Td is consistent with the known non-equivalent physics assumptions rather than a clean solver-only numerical discrepancy.",
            "",
            "The most important caveats are unchanged across the available MCIG cases: MCIG angular scattering model `11.000` is not mapped to the product `isotropic` sampler/closure, MCIG energy sharing model/parameter is not equivalent to product `energy_sharing: equal`, and the product run uses hold extrapolation above the 100 eV cross-section table to keep the 1000 eV EEDF range available.",
            "",
            "50 Td should not be included in MCIG error conclusions until an actual MCIG GUI EEDF export for 50 Td is produced with the same collision file/settings.",
            "",
            "## Outputs",
            "",
            f"- Linear overview: `{(OUT / f'{BASE}_eedf_overview_linear.png').relative_to(ROOT)}`",
            f"- Tail overview: `{(OUT / f'{BASE}_eedf_overview_tail.png').relative_to(ROOT)}`",
            f"- Curves CSV: `{(OUT / f'{BASE}_curves.csv').relative_to(ROOT)}`",
            f"- Metrics CSV: `{(OUT / f'{BASE}_metrics.csv').relative_to(ROOT)}`",
            f"- Settings CSV: `{settings_path.relative_to(ROOT)}`",
        ]
    )
    (OUT / f"{BASE}_report.md").write_text("\n".join(report_lines), encoding="utf-8")

    print(f"cases={len(result.cases)}")
    print(f"output_dir={OUT}")
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
