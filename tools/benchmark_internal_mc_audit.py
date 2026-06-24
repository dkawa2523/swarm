"""Write compact audit tables for the product internal Monte Carlo backend."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from electron_swarm import load_config, run
from electron_swarm.core.numerics import widths_from_centers
from tools.benchmark_common import case_audit_metadata, metadata_value


AUDIT_METADATA = [
    "mc_seed",
    "mc_particles",
    "mc_population_model",
    "mc_warmup_collisions",
    "mc_production_collisions",
    "mc_max_sampled_energy_eV",
    "mc_max_cross_section_energy_eV",
    "mc_energy_samples_above_xs_max",
    "mc_energy_samples_above_xs_max_fraction",
    "mc_histogram_samples",
    "mc_tail_histogram_samples",
    "mc_tail_histogram_sample_fraction",
    "transport_definition",
    "swarm_population_treatment",
    "eedf_sampling_clock",
    "nonconservative_growth_treatment",
    "mc_collision_count_null",
    "mc_collision_count_elastic",
    "mc_collision_count_excitation",
    "mc_collision_count_ionization",
    "mc_collision_count_attachment",
    "mc_collision_count_superelastic",
    "mc_secondary_electron_count",
    "mc_branching_resample_count",
    "mc_population_total_weight_final",
    "mc_population_log_growth_estimate_s_inv",
    "mc_population_weight_cv",
    "ionization_source_treatment",
    "ionization_branching_model",
    "secondary_electron_tracking",
    "inelastic_angular_model",
    "mc_energy_balance_status",
    "mc_tracked_energy_balance_residual_fraction",
    "mc_physical_branching_gap_eV",
    "mc_population_resampling_energy_adjustment_eV",
    "mc_null_collision_acceptance_fraction",
    "mc_max_collision_to_trial_ratio",
    "mc_tail_uncertainty_status",
    "mc_min_tail_bin_count",
    "mc_tail_effective_sample_count_min",
    "mc_tail_weak_probability_fraction",
    "mc_max_resolved_energy_eV",
    "mc_tail_comparison_status",
]


def _audit_rows(result) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case in result.cases:
        metadata = case_audit_metadata(case)
        if (
            case.solver != "monte_carlo"
            or metadata_value(metadata, "monte_carlo_backend") != "internal"
        ):
            continue
        row = {
            "case_id": case.case_id,
            "E_over_N_Td": case.e_over_n_Td,
            "mean_energy_eV": case.mean_energy_eV,
            "drift_velocity_m_s": case.drift_velocity_m_s,
            "net_ionization_frequency_s": case.net_ionization_frequency_s,
        }
        row.update({key: metadata_value(metadata, key) for key in AUDIT_METADATA})
        rows.append(row)
    return rows


def _tail_rows(result) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case in result.cases:
        if case.solver != "monte_carlo" or case.eedf_counts is None:
            continue
        effective = case.eedf_effective_counts
        widths = case.energy_widths_eV
        if widths is None or len(widths) != len(case.energy_eV):
            widths = widths_from_centers(case.energy_eV)
        widths = np.asarray(widths, dtype=float)
        energy_values = np.asarray(case.energy_eV, dtype=float)
        eedf_values = np.asarray(case.eedf, dtype=float)
        probabilities = np.clip(eedf_values, 0.0, None) * widths
        total_probability = float(np.sum(probabilities))
        cumulative = 0.0
        for index, (energy, eedf, count) in enumerate(zip(
            case.energy_eV, case.eedf, case.eedf_counts, strict=False
        )):
            count_int = int(count)
            effective_count = (
                float(effective[index])
                if effective is not None and index < len(effective)
                else float(count_int)
            )
            width = (
                float(widths[index])
                if index < len(widths) and np.isfinite(widths[index]) and widths[index] > 0.0
                else float("nan")
            )
            probability = (
                float(probabilities[index])
                if index < len(probabilities) and np.isfinite(probabilities[index])
                else float("nan")
            )
            survival = (
                max(total_probability - cumulative, 0.0)
                if np.isfinite(probability)
                else float("nan")
            )
            if np.isfinite(probability):
                cumulative += probability
            rows.append(
                {
                    "case_id": case.case_id,
                    "E_over_N_Td": case.e_over_n_Td,
                    "energy_eV": float(energy),
                    "energy_left_eV": float(energy_values[index] - 0.5 * width),
                    "energy_right_eV": float(energy_values[index] + 0.5 * width),
                    "energy_width_eV": width,
                    "eedf": float(eedf),
                    "probability_mass": probability,
                    "cumulative_probability": (
                        float(cumulative) if np.isfinite(probability) else ""
                    ),
                    "survival_probability": survival,
                    "sample_count": count_int,
                    "effective_sample_count": effective_count,
                    "relative_standard_error": (
                        float(1.0 / effective_count**0.5)
                        if effective_count > 0.0
                        else ""
                    ),
                }
            )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/benchmarks/internal_mc_audit"))
    args = parser.parse_args()

    cfg = load_config(args.config)
    result = run(cfg, write=False, collect_diagnostics=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    audit = pd.DataFrame(_audit_rows(result))
    tail = pd.DataFrame(_tail_rows(result))
    audit_path = args.output_dir / "internal_mc_audit_summary.csv"
    tail_path = args.output_dir / "internal_mc_tail_uncertainty.csv"
    report_path = args.output_dir / "internal_mc_audit_report.md"
    audit.to_csv(audit_path, index=False)
    tail.to_csv(tail_path, index=False)

    if audit.empty:
        body = "No internal monte_carlo cases were runnable in this config.\n"
    else:
        body = audit.to_csv(index=False)
    report_path.write_text(
        "# Internal MC Audit\n\n"
        "This benchmark audits the product internal MC model; it does not tune "
        "or fit MCIG.\n\n"
        f"- Audit CSV: `{audit_path.resolve()}`\n"
        f"- Tail uncertainty CSV: `{tail_path.resolve()}`\n\n"
        "```csv\n"
        f"{body}"
        "```\n",
        encoding="utf-8",
    )
    print(f"audit_csv: {audit_path}")
    print(f"tail_uncertainty_csv: {tail_path}")
    print(f"report: {report_path}")


if __name__ == "__main__":
    main()
