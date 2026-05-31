"""Write compact audit tables for the product internal Monte Carlo backend."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from electron_swarm import load_config, run


AUDIT_METADATA = [
    "mc_population_model",
    "mc_warmup_collisions",
    "mc_production_collisions",
    "ionization_source_treatment",
    "ionization_branching_model",
    "secondary_electron_tracking",
    "inelastic_angular_model",
    "mc_energy_balance_status",
    "mc_tracked_energy_balance_residual_fraction",
    "mc_physical_branching_gap_eV",
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
        if case.solver != "monte_carlo" or case.metadata.get("monte_carlo_backend") != "internal":
            continue
        row = {
            "case_id": case.case_id,
            "E_over_N_Td": case.e_over_n_Td,
            "mean_energy_eV": case.mean_energy_eV,
            "drift_velocity_m_s": case.drift_velocity_m_s,
            "net_ionization_frequency_s": case.net_ionization_frequency_s,
        }
        row.update({key: case.metadata.get(key, "") for key in AUDIT_METADATA})
        rows.append(row)
    return rows


def _tail_rows(result) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case in result.cases:
        if case.solver != "monte_carlo" or case.eedf_counts is None:
            continue
        effective = case.eedf_effective_counts
        for index, (energy, eedf, count) in enumerate(zip(
            case.energy_eV, case.eedf, case.eedf_counts, strict=False
        )):
            count_int = int(count)
            effective_count = (
                float(effective[index])
                if effective is not None and index < len(effective)
                else float(count_int)
            )
            rows.append(
                {
                    "case_id": case.case_id,
                    "E_over_N_Td": case.e_over_n_Td,
                    "energy_eV": float(energy),
                    "eedf": float(eedf),
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
    result = run(cfg, write=False)
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
