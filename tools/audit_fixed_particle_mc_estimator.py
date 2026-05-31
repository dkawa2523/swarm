"""Development audit for fixed-particle internal MC estimator stability.

This tool intentionally lives outside the product runner. It varies seed and
collision count, then checks whether EEDF statistics, energy bookkeeping, and
null-collision bookkeeping are stable enough to trust a benchmark comparison.
"""

from __future__ import annotations

import argparse
import copy
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from electron_swarm import load_config, run


RUN_COLUMNS = [
    "E_over_N_Td",
    "seed",
    "warmup_collisions",
    "max_collisions",
    "particles",
    "mean_energy_eV",
    "drift_velocity_m_s",
    "net_ionization_frequency_s",
    "mc_energy_balance_status",
    "mc_tracked_energy_balance_residual_fraction",
    "mc_tail_comparison_status",
    "mc_tail_weak_probability_fraction",
    "mc_tail_effective_sample_count_min",
    "mc_max_resolved_energy_eV",
    "mc_null_collision_acceptance_fraction",
    "mc_max_collision_to_trial_ratio",
    "mc_samples",
    "mc_nonzero_bins",
]


def _parse_ints(text: str) -> list[int]:
    return [int(item.strip()) for item in text.split(",") if item.strip()]


def _parse_floats(text: str) -> list[float]:
    return [float(item.strip()) for item in text.split(",") if item.strip()]


def _load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle)
    if not isinstance(loaded, dict):
        raise ValueError(f"{path} is not a mapping config")
    return loaded


def _absolutize_reference_paths(data: dict[str, Any], base_dir: Path) -> None:
    cross_sections = data.get("cross_sections")
    if not isinstance(cross_sections, dict):
        return
    files = cross_sections.get("files", [])
    if not isinstance(files, list):
        return
    for item in files:
        if not isinstance(item, dict) or "path" not in item:
            continue
        path = Path(str(item["path"]))
        if not path.is_absolute():
            item["path"] = str((base_dir / path).resolve())


def _prepare_run_config(
    base: dict[str, Any],
    *,
    seed: int,
    max_collisions: int,
    warmup_collisions: int,
    particles: int | None,
    e_over_n: list[float],
    output_dir: Path,
) -> dict[str, Any]:
    data = copy.deepcopy(base)
    run_block = data.setdefault("run", {})
    run_block["solvers"] = [{"id": "monte_carlo"}]
    run_block["e_over_n_Td"] = e_over_n
    run_block["case_prefix"] = f"fixed_mc_audit_s{seed}_c{max_collisions}"

    solver = data.setdefault("solvers", {}).setdefault("monte_carlo", {})
    solver["backend"] = "internal"
    solver["angular_scattering"] = "same_as_physics"
    solver["population_model"] = "fixed_particle_single_daughter"
    solver["seed"] = seed
    solver["warmup_collisions"] = warmup_collisions
    solver["max_collisions"] = max_collisions
    if particles is not None:
        solver["particles"] = particles

    output = data.setdefault("output", {})
    output["directory"] = str(output_dir.resolve())
    output["base_name"] = f"fixed_mc_audit_s{seed}_c{max_collisions}"
    return data


def _write_config(data: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(data, handle, sort_keys=False)


def _run_case_grid(
    base: dict[str, Any],
    *,
    seeds: list[int],
    collision_counts: list[int],
    warmup_collisions: int,
    particles: int | None,
    e_over_n: list[float],
    output_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    run_rows: list[dict[str, Any]] = []
    eedf_rows: list[dict[str, Any]] = []
    config_dir = output_dir / "_configs"
    for seed in seeds:
        for collisions in collision_counts:
            config_path = config_dir / f"fixed_mc_audit_s{seed}_c{collisions}.yaml"
            run_config = _prepare_run_config(
                base,
                seed=seed,
                max_collisions=collisions,
                warmup_collisions=warmup_collisions,
                particles=particles,
                e_over_n=e_over_n,
                output_dir=output_dir,
            )
            _write_config(run_config, config_path)
            result = run(load_config(config_path), write=False)
            for case in result.cases:
                meta = case.metadata
                row = {
                    "E_over_N_Td": case.e_over_n_Td,
                    "seed": seed,
                    "warmup_collisions": meta.get("mc_warmup_collisions", warmup_collisions),
                    "max_collisions": collisions,
                    "particles": meta.get("particles", particles),
                    "mean_energy_eV": case.mean_energy_eV,
                    "drift_velocity_m_s": case.drift_velocity_m_s,
                    "net_ionization_frequency_s": case.net_ionization_frequency_s,
                    "mc_energy_balance_status": meta.get("mc_energy_balance_status", ""),
                    "mc_tracked_energy_balance_residual_fraction": meta.get(
                        "mc_tracked_energy_balance_residual_fraction", np.nan
                    ),
                    "mc_tail_comparison_status": meta.get(
                        "mc_tail_comparison_status", ""
                    ),
                    "mc_tail_weak_probability_fraction": meta.get(
                        "mc_tail_weak_probability_fraction", np.nan
                    ),
                    "mc_tail_effective_sample_count_min": meta.get(
                        "mc_tail_effective_sample_count_min", np.nan
                    ),
                    "mc_max_resolved_energy_eV": meta.get(
                        "mc_max_resolved_energy_eV", np.nan
                    ),
                    "mc_null_collision_acceptance_fraction": meta.get(
                        "mc_null_collision_acceptance_fraction", np.nan
                    ),
                    "mc_max_collision_to_trial_ratio": meta.get(
                        "mc_max_collision_to_trial_ratio", np.nan
                    ),
                    "mc_samples": meta.get("mc_samples", np.nan),
                    "mc_nonzero_bins": meta.get("mc_nonzero_bins", np.nan),
                }
                run_rows.append(row)
                effective = case.eedf_effective_counts
                for index, (energy, eedf) in enumerate(
                    zip(case.energy_eV, case.eedf, strict=False)
                ):
                    eedf_rows.append(
                        {
                            "E_over_N_Td": case.e_over_n_Td,
                            "seed": seed,
                            "max_collisions": collisions,
                            "energy_eV": float(energy),
                            "eedf_eV_inv": float(eedf),
                            "effective_sample_count": (
                                float(effective[index])
                                if effective is not None and index < len(effective)
                                else np.nan
                            ),
                        }
                    )
    return pd.DataFrame(run_rows), pd.DataFrame(eedf_rows)


def _relative_l1(reference: pd.DataFrame, candidate: pd.DataFrame) -> float:
    ref = reference.sort_values("energy_eV")
    cand = candidate.sort_values("energy_eV")
    energy = ref["energy_eV"].to_numpy(float)
    ref_f = ref["eedf_eV_inv"].to_numpy(float)
    cand_f = np.interp(
        energy,
        cand["energy_eV"].to_numpy(float),
        cand["eedf_eV_inv"].to_numpy(float),
        left=np.nan,
        right=np.nan,
    )
    mask = np.isfinite(ref_f) & np.isfinite(cand_f)
    if int(np.count_nonzero(mask)) < 3:
        return float("nan")
    denom = np.trapezoid(np.abs(ref_f[mask]), energy[mask])
    if denom <= 0.0:
        return float("nan")
    return float(np.trapezoid(np.abs(ref_f[mask] - cand_f[mask]), energy[mask]) / denom)


def _convergence_rows(run_rows: pd.DataFrame, eedf_rows: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    if run_rows.empty:
        return pd.DataFrame(rows)
    max_collision = int(run_rows["max_collisions"].max())
    for (field, seed), group in run_rows.groupby(["E_over_N_Td", "seed"]):
        ref_group = group[group["max_collisions"] == max_collision]
        if ref_group.empty:
            continue
        ref_run = ref_group.iloc[0]
        ref_curve = eedf_rows[
            (eedf_rows["E_over_N_Td"] == field)
            & (eedf_rows["seed"] == seed)
            & (eedf_rows["max_collisions"] == max_collision)
        ]
        for _, row in group.iterrows():
            cand_curve = eedf_rows[
                (eedf_rows["E_over_N_Td"] == field)
                & (eedf_rows["seed"] == seed)
                & (eedf_rows["max_collisions"] == row["max_collisions"])
            ]
            rows.append(
                {
                    "E_over_N_Td": field,
                    "seed": seed,
                    "max_collisions": int(row["max_collisions"]),
                    "reference_max_collisions": max_collision,
                    "mean_energy_relative_delta_to_longest": abs(
                        float(row["mean_energy_eV"]) - float(ref_run["mean_energy_eV"])
                    )
                    / max(abs(float(ref_run["mean_energy_eV"])), 1.0e-300),
                    "eedf_relative_l1_to_longest": _relative_l1(ref_curve, cand_curve),
                }
            )
    return pd.DataFrame(rows)


def _seed_stats(run_rows: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for (field, collisions), group in run_rows.groupby(["E_over_N_Td", "max_collisions"]):
        mean_values = group["mean_energy_eV"].astype(float)
        rows.append(
            {
                "E_over_N_Td": field,
                "max_collisions": collisions,
                "seed_count": len(group),
                "mean_energy_mean_eV": float(mean_values.mean()),
                "mean_energy_std_eV": float(mean_values.std(ddof=1))
                if len(group) > 1
                else 0.0,
                "mean_energy_relative_seed_std": (
                    float(mean_values.std(ddof=1) / mean_values.mean())
                    if len(group) > 1 and mean_values.mean() != 0.0
                    else 0.0
                ),
                "max_energy_balance_residual_fraction": float(
                    group["mc_tracked_energy_balance_residual_fraction"]
                    .astype(float)
                    .max()
                ),
                "max_collision_to_trial_ratio": float(
                    group["mc_max_collision_to_trial_ratio"].astype(float).max()
                ),
                "tail_statuses": ",".join(
                    sorted(set(str(item) for item in group["mc_tail_comparison_status"]))
                ),
            }
        )
    return pd.DataFrame(rows)


def _write_report(
    output_dir: Path,
    run_rows: pd.DataFrame,
    seed_stats: pd.DataFrame,
    convergence: pd.DataFrame,
) -> None:
    max_balance = (
        float(run_rows["mc_tracked_energy_balance_residual_fraction"].astype(float).max())
        if not run_rows.empty
        else float("nan")
    )
    max_trial_ratio = (
        float(run_rows["mc_max_collision_to_trial_ratio"].astype(float).max())
        if not run_rows.empty
        else float("nan")
    )
    tail_statuses = (
        sorted(set(str(item) for item in run_rows["mc_tail_comparison_status"]))
        if not run_rows.empty
        else []
    )
    report = output_dir / "fixed_particle_mc_estimator_audit_report.md"
    report.write_text(
        "# Fixed-particle MC estimator/warmup audit\n\n"
        "Development-only audit. This tool does not add product schema or solver "
        "options; it reruns the existing fixed-particle internal MC across seeds "
        "and collision counts.\n\n"
        "## Headline\n\n"
        f"- Tail statuses: `{', '.join(tail_statuses)}`\n"
        f"- Max tracked energy-balance residual fraction: `{max_balance:.6g}`\n"
        f"- Max null-collision collision/trial ratio: `{max_trial_ratio:.6g}`\n\n"
        "## Files\n\n"
        "- `fixed_particle_mc_estimator_runs.csv`\n"
        "- `fixed_particle_mc_estimator_seed_stats.csv`\n"
        "- `fixed_particle_mc_estimator_convergence.csv`\n"
        "- `fixed_particle_mc_estimator_eedf_curves.csv`\n\n"
        "## Interpretation\n\n"
        "- Large seed standard deviation points to ordinary stochastic uncertainty.\n"
        "- Large delta to the longest collision count points to warmup/run-length "
        "sensitivity because the current estimator samples from the first flight.\n"
        "- Large energy-balance residual or collision/trial ratio above 1 points to "
        "bookkeeping or null-collision majorant problems.\n",
        encoding="utf-8",
    )


def _plot(output_dir: Path, seed_stats: pd.DataFrame, convergence: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt

    if seed_stats.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for field, group in seed_stats.groupby("E_over_N_Td"):
        group = group.sort_values("max_collisions")
        axes[0].errorbar(
            group["max_collisions"],
            group["mean_energy_mean_eV"],
            yerr=group["mean_energy_std_eV"],
            marker="o",
            label=f"{field:g} Td",
        )
    axes[0].set_xlabel("max_collisions")
    axes[0].set_ylabel("mean energy (eV)")
    axes[0].set_title("Seed ensemble mean energy")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend()

    if not convergence.empty:
        for field, group in convergence.groupby("E_over_N_Td"):
            summary = (
                group.groupby("max_collisions")["eedf_relative_l1_to_longest"]
                .mean()
                .reset_index()
                .sort_values("max_collisions")
            )
            axes[1].plot(
                summary["max_collisions"],
                summary["eedf_relative_l1_to_longest"],
                marker="o",
                label=f"{field:g} Td",
            )
    axes[1].set_xlabel("max_collisions")
    axes[1].set_ylabel("EEDF L1 to longest run")
    axes[1].set_title("Run-length convergence proxy")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(output_dir / "fixed_particle_mc_estimator_audit.png", dpi=180)
    fig.savefig(output_dir / "fixed_particle_mc_estimator_audit.svg")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/benchmarks/fixed_particle_mc_estimator_audit"),
    )
    parser.add_argument("--e-over-n", default="")
    parser.add_argument("--seeds", default="20260601,20260602,20260603")
    parser.add_argument("--collision-counts", default="300,600,1200")
    parser.add_argument("--warmup-collisions", type=int, default=0)
    parser.add_argument("--particles", type=int, default=None)
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()

    base = _load_yaml(args.config)
    _absolutize_reference_paths(base, args.config.resolve().parent)
    e_over_n = (
        _parse_floats(args.e_over_n)
        if args.e_over_n.strip()
        else [float(item) for item in base["run"]["e_over_n_Td"]]
    )
    seeds = _parse_ints(args.seeds)
    collision_counts = _parse_ints(args.collision_counts)
    if not seeds or not collision_counts or not e_over_n:
        raise ValueError("seeds, collision-counts, and e-over-n must be nonempty")
    if args.warmup_collisions < 0:
        raise ValueError("warmup-collisions must be nonnegative")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    runs, eedfs = _run_case_grid(
        base,
        seeds=seeds,
        collision_counts=collision_counts,
        warmup_collisions=args.warmup_collisions,
        particles=args.particles,
        e_over_n=e_over_n,
        output_dir=args.output_dir,
    )
    seed_stats = _seed_stats(runs)
    convergence = _convergence_rows(runs, eedfs)
    runs.to_csv(args.output_dir / "fixed_particle_mc_estimator_runs.csv", index=False)
    eedfs.to_csv(args.output_dir / "fixed_particle_mc_estimator_eedf_curves.csv", index=False)
    seed_stats.to_csv(
        args.output_dir / "fixed_particle_mc_estimator_seed_stats.csv", index=False
    )
    convergence.to_csv(
        args.output_dir / "fixed_particle_mc_estimator_convergence.csv", index=False
    )
    _write_report(args.output_dir, runs, seed_stats, convergence)
    if args.plot:
        _plot(args.output_dir, seed_stats, convergence)
    print(f"runs_csv: {args.output_dir / 'fixed_particle_mc_estimator_runs.csv'}")
    print(f"seed_stats_csv: {args.output_dir / 'fixed_particle_mc_estimator_seed_stats.csv'}")
    print(f"convergence_csv: {args.output_dir / 'fixed_particle_mc_estimator_convergence.csv'}")
    print(f"eedf_curves_csv: {args.output_dir / 'fixed_particle_mc_estimator_eedf_curves.csv'}")
    print(f"report: {args.output_dir / 'fixed_particle_mc_estimator_audit_report.md'}")


if __name__ == "__main__":
    main()
