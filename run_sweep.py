"""
Standalone runner for swarm calculations using the improved swarm_mc modules.

Usage:
    python3 run_sweep.py --config configs/Ar_N2_sweep.yaml
"""

import argparse
import sys
from swarm_mc import SimulationRunner
from swarm_mc.visualization import (
    plot_summary,
    plot_summary_metrics,
    plot_time_series,
    plot_energy_distribution,
    plot_energy_distribution_combined,
)
from pathlib import Path


def _resolve_comsol_export() -> object:
    try:
        from swarm_comsol_export.hook import maybe_export_comsol
        return maybe_export_comsol
    except ModuleNotFoundError:
        repo_root = Path(__file__).resolve().parent
        candidate = repo_root / "swarm_comsol_exporter"
        if (candidate / "swarm_comsol_export").exists():
            sys.path.insert(0, str(candidate))
            from swarm_comsol_export.hook import maybe_export_comsol
            return maybe_export_comsol
    return None


def main():
    parser = argparse.ArgumentParser(description="Run swarm simulations from a YAML config")
    parser.add_argument("--config", "-c", default="configs/Ar_N2_sweep.yaml",
                        help="Path to YAML experiment config")
    parser.add_argument("--plot", action="store_true",
                        help="Plot summary and first-run outputs after execution if available")
    parser.add_argument("--save-plots", action="store_true",
                        help="Save summary/time-series/EEDF plots to files (png)")
    parser.add_argument("--profile-per-step", action="store_true",
                        help="Record per-step timing breakdown to CSV (opt-in)")
    args = parser.parse_args()

    runner = SimulationRunner.from_yaml(args.config)
    if args.profile_per_step:
        sim_cfg = runner.experiment.base_config.setdefault("simulation_settings", {})
        sim_cfg["profile_per_step"] = True

    runs = runner.experiment.build_runs()
    sims = runner.run_all()

    maybe_export_comsol = _resolve_comsol_export()
    if maybe_export_comsol is not None:
        run_dir = Path(runner.experiment.output_root) / runner.experiment.name
        out_dir = maybe_export_comsol(args.config, run_dir=run_dir)
        if out_dir:
            print(f"[swarm-comsol-export] Output saved to: {out_dir}")

    if runner.run_stats:
        print("\nRun timing (s):")
        for stat in runner.run_stats:
            en = stat.get("E/N (Td)")
            en_text = f" (E/N={en} Td)" if en is not None else ""
            print(f"  {stat['run_label']}{en_text}: {stat['run_time (s)']:.3f} s")
        if runner.total_runtime_s is not None:
            print(f"Total sweep time: {runner.total_runtime_s:.3f} s")
        if runner.log_path:
            print(f"Execution log saved to: {runner.log_path}")

    plotting_cfg = runner.experiment.plotting or {}
    do_plot_summary = args.plot or plotting_cfg.get("summary", False)
    do_plot_first = args.plot or plotting_cfg.get("first_run", False)
    do_save_plots = args.save_plots or plotting_cfg.get("save", False)
    save_dir = Path(plotting_cfg.get("save_dir",
                                     Path(runner.experiment.output_root) /
                                     runner.experiment.name / "plots"))
    if do_save_plots:
        save_dir.mkdir(parents=True, exist_ok=True)

    if do_plot_summary:
        summary_path = Path(runner.experiment.output_root) / runner.experiment.name / "summary.csv"
        if summary_path.exists():
            plot_summary(str(summary_path),
                         show=False,
                         save_path=(save_dir / "summary.png" if do_save_plots else None))
            metrics = [
                ("flux drift velocity (m.s-1)", "flux_drift_velocity.png", "flux drift velocity (m/s)"),
                ("bulk drift velocity (m.s-1)", "bulk_drift_velocity.png", "bulk drift velocity (m/s)"),
                ("flux L diffusion coeff. * N (m-1.s-1)", "flux_diffusion_long.png", "flux longitudinal diffusion * N"),
                ("flux T diffusion coeff. * N (m-1.s-1)", "flux_diffusion_trans.png", "flux transverse diffusion * N"),
                ("bulk L diffusion coeff. * N (m-1.s-1)", "bulk_diffusion_long.png", "bulk longitudinal diffusion * N"),
                ("bulk T diffusion coeff. * N (m-1.s-1)", "bulk_diffusion_trans.png", "bulk transverse diffusion * N"),
            ]
            if do_save_plots:
                plot_summary_metrics(
                    str(summary_path), metrics, show=False, save_dir=save_dir
                )
    if do_plot_first and runs:
        # plot first run outputs if present
        first_run_dir = Path(runs[0].config["output"]["output_directory"])
        base = runs[0].config["output"]["base_name"]
        te = first_run_dir / f"{base}_temporal_evolution.csv"
        ed = first_run_dir / f"{base}_energy_distribution.csv"
        if te.exists():
            plot_time_series(str(te),
                             show=False,
                             save_path=(save_dir / f"{base}_temporal.png" if do_save_plots else None))
        if ed.exists():
            plot_energy_distribution(str(ed),
                                     show=False,
                                     save_path=(save_dir / f"{base}_energy.png" if do_save_plots else None))
            plot_energy_distribution_combined(
                str(ed),
                show=False,
                save_path=(save_dir / f"{base}_energy_combined.png" if do_save_plots else None),
            )

    if do_save_plots:
        for planned in runs:
            run_dir = Path(planned.config["output"]["output_directory"])
            base = planned.config["output"]["base_name"]
            te = run_dir / f"{base}_temporal_evolution.csv"
            ed = run_dir / f"{base}_energy_distribution.csv"
            if te.exists():
                plot_time_series(str(te), show=False,
                                 save_path=run_dir / f"{base}_temporal.png")
            if ed.exists():
                plot_energy_distribution(str(ed), show=False,
                                         save_path=run_dir / f"{base}_energy.png")
                plot_energy_distribution_combined(
                    str(ed), show=False,
                    save_path=run_dir / f"{base}_energy_combined.png",
                )


if __name__ == "__main__":
    main()
