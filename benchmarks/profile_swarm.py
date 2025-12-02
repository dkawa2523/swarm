"""
Profile helper for swarm_mc runs.

Examples:
  python3 benchmarks/profile_swarm.py --num-e 15000 --num-col-max 4000 --seed 42
  python3 benchmarks/profile_swarm.py --matrix --config configs/Ar_N2_sweep.yaml
"""

from __future__ import annotations

import argparse
import cProfile
import os
from copy import deepcopy
from pathlib import Path
import time
import pstats

from swarm_mc import Simulation
from swarm_mc.caching import GasMixtureCache
from swarm_mc.runner import SimulationRunner


def _profile_single(config_path: str, num_e: int, num_col_max: int, seed: int,
                    label: str, output_root: Path,
                    gas_cache: GasMixtureCache) -> dict:
    runner = SimulationRunner.from_yaml(config_path)
    base_cfg = deepcopy(runner.experiment.base_config)
    base_cfg["initial_state"]["num_e_initial"] = num_e
    base_cfg["end_conditions"]["num_col_max"] = num_col_max
    base_cfg["simulation_settings"]["seed"] = seed
    base_cfg["output"]["base_name"] = label
    base_cfg["output"]["output_directory"] = str(output_root / label) + os.sep
    base_cfg["output"]["save_simulation_pickle"] = False

    sim = Simulation(base_cfg, gas_mixture_cache=gas_cache, run_label=label)
    profile_path = output_root / f"profile_{label}.prof"
    text_path = output_root / f"profile_{label}.txt"

    profile = cProfile.Profile()
    start = time.perf_counter()
    profile.runcall(sim.run)
    elapsed = time.perf_counter() - start
    profile.dump_stats(profile_path)
    with open(text_path, "w") as fh:
        stats = pstats.Stats(profile, stream=fh).sort_stats("cumulative")
        stats.print_stats(30)

    return {
        "label": label,
        "num_e": num_e,
        "num_col_max": num_col_max,
        "elapsed_time_s": elapsed,
        "profile_path": str(profile_path),
        "text_path": str(text_path),
    }


def _matrix_cases(base_label: str):
    return [
        (f"{base_label}_s", 8000, 2000, 42),
        (f"{base_label}_m", 15000, 4000, 42),
        (f"{base_label}_l", 25000, 5000, 42),
    ]


def main():
    parser = argparse.ArgumentParser(description="Profile a single swarm_mc run.")
    parser.add_argument("--config", default="configs/Ar_N2_sweep.yaml",
                        help="Path to YAML experiment config.")
    parser.add_argument("--num-e", type=int, default=15000,
                        help="Initial electrons for profiling.")
    parser.add_argument("--num-col-max", type=int, default=4000,
                        help="Collision cap for profiling.")
    parser.add_argument("--seed", type=int, default=42,
                        help="RNG seed for reproducibility.")
    parser.add_argument("--label", default=None,
                        help="Label for profile artifacts.")
    parser.add_argument("--matrix", action="store_true",
                        help="Run S/M/L matrix instead of a single case.")
    args = parser.parse_args()

    output_root = Path("outputs/profiles")
    output_root.mkdir(parents=True, exist_ok=True)
    gas_cache = GasMixtureCache()

    if args.matrix:
        base_label = args.label or Path(args.config).stem
        results = []
        for label, num_e, num_col_max, seed in _matrix_cases(base_label):
            results.append(_profile_single(args.config, num_e, num_col_max, seed,
                                           label, output_root, gas_cache))
        summary_path = output_root / "profile_summary.csv"
        header = "label,num_e,num_col_max,elapsed_time_s,profile_path,text_path\n"
        with open(summary_path, "w") as fh:
            fh.write(header)
            for row in results:
                fh.write(",".join([
                    row["label"],
                    str(row["num_e"]),
                    str(row["num_col_max"]),
                    f"{row['elapsed_time_s']:.4f}",
                    row["profile_path"],
                    row["text_path"],
                ]) + "\n")
        print(f"Wrote matrix summary to {summary_path}")
        return

    label = args.label or Path(args.config).stem
    result = _profile_single(args.config, args.num_e, args.num_col_max,
                             args.seed, label, output_root, gas_cache)
    print(f"Profiled '{label}' in {result['elapsed_time_s']:.3f} s")
    print(f"Stats: {result['text_path']}")
    print(f"cProfile: {result['profile_path']}")


if __name__ == "__main__":
    main()
