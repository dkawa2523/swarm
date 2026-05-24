"""Run lightweight schema v2 product benchmarks.

This script measures wall-clock time around product runner executions. It does
not write canonical solver outputs; only the compact benchmark summary CSV is
produced.
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from electron_swarm import load_config, run  # noqa: E402


def _configs_from_args(paths: list[Path]) -> list[Path]:
    if paths:
        return [path.resolve() for path in paths]
    return sorted((ROOT / "configs" / "benchmarks").glob("*.yaml"))


def _display_path(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def _run_one(config_path: Path) -> dict[str, object]:
    start = time.perf_counter()
    try:
        result = run(load_config(config_path), write=False)
        elapsed = time.perf_counter() - start
        solvers = sorted({case.solver for case in result.cases})
        mean_energies = [case.mean_energy_eV for case in result.cases]
        drifts = [abs(case.drift_velocity_m_s) for case in result.cases]
        return {
            "config": _display_path(config_path),
            "status": "ok",
            "cases": len(result.cases),
            "solvers": ";".join(solvers),
            "elapsed_s": f"{elapsed:.6g}",
            "max_mean_energy_eV": f"{max(mean_energies, default=float('nan')):.6g}",
            "max_abs_drift_velocity_m_s": f"{max(drifts, default=float('nan')):.6g}",
            "error": "",
        }
    except Exception as exc:  # benchmark summaries should capture failures too
        elapsed = time.perf_counter() - start
        return {
            "config": _display_path(config_path),
            "status": "error",
            "cases": 0,
            "solvers": "",
            "elapsed_s": f"{elapsed:.6g}",
            "max_mean_energy_eV": "",
            "max_abs_drift_velocity_m_s": "",
            "error": f"{type(exc).__name__}: {exc}",
        }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run schema v2 product benchmark configs."
    )
    parser.add_argument(
        "--config",
        action="append",
        type=Path,
        default=[],
        help="Benchmark config to run. May be provided multiple times.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "outputs" / "benchmarks" / "benchmark_summary.csv",
        help="CSV summary path.",
    )
    parser.add_argument(
        "--fail-fast",
        action="store_true",
        help="Return nonzero when any benchmark config fails.",
    )
    args = parser.parse_args(argv)

    configs = _configs_from_args(args.config)
    if not configs:
        raise SystemExit("no benchmark configs found")

    rows = [_run_one(path) for path in configs]
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as fp:
        writer = csv.DictWriter(
            fp,
            fieldnames=[
                "config",
                "status",
                "cases",
                "solvers",
                "elapsed_s",
                "max_mean_energy_eV",
                "max_abs_drift_velocity_m_s",
                "error",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"{row['status']:>5} {row['elapsed_s']}s {row['config']} "
            f"cases={row['cases']} solvers={row['solvers']}"
        )
        if row["error"]:
            print(f"      {row['error']}")
    print(f"wrote {output}")
    return 1 if args.fail_fast and any(row["status"] != "ok" for row in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
