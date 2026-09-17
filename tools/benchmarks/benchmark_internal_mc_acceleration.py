"""Bounded end-to-end benchmark for the internal Monte Carlo numeric kernels."""

from __future__ import annotations

import argparse
from copy import deepcopy
from datetime import datetime, timezone
import json
from pathlib import Path
import platform
from time import perf_counter

import numba
import numpy as np

from electron_swarm import load_config, run


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "examples" / "argon_gec_ccp_monte_carlo_weighted_branching.yaml"


def _run_case(config, kernel: str):  # type: ignore[no-untyped-def]
    candidate = deepcopy(config)
    candidate.solvers.monte_carlo.numeric_kernel = kernel
    started = perf_counter()
    case = run(candidate, write=False).cases[0]
    elapsed = perf_counter() - started
    return case, elapsed


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--fields", type=float, nargs="+", default=[30.0, 1000.0, 4000.0]
    )
    parser.add_argument("--particles", type=int, default=64)
    parser.add_argument("--warmup", type=int, default=16)
    parser.add_argument("--production", type=int, default=64)
    parser.add_argument("--tail", type=int, default=64)
    parser.add_argument("--lag", type=int, default=8)
    parser.add_argument("--seed", type=int, default=444)
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args()
    if args.repeats <= 0:
        parser.error("--repeats must be positive")

    base = load_config(args.config.resolve())
    mc = base.solvers.monte_carlo
    mc.seed = args.seed
    mc.particles = args.particles
    mc.warmup_collisions = args.warmup
    mc.max_collisions = args.production
    mc.tail_max_collisions = args.tail
    mc.transport_correlation_lag_barriers = args.lag

    # Compile outside the timed region.  Persistent worker processes pay this
    # once; the benchmark reports steady-state solver cost explicitly.
    warm = deepcopy(base)
    warm.run.e_over_n_Td = [float(args.fields[0])]
    warm.solvers.monte_carlo.particles = 2
    warm.solvers.monte_carlo.warmup_collisions = 1
    warm.solvers.monte_carlo.max_collisions = 1
    warm.solvers.monte_carlo.tail_max_collisions = 1
    warm.solvers.monte_carlo.transport_correlation_lag_barriers = 4
    warm.solvers.monte_carlo.numeric_kernel = "numba"
    run(warm, write=False)

    rows: list[dict[str, object]] = []
    for field in args.fields:
        case_config = deepcopy(base)
        case_config.run.e_over_n_Td = [float(field)]
        reference_samples: list[float] = []
        compiled_samples: list[float] = []
        for _ in range(args.repeats):
            reference, reference_s = _run_case(case_config, "python")
            compiled, compiled_s = _run_case(case_config, "numba")
            reference_samples.append(reference_s)
            compiled_samples.append(compiled_s)
        reference_s = float(np.median(reference_samples))
        compiled_s = float(np.median(compiled_samples))
        reference_rates = np.asarray(
            [item.rate_coefficient_m3_s for item in reference.rates], dtype=float
        )
        compiled_rates = np.asarray(
            [item.rate_coefficient_m3_s for item in compiled.rates], dtype=float
        )
        reference_provenance = reference.diagnostics["internal_monte_carlo_transport"][
            "mc_run_provenance"
        ]
        compiled_provenance = compiled.diagnostics["internal_monte_carlo_transport"][
            "mc_run_provenance"
        ]
        rows.append(
            {
                "e_over_n_Td": float(field),
                "reference_seconds": reference_s,
                "compiled_seconds": compiled_s,
                "reference_seconds_samples": reference_samples,
                "compiled_seconds_samples": compiled_samples,
                "steady_state_speedup": reference_s / compiled_s,
                "reference_kernel_used": reference_provenance["numeric_kernel_used"],
                "compiled_kernel_used": compiled_provenance["numeric_kernel_used"],
                "reference_tail_collisions_executed": reference_provenance[
                    "tail_collisions_executed"
                ],
                "compiled_tail_collisions_executed": compiled_provenance[
                    "tail_collisions_executed"
                ],
                "eedf_counts_equal": bool(
                    np.array_equal(reference.eedf_counts, compiled.eedf_counts)
                ),
                "eedf_equal": bool(np.array_equal(reference.eedf, compiled.eedf)),
                "mean_energy_relative_difference": abs(
                    compiled.mean_energy_eV - reference.mean_energy_eV
                )
                / max(abs(reference.mean_energy_eV), 1.0e-300),
                "mobility_relative_difference": abs(
                    compiled.transport.mobility_m2_V_s
                    - reference.transport.mobility_m2_V_s
                )
                / max(abs(reference.transport.mobility_m2_V_s), 1.0e-300),
                "rate_peak_normalized_max_difference": float(
                    np.max(np.abs(compiled_rates - reference_rates))
                    / max(float(np.max(np.abs(reference_rates))), 1.0e-300)
                ),
            }
        )

    payload = {
        "evaluated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "bounded end-to-end run(config, write=False), compilation warmed",
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "numba": numba.__version__,
        },
        "settings": {
            "particles": args.particles,
            "warmup_collisions": args.warmup,
            "max_collisions": args.production,
            "tail_max_collisions": args.tail,
            "transport_correlation_lag_barriers": args.lag,
            "seed": args.seed,
            "repeats": args.repeats,
        },
        "cases": rows,
        "geometric_mean_speedup": float(
            np.exp(
                np.mean(np.log([float(row["steady_state_speedup"]) for row in rows]))
            )
        ),
        "production_scale_speedup_not_claimed": True,
    }
    encoded = json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    if args.output is not None:
        output = args.output.resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(encoded, encoding="utf-8")
    print(encoded, end="")


if __name__ == "__main__":
    main()
