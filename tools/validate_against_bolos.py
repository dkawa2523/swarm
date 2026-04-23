#!/usr/bin/env python3
"""Compare native_bolsig and optional BOLOS backends for one YAML input.

The script exits with code 77 when BOLOS is not installed, which lets CI mark the
comparison as skipped while still running native tests everywhere.
"""

from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from electron_swarm import load_config, run  # noqa: E402


def relerr(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1.0e-300)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", type=Path)
    parser.add_argument("--rtol", type=float, default=0.15, help="relative tolerance for scalar transport checks")
    args = parser.parse_args(argv)

    try:
        import bolos  # noqa: F401
    except Exception:
        print("SKIP: optional bolos package is not installed")
        return 77

    base = load_config(args.config)
    native = copy.deepcopy(base)
    native.boltzmann_two_term.backend = "native_bolsig"
    bolos_cfg = copy.deepcopy(base)
    bolos_cfg.boltzmann_two_term.backend = "bolos"

    native_res = run(native, write=False)
    bolos_res = run(bolos_cfg, write=False)
    rows = []
    failures = []
    fields = [
        "mean_energy_eV",
        "reduced_mobility_m2_V_s_m3",
        "reduced_diffusion_L_m2_s_m3",
        "drift_velocity_m_s",
    ]
    for n_case, b_case in zip(native_res.cases, bolos_res.cases, strict=True):
        for field in fields:
            a = float(getattr(n_case, field))
            b = float(getattr(b_case, field))
            err = relerr(a, b)
            rows.append({"case_id": n_case.case_id, "field": field, "native": a, "bolos": b, "relerr": err})
            if err > args.rtol:
                failures.append(rows[-1])
    frame = pd.DataFrame(rows)
    print(frame.to_string(index=False))
    if failures:
        print(f"FAIL: {len(failures)} comparisons exceeded rtol={args.rtol}")
        return 1
    print("PASS: native_bolsig agrees with BOLOS within configured tolerance")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
