#!/usr/bin/env python3
"""Lightweight smoke sweep for the experimental multi-term operator backend.

This tool intentionally stays outside the solver package. It exercises a small
set of representative cases and reports lmax trends. Passing this tool means
the lmax>1 reference-anchored closure is numerically executable; it is not a
full multi-term physics validation. Use tools/benchmark_operator_gate.py for
two-term, BOLOS, and MC reference comparisons.
"""

from __future__ import annotations

import argparse
import copy
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from electron_swarm import load_config  # noqa: E402
from electron_swarm.core.config import SwarmConfig  # noqa: E402
from electron_swarm.core.cross_sections import (  # noqa: E402
    CrossSectionProcess,
    CrossSectionSet,
    ProcessType,
    load_cross_sections,
)
from electron_swarm.solvers.multiterm_boltzmann import MultiTermBoltzmannSolver  # noqa: E402


@dataclass(frozen=True, slots=True)
class Scenario:
    name: str
    config: SwarmConfig
    cross_sections: CrossSectionSet
    e_over_n_Td: float


def _relerr(a: float, b: float) -> float:
    return abs(a - b) / max(abs(a), abs(b), 1.0e-300)


def _base_operator_config(*, quick: bool) -> SwarmConfig:
    cfg = load_config(
        ROOT / "configs" / "unified" / "multiterm_operator_experimental.yaml"
    )
    cfg.output.write_plots = False
    cfg.multiterm_boltzmann.method = "operator"
    cfg.multiterm_boltzmann.allow_experimental_operator = True
    cfg.multiterm_boltzmann.hydrodynamic = False
    cfg.multiterm_boltzmann.energy_grid.n = 56 if quick else 80
    cfg.multiterm_boltzmann.energy_grid.max_eV = 60.0 if quick else 80.0
    return cfg


def _ar_n2_config(*, quick: bool) -> SwarmConfig:
    cfg = load_config(ROOT / "configs" / "unified" / "ar_n2_both.yaml")
    cfg.output.write_plots = False
    cfg.conditions.pressure_Pa = 100.0
    cfg.multiterm_boltzmann.method = "operator"
    cfg.multiterm_boltzmann.allow_experimental_operator = True
    cfg.multiterm_boltzmann.hydrodynamic = False
    cfg.multiterm_boltzmann.energy_grid.n = 56 if quick else 80
    cfg.multiterm_boltzmann.energy_grid.max_eV = 70.0 if quick else 100.0
    return cfg


def _synthetic_processes(kind: str) -> CrossSectionSet:
    energy = np.array([0.0, 3.0, 10.0, 40.0])
    momentum = CrossSectionProcess(
        species="Ar",
        process="momentum",
        process_type=ProcessType.MOMENTUM,
        threshold_eV=None,
        mass_amu=39.948,
        energy_eV=energy,
        cross_section_m2=np.array([1.0e-20, 1.0e-20, 1.2e-20, 1.2e-20]),
    )
    if kind == "attachment":
        extra = CrossSectionProcess(
            species="Ar",
            process="attachment",
            process_type=ProcessType.ATTACHMENT,
            threshold_eV=0.0,
            energy_eV=energy,
            cross_section_m2=np.array([0.0, 1.0e-22, 1.0e-22, 3.0e-23]),
        )
    elif kind == "superelastic":
        excitation = CrossSectionProcess(
            species="Ar",
            process="excitation",
            process_type=ProcessType.EXCITATION,
            threshold_eV=1.5,
            energy_eV=energy,
            cross_section_m2=np.array([0.0, 1.0e-22, 1.0e-22, 1.0e-22]),
        )
        superelastic = CrossSectionProcess(
            species="Ar",
            process="deexcitation",
            process_type=ProcessType.SUPERELASTIC,
            threshold_eV=1.5,
            energy_eV=energy,
            cross_section_m2=np.array([5.0e-26, 5.0e-26, 5.0e-26, 5.0e-26]),
        )
        return CrossSectionSet([momentum, excitation, superelastic])
    else:
        raise ValueError(kind)
    return CrossSectionSet([momentum, extra])


def _scenarios(*, quick: bool) -> list[Scenario]:
    ar = _base_operator_config(quick=quick)
    ar_xs = load_cross_sections(ar.cross_sections, ar.conditions)

    ar_n2 = _ar_n2_config(quick=quick)
    ar_n2_xs = load_cross_sections(ar_n2.cross_sections, ar_n2.conditions)

    attach_cfg = copy.deepcopy(ar)
    attach_cfg.multiterm_boltzmann.energy_grid.max_eV = 40.0
    super_cfg = copy.deepcopy(attach_cfg)

    return [
        Scenario("Ar", ar, ar_xs, 50.0),
        Scenario("Ar_N2", ar_n2, ar_n2_xs, 80.0),
        Scenario("Ar_high_EN_ionization", ar, ar_xs, 300.0),
        Scenario(
            "synthetic_attachment",
            attach_cfg,
            _synthetic_processes("attachment"),
            40.0,
        ),
        Scenario(
            "synthetic_superelastic",
            super_cfg,
            _synthetic_processes("superelastic"),
            40.0,
        ),
    ]


def _solve_operator(scenario: Scenario, *, lmax: int):
    cfg = copy.deepcopy(scenario.config)
    cfg.multiterm_boltzmann.method = "operator"
    cfg.multiterm_boltzmann.allow_experimental_operator = True
    cfg.multiterm_boltzmann.lmax = int(lmax)
    cfg.multiterm_boltzmann.hydrodynamic = False
    solver = MultiTermBoltzmannSolver(cfg, scenario.cross_sections)
    return solver.solve_case(
        scenario.e_over_n_Td,
        f"{scenario.name}_lmax{lmax}",
    )


def _operator_row(scenario: Scenario, lmax: int) -> tuple[dict, bool]:
    try:
        case = _solve_operator(scenario, lmax=lmax)
        norm = float(case.metadata.get("normalization_integral", np.nan))
        if not np.isfinite(norm) or abs(norm - 1.0) > 5.0e-6:
            raise RuntimeError(f"EEDF normalization failed ({norm:.6g})")
        fields = [
            case.mean_energy_eV,
            case.drift_velocity_m_s,
            case.net_ionization_frequency_s,
            float(case.metadata.get("operator_tail_rate_fraction", np.nan)),
            float(case.metadata.get("operator_highest_l_relative_l1", np.nan)),
        ]
        if not all(np.isfinite(v) for v in fields):
            raise RuntimeError("non-finite operator metric")
        row = {
            "scenario": scenario.name,
            "lmax": int(lmax),
            "status": "ok",
            "mean_energy_eV": case.mean_energy_eV,
            "drift_velocity_m_s": case.drift_velocity_m_s,
            "net_ionization_frequency_s": case.net_ionization_frequency_s,
            "attachment_rate_m3_s": case.metadata.get(
                "convolution_attachment_rate_coefficient_m3_s", np.nan
            ),
            "tail_rate_fraction": case.metadata.get(
                "operator_tail_rate_fraction", np.nan
            ),
            "highest_l_l1": case.metadata.get(
                "operator_highest_l_relative_l1", np.nan
            ),
            "diffusion_L_m2_s": case.diffusion_L_m2_s,
            "closure": case.metadata.get("physical_validity", ""),
        }
        return row, True
    except Exception as exc:
        return {
            "scenario": scenario.name,
            "lmax": int(lmax),
            "status": f"FAIL: {exc}",
        }, False


def _add_lmax_trends(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    for field in ["mean_energy_eV", "drift_velocity_m_s", "net_ionization_frequency_s"]:
        frame[f"{field}_rel_change_prev"] = np.nan
    for scenario, idxs in frame.groupby("scenario").groups.items():
        prev = None
        for idx in sorted(idxs, key=lambda i: frame.at[i, "lmax"]):
            if prev is not None and frame.at[idx, "status"] == "ok":
                for field in [
                    "mean_energy_eV",
                    "drift_velocity_m_s",
                    "net_ionization_frequency_s",
                ]:
                    frame.at[idx, f"{field}_rel_change_prev"] = _relerr(
                        float(frame.at[idx, field]), float(frame.at[prev, field])
                    )
            if frame.at[idx, "status"] == "ok":
                prev = idx
    return frame


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="use smaller grids for fast checks")
    args = parser.parse_args(argv)

    scenarios = _scenarios(quick=args.quick)
    rows = []
    ok = True
    for scenario in scenarios:
        for lmax in (2, 3, 4):
            row, passed = _operator_row(scenario, lmax)
            rows.append(row)
            ok = ok and passed

    frame = _add_lmax_trends(pd.DataFrame(rows))
    print("\nMULTI-TERM OPERATOR SMOKE SWEEP")
    print("NOTE: lmax>1 is reference-anchored; this is not a full multi-term validation.")
    print(frame.to_string(index=False))

    if not ok:
        print("\nFAIL: smoke checks found one or more hard failures")
        return 1
    print("\nPASS: smoke checks completed; benchmark validation is still required")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
