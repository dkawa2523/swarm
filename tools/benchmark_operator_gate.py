#!/usr/bin/env python3
"""Small benchmark gate for the experimental multi-term operator.

The gate keeps three reference layers separate:

* native two-term is the always-available reference;
* BOLOS is optional unless ``--require-bolos`` is set;
* MC is opt-in because it is stochastic and slower.

The lmax>1 comparison is report-only because that path is a reference-anchored
integral-cross-section closure. Hard failures are reserved for non-finite
values, normalization failures, and reference paths that should be deterministic.
"""

from __future__ import annotations

import argparse
import contextlib
import copy
import io
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from electron_swarm import load_config  # noqa: E402
from electron_swarm.core.config import SwarmConfig  # noqa: E402
from electron_swarm.core.cross_sections import CrossSectionSet, load_cross_sections  # noqa: E402
from electron_swarm.solvers.boltzmann_two_term import BoltzmannTwoTermSolver  # noqa: E402
from electron_swarm.solvers.monte_carlo_adapter import MonteCarloAdapter  # noqa: E402
from electron_swarm.solvers.multiterm_boltzmann import MultiTermBoltzmannSolver  # noqa: E402


REPORT_COLUMNS = (
    "scenario",
    "candidate",
    "reference",
    "metric",
    "candidate_value",
    "reference_value",
    "relerr",
    "tolerance",
    "status",
    "note",
)

LMAX1_TOLERANCE = 1.0e-10
BOLOS_TOLERANCE = 0.20
MC_BASE_TOLERANCES = {
    "mean_energy_eV": 0.40,
    "drift_velocity_m_s": 0.40,
    "reduced_mobility_m2_V_s_m3": 0.40,
    "reduced_diffusion_L_m2_s_m3": 0.60,
    "net_ionization_frequency_s": 0.75,
}
REFERENCE_FIELDS = (
    "mean_energy_eV",
    "drift_velocity_m_s",
    "reduced_mobility_m2_V_s_m3",
    "reduced_diffusion_L_m2_s_m3",
    "net_ionization_frequency_s",
)


@dataclass(frozen=True, slots=True)
class ScenarioSpec:
    name: str
    config_path: Path
    e_over_n_Td: float
    allow_bolos: bool = True
    allow_mc: bool = False


@dataclass(frozen=True, slots=True)
class LoadedScenario:
    name: str
    config: SwarmConfig
    cross_sections: CrossSectionSet
    e_over_n_Td: float
    allow_bolos: bool = True
    allow_mc: bool = False


def _scenario_specs() -> tuple[ScenarioSpec, ...]:
    return (
        ScenarioSpec("Ar_50Td", ROOT / "configs" / "unified" / "both_template.yaml", 50.0, allow_mc=True),
        ScenarioSpec("Ar_100Td", ROOT / "configs" / "unified" / "both_template.yaml", 100.0),
        ScenarioSpec("Ar_300Td", ROOT / "configs" / "unified" / "both_template.yaml", 300.0),
        ScenarioSpec("Ar_N2_80Td", ROOT / "configs" / "unified" / "ar_n2_both.yaml", 80.0),
    )


def _set_nested(mapping: dict, *keys: str, value) -> None:
    current = mapping
    for key in keys[:-1]:
        current = current.setdefault(key, {})
    current[keys[-1]] = value


def _prepare_config(config: SwarmConfig, e_over_n_Td: float, *, quick: bool) -> SwarmConfig:
    cfg = copy.deepcopy(config)
    cfg.run.e_over_n_Td = [float(e_over_n_Td)]
    cfg.output.write_plots = False
    cfg.output.write_eedf = False
    cfg.output.write_rates = False

    boltz = cfg.boltzmann_two_term
    boltz.backend = "native_bolsig"
    boltz.adaptive_grid.enabled = False
    boltz.energy_grid.n = 160 if quick else max(300, boltz.energy_grid.n)
    boltz.energy_grid.max_eV = max(float(boltz.energy_grid.max_eV), 120.0)
    boltz.convergence.max_iterations = min(int(boltz.convergence.max_iterations), 120)

    mt = cfg.multiterm_boltzmann
    mt.enabled = True
    mt.method = "operator"
    mt.lmax = 1
    mt.hydrodynamic = False
    mt.allow_experimental_operator = False
    mt.energy_grid.n = 56 if quick else max(80, mt.energy_grid.n)
    mt.energy_grid.max_eV = max(float(mt.energy_grid.max_eV), 80.0)
    return cfg


def _prepare_mc_config(config: SwarmConfig, scenario_name: str, *, quick: bool) -> SwarmConfig:
    cfg = copy.deepcopy(config)
    cfg.run.case_prefix = scenario_name
    cfg.monte_carlo.enabled = True
    cfg.monte_carlo.python_api = "swarm_mc.unified_adapter:run_swarm"
    base = cfg.monte_carlo.passthrough.setdefault("base_config", {})
    scale = 1 if quick else 4
    _set_nested(base, "initial_state", "num_e_initial", value=1200 * scale)
    _set_nested(base, "simulation_settings", "num_energy_bins", value=128 if quick else 384)
    _set_nested(base, "simulation_settings", "num_e_max", value=60000 * scale)
    _set_nested(base, "simulation_settings", "seed", value=7)
    _set_nested(base, "end_conditions", "num_col_max", value=10000 * scale)
    _set_nested(base, "end_conditions", "w_tol", value=0.18 if quick else 0.08)
    _set_nested(base, "end_conditions", "DN_tol", value=0.18 if quick else 0.08)
    _set_nested(base, "output", "save_simulation_pickle", value=False)
    _set_nested(base, "output", "save_temporal_evolution", value=False)
    _set_nested(base, "output", "save_swarm_parameters", value=False)
    _set_nested(base, "output", "save_energy_distribution", value=False)
    return cfg


def _load_scenarios(*, quick: bool) -> list[LoadedScenario]:
    scenarios: list[LoadedScenario] = []
    for spec in _scenario_specs():
        cfg = _prepare_config(load_config(spec.config_path), spec.e_over_n_Td, quick=quick)
        cross_sections = load_cross_sections(cfg.cross_sections, cfg.conditions)
        scenarios.append(
            LoadedScenario(
                spec.name,
                cfg,
                cross_sections,
                spec.e_over_n_Td,
                allow_bolos=spec.allow_bolos,
                allow_mc=spec.allow_mc,
            )
        )
    return scenarios


def _relerr(candidate: float, reference: float) -> float:
    if not np.isfinite(candidate) or not np.isfinite(reference):
        return float("inf")
    return float(abs(candidate - reference) / max(abs(candidate), abs(reference), 1.0e-300))


def _row(
    scenario: str,
    candidate: str,
    reference: str,
    metric: str,
    candidate_value: float,
    reference_value: float,
    tolerance: float,
    status: str,
    note: str = "",
) -> dict[str, object]:
    finite_pair = np.isfinite(candidate_value) and np.isfinite(reference_value)
    return {
        "scenario": scenario,
        "candidate": candidate,
        "reference": reference,
        "metric": metric,
        "candidate_value": float(candidate_value) if np.isfinite(candidate_value) else None,
        "reference_value": float(reference_value) if np.isfinite(reference_value) else None,
        "relerr": _relerr(candidate_value, reference_value) if finite_pair else None,
        "tolerance": float(tolerance) if np.isfinite(tolerance) else None,
        "status": status,
        "note": note,
    }


def _comparison_rows(
    scenario: str,
    candidate_name: str,
    reference_name: str,
    candidate_case,
    reference_case,
    fields: Iterable[str],
    tolerance: float,
    *,
    report_only: bool = False,
    note: str = "",
) -> list[dict[str, object]]:
    rows = []
    for field in fields:
        candidate_value = float(getattr(candidate_case, field))
        reference_value = float(getattr(reference_case, field))
        relerr = _relerr(candidate_value, reference_value)
        passed = np.isfinite(relerr) and (report_only or relerr <= tolerance)
        status = "report" if report_only and passed else "ok" if passed else "FAIL"
        rows.append(
            _row(
                scenario,
                candidate_name,
                reference_name,
                field,
                candidate_value,
                reference_value,
                tolerance,
                status,
                note,
            )
        )
    return rows


def _native_reference(scenario: LoadedScenario):
    cfg = copy.deepcopy(scenario.config)
    cfg.boltzmann_two_term.backend = "native_bolsig"
    return BoltzmannTwoTermSolver(cfg, scenario.cross_sections).solve_native_reference_case(
        scenario.e_over_n_Td,
        f"{scenario.name}_native_two_term",
    )


def run_lmax1_gate(scenarios: Iterable[LoadedScenario]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for scenario in scenarios:
        reference = _native_reference(scenario)
        cfg = copy.deepcopy(scenario.config)
        cfg.multiterm_boltzmann.method = "operator"
        cfg.multiterm_boltzmann.lmax = 1
        cfg.multiterm_boltzmann.hydrodynamic = False
        candidate = MultiTermBoltzmannSolver(cfg, scenario.cross_sections).solve_case(
            scenario.e_over_n_Td,
            f"{scenario.name}_operator_lmax1",
        )
        rows.extend(
            _comparison_rows(
                scenario.name,
                "operator_lmax1",
                "native_two_term",
                candidate,
                reference,
                REFERENCE_FIELDS,
                LMAX1_TOLERANCE,
                note="deterministic lmax=1 reference gate",
            )
        )
    return rows


def run_lmax_gt1_report(scenarios: Iterable[LoadedScenario]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for scenario in scenarios:
        reference = _native_reference(scenario)
        cfg = copy.deepcopy(scenario.config)
        cfg.multiterm_boltzmann.method = "operator"
        cfg.multiterm_boltzmann.lmax = 2
        cfg.multiterm_boltzmann.allow_experimental_operator = True
        cfg.multiterm_boltzmann.hydrodynamic = False
        try:
            candidate = MultiTermBoltzmannSolver(cfg, scenario.cross_sections).solve_case(
                scenario.e_over_n_Td,
                f"{scenario.name}_operator_lmax2",
            )
        except Exception as exc:
            rows.append(
                _row(
                    scenario.name,
                    "operator_lmax2",
                    "native_two_term",
                    "all",
                    float("nan"),
                    float("nan"),
                    float("nan"),
                    "FAIL",
                    f"lmax>1 execution failed: {exc}",
                )
            )
            continue
        norm = float(candidate.metadata.get("normalization_integral", float("nan")))
        if not np.isfinite(norm) or abs(norm - 1.0) > 5.0e-6:
            rows.append(
                _row(
                    scenario.name,
                    "operator_lmax2",
                    "native_two_term",
                    "normalization_integral",
                    norm,
                    1.0,
                    5.0e-6,
                    "FAIL",
                    "lmax>1 operator normalization failed",
                )
            )
            continue
        rows.extend(
            _comparison_rows(
                scenario.name,
                "operator_lmax2",
                "native_two_term",
                candidate,
                reference,
                ("mean_energy_eV", "drift_velocity_m_s", "net_ionization_frequency_s"),
                float("nan"),
                report_only=True,
                note="reference-anchored lmax>1 report-only comparison",
            )
        )
    return rows


def run_bolos_gate(
    scenarios: Iterable[LoadedScenario],
    *,
    require_bolos: bool,
    rtol: float = BOLOS_TOLERANCE,
) -> list[dict[str, object]]:
    scenarios = [scenario for scenario in scenarios if scenario.allow_bolos]
    try:
        import bolos  # noqa: F401
    except Exception:
        status = "FAIL" if require_bolos else "skipped_optional"
        return [
            _row(
                scenario.name,
                "bolos",
                "native_two_term",
                "all",
                float("nan"),
                float("nan"),
                rtol,
                status,
                "optional bolos package is not installed",
            )
            for scenario in scenarios
        ]

    rows: list[dict[str, object]] = []
    for scenario in scenarios:
        reference = _native_reference(scenario)
        cfg = copy.deepcopy(scenario.config)
        cfg.boltzmann_two_term.backend = "bolos"
        try:
            candidate = BoltzmannTwoTermSolver(cfg, scenario.cross_sections).solve_case(
                scenario.e_over_n_Td,
                f"{scenario.name}_bolos",
            )
        except Exception as exc:
            rows.append(
                _row(
                    scenario.name,
                    "bolos",
                    "native_two_term",
                    "all",
                    float("nan"),
                    float("nan"),
                    rtol,
                    "FAIL",
                    f"BOLOS comparison failed: {exc}",
                )
            )
            continue
        rows.extend(
            _comparison_rows(
                scenario.name,
                "bolos",
                "native_two_term",
                candidate,
                reference,
                REFERENCE_FIELDS,
                rtol,
                note="optional BOLOS reference comparison",
            )
        )
    return rows


def _mc_absolute_error(case, metric: str) -> float | None:
    meta = case.metadata
    if metric == "mean_energy_eV":
        return meta.get("mean_energy_error_eV")
    if metric == "drift_velocity_m_s":
        return meta.get("flux_drift_velocity_error_m_s")
    if metric == "reduced_mobility_m2_V_s_m3":
        drift_err = meta.get("flux_drift_velocity_error_m_s")
        electric_field = meta.get("electric_field_V_m")
        density = meta.get("gas_number_density_m-3")
        if drift_err is not None and electric_field and density:
            return abs(float(drift_err) / float(electric_field) * float(density))
    if metric == "reduced_diffusion_L_m2_s_m3":
        return meta.get("flux_reduced_diffusion_L_error_m-1_s-1")
    if metric == "net_ionization_frequency_s":
        rate_err = meta.get("counted_effective_rate_error_m3_s")
        density = meta.get("gas_number_density_m-3")
        if rate_err is not None and density:
            return abs(float(rate_err) * float(density))
    return None


def _mc_tolerance(metric: str, candidate_case, reference_value: float) -> float:
    base = MC_BASE_TOLERANCES[metric]
    err = _mc_absolute_error(candidate_case, metric)
    if err is None or not np.isfinite(float(err)):
        return base
    candidate_value = float(getattr(candidate_case, metric))
    rel_sigma = abs(float(err)) / max(abs(candidate_value), abs(reference_value), 1.0e-300)
    return float(max(base, 3.0 * rel_sigma))


def run_mc_gate(
    scenarios: Iterable[LoadedScenario],
    *,
    quick: bool,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for scenario in scenarios:
        if not scenario.allow_mc:
            continue
        reference = _native_reference(scenario)
        cfg = _prepare_mc_config(scenario.config, scenario.name, quick=quick)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                [candidate] = MonteCarloAdapter(cfg, scenario.cross_sections).solve_all()
        except Exception as exc:
            rows.append(
                _row(
                    scenario.name,
                    "monte_carlo",
                    "native_two_term",
                    "all",
                    float("nan"),
                    float("nan"),
                    float("nan"),
                    "FAIL",
                    f"MC comparison failed: {exc}",
                )
            )
            continue
        for metric in MC_BASE_TOLERANCES:
            reference_value = float(getattr(reference, metric))
            tolerance = _mc_tolerance(metric, candidate, reference_value)
            rows.extend(
                _comparison_rows(
                    scenario.name,
                    "monte_carlo",
                    "native_two_term",
                    candidate,
                    reference,
                    (metric,),
                    tolerance,
                    note="optional stochastic MC comparison",
                )
            )
    return rows


def _is_failure(row: dict[str, object]) -> bool:
    return str(row.get("status", "")).startswith("FAIL")


def _as_jsonable(value):
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    return value


def write_json(rows: list[dict[str, object]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data = [{key: _as_jsonable(row.get(key)) for key in REPORT_COLUMNS} for row in rows]
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def _format_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def print_report(rows: list[dict[str, object]]) -> None:
    widths = {
        column: max(
            len(column),
            *(len(_format_value(row.get(column))) for row in rows),
        )
        for column in REPORT_COLUMNS
    }
    header = "  ".join(column.ljust(widths[column]) for column in REPORT_COLUMNS)
    print(header)
    print("  ".join("-" * widths[column] for column in REPORT_COLUMNS))
    for row in rows:
        print(
            "  ".join(
                _format_value(row.get(column)).ljust(widths[column])
                for column in REPORT_COLUMNS
            )
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="use small grids and short optional MC runs")
    parser.add_argument("--with-mc", action="store_true", help="include stochastic MC comparisons")
    parser.add_argument("--require-bolos", action="store_true", help="fail if BOLOS is unavailable")
    parser.add_argument("--bolos-rtol", type=float, default=BOLOS_TOLERANCE)
    parser.add_argument("--json", type=Path, help="write the gate report as JSON")
    args = parser.parse_args(argv)

    scenarios = _load_scenarios(quick=args.quick)
    rows: list[dict[str, object]] = []
    rows.extend(run_lmax1_gate(scenarios))
    rows.extend(run_bolos_gate(scenarios, require_bolos=args.require_bolos, rtol=args.bolos_rtol))
    rows.extend(run_lmax_gt1_report(scenarios))
    if args.with_mc:
        rows.extend(run_mc_gate(scenarios, quick=args.quick))

    print("\nMULTI-TERM OPERATOR BENCHMARK GATE")
    print("NOTE: lmax>1 is reference-anchored; comparisons are report-only.")
    print_report(rows)
    if args.json is not None:
        write_json(rows, args.json)
        print(f"\nWROTE: {args.json}")

    if any(_is_failure(row) for row in rows):
        print("\nFAIL: benchmark gate found hard failures")
        return 1
    print("\nPASS: benchmark gate completed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
