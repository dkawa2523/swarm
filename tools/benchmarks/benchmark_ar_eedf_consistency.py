"""Run the Ar product EEDF consistency benchmark."""

# The script must make the repository package importable before local imports.
# ruff: noqa: E402

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from electron_swarm import run
from electron_swarm.core.config import SwarmConfig
from tools.benchmarks.eedf_compare import compare_eedf_cases
from tools.benchmarks.references import load_reference_cases
from tools.benchmarks.references.common import (
    ExternalReferenceConfig,
    reference_comparison_metrics,
)
from tools.benchmarks.benchmark_ar_external_references import load_benchmark_config
from tools.benchmarks.common import variant_config, write_csv


PN_RESIDUAL_LIMIT = 1.0e-6
PN_NEGATIVE_MASS_LIMIT = 1.0e-8
PN_LMAX4_TO_6_EEDF_L1_LIMIT = 0.01

REFERENCE_COMMON_FIELDS = [
    "reference_id",
    "reference_case_id",
    "candidate_solver",
    "candidate_method",
    "candidate_lmax",
    "E_over_N_Td",
    "reference_format",
    "reference_role",
    "uncertainty_unavailable",
    "reported_mean_energy_difference_eV",
]

REFERENCE_SUMMARY_FIELDS = [
    *REFERENCE_COMMON_FIELDS,
    "mean_energy_relative_difference",
    "drift_velocity_relative_difference",
    "mobility_relative_difference",
    "diffusion_L_relative_difference",
    "diffusion_T_relative_difference",
    "major_rate_relative_difference",
    "eedf_relative_l1",
]

REFERENCE_METRIC_FIELDS = [
    *REFERENCE_COMMON_FIELDS,
    "eedf_l1",
    "eedf_relative_l1",
    "eedf_l2_weighted",
    "log_tail_error",
    "mean_energy_relative_difference",
    "eedf_mean_energy_relative_difference",
    "tail_probability_difference",
    "rate_weighted_eedf_error",
    "drift_velocity_relative_difference",
    "major_rate_relative_difference",
    "normalization_error_reference",
    "normalization_error_candidate",
    "E50_difference_eV",
    "E90_difference_eV",
    "E99_difference_eV",
    "mobility_relative_difference",
    "diffusion_L_relative_difference",
    "diffusion_T_relative_difference",
]


def _run_single(cfg: SwarmConfig):
    result = run(cfg, write=False)
    return result.cases[0] if result.cases else None


def _reference_rows(
    cfg: SwarmConfig,
    reference_configs: list[ExternalReferenceConfig],
    solver_cases,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    summary_rows: list[dict[str, object]] = []
    metric_rows: list[dict[str, object]] = []
    failure_rows: list[dict[str, object]] = []
    if not reference_configs:
        return summary_rows, metric_rows, failure_rows
    for ref_config in reference_configs:
        try:
            reference_cases = load_reference_cases(ref_config)
        except FileNotFoundError as exc:
            failure_rows.append(
                {
                    "category": "missing_external_reference",
                    "evidence": str(exc),
                    "affected_E_over_N_Td": "",
                    "affected_solver": f"reference:{ref_config.id}",
                    "suspected_code_area": "external reference ingest",
                    "recommended_fix": "provide the configured BOLSIG+/MCIG output file or remove the reference entry",
                    "severity": "medium",
                }
            )
            continue
        except ValueError as exc:
            failure_rows.append(
                {
                    "category": "unsupported_external_reference_format",
                    "evidence": str(exc),
                    "affected_E_over_N_Td": "",
                    "affected_solver": f"reference:{ref_config.id}",
                    "suspected_code_area": "external reference ingest",
                    "recommended_fix": "convert the external output to electron_swarm_reference_csv",
                    "severity": "high",
                }
            )
            continue
        for reference in reference_cases:
            matched = [
                case
                for case in solver_cases
                if abs(float(case.e_over_n_Td) - float(reference.e_over_n_Td)) < 1.0e-9
            ]
            if not matched:
                failure_rows.append(
                    {
                        "category": "missing_solver_case_for_reference",
                        "evidence": f"no solver result at E/N={reference.e_over_n_Td}",
                        "affected_E_over_N_Td": reference.e_over_n_Td,
                        "affected_solver": f"reference:{reference.reference_id}",
                        "suspected_code_area": "benchmark reference matching",
                        "recommended_fix": "include the reference E/N in run.e_over_n_Td",
                        "severity": "medium",
                    }
                )
                continue
            for case in matched:
                metrics = reference_comparison_metrics(reference, case)
                common = {
                    "reference_id": reference.reference_id,
                    "reference_case_id": reference.case_id,
                    "candidate_solver": case.solver,
                    "candidate_method": case.metadata.get("solver_method", ""),
                    "candidate_lmax": case.metadata.get("lmax", ""),
                    "E_over_N_Td": case.e_over_n_Td,
                    "reference_format": reference.metadata.get("reference_format", ""),
                    "reference_role": reference.metadata.get("reference_role", ""),
                    "uncertainty_unavailable": reference.metadata.get(
                        "uncertainty_unavailable", ""
                    ),
                    "reported_mean_energy_difference_eV": reference.metadata.get(
                        "reported_mean_energy_difference_eV", ""
                    ),
                }
                summary_rows.append(
                    {
                        **common,
                        "mean_energy_relative_difference": metrics.get(
                            "mean_energy_relative_difference", ""
                        ),
                        "drift_velocity_relative_difference": metrics.get(
                            "drift_velocity_relative_difference", ""
                        ),
                        "mobility_relative_difference": metrics.get(
                            "mobility_relative_difference", ""
                        ),
                        "diffusion_L_relative_difference": metrics.get(
                            "diffusion_L_relative_difference", ""
                        ),
                        "diffusion_T_relative_difference": metrics.get(
                            "diffusion_T_relative_difference", ""
                        ),
                        "major_rate_relative_difference": metrics.get(
                            "major_rate_relative_difference", ""
                        ),
                        "eedf_relative_l1": metrics.get("eedf_relative_l1", ""),
                    }
                )
                metric_rows.append({**common, **metrics})
    return summary_rows, metric_rows, failure_rows


def run_benchmark(config_path: Path) -> tuple[Path, ...]:
    cfg, reference_configs = load_benchmark_config(config_path)
    references = [
        _run_single(variant_config(cfg, "two_term")) for _ in cfg.run.e_over_n_Td[:1]
    ]
    reference = references[0]
    if reference is None:
        raise RuntimeError("Ar consistency benchmark requires a two_term reference")
    candidates = [
        ("multi_term", "pn_closure_direct", 1),
        ("multi_term", "pn_closure_direct", 2),
        ("multi_term", "pn_closure_direct", 4),
        ("multi_term", "pn_closure_direct", 6),
    ]
    if any(item.id == "monte_carlo" for item in cfg.run.solvers):
        candidates.append(("monte_carlo", None, None))
    metric_rows: list[dict[str, object]] = []
    failure_rows: list[dict[str, object]] = []
    solver_cases = [reference]
    direct_cases: dict[int, object] = {}
    previous_direct_case = None
    for solver, method, lmax in candidates:
        try:
            case = _run_single(variant_config(cfg, solver, method=method, lmax=lmax))
        except NotImplementedError as exc:
            category = (
                "higher_l_field_coupling_error"
                if solver == "multi_term"
                and method == "pn_closure_direct"
                and (lmax or 1) > 1
                else "collision_moment_or_angular_closure_error"
            )
            failure_rows.append(
                {
                    "category": category,
                    "evidence": str(exc),
                    "affected_E_over_N_Td": cfg.run.e_over_n_Td[0],
                    "affected_solver": solver,
                    "suspected_code_area": "multi_term direct PN higher-l closure",
                    "recommended_fix": (
                        "derive coefficient-space field coupling and validate "
                        "higher-l damping, source/sink, and boundary conditions"
                    ),
                    "severity": "medium",
                }
            )
            continue
        if case is None:
            continue
        solver_cases.append(case)
        comparison = compare_eedf_cases(reference, case)
        previous_lmax_l1 = ""
        if solver == "multi_term" and method == "pn_closure_direct":
            direct_cases[int(lmax or 1)] = case
            if previous_direct_case is not None:
                previous_lmax_l1 = compare_eedf_cases(
                    previous_direct_case,
                    case,
                ).metrics["eedf_relative_l1"]
            previous_direct_case = case
        metric_rows.append(
            {
                "reference_solver": "two_term",
                "candidate_solver": case.solver,
                "candidate_method": case.metadata.get("solver_method", method or ""),
                "candidate_lmax": case.metadata.get("lmax", lmax or ""),
                "E_over_N_Td": case.e_over_n_Td,
                "eedf_relative_l1_vs_previous_lmax": previous_lmax_l1,
                "pn_residual": case.diagnostics.get("multi_term", {}).get(
                    "pn_full_relative_residual", ""
                ),
                "negative_mass_fraction": case.diagnostics.get(
                    "multi_term", {}
                ).get(
                    "negative_mass_fraction", ""
                ),
                **comparison.metrics,
            }
        )
        failure_rows.extend(comparison.failures)

    out_dir = cfg.output.directory
    metrics_path = out_dir / "ar_eedf_consistency_eedf_metrics.csv"
    failures_path = out_dir / "ar_eedf_consistency_failure_analysis.csv"
    report_path = out_dir / "ar_eedf_consistency_report.md"
    ref_summary_path = out_dir / "ar_reference_comparison_summary.csv"
    ref_metrics_path = out_dir / "ar_reference_eedf_metrics.csv"
    ref_failures_path = out_dir / "ar_reference_failure_analysis.csv"
    ref_report_path = out_dir / "ar_reference_report.md"
    ref_summary_rows, ref_metric_rows, ref_failure_rows = _reference_rows(
        cfg,
        reference_configs,
        solver_cases,
    )
    write_csv(metrics_path, metric_rows)
    write_csv(
        failures_path,
        failure_rows,
        fields=[
            "category",
            "evidence",
            "affected_E_over_N_Td",
            "affected_solver",
            "suspected_code_area",
            "recommended_fix",
            "severity",
        ],
    )
    write_csv(ref_summary_path, ref_summary_rows, REFERENCE_SUMMARY_FIELDS)
    write_csv(ref_metrics_path, ref_metric_rows, REFERENCE_METRIC_FIELDS)
    write_csv(
        ref_failures_path,
        ref_failure_rows,
        fields=[
            "category",
            "evidence",
            "affected_E_over_N_Td",
            "affected_solver",
            "suspected_code_area",
            "recommended_fix",
            "severity",
        ],
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)
    pn_gate_failures = []
    for lmax in (1, 2, 4, 6):
        case = direct_cases.get(lmax)
        if case is None:
            pn_gate_failures.append(f"pn_closure_direct lmax={lmax} did not run")
            continue
        diagnostics = case.diagnostics.get("multi_term", {})
        residual = float(diagnostics.get("pn_full_relative_residual", float("inf")))
        negative_mass = float(
            diagnostics.get("negative_mass_fraction", float("inf"))
        )
        if not diagnostics.get("pn_converged", False):
            pn_gate_failures.append(f"lmax={lmax} did not report convergence")
        if residual > PN_RESIDUAL_LIMIT:
            pn_gate_failures.append(
                f"lmax={lmax} residual={residual:.3e} > {PN_RESIDUAL_LIMIT:.3e}"
            )
        if negative_mass > PN_NEGATIVE_MASS_LIMIT:
            pn_gate_failures.append(
                f"lmax={lmax} negative_mass={negative_mass:.3e} > "
                f"{PN_NEGATIVE_MASS_LIMIT:.3e}"
            )
    if 4 in direct_cases and 6 in direct_cases:
        lmax_change = compare_eedf_cases(
            direct_cases[4], direct_cases[6]
        ).metrics["eedf_relative_l1"]
        if lmax_change > PN_LMAX4_TO_6_EEDF_L1_LIMIT:
            pn_gate_failures.append(
                f"lmax4_to_6_eedf_l1={lmax_change:.3e} > "
                f"{PN_LMAX4_TO_6_EEDF_L1_LIMIT:.3e}"
            )
    status = "PASS" if not pn_gate_failures else "FAIL"
    report_path.write_text(
        "\n".join(
            [
                "# Ar EEDF Consistency",
                "",
                f"Config: `{config_path}`",
                f"Status: **{status}**",
                "Full PN numerical gate: "
                + ("PASS" if not pn_gate_failures else "FAIL"),
                *[f"- {item}" for item in pn_gate_failures],
                f"Metric rows: {len(metric_rows)}",
                f"Diagnostic failure rows: {len(failure_rows)}",
                "",
                "The comparison uses normalized EEDF `F(E)` in `1/eV`, not EEPF.",
                "The two-term comparison is diagnostic. Every multi-term order "
                "uses the independent staggered PN operator and is accepted from "
                "its full-system residual, negative mass, and lmax convergence.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    ref_report_path.write_text(
        "\n".join(
            [
                "# Ar External Reference Comparison",
                "",
                f"Configured external references: {len(reference_configs)}",
                f"Reference comparison rows: {len(ref_summary_rows)}",
                f"Reference metric rows: {len(ref_metric_rows)}",
                f"Reference failure rows: {len(ref_failure_rows)}",
                "",
                "External references are benchmark inputs, not solver modes.",
                "All external EEDF/EEPF data are converted to normalized EEDF "
                "`F(E)` in `1/eV` before comparison.",
                "Missing BOLSIG+/MCIG output files are reported here and do not "
                "fabricate reference data.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return (
        metrics_path,
        failures_path,
        report_path,
        ref_summary_path,
        ref_metrics_path,
        ref_failures_path,
        ref_report_path,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/benchmarks/ar_bolsig_eedf_consistency.yaml"),
    )
    args = parser.parse_args()
    paths = run_benchmark(args.config)
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
