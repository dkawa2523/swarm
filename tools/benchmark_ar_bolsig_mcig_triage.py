"""Triage Ar comparisons across product solvers, BOLSIG+, and MCIG."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.eedf_compare import compare_eedf_cases
from tools.references.common import (
    ReferenceCaseResult,
    reference_comparison_metrics,
    reference_to_swarm_case,
)
from tools.benchmark_ar_external_references import (
    BOLSIG_THRESHOLDS,
    MCIG_THRESHOLDS,
    _angular_status,
    _failure,
    _mc_confidence_status,
    _write_plot,
    load_benchmark_config,
    load_selected_references,
    run_requested_external_reference_commands,
    run_solver_variants,
)
from tools.benchmark_common import write_csv


TRIAGE_FIELDS = [
    "E_over_N_Td",
    "overall_status",
    "two_term_vs_bolsig",
    "direct_lmax1_vs_two_term",
    "mcig_vs_bolsig",
    "multi_lmax_gt1_vs_mcig",
    "internal_mc_vs_mcig",
    "likely_cause",
    "implementation_bug",
    "physics_model_limitation",
    "mc_uncertainty_limited",
    "recommended_next_fix",
]

METRIC_FIELDS = [
    "E_over_N_Td",
    "reference",
    "candidate",
    "candidate_method",
    "candidate_lmax",
    "angular_model_status",
    "confidence_status",
    "eedf_relative_l1",
    "log_tail_error",
    "tail_probability_difference",
    "mean_energy_relative_difference",
    "drift_velocity_relative_difference",
    "mobility_relative_difference",
    "diffusion_L_relative_difference",
    "diffusion_T_relative_difference",
    "major_rate_relative_difference",
    "rate_weighted_eedf_error",
]

FAILURE_FIELDS = [
    "category",
    "evidence",
    "affected_E_over_N_Td",
    "affected_solver",
    "suspected_code_area",
    "recommended_fix",
    "severity",
]


def _solver_label(case) -> str:
    if case.solver == "multi_term":
        return f"multi_term_lmax{case.metadata.get('lmax', '')}"
    return str(case.solver)


def _metric_row(
    eover: float,
    reference: str,
    candidate,
    metrics: dict[str, float],
    *,
    angular_status: str = "",
    confidence: str = "",
) -> dict[str, object]:
    return {
        "E_over_N_Td": eover,
        "reference": reference,
        "candidate": _solver_label(candidate) if not isinstance(candidate, str) else candidate,
        "candidate_method": "" if isinstance(candidate, str) else candidate.metadata.get("solver_method", ""),
        "candidate_lmax": "" if isinstance(candidate, str) else candidate.metadata.get("lmax", ""),
        "angular_model_status": angular_status,
        "confidence_status": confidence,
        **metrics,
    }


def _is_bad(metrics: dict[str, float], thresholds: dict[str, float]) -> bool:
    return any(metrics.get(key, 0.0) > limit for key, limit in thresholds.items())


def _find_solver(cases, eover: float, solver: str, lmax: int | None = None):
    for case in cases:
        if case.solver != solver or abs(case.e_over_n_Td - eover) > 1.0e-9:
            continue
        if lmax is not None and str(case.metadata.get("lmax", "")) != str(lmax):
            continue
        return case
    return None


def _find_reference(cases: list[ReferenceCaseResult], eover: float, ref_id: str) -> ReferenceCaseResult | None:
    for case in cases:
        if case.reference_id == ref_id and abs(case.e_over_n_Td - eover) < 1.0e-9:
            return case
    return None


def _add_failure(rows: list[dict[str, object]], category: str, evidence: str, eover: float, pair: str, area: str, fix: str, severity: str) -> None:
    rows.append(_failure(category, evidence, eover, pair, area, fix, severity))


def run_triage(
    config_path: Path,
    *,
    bolsig_output: Path | None = None,
    mcig_output: Path | None = None,
    bolsig_command: str | None = None,
    mcig_command: str | None = None,
    bolsig_input: Path | None = None,
    mcig_input: Path | None = None,
    external_working_directory: Path | None = None,
    external_timeout_s: float | None = None,
    require_bolsig: bool = False,
    require_mcig: bool = False,
    fail_on_code_regression: bool = False,
    fail_on_physics_mismatch: bool = False,
    plot: bool = False,
) -> tuple[Path, ...]:
    cfg, reference_configs = load_benchmark_config(config_path)
    out_dir = cfg.output.directory
    run_requested_external_reference_commands(
        cfg,
        reference_configs=reference_configs,
        config_path=config_path,
        bolsig_command=bolsig_command,
        mcig_command=mcig_command,
        bolsig_output=bolsig_output,
        mcig_output=mcig_output,
        bolsig_input=bolsig_input,
        mcig_input=mcig_input,
        working_directory=external_working_directory,
        timeout_s=external_timeout_s,
    )
    references, ref_failures = load_selected_references(
        cfg,
        reference_configs=reference_configs,
        bolsig_output=bolsig_output,
        mcig_output=mcig_output,
        require_bolsig=require_bolsig,
        require_mcig=require_mcig,
    )
    solver_cases, solver_failures = run_solver_variants(cfg)

    matrix_rows: list[dict[str, object]] = []
    eedf_rows: list[dict[str, object]] = []
    transport_rows: list[dict[str, object]] = []
    rate_rows: list[dict[str, object]] = []
    failure_rows: list[dict[str, object]] = [*ref_failures, *solver_failures]

    for eover in cfg.run.e_over_n_Td:
        eover = float(eover)
        likely: list[str] = []
        implementation_bug = False
        physics_limit = False
        mc_uncertain = False
        statuses: dict[str, str] = {}

        two = _find_solver(solver_cases, eover, "two_term")
        direct1 = _find_solver(solver_cases, eover, "multi_term", 1)
        bolsig = _find_reference(references, eover, "bolsig_plus")
        mcig = _find_reference(references, eover, "mcig")

        if two is not None and direct1 is not None:
            metrics = compare_eedf_cases(two, direct1).metrics
            eedf_rows.append(_metric_row(eover, "two_term", direct1, metrics))
            statuses["direct_lmax1_vs_two_term"] = "PASS"
            if metrics.get("eedf_relative_l1", 0.0) > 1.0e-8:
                statuses["direct_lmax1_vs_two_term"] = "FAIL"
                implementation_bug = True
                likely.append("code_regression_multi_term_lmax1")
                _add_failure(
                    failure_rows,
                    "code_regression_multi_term_lmax1",
                    f"eedf_relative_l1={metrics['eedf_relative_l1']:.3e}",
                    eover,
                    "two_term->multi_term_lmax1",
                    "multi_term direct lmax1 reduction",
                    "inspect shared projection/source-sink/rate convolution",
                    "critical",
                )
        else:
            statuses["direct_lmax1_vs_two_term"] = "missing"
            if two is not None:
                implementation_bug = True
                likely.append("code_regression_multi_term_lmax1")
                _add_failure(
                    failure_rows,
                    "code_regression_multi_term_lmax1",
                    "multi_term pn_closure_direct lmax=1 result is missing",
                    eover,
                    "two_term->multi_term_lmax1",
                    "multi_term direct lmax1 execution",
                    "inspect solver_execution_failure rows and direct lmax=1 gate",
                    "critical",
                )

        if bolsig is not None and two is not None:
            metrics = reference_comparison_metrics(bolsig, two)
            eedf_rows.append(_metric_row(eover, "bolsig_plus", two, metrics))
            transport_rows.append(_metric_row(eover, "bolsig_plus", two, metrics))
            rate_rows.append(_metric_row(eover, "bolsig_plus", two, metrics))
            statuses["two_term_vs_bolsig"] = "PASS" if not _is_bad(metrics, BOLSIG_THRESHOLDS) else "FAIL"
            if statuses["two_term_vs_bolsig"] == "FAIL":
                implementation_bug = True
                likely.append("code_or_bolsig_input_mismatch")
                _add_failure(
                    failure_rows,
                    "code_or_bolsig_input_mismatch",
                    f"eedf_l1={metrics.get('eedf_relative_l1', 0.0):.3e}; mean={metrics.get('mean_energy_relative_difference', 0.0):.3e}",
                    eover,
                    "bolsig_plus->two_term",
                    "BOLSIG+ ingest / two_term projection",
                    "check EEDF convention, cross sections, grid, rates",
                    "high",
                )
        else:
            statuses["two_term_vs_bolsig"] = "missing"

        if bolsig is not None and mcig is not None:
            bolsig_case = reference_to_swarm_case(bolsig)
            metrics = reference_comparison_metrics(mcig, bolsig_case)
            eedf_rows.append(_metric_row(eover, "mcig", "reference:bolsig_plus", metrics))
            statuses["mcig_vs_bolsig"] = "PASS" if not _is_bad(metrics, MCIG_THRESHOLDS) else "DEGRADED"
            if statuses["mcig_vs_bolsig"] == "DEGRADED":
                physics_limit = True
                likely.append("possible_two_term_approximation_limit_or_angular_model_difference")
                _add_failure(
                    failure_rows,
                    "possible_two_term_approximation_limit_or_angular_model_difference",
                    f"eedf_l1={metrics.get('eedf_relative_l1', 0.0):.3e}",
                    eover,
                    "mcig->bolsig_plus",
                    "two-term approximation / angular model",
                    "compare same-angular MCIG metadata and BOLSIG+ assumptions",
                    "medium",
                )
        else:
            statuses["mcig_vs_bolsig"] = "missing"

        if mcig is not None:
            for candidate in [case for case in solver_cases if abs(case.e_over_n_Td - eover) < 1.0e-9]:
                metrics = reference_comparison_metrics(mcig, candidate)
                angular_status, angular_evidence = _angular_status(mcig, candidate)
                confidence, _ = _mc_confidence_status(mcig, candidate, metrics)
                eedf_rows.append(_metric_row(eover, "mcig", candidate, metrics, angular_status=angular_status, confidence=confidence))
                transport_rows.append(_metric_row(eover, "mcig", candidate, metrics, angular_status=angular_status, confidence=confidence))
                rate_rows.append(_metric_row(eover, "mcig", candidate, metrics, angular_status=angular_status, confidence=confidence))
                if confidence in {"unknown", "pass_within_mc_uncertainty"}:
                    mc_uncertain = True
                if confidence == "pass_within_mc_uncertainty":
                    likely.append("mc_uncertainty_limited")
                    _add_failure(
                        failure_rows,
                        "mc_uncertainty_limited",
                        "difference is within reported MCIG confidence interval",
                        eover,
                        f"mcig->{_solver_label(candidate)}",
                        "MCIG uncertainty",
                        "increase MCIG statistics before treating this as a solver mismatch",
                        "low",
                    )
                    continue
                if angular_status != "match":
                    likely.append("angular_model_mismatch")
                    _add_failure(
                        failure_rows,
                        "angular_model_mismatch",
                        angular_evidence or f"angular_model_status={angular_status}",
                        eover,
                        f"mcig->{_solver_label(candidate)}",
                        "angular scattering metadata",
                        "use same angular model before treating mismatch as code failure",
                        "medium",
                    )
                    continue
                if _is_bad(metrics, MCIG_THRESHOLDS):
                    if candidate.solver == "multi_term" and str(candidate.metadata.get("lmax", "")) not in {"", "1"}:
                        likely.append("multi_term_higher_l_model_issue")
                        _add_failure(
                            failure_rows,
                            "multi_term_higher_l_model_issue",
                            f"eedf_l1={metrics.get('eedf_relative_l1', 0.0):.3e}",
                            eover,
                            f"mcig->{_solver_label(candidate)}",
                            "higher-l direct PN angular closure",
                            "inspect higher-l damping/source-sink/tail behavior",
                            "medium",
                        )
                    elif candidate.solver == "monte_carlo":
                        likely.append("mc_uncertainty_limited")
                    elif candidate.solver == "two_term":
                        likely.append("BOLSIG_two_term_approximation_difference")
                        physics_limit = True

                if metrics.get("major_rate_relative_difference", 0.0) > MCIG_THRESHOLDS["major_rate_relative_difference"] and metrics.get("eedf_relative_l1", 0.0) <= MCIG_THRESHOLDS["eedf_relative_l1"]:
                    likely.append("rate_convolution_or_cross_section_projection_issue")
                    _add_failure(
                        failure_rows,
                        "rate_convolution_or_cross_section_projection_issue",
                        f"rate_diff={metrics.get('major_rate_relative_difference', 0.0):.3e}; eedf_l1={metrics.get('eedf_relative_l1', 0.0):.3e}",
                        eover,
                        f"mcig->{_solver_label(candidate)}",
                        "rate convolution / cross-section projection",
                        "compare rate labels, thresholds, and mixture weighting",
                        "medium",
                    )
                tail_mismatch = (
                    metrics.get("log_tail_error", 0.0) > 2.0
                    or metrics.get("tail_probability_difference", 0.0) > 0.02
                    or metrics.get("E99_difference_eV", 0.0) > 0.5
                )
                if tail_mismatch:
                    likely.append("tail_boundary_or_sampling_issue")
                    _add_failure(
                        failure_rows,
                        "tail_boundary_or_sampling_issue",
                        "log_tail_error="
                        f"{metrics.get('log_tail_error', 0.0):.3e}; "
                        "tail_probability_difference="
                        f"{metrics.get('tail_probability_difference', 0.0):.3e}; "
                        f"E99_difference_eV={metrics.get('E99_difference_eV', 0.0):.3e}",
                        eover,
                        f"mcig->{_solver_label(candidate)}",
                        "tail boundary / MC sampling",
                        "inspect high-energy cutoff and MCIG tail statistics",
                        "medium",
                    )

            high_l_rows = [
                row for row in eedf_rows
                if row["E_over_N_Td"] == eover
                and row["reference"] == "mcig"
                and str(row["candidate"]).startswith("multi_term_lmax")
                and row["candidate"] != "multi_term_lmax1"
            ]
            statuses["multi_lmax_gt1_vs_mcig"] = (
                "PASS" if high_l_rows and all(float(row.get("eedf_relative_l1", 0.0)) <= MCIG_THRESHOLDS["eedf_relative_l1"] for row in high_l_rows)
                else "DEGRADED" if high_l_rows else "missing"
            )
            mc_rows = [
                row for row in eedf_rows
                if row["E_over_N_Td"] == eover and row["reference"] == "mcig" and row["candidate"] == "monte_carlo"
            ]
            statuses["internal_mc_vs_mcig"] = (
                "PASS" if mc_rows and all(float(row.get("eedf_relative_l1", 0.0)) <= MCIG_THRESHOLDS["eedf_relative_l1"] for row in mc_rows)
                else "DEGRADED" if mc_rows else "missing"
            )
        else:
            statuses["multi_lmax_gt1_vs_mcig"] = "missing"
            statuses["internal_mc_vs_mcig"] = "missing"

        if not likely:
            likely.append("no_external_reference_available" if bolsig is None and mcig is None else "pass")
        overall = "FAIL" if implementation_bug else "DEGRADED" if physics_limit or mc_uncertain or any("mismatch" in item for item in likely) else "PASS"
        matrix_rows.append(
            {
                "E_over_N_Td": eover,
                "overall_status": overall,
                "two_term_vs_bolsig": statuses.get("two_term_vs_bolsig", "missing"),
                "direct_lmax1_vs_two_term": statuses.get("direct_lmax1_vs_two_term", "missing"),
                "mcig_vs_bolsig": statuses.get("mcig_vs_bolsig", "missing"),
                "multi_lmax_gt1_vs_mcig": statuses.get("multi_lmax_gt1_vs_mcig", "missing"),
                "internal_mc_vs_mcig": statuses.get("internal_mc_vs_mcig", "missing"),
                "likely_cause": ";".join(dict.fromkeys(likely)),
                "implementation_bug": implementation_bug,
                "physics_model_limitation": physics_limit,
                "mc_uncertainty_limited": mc_uncertain,
                "recommended_next_fix": _recommend(likely),
            }
        )

    paths = {
        "matrix": out_dir / "ar_triage_matrix.csv",
        "eedf": out_dir / "ar_triage_eedf_metrics.csv",
        "transport": out_dir / "ar_triage_transport_metrics.csv",
        "rates": out_dir / "ar_triage_rate_metrics.csv",
        "failures": out_dir / "ar_triage_failure_analysis.csv",
        "report": out_dir / "ar_triage_report.md",
    }
    write_csv(paths["matrix"], matrix_rows, TRIAGE_FIELDS)
    write_csv(paths["eedf"], eedf_rows, METRIC_FIELDS)
    write_csv(paths["transport"], transport_rows, METRIC_FIELDS)
    write_csv(paths["rates"], rate_rows, METRIC_FIELDS)
    write_csv(paths["failures"], failure_rows, FAILURE_FIELDS)
    status = "PASS"
    if any(row["overall_status"] == "FAIL" for row in matrix_rows):
        status = "FAIL"
    elif any(row["overall_status"] == "DEGRADED" for row in matrix_rows):
        status = "DEGRADED"
    if all(row["likely_cause"] == "no_external_reference_available" for row in matrix_rows):
        status = "SKIP_EXTERNAL_REFERENCE"
    written_plot = _write_plot(out_dir, "ar_triage", references, solver_cases) if plot else None
    paths["report"].write_text(
        "\n".join(
            [
                "# Ar BOLSIG+ / MCIG Triage",
                "",
                f"Config: `{config_path}`",
                f"Status: **{status}**",
                f"Rows: {len(matrix_rows)}",
                f"Failure rows: {len(failure_rows)}",
                f"Plot: `{written_plot}`" if written_plot is not None else "Plot: not requested or no reference cases",
                "",
                "This benchmark classifies code regressions, physics-model limits, "
                "angular-model mismatch, and MC uncertainty without fabricating "
                "missing external references.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    if fail_on_code_regression and any(row["implementation_bug"] for row in matrix_rows):
        raise RuntimeError("triage detected code regression")
    if fail_on_physics_mismatch and any(row["physics_model_limitation"] for row in matrix_rows):
        raise RuntimeError("triage detected physics-model mismatch")
    output_paths = list(paths.values())
    if written_plot is not None:
        output_paths.append(written_plot)
    return tuple(output_paths)


def _recommend(likely: list[str]) -> str:
    causes = set(likely)
    if "code_regression_multi_term_lmax1" in causes:
        return "fix direct lmax=1 shared projection/regression first"
    if "code_or_bolsig_input_mismatch" in causes:
        return "inspect BOLSIG+ EEDF convention, cross sections, grid, and rates"
    if "multi_term_higher_l_model_issue" in causes:
        return "inspect higher-l damping/source-sink/tail model"
    if "rate_convolution_or_cross_section_projection_issue" in causes:
        return "compare rate convolution and cross-section projection"
    if "tail_boundary_or_sampling_issue" in causes:
        return "inspect high-energy boundary and MC tail statistics"
    if "angular_model_mismatch" in causes:
        return "rerun with same angular scattering metadata"
    if "mc_uncertainty_limited" in causes:
        return "increase MCIG/internal MC statistics and provide confidence intervals"
    return "no code fix indicated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("configs/benchmarks/ar_bolsig_mcig_triage.yaml"))
    parser.add_argument("--bolsig-output", type=Path, default=None)
    parser.add_argument("--mcig-output", type=Path, default=None)
    parser.add_argument("--run-bolsig", default=None, help="external BOLSIG+ command template")
    parser.add_argument("--run-mcig", default=None, help="external MCIG command template")
    parser.add_argument("--bolsig-input", type=Path, default=None)
    parser.add_argument("--mcig-input", type=Path, default=None)
    parser.add_argument("--external-working-directory", type=Path, default=None)
    parser.add_argument("--external-timeout-s", type=float, default=None)
    parser.add_argument("--require-bolsig", action="store_true")
    parser.add_argument("--require-mcig", action="store_true")
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--fail-on-code-regression", action="store_true")
    parser.add_argument("--fail-on-physics-mismatch", action="store_true")
    args = parser.parse_args()
    for path in run_triage(
        args.config,
        bolsig_output=args.bolsig_output,
        mcig_output=args.mcig_output,
        bolsig_command=args.run_bolsig,
        mcig_command=args.run_mcig,
        bolsig_input=args.bolsig_input,
        mcig_input=args.mcig_input,
        external_working_directory=args.external_working_directory,
        external_timeout_s=args.external_timeout_s,
        require_bolsig=args.require_bolsig,
        require_mcig=args.require_mcig,
        fail_on_code_regression=args.fail_on_code_regression,
        fail_on_physics_mismatch=args.fail_on_physics_mismatch,
        plot=args.plot,
    ):
        print(path)


if __name__ == "__main__":
    main()
