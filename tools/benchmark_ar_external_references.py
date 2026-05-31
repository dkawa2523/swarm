"""Run Ar external-reference equivalence benchmarks."""

from __future__ import annotations

import argparse
import copy
import csv
from pathlib import Path

from electron_swarm import load_config, run
from electron_swarm.core.config import (
    ExternalReferenceConfig,
    RequestedSolverConfig,
    SwarmConfig,
)
from electron_swarm.diagnostics.eedf_compare import compare_eedf_cases
from electron_swarm.references import load_reference_cases
from electron_swarm.references.common import (
    ReferenceCaseResult,
    reference_comparison_metrics,
)
from electron_swarm.references.runner import run_external_reference_command


BOLSIG_THRESHOLDS = {
    "mean_energy_relative_difference": 0.01,
    "major_rate_relative_difference": 0.03,
    "eedf_relative_l1": 0.05,
}
MCIG_THRESHOLDS = {
    "mean_energy_relative_difference": 0.05,
    "major_rate_relative_difference": 0.10,
    "eedf_relative_l1": 0.10,
}

DIRECT_LMAX1_THRESHOLDS = {
    "eedf_relative_l1": 1.0e-8,
    "mean_energy_relative_difference": 0.005,
    "drift_velocity_relative_difference": 0.01,
    "major_rate_relative_difference": 0.02,
    "normalization_error_candidate": 1.0e-8,
}

SUMMARY_FIELDS = [
    "reference_id",
    "reference_case_id",
    "candidate_solver",
    "candidate_method",
    "candidate_lmax",
    "E_over_N_Td",
    "status",
    "confidence_status",
    "same_angular_model",
    "angular_model_reference",
    "angular_model_candidate",
    "angular_model_mismatch_reason",
    "mean_energy_relative_difference",
    "drift_velocity_relative_difference",
    "mobility_relative_difference",
    "diffusion_L_relative_difference",
    "diffusion_T_relative_difference",
    "major_rate_relative_difference",
    "ionization_rate_relative_difference",
    "excitation_rate_relative_difference",
    "eedf_relative_l1",
    "log_tail_error",
    "tail_probability_difference",
    "mc_confidence_interval_coverage",
]

METRIC_FIELDS = [
    *SUMMARY_FIELDS,
    "eedf_l1",
    "eedf_l2_weighted",
    "eedf_mean_energy_relative_difference",
    "rate_weighted_eedf_error",
    "normalization_error_reference",
    "normalization_error_candidate",
    "E50_difference_eV",
    "E90_difference_eV",
    "E99_difference_eV",
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


def _write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _variant_config(
    cfg: SwarmConfig,
    solver: str,
    *,
    method: str | None = None,
    lmax: int | None = None,
    e_over_n_Td: float | None = None,
) -> SwarmConfig:
    variant = copy.deepcopy(cfg)
    variant.run.solvers = [RequestedSolverConfig(id=solver)]  # type: ignore[arg-type]
    if e_over_n_Td is not None:
        variant.run.e_over_n_Td = [float(e_over_n_Td)]
    if solver == "multi_term":
        variant.solvers.multi_term.method = method or "pn_closure_direct"  # type: ignore[assignment]
        variant.internal.multi_term.product_method = variant.solvers.multi_term.method
        if lmax is not None:
            variant.solvers.multi_term.lmax = int(lmax)
            variant.internal.multi_term.lmax = int(lmax)
    return variant


def _run_one(
    cfg: SwarmConfig,
    solver: str,
    *,
    method: str | None = None,
    lmax: int | None = None,
    e_over_n_Td: float | None = None,
):
    result = run(
        _variant_config(
            cfg,
            solver,
            method=method,
            lmax=lmax,
            e_over_n_Td=e_over_n_Td,
        ),
        write=False,
    )
    return result.cases


def _reference_key(name: str) -> str:
    if name in {"bolsig", "bolsig_plus"}:
        return "bolsig_plus"
    if name == "mcig":
        return "mcig"
    if name == "all":
        return "all"
    raise ValueError("--reference must be bolsig, mcig, or all")


def _matching_references(
    cfg: SwarmConfig,
    *,
    reference: str,
    bolsig_output: Path | None,
    mcig_output: Path | None,
    require_bolsig: bool,
    require_mcig: bool,
    angular_model: str | None,
) -> tuple[list[ReferenceCaseResult], list[dict[str, object]]]:
    desired = _reference_key(reference)
    configs = [
        ref
        for ref in cfg.references.external
        if desired == "all" or ref.id == desired or (desired == "mcig" and ref.id == "bolsig_plus")
    ]
    if bolsig_output is not None:
        configs = [
            ref for ref in configs if ref.id != "bolsig_plus"
        ]
        configs.append(
            ExternalReferenceConfig(
                id="bolsig_plus",
                path=bolsig_output,
                format="electron_swarm_reference_csv",
                eedf_convention="eedf",
            )
        )
    if mcig_output is not None:
        configs = [ref for ref in configs if ref.id != "mcig"]
        configs.append(
            ExternalReferenceConfig(
                id="mcig",
                path=mcig_output,
                format="electron_swarm_reference_csv",
                eedf_convention="eedf",
                angular_model=angular_model or "unknown",
            )
        )
    failures: list[dict[str, object]] = []
    cases: list[ReferenceCaseResult] = []
    if not configs and (desired in {"all", "bolsig_plus"}):
        failures.append(
            _failure(
                "missing_external_reference",
                "no BOLSIG+ reference path configured",
                "",
                "reference:bolsig_plus",
                "benchmark reference config",
                "provide --bolsig-output or references.external entry",
                "high" if require_bolsig else "medium",
            )
        )
        if require_bolsig:
            raise FileNotFoundError("BOLSIG+ reference is required but not configured")
    if not any(ref.id == "mcig" for ref in configs) and desired in {"all", "mcig"}:
        failures.append(
            _failure(
                "missing_external_reference",
                "no MCIG reference path configured",
                "",
                "reference:mcig",
                "benchmark reference config",
                "provide --mcig-output or references.external entry",
                "high" if require_mcig else "medium",
            )
        )
        if require_mcig:
            raise FileNotFoundError("MCIG reference is required but not configured")
    for ref_config in configs:
        if ref_config.id == "mcig" and angular_model is not None:
            ref_config = ExternalReferenceConfig(
                id=ref_config.id,
                path=ref_config.path,
                format=ref_config.format,
                eedf_convention=ref_config.eedf_convention,
                uncertainty=ref_config.uncertainty,
                angular_model=angular_model,  # type: ignore[arg-type]
            )
        try:
            cases.extend(load_reference_cases(ref_config))
        except FileNotFoundError as exc:
            failures.append(
                _failure(
                    "missing_external_reference",
                    str(exc),
                    "",
                    f"reference:{ref_config.id}",
                    "external reference ingest",
                    "provide the configured BOLSIG+/MCIG output file",
                    "high"
                    if (
                        (require_bolsig and ref_config.id == "bolsig_plus")
                        or (require_mcig and ref_config.id == "mcig")
                    )
                    else "medium",
                )
            )
            if require_bolsig and ref_config.id == "bolsig_plus":
                raise
            if require_mcig and ref_config.id == "mcig":
                raise
        except ValueError as exc:
            failures.append(
                _failure(
                    "unsupported_external_reference_format",
                    str(exc),
                    "",
                    f"reference:{ref_config.id}",
                    "external reference ingest",
                    "convert external output to electron_swarm_reference_csv",
                    "high",
                )
            )
            if require_bolsig and ref_config.id == "bolsig_plus":
                raise
            if require_mcig and ref_config.id == "mcig":
                raise
    return cases, failures


def load_selected_references(
    cfg: SwarmConfig,
    *,
    bolsig_output: Path | None = None,
    mcig_output: Path | None = None,
    require_bolsig: bool = False,
    require_mcig: bool = False,
    angular_model: str | None = None,
) -> tuple[list[ReferenceCaseResult], list[dict[str, object]]]:
    return _matching_references(
        cfg,
        reference="all",
        bolsig_output=bolsig_output,
        mcig_output=mcig_output,
        require_bolsig=require_bolsig,
        require_mcig=require_mcig,
        angular_model=angular_model,
    )


def _configured_reference_path(cfg: SwarmConfig, reference_id: str) -> Path | None:
    for ref in cfg.references.external:
        if ref.id == reference_id:
            return ref.path
    return None


def run_requested_external_reference_commands(
    cfg: SwarmConfig,
    *,
    config_path: Path,
    bolsig_command: str | None = None,
    mcig_command: str | None = None,
    bolsig_output: Path | None = None,
    mcig_output: Path | None = None,
    bolsig_input: Path | None = None,
    mcig_input: Path | None = None,
    working_directory: Path | None = None,
    timeout_s: float | None = None,
) -> None:
    if bolsig_command is not None:
        output = bolsig_output or _configured_reference_path(cfg, "bolsig_plus")
        if output is None:
            raise ValueError(
                "--run-bolsig requires --bolsig-output or a bolsig_plus "
                "references.external path"
            )
        run_external_reference_command(
            reference_id="bolsig_plus",
            command=bolsig_command,
            output=output,
            input_path=bolsig_input,
            config_path=config_path,
            working_directory=working_directory,
            timeout_s=timeout_s,
        )
    if mcig_command is not None:
        output = mcig_output or _configured_reference_path(cfg, "mcig")
        if output is None:
            raise ValueError(
                "--run-mcig requires --mcig-output or an mcig references.external path"
            )
        run_external_reference_command(
            reference_id="mcig",
            command=mcig_command,
            output=output,
            input_path=mcig_input,
            config_path=config_path,
            working_directory=working_directory,
            timeout_s=timeout_s,
        )


def run_solver_variants(cfg: SwarmConfig) -> tuple[list[object], list[dict[str, object]]]:
    failures: list[dict[str, object]] = []
    solver_cases = []
    solver_variants = [
        ("two_term", None, None),
        ("multi_term", "pn_closure_direct", 1),
        ("multi_term", "pn_closure_direct", 2),
        ("multi_term", "pn_closure_direct", 4),
        ("multi_term", "pn_closure_direct", 6),
    ]
    if any(item.id == "monte_carlo" for item in cfg.run.solvers):
        solver_variants.append(("monte_carlo", None, None))
    for solver, method, lmax in solver_variants:
        for eover in cfg.run.e_over_n_Td:
            try:
                solver_cases.extend(
                    _run_one(
                        cfg,
                        solver,
                        method=method,
                        lmax=lmax,
                        e_over_n_Td=eover,
                    )
                )
            except Exception as exc:
                affected_solver = f"{solver}_lmax{lmax}" if solver == "multi_term" and lmax is not None else solver
                failures.append(
                    _failure(
                        "solver_execution_failure",
                        str(exc),
                        eover,
                        affected_solver,
                        "benchmark solver execution",
                        "inspect solver config, grid range, and high_energy_extrapolation",
                        "high",
                    )
                )
    return solver_cases, failures


def _failure(
    category: str,
    evidence: str,
    eover: object,
    solver: str,
    area: str,
    fix: str,
    severity: str,
) -> dict[str, object]:
    return {
        "category": category,
        "evidence": evidence,
        "affected_E_over_N_Td": eover,
        "affected_solver": solver,
        "suspected_code_area": area,
        "recommended_fix": fix,
        "severity": severity,
    }


def _rate_metric(metrics: dict[str, float], rate_name: str) -> float | str:
    value = metrics.get("major_rate_relative_difference", "")
    return value if rate_name else value


def _angular_status(reference_case: ReferenceCaseResult, candidate) -> tuple[str, object, object, str]:
    ref_model = reference_case.metadata.get("angular_model", "unknown")
    cand_model = candidate.metadata.get("angular_model", "unknown")
    if ref_model in {"", None, "unknown"} or cand_model in {"", None, "unknown"}:
        return "unknown", ref_model, cand_model, "angular_model_unknown"
    if ref_model == cand_model:
        return "true", ref_model, cand_model, ""
    return "false", ref_model, cand_model, f"{ref_model}!={cand_model}"


def _mc_confidence_status(
    reference_case: ReferenceCaseResult,
    candidate,
    metrics: dict[str, float],
) -> tuple[str, float | str]:
    ci = reference_case.metadata.get("scalar_ci95")
    if not isinstance(ci, dict) or not ci:
        return "unknown", ""
    covered = 0
    total = 0
    for scalar, attr in (
        ("mean_energy_eV", "mean_energy_eV"),
        ("drift_velocity_m_s", "drift_velocity_m_s"),
        ("mobility_m2_V_s", "mobility_m2_V_s"),
        ("diffusion_L_m2_s", "diffusion_L_m2_s"),
        ("diffusion_T_m2_s", "diffusion_T_m2_s"),
    ):
        if scalar not in ci or scalar not in reference_case.scalars:
            continue
        total += 1
        if abs(float(getattr(candidate, attr)) - reference_case.scalars[scalar]) <= float(ci[scalar]):
            covered += 1
    if total == 0:
        return "unknown", ""
    coverage = covered / total
    return ("pass_within_mc_uncertainty" if covered == total else "outside_mc_uncertainty"), coverage


def _classify_bolsig_mismatch(
    metrics: dict[str, float],
    *,
    candidate_solver: str,
    eover: float,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    if metrics.get("normalization_error_reference", 0.0) > 1.0e-8:
        rows.append(
            _failure(
                "normalization_mismatch",
                f"normalization_error_reference={metrics['normalization_error_reference']:.3e}",
                eover,
                candidate_solver,
                "external BOLSIG+ EEDF ingest",
                "check whether the BOLSIG+ file is EEDF or EEPF and normalize F(E)",
                "high",
            )
        )
    if metrics.get("eedf_relative_l1", 0.0) > BOLSIG_THRESHOLDS["eedf_relative_l1"]:
        rows.append(
            _failure(
                "energy_grid_mismatch",
                f"eedf_relative_l1={metrics['eedf_relative_l1']:.3e}",
                eover,
                candidate_solver,
                "EEDF interpolation/grid alignment",
                "compare BOLSIG+ grid, EEDF convention, and high-energy boundary",
                "medium",
            )
        )
    if metrics.get("mean_energy_relative_difference", 0.0) > BOLSIG_THRESHOLDS["mean_energy_relative_difference"]:
        rows.append(
            _failure(
                "bolsig_convention_mismatch",
                f"mean_energy_relative_difference={metrics['mean_energy_relative_difference']:.3e}",
                eover,
                candidate_solver,
                "transport scalar convention",
                "verify BOLSIG+ units and EEDF/EEPF conversion",
                "medium",
            )
        )
    if metrics.get("major_rate_relative_difference", 0.0) > BOLSIG_THRESHOLDS["major_rate_relative_difference"]:
        rows.append(
            _failure(
                "rate_convolution_mismatch",
                f"major_rate_relative_difference={metrics['major_rate_relative_difference']:.3e}",
                eover,
                candidate_solver,
                "rate convolution/source-sink",
                "compare threshold interpolation, source/sink, and rate weights",
                "medium",
            )
        )
    if metrics.get("log_tail_error", 0.0) > 2.0:
        rows.append(
            _failure(
                "high_energy_boundary_mismatch",
                f"log_tail_error={metrics['log_tail_error']:.3e}",
                eover,
                candidate_solver,
                "high-energy tail/boundary",
                "compare energy cutoff and high_energy_extrapolation behavior",
                "medium",
            )
        )
    return rows


def _classify_mcig_mismatch(
    metrics: dict[str, float],
    *,
    candidate_solver: str,
    candidate_method: object,
    reference_case: ReferenceCaseResult,
    eover: float,
    same_angular_model: str,
    confidence_status: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    if confidence_status == "pass_within_mc_uncertainty":
        rows.append(
            _failure(
                "MCIG_statistical_uncertainty",
                "difference is within reported MCIG confidence interval",
                eover,
                f"mcig->{candidate_solver}",
                "MCIG uncertainty",
                "treat this row as pass_within_mc_uncertainty",
                "low",
            )
        )
        return rows
    if same_angular_model != "true":
        rows.append(
            _failure(
                "angular_model_mismatch",
                f"same_angular_model={same_angular_model}",
                eover,
                f"mcig->{candidate_solver}",
                "angular scattering metadata",
                "compare only same-angular rows before assigning implementation mismatch",
                "medium",
            )
        )
        return rows
    if metrics.get("normalization_error_reference", 0.0) > 1.0e-8:
        rows.append(
            _failure(
                "EEDF_convention_mismatch",
                f"normalization_error_reference={metrics['normalization_error_reference']:.3e}",
                eover,
                f"mcig->{candidate_solver}",
                "MCIG EEDF ingest",
                "verify MCIG EEDF/EEPF convention and normalization",
                "high",
            )
        )
    if metrics.get("eedf_relative_l1", 0.0) > MCIG_THRESHOLDS["eedf_relative_l1"]:
        category = (
            "BOLSIG_two_term_approximation_difference"
            if candidate_solver == "two_term"
            else "multi_term_higher_l_model_issue"
            if candidate_solver == "multi_term" and str(candidate_method) == "pn_closure_direct"
            else "tail_sampling_uncertainty"
        )
        rows.append(
            _failure(
                category,
                f"eedf_relative_l1={metrics['eedf_relative_l1']:.3e}",
                eover,
                f"mcig->{candidate_solver}",
                "EEDF/tail comparison",
                "check MCIG confidence intervals, angular model, and high-energy tail sampling",
                "medium",
            )
        )
    if metrics.get("major_rate_relative_difference", 0.0) > MCIG_THRESHOLDS["major_rate_relative_difference"]:
        rows.append(
            _failure(
                "rate_convolution_mismatch",
                f"major_rate_relative_difference={metrics['major_rate_relative_difference']:.3e}",
                eover,
                f"mcig->{candidate_solver}",
                "rate convolution",
                "compare process labels, thresholds, and MCIG rate definitions",
                "medium",
            )
        )
    if metrics.get("mean_energy_relative_difference", 0.0) > MCIG_THRESHOLDS["mean_energy_relative_difference"] and candidate_solver == "two_term":
        rows.append(
            _failure(
                "code_two_term_projection_mismatch",
                f"mean_energy_relative_difference={metrics['mean_energy_relative_difference']:.3e}",
                eover,
                f"mcig->{candidate_solver}",
                "two_term projection or two-term approximation",
                "compare BOLSIG+ vs MCIG first to separate approximation error from code error",
                "medium",
            )
        )
    if metrics.get("log_tail_error", 0.0) > 2.0:
        rows.append(
            _failure(
                "tail_sampling_uncertainty",
                f"log_tail_error={metrics['log_tail_error']:.3e}",
                eover,
                f"mcig->{candidate_solver}",
                "MCIG tail statistics",
                "increase MCIG samples or inspect tail confidence intervals",
                "medium",
            )
        )
    if not rows and reference_case.metadata.get("uncertainty_unavailable", False):
        rows.append(
            _failure(
                "MCIG_statistical_uncertainty",
                "MCIG uncertainty unavailable; confidence is unknown",
                eover,
                f"mcig->{candidate_solver}",
                "MCIG uncertainty metadata",
                "provide MCIG confidence intervals to avoid overclaiming agreement",
                "low",
            )
        )
    return rows


def _write_plot(
    out_dir: Path,
    base_name: str,
    reference_cases: list[ReferenceCaseResult],
    solver_cases,
) -> Path | None:
    if not reference_cases:
        return None
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - depends on optional extra
        raise RuntimeError("--plot requires matplotlib; install the plot extra") from exc

    count = len(reference_cases)
    cols = 2 if count > 1 else 1
    rows = (count + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(6.0 * cols, 4.5 * rows), squeeze=False)
    for axis, reference_case in zip(axes.ravel(), reference_cases, strict=False):
        axis.semilogy(
            reference_case.energy_eV,
            reference_case.eedf_eV_inv,
            label=f"{reference_case.reference_id}",
            linewidth=2.0,
        )
        for case in solver_cases:
            if abs(case.e_over_n_Td - reference_case.e_over_n_Td) > 1.0e-9:
                continue
            method = case.metadata.get("solver_method", case.solver)
            lmax = case.metadata.get("lmax", "")
            suffix = f" lmax={lmax}" if lmax not in {"", None} else ""
            axis.semilogy(
                case.energy_eV,
                case.eedf,
                label=f"{case.solver} {method}{suffix}",
                alpha=0.85,
            )
        axis.set_title(f"E/N={reference_case.e_over_n_Td:g} Td")
        axis.set_xlabel("Energy [eV]")
        axis.set_ylabel("EEDF F(E) [1/eV]")
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(fontsize=8)
    for axis in axes.ravel()[count:]:
        axis.axis("off")
    fig.tight_layout()
    path = out_dir / f"{base_name}_eedf.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def run_benchmark(
    config_path: Path,
    *,
    reference: str = "bolsig",
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
    angular_model: str | None = None,
    fail_on_mismatch: bool = False,
    plot: bool = False,
) -> tuple[Path, ...]:
    cfg = load_config(config_path)
    out_dir = cfg.output.directory
    base_name = cfg.output.base_name
    summary_path = out_dir / f"{base_name}_summary.csv"
    metrics_path = out_dir / f"{base_name}_eedf_metrics.csv"
    failures_path = out_dir / f"{base_name}_failure_analysis.csv"
    report_path = out_dir / f"{base_name}_report.md"

    failure_rows: list[dict[str, object]] = []
    run_requested_external_reference_commands(
        cfg,
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
    reference_cases, reference_failures = _matching_references(
        cfg,
        reference=reference,
        bolsig_output=bolsig_output,
        mcig_output=mcig_output,
        require_bolsig=require_bolsig,
        require_mcig=require_mcig,
        angular_model=angular_model,
    )
    failure_rows.extend(reference_failures)

    solver_cases, solver_failures = run_solver_variants(cfg)
    failure_rows.extend(solver_failures)

    summary_rows: list[dict[str, object]] = []
    metric_rows: list[dict[str, object]] = []
    for reference_case in reference_cases:
        if reference_case.reference_id not in {_reference_key(reference), "bolsig_plus"} and _reference_key(reference) != "all":
            continue
        for candidate in [
            case
            for case in solver_cases
            if abs(case.e_over_n_Td - reference_case.e_over_n_Td) < 1.0e-9
        ]:
            metrics = reference_comparison_metrics(reference_case, candidate)
            same_angular_model, ref_angular, cand_angular, angular_reason = _angular_status(
                reference_case,
                candidate,
            )
            confidence_status, confidence_coverage = (
                _mc_confidence_status(reference_case, candidate, metrics)
                if reference_case.reference_id == "mcig"
                else ("not_applicable", "")
            )
            status = "PASS"
            thresholds = MCIG_THRESHOLDS if reference_case.reference_id == "mcig" else BOLSIG_THRESHOLDS
            for threshold_name, threshold in thresholds.items():
                if metrics.get(threshold_name, 0.0) > threshold:
                    status = "FAIL"
            if reference_case.reference_id == "mcig" and confidence_status == "pass_within_mc_uncertainty":
                status = "PASS_WITHIN_MC_UNCERTAINTY"
            if reference_case.reference_id == "mcig" and same_angular_model != "true":
                status = "DEGRADED"
            common = {
                "reference_id": reference_case.reference_id,
                "reference_case_id": reference_case.case_id,
                "candidate_solver": candidate.solver,
                "candidate_method": candidate.metadata.get("solver_method", ""),
                "candidate_lmax": candidate.metadata.get("lmax", ""),
                "E_over_N_Td": candidate.e_over_n_Td,
                "status": status,
                "confidence_status": confidence_status,
                "same_angular_model": same_angular_model,
                "angular_model_reference": ref_angular,
                "angular_model_candidate": cand_angular,
                "angular_model_mismatch_reason": angular_reason,
                "mc_confidence_interval_coverage": confidence_coverage,
            }
            summary_rows.append(
                {
                    **common,
                    "mean_energy_relative_difference": metrics.get("mean_energy_relative_difference", ""),
                    "drift_velocity_relative_difference": metrics.get("drift_velocity_relative_difference", ""),
                    "mobility_relative_difference": metrics.get("mobility_relative_difference", ""),
                    "diffusion_L_relative_difference": metrics.get("diffusion_L_relative_difference", ""),
                    "diffusion_T_relative_difference": metrics.get("diffusion_T_relative_difference", ""),
                    "ionization_rate_relative_difference": _rate_metric(metrics, "ionization"),
                    "excitation_rate_relative_difference": _rate_metric(metrics, "excitation"),
                    "eedf_relative_l1": metrics.get("eedf_relative_l1", ""),
                    "log_tail_error": metrics.get("log_tail_error", ""),
                    "tail_probability_difference": metrics.get("tail_probability_difference", ""),
                }
            )
            metric_rows.append({**common, **metrics})
            if reference_case.reference_id == "mcig":
                failure_rows.extend(
                    _classify_mcig_mismatch(
                        metrics,
                        candidate_solver=candidate.solver,
                        candidate_method=candidate.metadata.get("solver_method", ""),
                        reference_case=reference_case,
                        eover=candidate.e_over_n_Td,
                        same_angular_model=same_angular_model,
                        confidence_status=confidence_status,
                    )
                )
            else:
                failure_rows.extend(
                    _classify_bolsig_mismatch(
                        metrics,
                        candidate_solver=candidate.solver,
                        eover=candidate.e_over_n_Td,
                    )
                )

    if _reference_key(reference) in {"all", "mcig"}:
        mcig_cases = [case for case in reference_cases if case.reference_id == "mcig"]
        bolsig_cases = [case for case in reference_cases if case.reference_id == "bolsig_plus"]
        for mcig_case in mcig_cases:
            for bolsig_case in bolsig_cases:
                if abs(mcig_case.e_over_n_Td - bolsig_case.e_over_n_Td) > 1.0e-9:
                    continue
                from electron_swarm.references.common import reference_to_swarm_case

                bolsig_as_case = reference_to_swarm_case(bolsig_case)
                metrics = reference_comparison_metrics(mcig_case, bolsig_as_case)
                same_angular_model, ref_angular, cand_angular, angular_reason = _angular_status(
                    mcig_case,
                    bolsig_as_case,
                )
                common = {
                    "reference_id": "mcig",
                    "reference_case_id": mcig_case.case_id,
                    "candidate_solver": "reference:bolsig_plus",
                    "candidate_method": "bolsig_plus",
                    "candidate_lmax": "",
                    "E_over_N_Td": mcig_case.e_over_n_Td,
                    "status": "DEGRADED" if same_angular_model != "true" else "PASS",
                    "confidence_status": "unknown" if mcig_case.metadata.get("uncertainty_unavailable", False) else "reported",
                    "same_angular_model": same_angular_model,
                    "angular_model_reference": ref_angular,
                    "angular_model_candidate": cand_angular,
                    "angular_model_mismatch_reason": angular_reason,
                    "mc_confidence_interval_coverage": "",
                }
                summary_rows.append(
                    {
                        **common,
                        "mean_energy_relative_difference": metrics.get("mean_energy_relative_difference", ""),
                        "drift_velocity_relative_difference": metrics.get("drift_velocity_relative_difference", ""),
                        "mobility_relative_difference": metrics.get("mobility_relative_difference", ""),
                        "diffusion_L_relative_difference": metrics.get("diffusion_L_relative_difference", ""),
                        "diffusion_T_relative_difference": metrics.get("diffusion_T_relative_difference", ""),
                        "major_rate_relative_difference": metrics.get("major_rate_relative_difference", ""),
                        "ionization_rate_relative_difference": _rate_metric(metrics, "ionization"),
                        "excitation_rate_relative_difference": _rate_metric(metrics, "excitation"),
                        "eedf_relative_l1": metrics.get("eedf_relative_l1", ""),
                        "log_tail_error": metrics.get("log_tail_error", ""),
                        "tail_probability_difference": metrics.get("tail_probability_difference", ""),
                    }
                )
                metric_rows.append({**common, **metrics})

    for two in [case for case in solver_cases if case.solver == "two_term"]:
        direct = next(
            (
                case
                for case in solver_cases
                if case.solver == "multi_term"
                and case.e_over_n_Td == two.e_over_n_Td
                and str(case.metadata.get("lmax", "")) == "1"
            ),
            None,
        )
        if direct is None:
            continue
        gate = compare_eedf_cases(two, direct).metrics
        for metric, threshold in DIRECT_LMAX1_THRESHOLDS.items():
            if gate.get(metric, 0.0) > threshold:
                failure_rows.append(
                    _failure(
                        "lmax1_reduction_failure",
                        f"{metric}={gate.get(metric, 0.0):.3e}",
                        two.e_over_n_Td,
                        "multi_term",
                        "multi_term direct PN lmax=1 gate",
                        "inspect shared projection/source-sink/rate convolution",
                        "high",
                    )
                )

    _write_csv(summary_path, summary_rows, SUMMARY_FIELDS)
    _write_csv(metrics_path, metric_rows, METRIC_FIELDS)
    _write_csv(failures_path, failure_rows, FAILURE_FIELDS)
    written_plot = _write_plot(out_dir, base_name, reference_cases, solver_cases) if plot else None
    status = "PASS"
    if [row for row in summary_rows if row.get("status") == "FAIL"]:
        status = "FAIL"
    elif [row for row in summary_rows if row.get("status") == "DEGRADED"]:
        status = "DEGRADED"
    if not reference_cases:
        status = "SKIP_EXTERNAL_REFERENCE"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    selected_thresholds = MCIG_THRESHOLDS if _reference_key(reference) == "mcig" else BOLSIG_THRESHOLDS
    report_path.write_text(
        "\n".join(
            [
                "# Ar External Reference Benchmark",
                "",
                f"Config: `{config_path}`",
                f"Reference selector: `{reference}`",
                f"Status: **{status}**",
                f"Reference cases: {len(reference_cases)}",
                f"Solver cases: {len(solver_cases)}",
                f"Comparison rows: {len(summary_rows)}",
                f"Failure rows: {len(failure_rows)}",
                f"Plot: `{written_plot}`" if written_plot is not None else "Plot: not requested or no reference cases",
                "",
                "Pass thresholds:",
                f"- mean energy relative difference < {selected_thresholds['mean_energy_relative_difference']:.0%}",
                f"- major rates relative difference < {selected_thresholds['major_rate_relative_difference']:.0%}",
                f"- EEDF relative L1 < {selected_thresholds['eedf_relative_l1']:.0%}",
                "",
                "For MCIG rows, angular-model mismatches are marked degraded and "
                "reported MC confidence intervals are used when present.",
                "The lmax=1 direct PN vs two_term strict gate is also checked.",
                "External BOLSIG+/MCIG binaries are not run by this benchmark; provide output files.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    if fail_on_mismatch and status == "FAIL":
        raise RuntimeError("Ar BOLSIG+ equivalence benchmark mismatch")
    paths: list[Path] = [summary_path, metrics_path, failures_path, report_path]
    if written_plot is not None:
        paths.append(written_plot)
    return tuple(paths)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/benchmarks/ar_bolsig_plus_equivalence.yaml"),
    )
    parser.add_argument("--reference", default="bolsig", choices=["bolsig", "bolsig_plus", "mcig", "all"])
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
    parser.add_argument("--angular-model", choices=["isotropic", "mcig_default", "unknown"], default=None)
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("--fail-on-mismatch", action="store_true")
    args = parser.parse_args()
    paths = run_benchmark(
        args.config,
        reference=args.reference,
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
        angular_model=args.angular_model,
        fail_on_mismatch=args.fail_on_mismatch,
        plot=args.plot,
    )
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
