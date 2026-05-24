"""Run the Ar product EEDF consistency benchmark."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from electron_swarm import load_config, run
from electron_swarm.diagnostics.eedf_compare import compare_eedf_cases


def _write_csv(
    path: Path,
    rows: list[dict[str, object]],
    *,
    fields: list[str] | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0]) if rows else ["status"]
    with path.open("w", encoding="utf-8", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run_benchmark(config_path: Path) -> tuple[Path, Path, Path]:
    cfg = load_config(config_path)
    result = run(cfg, write=False)
    grouped: dict[tuple[str, float], object] = {
        (case.solver, case.e_over_n_Td): case for case in result.cases
    }
    metric_rows: list[dict[str, object]] = []
    failure_rows: list[dict[str, object]] = []
    for case in result.cases:
        if case.solver == "two_term":
            continue
        reference = grouped.get(("two_term", case.e_over_n_Td))
        if reference is None:
            continue
        comparison = compare_eedf_cases(reference, case)  # type: ignore[arg-type]
        metric_rows.append(
            {
                "reference_solver": "two_term",
                "candidate_solver": case.solver,
                "E_over_N_Td": case.e_over_n_Td,
                **comparison.metrics,
            }
        )
        failure_rows.extend(comparison.failures)

    out_dir = cfg.output.directory
    metrics_path = out_dir / "ar_eedf_consistency_eedf_metrics.csv"
    failures_path = out_dir / "ar_eedf_consistency_failure_analysis.csv"
    report_path = out_dir / "ar_eedf_consistency_report.md"
    _write_csv(metrics_path, metric_rows)
    _write_csv(
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
    report_path.parent.mkdir(parents=True, exist_ok=True)
    status = "PASS" if not failure_rows else "FAIL"
    report_path.write_text(
        "\n".join(
            [
                "# Ar EEDF Consistency",
                "",
                f"Config: `{config_path}`",
                f"Status: **{status}**",
                f"Metric rows: {len(metric_rows)}",
                f"Failure rows: {len(failure_rows)}",
                "",
                "The comparison uses normalized EEDF `F(E)` in `1/eV`, not EEPF.",
                "`pn_closure_direct` is intentionally excluded until a true "
                "coupled f0/f1 block PN implementation exists.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return metrics_path, failures_path, report_path


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
