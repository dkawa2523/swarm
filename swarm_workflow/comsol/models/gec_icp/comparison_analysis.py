"""Aggregate and plot provenance-bound GEC-ICP saved-solution exports.

The exact common-time rows are the only like-for-like transient comparison.
Each model's terminal row is retained as a separately labelled end-state
comparison because the saved solutions may end at different physical times.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import re
from typing import Mapping, Sequence

from swarm_workflow._io import write_csv, write_json

from .comparison_analysis_contracts import (
    DEFAULT_STATIONARITY_LIMITS,
    GecIcpComparisonAnalysisError,
    GecIcpComparisonAnalysisSummary,
    GecIcpComparisonCaseData,
    MAX_COMPARISON_CASES,
    METRIC_KEYS,
    MIN_COMPARISON_CASES,
)
from .comparison_analysis_data import load_comparison_case, sha256_file
from .comparison_analysis_render import (
    render_common_time_relative_change,
    render_time_history,
)


COMPARISON_COLUMNS = (
    "case_id",
    "label",
    "basis",
    "time_s",
    "is_baseline",
    "like_for_like_time",
    "terminal_stationarity_passed",
    "axisymmetric_volume_m3",
    *METRIC_KEYS,
    *(f"{key}_signed_change_from_baseline" for key in METRIC_KEYS),
)
STATIONARITY_COLUMNS = (
    "case_id",
    "label",
    "terminal_time_s",
    "window_start_time_s",
    "window_points",
    "metric",
    "maximum_adjacent_relative_change",
    "limit",
    "passed",
)


def generate_gec_icp_saved_solution_comparison(
    *,
    cases: Mapping[str, str | Path],
    output_directory: str | Path,
    baseline_case_id: str = "original",
    common_time_s: float = 1.0e-3,
    case_labels: Mapping[str, str] | None = None,
    stationarity_limits: Mapping[str, float] | None = None,
) -> GecIcpComparisonAnalysisSummary:
    """Validate, aggregate, and render two to four comparison cases."""

    limits = _validate_request(
        cases,
        output_directory,
        baseline_case_id,
        common_time_s,
        case_labels,
        stationarity_limits,
    )
    labels = dict(case_labels or {})
    loaded = [
        load_comparison_case(
            case_id,
            labels.get(case_id, _default_label(case_id)),
            directory,
            common_time_s=common_time_s,
            stationarity_limits=limits,
        )
        for case_id, directory in cases.items()
    ]
    loaded = [
        *[case for case in loaded if case.case_id == baseline_case_id],
        *[case for case in loaded if case.case_id != baseline_case_id],
    ]
    _validate_cross_case_contract(loaded)
    baseline = next(case for case in loaded if case.case_id == baseline_case_id)
    terminal_time_aligned = all(
        math.isclose(
            case.terminal.time_s,
            baseline.terminal.time_s,
            rel_tol=1.0e-10,
            abs_tol=1.0e-15,
        )
        for case in loaded
    )

    output = Path(output_directory).resolve()
    output.mkdir(parents=True, exist_ok=True)
    comparison_csv = output / "icp_comparison_summary.csv"
    stationarity_csv = output / "icp_stationarity_summary.csv"
    comparison_rows = _comparison_rows(
        loaded,
        baseline,
        terminal_time_aligned=terminal_time_aligned,
    )
    stationarity_rows = _stationarity_rows(loaded)
    write_csv(comparison_csv, COMPARISON_COLUMNS, comparison_rows)
    write_csv(stationarity_csv, STATIONARITY_COLUMNS, stationarity_rows)
    figures = (
        *render_time_history(output, loaded, common_time_s=common_time_s),
        *render_common_time_relative_change(
            output,
            loaded,
            baseline_case_id=baseline_case_id,
            common_time_s=common_time_s,
        ),
    )
    manifest_path = output / "icp_comparison_manifest.json"
    manifest = _manifest_payload(
        loaded,
        baseline,
        common_time_s=common_time_s,
        terminal_time_aligned=terminal_time_aligned,
        comparison_csv=comparison_csv,
        stationarity_csv=stationarity_csv,
        figures=figures,
    )
    write_json(manifest_path, manifest)
    return GecIcpComparisonAnalysisSummary(
        output_directory=output,
        comparison_csv=comparison_csv,
        stationarity_csv=stationarity_csv,
        figures=tuple(figures),
        manifest=manifest_path,
        cases=tuple(loaded),
    )


def _validate_request(
    cases: Mapping[str, str | Path],
    output_directory: str | Path,
    baseline_case_id: str,
    common_time_s: float,
    case_labels: Mapping[str, str] | None,
    stationarity_limits: Mapping[str, float] | None,
) -> dict[str, float]:
    if not isinstance(cases, Mapping) or not (
        MIN_COMPARISON_CASES <= len(cases) <= MAX_COMPARISON_CASES
    ):
        raise GecIcpComparisonAnalysisError(
            f"comparison requires {MIN_COMPARISON_CASES} to "
            f"{MAX_COMPARISON_CASES} cases"
        )
    for case_id in cases:
        if (
            not isinstance(case_id, str)
            or re.fullmatch(r"[a-z][a-z0-9_]*", case_id) is None
        ):
            raise GecIcpComparisonAnalysisError(
                "comparison case ids must match [a-z][a-z0-9_]*"
            )
    if baseline_case_id not in cases:
        raise GecIcpComparisonAnalysisError(
            f"baseline case is absent: {baseline_case_id}"
        )
    directories = [Path(path).resolve() for path in cases.values()]
    if len(set(directories)) != len(directories):
        raise GecIcpComparisonAnalysisError(
            "comparison cases must use distinct export directories"
        )
    output = Path(output_directory).resolve()
    if output in directories:
        raise GecIcpComparisonAnalysisError(
            "analysis output must be separate from every case export directory"
        )
    if not math.isfinite(common_time_s) or common_time_s < 0.0:
        raise GecIcpComparisonAnalysisError(
            "common comparison time must be finite and nonnegative"
        )
    if case_labels is not None:
        unknown = set(case_labels) - set(cases)
        if unknown:
            raise GecIcpComparisonAnalysisError(
                f"labels name unknown cases: {sorted(unknown)}"
            )
        if any(not str(value).strip() for value in case_labels.values()):
            raise GecIcpComparisonAnalysisError("case labels must be nonempty")
    limits = dict(DEFAULT_STATIONARITY_LIMITS)
    if stationarity_limits is not None:
        if set(stationarity_limits) != set(METRIC_KEYS):
            raise GecIcpComparisonAnalysisError(
                "stationarity limits must define every comparison metric exactly"
            )
        limits = {key: float(value) for key, value in stationarity_limits.items()}
    if any(not math.isfinite(value) or value <= 0.0 for value in limits.values()):
        raise GecIcpComparisonAnalysisError(
            "stationarity limits must be finite and positive"
        )
    return limits


def _validate_cross_case_contract(
    cases: Sequence[GecIcpComparisonCaseData],
) -> None:
    mapping_hashes = {case.model_mapping_sha256 for case in cases}
    if len(mapping_hashes) != 1:
        raise GecIcpComparisonAnalysisError(
            "comparison cases were not exported through the same model mapping"
        )
    reference_volume = cases[0].common.axisymmetric_volume_m3
    for case in cases:
        for state in (case.common, case.terminal):
            if not math.isclose(
                state.axisymmetric_volume_m3,
                reference_volume,
                rel_tol=1.0e-9,
                abs_tol=1.0e-15,
            ):
                raise GecIcpComparisonAnalysisError(
                    "comparison cases do not share one axisymmetric plasma volume"
                )


def _comparison_rows(
    cases: Sequence[GecIcpComparisonCaseData],
    baseline: GecIcpComparisonCaseData,
    *,
    terminal_time_aligned: bool,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for basis, baseline_state in (
        ("common_time_exact", baseline.common),
        ("terminal_model_specific", baseline.terminal),
    ):
        for case in cases:
            state = case.common if basis == "common_time_exact" else case.terminal
            row: dict[str, object] = {
                "case_id": case.case_id,
                "label": case.label,
                "basis": basis,
                "time_s": state.time_s,
                "is_baseline": int(case.case_id == baseline.case_id),
                "like_for_like_time": int(
                    basis == "common_time_exact" or terminal_time_aligned
                ),
                "terminal_stationarity_passed": int(case.converged),
                "axisymmetric_volume_m3": state.axisymmetric_volume_m3,
                **state.metrics(),
            }
            for key in METRIC_KEYS:
                reference = float(getattr(baseline_state, key))
                row[f"{key}_signed_change_from_baseline"] = (
                    float(getattr(state, key)) - reference
                ) / abs(reference)
            rows.append(row)
    return rows


def _stationarity_rows(
    cases: Sequence[GecIcpComparisonCaseData],
) -> list[dict[str, object]]:
    return [
        {
            "case_id": case.case_id,
            "label": case.label,
            "terminal_time_s": case.terminal.time_s,
            "window_start_time_s": case.time_series[-4].time_s,
            "window_points": 4,
            "metric": key,
            "maximum_adjacent_relative_change": case.stationarity_changes[key],
            "limit": case.stationarity_limits[key],
            "passed": int(case.stationarity_passed[key]),
        }
        for case in cases
        for key in METRIC_KEYS
    ]


def _manifest_payload(
    cases: Sequence[GecIcpComparisonCaseData],
    baseline: GecIcpComparisonCaseData,
    *,
    common_time_s: float,
    terminal_time_aligned: bool,
    comparison_csv: Path,
    stationarity_csv: Path,
    figures: Sequence[Path],
) -> dict[str, object]:
    failed_baseline_metrics = [
        key for key in METRIC_KEYS if not baseline.stationarity_passed[key]
    ]
    warnings: list[str] = []
    if failed_baseline_metrics:
        warnings.append(
            "baseline terminal state is not stationary under the final-four-point "
            "contract; failed metrics: " + ", ".join(failed_baseline_metrics)
        )
    if not terminal_time_aligned:
        warnings.append(
            "model-specific terminal states occur at different physical times and "
            "are not a like-for-like transient comparison"
        )
    outputs = (comparison_csv, stationarity_csv, *figures)
    return {
        "schema": "swarm.gec_icp_saved_solution_comparison.v1",
        "status": "completed",
        "baseline_case_id": baseline.case_id,
        "comparison_contract": {
            "common_time": {
                "time_s": common_time_s,
                "basis": "exact transient interpolation at one shared physical time",
                "like_for_like": True,
            },
            "terminal": {
                "basis": "each model's last saved state",
                "like_for_like_time": terminal_time_aligned,
                "terminal_times_s": {
                    case.case_id: case.terminal.time_s for case in cases
                },
            },
            "stationarity": {
                "window": "final four saved time points",
                "metric": "maximum symmetric adjacent relative change",
                "limits": dict(baseline.stationarity_limits),
            },
        },
        "cases": [
            {
                "case_id": case.case_id,
                "label": case.label,
                "terminal_time_s": case.terminal.time_s,
                "terminal_stationarity_passed": case.converged,
                "stationarity_changes": dict(case.stationarity_changes),
                "source": {
                    "export_directory": str(case.export_directory),
                    "comparison_export_manifest": {
                        "path": str(case.source_manifest),
                        "sha256": case.source_manifest_sha256,
                    },
                    "input_mph": {
                        "path": str(case.input_mph),
                        "sha256": case.input_mph_sha256,
                    },
                    "model_mapping": {
                        "path": str(case.model_mapping),
                        "sha256": case.model_mapping_sha256,
                    },
                    "csv_sha256": dict(case.source_hashes),
                },
            }
            for case in cases
        ],
        "assessment": {
            "baseline_terminal_stationarity_passed": baseline.converged,
            "baseline_failed_stationarity_metrics": failed_baseline_metrics,
            "warnings": warnings,
        },
        "outputs": [
            {
                "path": path.name,
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
            for path in outputs
        ],
    }


def _default_label(case_id: str) -> str:
    standard = {
        "original": "Original COMSOL closure",
        "two_term": "Two-term",
        "propagator": "Propagator P1",
        "composite": "MC + two-term fallback",
    }
    return standard.get(case_id, case_id.replace("_", " ").title())


def _parse_assignments(values: Sequence[str], name: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for value in values:
        key, separator, item = value.partition("=")
        if not separator or not key or not item or key in parsed:
            raise GecIcpComparisonAnalysisError(
                f"{name} entries must be unique ID=VALUE assignments"
            )
        parsed[key] = item
    return parsed


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case",
        action="append",
        required=True,
        metavar="ID=DIRECTORY",
        help="completed read-only export case; repeat two to four times",
    )
    parser.add_argument(
        "--label",
        action="append",
        default=[],
        metavar="ID=LABEL",
    )
    parser.add_argument("--baseline", default="original")
    parser.add_argument("--common-time-s", type=float, default=1.0e-3)
    parser.add_argument("--output-directory", required=True, type=Path)
    return parser


def main() -> None:
    args = _parser().parse_args()
    summary = generate_gec_icp_saved_solution_comparison(
        cases=_parse_assignments(args.case, "case"),
        output_directory=args.output_directory,
        baseline_case_id=args.baseline,
        common_time_s=args.common_time_s,
        case_labels=_parse_assignments(args.label, "label"),
    )
    print(summary.manifest)


if __name__ == "__main__":
    main()


__all__ = (
    "GecIcpComparisonAnalysisError",
    "GecIcpComparisonAnalysisSummary",
    "generate_gec_icp_saved_solution_comparison",
)
