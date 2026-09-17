"""Reproducible comparison of independently calculated solver EEDFs.

All distributions use the probability-density convention ``F(E)`` in
``eV^-1`` with ``integral F(E) dE = 1``. Comparisons conservatively rebin
probability mass onto the union of the supplied finite-volume cell edges;
point interpolation is deliberately not used for distribution metrics.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Mapping, Sequence

from swarm_workflow.plots.eedf_contracts import (
    DEFAULT_REPRESENTATIVE_FIELDS_TD,
    DEFAULT_TAIL_THRESHOLDS_EV,
    EedfCase,
    EedfComparisonError,
    EedfComparisonSummary,
    EedfDataset,
)
from swarm_workflow.plots.eedf_data import (
    load_mc_eedf_database,
    load_table_eedf,
    qualification_by_field,
    validate_compatible_physical_sources,
    validate_mc_qualification_manifest,
)
from swarm_workflow.plots.eedf_provenance import (
    sha256,
    validate_manifest_artifact,
)
from swarm_workflow.plots.eedf_metrics import compare_eedf_cases
from swarm_workflow.plots.eedf_render import plot_distributions, plot_metrics

__all__ = [
    "EedfCase",
    "EedfComparisonError",
    "EedfComparisonSummary",
    "EedfDataset",
    "compare_eedf_cases",
    "generate_eedf_solver_comparison",
    "load_mc_eedf_database",
    "load_table_eedf",
]


def _write_metrics(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        raise EedfComparisonError("no EEDF comparison metrics were generated")
    fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _metric_rows(
    fields: Sequence[float],
    two_term: EedfDataset,
    propagator: EedfDataset,
    monte_carlo: EedfDataset,
    tail_thresholds_eV: Sequence[float],
) -> list[dict[str, object]]:
    comparisons: list[
        tuple[str, Mapping[float, EedfCase], Mapping[float, EedfCase]]
    ] = [
        ("raw_mc_vs_two_term", monte_carlo.cases, two_term.cases),
        ("propagator_vs_two_term", propagator.cases, two_term.cases),
        ("raw_mc_vs_propagator", monte_carlo.cases, propagator.cases),
    ]
    rows: list[dict[str, object]] = []
    for field in fields:
        for name, candidates, references in comparisons:
            rows.append(
                {
                    "E_over_N_Td": field,
                    "comparison": name,
                    **compare_eedf_cases(
                        references[field],
                        candidates[field],
                        tail_thresholds_eV=tail_thresholds_eV,
                    ),
                }
            )
    return rows


def _write_mc_status(
    path: Path,
    fields: Sequence[float],
    statuses: Mapping[float, bool],
    qualification_column: str,
) -> None:
    """Record MC closure status without inferring downstream source selection."""

    not_passed = [field for field in fields if not statuses[field]]
    payload = {
        "format_version": 1,
        "scope": "MC active-closure status annotation only",
        "comparison_invariant": (
            "EEDF figures and metrics always use the independent raw output "
            "from each solver"
        ),
        "selection_policy": (
            "not inferred here; any substitution must be recorded by a separate "
            "downstream source-selection artifact"
        ),
        "qualification_column": qualification_column,
        "per_anchor": [
            {
                "E_over_N_Td": field,
                "mc_active_closure_gate": (
                    "passed" if statuses[field] else "not_passed"
                ),
            }
            for field in fields
        ],
        "passed_fields_Td": [field for field in fields if statuses[field]],
        "not_passed_fields_Td": not_passed,
    }
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def generate_eedf_solver_comparison(
    *,
    two_term_directory: str | Path,
    propagator_directory: str | Path,
    mc_database: str | Path,
    mc_qualification_csv: str | Path,
    output_directory: str | Path,
    mixture_id: int = 0,
    qualification_column: str = "active_closure_quality_passed",
    representative_fields_td: Sequence[float] = DEFAULT_REPRESENTATIVE_FIELDS_TD,
    tail_thresholds_eV: Sequence[float] = DEFAULT_TAIL_THRESHOLDS_EV,
) -> EedfComparisonSummary:
    """Compare independent raw 2T, P1, and MC EEDFs.

    The supplied MC qualification remains visible as downstream status context,
    but it never replaces a Monte Carlo curve or changes a pairwise metric.
    """

    output = Path(output_directory).resolve()
    qualification_path = Path(mc_qualification_csv).resolve()
    mc_manifest_path = qualification_path.parent / "manifest.json"
    two_term = load_table_eedf(
        two_term_directory,
        "two_term",
        require_manifest=True,
    )
    propagator = load_table_eedf(
        propagator_directory,
        "propagator",
        require_manifest=True,
    )
    monte_carlo = load_mc_eedf_database(
        mc_database,
        mixture_id=mixture_id,
        source_manifest_path=mc_manifest_path,
    )
    validate_manifest_artifact(
        mc_manifest_path,
        qualification_path,
        solver="monte_carlo",
    )
    validate_mc_qualification_manifest(mc_manifest_path, qualification_path)
    validate_compatible_physical_sources(two_term, propagator, monte_carlo)
    statuses = qualification_by_field(qualification_path, qualification_column)
    if set(statuses) != set(monte_carlo.cases):
        raise EedfComparisonError(
            "Monte Carlo qualification and aggregate EEDF anchor sets do not match"
        )
    common_fields = sorted(
        set(monte_carlo.cases) & set(two_term.cases) & set(propagator.cases)
    )
    if len(common_fields) < 2:
        raise EedfComparisonError(
            "EEDF comparison requires at least two common physical E/N anchors"
        )
    status_not_passed_fields = {
        field for field in common_fields if not statuses[field]
    }
    representatives = tuple(float(value) for value in representative_fields_td)
    if len(set(representatives)) != len(representatives):
        raise EedfComparisonError("representative E/N anchors must be unique")
    for field in representatives:
        if field not in common_fields:
            raise EedfComparisonError(
                f"representative E/N lacks common evidence: {field:g} Td"
            )

    output.mkdir(parents=True, exist_ok=True)

    metric_rows = _metric_rows(
        common_fields,
        two_term,
        propagator,
        monte_carlo,
        tail_thresholds_eV,
    )
    metrics_path = output / "eedf_comparison_metrics.csv"
    _write_metrics(metrics_path, metric_rows)
    status_path = output / "eedf_mc_status.json"
    _write_mc_status(
        status_path,
        common_fields,
        statuses,
        qualification_column,
    )
    figures = (
        *plot_distributions(
            output,
            representatives,
            statuses,
            monte_carlo,
            two_term,
            propagator,
            tail_thresholds_eV,
        ),
        *plot_metrics(output, metric_rows),
    )
    cases = [
        *two_term.cases.values(),
        *propagator.cases.values(),
        *monte_carlo.cases.values(),
    ]
    comparison_scope = {
        "mode": "independent_raw_solver_outputs",
        "solvers": ["two_term", "propagator", "monte_carlo"],
        "common_comparison_fields_Td": common_fields,
        "source_specific_support_fields_Td": {
            "two_term": sorted(set(two_term.cases) - set(common_fields)),
            "propagator": sorted(set(propagator.cases) - set(common_fields)),
            "monte_carlo": sorted(set(monte_carlo.cases) - set(common_fields)),
        },
        "metric_pairs": [
            "raw_mc_vs_two_term",
            "propagator_vs_two_term",
            "raw_mc_vs_propagator",
        ],
    }

    uncertainty = {
        "mc_band": "pointwise 95% Student-t interval across independent replicas",
        "scope": (
            "not a simultaneous confidence region and not a significance "
            "test against deterministic solvers"
        ),
    }
    qa = {
        "maximum_source_normalization_error": max(
            abs(case.source_normalization - 1.0) for case in cases
        ),
        "maximum_reported_vs_reconstructed_mean_relative_difference": max(
            abs(case.reported_mean_energy_eV - case.reconstructed_mean_energy_eV)
            / case.reported_mean_energy_eV
            for case in cases
        ),
        "metrics_use_only_raw_solver_pairs": {
            str(row["comparison"]) for row in metric_rows
        }
        == {
            "raw_mc_vs_two_term",
            "propagator_vs_two_term",
            "raw_mc_vs_propagator",
        },
        "mc_active_closure_gate_passed_anchors": sum(
            int(statuses[field]) for field in common_fields
        ),
    }

    manifest = {
        "format_version": 2,
        "comparison": "raw_two_term_propagator_and_raw_monte_carlo",
        "convention": {
            "quantity": "F(E) probability-density EEDF",
            "units": "eV^-1",
            "normalization": "integral F(E) dE = 1",
            "common_grid": (
                "union of exact cell edges with conservative probability-mass rebinning"
            ),
            "plot_energy_range": (
                "field-specific; shared by all solver curves within each panel"
            ),
            "physical_distance": "total variation on physical electron energy",
            "shape_distance": (
                "total variation after each EEDF is scaled by its reconstructed "
                "mean energy"
            ),
        },
        "comparison_scope": comparison_scope,
        "mc_status_annotation": {
            "qualification_column": qualification_column,
            "passed_fields_Td": [
                field for field in common_fields if statuses[field]
            ],
            "not_passed_fields_Td": sorted(status_not_passed_fields),
            "meaning": (
                "downstream active-closure status; not an EEDF curve replacement"
            ),
        },
        "mc_status": {
            "artifact": status_path.name,
            "scope": "annotation only; separate from solver comparison",
        },
        "representative_fields_Td": list(representatives),
        "tail_thresholds_eV": [float(value) for value in tail_thresholds_eV],
        "uncertainty": uncertainty,
        "sources": {
            "two_term": dict(two_term.provenance),
            "propagator": dict(propagator.provenance),
            "monte_carlo": dict(monte_carlo.provenance),
            "monte_carlo_qualification": {
                "path": str(qualification_path),
                "sha256": sha256(qualification_path),
            },
        },
        "qa": qa,
        "outputs": {
            "figures": [path.name for path in figures],
            "metrics_csv": metrics_path.name,
            "mc_status": status_path.name,
        },
    }
    manifest_path = output / "eedf_comparison_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return EedfComparisonSummary(
        output_directory=output,
        figures=tuple(figures),
        metrics_csv=metrics_path,
        manifest=manifest_path,
        mc_status=status_path,
        mc_status_not_passed_fields_td=tuple(sorted(status_not_passed_fields)),
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--two-term-directory", required=True, type=Path)
    parser.add_argument("--propagator-directory", required=True, type=Path)
    parser.add_argument("--mc-database", required=True, type=Path)
    parser.add_argument("--mc-qualification-csv", required=True, type=Path)
    parser.add_argument("--output-directory", required=True, type=Path)
    parser.add_argument("--mixture-id", type=int, default=0)
    parser.add_argument(
        "--qualification-column",
        default="active_closure_quality_passed",
    )
    parser.add_argument(
        "--representative-fields-td",
        type=float,
        nargs="+",
        default=list(DEFAULT_REPRESENTATIVE_FIELDS_TD),
    )
    parser.add_argument(
        "--tail-thresholds-eV",
        type=float,
        nargs="*",
        default=list(DEFAULT_TAIL_THRESHOLDS_EV),
    )
    return parser


def main() -> None:
    args = _parser().parse_args()
    summary = generate_eedf_solver_comparison(
        two_term_directory=args.two_term_directory,
        propagator_directory=args.propagator_directory,
        mc_database=args.mc_database,
        mc_qualification_csv=args.mc_qualification_csv,
        output_directory=args.output_directory,
        mixture_id=args.mixture_id,
        qualification_column=args.qualification_column,
        representative_fields_td=args.representative_fields_td,
        tail_thresholds_eV=args.tail_thresholds_eV,
    )
    print(summary.manifest)


if __name__ == "__main__":
    main()
