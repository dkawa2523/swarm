"""Reproducible static plots for the argon GEC CCP workflow."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from swarm_workflow.comsol.models.gec_ccp.plots import eedf as _eedf
from swarm_workflow.comsol.models.gec_ccp.plots import data as _data
from swarm_workflow.comsol.models.gec_ccp.plots import spatial as _spatial
from swarm_workflow.comsol.models.gec_ccp.plots.contracts import (
    CANONICAL_RESTRICTED_TRANSPORT_CLOSURE,
    COMSOL_DIFFERENCE_COLORS,
    COMSOL_RAINBOW_LIGHT_COLORS,
    CURRENT_SOLVER_COMPARISON_CLOSURES,
    FUNCTION_EEDF_REACTION_MODELS,
    GROUNDED_BOUNDARY_COLOR,
    GecCcpPlotError,
    GecCcpPlotSummary,
    HYBRID_FUNCTION_EEDF_REACTION_MODEL,
    POWERED_BOUNDARY_COLOR,
)
from swarm_workflow.comsol.models.gec_ccp.contracts import (
    GEC_TWO_TERM_HYBRID_TRANSPORT_CLOSURE,
)


__all__ = [
    "CANONICAL_RESTRICTED_TRANSPORT_CLOSURE",
    "COMSOL_DIFFERENCE_COLORS",
    "COMSOL_RAINBOW_LIGHT_COLORS",
    "CURRENT_SOLVER_COMPARISON_CLOSURES",
    "FUNCTION_EEDF_REACTION_MODELS",
    "GEC_TWO_TERM_HYBRID_TRANSPORT_CLOSURE",
    "GROUNDED_BOUNDARY_COLOR",
    "GecCcpPlotError",
    "GecCcpPlotSummary",
    "HYBRID_FUNCTION_EEDF_REACTION_MODEL",
    "POWERED_BOUNDARY_COLOR",
    "plot_gec_ccp_results",
    "plot_gec_ccp_solver_comparison",
]


def plot_gec_ccp_results(
    bundle_path: str | Path,
    *,
    output_dir: str | Path,
    comsol_results_dir: str | Path | None = None,
) -> GecCcpPlotSummary:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise GecCcpPlotError(
            "plot-gec-ccp requires the optional matplotlib dependency"
        ) from exc

    bundle = Path(bundle_path).resolve()
    bundle_manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    closure_source = str(bundle_manifest.get("source", "Swarm"))
    run_provenance: dict[str, Any] | None = None
    include_builtin_reference = False
    external_label = f"External {closure_source} closure"
    reaction_model: str | None = None
    function_eedf_table: str | None = None
    if comsol_results_dir is not None:
        run_provenance = _data._validated_plot_provenance(
            bundle,
            bundle_manifest,
            Path(comsol_results_dir).resolve(),
        )
        external_label = run_provenance["external_label"]
        include_builtin_reference = run_provenance["include_builtin_reference"]
        reaction_model = run_provenance["closure"]["reaction_model"]
        function_summary = run_provenance["bundle"].get("function_eedf")
        if isinstance(function_summary, dict):
            function_eedf_table = str(function_summary.get("table") or "") or None
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    figures: list[Path] = []
    skipped: list[str] = []
    comparison_metrics: dict[str, Any] = {}

    figures.extend(
        _eedf.plot_closure_overview(
            bundle,
            output,
            plt,
            closure_source=closure_source,
            reaction_model=reaction_model,
        )
    )
    if comsol_results_dir is None:
        skipped.append("COMSOL comparison: no result directory supplied")
    else:
        results = Path(comsol_results_dir).resolve()
        spatial_metrics, spatial_figures, spatial_skipped = (
            _plot_comsol_result_comparisons(
                results,
                output,
                plt,
                include_builtin_reference=include_builtin_reference,
                external_label=external_label,
            )
        )
        comparison_metrics.update(spatial_metrics)
        figures.extend(spatial_figures)
        skipped.extend(spatial_skipped)
        operating_metrics, operating_figures, operating_skipped = (
            _eedf.plot_operating_eedf_if_available(
                bundle,
                results / "swarm_tables",
                output,
                plt,
                reaction_model=reaction_model,
                function_eedf_table=function_eedf_table,
            )
        )
        comparison_metrics.update(operating_metrics)
        figures.extend(operating_figures)
        skipped.extend(operating_skipped)

    manifest = output / "plot_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "stage": "plot-gec-ccp",
                "status": "ok",
                "bundle": str(bundle),
                "run_provenance": run_provenance,
                "figures": [str(path) for path in figures],
                "skipped": skipped,
                "comparison_labels": {
                    **(
                        {"baseline": "COMSOL built-in Druyvesteyn"}
                        if include_builtin_reference
                        else {}
                    ),
                    "external": external_label,
                },
                "comparison_metrics": comparison_metrics,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return GecCcpPlotSummary(output, tuple(figures), manifest, tuple(skipped))


def _plot_comsol_result_comparisons(
    results: Path,
    output: Path,
    plt: Any,
    *,
    include_builtin_reference: bool,
    external_label: str,
) -> tuple[dict[str, Any], list[Path], list[str]]:
    """Compose the result-specific renderers owned by the spatial module."""

    baseline = results / "builtin_druyvesteyn"
    external = results / "swarm_tables"
    figures: list[Path] = []
    skipped: list[str] = []
    metrics: dict[str, Any] = {}
    external_domain = external / "domain_period_average.csv"
    if external_domain.exists():
        figures.extend(
            _spatial._plot_gec_geometry_2d(
                external_domain,
                output / "gec_ccp_geometry_2d",
                plt,
            )
        )
    else:
        skipped.append("2D-domain COMSOL result CSV is absent")
    domain_files = (baseline / "domain_period_average.csv", external_domain)
    if include_builtin_reference and all(path.exists() for path in domain_files):
        domain_summary, domain_figures = _spatial._plot_domain_2d_comparison(
            domain_files[0],
            domain_files[1],
            output / "domain_2d_period_average_comparison",
            plt,
            external_label=external_label,
        )
        difference_summary, difference_figures = _spatial._plot_domain_2d_difference(
            domain_files[0],
            domain_files[1],
            output / "domain_2d_period_average_difference",
            plt,
            external_label=external_label,
        )
        domain_summary["difference_maps"] = difference_summary
        metrics["domain_2d_period_average"] = domain_summary
        figures.extend((*domain_figures, *difference_figures))
    elif include_builtin_reference:
        skipped.append("2D-domain built-in-reference comparison was not available")
    axis_files = (
        baseline / "axis_period_average.csv",
        external / "axis_period_average.csv",
    )
    radial_files = (
        baseline / "radial_period_average.csv",
        external / "radial_period_average.csv",
    )
    for name, paths, skip_message in (
        (
            "axis_period_average",
            axis_files,
            "axis built-in-reference comparison was not available",
        ),
        (
            "radial_period_average",
            radial_files,
            "radial built-in-reference comparison was not available",
        ),
    ):
        if include_builtin_reference and all(path.exists() for path in paths):
            metrics[name] = _data._spatial_comparison_metrics(paths[0], paths[1])
        elif include_builtin_reference:
            skipped.append(skip_message)
    if include_builtin_reference and all(
        path.exists() for path in (*axis_files, *radial_files)
    ):
        figures.extend(
            _spatial._plot_spatial_cut_overview(
                axis_files,
                radial_files,
                output / "period_averaged_spatial_distributions",
                plt,
                external_label=external_label,
            )
        )
    phase_path = external / "axis_phase_resolved.csv"
    if phase_path.exists():
        phase_summary, phase_figures = _spatial._plot_phase_axial_distributions(
            phase_path,
            output / "rf_phase_axial_distributions",
            plt,
            external_label=external_label,
        )
        metrics["swarm_phase_resolved_axial"] = phase_summary
        figures.extend(phase_figures)
    else:
        skipped.append("phase-resolved axial COMSOL result CSV is absent")
    waveform_files = (
        baseline / "electrode_waveform.csv",
        external / "electrode_waveform.csv",
    )
    if include_builtin_reference and all(path.exists() for path in waveform_files):
        figures.extend(
            _spatial._plot_waveform_comparison(
                waveform_files[0],
                waveform_files[1],
                output / "electrode_waveform_comparison",
                plt,
                external_label=external_label,
            )
        )
        metrics["electrode_waveform"] = _data._waveform_comparison_metrics(
            waveform_files[0], waveform_files[1]
        )
    elif include_builtin_reference:
        skipped.append("waveform built-in-reference comparison was not available")
    return metrics, figures, skipped


def plot_gec_ccp_solver_comparison(
    two_term_bundle_path: str | Path,
    monte_carlo_bundle_path: str | Path,
    *,
    two_term_results_dir: str | Path,
    monte_carlo_results_dir: str | Path,
    output_dir: str | Path,
) -> GecCcpPlotSummary:
    """Plot accepted two-term and Monte Carlo physical targets together.

    Each result must pass its own physics audits and originate from the same
    input MPH.  A repeated COMSOL built-in reference is not part of this
    cross-solver comparison.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise GecCcpPlotError(
            "plot-gec-ccp-comparison requires the optional matplotlib dependency"
        ) from exc

    cases = {
        "two_term": (
            Path(two_term_bundle_path).resolve(),
            Path(two_term_results_dir).resolve(),
        ),
        "monte_carlo": (
            Path(monte_carlo_bundle_path).resolve(),
            Path(monte_carlo_results_dir).resolve(),
        ),
    }
    provenance: dict[str, dict[str, Any]] = {}
    for solver, (bundle, results) in cases.items():
        try:
            bundle_manifest = json.loads(
                (bundle / "manifest.json").read_text(encoding="utf-8")
            )
        except (OSError, json.JSONDecodeError) as exc:
            raise GecCcpPlotError(f"{solver} bundle manifest is unreadable") from exc
        case_provenance = _data._validated_plot_provenance(
            bundle, bundle_manifest, results
        )
        if case_provenance["bundle"].get("source") != solver:
            raise GecCcpPlotError(
                f"comparison slot {solver} contains "
                f"{case_provenance['bundle'].get('source')}"
            )
        if (
            case_provenance["planned_result_role"] != "physical_target"
            or case_provenance["result_role"] != "physical_target"
            or case_provenance["physical_target_accepted"] is not True
        ):
            raise GecCcpPlotError(
                f"{solver} comparison requires an accepted physical-target result"
            )
        closure = case_provenance["closure"]
        expected_closure = CURRENT_SOLVER_COMPARISON_CLOSURES[solver]
        if closure != expected_closure:
            raise GecCcpPlotError(
                f"{solver} comparison requires the current source-specific "
                f"closure: {expected_closure}"
            )
        provenance[solver] = case_provenance

    input_hashes = {
        solver: _data._comparison_input_mph_sha256(case_provenance)
        for solver, case_provenance in provenance.items()
    }
    if len(set(input_hashes.values())) != 1:
        raise GecCcpPlotError(
            "two_term and monte_carlo results do not share the same original MPH"
        )

    two_term_results = cases["two_term"][1]
    monte_carlo_results = cases["monte_carlo"][1]
    external_paths = (
        two_term_results / "swarm_tables" / "domain_period_average.csv",
        monte_carlo_results / "swarm_tables" / "domain_period_average.csv",
    )
    external_radial_paths = (
        two_term_results / "swarm_tables" / "radial_period_average.csv",
        monte_carlo_results / "swarm_tables" / "radial_period_average.csv",
    )
    required = (*external_paths, *external_radial_paths)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise GecCcpPlotError(
            "solver comparison requires all 2D domain and radial-midplane "
            "exports: " + ", ".join(missing)
        )
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    two_term_label = "Two-term Swarm\nFunction-EEDF + Einstein diffusion"
    monte_carlo_label = "Monte Carlo Swarm\nDirect rates + Einstein diffusion"
    metrics, figures = _spatial._plot_domain_2d_comparison(
        external_paths[0],
        external_paths[1],
        output / "domain_2d_solver_comparison",
        plt,
        baseline_label=two_term_label,
        external_label=monte_carlo_label,
        identify_electrodes=True,
        include_ionization_source=False,
        include_svg=True,
    )
    metrics["reference_solver"] = "two_term"
    metrics["comparison_solver"] = "monte_carlo"
    domain_coordinates, _ = _data._comsol_domain_data(external_paths[0])
    geometry = _spatial._gec_domain_geometry(domain_coordinates)
    radial_metrics, radial_figures = _spatial._plot_radial_midplane_solver_comparison(
        (
            external_radial_paths[0],
            external_radial_paths[1],
        ),
        (
            "two_term Swarm",
            "monte_carlo Swarm",
        ),
        output / "radial_midplane_solver_comparison",
        plt,
        geometry=geometry,
        case_ids=("two_term", "monte_carlo"),
        include_ionization_source=False,
        include_svg=True,
    )
    metrics["radial_midplane"] = radial_metrics
    figures.extend(radial_figures)

    manifest = output / "solver_comparison_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "stage": "plot-gec-ccp-comparison",
                "status": "ok",
                "acceptance_policy": (
                    "both solver-specific physical-target runs passed "
                    "physics quality; source-specific Function-EEDF inputs "
                    "share swarm_mobility_einstein "
                    "restricted-LMEA transport"
                ),
                "input_mph_sha256": next(iter(input_hashes.values())),
                "run_provenance": provenance,
                "comparison_labels": {
                    "two_term": provenance["two_term"]["external_label"],
                    "monte_carlo": provenance["monte_carlo"]["external_label"],
                },
                "figures": [str(path) for path in figures],
                "comparison_metrics": metrics,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return GecCcpPlotSummary(output, tuple(figures), manifest, ())
