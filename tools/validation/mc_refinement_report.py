"""Metrics, plots, and evidence reports for MC population refinement."""

from __future__ import annotations

import csv
from hashlib import sha256
import itertools
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from swarm_workflow.plots.eedf_metrics import (
    compare_eedf_cases,
    mass_on,
    tail_probability,
    union_edges,
)

from tools.validation.mc_refinement_contracts import (
    DEFAULT_POPULATIONS,
    DatabaseInput,
    HIGH_FIELD_TD,
    LOW_FIELD_TAIL_ANCHORS_TD,
    LoadedDatabase,
    QualificationStatus,
    RefinementEvaluationError,
    SCHEMA,
    TAIL_THRESHOLDS_EV,
)
from tools.validation.mc_refinement_database import (
    file_sha256,
    load_database,
    load_qualifications,
    validate_qualification_manifest,
)


def _tail_summaries(
    tail_rows: Sequence[dict[str, object]],
) -> list[dict[str, object]]:
    grouped: dict[tuple[str, int, float], list[dict[str, object]]] = {}
    for row in tail_rows:
        field = float(row["E_over_N_Td"])
        if field in LOW_FIELD_TAIL_ANCHORS_TD:
            grouped.setdefault(
                (str(row["mixture"]), int(row["database_population"]), field), []
            ).append(row)
    summaries: list[dict[str, object]] = []
    for (mixture, population, field), rows in sorted(grouped.items()):
        executed = sum(int(row["tail_executed"]) for row in rows)
        configured = sum(int(row["tail_configured"]) for row in rows)
        summaries.append(
            {
                "mixture": mixture,
                "database_population": population,
                "E_over_N_Td": field,
                "replicas": len(rows),
                "tail_configured_replicas": configured,
                "tail_executed_replicas": executed,
                "tail_execution_fraction": executed / len(rows),
            }
        )
    return summaries


def _replica_stability(
    loaded: Mapping[tuple[str, int], LoadedDatabase],
    *,
    final_population: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    pair_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    mixtures = sorted({key[0] for key in loaded})
    for mixture in mixtures:
        database = loaded[(mixture, final_population)]
        for field, cases in sorted(database.final_replica_cases.items()):
            particle_counts = {
                int(row["particles"])
                for row in database.tail_rows
                if float(row["E_over_N_Td"]) == field
            }
            if len(particle_counts) != 1:
                raise RefinementEvaluationError(
                    f"{mixture} {field:g} Td replicas use inconsistent particles"
                )
            values: list[float] = []
            for (left_replica, left), (right_replica, right) in itertools.combinations(
                sorted(cases), 2
            ):
                metrics = compare_eedf_cases(left, right, tail_thresholds_eV=())
                value = float(metrics["total_variation"])
                values.append(value)
                pair_rows.append(
                    {
                        "mixture": mixture,
                        "database_population": final_population,
                        "E_over_N_Td": field,
                        "replica_a": left_replica,
                        "replica_b": right_replica,
                        "total_variation": value,
                    }
                )
            if not values:
                raise RefinementEvaluationError(
                    f"{mixture} {field:g} Td final database needs at least two replicas"
                )
            summary_rows.append(
                {
                    "mixture": mixture,
                    "database_population": final_population,
                    "E_over_N_Td": field,
                    "anchor_particles": next(iter(particle_counts)),
                    "replicas": len(cases),
                    "pair_count": len(values),
                    "pairwise_total_variation_mean": float(np.mean(values)),
                    "pairwise_total_variation_max": float(np.max(values)),
                }
            )
    return pair_rows, summary_rows


def _relative_change(candidate: float, reference: float) -> float | None:
    if reference == 0.0:
        return None
    return (candidate - reference) / reference


def _population_refinement(
    loaded: Mapping[tuple[str, int], LoadedDatabase],
    *,
    populations: Sequence[int],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    mixtures = sorted({key[0] for key in loaded})
    for mixture in mixtures:
        for reference_population, candidate_population in itertools.pairwise(populations):
            reference = loaded[(mixture, reference_population)].aggregate_high_field
            candidate = loaded[(mixture, candidate_population)].aggregate_high_field
            metrics = compare_eedf_cases(
                reference,
                candidate,
                tail_thresholds_eV=TAIL_THRESHOLDS_EV,
            )
            edges = union_edges((reference, candidate))
            reference_mass = mass_on(reference, edges)
            candidate_mass = mass_on(candidate, edges)
            reference_mean = reference.reported_mean_energy_eV
            candidate_mean = candidate.reported_mean_energy_eV
            reference_reconstructed_mean = reference.reconstructed_mean_energy_eV
            candidate_reconstructed_mean = candidate.reconstructed_mean_energy_eV
            row: dict[str, object] = {
                "mixture": mixture,
                "E_over_N_Td": HIGH_FIELD_TD,
                "reference_population": reference_population,
                "candidate_population": candidate_population,
                "total_variation": float(metrics["total_variation"]),
                "reference_mean_energy_eV": reference_mean,
                "candidate_mean_energy_eV": candidate_mean,
                "mean_energy_relative_change": _relative_change(
                    candidate_mean, reference_mean
                ),
                "mean_energy_absolute_relative_difference": abs(
                    candidate_mean - reference_mean
                )
                / reference_mean,
                "reference_reconstructed_mean_energy_eV": (
                    reference_reconstructed_mean
                ),
                "candidate_reconstructed_mean_energy_eV": (
                    candidate_reconstructed_mean
                ),
                "reconstructed_mean_energy_absolute_relative_difference": abs(
                    candidate_reconstructed_mean - reference_reconstructed_mean
                )
                / reference_reconstructed_mean,
            }
            for threshold in TAIL_THRESHOLDS_EV:
                label = f"{threshold:g}"
                reference_tail = tail_probability(edges, reference_mass, threshold)
                candidate_tail = tail_probability(edges, candidate_mass, threshold)
                row[f"reference_tail_mass_gt_{label}_eV"] = reference_tail
                row[f"candidate_tail_mass_gt_{label}_eV"] = candidate_tail
                row[f"tail_mass_relative_change_gt_{label}_eV"] = _relative_change(
                    candidate_tail, reference_tail
                )
            rows.append(row)
    return rows


def _qualification_rows(
    qualifications: Mapping[str, Mapping[float, QualificationStatus]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for mixture, by_field in sorted(qualifications.items()):
        for field, status in sorted(by_field.items()):
            rows.append(
                {
                    "mixture": mixture,
                    "E_over_N_Td": field,
                    **status.csv_values(),
                }
            )
    return rows


def _replica_refinement_rows(
    loaded: Mapping[tuple[str, int], LoadedDatabase],
    *,
    mixtures: Sequence[str],
    populations: Sequence[int],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for mixture in mixtures:
        previous: LoadedDatabase | None = None
        for population in populations:
            database = loaded[(mixture, population)]
            if previous is not None:
                for name, value in database.campaign_limits.items():
                    if value < previous.campaign_limits[name]:
                        raise RefinementEvaluationError(
                            f"{mixture} MC campaign {name} decreases from "
                            f"p{previous.source.population} to p{population}"
                        )
            for field in sorted(database.replicas_by_anchor):
                replicas = database.replicas_by_anchor[field]
                previous_replicas = (
                    previous.replicas_by_anchor[field] if previous is not None else None
                )
                if previous_replicas is not None and replicas < previous_replicas:
                    raise RefinementEvaluationError(
                        f"{mixture} {field:g} Td replicas decrease from "
                        f"{previous_replicas} to {replicas} at p{population}"
                    )
                rows.append(
                    {
                        "mixture": mixture,
                        "database_population": population,
                        "E_over_N_Td": field,
                        "previous_database_population": (
                            "" if previous is None else previous.source.population
                        ),
                        "previous_replicas": (
                            "" if previous_replicas is None else previous_replicas
                        ),
                        "replicas": replicas,
                        "replica_increment": (
                            0 if previous_replicas is None else replicas - previous_replicas
                        ),
                        "final_population": int(population == populations[-1]),
                        "campaign_maximum_replicas": database.campaign_limits[
                            "maximum_replicas"
                        ],
                        "campaign_maximum_total_particle_barriers": (
                            database.campaign_limits[
                                "maximum_total_particle_barriers"
                            ]
                        ),
                    }
                )
            previous = database
    return rows


def _write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        raise RefinementEvaluationError(f"refusing to write empty evidence: {path}")
    fieldnames = list(rows[0])
    if any(list(row) != fieldnames for row in rows):
        raise RefinementEvaluationError(f"inconsistent CSV columns for {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _pretty_mixture(mixture: str) -> str:
    return mixture.upper().replace("_", "/").replace("CL2", "Cl2")


def _plot_raw_high_field(
    path_png: Path,
    path_svg: Path,
    loaded: Mapping[tuple[str, int], LoadedDatabase],
    *,
    mixtures: Sequence[str],
    populations: Sequence[int],
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RefinementEvaluationError(
            "matplotlib is required for MC refinement figures"
        ) from exc

    colors = {512: "#2563EB", 1024: "#D97706", 2048: "#4D7C0F"}
    fig, axes = plt.subplots(
        1,
        len(mixtures),
        figsize=(5.2 * len(mixtures), 4.6),
        squeeze=False,
        sharey=True,
    )
    for column, mixture in enumerate(mixtures):
        axis = axes[0, column]
        positive: list[float] = []
        for population in populations:
            case = loaded[(mixture, population)].aggregate_high_field
            axis.stairs(
                case.density_eV_inv,
                case.edges_eV,
                baseline=None,
                linewidth=1.5,
                color=colors.get(population),
                label=f"p{population}",
            )
            visible = case.density_eV_inv[
                (case.centers_eV <= 200.0) & (case.density_eV_inv > 0.0)
            ]
            positive.extend(visible.tolist())
        for threshold in TAIL_THRESHOLDS_EV:
            axis.axvline(threshold, color="#6B7280", linewidth=0.7, linestyle=":")
        axis.set_yscale("log")
        axis.set_xlim(0.0, 200.0)
        if positive:
            maximum = max(positive)
            minimum = min(positive)
            axis.set_ylim(max(minimum * 0.5, maximum * 1.0e-12), maximum * 2.0)
        axis.grid(True, which="both", linewidth=0.35, alpha=0.35)
        axis.set_title(_pretty_mixture(mixture))
        axis.set_xlabel("Electron energy [eV]")
        axis.legend(frameon=False)
    axes[0, 0].set_ylabel("Raw aggregate MC EEDF [eV$^{-1}$]")
    fig.suptitle("1000 Td Monte Carlo population refinement (no smoothing or fallback)")
    fig.tight_layout()
    fig.savefig(path_png, dpi=220, bbox_inches="tight")
    fig.savefig(path_svg, format="svg", bbox_inches="tight", metadata={"Date": None})
    plt.close(fig)


def _json_compatible(value: object) -> object:
    if isinstance(value, Mapping):
        return {str(key): _json_compatible(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_compatible(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        raise RefinementEvaluationError("non-finite value in evaluation output")
    return value


def evaluate_mc_refinement(
    database_inputs: Sequence[DatabaseInput],
    *,
    output_directory: str | Path,
    qualification_paths: Mapping[str, Path] | None = None,
    mixture_id: int = 0,
    populations: Sequence[int] = DEFAULT_POPULATIONS,
    overwrite: bool = False,
) -> dict[str, Path]:
    populations = tuple(int(value) for value in populations)
    if len(populations) < 2 or len(set(populations)) != len(populations):
        raise RefinementEvaluationError("populations must be distinct and ordered")
    if tuple(sorted(populations)) != populations:
        raise RefinementEvaluationError("populations must be strictly increasing")
    by_key: dict[tuple[str, int], DatabaseInput] = {}
    for item in database_inputs:
        key = (item.mixture, int(item.population))
        if key in by_key:
            raise RefinementEvaluationError(f"duplicate database input for {key}")
        by_key[key] = DatabaseInput(item.mixture, int(item.population), item.path)
    mixtures = tuple(sorted({key[0] for key in by_key}))
    if not mixtures:
        raise RefinementEvaluationError("at least one mixture is required")
    expected = {
        (mixture, population) for mixture in mixtures for population in populations
    }
    if set(by_key) != expected:
        missing = sorted(expected - set(by_key))
        extra = sorted(set(by_key) - expected)
        raise RefinementEvaluationError(
            f"database matrix is incomplete (missing={missing}, extra={extra})"
        )

    output = Path(output_directory).resolve()
    paths = {
        "tail_execution": output / "tail_execution.csv",
        "tail_summary": output / "tail_execution_summary.csv",
        "replica_pairs": output / "replica_pairwise_tv.csv",
        "replica_stability": output / "replica_stability.csv",
        "population_refinement": output / "population_refinement.csv",
        "qualification": output / "qualification_status.csv",
        "replica_refinement": output / "replica_refinement.csv",
        "evaluation": output / "mc_refinement_evaluation.json",
        "figure_png": output / "mc_refinement_1000Td.png",
        "figure_svg": output / "mc_refinement_1000Td.svg",
        "manifest": output / "manifest.json",
    }
    existing = [path for path in paths.values() if path.exists()]
    if existing and not overwrite:
        raise RefinementEvaluationError(
            f"refusing to overwrite existing outputs: {[str(path) for path in existing]}"
        )
    output.mkdir(parents=True, exist_ok=True)

    final_population = populations[-1]
    loaded: dict[tuple[str, int], LoadedDatabase] = {}
    for key in sorted(by_key):
        loaded[key] = load_database(
            by_key[key],
            mixture_id=mixture_id,
            load_replicas=key[1] == final_population,
        )
    for mixture in mixtures:
        reference_database = loaded[(mixture, populations[0])]
        reference_fractions = reference_database.fractions
        for population in populations[1:]:
            candidate_database = loaded[(mixture, population)]
            if candidate_database.fractions != reference_fractions:
                raise RefinementEvaluationError(
                    f"{mixture} mixture fractions change across refinement databases"
                )
            if (
                candidate_database.refinement_context_sha256
                != reference_database.refinement_context_sha256
            ):
                raise RefinementEvaluationError(
                    f"{mixture} physical, cross-section, estimator, seed, campaign, "
                    "or quality context changes across refinement databases"
                )
            if (
                candidate_database.sampling_contract_sha256
                != reference_database.sampling_contract_sha256
            ):
                raise RefinementEvaluationError(
                    f"{mixture} sampling plans differ by more than the 1000 Td "
                    "particle population and replica counts"
                )

    replica_refinement_rows = _replica_refinement_rows(
        loaded, mixtures=mixtures, populations=populations
    )

    qualifications: dict[str, dict[float, QualificationStatus]] = {}
    qualification_hashes: dict[str, dict[str, str]] = {}
    for mixture, raw_path in sorted((qualification_paths or {}).items()):
        if mixture not in mixtures:
            raise RefinementEvaluationError(
                f"qualification supplied for unknown mixture {mixture!r}"
            )
        path = Path(raw_path).resolve()
        if not path.is_file():
            raise RefinementEvaluationError(f"missing qualification CSV: {path}")
        statuses = load_qualifications(path)
        csv_hash = file_sha256(path)
        qualifications[mixture] = statuses
        qualification_hashes[mixture] = validate_qualification_manifest(
            path,
            csv_sha256=csv_hash,
            statuses=statuses,
            final_database=loaded[(mixture, final_population)],
        )

    tail_rows = [row for database in loaded.values() for row in database.tail_rows]
    tail_rows.sort(
        key=lambda row: (
            str(row["mixture"]),
            int(row["database_population"]),
            float(row["E_over_N_Td"]),
            int(row["replicate"]),
        )
    )
    tail_summary = _tail_summaries(tail_rows)
    pair_rows, stability_rows = _replica_stability(
        loaded,
        final_population=final_population,
    )
    refinement_rows = _population_refinement(loaded, populations=populations)
    qualification_rows = _qualification_rows(qualifications)

    _write_csv(paths["tail_execution"], tail_rows)
    _write_csv(paths["tail_summary"], tail_summary)
    _write_csv(paths["replica_pairs"], pair_rows)
    _write_csv(paths["replica_stability"], stability_rows)
    _write_csv(paths["population_refinement"], refinement_rows)
    _write_csv(paths["replica_refinement"], replica_refinement_rows)
    if qualification_rows:
        _write_csv(paths["qualification"], qualification_rows)
    else:
        with paths["qualification"].open("w", encoding="utf-8", newline="") as handle:
            handle.write(
                "mixture,E_over_N_Td,qualification_profile,"
                "qualification_overall_passed,qualification_aggregate_quality_passed,"
                "qualification_active_closure_passed,"
                "qualification_solver_transport_passed,"
                "qualification_explicit_eedf_passed,"
                "qualification_eedf_normalization_passed\n"
            )
    _plot_raw_high_field(
        paths["figure_png"],
        paths["figure_svg"],
        loaded,
        mixtures=mixtures,
        populations=populations,
    )

    low_rows = [
        row
        for row in tail_rows
        if float(row["E_over_N_Td"]) in LOW_FIELD_TAIL_ANCHORS_TD
    ]
    final_low_rows = [
        row
        for row in low_rows
        if int(row["database_population"]) == final_population
    ]
    unique_low_rows = list(
        {
            (
                str(row["mixture"]),
                float(row["E_over_N_Td"]),
                int(row["particles"]),
                int(row["seed"]),
                int(row["tail_configured_collisions"]),
            ): row
            for row in low_rows
        }.values()
    )
    final_low_executed = sum(int(row["tail_executed"]) for row in final_low_rows)
    unique_low_executed = sum(int(row["tail_executed"]) for row in unique_low_rows)
    evaluation = {
        "schema": SCHEMA,
        "comparison_scope": {
            "solver": "monte_carlo",
            "mode": "raw_population_and_replica_refinement",
            "high_field_E_over_N_Td": HIGH_FIELD_TD,
            "populations": list(populations),
            "tail_thresholds_eV": list(TAIL_THRESHOLDS_EV),
            "fallback_used": False,
            "smoothing_used": False,
            "distance_grid": "finite_volume_conservative_union",
        },
        "qualification_semantics": {
            "separation": (
                "qualification is reported in its own table and JSON section; "
                "it is never inserted into raw EEDF distance rows"
            ),
            "explicit_eedf_passed": (
                "reported only when the qualification CSV has a dedicated EEDF "
                "quality column; otherwise null"
            ),
            "active_closure_passed": (
                "downstream function-EEDF closure status; not relabeled as an "
                "EEDF-only or transport-only result"
            ),
            "solver_transport_passed": (
                "copied only from solver_transport_qualified; not inferred from "
                "EEDF TV or normalization"
            ),
            "eedf_normalization_passed": (
                "normalization gate only; it is not a shape-convergence gate"
            ),
        },
        "qualification_status": qualification_rows,
        "replica_refinement": replica_refinement_rows,
        "low_field_tail_execution": {
            "anchors_Td": list(LOW_FIELD_TAIL_ANCHORS_TD),
            "final_database": {
                "database_population": final_population,
                "executed_replicas": final_low_executed,
                "total_replicas": len(final_low_rows),
                "execution_fraction": final_low_executed / len(final_low_rows),
            },
            "unique_physical_cases": {
                "deduplication_key": [
                    "mixture",
                    "E_over_N_Td",
                    "particles",
                    "seed",
                    "tail_configured_collisions",
                ],
                "executed_replicas": unique_low_executed,
                "total_replicas": len(unique_low_rows),
                "execution_fraction": unique_low_executed / len(unique_low_rows),
            },
            "database_case_rows": {
                "note": "includes exact replicas reused into later population databases",
                "executed_replicas": sum(int(row["tail_executed"]) for row in low_rows),
                "total_replicas": len(low_rows),
            },
        },
        "tail_execution_summary": tail_summary,
        "final_population_replica_stability": stability_rows,
        "high_field_population_refinement": refinement_rows,
        "representational_projection_events": [
            event for database in loaded.values() for event in database.projection_events
        ],
    }
    paths["evaluation"].write_text(
        json.dumps(_json_compatible(evaluation), indent=2, sort_keys=True, allow_nan=False)
        + "\n",
        encoding="utf-8",
    )

    for database in loaded.values():
        if file_sha256(database.source.path) != database.source_sha256:
            raise RefinementEvaluationError(
                f"source database changed during evaluation: {database.source.path}"
            )
    for mixture, evidence in qualification_hashes.items():
        if file_sha256(Path(evidence["path"])) != evidence["sha256"]:
            raise RefinementEvaluationError(
                f"qualification CSV changed during evaluation: {mixture}"
            )
        if file_sha256(Path(evidence["manifest_path"])) != evidence["manifest_sha256"]:
            raise RefinementEvaluationError(
                f"qualification manifest changed during evaluation: {mixture}"
            )
    evaluator_path = Path(__file__).with_name("evaluate_mc_refinement.py").resolve()
    manifest = {
        "schema": SCHEMA,
        "inputs": [
            {
                "mixture": database.source.mixture,
                "population": database.source.population,
                "path": str(database.source.path),
                "sha256": database.source_sha256,
                "mixture_fractions": database.fractions,
                "refinement_context_sha256": database.refinement_context_sha256,
                "sampling_contract_excluding_1000Td_population_sha256": (
                    database.sampling_contract_sha256
                ),
                "workflow_config_sha256": database.metadata["workflow_config_sha256"],
                "mc_sampling_plan_sha256": sha256(
                    database.metadata["mc_sampling_plan_json"].encode("utf-8")
                ).hexdigest(),
                "mc_transport_estimator_schema": database.metadata[
                    "mc_transport_estimator_schema_version"
                ],
                "mc_eedf_estimator_schema": database.metadata[
                    "mc_eedf_estimator_schema_version"
                ],
                "mc_tail_estimator_schema": database.metadata[
                    "mc_tail_estimator_schema_version"
                ],
                "replicas_by_anchor": {
                    f"{field:g}": replicas
                    for field, replicas in sorted(
                        database.replicas_by_anchor.items()
                    )
                },
                "campaign_limits": database.campaign_limits,
            }
            for database in (loaded[key] for key in sorted(loaded))
        ],
        "qualification_inputs": qualification_hashes,
        "replica_refinement": replica_refinement_rows,
        "evaluator": {
            "path": str(evaluator_path),
            "sha256": file_sha256(evaluator_path),
        },
        "invariants": {
            "raw_monte_carlo_only": True,
            "two_term_fallback_used": False,
            "histogram_smoothing_used": False,
            "source_databases_modified": False,
        },
        "representational_projection": {
            "scope": "in-memory machine-resolution cells only",
            "events": evaluation["representational_projection_events"],
        },
        "artifacts": {
            name: {"path": str(path), "sha256": file_sha256(path)}
            for name, path in paths.items()
            if name != "manifest"
        },
    }
    paths["manifest"].write_text(
        json.dumps(_json_compatible(manifest), indent=2, sort_keys=True, allow_nan=False)
        + "\n",
        encoding="utf-8",
    )
    return paths


__all__ = ["evaluate_mc_refinement"]
