"""Read-only SQLite snapshots and provenance for MC refinement evidence."""

from __future__ import annotations

import json
import math
from hashlib import sha256
from pathlib import Path
import sqlite3
from typing import Mapping, Sequence

import numpy as np

from electron_swarm.solvers.monte_carlo.evidence import (
    MC_EEDF_ESTIMATOR_SCHEMA_VERSION,
)

from swarm_workflow.campaign.aggregate import stable_mc_seed
from swarm_workflow.campaign.statistics import energy_edges_from_cells
from swarm_workflow.campaign.store import (
    WorkflowSchemaError,
    mc_sampling_plan_provenance,
)
from swarm_workflow.plots.eedf_contracts import EedfCase, EedfComparisonError
from swarm_workflow.plots.eedf_data import csv_rows, finite_float, validated_case
from swarm_workflow.plots.eedf_metrics import mass_on, union_edges

from tools.validation.mc_refinement_contracts import (
    DatabaseInput,
    HIGH_FIELD_TD,
    LoadedDatabase,
    MAX_AGGREGATE_MEAN_RELATIVE_ERROR,
    MAX_AGGREGATE_RECONSTRUCTION_TV,
    MAX_REPRESENTATIONAL_MERGE_MASS,
    MAX_REPRESENTATIONAL_MOMENT_CHANGE_EV,
    QualificationStatus,
    REFINEMENT_CONTEXT_METADATA,
    REQUIRED_ANCHORS_TD,
    RefinementEvaluationError,
)


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _parse_bool(row: Mapping[str, str], name: str) -> bool | None:
    raw = str(row.get(name, "")).strip().lower()
    if not raw:
        return None
    if raw in {"1", "true"}:
        return True
    if raw in {"0", "false"}:
        return False
    raise RefinementEvaluationError(
        f"qualification column {name!r} must be boolean when present"
    )


def canonical_anchor(value: float) -> float:
    for anchor in REQUIRED_ANCHORS_TD:
        if math.isclose(value, anchor, rel_tol=1.0e-12, abs_tol=1.0e-12):
            return anchor
    return float(value)


def require_mapping(value: object, label: str) -> dict[str, object]:
    if not isinstance(value, dict):
        raise RefinementEvaluationError(f"{label} must be an object")
    return value


def require_int(value: object, label: str, *, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise RefinementEvaluationError(f"{label} must be an integer")
    result = int(value)
    if float(value) != float(result) or result < minimum:
        raise RefinementEvaluationError(
            f"{label} must be an integer greater than or equal to {minimum}"
        )
    return result


def _read_metadata(connection: sqlite3.Connection) -> dict[str, str]:
    try:
        rows = connection.execute("SELECT key, value FROM metadata").fetchall()
    except sqlite3.Error as exc:
        raise RefinementEvaluationError(
            "database metadata table is unavailable"
        ) from exc
    return {str(row[0]): str(row[1]) for row in rows}


def _refinement_provenance(
    metadata: Mapping[str, str],
    *,
    declared_population: int,
) -> tuple[
    str,
    str,
    dict[float, dict[str, object]],
    dict[float, int],
    dict[str, int],
]:
    missing = [
        name
        for name in (*REFINEMENT_CONTEXT_METADATA, "mc_campaign_json")
        if not metadata.get(name)
    ]
    if missing:
        raise RefinementEvaluationError(
            f"MC database lacks refinement provenance {missing}"
        )
    if (
        metadata["mc_eedf_estimator_schema_version"]
        != MC_EEDF_ESTIMATOR_SCHEMA_VERSION
    ):
        raise RefinementEvaluationError(
            "MC database uses a different finite-volume EEDF estimator schema"
        )
    try:
        campaign = require_mapping(
            json.loads(metadata["mc_campaign_json"]), "MC campaign"
        )
        policy = require_mapping(campaign.get("policy"), "MC campaign policy")
    except json.JSONDecodeError as exc:
        raise RefinementEvaluationError("MC campaign JSON is invalid") from exc
    campaign_limits = {
        name: require_int(policy.get(name), f"MC campaign policy {name}", minimum=1)
        for name in ("maximum_replicas", "maximum_total_particle_barriers")
    }
    normalized_campaign = dict(campaign)
    normalized_policy = dict(policy)
    for name in campaign_limits:
        normalized_policy[name] = "replica_refinement_variable"
    normalized_campaign["policy"] = normalized_policy
    context = {
        name: metadata[name] for name in REFINEMENT_CONTEXT_METADATA
    }
    context["mc_campaign_json"] = json.dumps(
        normalized_campaign, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    context_json = json.dumps(context, sort_keys=True, separators=(",", ":"))
    try:
        provenance = mc_sampling_plan_provenance(dict(metadata), required=True)
    except WorkflowSchemaError as exc:
        raise RefinementEvaluationError(str(exc)) from exc
    assert provenance is not None
    entries = provenance["entries"]
    if not isinstance(entries, list):
        raise RefinementEvaluationError("MC sampling plan entries must be a list")
    normalized: list[dict[str, object]] = []
    observed_fields: set[float] = set()
    plan_by_field: dict[float, dict[str, object]] = {}
    for raw in entries:
        if not isinstance(raw, dict):
            raise RefinementEvaluationError("MC sampling plan row must be an object")
        row = dict(raw)
        field = canonical_anchor(float(row["e_over_n_Td"]))
        if field in plan_by_field:
            raise RefinementEvaluationError(
                f"MC sampling plan has duplicate {field:g} Td rows"
            )
        observed_fields.add(field)
        plan_by_field[field] = dict(row)
        if field == HIGH_FIELD_TD:
            if int(row["particles"]) != declared_population:
                raise RefinementEvaluationError(
                    f"1000 Td sampling plan records {row['particles']} particles, "
                    f"not declared p{declared_population}"
                )
            row["particles"] = "population_refinement_variable"
        row["replicas"] = "replica_refinement_variable"
        normalized.append(row)
    if observed_fields != set(REQUIRED_ANCHORS_TD):
        raise RefinementEvaluationError(
            "MC sampling plan must contain exactly 1, 10, 100, and 1000 Td"
        )
    sampling_contract = json.dumps(
        normalized, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    replicas_by_anchor = {
        field: require_int(
            plan.get("replicas"), f"{field:g} Td planned replicas", minimum=2
        )
        for field, plan in plan_by_field.items()
    }
    if campaign_limits["maximum_replicas"] < max(replicas_by_anchor.values()):
        raise RefinementEvaluationError(
            "MC campaign maximum_replicas is below the sampling-plan replica count"
        )
    return (
        sha256(context_json.encode("utf-8")).hexdigest(),
        sha256(sampling_contract.encode("utf-8")).hexdigest(),
        plan_by_field,
        replicas_by_anchor,
        campaign_limits,
    )


def _read_fractions(
    connection: sqlite3.Connection, mixture_id: int
) -> dict[str, float]:
    try:
        row = connection.execute(
            "SELECT fractions_json FROM mixtures WHERE mixture_id = ?",
            (int(mixture_id),),
        ).fetchone()
    except sqlite3.Error as exc:
        raise RefinementEvaluationError("database mixtures table is unavailable") from exc
    if row is None:
        raise RefinementEvaluationError(f"mixture_id {mixture_id} is absent")
    try:
        values = json.loads(str(row[0]))
    except json.JSONDecodeError as exc:
        raise RefinementEvaluationError("mixture fractions JSON is invalid") from exc
    if not isinstance(values, dict) or not values:
        raise RefinementEvaluationError("mixture fractions must be a nonempty object")
    result = {str(name): float(fraction) for name, fraction in values.items()}
    if any(not math.isfinite(value) or value < 0.0 for value in result.values()):
        raise RefinementEvaluationError(
            "mixture fractions must be finite and nonnegative"
        )
    if not math.isclose(sum(result.values()), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
        raise RefinementEvaluationError("mixture fractions do not sum to one")
    return result


def project_machine_resolution_cells(
    *,
    energy_eV: Sequence[float],
    widths_eV: Sequence[float],
    density_eV_inv: Sequence[float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, object]]:
    """Merge only numerically unrepresentable cells, conserving their mass."""

    centers = np.asarray(energy_eV, dtype=float)
    widths = np.asarray(widths_eV, dtype=float)
    density = np.asarray(density_eV_inv, dtype=float)
    if (
        centers.ndim != 1
        or len(centers) == 0
        or widths.shape != centers.shape
        or density.shape != centers.shape
        or np.any(~np.isfinite(centers))
        or np.any(~np.isfinite(widths))
        or np.any(~np.isfinite(density))
        or np.any(widths <= 0.0)
        or np.any(density < 0.0)
    ):
        raise RefinementEvaluationError("invalid EEDF cells cannot be projected")
    left = centers - 0.5 * widths
    right = centers + 0.5 * widths
    scale = np.maximum.reduce((np.ones_like(centers), np.abs(centers), widths))
    tolerance = 256.0 * np.finfo(float).eps * scale
    if left[0] < -tolerance[0]:
        raise RefinementEvaluationError("EEDF grid extends below zero energy")
    left[0] = max(0.0, left[0])
    if len(centers) > 1:
        gaps = left[1:] - right[:-1]
        boundary_tolerance = np.maximum(tolerance[1:], tolerance[:-1])
        if np.any(np.abs(gaps) > boundary_tolerance):
            raise RefinementEvaluationError(
                "EEDF grid has a macroscopic gap or overlap; no projection is allowed"
            )
        shared = 0.5 * (right[:-1] + left[1:])
    else:
        shared = np.asarray([], dtype=float)
    nominal_edges = np.concatenate(([left[0]], shared, [right[-1]]))
    edge_scale = np.maximum(1.0, np.abs(nominal_edges))
    edge_tolerance = 256.0 * np.finfo(float).eps * edge_scale
    keep = np.ones(len(nominal_edges), dtype=bool)
    last = 0
    for index in range(1, len(nominal_edges) - 1):
        if nominal_edges[index] - nominal_edges[last] <= max(
            edge_tolerance[index], edge_tolerance[last]
        ):
            keep[index] = False
        else:
            last = index
    if nominal_edges[-1] - nominal_edges[last] <= max(
        edge_tolerance[-1], edge_tolerance[last]
    ):
        keep[last] = False
    keep[-1] = True
    projected_edges = nominal_edges[keep]
    if len(projected_edges) < 2 or np.any(np.diff(projected_edges) <= 0.0):
        raise RefinementEvaluationError("machine-resolution projection is invalid")
    if int(np.count_nonzero(~keep)) == 0:
        raise RefinementEvaluationError(
            "EEDF grid is invalid for a reason other than machine-resolution cells"
        )

    source_mass = density * widths
    nominal_widths = np.diff(nominal_edges)
    collapsed = nominal_widths <= np.maximum(
        edge_tolerance[:-1], edge_tolerance[1:]
    )
    collapsed_mass = float(np.sum(source_mass[collapsed]))
    target_mass = np.zeros(len(projected_edges) - 1, dtype=float)
    for index, mass in enumerate(source_mass):
        midpoint = 0.5 * (nominal_edges[index] + nominal_edges[index + 1])
        target = int(np.searchsorted(projected_edges, midpoint, side="right") - 1)
        target = max(0, min(target, len(target_mass) - 1))
        target_mass[target] += float(mass)
    projected_widths = np.diff(projected_edges)
    projected_centers = 0.5 * (projected_edges[:-1] + projected_edges[1:])
    projected_density = target_mass / projected_widths
    before_mass = float(np.sum(source_mass))
    after_mass = float(np.sum(target_mass))
    before_moment = float(np.sum(centers * source_mass))
    after_moment = float(np.sum(projected_centers * target_mass))
    moment_change = after_moment - before_moment
    if collapsed_mass > MAX_REPRESENTATIONAL_MERGE_MASS:
        raise RefinementEvaluationError(
            "machine-resolution EEDF cells carry material probability mass: "
            f"{collapsed_mass:.17g}"
        )
    if abs(moment_change) > MAX_REPRESENTATIONAL_MOMENT_CHANGE_EV:
        raise RefinementEvaluationError(
            "machine-resolution EEDF merge changes the first moment materially: "
            f"{moment_change:.17g} eV"
        )
    return projected_centers, projected_widths, projected_density, {
        "method": "machine_resolution_cell_merge_without_smoothing",
        "input_cells": int(len(centers)),
        "output_cells": int(len(projected_centers)),
        "removed_cells": int(len(centers) - len(projected_centers)),
        "collapsed_cell_probability_mass": collapsed_mass,
        "maximum_allowed_collapsed_mass": MAX_REPRESENTATIONAL_MERGE_MASS,
        "probability_mass_change": after_mass - before_mass,
        "first_moment_change_eV": moment_change,
        "maximum_allowed_first_moment_change_eV": (
            MAX_REPRESENTATIONAL_MOMENT_CHANGE_EV
        ),
    }


def _case_from_rows(
    *,
    solver: str,
    field: float,
    reported_mean: float,
    rows: Sequence[sqlite3.Row],
    projection_context: Mapping[str, object],
) -> tuple[EedfCase, dict[str, object] | None]:
    if not rows:
        raise RefinementEvaluationError(f"missing EEDF bins at {field:g} Td")
    energy = [float(row["energy_eV"]) for row in rows]
    widths = [float(row["energy_width_eV"]) for row in rows]
    density = [float(row["eedf"]) for row in rows]
    projection: dict[str, object] | None = None
    try:
        energy_edges_from_cells(energy, widths)
    except ValueError:
        energy_values, width_values, density_values, projection = (
            project_machine_resolution_cells(
                energy_eV=energy,
                widths_eV=widths,
                density_eV_inv=density,
            )
        )
        energy = energy_values.tolist()
        widths = width_values.tolist()
        density = density_values.tolist()
        projection = {**dict(projection_context), **projection}
    try:
        case = validated_case(
            solver=solver,
            e_over_n_td=field,
            energy_eV=energy,
            widths_eV=widths,
            density_eV_inv=density,
            reported_mean_energy_eV=reported_mean,
        )
    except EedfComparisonError as exc:
        raise RefinementEvaluationError(str(exc)) from exc
    return case, projection


def _tail_evidence(
    diagnostics_json: str,
    *,
    database_tail_schema: str | None,
    sampling_plan: Mapping[str, object],
    context: str,
) -> dict[str, object]:
    try:
        diagnostics = json.loads(diagnostics_json)
    except json.JSONDecodeError as exc:
        raise RefinementEvaluationError(f"{context} diagnostics JSON is invalid") from exc
    diagnostics = require_mapping(diagnostics, f"{context} diagnostics")
    transport = require_mapping(
        diagnostics.get("internal_monte_carlo_transport"),
        f"{context} internal_monte_carlo_transport",
    )
    run = require_mapping(transport.get("mc_run_provenance"), f"{context} MC run")
    tail = require_mapping(transport.get("tail_sampling"), f"{context} tail sampling")
    block_lag = require_mapping(
        transport.get("block_lag_sampling"), f"{context} block lag sampling"
    )
    configured = require_int(
        tail.get("configured_max_collisions"),
        f"{context} tail configured collisions",
        minimum=1,
    )
    configured_run = require_int(
        run.get("tail_max_collisions"),
        f"{context} run tail configured collisions",
        minimum=1,
    )
    executed = require_int(
        tail.get("executed_collisions"), f"{context} tail executed collisions"
    )
    executed_run = require_int(
        run.get("tail_collisions_executed"),
        f"{context} run tail executed collisions",
    )
    if configured != configured_run or executed != executed_run:
        raise RefinementEvaluationError(
            f"{context} tail sampling and run provenance disagree"
        )
    if executed > configured:
        raise RefinementEvaluationError(
            f"{context} executed tail budget exceeds configured budget"
        )
    model = str(tail.get("model", "")).strip()
    raw_schema = tail.get("estimator_schema_version")
    estimator_schema = None if raw_schema is None else str(raw_schema).strip() or None
    if executed not in {0, configured}:
        raise RefinementEvaluationError(
            f"{context} tail execution must be zero or the complete configured budget"
        )
    expected_model = (
        "reaction_kernel_weighted_ensemble"
        if executed > 0
        else "ordinary_trajectory_sampling_resolved"
    )
    if model != expected_model:
        raise RefinementEvaluationError(
            f"{context} tail model {model!r} does not match execution status"
        )
    if estimator_schema != database_tail_schema:
        raise RefinementEvaluationError(
            f"{context} tail estimator schema disagrees with database provenance"
        )
    strata = tail.get("strata_edges_eV")
    if not isinstance(strata, list) or len(strata) < 2:
        raise RefinementEvaluationError(f"{context} tail strata are missing")
    if tail.get("reported_eedf_includes_main_production") is not True:
        raise RefinementEvaluationError(
            f"{context} reported EEDF omits main production"
        )
    if tail.get("reported_rates_include_main_production") is not True:
        raise RefinementEvaluationError(
            f"{context} reported rates omit main production"
        )
    particles = require_int(run.get("particles"), f"{context} particles", minimum=1)
    warmup = require_int(run.get("warmup_collisions"), f"{context} warmup collisions")
    production = require_int(
        run.get("production_collisions"),
        f"{context} production collisions",
        minimum=1,
    )
    lag = require_int(
        block_lag.get("configured_correlation_lag_barriers"),
        f"{context} configured correlation lag",
        minimum=1,
    )
    transport_estimator = str(run.get("transport_estimator", "")).strip()
    expected = {
        "particles": particles,
        "warmup_collisions": warmup,
        "max_collisions": production,
        "tail_max_collisions": configured,
        "transport_correlation_lag_barriers": lag,
        "transport_estimator": transport_estimator,
    }
    for name, actual in expected.items():
        if sampling_plan.get(name) != actual:
            raise RefinementEvaluationError(
                f"{context} case provenance {name}={actual!r} disagrees with "
                f"sampling plan {sampling_plan.get(name)!r}"
            )
    return {
        "particles": particles,
        "warmup_collisions": warmup,
        "production_collisions": production,
        "transport_correlation_lag_barriers": lag,
        "transport_estimator": transport_estimator,
        "seed": require_int(run.get("seed"), f"{context} seed"),
        "tail_configured": 1,
        "tail_configured_collisions": configured,
        "tail_executed": int(executed > 0),
        "tail_executed_collisions": executed,
        "tail_model": model,
        "tail_estimator_schema": estimator_schema,
        "database_tail_estimator_schema": database_tail_schema,
    }


def load_database(
    source: DatabaseInput,
    *,
    mixture_id: int,
    load_replicas: bool,
) -> LoadedDatabase:
    path = source.path.resolve()
    if not path.is_file():
        raise RefinementEvaluationError(f"missing MC database: {path}")
    live_sidecars = [
        sidecar
        for sidecar in (
            path.with_name(path.name + "-wal"),
            path.with_name(path.name + "-journal"),
        )
        if sidecar.is_file() and sidecar.stat().st_size > 0
    ]
    if live_sidecars:
        raise RefinementEvaluationError(
            f"MC database has an active SQLite journal: {live_sidecars}"
        )
    source_hash = file_sha256(path)
    uri = f"file:{path.as_posix()}?mode=ro"
    try:
        connection = sqlite3.connect(uri, uri=True)
        connection.row_factory = sqlite3.Row
        integrity = connection.execute("PRAGMA quick_check").fetchone()
        if integrity is None or str(integrity[0]) != "ok":
            raise RefinementEvaluationError(
                f"MC database quick_check failed: {path}: {integrity}"
            )
        metadata = _read_metadata(connection)
        (
            refinement_context_hash,
            sampling_contract_hash,
            sampling_plan,
            replicas_by_anchor,
            campaign_limits,
        ) = _refinement_provenance(metadata, declared_population=source.population)
        fractions = _read_fractions(connection, mixture_id)
        case_rows = connection.execute(
            """
            SELECT e_over_n_Td, replicate, mean_energy_eV, diagnostics_json
            FROM cases
            WHERE mixture_id = ? AND solver = 'monte_carlo'
            ORDER BY e_over_n_Td, replicate
            """,
            (int(mixture_id),),
        ).fetchall()
        if not case_rows:
            raise RefinementEvaluationError(f"no Monte Carlo cases in {path}")
        fields = {canonical_anchor(float(row["e_over_n_Td"])) for row in case_rows}
        missing = [anchor for anchor in REQUIRED_ANCHORS_TD if anchor not in fields]
        extra = sorted(fields - set(REQUIRED_ANCHORS_TD))
        if missing or extra:
            raise RefinementEvaluationError(
                f"{path} anchor set differs from the refinement plan "
                f"(missing={missing}, extra={extra})"
            )
        database_tail_schema = metadata["mc_tail_estimator_schema_version"]
        tail_rows: list[dict[str, object]] = []
        final_cases: dict[float, list[tuple[int, EedfCase]]] = {}
        projection_events: list[dict[str, object]] = []
        seen_case_keys: set[tuple[float, int]] = set()
        replicas_by_field: dict[float, set[int]] = {}
        seeds: set[int] = set()
        for row in case_rows:
            field = canonical_anchor(float(row["e_over_n_Td"]))
            replicate = int(row["replicate"])
            key = (field, replicate)
            if key in seen_case_keys:
                raise RefinementEvaluationError(f"duplicate MC case {key} in {path}")
            seen_case_keys.add(key)
            replicas_by_field.setdefault(field, set()).add(replicate)
            context = (
                f"{source.mixture} p{source.population} {field:g} Td "
                f"replica {replicate}"
            )
            evidence = _tail_evidence(
                str(row["diagnostics_json"]),
                database_tail_schema=database_tail_schema,
                sampling_plan=sampling_plan[field],
                context=context,
            )
            seed = int(evidence["seed"])
            expected_seed = stable_mc_seed(
                base_seed=int(metadata["mc_base_seed"]),
                mixture_id=mixture_id,
                e_over_n_Td=field,
                replicate=replicate,
                solver="monte_carlo",
            )
            if seed != expected_seed:
                raise RefinementEvaluationError(
                    f"{context} seed {seed} does not match stable seed {expected_seed}"
                )
            if seed in seeds:
                raise RefinementEvaluationError(
                    f"{path} reuses MC seed {seed} across case replicas"
                )
            seeds.add(seed)
            if field == HIGH_FIELD_TD and evidence["particles"] != source.population:
                raise RefinementEvaluationError(
                    f"{context} records {evidence['particles']} particles, not "
                    f"the declared refinement population {source.population}"
                )
            tail_rows.append(
                {
                    "mixture": source.mixture,
                    "database_population": source.population,
                    "E_over_N_Td": field,
                    "replicate": replicate,
                    **evidence,
                    "source_database": str(path),
                    "source_database_sha256": source_hash,
                }
            )
            if load_replicas or field == HIGH_FIELD_TD:
                bins = connection.execute(
                    """
                    SELECT energy_eV, energy_width_eV, eedf
                    FROM eedf_bins
                    WHERE mixture_id = ? AND solver = 'monte_carlo'
                      AND e_over_n_Td = ? AND replicate = ?
                    ORDER BY bin_index
                    """,
                    (int(mixture_id), float(row["e_over_n_Td"]), replicate),
                ).fetchall()
                case, projection = _case_from_rows(
                    solver=f"monte_carlo_replica_{replicate}",
                    field=field,
                    reported_mean=float(row["mean_energy_eV"]),
                    rows=bins,
                    projection_context={
                        "mixture": source.mixture,
                        "database_population": source.population,
                        "E_over_N_Td": field,
                        "replicate": replicate,
                        "source": "eedf_bins",
                    },
                )
                final_cases.setdefault(field, []).append((replicate, case))
                if projection is not None:
                    projection_events.append(projection)

        for field, plan in sampling_plan.items():
            expected_replicas = require_int(
                plan.get("replicas"),
                f"{source.mixture} p{source.population} {field:g} Td planned replicas",
                minimum=2,
            )
            actual_replicas = replicas_by_field.get(field, set())
            if actual_replicas != set(range(expected_replicas)):
                raise RefinementEvaluationError(
                    f"{source.mixture} p{source.population} {field:g} Td replicas "
                    f"{sorted(actual_replicas)} do not match plan 0..{expected_replicas - 1}"
                )

        aggregate_mean_row = connection.execute(
            """
            SELECT mean
            FROM aggregate_scalars
            WHERE mixture_id = ? AND solver = 'monte_carlo'
              AND e_over_n_Td = ? AND scalar_group = 'case'
              AND scalar_name = 'mean_energy_eV'
            """,
            (int(mixture_id), HIGH_FIELD_TD),
        ).fetchall()
        if len(aggregate_mean_row) != 1 or aggregate_mean_row[0]["mean"] is None:
            raise RefinementEvaluationError(
                f"{path} lacks one aggregate mean-energy result at {HIGH_FIELD_TD:g} Td"
            )
        aggregate_bins = connection.execute(
            """
            SELECT energy_eV, energy_width_eV, eedf
            FROM aggregate_eedf_bins
            WHERE mixture_id = ? AND solver = 'monte_carlo' AND e_over_n_Td = ?
            ORDER BY bin_index
            """,
            (int(mixture_id), HIGH_FIELD_TD),
        ).fetchall()
        aggregate_case, projection = _case_from_rows(
            solver=f"monte_carlo_p{source.population}_aggregate",
            field=HIGH_FIELD_TD,
            reported_mean=float(aggregate_mean_row[0]["mean"]),
            rows=aggregate_bins,
            projection_context={
                "mixture": source.mixture,
                "database_population": source.population,
                "E_over_N_Td": HIGH_FIELD_TD,
                "replicate": None,
                "source": "aggregate_eedf_bins",
            },
        )
        if projection is not None:
            projection_events.append(projection)
        _validate_aggregate_against_replicas(
            aggregate_case,
            final_cases[HIGH_FIELD_TD],
            context=f"{source.mixture} p{source.population} {HIGH_FIELD_TD:g} Td",
        )
    except sqlite3.Error as exc:
        raise RefinementEvaluationError(f"could not read MC database {path}: {exc}") from exc
    finally:
        if "connection" in locals():
            connection.close()
    return LoadedDatabase(
        source=DatabaseInput(source.mixture, source.population, path),
        source_sha256=source_hash,
        fractions=fractions,
        tail_rows=tail_rows,
        aggregate_high_field=aggregate_case,
        final_replica_cases=final_cases,
        projection_events=projection_events,
        metadata=metadata,
        refinement_context_sha256=refinement_context_hash,
        sampling_contract_sha256=sampling_contract_hash,
        replicas_by_anchor=replicas_by_anchor,
        campaign_limits=campaign_limits,
    )


def _validate_aggregate_against_replicas(
    aggregate: EedfCase,
    replicas: Sequence[tuple[int, EedfCase]],
    *,
    context: str,
) -> None:
    if len(replicas) < 2:
        raise RefinementEvaluationError(
            f"{context} needs at least two replicas for aggregate reconstruction"
        )
    edges = union_edges([aggregate, *(case for _replica, case in replicas)])
    aggregate_mass = mass_on(aggregate, edges)
    replica_mass = np.asarray(
        [mass_on(case, edges) for _replica, case in replicas], dtype=float
    )
    reconstructed_mass = np.mean(replica_mass, axis=0)
    total_variation = 0.5 * float(np.sum(np.abs(aggregate_mass - reconstructed_mass)))
    if total_variation > MAX_AGGREGATE_RECONSTRUCTION_TV:
        raise RefinementEvaluationError(
            f"{context} aggregate EEDF does not reproduce the raw replica mean: "
            f"TV={total_variation:.17g}"
        )
    reconstructed_mean = float(
        np.mean([case.reported_mean_energy_eV for _replica, case in replicas])
    )
    relative_error = abs(aggregate.reported_mean_energy_eV - reconstructed_mean) / abs(
        reconstructed_mean
    )
    if relative_error > MAX_AGGREGATE_MEAN_RELATIVE_ERROR:
        raise RefinementEvaluationError(
            f"{context} aggregate mean energy does not reproduce the raw replica mean: "
            f"relative error={relative_error:.17g}"
        )


def load_qualifications(path: Path) -> dict[float, QualificationStatus]:
    statuses: dict[float, QualificationStatus] = {}
    for row in csv_rows(path):
        field = canonical_anchor(finite_float(row.get("E_over_N_Td"), "qualification E/N"))
        if field in statuses:
            raise RefinementEvaluationError(
                f"duplicate qualification status at {field:g} Td in {path}"
            )
        explicit_eedf: bool | None = None
        for name in (
            "eedf_quality_passed",
            "eedf_distribution_quality_passed",
            "eedf_shape_quality_passed",
        ):
            if name in row:
                explicit_eedf = _parse_bool(row, name)
                break
        normalization_passed: bool | None = None
        error_raw = str(row.get("eedf_normalization_error", "")).strip()
        thresholds_raw = str(row.get("thresholds_json", "")).strip()
        if error_raw and thresholds_raw:
            try:
                threshold_payload = json.loads(thresholds_raw)
                threshold = float(threshold_payload["eedf_normalization_error"])
            except (json.JSONDecodeError, KeyError, TypeError, ValueError) as exc:
                raise RefinementEvaluationError(
                    f"invalid EEDF normalization threshold in {path}"
                ) from exc
            normalization_passed = float(error_raw) <= threshold
        statuses[field] = QualificationStatus(
            profile=str(row.get("qualification_profile", "")).strip() or None,
            overall_passed=_parse_bool(row, "passed"),
            aggregate_quality_passed=_parse_bool(row, "aggregate_quality_passed"),
            active_closure_passed=_parse_bool(row, "active_closure_quality_passed"),
            solver_transport_passed=_parse_bool(row, "solver_transport_qualified"),
            explicit_eedf_passed=explicit_eedf,
            eedf_normalization_passed=normalization_passed,
        )
    if not statuses:
        raise RefinementEvaluationError(f"qualification CSV is empty: {path}")
    if set(statuses) != set(REQUIRED_ANCHORS_TD):
        raise RefinementEvaluationError(
            "qualification CSV must contain exactly 1, 10, 100, and 1000 Td"
        )
    return statuses


def validate_qualification_manifest(
    csv_path: Path,
    *,
    csv_sha256: str,
    statuses: Mapping[float, QualificationStatus],
    final_database: LoadedDatabase,
) -> dict[str, str]:
    manifest_path = csv_path.parent / "manifest.json"
    if not manifest_path.is_file():
        raise RefinementEvaluationError(
            f"qualification CSV requires its sibling manifest.json: {csv_path}"
        )
    try:
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RefinementEvaluationError(
            f"qualification manifest is unreadable: {manifest_path}"
        ) from exc
    payload = require_mapping(payload, "qualification manifest")
    if payload.get("source") != "monte_carlo":
        raise RefinementEvaluationError("qualification manifest is not Monte Carlo")
    tables = require_mapping(payload.get("tables"), "qualification manifest tables")
    table = require_mapping(
        tables.get(csv_path.name), f"qualification manifest table {csv_path.name}"
    )
    if table.get("sha256") != csv_sha256:
        raise RefinementEvaluationError(
            "qualification CSV hash does not match its table manifest"
        )
    hashes = require_mapping(payload.get("hashes"), "qualification manifest hashes")
    database_metadata = final_database.metadata
    for key in (
        "workflow_config_sha256",
        "base_config_sha256",
        "cross_sections_sha256",
        "mc_sampling_plan_json",
        "mc_transport_estimator_schema_version",
        "mc_eedf_estimator_schema_version",
        "mc_tail_estimator_schema_version",
    ):
        if hashes.get(key) != database_metadata.get(key):
            raise RefinementEvaluationError(
                f"qualification manifest {key} does not bind to the final MC database"
            )
    try:
        manifest_context = json.dumps(
            payload["physical_context"], sort_keys=True, separators=(",", ":")
        )
        database_context = json.dumps(
            json.loads(database_metadata["physical_context_json"]),
            sort_keys=True,
            separators=(",", ":"),
        )
    except (KeyError, TypeError, json.JSONDecodeError) as exc:
        raise RefinementEvaluationError(
            "qualification physical-context provenance is invalid"
        ) from exc
    if manifest_context != database_context:
        raise RefinementEvaluationError(
            "qualification physical context does not bind to the final MC database"
        )
    mixture = require_mapping(payload.get("mixture"), "qualification mixture")
    species_rows = mixture.get("species")
    if not isinstance(species_rows, list):
        raise RefinementEvaluationError("qualification mixture species are missing")
    try:
        manifest_fractions = {
            str(item["species"]): float(item["fraction"])
            for item in species_rows
            if isinstance(item, dict)
        }
    except (KeyError, TypeError, ValueError) as exc:
        raise RefinementEvaluationError(
            "qualification mixture fractions are invalid"
        ) from exc
    if manifest_fractions != final_database.fractions:
        raise RefinementEvaluationError(
            "qualification mixture does not bind to the final MC database"
        )
    profiles = {status.profile for status in statuses.values()}
    expected_profile = require_mapping(
        payload.get("source_policy"), "qualification source policy"
    ).get("qualification_profile")
    if profiles != {expected_profile}:
        raise RefinementEvaluationError(
            "qualification CSV profile does not match its manifest"
        )
    return {
        "path": str(csv_path),
        "sha256": csv_sha256,
        "manifest_path": str(manifest_path.resolve()),
        "manifest_sha256": file_sha256(manifest_path),
        "bound_final_database": str(final_database.source.path),
        "bound_workflow_config_sha256": database_metadata["workflow_config_sha256"],
    }


__all__ = [
    "canonical_anchor",
    "file_sha256",
    "load_database",
    "load_qualifications",
    "project_machine_resolution_cells",
    "require_mapping",
    "validate_qualification_manifest",
]
