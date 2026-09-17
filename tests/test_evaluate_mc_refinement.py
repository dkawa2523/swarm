from __future__ import annotations

import csv
from hashlib import sha256
import json
from pathlib import Path
import sqlite3

import numpy as np
import pytest

from swarm_workflow.campaign.aggregate import stable_mc_seed
from tools.validation.evaluate_mc_refinement import (
    DatabaseInput,
    RefinementEvaluationError,
    _project_machine_resolution_cells,
    evaluate_mc_refinement,
)


MIXTURES = {
    "ar_o2": {"Ar": 0.9, "O2": 0.1},
    "ar_n2": {"Ar": 0.9, "N2": 0.1},
    "ar_cl2": {"Ar": 0.9, "Cl2": 0.1},
}
POPULATIONS = (512, 1024, 2048)
def _population_grid(population: int, field: float) -> tuple[np.ndarray, np.ndarray]:
    if field != 1000.0:
        base_mass = np.asarray([0.90, 0.05, 0.03, 0.02])
    else:
        base_mass = {
        512: np.asarray([0.70, 0.12, 0.10, 0.08]),
        1024: np.asarray([0.72, 0.11, 0.095, 0.075]),
        2048: np.asarray([0.73, 0.105, 0.09, 0.075]),
        }[population]
    if population == 512:
        return np.asarray([0.0, 80.0, 100.0, 120.0, 160.0]), base_mass
    if population == 1024:
        return (
            np.asarray([0.0, 40.0, 80.0, 100.0, 120.0, 160.0]),
            np.concatenate(([0.5 * base_mass[0], 0.5 * base_mass[0]], base_mass[1:])),
        )
    return (
        np.asarray([0.0, 20.0, 60.0, 80.0, 100.0, 120.0, 160.0]),
        np.concatenate(
            (
                [0.25 * base_mass[0], 0.5 * base_mass[0], 0.25 * base_mass[0]],
                base_mass[1:],
            )
        ),
    )


def _write_database(
    path: Path,
    *,
    mixture: dict[str, float],
    population: int,
    replicas_by_field: dict[float, int] | None = None,
    campaign_maximum_replicas: int = 2,
    campaign_maximum_total_particle_barriers: int = 1_000_000,
) -> None:
    path.parent.mkdir(parents=True)
    connection = sqlite3.connect(path)
    connection.executescript(
        """
        CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);
        CREATE TABLE mixtures (mixture_id INTEGER PRIMARY KEY, fractions_json TEXT NOT NULL);
        CREATE TABLE cases (
            mixture_id INTEGER, solver TEXT, e_over_n_Td REAL, replicate INTEGER,
            mean_energy_eV REAL, diagnostics_json TEXT
        );
        CREATE TABLE eedf_bins (
            mixture_id INTEGER, solver TEXT, e_over_n_Td REAL, replicate INTEGER,
            bin_index INTEGER, energy_eV REAL, energy_width_eV REAL, eedf REAL
        );
        CREATE TABLE aggregate_scalars (
            mixture_id INTEGER, solver TEXT, e_over_n_Td REAL,
            scalar_group TEXT, scalar_name TEXT, mean REAL
        );
        CREATE TABLE aggregate_eedf_bins (
            mixture_id INTEGER, solver TEXT, e_over_n_Td REAL, bin_index INTEGER,
            energy_eV REAL, energy_width_eV REAL, eedf REAL
        );
        """
    )
    replica_counts = replicas_by_field or {
        field: 2 for field in (1.0, 10.0, 100.0, 1000.0)
    }
    sampling_plan = [
        {
            "e_over_n_Td": field,
            "particles": {1.0: 128, 10.0: 256, 100.0: 512}.get(
                field, population
            ),
            "warmup_collisions": 16,
            "max_collisions": 32,
            "tail_max_collisions": 65536 if field in (1.0, 10.0) else 16384,
            "replicas": replica_counts[field],
            "transport_correlation_lag_barriers": 4,
            "transport_estimator": "single_field",
        }
        for field in (1.0, 10.0, 100.0, 1000.0)
    ]
    metadata = {
        "base_config_sha256": "a" * 64,
        "cross_section_files_json": json.dumps(
            {"fixture.txt": "b" * 64}, sort_keys=True, separators=(",", ":")
        ),
        "cross_sections_sha256": "c" * 64,
        "mc_base_seed": "12345",
        "mc_campaign_json": json.dumps(
            {
                "policy": {
                    "initial_replicas": 2,
                    "replica_increment": 2,
                    "maximum_replicas": campaign_maximum_replicas,
                    "maximum_attempts": 1,
                    "time_extension_factor": 2,
                    "maximum_particle_barriers_per_replica": 100_000,
                    "maximum_total_particle_barriers": (
                        campaign_maximum_total_particle_barriers
                    ),
                },
                "previous_decision": None,
            },
            sort_keys=True,
            separators=(",", ":"),
        ),
        "mc_sampling_plan_json": json.dumps(
            sampling_plan, sort_keys=True, separators=(",", ":")
        ),
        "mc_tail_estimator_schema_version": (
            "reaction_kernel_weighted_ensemble_pooled_estimators.v4"
        ),
        "mc_eedf_estimator_schema_version": (
            "exact_bin_residence_piecewise_linear_rate_integrals.v2"
        ),
        "mc_transport_estimator_schema_version": "fixture_transport.v1",
        "physical_context_json": json.dumps(
            {"schema": "fixture_context.v1"},
            sort_keys=True,
            separators=(",", ":"),
        ),
        "quality_thresholds_json": json.dumps(
            {"eedf_normalization_error": 1.0e-6},
            sort_keys=True,
            separators=(",", ":"),
        ),
        "workflow_config_sha256": f"{population:064x}",
    }
    connection.executemany(
        "INSERT INTO metadata VALUES (?, ?)", sorted(metadata.items())
    )
    connection.execute(
        "INSERT INTO mixtures VALUES (?, ?)",
        (0, json.dumps(mixture, sort_keys=True)),
    )
    for field in (1.0, 10.0, 100.0, 1000.0):
        edges, aggregate_mass = _population_grid(population, field)
        centers = 0.5 * (edges[:-1] + edges[1:])
        widths = np.diff(edges)
        aggregate_mean = float(np.sum(centers * aggregate_mass))
        connection.execute(
            "INSERT INTO aggregate_scalars VALUES (?, ?, ?, ?, ?, ?)",
            (0, "monte_carlo", field, "case", "mean_energy_eV", aggregate_mean),
        )
        for index, (energy, width, mass) in enumerate(
            zip(centers, widths, aggregate_mass, strict=True)
        ):
            connection.execute(
                "INSERT INTO aggregate_eedf_bins VALUES (?, ?, ?, ?, ?, ?, ?)",
                (0, "monte_carlo", field, index, energy, width, mass / width),
            )
        for replicate in range(replica_counts[field]):
            sign = -1.0 if replicate % 2 == 0 else 1.0
            mass = aggregate_mass.copy()
            mass[0] += sign * 0.005
            mass[-1] -= sign * 0.005
            mean = float(np.sum(centers * mass))
            configured = 65536 if field in (1.0, 10.0) else 16384
            executed = configured if field == 10.0 or (field == 1.0 and replicate == 0) else 0
            tail = {
                "configured_max_collisions": configured,
                "executed_collisions": executed,
                "model": (
                    "reaction_kernel_weighted_ensemble"
                    if executed
                    else "ordinary_trajectory_sampling_resolved"
                ),
                "estimator_schema_version": (
                    "reaction_kernel_weighted_ensemble_pooled_estimators.v4"
                ),
                "strata_edges_eV": [0.0, 80.0, 160.0],
                "reported_eedf_includes_main_production": True,
                "reported_rates_include_main_production": True,
            }
            diagnostics = {
                "internal_monte_carlo_transport": {
                    "mc_run_provenance": {
                        "particles": {1.0: 128, 10.0: 256, 100.0: 512}.get(
                            field, population
                        ),
                        "seed": stable_mc_seed(
                            base_seed=12345,
                            mixture_id=0,
                            e_over_n_Td=field,
                            replicate=replicate,
                            solver="monte_carlo",
                        ),
                        "warmup_collisions": 16,
                        "production_collisions": 32,
                        "tail_max_collisions": configured,
                        "tail_collisions_executed": executed,
                        "transport_estimator": "single_field",
                    },
                    "tail_sampling": tail,
                    "block_lag_sampling": {
                        "configured_correlation_lag_barriers": 4
                    },
                }
            }
            connection.execute(
                "INSERT INTO cases VALUES (?, ?, ?, ?, ?, ?)",
                (0, "monte_carlo", field, replicate, mean, json.dumps(diagnostics)),
            )
            for index, (energy, width, probability) in enumerate(
                zip(centers, widths, mass, strict=True)
            ):
                connection.execute(
                    "INSERT INTO eedf_bins VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                    (
                        0,
                        "monte_carlo",
                        field,
                        replicate,
                        index,
                        energy,
                        width,
                        probability / width,
                    ),
                )
    connection.commit()
    connection.close()


def _write_qualification(
    path: Path,
    *,
    database: Path,
    mixture: dict[str, float],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "E_over_N_Td",
                "qualification_profile",
                "passed",
                "aggregate_quality_passed",
                "active_closure_quality_passed",
                "solver_transport_qualified",
                "eedf_quality_passed",
                "eedf_normalization_error",
                "thresholds_json",
            ),
            lineterminator="\n",
        )
        writer.writeheader()
        for field in (1.0, 10.0, 100.0, 1000.0):
            writer.writerow(
                {
                    "E_over_N_Td": field,
                    "qualification_profile": "function_eedf_restricted_lmea",
                    "passed": 1,
                    "aggregate_quality_passed": 1,
                    "active_closure_quality_passed": int(field != 1000.0),
                    "solver_transport_qualified": 1,
                    "eedf_quality_passed": int(field != 1000.0),
                    "eedf_normalization_error": 0.0,
                    "thresholds_json": json.dumps(
                        {"eedf_normalization_error": 1.0e-6}
                    ),
                }
            )
    csv_hash = sha256(path.read_bytes()).hexdigest()
    connection = sqlite3.connect(database)
    metadata = dict(connection.execute("SELECT key, value FROM metadata"))
    connection.close()
    manifest = {
        "source": "monte_carlo",
        "hashes": {
            key: metadata[key]
            for key in (
                "workflow_config_sha256",
                "base_config_sha256",
                "cross_sections_sha256",
                "mc_sampling_plan_json",
                "mc_transport_estimator_schema_version",
                "mc_eedf_estimator_schema_version",
                "mc_tail_estimator_schema_version",
            )
        },
        "physical_context": json.loads(metadata["physical_context_json"]),
        "mixture": {
            "species": [
                {"species": species, "fraction": fraction}
                for species, fraction in mixture.items()
            ]
        },
        "source_policy": {
            "qualification_profile": "function_eedf_restricted_lmea"
        },
        "tables": {path.name: {"sha256": csv_hash}},
    }
    (path.parent / "manifest.json").write_text(
        json.dumps(manifest), encoding="utf-8"
    )


def test_evaluator_keeps_raw_population_and_qualification_evidence_separate(
    tmp_path: Path,
) -> None:
    pytest.importorskip("matplotlib")
    inputs: list[DatabaseInput] = []
    qualifications: dict[str, Path] = {}
    for mixture_name, mixture in MIXTURES.items():
        for population in POPULATIONS:
            database = tmp_path / mixture_name / f"p{population}" / "monte_carlo.sqlite"
            _write_database(database, mixture=mixture, population=population)
            inputs.append(DatabaseInput(mixture_name, population, database))
        qualification = tmp_path / mixture_name / "mc_qualification.csv"
        _write_qualification(
            qualification,
            database=tmp_path / mixture_name / "p2048" / "monte_carlo.sqlite",
            mixture=mixture,
        )
        qualifications[mixture_name] = qualification

    paths = evaluate_mc_refinement(
        inputs,
        output_directory=tmp_path / "evaluation",
        qualification_paths=qualifications,
    )

    assert all(path.is_file() and path.stat().st_size > 0 for path in paths.values())
    tail_rows = list(csv.DictReader(paths["tail_execution"].open(encoding="utf-8")))
    assert len(tail_rows) == 3 * 3 * 4 * 2
    tail_summary = list(csv.DictReader(paths["tail_summary"].open(encoding="utf-8")))
    p512_o2_1td = next(
        row
        for row in tail_summary
        if row["mixture"] == "ar_o2"
        and row["database_population"] == "512"
        and row["E_over_N_Td"] == "1.0"
    )
    assert float(p512_o2_1td["tail_execution_fraction"]) == pytest.approx(0.5)

    stability = list(
        csv.DictReader(paths["replica_stability"].open(encoding="utf-8"))
    )
    assert len(stability) == 3 * 4
    high_o2 = next(
        row
        for row in stability
        if row["mixture"] == "ar_o2" and row["E_over_N_Td"] == "1000.0"
    )
    assert float(high_o2["pairwise_total_variation_max"]) == pytest.approx(0.01)
    assert "qualification_explicit_eedf_passed" not in high_o2

    refinement = list(
        csv.DictReader(paths["population_refinement"].open(encoding="utf-8"))
    )
    assert len(refinement) == 3 * 2
    first = next(
        row
        for row in refinement
        if row["mixture"] == "ar_o2" and row["reference_population"] == "512"
    )
    assert float(first["total_variation"]) == pytest.approx(0.02)
    assert float(first["reference_tail_mass_gt_80_eV"]) == pytest.approx(0.30)
    assert float(first["candidate_tail_mass_gt_80_eV"]) == pytest.approx(0.28)
    assert "qualification_solver_transport_passed" not in first

    qualification_rows = list(
        csv.DictReader(paths["qualification"].open(encoding="utf-8"))
    )
    high_qualification = next(
        row
        for row in qualification_rows
        if row["mixture"] == "ar_o2" and row["E_over_N_Td"] == "1000.0"
    )
    assert high_qualification["qualification_explicit_eedf_passed"] == "0"
    assert high_qualification["qualification_solver_transport_passed"] == "1"

    payload = json.loads(paths["evaluation"].read_text(encoding="utf-8"))
    assert payload["comparison_scope"]["fallback_used"] is False
    assert payload["comparison_scope"]["smoothing_used"] is False
    manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    assert len(manifest["inputs"]) == 9
    assert all(len(item["sha256"]) == 64 for item in manifest["inputs"])
    svg = paths["figure_svg"].read_text(encoding="utf-8").lower()
    assert "no smoothing or fallback" in svg
    assert "two_term" not in svg


def test_evaluator_rejects_nonpopulation_context_change(tmp_path: Path) -> None:
    pytest.importorskip("matplotlib")
    inputs: list[DatabaseInput] = []
    for mixture_name, mixture in MIXTURES.items():
        for population in POPULATIONS:
            database = tmp_path / mixture_name / f"p{population}" / "monte_carlo.sqlite"
            _write_database(database, mixture=mixture, population=population)
            inputs.append(DatabaseInput(mixture_name, population, database))
    changed = tmp_path / "ar_o2" / "p1024" / "monte_carlo.sqlite"
    connection = sqlite3.connect(changed)
    connection.execute(
        "UPDATE metadata SET value = ? WHERE key = 'cross_sections_sha256'",
        ("d" * 64,),
    )
    connection.commit()
    connection.close()

    with pytest.raises(RefinementEvaluationError, match="context changes"):
        evaluate_mc_refinement(inputs, output_directory=tmp_path / "evaluation")


def test_evaluator_accepts_targeted_replica_refinement_and_reports_it(
    tmp_path: Path,
) -> None:
    pytest.importorskip("matplotlib")
    inputs: list[DatabaseInput] = []
    for population in POPULATIONS:
        database = tmp_path / "ar_n2" / f"p{population}" / "monte_carlo.sqlite"
        replicas = {
            field: (4 if population == 2048 and field == 1.0 else 2)
            for field in (1.0, 10.0, 100.0, 1000.0)
        }
        _write_database(
            database,
            mixture=MIXTURES["ar_n2"],
            population=population,
            replicas_by_field=replicas,
            campaign_maximum_replicas=4 if population == 2048 else 2,
            campaign_maximum_total_particle_barriers=(
                2_000_000 if population == 2048 else 1_000_000
            ),
        )
        inputs.append(DatabaseInput("ar_n2", population, database))

    paths = evaluate_mc_refinement(
        inputs, output_directory=tmp_path / "evaluation"
    )

    rows = list(
        csv.DictReader(paths["replica_refinement"].open(encoding="utf-8"))
    )
    final = next(
        row
        for row in rows
        if row["database_population"] == "2048" and row["E_over_N_Td"] == "1.0"
    )
    assert final["previous_replicas"] == "2"
    assert final["replicas"] == "4"
    assert final["replica_increment"] == "2"
    assert final["campaign_maximum_replicas"] == "4"
    evaluation = json.loads(paths["evaluation"].read_text(encoding="utf-8"))
    manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    assert final in [
        {key: str(value) for key, value in row.items()}
        for row in evaluation["replica_refinement"]
    ]
    assert manifest["inputs"][-1]["replicas_by_anchor"]["1"] == 4
    assert manifest["replica_refinement"][-4]["replica_increment"] == 2


@pytest.mark.parametrize(
    ("change", "message"),
    (
        ("replica_decrease", "replicas decrease"),
        ("seed", "does not match stable seed"),
        ("sampling", "sampling plans differ"),
        (
            "campaign_maximum_replicas_decrease",
            "campaign maximum_replicas decreases",
        ),
        (
            "campaign_budget_decrease",
            "campaign maximum_total_particle_barriers decreases",
        ),
    ),
)
def test_evaluator_rejects_nonmonotone_or_nonreplica_refinement(
    tmp_path: Path,
    change: str,
    message: str,
) -> None:
    pytest.importorskip("matplotlib")
    inputs: list[DatabaseInput] = []
    for population in POPULATIONS:
        database = tmp_path / "ar_n2" / f"p{population}" / "monte_carlo.sqlite"
        replica_count = (
            4 if change == "replica_decrease" and population == 1024 else 2
        )
        total_budget = (
            2_000_000
            if change == "campaign_budget_decrease" and population == 512
            else 1_000_000
        )
        maximum_replicas = (
            4
            if change == "campaign_maximum_replicas_decrease" and population == 512
            else max(2, replica_count)
        )
        _write_database(
            database,
            mixture=MIXTURES["ar_n2"],
            population=population,
            replicas_by_field={
                field: replica_count
                for field in (1.0, 10.0, 100.0, 1000.0)
            },
            campaign_maximum_replicas=maximum_replicas,
            campaign_maximum_total_particle_barriers=total_budget,
        )
        inputs.append(DatabaseInput("ar_n2", population, database))

    changed = tmp_path / "ar_n2" / "p1024" / "monte_carlo.sqlite"
    if change == "seed":
        connection = sqlite3.connect(changed)
        row = connection.execute(
            "SELECT diagnostics_json FROM cases WHERE e_over_n_Td = 1 AND replicate = 0"
        ).fetchone()
        diagnostics = json.loads(row[0])
        diagnostics["internal_monte_carlo_transport"]["mc_run_provenance"][
            "seed"
        ] += 1
        connection.execute(
            "UPDATE cases SET diagnostics_json = ? "
            "WHERE e_over_n_Td = 1 AND replicate = 0",
            (json.dumps(diagnostics),),
        )
        connection.commit()
        connection.close()
    elif change == "sampling":
        connection = sqlite3.connect(changed)
        metadata = dict(connection.execute("SELECT key, value FROM metadata"))
        plan = json.loads(metadata["mc_sampling_plan_json"])
        next(row for row in plan if row["e_over_n_Td"] == 1.0)[
            "warmup_collisions"
        ] = 17
        connection.execute(
            "UPDATE metadata SET value = ? WHERE key = 'mc_sampling_plan_json'",
            (json.dumps(plan, sort_keys=True, separators=(",", ":")),),
        )
        rows = connection.execute(
            "SELECT replicate, diagnostics_json FROM cases WHERE e_over_n_Td = 1"
        ).fetchall()
        for replicate, raw in rows:
            diagnostics = json.loads(raw)
            diagnostics["internal_monte_carlo_transport"]["mc_run_provenance"][
                "warmup_collisions"
            ] = 17
            connection.execute(
                "UPDATE cases SET diagnostics_json = ? "
                "WHERE e_over_n_Td = 1 AND replicate = ?",
                (json.dumps(diagnostics), replicate),
            )
        connection.commit()
        connection.close()

    with pytest.raises(RefinementEvaluationError, match=message):
        evaluate_mc_refinement(inputs, output_directory=tmp_path / "evaluation")


def test_machine_resolution_merge_rejects_material_probability_mass() -> None:
    with pytest.raises(RefinementEvaluationError, match="material probability mass"):
        _project_machine_resolution_cells(
            energy_eV=(0.5, 1.0 + 5.0e-16, 1.5 + 1.0e-15),
            widths_eV=(1.0, 1.0e-15, 1.0),
            density_eV_inv=(0.45, 1.0e14, 0.45),
        )
