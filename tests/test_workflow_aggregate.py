from __future__ import annotations

import json
import sqlite3
from pathlib import Path

import pytest
import yaml

from swarm_workflow.aggregate import (
    QualityThresholds,
    aggregate_database,
    aggregate_monte_carlo_results,
    conservative_rebin_probability_mass,
    stable_mc_seed,
    summarize_replicates,
)
from swarm_workflow.store import WorkflowStore
from swarm_workflow.sweep import load_workflow
from product_helpers import (
    base_product_config,
    insert_workflow_case as _insert_case,
    insert_workflow_eedf as _insert_eedf,
    insert_workflow_rate as _insert_rate,
    write_config,
)


def _build_synthetic_db(
    path: Path,
    *,
    order: tuple[int, ...] = (0, 1, 2),
    e_over_n: float = 10.0,
    normalized_eedf: bool = True,
) -> sqlite3.Connection:
    store = WorkflowStore(path)
    conn = store.connection
    values = {
        0: (1.0, 4.0, 10.0, 20.0, 30.0, 0.9e-15, (0.4, 0.6)),
        1: (2.0, 5.0, 12.0, 22.0, 33.0, 1.0e-15, (0.5, 0.5)),
        2: (3.0, 6.0, 14.0, 24.0, 36.0, 1.1e-15, (0.45, 0.55)),
    }
    if not normalized_eedf:
        values[2] = (*values[2][:6], (0.3, 0.4))
    for replicate in order:
        mean_energy, drift, mobility, diff_l, diff_t, rate, masses = values[replicate]
        _insert_case(
            conn,
            replicate=replicate,
            e_over_n=e_over_n,
            mean_energy=mean_energy,
            drift=drift,
            reduced_mobility=mobility,
            reduced_diffusion_l=diff_l,
            reduced_diffusion_t=diff_t,
        )
        _insert_rate(conn, replicate=replicate, value=rate, e_over_n=e_over_n)
        _insert_eedf(conn, replicate=replicate, probability_masses=masses, e_over_n=e_over_n)
    conn.commit()
    return conn


def test_sha256_seed_is_stable_and_order_independent() -> None:
    seeds = [
        stable_mc_seed(
            base_seed=12345,
            mixture_id=0,
            e_over_n_Td=10.0,
            replicate=index,
            solver="monte_carlo",
        )
        for index in [0, 1, 2]
    ]

    assert seeds[0] == 2181922784
    assert seeds == [
        stable_mc_seed(
            base_seed=12345,
            mixture_id=0,
            e_over_n_Td=10.0,
            replicate=index,
            solver="monte_carlo",
        )
        for index in [0, 1, 2]
    ]
    assert sorted(seeds) == sorted(reversed(seeds))
    assert seeds[0] != stable_mc_seed(
        base_seed=12345,
        mixture_id=0,
        e_over_n_Td=10.0,
        replicate=0,
        solver="two_term",
    )
    assert stable_mc_seed(
        base_seed=None,
        mixture_id=0,
        e_over_n_Td=10.0,
        replicate=0,
        solver="monte_carlo",
    ) is None


def test_scalar_replicate_statistics_fix_mean_se_and_ci() -> None:
    stats = summarize_replicates([1.0, 2.0, 3.0])

    assert stats.mean == pytest.approx(2.0)
    assert stats.sample_stddev == pytest.approx(1.0)
    assert stats.standard_error == pytest.approx(1.0 / (3.0**0.5))
    assert stats.relative_standard_error == pytest.approx(
        (1.0 / (3.0**0.5)) / 2.0
    )
    assert stats.ci95_low == pytest.approx(2.0 - 1.96 / (3.0**0.5))
    assert stats.ci95_high == pytest.approx(2.0 + 1.96 / (3.0**0.5))
    assert stats.valid_replicates == 3
    assert stats.uncertainty_available is True


def test_single_replicate_uncertainty_is_unavailable(tmp_path: Path) -> None:
    with WorkflowStore(tmp_path / "one.sqlite") as store:
        _insert_case(store.connection, replicate=0, e_over_n=20.0, mean_energy=5.0)
        _insert_rate(store.connection, replicate=0, value=1.0e-15, e_over_n=20.0)
        _insert_eedf(store.connection, replicate=0, e_over_n=20.0)
        aggregate_monte_carlo_results(store.connection)
        row = store.connection.execute(
            """
            SELECT mean, sample_stddev, standard_error, relative_standard_error,
                   valid_replicates, uncertainty_available
            FROM aggregate_scalars
            WHERE scalar_name = 'mean_energy_eV'
            """
        ).fetchone()
        quality = store.connection.execute(
            "SELECT passed, uncertainty_available FROM aggregate_quality"
        ).fetchone()

    assert row["mean"] == pytest.approx(5.0)
    assert row["sample_stddev"] is None
    assert row["standard_error"] is None
    assert row["relative_standard_error"] is None
    assert row["valid_replicates"] == 1
    assert row["uncertainty_available"] == 0
    assert quality["passed"] == 0
    assert quality["uncertainty_available"] == 0


def test_conservative_rebin_preserves_probability_mass() -> None:
    rebinned = conservative_rebin_probability_mass(
        [0.0, 1.0, 3.0],
        [0.2, 0.8],
        [0.0, 2.0, 3.0],
    )

    assert rebinned == pytest.approx([0.6, 0.4])
    assert float(sum(rebinned)) == pytest.approx(1.0)


def test_aggregate_tables_are_independent_of_replicate_insert_order(
    tmp_path: Path,
) -> None:
    conn_a = _build_synthetic_db(tmp_path / "a.sqlite", order=(0, 1, 2))
    conn_b = _build_synthetic_db(tmp_path / "b.sqlite", order=(2, 0, 1))
    loose = QualityThresholds(
        mobility_rse=0.2,
        diffusion_rse=0.2,
        major_rate_rse=0.2,
    )
    aggregate_monte_carlo_results(conn_a, loose)
    aggregate_monte_carlo_results(conn_b, loose)

    query = """
        SELECT scalar_group, scalar_name, species, process, mean,
               standard_error, relative_standard_error, valid_replicates
        FROM aggregate_scalars
        ORDER BY scalar_group, scalar_name, species, process
    """
    rows_a = [tuple(row) for row in conn_a.execute(query).fetchall()]
    rows_b = [tuple(row) for row in conn_b.execute(query).fetchall()]
    quality = conn_a.execute("SELECT passed FROM aggregate_quality").fetchone()
    eedf_norm = conn_a.execute(
        "SELECT SUM(mean_probability_mass) FROM aggregate_eedf_bins"
    ).fetchone()[0]

    conn_a.close()
    conn_b.close()
    assert rows_a == rows_b
    assert quality["passed"] == 1
    assert eedf_norm == pytest.approx(1.0)


def test_aggregate_quality_records_comsol_visible_failure(tmp_path: Path) -> None:
    with WorkflowStore(tmp_path / "quality.sqlite") as store:
        for replicate, mobility in enumerate([1.0, 10.0, 20.0]):
            _insert_case(
                store.connection,
                replicate=replicate,
                reduced_mobility=mobility,
                reduced_diffusion_l=1.0,
                reduced_diffusion_t=1.0,
            )
            _insert_rate(store.connection, replicate=replicate, value=1.0e-15)
            _insert_eedf(
                store.connection,
                replicate=replicate,
                probability_masses=(0.3, 0.4),
            )
        aggregate_monte_carlo_results(store.connection)
        row = store.connection.execute(
            """
            SELECT passed, failure_reasons_json, mobility_rse,
                   eedf_normalization_error
            FROM aggregate_quality
            """
        ).fetchone()

    reasons = json.loads(row["failure_reasons_json"])
    assert row["passed"] == 0
    assert row["mobility_rse"] > 0.02
    assert row["eedf_normalization_error"] == pytest.approx(0.3)
    assert any("reduced_mobility" in reason for reason in reasons)
    assert "eedf_normalization_error_exceeds_threshold" in reasons


def test_workflow_quality_thresholds_are_loaded_from_yaml(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["two_term"])
    base_path = write_config(tmp_path, data, "base.yaml")
    workflow_path = tmp_path / "workflow.yaml"
    workflow_path.write_text(
        yaml.safe_dump(
            {
                "base_config": base_path.name,
                "database": "outputs/swarm.sqlite",
                "e_over_n_Td": [10.0],
                "mixtures": [{"Ar": 1.0}],
                "quality": {
                    "mobility_rse": 0.1,
                    "major_rate_fraction": 0.2,
                },
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    workflow = load_workflow(workflow_path)
    assert workflow.quality.mobility_rse == pytest.approx(0.1)
    assert workflow.quality.major_rate_fraction == pytest.approx(0.2)

    bad_path = tmp_path / "bad_workflow.yaml"
    bad_path.write_text(
        yaml.safe_dump(
            {
                "base_config": base_path.name,
                "database": "outputs/swarm.sqlite",
                "e_over_n_Td": [10.0],
                "mixtures": [{"Ar": 1.0}],
                "quality": {"major_rate_fraction": 2.0},
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="major_rate_fraction"):
        load_workflow(bad_path)


def test_aggregate_cli_output_is_csv_and_manifest_only(tmp_path: Path) -> None:
    db_path = tmp_path / "aggregate.sqlite"
    conn = _build_synthetic_db(db_path)
    conn.close()

    aggregate_database(db_path, tmp_path / "aggregate")

    assert {path.name for path in (tmp_path / "aggregate").iterdir()} == {
        "manifest.json",
        "aggregate_scalars.csv",
        "aggregate_eedf_bins.csv",
        "aggregate_quality.csv",
    }
