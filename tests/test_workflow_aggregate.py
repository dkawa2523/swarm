from __future__ import annotations

import csv
import json
import math
import sqlite3
from pathlib import Path

import pytest
import yaml

from swarm_workflow.campaign.aggregate import (
    aggregate_database,
    aggregate_workflow_results,
    stable_mc_seed,
)
from swarm_workflow.campaign.statistics import (
    conservative_rebin_probability_mass,
    summarize_replicates,
)
from swarm_workflow.quality.policy import QualityThresholds, RequiredRateRse
from swarm_workflow.campaign.config import load_workflow
from swarm_workflow.campaign.store import WorkflowSchemaError, WorkflowStore
from product_helpers import (
    base_product_config,
    ensure_workflow_provenance,
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
    quality: QualityThresholds | None = None,
) -> sqlite3.Connection:
    store = WorkflowStore(path)
    conn = store.connection
    ensure_workflow_provenance(conn, quality)
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

    assert seeds[0] == 1277222181
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
    assert stats.ci95_critical_value == pytest.approx(4.3026527299)
    assert stats.ci95_low == pytest.approx(2.0 - 4.3026527299 / (3.0**0.5))
    assert stats.ci95_high == pytest.approx(2.0 + 4.3026527299 / (3.0**0.5))
    assert stats.valid_replicates == 3
    assert stats.uncertainty_available is True


def test_student_t_ci_converges_to_normal_for_many_replicates() -> None:
    stats = summarize_replicates(range(1, 10_002))

    assert stats.ci95_critical_value == pytest.approx(1.96, rel=2.0e-4)


def test_all_zero_rate_replicates_are_censored_not_precise() -> None:
    stats = summarize_replicates([0.0, 0.0, 0.0], zero_is_censored=True)

    assert stats.mean == 0.0
    assert stats.sample_stddev == 0.0
    assert stats.standard_error is None
    assert stats.relative_standard_error is None
    assert stats.ci95_low is None
    assert stats.ci95_high is None
    assert stats.uncertainty_available is False
    assert stats.estimate_status == "censored_all_zero"

    exact_zero = summarize_replicates([0.0, 0.0, 0.0])
    assert exact_zero.relative_standard_error == 0.0
    assert exact_zero.estimate_status == "estimate"


def test_single_replicate_uncertainty_is_unavailable(tmp_path: Path) -> None:
    with WorkflowStore(tmp_path / "one.sqlite") as store:
        _insert_case(store.connection, replicate=0, e_over_n=20.0, mean_energy=5.0)
        _insert_rate(store.connection, replicate=0, value=1.0e-15, e_over_n=20.0)
        _insert_eedf(store.connection, replicate=0, e_over_n=20.0)
        aggregate_workflow_results(store.connection)
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


def test_monte_carlo_quality_allows_unreported_energy_transport(tmp_path: Path) -> None:
    thresholds = QualityThresholds(
        mobility_rse=0.1,
        diffusion_rse=0.1,
        major_rate_rse=0.1,
    )
    with WorkflowStore(tmp_path / "mc_optional_transport.sqlite") as store:
        ensure_workflow_provenance(store.connection, thresholds)
        for replicate, scale in enumerate((1.0, 1.01, 0.99)):
            _insert_case(
                store.connection,
                replicate=replicate,
                mean_energy=5.0 * scale,
                reduced_mobility=12.0 * scale,
                reduced_diffusion_l=22.0 * scale,
                reduced_diffusion_t=33.0 * scale,
                reduced_energy_mobility=None,
                reduced_energy_diffusion=None,
            )
            _insert_rate(store.connection, replicate=replicate, value=1.0e-15 * scale)
            _insert_eedf(store.connection, replicate=replicate)
        aggregate_workflow_results(
            store.connection,
            thresholds,
        )
        quality = store.connection.execute(
            "SELECT passed, failure_reasons_json FROM aggregate_quality"
        ).fetchone()

    assert quality["passed"] == 1
    assert json.loads(quality["failure_reasons_json"]) == []


def test_required_reaction_rse_is_independent_of_elastic_major_rate(
    tmp_path: Path,
) -> None:
    thresholds = QualityThresholds(
        mobility_rse=0.1,
        diffusion_rse=0.1,
        major_rate_rse=0.1,
        major_rate_fraction=0.01,
    )
    with WorkflowStore(tmp_path / "required_reaction.sqlite") as store:
        ensure_workflow_provenance(store.connection, thresholds)
        for replicate, ionization in enumerate((1.0e-18, 4.0e-18, 1.0e-18)):
            _insert_case(store.connection, replicate=replicate)
            _insert_rate(
                store.connection,
                replicate=replicate,
                rate_index=0,
                process="elastic momentum",
                process_type="elastic",
                threshold_eV=0.0,
                value=1.0e-12,
            )
            _insert_rate(
                store.connection,
                replicate=replicate,
                rate_index=1,
                process="rare ionization",
                process_type="ionization",
                value=ionization,
            )
            _insert_eedf(store.connection, replicate=replicate)
        aggregate_workflow_results(
            store.connection,
            thresholds,
        )
        quality = store.connection.execute(
            "SELECT passed, failure_reasons_json, max_major_rate_rse "
            "FROM aggregate_quality"
        ).fetchone()

    reasons = json.loads(quality["failure_reasons_json"])
    assert quality["passed"] == 0
    assert quality["max_major_rate_rse"] == pytest.approx(0.0)
    assert any(
        reason.startswith("required_rate_rse_exceeds_threshold:Ar:rare ionization")
        for reason in reasons
    )
    assert not any(
        reason.startswith("major_rate_rse_exceeds_threshold:Ar:rare ionization")
        for reason in reasons
    )


def test_required_rate_relevance_uses_same_process_peak_across_en(
    tmp_path: Path,
) -> None:
    thresholds = QualityThresholds(
        mobility_rse=0.1,
        diffusion_rse=0.1,
        major_rate_rse=0.2,
        major_rate_fraction=0.01,
        required_rate_min_process_peak_fraction=1.0e-5,
        required_rate_rse=RequiredRateRse(
            process_type={"excitation": 0.2, "ionization": 0.2}
        ),
    )

    def triplet(mean: float, rse: float) -> tuple[float, float, float]:
        delta = mean * rse * math.sqrt(3.0)
        return mean - delta, mean, mean + delta

    excitation = {
        20.0: triplet(1.20e-6, 0.36),
        35.0: triplet(2.0e-6, 0.05),
        50.0: triplet(1.0, 0.10),
    }
    ionization = {
        20.0: (0.0, 0.0, 0.0),
        35.0: triplet(1.74e-6, 0.285),
        50.0: triplet(1.0, 0.10),
    }
    database = tmp_path / "rate_relevance.sqlite"
    with WorkflowStore(database) as store:
        ensure_workflow_provenance(store.connection, thresholds)
        for e_over_n in (20.0, 35.0, 50.0):
            for replicate, scale in enumerate((0.99, 1.0, 1.01)):
                _insert_case(
                    store.connection,
                    e_over_n=e_over_n,
                    replicate=replicate,
                    reduced_mobility=12.0 * scale,
                    reduced_diffusion_l=22.0 * scale,
                    reduced_diffusion_t=33.0 * scale,
                )
                _insert_rate(
                    store.connection,
                    e_over_n=e_over_n,
                    replicate=replicate,
                    rate_index=0,
                    process="elastic momentum",
                    process_type="elastic",
                    threshold_eV=0.0,
                    value=1.0 * scale,
                )
                _insert_rate(
                    store.connection,
                    e_over_n=e_over_n,
                    replicate=replicate,
                    rate_index=1,
                    process="excitation",
                    process_type="excitation",
                    threshold_eV=11.5,
                    value=excitation[e_over_n][replicate],
                )
                _insert_rate(
                    store.connection,
                    e_over_n=e_over_n,
                    replicate=replicate,
                    rate_index=2,
                    process="ionization",
                    process_type="ionization",
                    threshold_eV=15.0,
                    value=ionization[e_over_n][replicate],
                )
                _insert_eedf(
                    store.connection,
                    e_over_n=e_over_n,
                    replicate=replicate,
                )
        store.connection.commit()
        aggregate_workflow_results(store.connection, thresholds)
        quality = {
            float(row["e_over_n_Td"]): row
            for row in store.connection.execute(
                """
                SELECT e_over_n_Td, passed, failure_reasons_json
                FROM aggregate_quality ORDER BY e_over_n_Td
                """
            )
        }
        low_rate_stats = {
            (float(row["e_over_n_Td"]), str(row["process"])): row
            for row in store.connection.execute(
                """
                SELECT e_over_n_Td, process, mean, relative_standard_error,
                       estimate_status
                FROM aggregate_scalars
                WHERE scalar_group = 'rate'
                  AND process IN ('excitation', 'ionization')
                """
            )
        }

        assert all(quality[e]["passed"] == 1 for e in (20.0, 35.0, 50.0))
        assert low_rate_stats[(20.0, "excitation")][
            "relative_standard_error"
        ] == pytest.approx(0.36)
        assert low_rate_stats[(20.0, "ionization")]["estimate_status"] == (
            "censored_all_zero"
        )
        assert low_rate_stats[(35.0, "ionization")][
            "relative_standard_error"
        ] == pytest.approx(0.285)

        for replicate, value in enumerate(triplet(1.0, 0.30)):
            store.connection.execute(
                """
                UPDATE rates SET rate_coefficient_m3_s = ?
                WHERE solver = 'monte_carlo' AND mixture_id = 0
                  AND e_over_n_Td = 50.0 AND replicate = ?
                  AND process = 'excitation'
                """,
                (value, replicate),
            )
        store.connection.commit()
        aggregate_workflow_results(store.connection, thresholds)
        reevaluated = {
            float(row["e_over_n_Td"]): row
            for row in store.connection.execute(
                """
                SELECT e_over_n_Td, passed, failure_reasons_json
                FROM aggregate_quality ORDER BY e_over_n_Td
                """
            )
        }

    assert reevaluated[20.0]["passed"] == 1
    assert reevaluated[35.0]["passed"] == 1
    assert reevaluated[50.0]["passed"] == 0
    assert any(
        reason.startswith("required_rate_rse_exceeds_threshold:Ar:excitation")
        for reason in json.loads(reevaluated[50.0]["failure_reasons_json"])
    )


def test_required_process_name_gate_and_all_zero_censor_metadata(
    tmp_path: Path,
) -> None:
    thresholds = QualityThresholds(
        mobility_rse=0.1,
        diffusion_rse=0.1,
        major_rate_rse=0.1,
        required_rate_rse=RequiredRateRse(
            process_type={}, process={"NAMED TAIL CHANNEL": 0.2}
        ),
    )
    with WorkflowStore(tmp_path / "censored_reaction.sqlite") as store:
        ensure_workflow_provenance(store.connection, thresholds)
        for replicate in range(3):
            _insert_case(store.connection, replicate=replicate)
            _insert_rate(
                store.connection,
                replicate=replicate,
                process="named tail channel",
                process_type="other",
                value=0.0,
            )
            _insert_eedf(store.connection, replicate=replicate)
        aggregate_workflow_results(
            store.connection,
            thresholds,
        )
        scalar = store.connection.execute(
            """
            SELECT mean, sample_stddev, standard_error,
                   relative_standard_error, ci95_low, ci95_high,
                   ci95_critical_value, estimate_status,
                   uncertainty_available
            FROM aggregate_scalars
            WHERE scalar_group = 'rate'
            """
        ).fetchone()
        quality = store.connection.execute(
            "SELECT passed, failure_reasons_json, uncertainty_available "
            "FROM aggregate_quality"
        ).fetchone()

    assert scalar["mean"] == 0.0
    assert scalar["sample_stddev"] == 0.0
    assert scalar["standard_error"] is None
    assert scalar["relative_standard_error"] is None
    assert scalar["ci95_low"] is None
    assert scalar["ci95_high"] is None
    assert scalar["ci95_critical_value"] is None
    assert scalar["estimate_status"] == "censored_all_zero"
    assert scalar["uncertainty_available"] == 0
    assert quality["passed"] == 0
    assert quality["uncertainty_available"] == 0
    assert any(
        reason.startswith("required_rate_censored_all_zero:Ar:named tail channel")
        for reason in json.loads(quality["failure_reasons_json"])
    )


def test_conservative_rebin_preserves_probability_mass() -> None:
    rebinned = conservative_rebin_probability_mass(
        [0.0, 1.0, 3.0],
        [0.2, 0.8],
        [0.0, 2.0, 3.0],
    )

    assert rebinned == pytest.approx([0.6, 0.4])
    assert float(sum(rebinned)) == pytest.approx(1.0)


def test_nonuniform_eedf_cells_keep_exact_widths_during_aggregation(
    tmp_path: Path,
) -> None:
    database = tmp_path / "nonuniform.sqlite"
    with WorkflowStore(database) as store:
        for replicate in range(3):
            _insert_case(store.connection, replicate=replicate)
            _insert_rate(
                store.connection,
                replicate=replicate,
                value=1.0e-15,
            )
            store.connection.executemany(
                """
                INSERT INTO eedf_bins(
                    mixture_id, solver, e_over_n_Td, replicate, bin_index,
                    energy_eV, energy_width_eV, eedf, eepf, sample_count,
                    effective_sample_count, relative_standard_error
                ) VALUES (0, 'monte_carlo', 10.0, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [
                    (replicate, 0, 0.5, 1.0, 0.5, 0.5, 100, 100.0, 0.1),
                    (
                        replicate,
                        1,
                        1.5,
                        1.0,
                        0.3,
                        0.3 / 1.5**0.5,
                        80,
                        80.0,
                        1.0 / 80.0**0.5,
                    ),
                    (
                        replicate,
                        2,
                        7.0,
                        10.0,
                        0.02,
                        0.02 / 7.0**0.5,
                        40,
                        40.0,
                        1.0 / 40.0**0.5,
                    ),
                ],
            )
        store.connection.commit()

        aggregate_workflow_results(store.connection)
        rows = store.connection.execute(
            """
            SELECT energy_left_eV, energy_right_eV, energy_eV,
                   energy_width_eV, mean_probability_mass, eedf,
                   pooled_effective_sample_count
            FROM aggregate_eedf_bins
            ORDER BY bin_index
            """
        ).fetchall()

    expected_cells = (
        (0.0, 1.0, 0.5, 1.0),
        (1.0, 2.0, 1.5, 1.0),
        (2.0, 12.0, 7.0, 10.0),
    )
    for row, expected in zip(rows, expected_cells, strict=True):
        assert tuple(row[:4]) == pytest.approx(expected)
    assert [row[4] for row in rows] == pytest.approx([0.5, 0.3, 0.2])
    assert [row[5] for row in rows] == pytest.approx([0.5, 0.3, 0.02])
    assert [row[6] for row in rows] == pytest.approx([300.0, 240.0, 120.0])


def test_narrow_high_energy_cell_roundtrip_uses_scale_aware_tolerance(
    tmp_path: Path,
) -> None:
    database = tmp_path / "narrow_high_energy.sqlite"
    cells = (
        (88.735, 177.47, 0.9 / 177.47),
        (177.475, 0.01, 10.0),
    )
    with WorkflowStore(database) as store:
        for replicate in range(3):
            _insert_case(store.connection, replicate=replicate)
            _insert_rate(store.connection, replicate=replicate, value=1.0e-15)
            store.connection.executemany(
                """
                INSERT INTO eedf_bins(
                    mixture_id, solver, e_over_n_Td, replicate, bin_index,
                    energy_eV, energy_width_eV, eedf, eepf, sample_count,
                    effective_sample_count, relative_standard_error
                ) VALUES (0, 'monte_carlo', 10.0, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [
                    (
                        replicate,
                        index,
                        center,
                        width,
                        density,
                        density / center**0.5,
                        100,
                        100.0,
                        0.1,
                    )
                    for index, (center, width, density) in enumerate(cells)
                ],
            )
        store.connection.commit()

        aggregate_workflow_results(store.connection)
        rows = store.connection.execute(
            """
            SELECT energy_eV, energy_width_eV, mean_probability_mass
            FROM aggregate_eedf_bins ORDER BY bin_index
            """
        ).fetchall()

    assert [row[0] for row in rows] == pytest.approx([88.735, 177.475])
    assert [row[1] for row in rows] == pytest.approx([177.47, 0.01])
    assert [row[2] for row in rows] == pytest.approx([0.9, 0.1])


def test_aggregate_tables_are_independent_of_replicate_insert_order(
    tmp_path: Path,
) -> None:
    loose = QualityThresholds(
        mobility_rse=0.2,
        diffusion_rse=0.2,
        major_rate_rse=0.2,
    )
    conn_a = _build_synthetic_db(
        tmp_path / "a.sqlite", order=(0, 1, 2), quality=loose
    )
    conn_b = _build_synthetic_db(
        tmp_path / "b.sqlite", order=(2, 0, 1), quality=loose
    )
    aggregate_workflow_results(conn_a, loose)
    aggregate_workflow_results(conn_b, loose)

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
        aggregate_workflow_results(store.connection)
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
                    "required_rate_min_process_peak_fraction": 1.0e-5,
                    "required_rate_rse": {
                        "process_type": {
                            "ionization": 0.15,
                            "excitation": None,
                        },
                        "process": {"Ar resonance": 0.25},
                    },
                },
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    workflow = load_workflow(workflow_path)
    assert workflow.quality.mobility_rse == pytest.approx(0.1)
    assert workflow.quality.major_rate_fraction == pytest.approx(0.2)
    assert (
        workflow.quality.required_rate_min_process_peak_fraction
        == pytest.approx(1.0e-5)
    )
    assert workflow.quality.required_rate_rse.process_type == {
        "ionization": pytest.approx(0.15),
        "excitation": None,
    }
    assert workflow.quality.required_rate_rse.process == {
        "Ar resonance": pytest.approx(0.25)
    }

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

    invalid_relevance_path = tmp_path / "bad_required_rate_relevance.yaml"
    invalid_relevance_path.write_text(
        yaml.safe_dump(
            {
                "base_config": base_path.name,
                "database": "outputs/swarm.sqlite",
                "e_over_n_Td": [10.0],
                "mixtures": [{"Ar": 1.0}],
                "quality": {
                    "required_rate_min_process_peak_fraction": 1.01,
                },
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    with pytest.raises(
        ValueError,
        match="required_rate_min_process_peak_fraction",
    ):
        load_workflow(invalid_relevance_path)

    invalid_gate_path = tmp_path / "bad_required_gate.yaml"
    invalid_gate_path.write_text(
        yaml.safe_dump(
            {
                "base_config": base_path.name,
                "database": "outputs/swarm.sqlite",
                "e_over_n_Td": [10.0],
                "mixtures": [{"Ar": 1.0}],
                "quality": {
                    "required_rate_rse": {
                        "process_type": {"ionization": -0.1}
                    }
                },
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="required_rate_rse"):
        load_workflow(invalid_gate_path)


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
    quality_rows = list(
        csv.DictReader(
            (tmp_path / "aggregate" / "aggregate_quality.csv").open(
                encoding="utf-8", newline=""
            )
        )
    )
    assert quality_rows
    assert {row["quality_scope"] for row in quality_rows} == {
        "raw_replica_statistical_screen"
    }
    manifest = json.loads(
        (tmp_path / "aggregate" / "manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["quality_scope"]["aggregate_quality.csv"] == (
        "raw_replica_statistical_screen"
    )


def test_standalone_aggregate_rejects_missing_or_changed_quality_provenance(
    tmp_path: Path,
) -> None:
    missing_path = tmp_path / "missing.sqlite"
    missing = _build_synthetic_db(missing_path)
    missing.execute(
        "DELETE FROM metadata WHERE key = 'quality_thresholds_json'"
    )
    missing.commit()
    missing.close()

    with pytest.raises(WorkflowSchemaError, match="missing quality_thresholds_json"):
        aggregate_database(missing_path)

    mismatch_path = tmp_path / "mismatch.sqlite"
    mismatch = _build_synthetic_db(mismatch_path)
    mismatch.close()
    with pytest.raises(WorkflowSchemaError, match="differ from workflow database"):
        aggregate_database(
            mismatch_path,
            thresholds=QualityThresholds(mobility_rse=0.1),
        )
