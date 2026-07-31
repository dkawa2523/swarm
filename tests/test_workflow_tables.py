from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import pytest

from swarm_workflow.tables import TableBuildError, build_tables
from swarm_workflow.store import WorkflowStore
from product_helpers import (
    insert_workflow_case,
    insert_workflow_eedf,
    insert_workflow_rate,
)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as fp:
        return list(csv.DictReader(fp))


def _insert_case_set(
    store: WorkflowStore,
    *,
    solver: str,
    mixture_id: int = 0,
    e_values: tuple[float, ...] = (10.0, 20.0),
    replicate: int = 0,
    mean_values: tuple[float, ...] = (1.0, 2.0),
    mobility_values: tuple[float, ...] = (10.0, 11.0),
) -> None:
    for index, e_over_n in enumerate(e_values):
        insert_workflow_case(
            store.connection,
            mixture_id=mixture_id,
            solver=solver,
            e_over_n=e_over_n,
            replicate=replicate,
            mean_energy=mean_values[index],
            drift=5.0,
            reduced_mobility=mobility_values[index],
            reduced_diffusion_l=20.0 + index,
            reduced_diffusion_t=30.0 + index,
            effective_townsend=1.0e-16,
        )
        insert_workflow_rate(
            store.connection,
            mixture_id=mixture_id,
            solver=solver,
            e_over_n=e_over_n,
            replicate=replicate,
            value=1.0e-15 * (index + 1),
            target_fraction=0.5,
            process="ionization",
            process_type="ionization",
            threshold_eV=15.0,
        )
        insert_workflow_eedf(
            store.connection,
            mixture_id=mixture_id,
            solver=solver,
            e_over_n=e_over_n,
            replicate=replicate,
        )


def test_build_tables_two_term_writes_two_mixture_directories(tmp_path: Path) -> None:
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        _insert_case_set(store, solver="two_term", mixture_id=0)
        _insert_case_set(store, solver="two_term", mixture_id=1)
        store.connection.commit()

    summary = build_tables(
        tmp_path / "workflow.sqlite",
        tmp_path / "tables",
        source="two_term",
    )

    assert summary.mixtures == 2
    root_manifest = json.loads((tmp_path / "tables" / "manifest.json").read_text())
    assert [item["mixture_id"] for item in root_manifest["mixtures"]] == [0, 1]
    for mixture_id in [0, 1]:
        mixture_dir = tmp_path / "tables" / f"mixture_{mixture_id:04d}"
        rows = _read_csv(mixture_dir / "transport_vs_en.csv")
        assert len(rows) == 2
        assert float(rows[0]["E_over_N_V_m2"]) == pytest.approx(
            float(rows[0]["E_over_N_Td"]) * 1.0e-21
        )
        for table_name in ["transport_vs_en.csv", "rates_vs_en.csv"]:
            for row in _read_csv(mixture_dir / table_name):
                for value in row.values():
                    if value:
                        number = float(value) if _looks_numeric(value) else None
                        if number is not None:
                            assert math.isfinite(number)
        manifest = json.loads((mixture_dir / "manifest.json").read_text())
        assert manifest["table_argument"]["secondary"] == "mean_energy_eV"


def test_two_term_quality_includes_and_enforces_solver_diagnostics(
    tmp_path: Path,
) -> None:
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        _insert_case_set(
            store,
            solver="two_term",
            e_values=(1.0,),
            mean_values=(1.0,),
            mobility_values=(10.0,),
        )
        diagnostics = {
            "two_term": {
                "converged": False,
                "iterations": 600,
                "residual_L1": 2.0e-7,
                "residual_tolerance": 1.0e-7,
                "tail_probability": 1.0e-20,
                "tail_probability_target": 1.0e-9,
                "edge_to_peak": 1.0e-30,
                "edge_to_peak_target": 1.0e-10,
                "grid_max_eV": 100.0,
                "grid_max_limit_eV": 20000.0,
            }
        }
        store.connection.execute(
            """
            UPDATE cases
            SET diagnostics_json = ?
            WHERE solver = 'two_term' AND e_over_n_Td = 1.0
            """,
            (json.dumps(diagnostics),),
        )
        store.connection.commit()

    build_tables(
        tmp_path / "workflow.sqlite",
        tmp_path / "tables",
        source="two_term",
    )
    [quality] = _read_csv(
        tmp_path / "tables" / "mixture_0000" / "quality.csv"
    )
    assert quality["passed"] == "0"
    assert quality["solver_diagnostics_available"] == "1"
    assert quality["solver_converged"] == "0"
    assert quality["solver_diagnostics_passed"] == "0"
    assert "solver_not_converged" in json.loads(
        quality["failure_reasons_json"]
    )
    manifest = json.loads(
        (
            tmp_path / "tables" / "mixture_0000" / "manifest.json"
        ).read_text()
    )
    assert manifest["quality_summary"]["passed"] is False


def test_build_tables_monte_carlo_uses_quality_passing_points_only(
    tmp_path: Path,
) -> None:
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        for replicate, mobility_offset in enumerate([0.0, 0.1]):
            _insert_case_set(
                store,
                solver="monte_carlo",
                e_values=(10.0,),
                replicate=replicate,
                mobility_values=(10.0 + mobility_offset,),
                mean_values=(1.0 + 0.01 * replicate,),
            )
        for replicate, mobility in enumerate([1.0, 30.0]):
            _insert_case_set(
                store,
                solver="monte_carlo",
                e_values=(20.0,),
                replicate=replicate,
                mobility_values=(mobility,),
                mean_values=(2.0,),
            )
        store.connection.commit()

    build_tables(tmp_path / "workflow.sqlite", tmp_path / "tables", source="monte_carlo")

    rows = _read_csv(tmp_path / "tables" / "mixture_0000" / "transport_vs_en.csv")
    assert [float(row["E_over_N_Td"]) for row in rows] == [10.0]


def test_build_tables_hybrid_log_correction_and_net_townsend(
    tmp_path: Path,
) -> None:
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        for e_over_n in [10.0, 20.0]:
            insert_workflow_case(
                store.connection,
                solver="two_term",
                e_over_n=e_over_n,
                mean_energy=e_over_n / 10.0,
                drift=10.0,
                reduced_mobility=10.0,
                effective_townsend=999.0,
            )
            insert_workflow_rate(
                store.connection,
                solver="two_term",
                e_over_n=e_over_n,
                value=2.0e-15,
                target_fraction=1.0,
                process="ionization",
                process_type="ionization",
                threshold_eV=15.0,
                rate_index=0,
            )
            insert_workflow_rate(
                store.connection,
                solver="two_term",
                e_over_n=e_over_n,
                value=1.0e-15,
                target_fraction=1.0,
                process="attachment",
                process_type="attachment",
                threshold_eV=0.0,
                rate_index=1,
            )
            insert_workflow_eedf(store.connection, solver="two_term", e_over_n=e_over_n)
            for replicate, scale in enumerate([1.0, 1.01]):
                insert_workflow_case(
                    store.connection,
                    solver="monte_carlo",
                    e_over_n=e_over_n,
                    replicate=replicate,
                    mean_energy=e_over_n / 10.0 * scale,
                    drift=10.0,
                    reduced_mobility=12.0 * scale,
                    effective_townsend=1.0e-16,
                )
                insert_workflow_rate(
                    store.connection,
                    solver="monte_carlo",
                    e_over_n=e_over_n,
                    replicate=replicate,
                    value=4.0e-15 * scale,
                    target_fraction=1.0,
                    process="ionization",
                    process_type="ionization",
                    threshold_eV=15.0,
                    rate_index=0,
                )
                insert_workflow_rate(
                    store.connection,
                    solver="monte_carlo",
                    e_over_n=e_over_n,
                    replicate=replicate,
                    value=0.5e-15 * scale,
                    target_fraction=1.0,
                    process="attachment",
                    process_type="attachment",
                    threshold_eV=0.0,
                    rate_index=1,
                )
                insert_workflow_eedf(
                    store.connection,
                    solver="monte_carlo",
                    e_over_n=e_over_n,
                    replicate=replicate,
                )
        store.connection.commit()

    build_tables(tmp_path / "workflow.sqlite", tmp_path / "tables", source="hybrid")

    mixture_dir = tmp_path / "tables" / "mixture_0000"
    for row in _read_csv(mixture_dir / "transport_vs_en.csv"):
        assert float(row["reduced_mobility_m2_V_s_m3"]) > 0.0
        assert math.isfinite(float(row["reduced_mobility_m2_V_s_m3"]))
    townsend = _read_csv(mixture_dir / "townsend_vs_en.csv")
    assert all(float(row["effective_townsend_m2"]) != 999.0 for row in townsend)
    assert all(float(row["effective_townsend_m2"]) > 0.0 for row in townsend)


def test_build_tables_nonmonotonic_mean_energy_skips_mean_energy_tables(
    tmp_path: Path,
) -> None:
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        _insert_case_set(
            store,
            solver="two_term",
            e_values=(10.0, 20.0, 30.0),
            mean_values=(1.0, 3.0, 2.0),
            mobility_values=(10.0, 11.0, 12.0),
        )
        store.connection.commit()

    build_tables(tmp_path / "workflow.sqlite", tmp_path / "tables", source="two_term")

    mixture_dir = tmp_path / "tables" / "mixture_0000"
    manifest = json.loads((mixture_dir / "manifest.json").read_text())
    assert manifest["table_argument"]["secondary"] is None
    assert manifest["monotonicity"]["reason"] == "mean_energy_not_strictly_monotonic"
    assert not (mixture_dir / "transport_vs_mean_energy.csv").exists()


def test_build_tables_monte_carlo_fails_without_quality_passing_points(
    tmp_path: Path,
) -> None:
    with WorkflowStore(tmp_path / "workflow.sqlite") as store:
        _insert_case_set(store, solver="monte_carlo", e_values=(10.0,))
        store.connection.commit()

    with pytest.raises(TableBuildError, match="no quality-passing"):
        build_tables(
            tmp_path / "workflow.sqlite",
            tmp_path / "tables",
            source="monte_carlo",
        )


def _looks_numeric(value: str) -> bool:
    try:
        float(value)
    except ValueError:
        return False
    return True
