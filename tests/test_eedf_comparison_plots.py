from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sqlite3

import numpy as np
import pytest

from swarm_workflow.plots.eedf_comparison import (
    EedfCase,
    EedfComparisonError,
    compare_eedf_cases,
    generate_eedf_solver_comparison,
    load_mc_eedf_database,
    load_table_eedf,
)


_PHYSICAL_CONTEXT = {
    "schema": "test.swarm_physical_context.v1",
    "field": {"type": "dc"},
}
_MIXTURE = {
    "mixture_id": 0,
    "species": [{"species": "Ar", "fraction": 1.0, "mass_amu": 39.948}],
}
_CROSS_SECTIONS_SHA256 = "a" * 64
_MC_IDENTITY = {
    "workflow_config_sha256": "b" * 64,
    "base_config_sha256": "c" * 64,
    "mc_sampling_plan_json": "[]",
    "mc_transport_estimator_schema_version": "test_transport.v1",
    "mc_eedf_estimator_schema_version": "test_eedf.v1",
    "mc_tail_estimator_schema_version": "test_tail.v1",
    "mc_seed_derivation_schema_version": "test_seed.v1",
    "mc_solver_source_sha256": "d" * 64,
}
_MC_QUALIFICATION_PROFILE = "function_eedf_restricted_lmea"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _case(solver: str, density: tuple[float, float]) -> EedfCase:
    return EedfCase(
        solver=solver,
        e_over_n_td=10.0,
        edges_eV=np.asarray([0.0, 1.0, 3.0]),
        density_eV_inv=np.asarray(density, dtype=float),
        reported_mean_energy_eV=1.0,
    )


def test_eedf_metrics_are_bounded_and_identical_case_is_zero() -> None:
    reference = _case("two_term", (0.5, 0.25))

    metrics = compare_eedf_cases(reference, reference, tail_thresholds_eV=(1.0,))

    assert metrics["total_variation"] == pytest.approx(0.0, abs=1.0e-15)
    assert metrics["shape_total_variation"] == pytest.approx(0.0, abs=1.0e-15)
    assert metrics["hellinger_distance"] == pytest.approx(0.0, abs=1.0e-15)
    assert metrics["jensen_shannon_divergence_nats"] == pytest.approx(0.0, abs=1.0e-15)
    assert metrics["wasserstein_1_eV"] == pytest.approx(0.0, abs=1.0e-15)


def test_eedf_metrics_use_conservative_union_grid() -> None:
    reference = EedfCase(
        solver="two_term",
        e_over_n_td=10.0,
        edges_eV=np.asarray([0.0, 1.0, 2.0]),
        density_eV_inv=np.asarray([1.0, 0.0]),
        reported_mean_energy_eV=0.5,
    )
    candidate = EedfCase(
        solver="propagator",
        e_over_n_td=10.0,
        edges_eV=np.asarray([0.0, 0.5, 1.0, 2.0]),
        density_eV_inv=np.asarray([1.0, 1.0, 0.0]),
        reported_mean_energy_eV=0.5,
    )

    metrics = compare_eedf_cases(reference, candidate, tail_thresholds_eV=())

    assert metrics["total_variation"] == pytest.approx(0.0, abs=1.0e-15)
    assert metrics["shape_total_variation"] == pytest.approx(0.0, abs=1.0e-15)


def test_eedf_jensen_shannon_keeps_subnormal_tail_mass_finite() -> None:
    subnormal = float(np.nextafter(0.0, 1.0))
    reference = EedfCase(
        solver="monte_carlo",
        e_over_n_td=1000.0,
        edges_eV=np.asarray([0.0, 1.0, 2.0]),
        density_eV_inv=np.asarray([subnormal, 1.0]),
        reported_mean_energy_eV=1.5,
    )
    candidate = EedfCase(
        solver="two_term",
        e_over_n_td=1000.0,
        edges_eV=np.asarray([0.0, 1.0, 2.0]),
        density_eV_inv=np.asarray([0.0, 1.0]),
        reported_mean_energy_eV=1.5,
    )

    divergence = compare_eedf_cases(
        reference,
        candidate,
        tail_thresholds_eV=(),
    )["jensen_shannon_divergence_nats"]

    assert np.isfinite(divergence)
    assert 0.0 <= divergence <= np.log(2.0)


def _write_table(
    root: Path,
    solver_shift: float = 0.0,
    *,
    solver: str = "two_term",
    physical_context: dict[str, object] | None = None,
) -> None:
    root.mkdir(parents=True)
    with (root / "mean_energy_vs_en.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("E_over_N_Td", "mean_energy_eV"))
        writer.writerow((1, 0.5 + solver_shift))
        writer.writerow((5, 0.5 + solver_shift))
    with (root / "eedf.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            (
                "electron_energy_eV",
                "energy_width_eV",
                "E_over_N_Td",
                "mean_energy_eV",
                "eedf",
            )
        )
        for field in (1, 5):
            writer.writerow((0.25, 0.5, field, 0.5 + solver_shift, 1.0))
            writer.writerow((0.75, 0.5, field, 0.5 + solver_shift, 1.0))
    manifest = {
        "source": solver,
        "physical_context": physical_context or _PHYSICAL_CONTEXT,
        "mixture": _MIXTURE,
        "hashes": {"cross_sections_sha256": _CROSS_SECTIONS_SHA256},
        "tables": {
            name: {"sha256": _sha256(root / name)}
            for name in ("eedf.csv", "mean_energy_vs_en.csv")
        },
    }
    (root / "manifest.json").write_text(
        json.dumps(manifest, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _write_mc_database(path: Path) -> None:
    connection = sqlite3.connect(path)
    connection.executescript(
        """
        CREATE TABLE aggregate_scalars (
            mixture_id INTEGER, solver TEXT, e_over_n_Td REAL,
            scalar_group TEXT, scalar_name TEXT, mean REAL
        );
        CREATE TABLE aggregate_eedf_bins (
            mixture_id INTEGER, solver TEXT, e_over_n_Td REAL,
            bin_index INTEGER, energy_eV REAL, energy_width_eV REAL,
            eedf REAL, probability_mass_standard_error REAL,
            valid_replicates INTEGER, uncertainty_available INTEGER
        );
        CREATE TABLE metadata (key TEXT, value TEXT);
        CREATE TABLE mixtures (mixture_id INTEGER, fractions_json TEXT);
        CREATE TABLE mixture_species (
            mixture_id INTEGER, species TEXT, fraction REAL, mass_amu REAL
        );
        """
    )
    connection.execute(
        "INSERT INTO metadata VALUES (?, ?)",
        ("physical_context_json", json.dumps(_PHYSICAL_CONTEXT, sort_keys=True)),
    )
    connection.execute(
        "INSERT INTO metadata VALUES (?, ?)",
        ("cross_sections_sha256", _CROSS_SECTIONS_SHA256),
    )
    connection.executemany(
        "INSERT INTO metadata VALUES (?, ?)",
        _MC_IDENTITY.items(),
    )
    connection.execute(
        "INSERT INTO mixtures VALUES (?, ?)",
        (0, json.dumps({"Ar": 1.0}, sort_keys=True)),
    )
    connection.execute(
        "INSERT INTO mixture_species VALUES (?, ?, ?, ?)",
        (0, "Ar", 1.0, 39.948),
    )
    for field in (1.0, 5.0):
        connection.execute(
            "INSERT INTO aggregate_scalars VALUES (?, ?, ?, ?, ?, ?)",
            (0, "monte_carlo", field, "case", "mean_energy_eV", 0.5),
        )
        density = (1.2, 0.8) if field == 1.0 else (1.0, 1.0)
        for index, (center, value) in enumerate(zip((0.25, 0.75), density)):
            connection.execute(
                "INSERT INTO aggregate_eedf_bins VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (0, "monte_carlo", field, index, center, 0.5, value, 0.01, 4, 1),
            )
    connection.commit()
    connection.close()


def _write_mc_manifest(root: Path, qualification: Path, database: Path) -> None:
    payload = {
        "source": "monte_carlo",
        "physical_context": _PHYSICAL_CONTEXT,
        "mixture": _MIXTURE,
        "hashes": {
            "cross_sections_sha256": _CROSS_SECTIONS_SHA256,
            "source_database_sha256": _sha256(database),
            **_MC_IDENTITY,
        },
        "source_policy": {
            "qualification_profile": _MC_QUALIFICATION_PROFILE,
        },
        "mc_qualification": {"file": qualification.name},
        "tables": {
            qualification.name: {"sha256": _sha256(qualification)},
        },
    }
    (root / "manifest.json").write_text(
        json.dumps(payload, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def test_table_and_mc_loaders_preserve_normalization_and_uncertainty(
    tmp_path: Path,
) -> None:
    table = tmp_path / "table"
    database = tmp_path / "mc.sqlite"
    _write_table(table)
    _write_mc_database(database)

    deterministic = load_table_eedf(table, "two_term")
    mc = load_mc_eedf_database(database)

    assert deterministic.cases[1.0].source_normalization == pytest.approx(1.0)
    assert mc.cases[1.0].source_normalization == pytest.approx(1.0)
    assert mc.cases[1.0].ci95_half_density_eV_inv is not None
    assert np.all(mc.cases[1.0].ci95_half_density_eV_inv > 0.0)


def test_mc_loader_rejects_database_changed_after_table_materialization(
    tmp_path: Path,
) -> None:
    database = tmp_path / "mc.sqlite"
    qualification = tmp_path / "qualification.csv"
    _write_mc_database(database)
    qualification.write_text(
        "E_over_N_Td,active_closure_quality_passed,qualification_profile\n"
        f"1,1,{_MC_QUALIFICATION_PROFILE}\n"
        f"5,1,{_MC_QUALIFICATION_PROFILE}\n",
        encoding="utf-8",
    )
    _write_mc_manifest(tmp_path, qualification, database)
    with sqlite3.connect(database) as connection:
        connection.execute(
            "UPDATE metadata SET value = ? WHERE key = ?",
            ("changed", "base_config_sha256"),
        )

    with pytest.raises(EedfComparisonError, match="does not match its source manifest"):
        load_mc_eedf_database(
            database,
            source_manifest_path=tmp_path / "manifest.json",
        )


def test_comparison_keeps_raw_mc_when_downstream_gate_does_not_pass(
    tmp_path: Path,
) -> None:
    pytest.importorskip("matplotlib")
    two_term = tmp_path / "two_term"
    propagator = tmp_path / "propagator"
    database = tmp_path / "mc.sqlite"
    qualification = tmp_path / "qualification.csv"
    output = tmp_path / "plots"
    _write_table(two_term)
    _write_table(propagator, solver="propagator")
    _write_mc_database(database)
    with qualification.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            (
                "E_over_N_Td",
                "active_closure_quality_passed",
                "qualification_profile",
            )
        )
        writer.writerow((1, 0, _MC_QUALIFICATION_PROFILE))
        writer.writerow((5, 1, _MC_QUALIFICATION_PROFILE))
    _write_mc_manifest(qualification.parent, qualification, database)

    summary = generate_eedf_solver_comparison(
        two_term_directory=two_term,
        propagator_directory=propagator,
        mc_database=database,
        mc_qualification_csv=qualification,
        output_directory=output,
        representative_fields_td=(1.0, 5.0),
        tail_thresholds_eV=(0.5,),
    )

    assert summary.mc_status_not_passed_fields_td == (1.0,)
    assert all(path.is_file() and path.stat().st_size > 0 for path in summary.figures)
    assert summary.metrics_csv.is_file()
    manifest = json.loads(summary.manifest.read_text(encoding="utf-8"))
    assert manifest["format_version"] == 2
    assert manifest["comparison"] == "raw_two_term_propagator_and_raw_monte_carlo"
    assert manifest["comparison_scope"]["mode"] == "independent_raw_solver_outputs"
    assert manifest["mc_status_annotation"]["not_passed_fields_Td"] == [1.0]
    status = json.loads(summary.mc_status.read_text(encoding="utf-8"))
    assert status["comparison_invariant"].startswith(
        "EEDF figures and metrics always use"
    )
    assert status["not_passed_fields_Td"] == [1.0]
    assert status["selection_policy"].startswith("not inferred here")
    assert "downstream_selected_source" not in status["per_anchor"][0]
    rows = list(csv.DictReader(summary.metrics_csv.open(encoding="utf-8")))
    assert {row["comparison"] for row in rows} == {
        "raw_mc_vs_two_term",
        "propagator_vs_two_term",
        "raw_mc_vs_propagator",
    }
    assert all("selected_source" not in row for row in rows)
    assert all("mc_qualified_for_active_closure" not in row for row in rows)
    raw_mc = next(
        row
        for row in rows
        if row["E_over_N_Td"] == "1.0"
        and row["comparison"] == "raw_mc_vs_two_term"
    )
    assert float(raw_mc["total_variation"]) == pytest.approx(0.1)
    for path in summary.figures:
        if path.suffix == ".svg":
            svg = path.read_text(encoding="utf-8").lower()
            assert "selected mc" not in svg
            assert "2t substitute" not in svg
            assert "raw unqualified mc" not in svg


def test_all_passed_report_still_compares_three_raw_solver_outputs(
    tmp_path: Path,
) -> None:
    pytest.importorskip("matplotlib")
    two_term = tmp_path / "two_term"
    propagator = tmp_path / "propagator"
    database = tmp_path / "mc.sqlite"
    qualification = tmp_path / "qualification.csv"
    output = tmp_path / "plots"
    _write_table(two_term)
    _write_table(propagator, solver="propagator")
    _write_mc_database(database)
    with qualification.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            (
                "E_over_N_Td",
                "active_closure_quality_passed",
                "qualification_profile",
            )
        )
        writer.writerow((1, 1, _MC_QUALIFICATION_PROFILE))
        writer.writerow((5, 1, _MC_QUALIFICATION_PROFILE))
    _write_mc_manifest(qualification.parent, qualification, database)

    summary = generate_eedf_solver_comparison(
        two_term_directory=two_term,
        propagator_directory=propagator,
        mc_database=database,
        mc_qualification_csv=qualification,
        output_directory=output,
        representative_fields_td=(1.0, 5.0),
        tail_thresholds_eV=(0.5,),
    )

    manifest = json.loads(summary.manifest.read_text(encoding="utf-8"))
    status = json.loads(summary.mc_status.read_text(encoding="utf-8"))
    assert summary.mc_status_not_passed_fields_td == ()
    assert manifest["comparison_scope"]["mode"] == "independent_raw_solver_outputs"
    assert status["not_passed_fields_Td"] == []
    assert status["passed_fields_Td"] == [1.0, 5.0]
    for path in summary.figures:
        if path.suffix == ".svg":
            svg = path.read_text(encoding="utf-8").lower()
            assert "fallback" not in svg
            assert "selected" not in svg


def test_comparison_rejects_mismatched_physical_context(tmp_path: Path) -> None:
    two_term = tmp_path / "two_term"
    propagator = tmp_path / "propagator"
    database = tmp_path / "mc.sqlite"
    qualification = tmp_path / "qualification.csv"
    _write_table(two_term)
    _write_table(
        propagator,
        solver="propagator",
        physical_context={
            "schema": "test.swarm_physical_context.v1",
            "field": {"type": "rf", "frequency_Hz": 13.56e6},
        },
    )
    _write_mc_database(database)
    with qualification.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            (
                "E_over_N_Td",
                "active_closure_quality_passed",
                "qualification_profile",
            )
        )
        writer.writerow((1, 1, _MC_QUALIFICATION_PROFILE))
        writer.writerow((5, 1, _MC_QUALIFICATION_PROFILE))
    _write_mc_manifest(qualification.parent, qualification, database)

    with pytest.raises(EedfComparisonError, match="mismatched physical_context"):
        generate_eedf_solver_comparison(
            two_term_directory=two_term,
            propagator_directory=propagator,
            mc_database=database,
            mc_qualification_csv=qualification,
            output_directory=tmp_path / "plots",
            representative_fields_td=(1.0, 5.0),
        )
    assert not (tmp_path / "plots").exists()
