from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest
import yaml

from swarm_workflow.sweep import load_workflow, run_sweep
from swarm_workflow.cli import main as workflow_cli_main
from product_helpers import base_product_config, write_config


def _write_mixture_xs(tmp_path: Path) -> Path:
    path = tmp_path / "ar_o2_minimal.csv"
    rows = [
        "species,process,type,threshold_eV,mass_amu,energy_eV,cross_section_m2",
    ]
    for species, mass, scale in [("Ar", 39.948, 1.0), ("O2", 31.998, 1.4)]:
        for energy in [0.0, 1.0, 10.0, 100.0, 1000.0]:
            rows.append(
                f"{species},momentum,momentum,0.0,{mass},{energy},{scale * 1.0e-20}"
            )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def _write_workflow_base(tmp_path: Path) -> Path:
    xs_path = _write_mixture_xs(tmp_path)
    data = base_product_config(tmp_path, ["two_term"])
    data["run"]["e_over_n_Td"] = [50.0]
    data["conditions"]["gas_mixture"] = [
        {"species": "Ar", "fraction": 1.0, "mass_amu": 39.948},
        {"species": "O2", "fraction": 0.0, "mass_amu": 31.998},
    ]
    data["cross_sections"]["files"] = [
        {"path": xs_path.as_posix(), "format": "csv"}
    ]
    return write_config(tmp_path, data, "base.yaml")


def _write_workflow(tmp_path: Path, base_path: Path) -> Path:
    path = tmp_path / "workflow.yaml"
    path.write_text(
        yaml.safe_dump(
            {
                "base_config": base_path.name,
                "database": "outputs/swarm.sqlite",
                "e_over_n_Td": [10.0, 20.0],
                "mixtures": [
                    {"Ar": 1.0, "O2": 0.0},
                    {"Ar": 0.99, "O2": 0.01},
                ],
                "mc": {"replicas": 1, "base_seed": 12345},
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    return path


def test_two_term_workflow_sweep_resumes_without_duplicate_rows(tmp_path: Path) -> None:
    base_path = _write_workflow_base(tmp_path)
    workflow_path = _write_workflow(tmp_path, base_path)

    first = run_sweep(workflow_path)
    second = run_sweep(workflow_path)
    workflow_cli_main(["sweep", str(workflow_path)])

    assert first.cases_written == 4
    assert second.cases_written == 4
    with sqlite3.connect(tmp_path / "outputs" / "swarm.sqlite") as conn:
        assert conn.execute("SELECT COUNT(*) FROM mixtures").fetchone()[0] == 2
        assert conn.execute("SELECT COUNT(*) FROM cases").fetchone()[0] == 4
        assert conn.execute("SELECT COUNT(*) FROM metadata").fetchone()[0] >= 3
        assert {
            row[0]
            for row in conn.execute(
                "SELECT DISTINCT solver FROM cases ORDER BY solver"
            ).fetchall()
        } == {"two_term"}
        norms = conn.execute(
            """
            SELECT mixture_id, solver, e_over_n_Td, replicate,
                   SUM(eedf * energy_width_eV)
            FROM eedf_bins
            GROUP BY mixture_id, solver, e_over_n_Td, replicate
            ORDER BY mixture_id, e_over_n_Td
            """
        ).fetchall()

    assert len(norms) == 4
    for *_case_key, norm in norms:
        assert norm == pytest.approx(1.0, rel=1.0e-8, abs=1.0e-10)


def test_workflow_rejects_species_not_present_in_base_config(tmp_path: Path) -> None:
    base_path = _write_workflow_base(tmp_path)
    workflow_path = tmp_path / "bad_workflow.yaml"
    workflow_path.write_text(
        yaml.safe_dump(
            {
                "base_config": base_path.name,
                "database": "outputs/swarm.sqlite",
                "e_over_n_Td": [10.0],
                "mixtures": [{"Ar": 1.0, "N2": 0.1}],
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="unknown species"):
        load_workflow(workflow_path)


def test_electron_swarm_does_not_import_workflow_or_sqlite() -> None:
    root = Path(__file__).resolve().parents[1]
    for path in (root / "electron_swarm").rglob("*.py"):
        text = path.read_text(encoding="utf-8")
        assert "swarm_workflow" not in text
        assert "sqlite3" not in text
