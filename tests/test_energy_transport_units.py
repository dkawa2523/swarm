from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from electron_swarm.core.transport import ElectronTransport
from electron_swarm.io.writers import SUMMARY_COLUMNS
from swarm_workflow.aggregate import aggregate_database
from swarm_workflow.store import WorkflowSchemaError, WorkflowStore
from swarm_workflow.tables import CASE_COLUMNS, UNITS, build_tables


CANONICAL_REDUCED_COLUMNS = {
    "reduced_electron_energy_mobility_m2_V_s_m3",
    "reduced_electron_energy_diffusion_m2_s_m3",
}
LEGACY_REDUCED_COLUMNS = {
    "reduced_electron_energy_mobility_eV_m2_V_s_m3",
    "reduced_electron_energy_diffusion_eV_m2_s_m3",
}
LEGACY_ACTUAL_PROPERTIES = {
    "electron_energy_mobility_eV_m2_V_s",
    "electron_energy_diffusion_eV_m2_s",
}


def test_energy_transport_fields_have_dimensionally_correct_names() -> None:
    density = 2.5e20
    transport = ElectronTransport.from_actual(
        definition="flux",
        gas_number_density_m3=density,
        drift_velocity_m_s=1.0,
        mobility_m2_V_s=2.0,
        diffusion_L_m2_s=3.0,
        diffusion_T_m2_s=4.0,
        electron_energy_mobility_m2_V_s=5.0,
        electron_energy_diffusion_m2_s=6.0,
    )

    assert transport.reduced_electron_energy_mobility_m2_V_s_m3 == pytest.approx(
        5.0 * density
    )
    assert transport.reduced_electron_energy_diffusion_m2_s_m3 == pytest.approx(
        6.0 * density
    )
    assert transport.electron_energy_mobility_m2_V_s == pytest.approx(5.0)
    assert transport.electron_energy_diffusion_m2_s == pytest.approx(6.0)
    for legacy in LEGACY_REDUCED_COLUMNS | LEGACY_ACTUAL_PROPERTIES:
        assert not hasattr(transport, legacy)


def test_public_columns_and_table_units_exclude_obsolete_ev_factor() -> None:
    assert CANONICAL_REDUCED_COLUMNS <= set(SUMMARY_COLUMNS)
    assert CANONICAL_REDUCED_COLUMNS <= set(CASE_COLUMNS)
    assert LEGACY_REDUCED_COLUMNS.isdisjoint(SUMMARY_COLUMNS)
    assert LEGACY_REDUCED_COLUMNS.isdisjoint(CASE_COLUMNS)
    assert (
        UNITS["reduced_electron_energy_mobility_m2_V_s_m3"]
        == "1/(V m s)"
    )
    assert UNITS["reduced_electron_energy_diffusion_m2_s_m3"] == "1/(m s)"
    assert UNITS["reduced_mobility_m2_V_s_m3"] == "1/(V m s)"
    assert UNITS["reduced_diffusion_L_m2_s_m3"] == "1/(m s)"
    assert UNITS["reduced_diffusion_T_m2_s_m3"] == "1/(m s)"
    assert all("eV" not in UNITS[column] for column in CANONICAL_REDUCED_COLUMNS)


def test_new_workflow_database_uses_only_canonical_columns(tmp_path: Path) -> None:
    path = tmp_path / "workflow.sqlite"
    with WorkflowStore(path) as store:
        columns = {
            str(row[1])
            for row in store.connection.execute("PRAGMA table_info(cases)")
        }

    assert CANONICAL_REDUCED_COLUMNS <= columns
    assert LEGACY_REDUCED_COLUMNS.isdisjoint(columns)


@pytest.mark.parametrize("entrypoint", ["store", "aggregate", "tables"])
def test_legacy_energy_transport_database_fails_fast(
    tmp_path: Path,
    entrypoint: str,
) -> None:
    path = tmp_path / "legacy.sqlite"
    connection = sqlite3.connect(path)
    connection.execute(
        """
        CREATE TABLE cases (
            reduced_electron_energy_mobility_eV_m2_V_s_m3 REAL,
            reduced_electron_energy_diffusion_eV_m2_s_m3 REAL
        )
        """
    )
    connection.commit()
    connection.close()

    with pytest.raises(
        WorkflowSchemaError,
        match="obsolete energy-transport columns.*regenerate",
    ):
        if entrypoint == "store":
            WorkflowStore(path)
        elif entrypoint == "aggregate":
            aggregate_database(path)
        else:
            build_tables(path, tmp_path / "tables", source="two_term")
