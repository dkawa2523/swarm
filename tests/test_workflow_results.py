"""Case data survives persistence, resume, and publication without RF loss."""

from copy import deepcopy
from dataclasses import asdict
import sqlite3

import numpy as np
import pandas as pd
import pytest

from electron_swarm import run
from electron_swarm.core.output import OutputConfig
from electron_swarm.core.results import RfPhaseResult, SwarmRunResult
from electron_swarm.io import read_product_run, write_outputs
from product_helpers import base_product_config, write_config
from swarm_workflow.campaign.document import load_case_template
from swarm_workflow.campaign.store import WorkflowSchemaError, WorkflowStore
from swarm_workflow.campaign.sweep import run_sweep


@pytest.fixture
def phase_case(tmp_path):
    config = load_case_template(write_config(tmp_path, base_product_config(tmp_path)))
    case = run(config, write=False).cases[0]
    case.rf_phase = RfPhaseResult(
        phase_fraction=np.array([0.0, 0.5]),
        phase_rad=np.array([0.0, np.pi]),
        instantaneous_e_over_n_Td=np.array([50.0, -50.0]),
        mean_energy_eV=np.array([case.mean_energy_eV, case.mean_energy_eV]),
        ionization_rate_coefficient_m3_s=np.array([1e-16, 2e-16]),
        eedf=np.stack([case.eedf, case.eedf]),
    )
    return case


def test_native_case_roundtrip_and_rf_csvs_are_lossless(tmp_path, phase_case):
    key = dict(
        mixture_id=0,
        solver=phase_case.solver,
        e_over_n_Td=phase_case.e_over_n_Td,
        replicate=0,
    )
    database = tmp_path / "cases.sqlite"
    with WorkflowStore(database) as store:
        store.write_mixture(0, [("Ar", 1.0, 39.948)])
        store.write_case(mixture_id=0, replicate=0, case=phase_case)
    with WorkflowStore(database) as store:
        restored = store.read_case(**key)
        assert store.has_case(**key)
    assert restored.physical_context == phase_case.physical_context
    assert restored.metadata == phase_case.metadata
    assert restored.transport == phase_case.transport
    assert restored.rates == phase_case.rates
    for name in asdict(phase_case.rf_phase):
        np.testing.assert_array_equal(
            getattr(restored.rf_phase, name), getattr(phase_case.rf_phase, name)
        )
    for label, case in (("direct", phase_case), ("restored", restored)):
        paths = write_outputs(SwarmRunResult([case]), OutputConfig(tmp_path / label))
        artifact = read_product_run(paths["manifest_json"])
        if label == "direct":
            expected = artifact.tables
        else:
            for name in (
                "summary_csv",
                "rates_csv",
                "eedf_csv",
                "rf_phase_csv",
                "rf_phase_eedf_csv",
            ):
                pd.testing.assert_frame_equal(artifact.tables[name], expected[name])


def test_import_preserves_rf_and_rejects_an_incomplete_source(tmp_path, phase_case):
    key = (0, phase_case.solver, phase_case.e_over_n_Td, 0)
    with (
        WorkflowStore(tmp_path / "source.sqlite") as source,
        WorkflowStore(tmp_path / "target.sqlite") as target,
    ):
        for store in (source, target):
            store.write_mixture(0, [("Ar", 1.0, 39.948)])
        source.write_case(mixture_id=0, replicate=0, case=phase_case)
        assert target.import_cases(source.connection, [key]) == 1
        restored = target.read_case(
            mixture_id=0, solver=key[1], e_over_n_Td=key[2], replicate=0
        )
        np.testing.assert_array_equal(restored.rf_phase.eedf, phase_case.rf_phase.eedf)
        source.connection.execute(
            "DELETE FROM rf_phase_eedf_bins WHERE phase_index=0 AND bin_index=0"
        )
        source.connection.commit()
        assert key not in source.existing_case_keys()
        with pytest.raises(WorkflowSchemaError, match="incomplete RF"):
            target.import_cases(source.connection, [key])
        with pytest.raises(WorkflowSchemaError, match="incomplete RF"):
            source.read_case(
                mixture_id=0, solver=key[1], e_over_n_Td=key[2], replicate=0
            )


def test_rf_child_failure_rolls_back_the_whole_case(tmp_path, phase_case):
    with WorkflowStore(tmp_path / "cases.sqlite") as store:
        store.write_mixture(0, [("Ar", 1.0, 39.948)])
        store.write_case(mixture_id=0, replicate=0, case=phase_case)
        store.connection.execute(
            "CREATE TRIGGER fail_phase BEFORE INSERT ON rf_phase_eedf_bins BEGIN SELECT RAISE(ABORT, 'test storage failure'); END"
        )
        store.connection.commit()
        changed = deepcopy(phase_case)
        changed.metadata["must_not_commit"] = True
        with pytest.raises(sqlite3.IntegrityError, match="storage failure"):
            store.write_case(mixture_id=0, replicate=0, case=changed)
        restored = store.read_case(
            mixture_id=0,
            solver=phase_case.solver,
            e_over_n_Td=phase_case.e_over_n_Td,
            replicate=0,
        )
        assert "must_not_commit" not in restored.metadata
        np.testing.assert_array_equal(restored.rf_phase.eedf, phase_case.rf_phase.eedf)
        changed.rf_phase.eedf[0, 0] = np.nan
        with pytest.raises(ValueError, match="finite"):
            store.write_case(mixture_id=0, replicate=0, case=changed)


def test_rf_campaign_resumes_and_repairs_missing_phase_data(tmp_path):
    data = base_product_config(tmp_path)
    data["physics"]["field"].update(
        type="time_dependent",
        time_dependent={
            "frequency_Hz": 1e6,
            "phase_steps": 8,
            "max_periods": 100,
            "periodic_tolerance": 1e-4,
        },
    )
    data["campaign"] = {"database": "cases.sqlite"}
    path = write_config(tmp_path, data)
    first = run_sweep(path)
    assert first.cases_written == 1
    assert run_sweep(path).cases_written == 0
    with WorkflowStore(first.database_path) as store:
        case = store.read_case(
            mixture_id=0, solver="two_term", e_over_n_Td=50, replicate=0
        )
        assert len(case.rf_phase.phase_fraction) == 8
        store.connection.execute(
            "DELETE FROM rf_phase_eedf_bins WHERE phase_index=0 AND bin_index=0"
        )
        store.connection.commit()
    assert run_sweep(path).cases_written == 1
    assert run_sweep(path).cases_written == 0


def test_old_storage_without_rf_tables_requires_recomputation(tmp_path):
    database = tmp_path / "old.sqlite"
    with WorkflowStore(database) as store:
        store.connection.execute("DROP TABLE rf_phase_eedf_bins")
        store.connection.commit()
    with pytest.raises(WorkflowSchemaError, match="RF phase persistence"):
        WorkflowStore(database)
