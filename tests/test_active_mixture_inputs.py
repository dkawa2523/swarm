from __future__ import annotations

import numpy as np
import pytest

from electron_swarm import load_config
from electron_swarm.core.config import ConditionsConfig, GasComponent
from electron_swarm.core.cross_sections import (
    ActiveMixtureInputs,
    CrossSectionProcess,
    CrossSectionSet,
    CrossSectionValidationError,
    ProcessType,
    prepare_active_mixture_inputs,
)
from electron_swarm.core.solver_configs import build_internal_solver_configs
from electron_swarm.core.solver_registry import build_solver
from electron_swarm.grids.energy import build_energy_grid
from electron_swarm.solvers.boltzmann_common.collisions import (
    build_effective_collision_data,
)
from electron_swarm.solvers.boltzmann_common.observables import (
    compute_rates_from_eedf,
)
from electron_swarm.solvers.monte_carlo.collisions import (
    PreparedCollisionSampler,
    _project_processes,
)
from swarm_workflow.tables.collision_kernels import (
    _combined_hash,
    _file_sha256,
    active_collision_kernel_rows,
)
from swarm_workflow.tables.contracts import TableBuildError
from tests.product_helpers import base_product_config, write_config


def _process(
    process_type: ProcessType,
    *,
    species: str = "Ar",
    threshold_eV: float | None = None,
) -> CrossSectionProcess:
    return CrossSectionProcess(
        species=species,
        process=f"{species}_{process_type.value}",
        process_type=process_type,
        threshold_eV=threshold_eV,
        mass_amu=39.948,
        energy_eV=np.asarray([0.0, 10.0]),
        cross_section_m2=np.asarray([1.0e-20, 1.0e-20]),
    )


@pytest.mark.parametrize(
    ("process_type", "message"),
    [
        (ProcessType.UNKNOWN, "Unknown cross-section process type"),
        (ProcessType.EXCITATION, "threshold_eV is required"),
        (ProcessType.IONIZATION, "threshold_eV is required"),
        (ProcessType.SUPERELASTIC, "threshold_eV is required"),
    ],
)
@pytest.mark.parametrize(
    "solver_id",
    ["two_term", "multi_term", "monte_carlo", "propagator"],
)
def test_all_solver_entrypoints_reject_invalid_active_processes(
    tmp_path,
    solver_id: str,
    process_type: ProcessType,
    message: str,
) -> None:
    config = load_config(
        write_config(
            tmp_path,
            base_product_config(tmp_path, [solver_id]),
        )
    )
    internal = build_internal_solver_configs(config.solvers, config.physics)

    with pytest.raises(CrossSectionValidationError, match=message):
        build_solver(
            config,
            CrossSectionSet([_process(process_type)]),
            internal,
            solver_id,
        )


def test_inactive_species_is_absent_from_numeric_inputs(tmp_path) -> None:
    conditions = ConditionsConfig(
        gas_mixture=[
            GasComponent("Ar", 1.0, 39.948),
            GasComponent("inactive", 0.0, 28.0),
        ]
    )
    active = _process(ProcessType.MOMENTUM)
    inactive = _process(
        ProcessType.EXCITATION,
        species="inactive",
        threshold_eV=5.25,
    )
    full_inputs = prepare_active_mixture_inputs(
        CrossSectionSet([active, inactive]),
        conditions,
    )
    active_only_inputs = prepare_active_mixture_inputs(
        CrossSectionSet([active]),
        conditions,
    )

    assert full_inputs.active_species == ("Ar",)
    assert full_inputs.processes == [active]
    grid_options = {
        "min_eV": 0.0,
        "max_eV": 10.0,
        "n": 16,
        "refine": True,
    }
    full_grid = build_energy_grid(
        cross_sections=full_inputs,
        **grid_options,
    )
    active_only_grid = build_energy_grid(
        cross_sections=active_only_inputs,
        **grid_options,
    )

    np.testing.assert_array_equal(
        full_grid.centers_eV,
        active_only_grid.centers_eV,
    )
    assert full_grid.metadata == active_only_grid.metadata

    data = base_product_config(tmp_path)
    data["conditions"]["gas_mixture"].append(
        {"species": "inactive", "fraction": 0.0, "mass_amu": 28.0}
    )
    config = load_config(write_config(tmp_path, data))
    internal = build_internal_solver_configs(config.solvers, config.physics)
    full_collisions = build_effective_collision_data(
        config,
        full_inputs,
        full_grid.centers_eV,
        1.0e22,
        internal.two_term,
    )
    active_collisions = build_effective_collision_data(
        config,
        active_only_inputs,
        active_only_grid.centers_eV,
        1.0e22,
        internal.two_term,
    )
    np.testing.assert_array_equal(full_collisions.nu_m, active_collisions.nu_m)
    assert {process.species for process in full_collisions.processes} == {"Ar"}

    eedf = np.ones_like(full_grid.centers_eV) / np.sum(full_grid.widths_eV)
    rates = compute_rates_from_eedf(
        config,
        full_inputs,
        full_grid.centers_eV,
        full_grid.widths_eV,
        eedf,
        case_id="active",
        e_over_n_Td=10.0,
        solver_name="two_term",
    )
    assert {rate.species for rate in rates.rates} == {"Ar"}

    projected = _project_processes(config, full_inputs)
    sampler = PreparedCollisionSampler.build(
        config,
        projected,
        angular_model_name="isotropic",
        max_energy_eV=10.0,
    )
    assert {process.species for process in sampler.table.processes} == {"Ar"}


def test_table_collision_kernel_cannot_materialize_an_inactive_rate(
    tmp_path,
) -> None:
    inactive_path = tmp_path / "inactive.csv"
    inactive_path.write_text(
        "species,process,type,threshold_eV,mass_amu,energy_eV,cross_section_m2\n"
        "inactive,inactive momentum,momentum,0.0,28.0,0.0,9e-20\n"
        "inactive,inactive momentum,momentum,0.0,28.0,10.0,9e-20\n",
        encoding="utf-8",
    )
    data = base_product_config(tmp_path)
    data["conditions"]["gas_mixture"].append(
        {"species": "inactive", "fraction": 0.0, "mass_amu": 28.0}
    )
    data["cross_sections"]["files"].append(
        {
            "path": inactive_path.as_posix(),
            "species": "inactive",
            "format": "csv",
        }
    )
    config_path = write_config(tmp_path, data)
    config = load_config(config_path)
    hashes = {
        str(item.path): _file_sha256(item.path)
        for item in config.cross_sections.files
    }
    metadata = {
        "base_config_path": str(config_path),
        "cross_sections_sha256": _combined_hash(hashes),
    }
    active_rate = {
        "species": "Ar",
        "process": "Ar momentum transfer",
        "process_type": "momentum",
        "threshold_eV": 0.0,
    }
    mixture = [{"species": "Ar", "fraction": 1.0, "mass_amu": 39.948}]
    _columns, rows = active_collision_kernel_rows(
        metadata,
        [active_rate],
        mixture=mixture,
    )
    assert {row["species"] for row in rows} == {"Ar"}

    inactive_rate = {
        "species": "inactive",
        "process": "inactive momentum",
        "process_type": "momentum",
        "threshold_eV": 0.0,
    }
    with pytest.raises(TableBuildError, match="does not map"):
        active_collision_kernel_rows(
            metadata,
            [active_rate, inactive_rate],
            mixture=mixture,
        )


@pytest.mark.parametrize(
    "process",
    [
        pytest.param(_process(ProcessType.UNKNOWN), id="unknown"),
        pytest.param(_process(ProcessType.EXCITATION), id="missing-threshold"),
    ],
)
def test_prebuilt_active_inventory_cannot_bypass_process_validation(
    process: CrossSectionProcess,
) -> None:
    conditions = ConditionsConfig(
        gas_mixture=[GasComponent("Ar", 1.0, 39.948)]
    )
    forged = ActiveMixtureInputs(
        processes=[process],
        components=tuple(conditions.gas_mixture),
    )

    with pytest.raises(CrossSectionValidationError):
        prepare_active_mixture_inputs(forged, conditions)


def test_revalidated_active_inventory_does_not_share_its_process_list() -> None:
    conditions = ConditionsConfig(
        gas_mixture=[GasComponent("Ar", 1.0, 39.948)]
    )
    process = _process(ProcessType.MOMENTUM)
    caller_owned = ActiveMixtureInputs(
        processes=[process],
        components=tuple(conditions.gas_mixture),
    )

    solver_owned = prepare_active_mixture_inputs(caller_owned, conditions)
    caller_owned.processes.clear()

    assert solver_owned.processes == [process]
