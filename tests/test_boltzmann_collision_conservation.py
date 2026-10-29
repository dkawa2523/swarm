from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from electron_swarm import load_config, run
from electron_swarm.core.cross_sections import (
    CrossSectionProcess,
    CrossSectionSet,
    ProcessType,
    ScatteringRole,
)
from electron_swarm.core.capabilities import SupportLevel
from electron_swarm.core.solver_configs import build_internal_solver_configs
from electron_swarm.physics.electron_neutral import (
    first_order_elastic_recoil_energy_eV,
)
from electron_swarm.core.gas import gas_number_density
from electron_swarm.orchestration.plan import build_solve_plan
from electron_swarm.solvers.boltzmann_common.collisions import (
    build_effective_collision_data,
    build_inelastic_collision_data,
)
from electron_swarm.solvers.boltzmann_common.grid import electron_speed_m_s
from electron_swarm.solvers.multi_term.case import build_multiterm_case
from electron_swarm.solvers.multi_term.operator import assemble_pn_operator

from product_helpers import base_product_config, write_config, write_moment_table


def _elastic_process(role: ScatteringRole, sigma: float) -> CrossSectionProcess:
    kind = {
        ScatteringRole.ELASTIC_TOTAL: ProcessType.ELASTIC,
        ScatteringRole.ELASTIC_MOMENTUM_TRANSFER: ProcessType.MOMENTUM,
        ScatteringRole.EFFECTIVE_MOMENTUM_TRANSFER: ProcessType.EFFECTIVE,
    }[role]
    return CrossSectionProcess(
        species="Ar",
        process=role.value,
        process_type=kind,
        scattering_role=role,
        mass_amu=39.948,
        energy_eV=np.array([0.0, 1000.0]),
        cross_section_m2=np.full(2, sigma),
    )


@pytest.mark.parametrize(
    ("sharing", "nonconservative", "number_change"),
    [
        ("equal", "growth", 1.0),
        ("primary_secondary", "growth", 1.0),
        ("loss_only", "growth", 0.0),
        ("primary_secondary", "ignore", 0.0),
    ],
)
@pytest.mark.parametrize("incident", [16.0, 20.0])
def test_ionization_transfer_preserves_threshold_energy_loss_and_particle_number(
    tmp_path: Path,
    sharing: str,
    nonconservative: str,
    number_change: float,
    incident: float,
) -> None:
    data = base_product_config(tmp_path)
    data["physics"]["ionization"] = {
        "energy_sharing": sharing,
        "secondary_electron_energy_eV": 2.0 if sharing == "primary_secondary" else 0.0,
    }
    config = load_config(write_config(tmp_path, data))
    energy = np.array([0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 16.0, 20.0])
    widths = np.arange(1.0, 9.0)
    process = CrossSectionProcess(
        species="Ar", process="ionization", process_type=ProcessType.IONIZATION,
        threshold_eV=15.0, energy_eV=np.array([0.0, 100.0]),
        cross_section_m2=np.full(2, 1.0e-20),
    )
    density = gas_number_density(config)
    internal = build_internal_solver_configs(config.solvers, config.physics)
    internal.two_term.nonconservative_model = nonconservative
    collisions = build_inelastic_collision_data(
        config, CrossSectionSet([process]), energy, widths, density, internal.two_term,
    )
    source = int(np.flatnonzero(energy == incident)[0])
    population = np.zeros_like(energy)
    population[source] = 1.0 / widths[source]
    frequency = density * electron_speed_m_s(incident) * 1.0e-20
    change_per_event = collisions.matrix @ population / frequency

    assert np.sum(change_per_event * widths) == pytest.approx(number_change)
    assert collisions.net_creation_frequency_s_inv[source] / frequency == pytest.approx(
        number_change
    )
    assert -np.sum(energy * change_per_event * widths) == pytest.approx(15.0)
    daughters = change_per_event + population
    assert np.all(daughters >= -1.0e-14)
    assert np.sum(energy * daughters * widths) == pytest.approx(incident - 15.0)


@pytest.mark.parametrize("sharing", ["equal", "primary_secondary", "loss_only"])
@pytest.mark.parametrize("nonconservative", ["growth", "ignore"])
def test_discrete_number_source_matches_ionization_attachment_and_conservative_transfers(
    tmp_path: Path, sharing: str, nonconservative: str,
) -> None:
    data = base_product_config(tmp_path)
    data["physics"]["ionization"] = {
        "energy_sharing": sharing,
        "secondary_electron_energy_eV": 2.0 if sharing == "primary_secondary" else 0.0,
    }
    config = load_config(write_config(tmp_path, data))
    energy = np.array([0.0, 1.0, 2.0, 4.0, 8.0, 12.0])
    widths = np.arange(1.0, 7.0)
    processes = [
        CrossSectionProcess(
            species="Ar", process=kind.value, process_type=kind,
            threshold_eV=threshold, energy_eV=np.array([0.0, 100.0]),
            cross_section_m2=np.full(2, sigma),
        )
        for kind, threshold, sigma in [
            (ProcessType.EXCITATION, 1.0, 3.0e-20),
            (ProcessType.IONIZATION, 4.0, 2.0e-20),
            (ProcessType.ATTACHMENT, 0.0, 1.0e-20),
        ]
    ]
    density = gas_number_density(config)
    internal = build_internal_solver_configs(config.solvers, config.physics)
    internal.two_term.nonconservative_model = nonconservative
    collisions = build_inelastic_collision_data(
        config, CrossSectionSet(processes), energy, widths, density, internal.two_term,
    )
    expected = np.zeros_like(energy)
    if nonconservative == "growth":
        expected -= processes[2].sigma(energy) * electron_speed_m_s(energy) * density
        if sharing != "loss_only":
            expected += processes[1].sigma(energy) * electron_speed_m_s(energy) * density
    assert collisions.net_creation_frequency_s_inv == pytest.approx(expected)
    distribution = np.exp(-energy / 3.0)
    distribution /= np.sum(distribution * widths)
    measured_change = np.sum((collisions.matrix @ distribution) * widths)
    expected_change = np.sum(expected * distribution * widths)
    assert measured_change == pytest.approx(expected_change, rel=1.0e-12, abs=1.0e-6)


@pytest.mark.parametrize("angular", ["isotropic", "maxent_p1", "momentum_power"])
def test_elastic_recoil_uses_explicit_momentum_without_effective_or_inelastic_loss(
    tmp_path: Path, angular: str,
) -> None:
    data = base_product_config(tmp_path)
    data["physics"]["angular_scattering"] = {
        "model": angular,
        "higher_moment_closure": {
            "isotropic": "zero", "maxent_p1": "maxent", "momentum_power": "power",
        }[angular],
    }
    config = load_config(write_config(tmp_path, data))
    energy = np.array([1.0, 5.0, 10.0])
    density = gas_number_density(config)
    internal = build_internal_solver_configs(config.solvers, config.physics)
    momentum = _elastic_process(ScatteringRole.ELASTIC_MOMENTUM_TRANSFER, 1.0e-20)
    reference = build_effective_collision_data(
        config, CrossSectionSet([momentum]), energy, density, internal.two_term,
    )
    excitation = CrossSectionProcess(
        species="Ar", process="excitation", process_type=ProcessType.EXCITATION,
        threshold_eV=0.5, energy_eV=np.array([0.0, 1000.0]),
        cross_section_m2=np.full(2, 4.0e-20),
    )
    augmented = build_effective_collision_data(
        config,
        CrossSectionSet([
            momentum,
            _elastic_process(ScatteringRole.ELASTIC_TOTAL, 10.0e-20),
            _elastic_process(ScatteringRole.EFFECTIVE_MOMENTUM_TRANSFER, 20.0e-20),
            excitation,
        ]),
        energy, density, internal.two_term,
    )
    assert augmented.elastic_A_eV_s == pytest.approx(reference.elastic_A_eV_s)
    assert augmented.elastic_D_eV2_s == pytest.approx(reference.elastic_D_eV2_s)
    assert augmented.nu_m == pytest.approx(reference.nu_m * 20.0)
    assert np.all(augmented.inelastic_loss_frequency_s_inv > 0.0)

    # The cold component equals the angular average of the independent
    # first-mass-ratio collision kinematics, not the total collision count.
    cold_loss = (
        2.0 * energy - first_order_elastic_recoil_energy_eV(energy, -1.0, 39.948)
        - energy
    ) / 2.0
    cold_rate = density * electron_speed_m_s(energy) * 1.0e-20 * cold_loss
    temperature_free_drift = (
        -reference.elastic_A_eV_s
        + reference.elastic_D_eV2_s / (2.0 * energy)
    )
    assert temperature_free_drift == pytest.approx(cold_rate, rel=1.0e-10)


def test_effective_only_elastic_energy_relaxation_fails_explicitly(tmp_path: Path) -> None:
    config = load_config(write_config(tmp_path, base_product_config(tmp_path)))
    internal = build_internal_solver_configs(config.solvers, config.physics)
    with pytest.raises(NotImplementedError, match="cannot define elastic recoil"):
        build_effective_collision_data(
            config,
            CrossSectionSet([
                _elastic_process(ScatteringRole.EFFECTIVE_MOMENTUM_TRANSFER, 1.0e-20),
            ]),
            np.array([1.0, 2.0]), gas_number_density(config), internal.two_term,
        )


@pytest.mark.parametrize("solver", ["two_term", "multi_term"])
@pytest.mark.parametrize("degraded", ["record", "fail"])
@pytest.mark.parametrize("unsupported", ["fail", "skip_solver"])
def test_missing_elastic_recoil_inputs_follow_feature_policy_before_execution(
    tmp_path: Path, solver: str, degraded: str, unsupported: str,
) -> None:
    data = base_product_config(tmp_path, [solver])
    data["feature_policy"] = {"unsupported": unsupported, "degraded": degraded}
    config = load_config(write_config(tmp_path, data))
    cross_sections = CrossSectionSet([
        _elastic_process(ScatteringRole.EFFECTIVE_MOMENTUM_TRANSFER, 1.0e-20),
    ])
    if unsupported == "fail":
        with pytest.raises(ValueError, match="cannot define elastic recoil"):
            build_solve_plan(config, cross_sections)
    else:
        [item] = build_solve_plan(config, cross_sections)
        assert item.skipped is True
        assert item.runnable is False
        angular = item.feature_treatments["angular_scattering"]
        assert angular.support == SupportLevel.UNSUPPORTED
        assert angular.treatment == "isotropic:missing_elastic_recoil_inputs"


@pytest.mark.parametrize("unsupported", ["fail", "skip_solver"])
def test_isotropic_pn_with_unequal_total_and_momentum_requires_explicit_closure(
    tmp_path: Path, unsupported: str,
) -> None:
    data = base_product_config(tmp_path, ["multi_term"])
    data["feature_policy"]["unsupported"] = unsupported
    config = load_config(write_config(tmp_path, data))
    cross_sections = CrossSectionSet([
        _elastic_process(ScatteringRole.ELASTIC_TOTAL, 10.0e-20),
        _elastic_process(ScatteringRole.ELASTIC_MOMENTUM_TRANSFER, 1.0e-20),
    ])
    if unsupported == "fail":
        with pytest.raises(ValueError, match="use maxent_p1"):
            build_solve_plan(config, cross_sections)
    else:
        [item] = build_solve_plan(config, cross_sections)
        assert item.skipped is True
        assert "use maxent_p1" in item.skip_reason
    internal = build_internal_solver_configs(config.solvers, config.physics)
    case = build_multiterm_case(config, cross_sections, internal.multi_term, 1.0)
    with pytest.raises(ValueError, match="use maxent_p1"):
        assemble_pn_operator(case)
    config.physics.angular_scattering.model = "maxent_p1"
    config.physics.angular_scattering.higher_moment_closure = "maxent"
    [item] = build_solve_plan(config, cross_sections)
    assert item.runnable is True


@pytest.mark.parametrize("total_scale", [1.0, 10.0])
def test_two_term_first_moment_assumption_exposes_isotropic_particle_mismatch(
    tmp_path: Path, total_scale: float,
) -> None:
    config = load_config(write_config(
        tmp_path, base_product_config(tmp_path, ["two_term", "monte_carlo"]),
    ))
    cross_sections = CrossSectionSet([
        _elastic_process(ScatteringRole.ELASTIC_TOTAL, total_scale * 1.0e-20),
        _elastic_process(ScatteringRole.ELASTIC_MOMENTUM_TRANSFER, 1.0e-20),
    ])
    two_term, monte_carlo = build_solve_plan(config, cross_sections)
    assert two_term.runnable and monte_carlo.runnable
    two_term_assumption = two_term.assumption("angular_scattering")
    particle_assumption = monte_carlo.assumption("angular_scattering")
    if total_scale == 1.0:
        assert two_term_assumption == particle_assumption == "isotropic:zero"
    else:
        assert two_term_assumption == "isotropic:zero;first_moment=1-sigma_m/sigma_0:Ar"
        assert two_term_assumption != particle_assumption


@pytest.mark.parametrize(
    ("total_tail", "momentum_values", "ceiling", "match"),
    [
        ("zero", [2.0e-20, 2.0e-20, 2.0e-20], 30.0, False),
        ("zero", [2.0e-20, 2.0e-20, 0.0], 30.0, False),
        ("hold", [2.0e-20, 2.0e-20, 2.0e-20], 30.0, True),
        ("zero", [2.0e-20, 2.0e-20, 2.0e-20], 5.0, True),
        ("error", [2.0e-20, 2.0e-20, 2.0e-20], 5.0, True),
    ],
)
def test_isotropic_assumptions_respect_unequal_table_support_and_extrapolation(
    tmp_path: Path,
    total_tail: str,
    momentum_values: list[float],
    ceiling: float,
    match: bool,
) -> None:
    data = base_product_config(tmp_path, ["two_term", "monte_carlo"])
    data["physics"]["energy_grid_policy"]["max_eV_limit"] = ceiling
    config = load_config(write_config(tmp_path, data))
    total = CrossSectionProcess(
        species="Ar", process="total", process_type=ProcessType.ELASTIC,
        energy_eV=np.array([0.0, 10.0]), cross_section_m2=np.full(2, 2.0e-20),
        metadata={"high_energy_extrapolation": total_tail},
    )
    momentum = CrossSectionProcess(
        species="Ar", process="momentum", process_type=ProcessType.MOMENTUM,
        energy_eV=np.array([0.0, 10.0, 20.0]),
        cross_section_m2=np.array(momentum_values),
        metadata={"high_energy_extrapolation": "zero" if total_tail == "zero" else "hold"},
    )
    two_term, monte_carlo = build_solve_plan(config, CrossSectionSet([total, momentum]))
    same_assumption = (
        two_term.assumption("angular_scattering")
        == monte_carlo.assumption("angular_scattering")
    )
    assert same_assumption is match
    data["run"]["solvers"] = [{"id": "multi_term"}]
    config = load_config(write_config(tmp_path, data))
    if match:
        [pn] = build_solve_plan(config, CrossSectionSet([total, momentum]))
        assert pn.runnable is True
    else:
        with pytest.raises(ValueError, match="use maxent_p1"):
            build_solve_plan(config, CrossSectionSet([total, momentum]))


@pytest.mark.parametrize("source", ["total", "momentum", "both"])
def test_dcs_first_moment_owns_elastic_recoil_and_pn_particle_conservation(
    tmp_path: Path, source: str,
) -> None:
    data = base_product_config(tmp_path, ["multi_term"])
    data["solvers"]["multi_term"]["method"] = "pn_dcs"
    path = write_moment_table(tmp_path, m1=0.9)
    data["physics"]["angular_scattering"] = {
        "model": "moment_table", "higher_moment_closure": "table",
        "moment_table": {"path": path.as_posix(), "provenance": "dcs_derived"},
    }
    config = load_config(write_config(tmp_path, data))
    internal = build_internal_solver_configs(config.solvers, config.physics)
    momentum = _elastic_process(ScatteringRole.ELASTIC_MOMENTUM_TRANSFER, 1.0e-20)
    total = _elastic_process(ScatteringRole.ELASTIC_TOTAL, 10.0e-20)
    processes = {"total": [total], "momentum": [momentum], "both": [total, momentum]}[source]
    case = build_multiterm_case(config, CrossSectionSet(processes), internal.multi_term, 1.0)
    operator = assemble_pn_operator(case)
    reference_case = build_multiterm_case(
        config, CrossSectionSet([momentum]), internal.multi_term, 1.0,
    )
    reference = assemble_pn_operator(reference_case)
    difference = operator.matrix - reference.matrix
    assert np.max(np.abs(difference.data), initial=0.0) <= (
        1.0e-12 * np.max(np.abs(reference.matrix.data))
    )
    particle_weights = np.zeros(operator.matrix.shape[0])
    particle_weights[operator.moment_slices[0]] = operator.widths_eV
    column_balance = particle_weights @ operator.matrix
    scale = np.asarray(abs(operator.matrix).T @ particle_weights).reshape(-1)
    assert np.max(np.abs(column_balance) / np.maximum(scale, 1.0)) < 1.0e-12


def test_inconsistent_dcs_total_momentum_and_first_moment_are_rejected(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["multi_term"])
    path = write_moment_table(tmp_path, m1=0.9)
    data["solvers"]["multi_term"]["method"] = "pn_dcs"
    data["physics"]["angular_scattering"] = {
        "model": "moment_table", "higher_moment_closure": "table",
        "moment_table": {"path": path.as_posix(), "provenance": "dcs_derived"},
    }
    config = load_config(write_config(tmp_path, data))
    internal = build_internal_solver_configs(config.solvers, config.physics)
    with pytest.raises(ValueError, match="disagree.*supplied first angular moment"):
        build_effective_collision_data(
            config,
            CrossSectionSet([
                _elastic_process(ScatteringRole.ELASTIC_TOTAL, 10.0e-20),
                _elastic_process(ScatteringRole.ELASTIC_MOMENTUM_TRANSFER, 2.0e-20),
            ]),
            np.array([1.0, 2.0]), gas_number_density(config), internal.multi_term,
        )


def test_two_term_maxent_cooling_is_independent_of_total_at_fixed_momentum(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path)
    data["run"]["e_over_n_Td"] = [1.0]
    data["physics"]["energy_grid_policy"]["adaptive"] = False
    data["physics"]["angular_scattering"] = {
        "model": "maxent_p1", "higher_moment_closure": "maxent",
    }
    means = []
    for scale in [1.0, 10.0]:
        path = tmp_path / f"elastic_{int(scale)}.csv"
        path.write_text(
            "species,process,type,threshold_eV,mass_amu,energy_eV,cross_section_m2\n"
            f"Ar,total,elastic_total,0,39.948,0,{scale * 1e-20}\n"
            f"Ar,total,elastic_total,0,39.948,1000,{scale * 1e-20}\n"
            "Ar,momentum,elastic_momentum_transfer,0,39.948,0,1e-20\n"
            "Ar,momentum,elastic_momentum_transfer,0,39.948,1000,1e-20\n",
            encoding="utf-8",
        )
        data["cross_sections"]["files"] = [{"path": path.as_posix(), "format": "csv"}]
        config = load_config(write_config(tmp_path, data))
        [case] = run(config, write=False).cases
        assert case.diagnostics["two_term"]["converged"] is True
        assert case.diagnostics["two_term"]["elastic_energy_loss"]["status"] == "available"
        means.append(case.mean_energy_eV)
    assert means[0] == pytest.approx(means[1], rel=1.0e-12)


def test_anisotropic_mixture_generator_is_the_sum_of_species_generators(tmp_path: Path) -> None:
    from copy import deepcopy
    from electron_swarm.core.config import GasComponent

    data = base_product_config(tmp_path, ["multi_term"])
    data["physics"]["angular_scattering"] = {"model": "maxent_p1", "higher_moment_closure": "maxent"}
    config = load_config(write_config(tmp_path, data))
    config.conditions.gas_mixture = [GasComponent("Ar", 0.5, 39.948), GasComponent("He", 0.5, 4.0)]
    processes = []
    for species, mass, first in [("Ar", 39.948, 0.0), ("He", 4.0, 0.9)]:
        for role, value in [(ScatteringRole.ELASTIC_TOTAL, 4e-20),
                            (ScatteringRole.ELASTIC_MOMENTUM_TRANSFER, 4e-20 * (1 - first))]:
            process = _elastic_process(role, value)
            process.species, process.mass_amu = species, mass
            processes.append(process)
    xs = CrossSectionSet(processes)
    internal = build_internal_solver_configs(config.solvers, config.physics)

    def matrix(fraction: float):
        variant = deepcopy(config)
        variant.conditions.gas_mixture[0].fraction = fraction
        variant.conditions.gas_mixture[1].fraction = 1 - fraction
        return assemble_pn_operator(build_multiterm_case(variant, xs, internal.multi_term, 50)).matrix

    mixture, argon, helium = matrix(0.5), matrix(1.0), matrix(0.0)
    delta = mixture - 0.5 * (argon + helium)
    assert np.max(np.abs(delta.data), initial=0.0) <= 1e-12 * np.max(np.abs(mixture.data))


def test_momentum_regularization_is_shared_by_collision_frequency_and_transport(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["two_term"])
    config = load_config(write_config(tmp_path, data))
    internal = build_internal_solver_configs(config.solvers, config.physics).two_term
    internal.min_momentum_cross_section_m2 = 1e-18
    energy = np.array([0.1, 1.0, 10.0])
    density = gas_number_density(config)
    collision = build_effective_collision_data(
        config, CrossSectionSet([_elastic_process(ScatteringRole.ELASTIC_MOMENTUM_TRANSFER, 1e-20)]),
        energy, density, internal,
    )
    assert collision.sigma_m == pytest.approx(np.full(3, 1e-18))
    assert collision.nu_m == pytest.approx(density * electron_speed_m_s(energy) * collision.sigma_m)
