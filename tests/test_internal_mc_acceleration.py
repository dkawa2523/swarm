from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from electron_swarm import load_config, run
from electron_swarm.core.cross_sections import (
    CrossSectionProcess,
    ProcessType,
    load_active_mixture_inputs,
)
from electron_swarm.solvers.monte_carlo.compiled_kernel import (
    COMPILED_MC_KERNEL_SCHEMA_VERSION,
)
from electron_swarm.solvers.monte_carlo.case_result import (
    _tail_refinement_treatment,
)
from electron_swarm.solvers.monte_carlo.cross_section_table import (
    NUMBA_AVAILABLE,
    PreparedCrossSectionTable,
)
from electron_swarm.solvers.monte_carlo.population import _ParticleEnsemble
from electron_swarm.solvers.monte_carlo.setup import _prepare_monte_carlo_run


ROOT = Path(__file__).resolve().parents[1]


def test_tail_refinement_treatment_has_three_unambiguous_execution_states() -> None:
    disabled = SimpleNamespace(tail_collisions=0)
    configured = SimpleNamespace(tail_collisions=8)

    assert _tail_refinement_treatment(
        disabled, SimpleNamespace(tail_sampling_used=False)
    ) == "disabled"
    assert _tail_refinement_treatment(
        configured, SimpleNamespace(tail_sampling_used=False)
    ) == "configured_not_triggered"
    assert _tail_refinement_treatment(
        configured, SimpleNamespace(tail_sampling_used=True)
    ) == "executed"


def test_explicit_mc_tail_budget_is_independent_of_deterministic_grid_policy() -> None:
    config = load_config(
        ROOT / "examples" / "argon_gec_ccp_monte_carlo_weighted_branching.yaml"
    )
    config.physics.energy_grid_policy.threshold_refinement = False
    config.solvers.monte_carlo.tail_max_collisions = 8
    cross_sections = load_active_mixture_inputs(
        config.cross_sections,
        config.conditions,
    )

    setup = _prepare_monte_carlo_run(
        config,
        cross_sections,
        config.solvers.monte_carlo,
        collect_audit=False,
    )

    assert setup.tail_collisions == 8
    assert setup.tail_strata_edges is not None
    assert setup.tail_reaction_importance is not None


def test_omitted_mc_tail_budget_disables_tail_phase_and_records_zero() -> None:
    config = load_config(
        ROOT / "examples" / "argon_gec_ccp_monte_carlo_weighted_branching.yaml"
    )
    config.run.e_over_n_Td = [30.0]
    config.solvers.monte_carlo.seed = 11
    config.solvers.monte_carlo.particles = 4
    config.solvers.monte_carlo.warmup_collisions = 1
    config.solvers.monte_carlo.max_collisions = 4
    config.solvers.monte_carlo.tail_max_collisions = None
    config.solvers.monte_carlo.transport_correlation_lag_barriers = 4
    config.solvers.monte_carlo.numeric_kernel = "python"
    cross_sections = load_active_mixture_inputs(
        config.cross_sections,
        config.conditions,
    )

    setup = _prepare_monte_carlo_run(
        config,
        cross_sections,
        config.solvers.monte_carlo,
        collect_audit=False,
    )
    case = run(config, write=False).cases[0]
    provenance = case.diagnostics["internal_monte_carlo_transport"][
        "mc_run_provenance"
    ]

    assert setup.tail_collisions == 0
    assert setup.tail_strata_edges is None
    assert setup.tail_reaction_importance is None
    assert provenance["tail_max_collisions"] == 0
    assert provenance["tail_collisions_executed"] == 0
    assert provenance["tail_refinement_treatment"] == "disabled"
    assert case.metadata["tail_refinement_treatment"] == "disabled"
    assert case.diagnostics["internal_monte_carlo_transport"]["tail_sampling"][
        "state"
    ] == "disabled"
    assert case.diagnostics["internal_monte_carlo_transport"]["tail_sampling"][
        "reported_eedf_includes_main_production"
    ] is True
    assert case.diagnostics["internal_monte_carlo_transport"]["tail_sampling"][
        "reported_rates_include_main_production"
    ] is True
    assert case.diagnostics["internal_monte_carlo_reaction_rates"][
        "sampling_model"
    ] == "ordinary_trajectory_sampling"


def test_prepared_cross_sections_preserve_process_grids_and_endpoint_policies() -> None:
    hold = CrossSectionProcess(
        species="Ar",
        process="hold",
        process_type=ProcessType.ELASTIC,
        energy_eV=np.array([1.0, 3.0]),
        cross_section_m2=np.array([2.0, 6.0]),
        metadata={"high_energy_extrapolation": "hold"},
    )
    zero = CrossSectionProcess(
        species="Ar",
        process="zero",
        process_type=ProcessType.EXCITATION,
        energy_eV=np.array([0.0, 2.0, 4.0]),
        cross_section_m2=np.array([0.0, 8.0, 4.0]),
        metadata={"high_energy_extrapolation": "zero"},
    )
    multipliers = np.array([0.25, 0.75])
    table = PreparedCrossSectionTable.build(
        (hold, zero),
        multipliers=multipliers,
    )

    for energy in (-1.0, 0.0, 0.5, 1.0, 1.5, 3.0, 3.5, 4.0, 5.0):
        expected = np.array(
            [
                hold.sigma(np.array([energy]))[0] * multipliers[0],
                zero.sigma(np.array([energy]))[0] * multipliers[1],
            ]
        )
        assert table.evaluate(energy) == pytest.approx(expected, abs=0.0)

    error = CrossSectionProcess(
        species="Ar",
        process="error",
        process_type=ProcessType.IONIZATION,
        energy_eV=np.array([0.0, 2.0]),
        cross_section_m2=np.array([0.0, 1.0]),
        metadata={"high_energy_extrapolation": "error"},
    )
    error_table = PreparedCrossSectionTable.build((error, zero))
    with pytest.raises(ValueError, match="Ar:error"):
        error_table.evaluate(3.0)


def test_particle_ensemble_grows_geometrically_without_per_daughter_copy() -> None:
    rng = np.random.default_rng(7)
    ensemble = _ParticleEnsemble.initialize(2, 1.0, rng)
    positions, velocities, times, weights, lineages, _ = ensemble.compiled_storage()
    original_ids = tuple(map(id, (positions, velocities, times, weights, lineages)))

    for index in range(2):
        ensemble.append_particle(
            position=np.array([index + 1.0, 0.0, 0.0]),
            velocity=np.array([0.0, index + 2.0, 0.0]),
            time_s=0.5 + index,
            weight=0.25,
            lineage=10 + index,
        )

    storage = ensemble.compiled_storage()
    assert len(ensemble) == 4
    assert tuple(map(id, storage[:5])) == original_ids
    assert ensemble.lineages[-1] == 11

    ensemble.append_particle(
        position=np.zeros(3),
        velocity=np.ones(3),
        time_s=3.0,
        weight=0.5,
        lineage=12,
    )
    assert len(ensemble) == 5
    assert id(ensemble.compiled_storage()[0]) != original_ids[0]


@pytest.mark.mc
@pytest.mark.skipif(not NUMBA_AVAILABLE, reason="Numba is not installed")
def test_compiled_weighted_kernel_is_statistically_equivalent_to_reference() -> None:
    base = load_config(
        ROOT / "examples" / "argon_gec_ccp_monte_carlo_weighted_branching.yaml"
    )
    base.run.e_over_n_Td = [1000.0]
    base.solvers.monte_carlo.particles = 96
    base.solvers.monte_carlo.warmup_collisions = 24
    base.solvers.monte_carlo.max_collisions = 96
    base.solvers.monte_carlo.tail_max_collisions = 1
    base.solvers.monte_carlo.tail_rate_rse_trigger = 1.0
    base.solvers.monte_carlo.transport_correlation_lag_barriers = 8

    observations: dict[str, list[np.ndarray]] = {"python": [], "numba": []}
    cases = {}
    for seed in range(600, 612):
        for kernel in observations:
            config = deepcopy(base)
            config.solvers.monte_carlo.seed = seed
            config.solvers.monte_carlo.numeric_kernel = kernel
            case = run(config, write=False).cases[0]
            cases[kernel] = case
            observations[kernel].append(
                np.asarray(
                    [
                        case.mean_energy_eV,
                        case.drift_velocity_m_s,
                        *(item.rate_coefficient_m3_s for item in case.rates),
                    ],
                    dtype=float,
                )
            )

    reference = np.stack(observations["python"])
    compiled_values = np.stack(observations["numba"])
    reference_mean = np.mean(reference, axis=0)
    compiled_mean = np.mean(compiled_values, axis=0)
    combined_standard_error = np.sqrt(
        np.var(reference, axis=0, ddof=1) / reference.shape[0]
        + np.var(compiled_values, axis=0, ddof=1) / compiled_values.shape[0]
    )
    tolerance = np.maximum(
        0.02 * np.abs(reference_mean),
        2.0 * combined_standard_error,
    )
    assert np.all(np.abs(compiled_mean - reference_mean) <= tolerance)

    compiled = cases["numba"]
    provenance = compiled.diagnostics["internal_monte_carlo_transport"][
        "mc_run_provenance"
    ]
    assert provenance["numeric_kernel_requested"] == "numba"
    assert provenance["numeric_kernel_used"] == "numba"
    assert provenance["numeric_kernel_schema_version"] == (
        COMPILED_MC_KERNEL_SCHEMA_VERSION
    )
    assert provenance["tail_max_collisions"] == 1
    assert provenance["tail_collisions_executed"] == 0
    assert provenance["tail_refinement_treatment"] == "configured_not_triggered"
    assert compiled.metadata["tail_refinement_treatment"] == (
        "configured_not_triggered"
    )
    assert compiled.diagnostics["internal_monte_carlo_transport"]["tail_sampling"][
        "reported_eedf_includes_main_production"
    ] is True
    assert compiled.diagnostics["internal_monte_carlo_transport"]["tail_sampling"][
        "reported_rates_include_main_production"
    ] is True


@pytest.mark.mc
@pytest.mark.skipif(not NUMBA_AVAILABLE, reason="Numba is not installed")
def test_explicit_compiled_kernel_rejects_magnetic_orbit() -> None:
    config = load_config(
        ROOT / "examples" / "argon_gec_ccp_monte_carlo_weighted_branching.yaml"
    )
    config.run.e_over_n_Td = [30.0]
    config.solvers.monte_carlo.seed = 11
    config.solvers.monte_carlo.particles = 4
    config.solvers.monte_carlo.warmup_collisions = 0
    config.solvers.monte_carlo.max_collisions = 1
    config.solvers.monte_carlo.tail_max_collisions = 1
    config.physics.field.magnetic_field.enabled = True
    config.physics.field.magnetic_field.B_T = 0.01
    config.solvers.monte_carlo.numeric_kernel = "numba"

    with pytest.raises(
        NotImplementedError,
        match="magnetic_orbit",
    ):
        run(config, write=False)


@pytest.mark.mc
def test_detailed_event_audit_uses_reference_kernel_with_explicit_provenance() -> None:
    config = load_config(
        ROOT / "examples" / "argon_gec_ccp_monte_carlo_weighted_branching.yaml"
    )
    config.run.e_over_n_Td = [30.0]
    config.solvers.monte_carlo.seed = 9
    config.solvers.monte_carlo.particles = 4
    config.solvers.monte_carlo.warmup_collisions = 1
    config.solvers.monte_carlo.max_collisions = 4
    config.solvers.monte_carlo.tail_max_collisions = 2
    config.solvers.monte_carlo.transport_correlation_lag_barriers = 4
    config.solvers.monte_carlo.numeric_kernel = "numba"

    case = run(config, write=False, collect_diagnostics=True).cases[0]
    provenance = case.diagnostics["internal_monte_carlo_transport"][
        "mc_run_provenance"
    ]

    assert provenance["numeric_kernel_requested"] == "numba"
    assert provenance["numeric_kernel_used"] == "python"
    assert "detailed_event_audit" in provenance["numeric_kernel_fallback_reason"]
    assert provenance["tail_refinement_treatment"] == "executed"
    assert case.metadata["tail_refinement_treatment"] == "executed"
    assert case.diagnostics["internal_monte_carlo_transport"]["tail_sampling"][
        "state"
    ] == "executed"
    assert "internal_monte_carlo_audit" in case.diagnostics
