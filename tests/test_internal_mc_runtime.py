from __future__ import annotations

import pickle
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import electron_swarm.solvers.monte_carlo.batch as mc_batch
from electron_swarm import load_config, run
from electron_swarm.core.constants import E_CHARGE_C, ELECTRON_MASS_KG
from electron_swarm.solvers.monte_carlo.kinematics import boris_push
from electron_swarm.solvers.monte_carlo.seeding import derive_case_seed
from electron_swarm.solvers.monte_carlo.evidence import MC_CASE_SEED_DERIVATION


ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize(
    ("velocity", "electric_field", "dt_s"),
    [
        ([1.25e6, -2.5e5, 8.0e5], [0.0, 0.0, 3.5e4], 2.0e-12),
        ([-3.0, 2.0, -1.0], [4.0, -5.0, 6.0], 7.0e-10),
    ],
)
def test_zero_b_boris_fast_path_matches_reference_values(
    velocity: list[float],
    electric_field: list[float],
    dt_s: float,
) -> None:
    v = np.asarray(velocity, dtype=float)
    field = np.asarray(electric_field, dtype=float)
    zero_b = np.array([0.0, -0.0, 0.0])

    qmdt2 = -E_CHARGE_C / ELECTRON_MASS_KG * dt_s * 0.5
    v_minus = v + qmdt2 * field
    t = qmdt2 * zero_b
    s = 2.0 * t / (1.0 + float(np.dot(t, t)))
    v_prime = v_minus + np.cross(v_minus, t)
    v_plus = v_minus + np.cross(v_prime, s)
    reference = v_plus + qmdt2 * field

    actual = boris_push(v, field, zero_b, dt_s)

    np.testing.assert_array_equal(actual, reference)


def test_monte_carlo_batch_derives_one_rng_per_field_independent_of_order(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    prepared = object()
    config = SimpleNamespace(
        run=SimpleNamespace(e_over_n_Td=[75.0, 25.0], case_prefix="ordered")
    )
    solver_config = SimpleNamespace(seed=41)
    cross_sections = object()
    events: list[tuple[object, ...]] = []

    rng_by_seed: dict[int, object] = {}

    def fake_default_rng(seed: int) -> object:
        events.append(("rng", seed))
        return rng_by_seed.setdefault(seed, object())

    def fake_prepare(
        actual_config: object,
        actual_cross_sections: object,
        actual_solver_config: object,
        *,
        collect_audit: bool,
    ) -> object:
        events.append(
            (
                "prepare",
                actual_config,
                actual_cross_sections,
                actual_solver_config,
                collect_audit,
            )
        )
        return prepared

    def fake_case(
        actual_setup: object,
        actual_rng: object,
        *,
        e_over_n_Td: float,
        case_id: str,
        case_seed: int,
    ) -> str:
        events.append(
            (
                "case",
                actual_setup,
                actual_rng,
                e_over_n_Td,
                case_id,
                case_seed,
            )
        )
        return case_id

    monkeypatch.setattr(mc_batch.np.random, "default_rng", fake_default_rng)
    monkeypatch.setattr(mc_batch._setup, "_prepare_monte_carlo_run", fake_prepare)
    monkeypatch.setattr(mc_batch._case, "_run_monte_carlo_case", fake_case)

    results = mc_batch.run_internal_monte_carlo(
        config,  # type: ignore[arg-type]
        cross_sections,  # type: ignore[arg-type]
        solver_config,  # type: ignore[arg-type]
        collect_audit=True,
    )

    assert results == ["ordered_0000", "ordered_0001"]
    seed_75 = derive_case_seed(base_seed=41, e_over_n_Td=75.0)
    seed_25 = derive_case_seed(base_seed=41, e_over_n_Td=25.0)
    assert seed_75 != seed_25
    assert events == [
        ("prepare", config, cross_sections, solver_config, True),
        ("rng", seed_75),
        (
            "case",
            prepared,
            rng_by_seed[seed_75],
            75.0,
            "ordered_0000",
            seed_75,
        ),
        ("rng", seed_25),
        (
            "case",
            prepared,
            rng_by_seed[seed_25],
            25.0,
            "ordered_0001",
            seed_25,
        ),
    ]


@pytest.mark.mc
def test_fixed_base_seed_is_exact_and_each_field_matches_standalone() -> None:
    config = load_config(
        ROOT / "examples" / "argon_gec_ccp_monte_carlo_weighted_branching.yaml"
    )
    config.run.e_over_n_Td = [25.0, 75.0]
    mc = config.solvers.monte_carlo
    mc.population_model = "fixed_particle_single_daughter"
    mc.particles = 4
    mc.warmup_collisions = 1
    mc.max_collisions = 3
    mc.tail_max_collisions = None
    mc.seed = 123
    mc.numeric_kernel = "python"

    first = run(deepcopy(config), write=False).cases
    repeated = run(deepcopy(config), write=False).cases
    assert pickle.dumps(first, protocol=5) == pickle.dumps(repeated, protocol=5)

    standalone_config = deepcopy(config)
    standalone_config.run.e_over_n_Td = [75.0]
    standalone = run(standalone_config, write=False).cases[0]
    assert first[1].mean_energy_eV == standalone.mean_energy_eV
    np.testing.assert_array_equal(first[1].eedf_counts, standalone.eedf_counts)
    np.testing.assert_array_equal(first[1].eedf, standalone.eedf)
    assert first[1].metadata["monte_carlo_case_seed"] == standalone.metadata[
        "monte_carlo_case_seed"
    ]
    expected_case_seed = derive_case_seed(base_seed=123, e_over_n_Td=75.0)
    assert standalone.metadata["monte_carlo_base_seed"] == 123
    assert standalone.metadata["monte_carlo_case_seed"] == expected_case_seed
    provenance = standalone.diagnostics["internal_monte_carlo_transport"][
        "mc_run_provenance"
    ]
    assert provenance["seed"] == 123
    assert provenance["case_seed"] == expected_case_seed
    assert provenance["case_seed_derivation"] == MC_CASE_SEED_DERIVATION

    reordered_config = deepcopy(config)
    reordered_config.run.e_over_n_Td = [75.0, 25.0]
    reordered = run(reordered_config, write=False).cases
    by_field = {case.e_over_n_Td: case for case in first}
    reordered_by_field = {case.e_over_n_Td: case for case in reordered}
    for field in by_field:
        np.testing.assert_array_equal(
            by_field[field].eedf_counts,
            reordered_by_field[field].eedf_counts,
        )
        assert by_field[field].mean_energy_eV == reordered_by_field[field].mean_energy_eV


def test_direct_monte_carlo_requires_explicit_seed() -> None:
    config = load_config(
        ROOT / "examples" / "argon_gec_ccp_monte_carlo_weighted_branching.yaml"
    )
    config.run.e_over_n_Td = [25.0]
    config.solvers.monte_carlo.particles = 4
    config.solvers.monte_carlo.warmup_collisions = 1
    config.solvers.monte_carlo.max_collisions = 2
    config.solvers.monte_carlo.tail_max_collisions = None
    config.solvers.monte_carlo.seed = None

    with pytest.raises(ValueError, match=r"monte_carlo\.seed is required"):
        run(config, write=False)
