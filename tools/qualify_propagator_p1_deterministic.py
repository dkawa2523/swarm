"""Run the Propagator-only deterministic P0/P1 qualification matrix.

This command deliberately does not instantiate two-term, multi-term, Monte
Carlo, or COMSOL paths.  It qualifies the bounded homogeneous-DC P1 numerical
core by operator identities, grid refinement, energy-ceiling independence,
and a mixed-gas angular-refinement case.
"""

from __future__ import annotations

import argparse
import copy
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
from time import perf_counter
from typing import Any

import numpy as np
import scipy


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from electron_swarm import load_config  # noqa: E402
from electron_swarm.core.config import RequestedSolverConfig, SwarmConfig  # noqa: E402
from electron_swarm.core.cross_sections import (  # noqa: E402
    load_active_mixture_inputs,
)
from electron_swarm.core.solver_configs import (  # noqa: E402
    build_internal_solver_configs,
)
from electron_swarm.solvers.propagator.collisions import (  # noqa: E402
    build_collision_operator,
)
from electron_swarm.solvers.propagator.grid import (  # noqa: E402
    build_propagator_grid,
)
from electron_swarm.solvers.propagator.shell_response import (  # noqa: E402
    SHELL_COEFFICIENT_ERROR_TOLERANCE,
    SHELL_MAX_COEFFICIENT_SEGMENTS,
    SHELL_RATIONAL_ERROR_TOLERANCE,
)
from electron_swarm.solvers.propagator.solver import PropagatorSolver  # noqa: E402
from electron_swarm.solvers.propagator.steady import (  # noqa: E402
    PropagatorConvergenceError,
)
from swarm_workflow.quality.propagator_source import (  # noqa: E402
    propagator_qualification_source_fingerprint,
)


LOW_FIELD_TD = (0.05, 0.1, 1.0)
INELASTIC_FIELD_TD = (10.0, 100.0, 500.0)
GRID_PROFILES = {"medium": (300, 36), "fine": (600, 72)}
MIXTURE_PROFILES = {"medium": (64, 36), "fine": (64, 72)}
ENERGY_CEILINGS_EV = (7_500.0, 15_000.0, 30_000.0)
SCALAR_REFINEMENT_LIMIT = 0.01
EEDF_L1_REFINEMENT_LIMIT = 0.02
BOUNDARY_SCALAR_LIMIT = 1.0e-5
BOUNDARY_EEDF_L1_LIMIT = 1.0e-5
QUALIFICATION_INPUT_PATHS = (
    ROOT / "examples" / "argon_propagator.yaml",
    ROOT / "configs" / "benchmarks" / "argon_propagator_maxent_p1.yaml",
    ROOT / "configs" / "benchmarks" / "ar_o2_propagator_generalization.yaml",
    ROOT / "examples" / "cross_sections" / "argon_application_library.csv",
    ROOT
    / "examples"
    / "cross_sections"
    / "argon_magboltz_11_17_gas2_elastic.csv",
    ROOT
    / "examples"
    / "cross_sections"
    / "argon_magboltz_11_17_gas2_elastic.provenance.json",
    ROOT / "examples" / "cross_sections" / "ar_o2_minimal.csv",
)
FULL_QUALIFICATION_OUTPUT = (
    ROOT
    / "docs"
    / "dev"
    / "results"
    / "propagator_p1_deterministic_qualification_20260908.json"
)
QUICK_QUALIFICATION_OUTPUT = (
    ROOT / "outputs" / "benchmarks" / "propagator_p1_deterministic_quick.json"
)


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_value(*args: str) -> str | None:
    result = subprocess.run(
        ["git", *args],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    return result.stdout.strip() or None


def _qualification_source_state(
    input_paths: tuple[Path, ...] = QUALIFICATION_INPUT_PATHS,
) -> dict[str, Any]:
    return {
        "implementation_fingerprint": (
            propagator_qualification_source_fingerprint(ROOT)
        ),
        "input_sha256": {
            path.relative_to(ROOT).as_posix(): _sha256(path)
            for path in input_paths
        },
    }


def _require_unchanged_qualification_sources(
    starting_state: dict[str, Any],
    input_paths: tuple[Path, ...] = QUALIFICATION_INPUT_PATHS,
) -> None:
    try:
        current_state = _qualification_source_state(input_paths)
    except OSError as exc:
        raise RuntimeError(
            "Propagator qualification sources or inputs changed during execution; "
            "discard the results and rerun"
        ) from exc
    if current_state != starting_state:
        raise RuntimeError(
            "Propagator qualification sources or inputs changed during execution; "
            "discard the results and rerun"
        )


def _variant(
    base: SwarmConfig,
    field_Td: float,
    grid: tuple[int, int],
    *,
    maximum_energy_eV: float | None = None,
) -> SwarmConfig:
    config = copy.deepcopy(base)
    config.run.solvers = [RequestedSolverConfig(id="propagator")]
    config.run.e_over_n_Td = [float(field_Td)]
    config.run.case_prefix = "propagator_deterministic_qualification"
    config.solvers.propagator.energy_cells = int(grid[0])
    config.solvers.propagator.polar_cells = int(grid[1])
    if maximum_energy_eV is not None:
        config.physics.energy_grid_policy.max_eV_limit = float(
            maximum_energy_eV
        )
    return config


def _run_case(
    config: SwarmConfig,
    *,
    family: str,
    profile: str,
) -> tuple[dict[str, Any], Any | None]:
    started = perf_counter()
    try:
        active_inputs = load_active_mixture_inputs(
            config.cross_sections,
            config.conditions,
        )
        solver_config = build_internal_solver_configs(
            config.solvers,
            config.physics,
        ).propagator
        [case] = PropagatorSolver(
            config,
            active_inputs,
            solver_config,
        ).solve_all()
        diagnostic = case.diagnostics["propagator"]
        timing = diagnostic["timings_s"]
        quality = bool(
            diagnostic["operator_residual_L1"] <= 1.0e-8
            and diagnostic["number_balance_residual"] <= 1.0e-12
            and diagnostic["negative_population_mass"] <= 1.0e-14
            and diagnostic["tail_probability"]
            <= diagnostic["tail_probability_target"]
            and diagnostic["outer_acceleration_flux_fraction"]
            <= max(diagnostic["tail_probability_target"], 1.0e-10)
            and timing["response_shell_rational_error_estimate"]
            <= SHELL_RATIONAL_ERROR_TOLERANCE
            and timing["response_shell_coefficient_error_estimate"]
            <= SHELL_COEFFICIENT_ERROR_TOLERANCE
            and timing["response_shell_maximum_coefficient_segments"]
            <= SHELL_MAX_COEFFICIENT_SEGMENTS
            and timing["response_cone_compatible_ritz_pairs"] >= 1.0
            and diagnostic["estimated_peak_memory_bytes"]
            < config.solvers.propagator.max_memory_mb * 1024 * 1024
        )
        return (
            {
                "family": family,
                "profile": profile,
                "E_over_N_Td": float(case.e_over_n_Td),
                "status": "ok",
                "elapsed_s": perf_counter() - started,
                "energy_cells": int(diagnostic["energy_cells"]),
                "polar_cells": int(diagnostic["polar_cells"]),
                "energy_max_eV": float(diagnostic["energy_max_eV"]),
                "mean_energy_eV": float(case.mean_energy_eV),
                "drift_velocity_m_s": float(case.drift_velocity_m_s),
                "growth_frequency_s_inv": float(
                    diagnostic["growth_frequency_s_inv"]
                ),
                "response_applications": int(diagnostic["iterations"]),
                "operator_residual_L1": float(
                    diagnostic["operator_residual_L1"]
                ),
                "number_balance_residual": float(
                    diagnostic["number_balance_residual"]
                ),
                "negative_population_mass": float(
                    diagnostic["negative_population_mass"]
                ),
                "tail_probability": float(diagnostic["tail_probability"]),
                "outer_acceleration_flux_fraction": float(
                    diagnostic["outer_acceleration_flux_fraction"]
                ),
                "estimated_peak_memory_bytes": int(
                    diagnostic["estimated_peak_memory_bytes"]
                ),
                "growth_root_iterations": int(
                    timing["growth_root_iterations"]
                ),
                "growth_response_system_builds": int(
                    timing["growth_response_system_builds"]
                ),
                "growth_secant_steps": int(
                    timing["growth_secant_steps"]
                ),
                "response_system_build_s": float(
                    timing["response_system_build_s"]
                ),
                "response_build_blas_threads": int(
                    timing["response_build_blas_threads"]
                ),
                "perron_ritz_pairs": int(timing["response_ritz_pairs"]),
                "perron_cone_compatible_pairs": int(
                    timing["response_cone_compatible_ritz_pairs"]
                ),
                "quality_gates_passed": quality,
            },
            case,
        )
    except PropagatorConvergenceError as exc:
        diagnostic = exc.diagnostics.as_dict()
        return (
            {
                "family": family,
                "profile": profile,
                "E_over_N_Td": float(config.run.e_over_n_Td[0]),
                "status": "not_converged",
                "elapsed_s": perf_counter() - started,
                "stop_reason": diagnostic["stop_reason"],
                "response_applications": int(diagnostic["iterations"]),
                "quality_gates_passed": False,
                "error": str(exc),
            },
            None,
        )
    except Exception as exc:
        return (
            {
                "family": family,
                "profile": profile,
                "E_over_N_Td": float(config.run.e_over_n_Td[0]),
                "status": "error",
                "elapsed_s": perf_counter() - started,
                "quality_gates_passed": False,
                "error": f"{type(exc).__name__}: {exc}",
            },
            None,
        )


def _relative(left: float, right: float) -> float:
    return abs(float(left) - float(right)) / max(
        abs(float(left)),
        abs(float(right)),
        1.0e-300,
    )


def _report_progress(record: dict[str, Any]) -> None:
    print(
        f"{record['family']}/{record['profile']}/"
        f"{record['E_over_N_Td']:g}Td: {record['status']} "
        f"({record['elapsed_s']:.2f}s)",
        flush=True,
    )


def _eedf_l1(left: Any, right: Any) -> float:
    energy = np.asarray(right.energy_eV, dtype=float)
    widths = np.asarray(right.energy_widths_eV, dtype=float)
    projected = np.interp(
        energy,
        np.asarray(left.energy_eV, dtype=float),
        np.asarray(left.eedf, dtype=float),
        left=0.0,
        right=0.0,
    )
    return float(np.sum(np.abs(projected - right.eedf) * widths))


def _comparison(
    left: Any | None,
    right: Any | None,
    *,
    scalar_limit: float,
    eedf_limit: float,
    gate_growth: bool,
) -> dict[str, Any]:
    if left is None or right is None:
        return {"status": "unavailable", "passed": False}
    left_growth = float(
        left.diagnostics["propagator"]["growth_frequency_s_inv"]
    )
    right_growth = float(
        right.diagnostics["propagator"]["growth_frequency_s_inv"]
    )
    values = {
        "mean_energy": _relative(left.mean_energy_eV, right.mean_energy_eV),
        "drift_velocity": _relative(
            left.drift_velocity_m_s,
            right.drift_velocity_m_s,
        ),
        "growth_frequency": _relative(left_growth, right_growth),
    }
    gated = ["mean_energy", "drift_velocity"]
    if gate_growth:
        gated.append("growth_frequency")
    eedf = _eedf_l1(left, right)
    return {
        "status": "available",
        "relative_differences": values,
        "gated_scalars": gated,
        "eedf_weighted_L1": eedf,
        "scalar_limit": scalar_limit,
        "eedf_limit": eedf_limit,
        "passed": bool(
            max(values[name] for name in gated) <= scalar_limit
            and eedf <= eedf_limit
        ),
    }


def _inelastic_weak_transfer_refinement(base: SwarmConfig) -> dict[str, Any]:
    cross_sections = load_active_mixture_inputs(
        base.cross_sections,
        base.conditions,
    )
    internal = build_internal_solver_configs(
        base.solvers,
        base.physics,
    ).propagator
    rows: list[dict[str, Any]] = []
    errors: dict[str, list[float]] = {}
    for energy_cells in (64, 128, 256):
        config = copy.deepcopy(internal)
        config.energy_cells = energy_cells
        config.polar_cells = 8
        grid = build_propagator_grid(config, cross_sections)
        collisions = build_collision_operator(
            base,
            cross_sections,
            grid,
        )
        for transfer in collisions.inelastic:
            if transfer.daughter_count == 0:
                continue
            number = np.asarray(
                transfer.gain_matrix_s_inv.sum(axis=0)
            ).reshape(-1)
            expected_number = (
                transfer.daughter_count * transfer.frequency_s_inv
            )
            number_scale = max(float(np.max(expected_number)), 1.0)
            number_error = float(
                np.max(np.abs(number - expected_number)) / number_scale
            )
            kinematic = (
                transfer.incident_energy_rate_eV_s_inv
                - transfer.daughter_energy_rate_eV_s_inv
                - transfer.threshold_eV * transfer.frequency_s_inv
            )
            energy_scale = max(
                float(
                    np.max(
                        np.abs(transfer.incident_energy_rate_eV_s_inv)
                    )
                ),
                1.0,
            )
            kinematic_error = float(
                np.max(np.abs(kinematic)) / energy_scale
            )
            deposited = np.asarray(
                grid.energy_centers_eV @ transfer.gain_matrix_s_inv
            ).reshape(-1)
            discrete_change = (
                deposited
                - grid.energy_centers_eV * transfer.frequency_s_inv
            )
            threshold_change = (
                -transfer.threshold_eV * transfer.frequency_s_inv
            )
            discrete_error = float(
                np.sum(np.abs(discrete_change - threshold_change))
                / max(float(np.sum(np.abs(threshold_change))), 1.0e-300)
            )
            errors.setdefault(transfer.process_type, []).append(
                discrete_error
            )
            rows.append(
                {
                    "energy_cells": energy_cells,
                    "process_type": transfer.process_type,
                    "daughter_number_relative_error": number_error,
                    "kinematic_energy_relative_error": kinematic_error,
                    "discrete_energy_balance_relative_error": discrete_error,
                }
            )
    contraction = {
        process: [values[index + 1] / values[index] for index in range(2)]
        for process, values in errors.items()
    }
    passed = bool(
        all(row["daughter_number_relative_error"] <= 3.0e-13 for row in rows)
        and all(row["kinematic_energy_relative_error"] <= 5.0e-13 for row in rows)
        and all(ratio <= 0.27 for values in contraction.values() for ratio in values)
    )
    return {"rows": rows, "refinement_ratios": contraction, "passed": passed}


def run_qualification(*, quick: bool) -> dict[str, Any]:
    starting_source_state = _qualification_source_state()
    isotropic_path, maxent_path, mixture_path = QUALIFICATION_INPUT_PATHS[:3]
    isotropic = load_config(isotropic_path)
    maxent = load_config(maxent_path)
    mixture = load_config(mixture_path)
    records: list[dict[str, Any]] = []
    cases: dict[tuple[str, str, float], Any | None] = {}

    profiles = {"smoke": (64, 8)} if quick else GRID_PROFILES
    maxent_fields = (0.1,) if quick else LOW_FIELD_TD
    inelastic_fields = (100.0,) if quick else INELASTIC_FIELD_TD
    for family, base, fields in (
        ("maxent_p1_explicit_integrals", maxent, maxent_fields),
        ("argon_inelastic", isotropic, inelastic_fields),
    ):
        for profile, grid in profiles.items():
            for field in fields:
                record, case = _run_case(
                    _variant(base, field, grid),
                    family=family,
                    profile=profile,
                )
                records.append(record)
                _report_progress(record)
                cases[(family, profile, field)] = case

    refinement: list[dict[str, Any]] = []
    if not quick:
        inelastic_growth_samples = [
            abs(
                float(
                    cases[("argon_inelastic", profile, field)].diagnostics[
                        "propagator"
                    ]["growth_frequency_s_inv"]
                )
            )
            for profile in GRID_PROFILES
            for field in INELASTIC_FIELD_TD
            if cases[("argon_inelastic", profile, field)] is not None
        ]
        inelastic_growth_peak = max(inelastic_growth_samples, default=0.0)
        for family, fields in (
            ("maxent_p1_explicit_integrals", LOW_FIELD_TD),
            ("argon_inelastic", INELASTIC_FIELD_TD),
        ):
            for field in fields:
                field_cases = (
                    cases[(family, "medium", field)],
                    cases[(family, "fine", field)],
                )
                gate_growth = bool(
                    family == "argon_inelastic"
                    and all(case is not None for case in field_cases)
                    and max(
                        abs(
                            float(
                                case.diagnostics["propagator"][
                                    "growth_frequency_s_inv"
                                ]
                            )
                        )
                        for case in field_cases
                    )
                    >= 0.01 * inelastic_growth_peak
                )
                refinement.append(
                    {
                        "family": family,
                        "E_over_N_Td": field,
                        **_comparison(
                            cases[(family, "medium", field)],
                            cases[(family, "fine", field)],
                            scalar_limit=SCALAR_REFINEMENT_LIMIT,
                            eedf_limit=EEDF_L1_REFINEMENT_LIMIT,
                            gate_growth=gate_growth,
                        ),
                        "growth_major_peak_fraction": 0.01,
                    }
                )

    mixture_rows: list[dict[str, Any]] = []
    mixture_cases: dict[str, Any | None] = {}
    if not quick:
        for profile, grid in MIXTURE_PROFILES.items():
            record, case = _run_case(
                _variant(mixture, 100.0, grid),
                family="ar_o2_mixture",
                profile=profile,
            )
            records.append(record)
            _report_progress(record)
            mixture_cases[profile] = case
        mixture_rows.append(
            {
                "E_over_N_Td": 100.0,
                **_comparison(
                    mixture_cases["medium"],
                    mixture_cases["fine"],
                    scalar_limit=SCALAR_REFINEMENT_LIMIT,
                    eedf_limit=EEDF_L1_REFINEMENT_LIMIT,
                    gate_growth=True,
                ),
            }
        )

    boundary_rows: list[dict[str, Any]] = []
    if not quick:
        boundary_cases: dict[tuple[float, float], Any | None] = {}
        for field in INELASTIC_FIELD_TD:
            for ceiling in ENERGY_CEILINGS_EV:
                profile = f"ceiling_{ceiling:g}eV"
                record, case = _run_case(
                    _variant(
                        isotropic,
                        field,
                        (64, 8),
                        maximum_energy_eV=ceiling,
                    ),
                    family="energy_ceiling",
                    profile=profile,
                )
                records.append(record)
                _report_progress(record)
                boundary_cases[(field, ceiling)] = case
            reference = boundary_cases[(field, ENERGY_CEILINGS_EV[-1])]
            for ceiling in ENERGY_CEILINGS_EV[:-1]:
                boundary_rows.append(
                    {
                        "E_over_N_Td": field,
                        "ceiling_eV": ceiling,
                        "reference_ceiling_eV": ENERGY_CEILINGS_EV[-1],
                        **_comparison(
                            boundary_cases[(field, ceiling)],
                            reference,
                            scalar_limit=BOUNDARY_SCALAR_LIMIT,
                            eedf_limit=BOUNDARY_EEDF_L1_LIMIT,
                            gate_growth=True,
                        ),
                    }
                )

    weak_transfer = _inelastic_weak_transfer_refinement(isotropic)
    run_gate = bool(records) and all(
        row["quality_gates_passed"] for row in records
    )
    refinement_gate = bool(quick) or (
        bool(refinement) and all(row["passed"] for row in refinement)
    )
    mixture_gate = bool(quick) or (
        bool(mixture_rows) and all(row["passed"] for row in mixture_rows)
    )
    boundary_gate = bool(quick) or (
        bool(boundary_rows) and all(row["passed"] for row in boundary_rows)
    )
    passed = bool(
        run_gate
        and refinement_gate
        and mixture_gate
        and boundary_gate
        and weak_transfer["passed"]
    )
    blocking = []
    for name, value in (
        ("case_quality", run_gate),
        ("medium_fine_refinement", refinement_gate),
        ("mixed_gas_angular_refinement", mixture_gate),
        ("energy_ceiling_independence", boundary_gate),
        ("inelastic_weak_transfer", weak_transfer["passed"]),
    ):
        if not value:
            blocking.append(name)

    _require_unchanged_qualification_sources(starting_source_state)
    return {
        "schema": "swarm.propagator_p1_deterministic_qualification.v1",
        "scope": {
            "phase": "P0_P1",
            "quick": quick,
            "included_solver": "propagator",
            "excluded": [
                "two_term",
                "multi_term",
                "monte_carlo",
                "COMSOL",
                "P2_bulk_transport_and_diffusion",
            ],
            "field_model": "homogeneous_dc_B0",
            "grid_profiles": GRID_PROFILES,
            "scalar_refinement_limit": SCALAR_REFINEMENT_LIMIT,
            "eedf_weighted_L1_limit": EEDF_L1_REFINEMENT_LIMIT,
        },
        "environment": {
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "git_head": _git_value("rev-parse", "HEAD"),
            "git_dirty": bool(_git_value("status", "--porcelain")),
            "implementation_fingerprint": starting_source_state[
                "implementation_fingerprint"
            ],
            "python": sys.version,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
            "processor": platform.processor()
            or os.environ.get("PROCESSOR_IDENTIFIER", ""),
            "logical_cpu_count": os.cpu_count(),
            "input_sha256": starting_source_state["input_sha256"],
        },
        "runs": records,
        "medium_fine_refinement": refinement,
        "mixed_gas_angular_refinement": mixture_rows,
        "energy_ceiling_independence": boundary_rows,
        "inelastic_weak_transfer_refinement": weak_transfer,
        "structural_verification_tests": [
            "test_fused_redheffer_composition_matches_literal_block_equations",
            "test_shell_error_control_bounds_independent_continuous_references",
            "test_characteristic_origin_is_stable_under_refinement_for_maxent_argon",
            "test_characteristic_origin_obeys_collisionless_and_rate_scaling_limits",
            "test_perron_branch_matches_dense_reference_and_is_seed_independent",
            "test_propagator_mixture_is_density_scaled_and_order_independent",
        ],
        "decision": {
            "case_quality_passed": run_gate,
            "medium_fine_refinement_passed": refinement_gate,
            "mixed_gas_angular_refinement_passed": mixture_gate,
            "energy_ceiling_independence_passed": boundary_gate,
            "inelastic_weak_transfer_passed": weak_transfer["passed"],
            "p1_deterministic_core_qualified": passed,
            "blocking_gates": blocking,
            "release_claim": (
                "deterministic_P0_P1_scope_only; no stochastic or P2 claim"
            ),
        },
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run Propagator-only deterministic P0/P1 qualification."
    )
    parser.add_argument("--quick", action="store_true")
    parser.add_argument(
        "--output",
        type=Path,
    )
    args = parser.parse_args(argv)
    output = (
        args.output.resolve()
        if args.output is not None
        else (
            QUICK_QUALIFICATION_OUTPUT
            if args.quick
            else FULL_QUALIFICATION_OUTPUT
        ).resolve()
    )
    if args.quick and output == FULL_QUALIFICATION_OUTPUT.resolve():
        parser.error("--quick cannot overwrite the formal qualification artifact")
    payload = run_qualification(quick=bool(args.quick))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(output)
    print(json.dumps(payload["decision"], indent=2, sort_keys=True))
    return 0 if payload["decision"]["p1_deterministic_core_qualified"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
