"""Bounded continuation of deterministic solver mean-energy support."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import math
from typing import Sequence, cast

import numpy as np

from electron_swarm.core.config import RequestedSolverConfig, SolverId, SwarmConfig
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.runner import run as _run_product_config


_DETERMINISTIC_SOLVERS = frozenset({"two_term", "multi_term", "propagator"})


@dataclass(frozen=True, slots=True)
class MeanEnergyExtensionRequest:
    """One field point proposed from the calculated upper support."""

    current_max_mean_energy_eV: float
    required_max_mean_energy_eV: float
    guarded_target_mean_energy_eV: float
    suggested_E_over_N_Td: float


@dataclass(frozen=True, slots=True)
class MeanEnergyContinuationResult:
    """Calculated extension cases and whether they reached the guarded target."""

    solver_id: SolverId
    guarded_target_mean_energy_eV: float
    initial_max_mean_energy_eV: float
    final_max_mean_energy_eV: float
    reached_target: bool
    extension_cases: tuple[SwarmCaseResult, ...]


def propose_upper_mean_energy_extension(
    e_over_n_Td: Sequence[float] | np.ndarray,
    mean_energy_eV: Sequence[float] | np.ndarray,
    required_max_mean_energy_eV: float,
    *,
    relative_guard: float = 0.1,
    maximum_field_step_factor: float = 4.0,
) -> MeanEnergyExtensionRequest | None:
    """Propose one bounded upper E/N continuation step in log-log space.

    Only calculated ``(E/N, mean energy)`` pairs determine the local inverse.
    The result is a new solver input, never a synthetic EEDF or coefficient.
    """

    fields = np.asarray(e_over_n_Td, dtype=float)
    means = np.asarray(mean_energy_eV, dtype=float)
    required = float(required_max_mean_energy_eV)
    guard = float(relative_guard)
    step_cap = float(maximum_field_step_factor)
    if (
        fields.ndim != 1
        or means.shape != fields.shape
        or len(fields) < 2
        or np.any(~np.isfinite(fields))
        or np.any(~np.isfinite(means))
        or np.any(fields <= 0.0)
        or np.any(means <= 0.0)
        or np.any(np.diff(fields) <= 0.0)
        or not math.isfinite(required)
        or required <= 0.0
        or not math.isfinite(guard)
        or guard < 0.0
        or not math.isfinite(step_cap)
        or step_cap <= 1.0
    ):
        raise ValueError(
            "mean-energy continuation requires positive finite values and "
            "strictly increasing E/N"
        )

    guarded = required * (1.0 + guard)
    if not math.isfinite(guarded):
        raise ValueError("guarded mean-energy target must be finite")
    envelope_indices: list[int] = []
    record = -math.inf
    for index, value in enumerate(means):
        if float(value) > record:
            envelope_indices.append(index)
            record = float(value)
    current_max = float(means[envelope_indices[-1]])
    if guarded <= current_max:
        return None
    if len(envelope_indices) < 2:
        raise ValueError("upper mean-energy envelope is not locally invertible")
    if envelope_indices[-1] != len(fields) - 1:
        raise ValueError(
            "upper mean-energy endpoint is not locally invertible; calculate a "
            "new upper E/N anchor before extrapolating"
        )

    previous_index, current_index = envelope_indices[-2:]
    slope = math.log(means[current_index] / means[previous_index]) / math.log(
        fields[current_index] / fields[previous_index]
    )
    if not math.isfinite(slope) or slope <= 0.0:
        raise ValueError("upper mean-energy response is not locally invertible")

    # Compare logarithmic steps before exponentiation so a nearly flat local
    # response cannot overflow before the explicit field-step bound is applied.
    predicted_log_step = math.log(guarded / current_max) / slope
    bounded_log_step = min(predicted_log_step, math.log(step_cap))
    suggested = float(fields[-1]) * math.exp(bounded_log_step)
    suggested = max(suggested, math.nextafter(float(fields[-1]), math.inf))
    return MeanEnergyExtensionRequest(
        current_max_mean_energy_eV=current_max,
        required_max_mean_energy_eV=required,
        guarded_target_mean_energy_eV=guarded,
        suggested_E_over_N_Td=suggested,
    )


def continue_deterministic_mean_energy_support(
    config: SwarmConfig,
    solver_id: SolverId,
    e_over_n_Td: Sequence[float] | np.ndarray,
    mean_energy_eV: Sequence[float] | np.ndarray,
    required_max_mean_energy_eV: float,
    *,
    relative_guard: float = 0.1,
    maximum_field_step_factor: float = 4.0,
    maximum_steps: int = 2,
    maximum_e_over_n_Td: float | None = None,
) -> MeanEnergyContinuationResult:
    """Run bounded continuation cases through the selected product solver.

    The input config is copied and only its selected solver, E/N point, and
    case prefix are changed.  Monte Carlo is deliberately excluded because an
    isolated trajectory run cannot replace its configured replica campaign.
    """

    solver_name = str(solver_id)
    if solver_name == "monte_carlo":
        raise ValueError(
            "monte_carlo mean-energy support must be extended by its statistical "
            "campaign"
        )
    if solver_name not in _DETERMINISTIC_SOLVERS:
        raise ValueError(f"unsupported deterministic solver {solver_name!r}")
    if isinstance(maximum_steps, bool) or not isinstance(maximum_steps, int):
        raise ValueError("maximum_steps must be a positive integer")
    if maximum_steps <= 0:
        raise ValueError("maximum_steps must be a positive integer")

    fields = [float(value) for value in e_over_n_Td]
    means = [float(value) for value in mean_energy_eV]
    # Validate even when the existing support already satisfies the target.
    initial_request = propose_upper_mean_energy_extension(
        fields,
        means,
        required_max_mean_energy_eV,
        relative_guard=relative_guard,
        maximum_field_step_factor=maximum_field_step_factor,
    )
    initial_max = max(means)
    guarded_target = (
        initial_request.guarded_target_mean_energy_eV
        if initial_request is not None
        else float(required_max_mean_energy_eV) * (1.0 + float(relative_guard))
    )

    field_limit = math.inf
    if maximum_e_over_n_Td is not None:
        field_limit = float(maximum_e_over_n_Td)
        if not math.isfinite(field_limit) or field_limit <= fields[-1]:
            raise ValueError(
                "maximum_e_over_n_Td must exceed the current upper E/N anchor"
            )

    extension_cases: list[SwarmCaseResult] = []
    request = initial_request
    for step in range(maximum_steps):
        if request is None:
            break
        next_field = min(request.suggested_E_over_N_Td, field_limit)
        if next_field <= fields[-1]:
            break

        step_config = deepcopy(config)
        step_config.run.solvers = [
            RequestedSolverConfig(id=cast(SolverId, solver_name), enabled=True)
        ]
        step_config.run.e_over_n_Td = [next_field]
        step_config.run.case_prefix = (
            f"{config.run.case_prefix}_mean_energy_support_{step:02d}"
        )
        calculated = _run_product_config(step_config, write=False).cases
        if len(calculated) != 1:
            raise RuntimeError(
                f"{solver_name} support continuation returned "
                f"{len(calculated)} cases for one E/N point"
            )
        case = calculated[0]
        if (
            case.solver != solver_name
            or not math.isclose(
                float(case.e_over_n_Td), next_field, rel_tol=0.0, abs_tol=0.0
            )
            or not math.isfinite(float(case.mean_energy_eV))
            or float(case.mean_energy_eV) <= 0.0
        ):
            raise RuntimeError(
                f"{solver_name} support continuation did not return a positive "
                "finite mean-energy point"
            )
        extension_cases.append(case)
        fields.append(next_field)
        means.append(float(case.mean_energy_eV))
        request = propose_upper_mean_energy_extension(
            fields,
            means,
            required_max_mean_energy_eV,
            relative_guard=relative_guard,
            maximum_field_step_factor=maximum_field_step_factor,
        )

    final_max = max(means)
    return MeanEnergyContinuationResult(
        solver_id=cast(SolverId, solver_name),
        guarded_target_mean_energy_eV=guarded_target,
        initial_max_mean_energy_eV=initial_max,
        final_max_mean_energy_eV=final_max,
        reached_target=final_max >= guarded_target,
        extension_cases=tuple(extension_cases),
    )


__all__ = [
    "MeanEnergyContinuationResult",
    "MeanEnergyExtensionRequest",
    "continue_deterministic_mean_energy_support",
    "propose_upper_mean_energy_extension",
]
