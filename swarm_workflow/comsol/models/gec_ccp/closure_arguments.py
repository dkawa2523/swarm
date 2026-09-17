"""Mean-energy argument contract shared by Java generation and audits."""

from __future__ import annotations

from swarm_workflow.comsol.input.mean_energy import MeanEnergyArgument
from swarm_workflow.comsol.models.gec_ccp.contracts import (
    GecCcpMapping,
    GecCcpWorkflowError,
)


def smooth_log_energy_argument(
    log_mean_energy: str,
    *,
    minimum_eV: float,
    maximum_eV: float,
) -> str:
    """Map any Newton trial energy smoothly inside a log-energy table."""

    try:
        return MeanEnergyArgument(minimum_eV, maximum_eV).expression(log_mean_energy)
    except ValueError as exc:
        raise GecCcpWorkflowError(str(exc)) from exc


def closure_argument_range(
    mapping: GecCcpMapping,
    support: list[float],
) -> tuple[float, float]:
    """Return the active coefficient range, including an explicit test floor."""

    minimum_eV, maximum_eV = (float(value) for value in support)
    floor_eV = mapping.run.validation_mean_energy_floor_eV
    if floor_eV is not None:
        minimum_eV = max(minimum_eV, floor_eV)
    if not 0.0 < minimum_eV < maximum_eV:
        raise GecCcpWorkflowError(
            "validation mean-energy floor must remain inside every active "
            "coefficient support"
        )
    return minimum_eV, maximum_eV
