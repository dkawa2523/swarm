"""Finite-volume refinement comparisons shared by qualification runners."""

from __future__ import annotations

from typing import Any

from electron_swarm.grids.energy import cell_density_l1, reported_cell_edges


def compare_propagator_cases(
    left: Any | None,
    right: Any | None,
    *,
    scalar_limit: float,
    eedf_limit: float,
    gate_growth: bool,
) -> dict[str, Any]:
    if left is None or right is None:
        return {"status": "unavailable", "passed": False}

    def relative(a: float, b: float) -> float:
        return abs(float(a) - float(b)) / max(abs(float(a)), abs(float(b)), 1e-300)

    values = {
        "mean_energy": relative(left.mean_energy_eV, right.mean_energy_eV),
        "drift_velocity": relative(left.drift_velocity_m_s, right.drift_velocity_m_s),
        "growth_frequency": relative(
            left.diagnostics["propagator"]["growth_frequency_s_inv"],
            right.diagnostics["propagator"]["growth_frequency_s_inv"],
        ),
    }
    gated = ["mean_energy", "drift_velocity"]
    if gate_growth:
        gated.append("growth_frequency")
    eedf = cell_density_l1(
        reported_cell_edges(left.energy_eV, left.energy_widths_eV), left.eedf,
        reported_cell_edges(right.energy_eV, right.energy_widths_eV), right.eedf,
    )
    return {
        "status": "available",
        "relative_differences": values,
        "gated_scalars": gated,
        "eedf_weighted_L1": eedf,
        "scalar_limit": scalar_limit,
        "eedf_limit": eedf_limit,
        "passed": bool(max(values[name] for name in gated) <= scalar_limit
                       and eedf <= eedf_limit),
    }
