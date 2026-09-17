from __future__ import annotations

import numpy as np
import pytest

from electron_swarm.core.config import ComparisonConfig
from electron_swarm.core.results import SwarmCaseResult, SwarmRunResult
from electron_swarm.core.transport import ElectronTransport
from electron_swarm.orchestration.comparison import build_comparison_summary


def _case(
    solver: str,
    *,
    energy_eV: list[float],
    widths_eV: list[float],
    eedf: list[float],
) -> SwarmCaseResult:
    return SwarmCaseResult(
        solver=solver,
        case_id="case_0000",
        e_over_n_Td=50.0,
        mean_energy_eV=1.0,
        net_ionization_frequency_s=0.0,
        effective_townsend_m2=0.0,
        transport=ElectronTransport(
            definition="test",
            gas_number_density_m3=1.0,
            drift_velocity_m_s=1.0,
            reduced_mobility_m2_V_s_m3=1.0,
            reduced_diffusion_L_m2_s_m3=1.0,
            reduced_diffusion_T_m2_s_m3=1.0,
        ),
        energy_eV=np.asarray(energy_eV),
        energy_widths_eV=np.asarray(widths_eV),
        eedf=np.asarray(eedf),
    )


def test_required_comparison_rejects_missing_candidate_case() -> None:
    result = SwarmRunResult(
        cases=[
            _case(
                "two_term",
                energy_eV=[0.5, 1.5],
                widths_eV=[1.0, 1.0],
                eedf=[0.5, 0.5],
            )
        ]
    )
    comparison = ComparisonConfig(
        enabled=True,
        reference_solver="two_term",
        candidate_solvers=["multi_term"],
        required=True,
    )

    with pytest.raises(ValueError, match="missing candidate results"):
        build_comparison_summary(result, comparison)


def test_eedf_l1_uses_reported_cells_conservatively() -> None:
    result = SwarmRunResult(
        cases=[
            _case(
                "two_term",
                energy_eV=[0.5, 1.5],
                widths_eV=[1.0, 1.0],
                eedf=[0.5, 0.5],
            ),
            _case(
                "multi_term",
                energy_eV=[0.25, 0.75, 1.25, 1.75],
                widths_eV=[0.5, 0.5, 0.5, 0.5],
                eedf=[0.5, 0.5, 0.5, 0.5],
            ),
        ]
    )
    comparison = ComparisonConfig(
        enabled=True,
        reference_solver="two_term",
        candidate_solvers=["multi_term"],
        required=True,
    )

    [row] = build_comparison_summary(result, comparison)

    assert row["eedf_l1_error"] == pytest.approx(0.0)
