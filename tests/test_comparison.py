from __future__ import annotations

import numpy as np
import pytest

from electron_swarm.core.config import ComparisonConfig
from electron_swarm.core.results import SwarmCaseResult, SwarmRunResult
from electron_swarm.core.transport import ElectronTransport
from electron_swarm.orchestration.comparison import build_comparison_summary
from electron_swarm.core.config_parser import load_config_from_raw
from electron_swarm.core.physical_context import physical_context
from product_helpers import base_product_config
from pathlib import Path


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
            definition="flux",
            gas_number_density_m3=1.0,
            drift_velocity_m_s=1.0,
            reduced_mobility_m2_V_s_m3=1.0,
            reduced_diffusion_L_m2_s_m3=1.0,
            reduced_diffusion_T_m2_s_m3=1.0,
        ),
        energy_eV=np.asarray(energy_eV),
        energy_widths_eV=np.asarray(widths_eV),
        eedf=np.asarray(eedf),
        physical_context=physical_context(
            load_config_from_raw(
                base_product_config(Path.cwd()), Path.cwd() / "synthetic.yaml"
            ),
            "synthetic",
            50.0,
        ),
        metadata={
            "angular_model": "isotropic",
            "angular_moment_source": "ordinary_integral_xs",
            "exact_dcs_based": False,
            "ordinary_integral_xs_closure": True,
            "angular_scattering_assumption": "isotropic:zero",
        },
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


def test_comparison_identifies_different_first_moment_assumptions() -> None:
    cases = [
        _case(
            solver,
            energy_eV=[0.5, 1.5],
            widths_eV=[1.0, 1.0],
            eedf=[0.5, 0.5],
        )
        for solver in ("two_term", "monte_carlo")
    ]
    for case in cases:
        case.metadata.update(
            {
                "angular_model": "isotropic",
                "angular_moment_source": "ordinary_integral_xs",
                "exact_dcs_based": False,
                "ordinary_integral_xs_closure": True,
                "angular_scattering_assumption": "isotropic:zero",
            }
        )
    cases[0].metadata["angular_scattering_assumption"] += (
        ";first_moment=1-sigma_m/sigma_0:Ar"
    )
    comparison = ComparisonConfig(
        enabled=True,
        reference_solver="monte_carlo",
        candidate_solvers=["two_term"],
        required=False,
    )
    result = SwarmRunResult(cases=cases)
    [row] = build_comparison_summary(result, comparison)
    assert row["angular_model_status"] == "mismatch"
    comparison.required = True
    with pytest.raises(ValueError, match="angular_model_status='mismatch'"):
        build_comparison_summary(result, comparison)


def test_comparison_does_not_report_stale_transport_as_current() -> None:
    cases = [
        _case(solver, energy_eV=[0.5, 1.5], widths_eV=[1.0, 1.0], eedf=[0.5, 0.5])
        for solver in ("two_term", "multi_term")
    ]
    cases[1].metadata["electron_electron_transport_stale"] = True
    [row] = build_comparison_summary(
        SwarmRunResult(cases),
        ComparisonConfig(
            enabled=True, reference_solver="two_term", candidate_solvers=["multi_term"]
        ),
    )
    assert row["mobility_relative_difference"] is None
    assert row["mobility_m2_V_s_status"] == "candidate_transport_stale"
    assert row["mean_energy_eV_relative_difference"] == pytest.approx(0.0)
    assert row["candidate_numerical_status"] == "postprocessed_transport_stale"


def test_comparison_defaults_choose_a_present_deterministic_reference() -> None:
    cases = [
        _case(solver, energy_eV=[0.5, 1.5], widths_eV=[1.0, 1.0], eedf=[0.5, 0.5])
        for solver in ("two_term", "multi_term")
    ]
    result = SwarmRunResult(cases)
    [row] = build_comparison_summary(result, ComparisonConfig(enabled=True))
    assert row["reference_solver"] == "two_term"
    assert row["candidate_solver"] == "multi_term"
    assert result.metadata["comparison"]["available_pairs"] == 1


def test_optional_comparison_reports_missing_reference_without_switching_it() -> None:
    result = SwarmRunResult(
        [
            _case(
                "multi_term",
                energy_eV=[0.5, 1.5],
                widths_eV=[1.0, 1.0],
                eedf=[0.5, 0.5],
            )
        ]
    )
    result.metadata["configuration"] = {
        "run": {
            "solvers": [
                {"id": solver, "enabled": True} for solver in ("two_term", "multi_term")
            ],
            "case_prefix": "case",
            "e_over_n_Td": [50.0],
        }
    }
    [row] = build_comparison_summary(result, ComparisonConfig(enabled=True))
    assert row["reference_solver"] == "two_term"
    assert row["comparison_status"] == "reference_missing"
    assert row["mean_energy_eV_relative_difference"] is None
    assert result.metadata["comparison"]["status"] == "unavailable"


def test_required_comparison_covers_candidate_only_cases() -> None:
    cases = [
        _case(solver, energy_eV=[0.5, 1.5], widths_eV=[1.0, 1.0], eedf=[0.5, 0.5])
        for solver in ("two_term", "multi_term", "multi_term")
    ]
    cases[-1].case_id = "case_0001"
    result = SwarmRunResult(cases)
    with pytest.raises(ValueError, match="reference results"):
        build_comparison_summary(
            result,
            ComparisonConfig(
                enabled=True,
                reference_solver="two_term",
                candidate_solvers=["multi_term"],
                required=True,
            ),
        )
    assert result.metadata["comparison"]["status"] == "partial"
    assert result.metadata["comparison"]["expected_pairs"] == 2


def test_nonzero_against_zero_reference_has_no_finite_relative_difference() -> None:
    cases = [
        _case(solver, energy_eV=[0.5, 1.5], widths_eV=[1.0, 1.0], eedf=[0.5, 0.5])
        for solver in ("two_term", "multi_term")
    ]
    cases[1].net_ionization_frequency_s = 1.0
    [row] = build_comparison_summary(
        SwarmRunResult(cases), ComparisonConfig(enabled=True)
    )
    assert row["net_ionization_frequency_relative_difference"] is None
    assert row["net_ionization_frequency_s_status"] == "reference_zero"
