from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import textwrap

import numpy as np
import yaml

from electron_swarm import load_config, run_from_config
from electron_swarm.collisions.electron_electron import mean_energy_eV, relaxation_target
from electron_swarm.collisions.states import augment_cross_sections_from_config
from electron_swarm.core.cross_sections import CrossSectionProcess, CrossSectionSet, ProcessType
from electron_swarm.diagnostics.common import eedf_quality_metrics
from electron_swarm.grids.energy import build_energy_grid


ROOT = Path(__file__).resolve().parents[1]
MINIMAL_XS = ROOT / "examples" / "cross_sections" / "argon_minimal.csv"


def test_eedf_quality_metrics_are_scalar_and_normalized() -> None:
    energy = np.array([0.5, 1.5, 2.5])
    widths = np.array([1.0, 1.0, 1.0])
    eedf = np.array([0.2, 0.5, 0.3])

    quality = eedf_quality_metrics(energy, eedf, widths)

    assert quality.normalization_error < 1.0e-12
    assert quality.negative_fraction == 0.0
    assert quality.tail_probability >= 0.0


def test_threshold_aware_energy_grid_adds_points_near_process_threshold() -> None:
    process = CrossSectionProcess(
        species="Ar",
        process="excitation_4s",
        process_type=ProcessType.EXCITATION,
        threshold_eV=2.0,
        energy_eV=np.array([0.0, 1.0, 2.0, 5.0]),
        cross_section_m2=np.array([0.0, 0.0, 1.0e-20, 2.0e-20]),
    )
    xs = CrossSectionSet([process])

    grid = build_energy_grid(
        min_eV=0.1,
        max_eV=5.0,
        n=20,
        spacing="linear",
        cross_sections=xs,
        refine=True,
        threshold_padding_eV=0.1,
        points_per_threshold=5,
    )

    assert len(grid.centers_eV) > 20
    assert np.min(np.abs(grid.centers_eV - 2.0)) < 1.0e-8
    assert grid.metadata["threshold_refined"] is True


def test_state_resolved_config_generates_superelastic_process(tmp_path) -> None:
    process = CrossSectionProcess(
        species="Ar",
        process="excitation_4s",
        process_type=ProcessType.EXCITATION,
        threshold_eV=11.55,
        energy_eV=np.linspace(0.1, 30.0, 64),
        cross_section_m2=np.full(64, 1.0e-20),
    )
    xs = CrossSectionSet([process])
    cfg_path = tmp_path / "case.yaml"
    cfg_path.write_text(
        yaml.safe_dump(
            {
                "state_resolved": {
                    "enabled": True,
                    "generate_superelastic": True,
                    "detailed_balance": "simple",
                    "transitions": [
                        {
                            "species": "Ar",
                            "process": "excitation_4s",
                            "initial_state": "ground",
                            "final_state": "Ar_4s",
                        }
                    ],
                },
                "species_states": {
                    "Ar": {
                        "ground": {"energy_eV": 0.0, "population": 1.0},
                        "Ar_4s": {"energy_eV": 11.55, "population": 1.0e-4},
                    }
                },
            }
        ),
        encoding="utf-8",
    )

    augmented = augment_cross_sections_from_config(xs, SimpleNamespace(source_path=cfg_path))

    generated = [p for p in augmented.processes if p.process_type == ProcessType.SUPERELASTIC]
    assert len(generated) == 1
    assert generated[0].metadata["generated_by"] == "state_resolved_superelastic"
    assert generated[0].metadata["regularized_low_energy_singularity"] is True
    assert np.max(generated[0].cross_section_m2) < 2.0e-22


def test_electron_electron_relaxation_target_is_normalized_and_mean_preserving() -> None:
    energy = np.linspace(0.01, 20.0, 200)
    widths = np.gradient(energy)
    eedf = np.exp(-energy / 3.0)
    eedf = eedf / np.sum(eedf * widths)
    before = mean_energy_eV(energy, widths, eedf)

    target = relaxation_target(energy, widths, eedf, conserve_mean_energy=True)
    after = mean_energy_eV(energy, widths, target)

    assert abs(float(np.sum(target * widths)) - 1.0) < 1.0e-12
    assert abs(after - before) / before < 0.05


def _extension_config_text(
    *, refine: bool, electron_electron: bool, mc_api: str | None = None
) -> str:
    data = {
        "run": {"mode": "all", "e_over_n_Td": [50], "case_prefix": "Ext"},
        "conditions": {
            "gas_temperature_K": 300.0,
            "pressure_Pa": 100.0,
            "gas_mixture": [
                {"species": "Ar", "fraction": 1.0, "mass_amu": 39.948}
            ],
        },
        "cross_sections": {
            "format": "csv",
            "high_energy_extrapolation": "zero",
            "files": [
                {
                    "path": MINIMAL_XS.as_posix(),
                    "species": "Ar",
                    "format": "csv",
                }
            ],
        },
        "boltzmann_two_term": {
            "enabled": True,
            "backend": "native_bolsig",
            "energy_grid": {
                "min_eV": 0.0001,
                "max_eV": 40.0,
                "n": 80,
                "spacing": "quadratic",
            },
            "adaptive_grid": {"enabled": False},
            "convergence": {"max_iterations": 80},
        },
        "multiterm_boltzmann": {
            "enabled": True,
            "method": "moment_closure",
            "lmax": 3,
            "energy_grid": {
                "min_eV": 0.001,
                "max_eV": 50.0,
                "n": 48,
                "spacing": "log_linear",
                "linear_until_eV": 2.0,
            },
        },
        "monte_carlo": {"enabled": bool(mc_api)},
        "output": {
            "directory": "outputs",
            "base_name": "ext",
            "write_plots": False,
        },
    }
    if mc_api:
        data["monte_carlo"]["python_api"] = mc_api
    if refine:
        refinement = {
            "enabled": True,
            "threshold_padding_eV": 0.10,
            "points_per_threshold": 5,
        }
        data["boltzmann_two_term"]["energy_grid"]["refine"] = refinement
        data["multiterm_boltzmann"]["energy_grid"]["refine"] = dict(refinement)
    if electron_electron:
        data["electron_electron"] = {
            "enabled": True,
            "model": "relaxation",
            "relaxation_fraction": 0.07,
            "conserve_mean_energy": True,
        }
    return yaml.safe_dump(data, sort_keys=False)


def test_electron_electron_relaxation_runs_from_config_and_skips_mc(
    tmp_path, monkeypatch
) -> None:
    api_dir = tmp_path / "api"
    api_dir.mkdir()
    (api_dir / "dummy_mc_api.py").write_text(
        textwrap.dedent(
            """
            import numpy as np
            from electron_swarm.core.results import SwarmCaseResult

            def run_swarm(config, cross_sections, **kwargs):
                return [
                    SwarmCaseResult(
                        solver="monte_carlo",
                        case_id="mc_0000",
                        e_over_n_Td=float(config.run.e_over_n_Td[0]),
                        mean_energy_eV=1.0,
                        drift_velocity_m_s=10.0,
                        mobility_m2_V_s=1.0,
                        reduced_mobility_m2_V_s_m3=1.0,
                        diffusion_L_m2_s=1.0,
                        diffusion_T_m2_s=1.0,
                        reduced_diffusion_L_m2_s_m3=1.0,
                        reduced_diffusion_T_m2_s_m3=1.0,
                        net_ionization_frequency_s=0.0,
                        effective_townsend_m2=0.0,
                        energy_eV=np.array([0.0, 1.0]),
                        eedf=np.array([1.0, 0.0]),
                        eepf=np.array([1.0, 0.0]),
                    )
                ]
            """
        ),
        encoding="utf-8",
    )
    monkeypatch.syspath_prepend(str(api_dir))
    cfg_path = tmp_path / "ee.yaml"
    cfg_path.write_text(
        _extension_config_text(
            refine=False,
            electron_electron=True,
            mc_api="dummy_mc_api:run_swarm",
        ),
        encoding="utf-8",
    )

    result = run_from_config(cfg_path, write=False)

    boltzmann_cases = [
        c
        for c in result.cases
        if c.solver in {"boltzmann_two_term", "multiterm_boltzmann"}
    ]
    assert {c.solver for c in boltzmann_cases} == {
        "boltzmann_two_term",
        "multiterm_boltzmann",
    }
    for case in boltzmann_cases:
        assert case.metadata["ee_enabled"] is True
        assert case.metadata["ee_model"] == "relaxation"
        assert case.metadata["ee_model_implementation"] == "relaxation_postprocess"
        assert case.metadata["ee_relaxation_fraction"] == 0.07
        assert case.metadata["ee_transport_recomputed"] is False
        assert case.metadata["ee_transport_stale_after_relaxation"] is True
        assert case.metadata["ee_rates_recomputed"] is True
        assert case.mean_energy_eV == case.metadata["ee_mean_energy_after_eV"]

    [mc_case] = [c for c in result.cases if c.solver == "monte_carlo"]
    assert "ee_enabled" not in mc_case.metadata
    assert mc_case.metadata["ee_relaxation_skipped"] == "non_boltzmann_solver"


def test_energy_grid_refinement_is_optional_and_reported(tmp_path) -> None:
    default_path = tmp_path / "default.yaml"
    default_path.write_text(
        _extension_config_text(refine=False, electron_electron=False),
        encoding="utf-8",
    )
    refined_path = tmp_path / "refined.yaml"
    refined_path.write_text(
        _extension_config_text(refine=True, electron_electron=False),
        encoding="utf-8",
    )

    default = run_from_config(default_path, write=False)
    refined = run_from_config(refined_path, write=False)

    for case in default.cases:
        expected_n = 80 if case.solver == "boltzmann_two_term" else 48
        assert case.metadata["grid_n_cells"] == expected_n
        assert case.metadata["threshold_refined"] is False
        assert "ee_enabled" not in case.metadata

    for case in refined.cases:
        expected_n = 80 if case.solver == "boltzmann_two_term" else 48
        assert case.metadata["grid_n_cells"] > expected_n
        assert case.metadata["threshold_refined"] is True
        assert case.metadata["grid_min_eV"] >= 0.0
        assert case.metadata["grid_max_eV"] > case.metadata["grid_min_eV"]

    cfg = load_config(refined_path)
    assert cfg.boltzmann_two_term.energy_grid.refine.enabled is True
    assert cfg.multiterm_boltzmann.energy_grid.refine.enabled is True
