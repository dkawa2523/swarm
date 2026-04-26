from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import yaml

from electron_swarm.collisions.electron_electron import mean_energy_eV, relaxation_target
from electron_swarm.collisions.states import augment_cross_sections_from_config
from electron_swarm.core.cross_sections import CrossSectionProcess, CrossSectionSet, ProcessType
from electron_swarm.diagnostics.common import eedf_quality_metrics
from electron_swarm.grids.energy import build_energy_grid


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
