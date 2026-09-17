from __future__ import annotations

import numpy as np
import pytest

import electron_swarm.solvers.two_term.transport as kinetic
from electron_swarm import load_config, run
from electron_swarm.core.cross_sections import ProcessType, load_cross_sections
from electron_swarm.core.solver_configs import build_internal_solver_configs
from electron_swarm.physics.kinetics import gas_number_density
from electron_swarm.solvers.boltzmann_common.collisions import (
    build_effective_collision_data,
)
from electron_swarm.solvers.boltzmann_common.grid import cell_edges_from_centers
from electron_swarm.solvers.two_term.steady import assemble_native_operator_block
from electron_swarm.solvers.two_term.transport import (
    temporal_growth_effective_momentum_frequency,
)

from product_helpers import ROOT, base_product_config, write_config


def test_pt_effective_momentum_frequency_adds_growth_and_fails_closed() -> None:
    base = np.array([10.0, 20.0])
    assert temporal_growth_effective_momentum_frequency(
        base, 3.0
    ) == pytest.approx([13.0, 23.0])

    with pytest.raises(FloatingPointError, match="must remain finite and positive"):
        temporal_growth_effective_momentum_frequency(base, -10.0)
    with pytest.raises(FloatingPointError, match="finite growth frequency"):
        temporal_growth_effective_momentum_frequency(base, float("nan"))
    with pytest.raises(FloatingPointError, match="finite positive nu_m"):
        temporal_growth_effective_momentum_frequency(
            np.array([0.0, 20.0]), 1.0
        )


def test_dc_growth_uses_converged_pt_frequency_for_flux_transport(
    tmp_path,
) -> None:
    data = base_product_config(tmp_path, ["two_term"])
    data["run"]["e_over_n_Td"] = [500.0]
    data["conditions"]["pressure_Pa"] = 13.3322
    data["cross_sections"]["files"] = [
        {
            "path": (
                ROOT / "examples" / "cross_sections" / "argon_application_library.csv"
            ).as_posix(),
            "species": "Ar",
            "format": "csv",
        }
    ]
    data["physics"]["energy_grid_policy"]["threshold_refinement"] = True
    data["physics"]["energy_grid_policy"]["tail_probability_target"] = 1.0e-9
    data["physics"]["energy_grid_policy"]["max_eV_limit"] = 30000.0
    cfg = load_config(write_config(tmp_path, data, "growth.yaml"))
    [case] = run(cfg, write=False).cases

    diagnostic = case.diagnostics["two_term"]
    correction = diagnostic["temporal_growth_momentum_correction"]
    growth = diagnostic["growth_frequency_s-1"]
    assert correction["applied"] is True
    assert correction["growth_frequency_s-1"] == pytest.approx(growth)
    assert correction["model"] == (
        "nu_m_eff(epsilon)=nu_m(epsilon)+growth_frequency"
    )
    assert growth > 0.0

    # BOLSIG+ 07/2024, temporal growth, equal sharing, same Ar cross sections.
    reference = {
        "mean_energy_eV": 10.3854,
        "muN": 6.63936e23,
        "DeN": 6.24742e24,
        "muenN": 1.02941e24,
        "DenN": 6.64794e24,
    }
    observed = {
        "mean_energy_eV": case.mean_energy_eV,
        "muN": case.reduced_mobility_m2_V_s_m3,
        "DeN": case.reduced_diffusion_L_m2_s_m3,
        "muenN": case.reduced_electron_energy_mobility_m2_V_s_m3,
        "DenN": case.reduced_electron_energy_diffusion_m2_s_m3,
    }
    assert max(
        abs(observed[name] - value) / value
        for name, value in reference.items()
    ) < 2.0e-3

    cross_sections = load_cross_sections(cfg.cross_sections, cfg.conditions)
    internal = build_internal_solver_configs(cfg.solvers, cfg.physics)
    probe_energy = np.array([1.0])
    probe_collisions = build_effective_collision_data(
        cfg,
        cross_sections,
        probe_energy,
        gas_number_density(cfg),
        internal.two_term,
    )
    raw_momentum_sum = sum(
        process.sigma(probe_energy)
        for process in cross_sections.processes
        if process.process_type
        in {
            ProcessType.MOMENTUM,
            ProcessType.EFFECTIVE,
            ProcessType.ELASTIC,
            ProcessType.EXCITATION,
            ProcessType.IONIZATION,
            ProcessType.ATTACHMENT,
            ProcessType.SUPERELASTIC,
        }
    )
    assert probe_collisions.sigma_m == pytest.approx(raw_momentum_sum)
    edges, widths = cell_edges_from_centers(case.energy_eV)
    block = assemble_native_operator_block(
        cfg,
        cross_sections,
        internal.two_term,
        case.e_over_n_Td,
        case.energy_eV,
        edges,
        widths,
    )
    without_pt = kinetic.transport_from_eedf(
        cfg,
        block,
        case.eedf,
        case.e_over_n_Td,
        internal.two_term,
    )
    assert case.reduced_diffusion_L_m2_s_m3 < (
        without_pt.reduced_diffusion_L_m2_s_m3
    )
    assert case.reduced_electron_energy_diffusion_m2_s_m3 < (
        without_pt.reduced_electron_energy_diffusion_m2_s_m3
    )


def test_ignore_mode_does_not_apply_pt_momentum_correction(tmp_path) -> None:
    data = base_product_config(tmp_path, ["two_term"])
    data["solvers"]["two_term"]["nonconservative_model"] = "ignore"
    cfg = load_config(write_config(tmp_path, data, "ignore.yaml"))
    [case] = run(cfg, write=False).cases

    correction = case.diagnostics["two_term"][
        "temporal_growth_momentum_correction"
    ]
    assert correction["applied"] is False


def test_pt_transport_keyword_rejects_non_dc_use(tmp_path) -> None:
    data = base_product_config(tmp_path, ["two_term"])
    data["physics"]["field"] = {
        "type": "time_dependent",
        "magnetic_field": {
            "enabled": False,
            "B_T": 0.0,
            "angle_EB_deg": 0.0,
        },
        "time_dependent": {
            "waveform": "sinusoidal",
            "frequency_Hz": 13.56e6,
            "amplitude_definition": "rms",
            "momentum_response": "instantaneous",
            "phase_steps": 8,
            "max_periods": 4,
            "periodic_tolerance": 0.5,
        },
    }
    dc_data = base_product_config(tmp_path, ["two_term"])
    dc_cfg = load_config(write_config(tmp_path, dc_data, "dc.yaml"))
    [dc_case] = run(dc_cfg, write=False).cases
    rf_cfg = load_config(write_config(tmp_path, data, "rf.yaml"))
    cross_sections = load_cross_sections(rf_cfg.cross_sections, rf_cfg.conditions)
    internal = build_internal_solver_configs(rf_cfg.solvers, rf_cfg.physics)
    edges, widths = cell_edges_from_centers(dc_case.energy_eV)
    block = assemble_native_operator_block(
        rf_cfg,
        cross_sections,
        internal.two_term,
        dc_case.e_over_n_Td,
        dc_case.energy_eV,
        edges,
        widths,
    )

    with pytest.raises(ValueError, match="limited to DC two_term growth mode"):
        kinetic.transport_from_eedf(
            rf_cfg,
            block,
            dc_case.eedf,
            dc_case.e_over_n_Td,
            internal.two_term,
            temporal_growth_frequency_s_inv=1.0,
        )
