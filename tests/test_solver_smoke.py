from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from electron_swarm import load_config, run

from product_helpers import base_product_config, write_config


def test_two_term_smoke_run(tmp_path: Path) -> None:
    cfg = load_config(write_config(tmp_path, base_product_config(tmp_path, ["two_term"])))
    result = run(cfg, write=False)
    [case] = result.cases
    assert case.solver == "two_term"
    assert case.schema_version == "2"
    assert case.mean_energy_eV > 0.0
    assert case.metadata["energy_grid_tail_status"] in {"ok", "warning", "insufficient"}
    assert case.metadata["tail_probability"] >= 0.0
    assert np.all(np.isfinite(case.eedf))


def test_multi_term_surrogate_smoke_and_minimum_metadata(tmp_path: Path) -> None:
    cfg = load_config(write_config(tmp_path, base_product_config(tmp_path, ["multi_term"])))
    result = run(cfg, write=False)
    [case] = result.cases
    assert case.solver == "multi_term"
    assert case.metadata["physics_level"] == "surrogate"
    assert case.metadata["angular_model"] == "isotropic"
    assert case.metadata["angular_moment_source"] == "isotropic_closure"
    assert case.metadata["angular_closure_assumption"] == "isotropic_scattering"
    assert case.metadata["solver_method"] == "pn_closure_surrogate"
    assert case.metadata["lmax"] == 3
    assert case.metadata["exact_dcs_based"] is False
    assert case.metadata["ordinary_integral_xs_closure"] is True
    assert case.metadata["direct_pn_operator"] is False
    assert case.metadata["energy_grid_tail_status"] in {"ok", "warning", "insufficient"}
    assert case.mean_energy_eV > 0.0


@pytest.mark.parametrize(
    ("model", "closure", "assumption"),
    [
        (
            "momentum_power",
            "power",
            "m1_from_total_and_momentum_xs_power_closure",
        ),
        ("maxent_p1", "maxent", "maximum_entropy_p1_closure"),
    ],
)
def test_multi_term_uses_selected_angular_closure(
    tmp_path: Path, model: str, closure: str, assumption: str
) -> None:
    data = base_product_config(tmp_path, ["multi_term"])
    data["physics"]["angular_scattering"] = {
        "model": model,
        "higher_moment_closure": closure,
    }
    cfg = load_config(write_config(tmp_path, data))
    result = run(cfg, write=False)
    [case] = result.cases

    assert case.metadata["angular_model"] == model
    assert case.metadata["angular_moment_source"] == "ordinary_integral_xs_closure"
    assert case.metadata["angular_closure_assumption"] == assumption
    assert case.metadata["lmax"] == 3
    assert case.metadata["direct_pn_operator"] is False
    assert case.metadata["exact_dcs_based"] is False
    assert case.metadata["ordinary_integral_xs_closure"] is True


def test_electron_electron_postprocess_marks_transport_stale(tmp_path: Path) -> None:
    baseline_cfg = load_config(
        write_config(
            tmp_path,
            base_product_config(tmp_path, ["two_term"]),
            name="baseline.yaml",
        )
    )
    baseline = run(baseline_cfg, write=False).cases[0]

    data = base_product_config(tmp_path, ["two_term"])
    data["physics"]["electron_electron"] = {
        "enabled": True,
        "model": "relaxation_postprocess",
        "relaxation_fraction": 0.5,
        "conserve_mean_energy": True,
    }
    cfg = load_config(write_config(tmp_path, data, name="ee.yaml"))
    result = run(cfg, write=False)
    [case] = result.cases
    assert not np.allclose(case.eedf, baseline.eedf)
    assert case.metadata["electron_electron_treatment"] == "relaxation_postprocess"
    assert case.metadata["electron_electron_affects_eedf"] is True
    assert case.metadata["electron_electron_affects_rates"] is True
    assert case.metadata["electron_electron_affects_transport"] is False
    assert case.metadata["electron_electron_transport_stale"] is True
