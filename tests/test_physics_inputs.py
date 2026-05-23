from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from electron_swarm import load_config
from electron_swarm.physics.angular_scattering import (
    IsotropicAngularModel,
    MaxEntP1AngularModel,
    MomentumPowerAngularModel,
    build_angular_model,
)

from product_helpers import base_product_config, write_config


def test_isotropic_angular_model_moments() -> None:
    model = IsotropicAngularModel()
    moments = model.moments(np.array([1.0, 2.0]), 3)
    assert moments.shape == (4, 2)
    assert np.allclose(moments[0], 1.0)
    assert np.allclose(moments[1:], 0.0)
    assert model.metadata()["angular_closure_assumption"] == "isotropic_scattering"


def test_momentum_power_angular_model_power_closure() -> None:
    energy = np.array([1.0, 2.0])
    sigma_total = np.array([4.0, 10.0])
    sigma_momentum = np.array([1.0, 12.0])
    model = MomentumPowerAngularModel()

    moments = model.moments(
        energy,
        3,
        sigma_total=sigma_total,
        sigma_momentum=sigma_momentum,
    )

    expected_m1 = np.array([0.75, -0.2])
    assert np.allclose(moments[0], 1.0)
    assert np.allclose(moments[1], expected_m1)
    assert np.allclose(moments[2], expected_m1**2)
    assert np.allclose(moments[3], expected_m1**3)
    assert model.metadata()["angular_higher_moment_closure"] == "power"


def test_maxent_p1_isotropic_limit_and_finite_shape() -> None:
    energy = np.array([1.0, 2.0])
    model = MaxEntP1AngularModel()

    isotropic = model.moments(
        energy,
        4,
        sigma_total=np.array([2.0, 2.0]),
        sigma_momentum=np.array([2.0, 2.0]),
    )
    assert isotropic.shape == (5, 2)
    assert np.allclose(isotropic[0], 1.0)
    assert np.allclose(isotropic[1:], 0.0)

    moments = model.moments(
        energy,
        4,
        sigma_total=np.array([4.0, 5.0]),
        sigma_momentum=np.array([2.0, 6.0]),
    )
    assert moments.shape == (5, 2)
    assert np.all(np.isfinite(moments))
    assert np.all(moments <= 1.0)
    assert np.all(moments >= -1.0)
    assert np.allclose(moments[0], 1.0)
    assert np.allclose(moments[1], np.array([0.5, -0.2]), atol=1.0e-10)


def test_angular_models_reject_invalid_cross_sections() -> None:
    model = MomentumPowerAngularModel()
    energy = np.array([1.0, 2.0])
    with pytest.raises(ValueError, match="sigma_total is required"):
        model.moments(energy, 2, sigma_momentum=np.ones(2))
    with pytest.raises(ValueError, match="sigma_momentum is required"):
        model.moments(energy, 2, sigma_total=np.ones(2))
    with pytest.raises(ValueError, match="sigma_total must be positive"):
        model.moments(
            energy,
            2,
            sigma_total=np.array([1.0, 0.0]),
            sigma_momentum=np.ones(2),
        )
    with pytest.raises(ValueError, match="must match energy_eV shape"):
        model.moments(
            energy,
            2,
            sigma_total=np.ones(1),
            sigma_momentum=np.ones(1),
        )


def test_build_angular_model_from_config_and_rejects_invalid_schema(
    tmp_path: Path,
) -> None:
    data = base_product_config(tmp_path)
    data["physics"]["angular_scattering"] = {
        "model": "momentum_power",
        "higher_moment_closure": "power",
    }
    cfg = load_config(write_config(tmp_path, data))
    assert isinstance(build_angular_model(cfg), MomentumPowerAngularModel)

    data["physics"]["angular_scattering"] = {
        "model": "maxent_p1",
        "higher_moment_closure": "maxent",
    }
    cfg = load_config(write_config(tmp_path, data))
    assert isinstance(build_angular_model(cfg), MaxEntP1AngularModel)

    data["physics"]["angular_scattering"] = {
        "model": "momentum_only",
        "higher_moment_closure": "power",
    }
    with pytest.raises(ValueError, match="momentum_power"):
        load_config(write_config(tmp_path, data))

    data["physics"]["angular_scattering"] = {
        "model": "maxent_p1",
        "higher_moment_closure": "power",
    }
    with pytest.raises(ValueError, match="model and higher_moment_closure mismatch"):
        load_config(write_config(tmp_path, data))

    data["physics"]["angular_scattering"] = {
        "model": "momentum_power",
        "higher_moment_closure": "power",
        "uncertainty_ensemble": False,
    }
    with pytest.raises(ValueError, match="Unsupported physics.angular_scattering"):
        load_config(write_config(tmp_path, data))


def test_monte_carlo_uses_shared_angular_model_config(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {"angular_scattering": "isotropic"}
    with pytest.raises(ValueError, match="Unsupported solvers.monte_carlo fields"):
        load_config(write_config(tmp_path, data))


def test_unimplemented_state_resolved_schema_is_rejected(tmp_path: Path) -> None:
    data = base_product_config(tmp_path)
    data["physics"]["state_resolved"] = {"enabled": True}
    with pytest.raises(ValueError, match="Unsupported physics fields"):
        load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["physics"]["superelastic"] = {"enabled": True}
    with pytest.raises(ValueError, match="Unsupported physics fields"):
        load_config(write_config(tmp_path, data))
