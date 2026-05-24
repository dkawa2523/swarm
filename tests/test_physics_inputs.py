from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from electron_swarm import load_config
from electron_swarm.collisions.ee_fp_energy import apply_fp_energy_operator
from electron_swarm.collisions.electron_electron import mean_energy_eV
from electron_swarm.physics.angular_scattering import (
    IsotropicAngularModel,
    MaxEntP1AngularModel,
    MomentTableAngularModel,
    MomentumPowerAngularModel,
    build_angular_model,
)
from electron_swarm.solvers.internal_monte_carlo import (
    boris_push,
    magnetic_field_vector,
)

from product_helpers import base_product_config, write_config


def write_moment_table(path: Path, *, bad: str | None = None) -> Path:
    if bad == "missing_m0":
        path.write_text("energy_eV,m1\n0,0\n1,0.1\n", encoding="utf-8")
        return path
    if bad == "energy_order":
        path.write_text("energy_eV,m0,m1\n1,1,0.1\n0,1,0\n", encoding="utf-8")
        return path
    if bad == "moment_range":
        path.write_text("energy_eV,m0,m1\n0,1,0\n1,1,1.5\n", encoding="utf-8")
        return path
    if bad == "sigma_columns":
        path.write_text("energy_eV,sigma0,sigma1\n0,1,0\n1,1,0.1\n", encoding="utf-8")
        return path
    path.write_text(
        "\n".join(
            [
                "energy_eV,m0,m1,m2,m3",
                "0,1,0.0,0.0,0.0",
                "10,1,0.2,0.04,0.008",
                "100,1,0.4,0.16,0.064",
                "1000,1,0.5,0.25,0.125",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return path


def test_isotropic_angular_model_moments() -> None:
    model = IsotropicAngularModel()
    moments = model.moments(np.array([1.0, 2.0]), 3)
    assert moments.shape == (4, 2)
    assert np.allclose(moments[0], 1.0)
    assert np.allclose(moments[1:], 0.0)
    assert model.metadata() == {
        "angular_model": "isotropic",
        "angular_moment_source": "isotropic_closure",
    }


def test_fp_energy_operator_preserves_normalized_finite_eedf() -> None:
    energy = np.linspace(0.05, 20.0, 120)
    widths = np.full_like(energy, energy[1] - energy[0])
    eedf = np.exp(-energy / 1.2) * (1.0 + 0.4 * np.sin(energy))
    eedf = np.clip(eedf, 0.0, None)
    eedf /= np.sum(eedf * widths)
    before = mean_energy_eV(energy, widths, eedf)

    updated, metadata = apply_fp_energy_operator(
        energy,
        widths,
        eedf,
        relaxation_fraction=0.35,
        conserve_mean_energy=True,
    )

    assert np.sum(updated * widths) == pytest.approx(1.0, abs=1.0e-12)
    assert np.all(np.isfinite(updated))
    assert np.all(updated >= 0.0)
    assert mean_energy_eV(energy, widths, updated) == pytest.approx(before, rel=0.05)
    assert metadata["electron_electron_operator_scope"] == "f0_energy_only"


def test_magnetic_field_vector_geometry() -> None:
    assert np.allclose(magnetic_field_vector(0.2, 0.0), [0.0, 0.0, 0.2])
    assert np.allclose(magnetic_field_vector(0.2, 90.0), [0.2, 0.0, 0.0], atol=1e-14)
    assert np.allclose(magnetic_field_vector(0.2, 180.0), [0.0, 0.0, -0.2], atol=1e-14)


def test_boris_push_pure_b_conserves_speed_and_zero_b_matches_e_acceleration() -> None:
    v = np.array([1.0e5, 2.0e5, 3.0e5])
    B = np.array([0.0, 0.0, 0.01])
    E = np.zeros(3)
    pushed = boris_push(v, E, B, 1.0e-11)
    assert np.linalg.norm(pushed) == pytest.approx(np.linalg.norm(v), rel=1.0e-12)

    E = np.array([0.0, 0.0, 100.0])
    dt = 2.0e-12
    pushed = boris_push(v, E, np.zeros(3), dt)
    expected = v.copy()
    expected[2] += -1.602176634e-19 / 9.1093837015e-31 * E[2] * dt
    assert np.allclose(pushed, expected)


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
    assert model.metadata() == {
        "angular_model": "momentum_power",
        "angular_moment_source": "ordinary_integral_xs_closure",
    }


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


def test_angular_model_sampling_contracts(tmp_path: Path) -> None:
    rng = np.random.default_rng(1234)
    isotropic = IsotropicAngularModel()
    samples = isotropic.sample_mu(1.0, rng, size=1000)
    assert samples.shape == (1000,)
    assert np.all(np.isfinite(samples))
    assert np.all(samples >= -1.0)
    assert np.all(samples <= 1.0)

    maxent = MaxEntP1AngularModel()
    maxent_samples = maxent.sample_mu(
        1.0,
        rng,
        size=20000,
        sigma_total=4.0,
        sigma_momentum=2.0,
    )
    assert maxent_samples.shape == (20000,)
    assert np.all(np.isfinite(maxent_samples))
    assert np.all(maxent_samples >= -1.0)
    assert np.all(maxent_samples <= 1.0)
    assert np.mean(maxent_samples) == pytest.approx(0.5, abs=0.03)

    with pytest.raises(NotImplementedError, match="does not define a unique MC sampler"):
        MomentumPowerAngularModel().sample_mu(1.0, rng, size=1)

    table_path = write_moment_table(tmp_path / "moments.csv")
    with pytest.raises(NotImplementedError, match="do not define a unique MC sampler"):
        MomentTableAngularModel(table_path).sample_mu(1.0, rng, size=1)


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

    table_path = write_moment_table(tmp_path / "moments.csv")
    data["physics"]["angular_scattering"] = {
        "model": "moment_table",
        "moment_table": {
            "path": table_path.as_posix(),
            "format": "csv",
            "provenance": "precomputed_moments",
            "extrapolation": "error",
        },
    }
    cfg = load_config(write_config(tmp_path, data))
    assert isinstance(build_angular_model(cfg), MomentTableAngularModel)

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


def test_moment_table_model_interpolation_and_validation(tmp_path: Path) -> None:
    table_path = write_moment_table(tmp_path / "moments.csv")
    model = MomentTableAngularModel(table_path)
    moments = model.moments(np.array([5.0, 55.0]), 3)
    assert moments.shape == (4, 2)
    assert np.allclose(moments[0], 1.0)
    assert np.all(np.isfinite(moments))
    assert np.all(moments <= 1.0)
    assert np.all(moments >= -1.0)
    assert model.metadata()["angular_moment_source"] == "moment_table"

    with pytest.raises(ValueError, match="does not contain enough"):
        model.moments(np.array([5.0]), 4)
    with pytest.raises(ValueError, match="does not cover"):
        model.moments(np.array([2000.0]), 1)
    with pytest.raises(ValueError, match="requires m0"):
        MomentTableAngularModel(
            write_moment_table(tmp_path / "bad1.csv", bad="missing_m0")
        ).moments(np.array([0.5]), 1)
    with pytest.raises(ValueError, match="strictly increasing"):
        MomentTableAngularModel(
            write_moment_table(tmp_path / "bad2.csv", bad="energy_order")
        ).moments(np.array([0.5]), 1)
    with pytest.raises(ValueError, match=r"\[-1, 1\]"):
        MomentTableAngularModel(
            write_moment_table(tmp_path / "bad3.csv", bad="moment_range")
        ).moments(np.array([0.5]), 1)
    with pytest.raises(ValueError, match="requires m0"):
        MomentTableAngularModel(
            write_moment_table(tmp_path / "bad4.csv", bad="sigma_columns")
        ).moments(np.array([0.5]), 1)


def test_moment_table_schema_validation(tmp_path: Path) -> None:
    table_path = write_moment_table(tmp_path / "moments.csv")
    data = base_product_config(tmp_path)
    data["physics"]["angular_scattering"] = {
        "model": "moment_table",
        "higher_moment_closure": "table",
        "moment_table": {"path": table_path.as_posix()},
    }
    cfg = load_config(write_config(tmp_path, data))
    assert cfg.physics.angular_scattering.moment_table is not None

    data["physics"]["angular_scattering"]["moment_table"] = {}
    with pytest.raises(ValueError, match="moment_table.path"):
        load_config(write_config(tmp_path, data))

    data["physics"]["angular_scattering"]["moment_table"] = {
        "path": table_path.as_posix(),
        "format": "json",
    }
    with pytest.raises(ValueError, match="moment_table.format"):
        load_config(write_config(tmp_path, data))

    data["physics"]["angular_scattering"]["moment_table"] = {
        "path": table_path.as_posix(),
        "provenance": "unknown",
    }
    with pytest.raises(ValueError, match="moment_table.provenance"):
        load_config(write_config(tmp_path, data))

    data["physics"]["angular_scattering"]["moment_table"] = {
        "path": table_path.as_posix(),
        "extrapolation": "hold",
    }
    with pytest.raises(ValueError, match="moment_table.extrapolation"):
        load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["physics"]["angular_scattering"] = {
        "model": "dcs_table",
        "dcs_table": {
            "path": "dcs.csv",
            "energy_column": "energy_eV",
            "mu_column": "mu",
            "value_column": "dcs_m2_sr",
        },
    }
    with pytest.raises(ValueError, match="dcs_table input is not implemented"):
        load_config(write_config(tmp_path, data))

    data = base_product_config(tmp_path)
    data["physics"]["angular_scattering"] = {
        "model": "moment_table",
        "dcs_table": {"path": "dcs.csv"},
    }
    with pytest.raises(ValueError, match="dcs_table input is not implemented"):
        load_config(write_config(tmp_path, data))


def test_monte_carlo_uses_shared_angular_model_config(tmp_path: Path) -> None:
    data = base_product_config(tmp_path, ["monte_carlo"])
    data["solvers"]["monte_carlo"] = {"angular_scattering": "isotropic"}
    with pytest.raises(ValueError, match="solvers.monte_carlo.angular_scattering"):
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
