from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from electron_swarm import load_config
from electron_swarm.core.cross_sections import (
    CrossSectionProcess,
    CrossSectionSet,
    ProcessType,
    ScatteringRole,
    load_cross_sections,
    normalize_scattering_role,
)
from electron_swarm.core.scattering import (
    resolve_particle_scattering,
    transport_momentum_processes,
)
from electron_swarm.physics.electron_neutral import (
    energy_from_speed_eV,
    first_order_elastic_recoil_energy_eV,
    ionization_daughters,
    speed_from_energy_m_s,
)

from product_helpers import base_product_config, write_config


def _write_long_csv(tmp_path: Path, rows: list[str]) -> Path:
    path = tmp_path / "roles.csv"
    path.write_text(
        "species,process,type,threshold_eV,mass_amu,energy_eV,"
        "cross_section_m2\n" + "\n".join(rows) + "\n",
        encoding="utf-8",
    )
    return path


def _load_from_path(tmp_path: Path, path: Path) -> CrossSectionSet:
    data = base_product_config(tmp_path)
    data["cross_sections"]["files"] = [
        {"path": path.as_posix(), "species": "Ar", "format": "csv"}
    ]
    config = load_config(write_config(tmp_path, data))
    return load_cross_sections(config.cross_sections, config.conditions)


def test_canonical_csv_assigns_explicit_scattering_roles(tmp_path: Path) -> None:
    rows: list[str] = []
    for process, role, scale in (
        ("total", "elastic_total", 3.0),
        ("momentum", "elastic_momentum_transfer", 2.0),
        ("effective", "effective_momentum_transfer", 4.0),
    ):
        rows.extend(
            [
                f"Ar,{process},{role},0,39.948,0,{scale}e-20",
                f"Ar,{process},{role},0,39.948,10,{scale}e-20",
            ]
        )
    cross_sections = _load_from_path(tmp_path, _write_long_csv(tmp_path, rows))
    assert {item.process: item.scattering_role for item in cross_sections.processes} == {
        "total": ScatteringRole.ELASTIC_TOTAL,
        "momentum": ScatteringRole.ELASTIC_MOMENTUM_TRANSFER,
        "effective": ScatteringRole.EFFECTIVE_MOMENTUM_TRANSFER,
    }
    assert normalize_scattering_role(
        "ELASTIC", source_format="lxcat"
    ) == ScatteringRole.ELASTIC_MOMENTUM_TRANSFER
    assert normalize_scattering_role(
        "EFFECTIVE", source_format="bolsig"
    ) == ScatteringRole.EFFECTIVE_MOMENTUM_TRANSFER


@pytest.mark.parametrize("invalid", ["-1e-20", "nan", "inf"])
def test_invalid_long_csv_cross_section_is_not_masked_by_wide_fallback(
    tmp_path: Path,
    invalid: str,
) -> None:
    path = _write_long_csv(
        tmp_path,
        [
            "Ar,total,elastic_total,0,39.948,0,1e-20",
            f"Ar,total,elastic_total,0,39.948,10,{invalid}",
        ],
    )
    with pytest.raises(ValueError, match="Invalid cross section at row 1"):
        _load_from_path(tmp_path, path)


def test_particle_and_transport_resolvers_do_not_double_count_integrals() -> None:
    energy = np.array([0.0, 10.0])

    def process(
        name: str,
        process_type: ProcessType,
        role: ScatteringRole,
    ) -> CrossSectionProcess:
        return CrossSectionProcess(
            species="Ar",
            process=name,
            process_type=process_type,
            scattering_role=role,
            energy_eV=energy,
            cross_section_m2=np.full(2, 1.0e-20),
            mass_amu=39.948,
        )

    total = process("total", ProcessType.ELASTIC, ScatteringRole.ELASTIC_TOTAL)
    momentum = process(
        "momentum",
        ProcessType.MOMENTUM,
        ScatteringRole.ELASTIC_MOMENTUM_TRANSFER,
    )
    effective = process(
        "effective",
        ProcessType.EFFECTIVE,
        ScatteringRole.EFFECTIVE_MOMENTUM_TRANSFER,
    )
    cross_sections = CrossSectionSet([total, momentum, effective])

    [particle] = resolve_particle_scattering(
        cross_sections, angular_model="maxent_p1"
    )
    assert particle.collision_total == (total,)
    assert particle.momentum_transfer == (momentum,)
    selected, includes_inelastic, model = transport_momentum_processes(
        cross_sections, "Ar"
    )
    assert selected == (effective,)
    assert includes_inelastic is True
    assert model == "effective_momentum_includes_inelastic"

    with pytest.raises(NotImplementedError, match="do not define a particle"):
        resolve_particle_scattering(
            CrossSectionSet([effective]), angular_model="isotropic"
        )

    inactive_effective = CrossSectionProcess(
        species="O2",
        process="inactive effective",
        process_type=ProcessType.EFFECTIVE,
        scattering_role=ScatteringRole.EFFECTIVE_MOMENTUM_TRANSFER,
        energy_eV=energy,
        cross_section_m2=np.full(2, 2.0e-20),
        mass_amu=31.998,
    )
    [active] = resolve_particle_scattering(
        CrossSectionSet([total, inactive_effective]),
        angular_model="isotropic",
        active_species=["Ar"],
    )
    assert active.species == "Ar"
    with pytest.raises(ValueError, match="Active gas species.*O2"):
        resolve_particle_scattering(
            CrossSectionSet([total]),
            angular_model="isotropic",
            active_species=["Ar", "O2"],
        )


def test_shared_electron_neutral_kinematics_preserve_declared_energy() -> None:
    energy = np.array([0.0, 0.1, 5.0, 100.0])
    assert energy_from_speed_eV(speed_from_energy_m_s(energy)) == pytest.approx(
        energy
    )
    recoil = first_order_elastic_recoil_energy_eV(10.0, 0.25, 39.948)
    assert 0.0 < recoil < 10.0

    equal = ionization_daughters(30.0, 15.0, model="equal")
    assert equal.energies_eV == pytest.approx((7.5, 7.5))
    primary_secondary = ionization_daughters(
        30.0,
        15.0,
        model="primary_secondary",
        secondary_electron_energy_eV=2.0,
    )
    assert primary_secondary.energies_eV == pytest.approx((13.0, 2.0))
    loss_only = ionization_daughters(30.0, 15.0, model="loss_only")
    assert loss_only.energies_eV == pytest.approx((15.0,))
    for daughters in (equal, primary_secondary, loss_only):
        assert sum(daughters.energies_eV) == pytest.approx(15.0)


@pytest.mark.parametrize(
    "process_type",
    [ProcessType.ATTACHMENT, ProcessType.EXCITATION, ProcessType.IONIZATION],
)
def test_incident_cross_sections_are_zero_below_physical_threshold(
    process_type: ProcessType,
) -> None:
    process = CrossSectionProcess(
        species="Ar",
        process=process_type.value,
        process_type=process_type,
        threshold_eV=1.0,
        energy_eV=np.array([0.0, 2.0]),
        cross_section_m2=np.array([2.0e-20, 2.0e-20]),
    )

    assert process.sigma(np.array([0.0, 0.999, 1.0, 1.5])) == pytest.approx(
        [0.0, 0.0, 2.0e-20, 2.0e-20]
    )


def test_superelastic_threshold_is_energy_release_not_incident_support() -> None:
    process = CrossSectionProcess(
        species="Ar",
        process="superelastic",
        process_type=ProcessType.SUPERELASTIC,
        threshold_eV=1.0,
        energy_eV=np.array([0.0, 2.0]),
        cross_section_m2=np.array([2.0e-20, 2.0e-20]),
    )

    assert process.sigma(np.array([0.0, 0.999, 1.0])) == pytest.approx(
        [2.0e-20, 2.0e-20, 2.0e-20]
    )
