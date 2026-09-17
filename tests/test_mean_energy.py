from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from electron_swarm.core.constants import E_CHARGE_C
from swarm_workflow._io import write_csv
from swarm_workflow.comsol.input.function_eedf import (
    read_comsol_function_eedf_grid,
)
from swarm_workflow.comsol.models.gec_ccp.validation.bundle_guards import (
    GEC_ELASTIC_ENERGY_LOSS_COLUMN,
    GEC_ELASTIC_ENERGY_LOSS_TABLE,
)
from swarm_workflow.comsol.models.gec_ccp.validation.eedf import (
    _integrate_native_function_eedf_rate,
)
from swarm_workflow.comsol.models.gec_ccp.contracts import GecReactionSpec
from swarm_workflow.comsol.models.gec_ccp.mapping import load_gec_ccp_mapping
from swarm_workflow.comsol.models.gec_ccp.audits.function_eedf import (
    _audit_function_eedf_rates,
)
from swarm_workflow.comsol.models.gec_ccp.audits.transport_values import (
    _independent_elastic_energy_loss_value_audit,
)
from swarm_workflow.comsol.input.mean_energy import MeanEnergyArgument
from test_comsol_gec_ccp import _write_fixture
from test_comsol_eedf_audit import _write_eedf


def test_mean_energy_argument_preserves_interior_and_rejects_outside_support() -> None:
    argument = MeanEnergyArgument(3.0, 10.0)
    assert argument.values(np.array([3.0, 5.0, 10.0])) == pytest.approx(
        [3.0, 5.0, 10.0], rel=1e-10
    )
    with pytest.raises(ValueError, match="outside qualified support"):
        argument.values(np.array([0.3, 5.0, 100.0]))
    with pytest.raises(ValueError, match="positive and finite"):
        argument.values(np.array([0.0, np.nan]))


@pytest.mark.parametrize("limits", [(0.0, 2.0), (2.0, 2.0), (2.0, float("inf"))])
def test_mean_energy_argument_rejects_invalid_support(
    limits: tuple[float, float],
) -> None:
    with pytest.raises(ValueError, match="positive increasing finite"):
        MeanEnergyArgument(*limits)


def test_elastic_audit_rejects_values_below_declared_argument_support(
    tmp_path: Path,
) -> None:
    mapping_path, _ = _write_fixture(tmp_path)
    mapping = load_gec_ccp_mapping(mapping_path)
    mapping = replace(
        mapping,
        bundle=replace(
            mapping.bundle,
            path=tmp_path,
            expected_source="monte_carlo",
        ),
        run=replace(mapping.run, validation_mean_energy_floor_eV=3.0),
        closure=replace(
            mapping.closure,
            elastic_energy_loss_model="external_solver_native",
        ),
        results=replace(mapping.results, role="historical_validation"),
    )
    write_csv(
        tmp_path / GEC_ELASTIC_ENERGY_LOSS_TABLE,
        ("mean_energy_eV", GEC_ELASTIC_ENERGY_LOSS_COLUMN),
        [
            {
                "mean_energy_eV": energy,
                GEC_ELASTIC_ENERGY_LOSS_COLUMN: 1e-16 * energy**2,
            }
            for energy in (1.0, 2.0, 4.0, 8.0)
        ],
    )
    expression = "-ptp.ne*ptp.n_wAr*exp(sw_logKel(sw_logeps_el))*1[eV*m^3/s]"
    headers = [
        f"{name} @ t=0" for name in ("ptp.ebar", "ptp.ne", "ptp.n_wAr", expression)
    ]
    values = np.array(
        [
            [1.0, 1e14, 1e20, -1e34 * 9e-16 * E_CHARGE_C],
            [5.0, 1e14, 1e20, -1e34 * 25e-16 * E_CHARGE_C],
        ]
    )
    support = {"elastic_energy_loss_support_mean_energy_eV": [1.0, 8.0]}
    with pytest.raises(ValueError, match="outside qualified support"):
        _independent_elastic_energy_loss_value_audit(mapping, headers, values, support)


def test_preintegrated_rate_audit_uses_raw_mean_and_rejects_explicitly_unsupported_mean(
    tmp_path: Path,
) -> None:
    table = tmp_path / "eedf.txt"
    _write_eedf(table, structure="spreadsheet")
    grid = read_comsol_function_eedf_grid(table)
    argument = MeanEnergyArgument(3.0, 7.0)
    raw_means = np.array([[1.0], [5.0]])
    energies = np.array([0.0, 1.0, 2.0, 40.0])
    sigma = np.array([0.0, 0.0, 1e-20, 1e-20])
    reaction = GecReactionSpec("excitation", "eir2", "excitation")
    rates = np.array(
        [
            [
                _integrate_native_function_eedf_rate(grid, mean, energies, sigma)
                * 6.02214076e23
            ]
            for mean in raw_means.reshape(-1)
        ]
    )
    args = (
        grid,
        {"excitation": (energies, sigma)},
        (reaction,),
        raw_means,
        rates,
        {"eir2": [0]},
    )
    assert _audit_function_eedf_rates(*args)["passed"] is True
    with pytest.raises(ValueError, match="outside qualified support"):
        _audit_function_eedf_rates(*args, rate_arguments={"excitation": argument})
