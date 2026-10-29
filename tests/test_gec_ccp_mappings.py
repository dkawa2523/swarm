from dataclasses import asdict, replace
from pathlib import Path
import pytest

from swarm_workflow.comsol.models.gec_ccp import contracts
from swarm_workflow.comsol.models.gec_ccp import mapping as mapping_parser
from swarm_workflow.comsol.models.gec_ccp import mph as mph_contract
from swarm_workflow.comsol.models.gec_ccp import closure as closure_contract

PRODUCTION_MAPS = (
    "argon_gec_ccp_two_term_function_eedf.yaml",
    "argon_gec_ccp_monte_carlo_function_eedf_restricted_lmea.yaml",
    "argon_gec_ccp_propagator_function_eedf.yaml",
)

@pytest.mark.parametrize("name", PRODUCTION_MAPS)
def test_production_mapping_contract_remains_stable(name: str) -> None:
    root = Path(__file__).resolve().parents[1]
    path = root / "comsol_modes" / "maps" / name

    mapping = mapping_parser.load_gec_ccp_mapping(path)
    payload = asdict(mapping)
    contract = mph_contract.inspect_gec_ccp_mph(mapping.model.input_mph)
    mph_contract._validate_contract(mapping, contract)
    transport_contract = closure_contract._transport_input_contract(
        mapping.closure,
        source=mapping.bundle.expected_source,
    )

    assert type(mapping) is contracts.GecCcpMapping
    assert type(contract) is contracts.MphContract
    assert transport_contract["muN"]["source"] == "external_swarm_table"
    assert set(payload) == {
        "path",
        "root",
        "model",
        "bundle",
        "reactions",
        "closure",
        "run",
        "results",
        "logs",
    }


def test_mapping_error_keeps_contract_exception_identity(tmp_path: Path) -> None:
    path = tmp_path / "mapping.yaml"
    path.write_text("schema_version: 1\n", encoding="utf-8")

    with pytest.raises(
        contracts.GecCcpWorkflowError,
        match=(
            "GEC CCP mapping requires schema_version: 2; migrate "
            "run.powers_W to the single direct-solve setting run.power_W"
        ),
    ):
        mapping_parser.load_gec_ccp_mapping(path)


def test_invalid_mph_contract_keeps_validation_error() -> None:
    root = Path(__file__).resolve().parents[1]
    mapping = mapping_parser.load_gec_ccp_mapping(
        root / "comsol_modes" / "maps" / "argon_gec_ccp_two_term_function_eedf.yaml"
    )
    contract = mph_contract.inspect_gec_ccp_mph(mapping.model.input_mph)

    with pytest.raises(
        contracts.GecCcpWorkflowError,
        match=("local MPH does not match the GEC CCP mapping: physics tag=invalid"),
    ):
        mph_contract._validate_contract(
            mapping,
            replace(contract, physics_tag="invalid"),
        )
