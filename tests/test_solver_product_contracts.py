from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from electron_swarm import load_config
from electron_swarm.core.results import SwarmCaseResult
from electron_swarm.core.solver_ids import CANONICAL_SOLVER_IDS
from electron_swarm.core.solver_registry import SOLVER_DESCRIPTORS
from electron_swarm.core.transport import ElectronTransport
from electron_swarm.orchestration.metadata import attach_product_metadata
from electron_swarm.orchestration.plan import SolverPlanItem, build_solve_plan
from electron_swarm.orchestration.solver_product_contracts import (
    ProductMetadataContractError,
    SOLVER_PRODUCT_CONTRACTS,
    TransportDefinition,
    TransportEvidence,
)

from product_helpers import base_product_config, write_config


TRANSPORT_BY_SOLVER = {
    "two_term": TransportDefinition.FLUX.value,
    "multi_term": TransportDefinition.PN_F1_DRIFT_F0_DIFFUSION.value,
    "monte_carlo": TransportDefinition.MC_FIXED_POPULATION.value,
    "propagator": TransportDefinition.PROPAGATOR_FLUX.value,
}


def _config_and_item(tmp_path: Path, solver: str) -> tuple[object, SolverPlanItem]:
    data = base_product_config(tmp_path, [solver])
    if solver == "monte_carlo":
        data["solvers"]["monte_carlo"] = {"seed": 17}
    config = load_config(
        write_config(
            tmp_path,
            data,
            name=f"{solver}.yaml",
        )
    )
    [item] = build_solve_plan(config)
    return config, item


def _solver_evidence(solver: str) -> dict[str, object]:
    definition = TRANSPORT_BY_SOLVER[solver]
    if solver == "two_term":
        return {}
    if solver == "multi_term":
        return {
            "solver_method": "pn_closure_direct",
            "angular_model": "isotropic",
            "direct_pn_operator": True,
            "exact_dcs_based": False,
            "ordinary_integral_xs_closure": True,
            "lmax": 3,
            "transport_definition": definition,
        }
    if solver == "monte_carlo":
        return {
            "magnetic_field_treatment": "none",
            "ionization_source_model": "equal",
            "ionization_source_treatment": "equal",
            "angular_model": "isotropic",
            "angular_moment_source": "isotropic_closure",
            "moment_table_provenance": "",
            "exact_dcs_based": False,
            "ordinary_integral_xs_closure": True,
            "tail_refinement_treatment": "disabled",
            "monte_carlo_base_seed": 17,
            "monte_carlo_case_seed": 123456,
            "transport_definition": definition,
        }
    if solver == "propagator":
        return {
            "velocity_space_representation": "axisymmetric_energy_theta_cells",
            "transport_definition": definition,
            "transport_components": "drift_velocity;mobility",
            "inelastic_angular_model": "isotropic_integral_xs_closure",
            "inelastic_radial_transfer": (
                "xs_knot_threshold_partitioned_positive_weak_projection"
            ),
            "elastic_recoil_model": (
                "momentum_transfer_driven_finite_temperature_"
                "first_mass_ratio_reversible_sg"
            ),
            "elastic_collision_xs_role": "Ar:elastic:total_and_momentum",
            "elastic_total_xs_process_ids": "Ar:elastic",
            "elastic_momentum_xs_process_ids": "Ar:elastic",
        }
    raise AssertionError(f"unexpected test solver {solver!r}")


def _case(solver: str, metadata: dict[str, object] | None = None) -> SwarmCaseResult:
    return SwarmCaseResult(
        solver=solver,
        case_id="contract_0000",
        e_over_n_Td=50.0,
        mean_energy_eV=2.0,
        net_ionization_frequency_s=0.0,
        effective_townsend_m2=0.0,
        transport=ElectronTransport(
            definition=TRANSPORT_BY_SOLVER[solver],
            gas_number_density_m3=1.0,
            drift_velocity_m_s=1.0,
            reduced_mobility_m2_V_s_m3=1.0,
            reduced_diffusion_L_m2_s_m3=None,
            reduced_diffusion_T_m2_s_m3=None,
        ),
        energy_eV=np.array([1.0]),
        eedf=np.array([1.0]),
        energy_widths_eV=np.array([1.0]),
        metadata=dict(_solver_evidence(solver) if metadata is None else metadata),
    )


def test_contract_catalog_is_typed_and_complete() -> None:
    assert set(SOLVER_PRODUCT_CONTRACTS) == {
        "two_term",
        "multi_term",
        "monte_carlo",
        "propagator",
    }
    assert (
        SOLVER_PRODUCT_CONTRACTS["two_term"].transport_evidence
        is TransportEvidence.RESULT
    )
    for solver in ("multi_term", "monte_carlo", "propagator"):
        assert (
            SOLVER_PRODUCT_CONTRACTS[solver].transport_evidence
            is TransportEvidence.RESULT_AND_SOLVER_METADATA
        )


def test_solver_catalogs_have_identical_canonical_keys_and_identities() -> None:
    canonical_ids = set(CANONICAL_SOLVER_IDS)
    assert set(SOLVER_DESCRIPTORS) == canonical_ids
    assert set(SOLVER_PRODUCT_CONTRACTS) == canonical_ids

    for solver_id in CANONICAL_SOLVER_IDS:
        descriptor = SOLVER_DESCRIPTORS[solver_id]
        contract = SOLVER_PRODUCT_CONTRACTS[solver_id]
        assert descriptor.id == solver_id
        assert descriptor.capabilities.solver == solver_id
        assert contract.solver == solver_id


@pytest.mark.parametrize(
    ("solver", "expected"),
    [
        (
            "two_term",
            {
                "solver_method": "native_sg",
                "transport_definition": "flux",
                "direct_pn_operator": False,
            },
        ),
        (
            "multi_term",
            {
                "solver_method": "pn_closure_direct",
                "transport_definition": (
                    "pn_f1_flux_drift_f0_gradient_diffusion"
                ),
                "direct_pn_operator": True,
                "lmax": 3,
                "exact_dcs_based": False,
                "ordinary_integral_xs_closure": True,
            },
        ),
        (
            "monte_carlo",
            {
                "solver_method": "internal",
                "transport_definition": (
                    "mc_flux_particle_tracking_fixed_population"
                ),
                "magnetic_field_treatment": "none",
                "ionization_source_treatment": "equal",
                "tail_refinement_treatment": "disabled",
                "monte_carlo_base_seed": 17,
                "monte_carlo_case_seed": 123456,
            },
        ),
        (
            "propagator",
            {
                "solver_method": "stationary_response",
                "transport_definition": "propagator_flux_velocity_moment",
                "velocity_space_representation": (
                    "axisymmetric_energy_theta_cells"
                ),
                "transport_components": "drift_velocity;mobility",
                "inelastic_angular_model": "isotropic_integral_xs_closure",
                "elastic_collision_xs_role": "Ar:elastic:total_and_momentum",
            },
        ),
    ],
)
def test_contracts_preserve_canonical_metadata(
    tmp_path: Path,
    solver: str,
    expected: dict[str, object],
) -> None:
    config, item = _config_and_item(tmp_path, solver)

    case = attach_product_metadata(_case(solver), config=config, item=item)

    assert {key: case.metadata[key] for key in expected} == expected
    assert case.metadata["angular_model"] == "isotropic"
    assert case.metadata["angular_moment_source"] == "isotropic_closure"
    assert case.metadata["exact_dcs_based"] is False
    assert case.metadata["ordinary_integral_xs_closure"] is True


@pytest.mark.parametrize(
    ("solver", "missing_key"),
    [
        ("multi_term", "direct_pn_operator"),
        ("monte_carlo", "magnetic_field_treatment"),
        ("propagator", "elastic_recoil_model"),
    ],
)
def test_contracts_reject_missing_solver_evidence(
    tmp_path: Path,
    solver: str,
    missing_key: str,
) -> None:
    config, item = _config_and_item(tmp_path, solver)
    metadata = _solver_evidence(solver)
    del metadata[missing_key]

    with pytest.raises(ProductMetadataContractError, match=missing_key):
        attach_product_metadata(
            _case(solver, metadata),
            config=config,
            item=item,
        )


def test_contract_rejects_transport_metadata_result_mismatch(tmp_path: Path) -> None:
    config, item = _config_and_item(tmp_path, "monte_carlo")
    metadata = _solver_evidence("monte_carlo")
    metadata["transport_definition"] = TransportDefinition.MC_WEIGHTED_GROWTH.value

    with pytest.raises(ProductMetadataContractError, match="does not match"):
        attach_product_metadata(
            _case("monte_carlo", metadata),
            config=config,
            item=item,
        )


def test_contract_rejects_missing_transport_with_domain_error(tmp_path: Path) -> None:
    config, item = _config_and_item(tmp_path, "two_term")
    case = _case("two_term")
    case.transport = None  # type: ignore[assignment]

    with pytest.raises(ProductMetadataContractError, match="required transport"):
        attach_product_metadata(case, config=config, item=item)


def test_contract_rejects_non_string_transport_definition(tmp_path: Path) -> None:
    config, item = _config_and_item(tmp_path, "two_term")
    case = _case("two_term")
    object.__setattr__(case.transport, "definition", 42)

    with pytest.raises(ProductMetadataContractError, match="must be a string"):
        attach_product_metadata(case, config=config, item=item)


def test_contract_rejects_false_exact_dcs_claim(tmp_path: Path) -> None:
    config, item = _config_and_item(tmp_path, "multi_term")
    metadata = _solver_evidence("multi_term")
    metadata["exact_dcs_based"] = True

    with pytest.raises(ProductMetadataContractError, match="exact_dcs_based"):
        attach_product_metadata(
            _case("multi_term", metadata),
            config=config,
            item=item,
        )
