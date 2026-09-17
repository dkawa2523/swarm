"""Typed product-metadata contracts for canonical solver results.

The solve plan describes what a solver is expected to do.  These contracts
validate the evidence returned by the solver before that evidence is exposed
as canonical product metadata.  They deliberately do not provide fallback
values for solver-specific execution facts.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from types import MappingProxyType
from typing import TYPE_CHECKING, Protocol

from electron_swarm.core.result_metadata import (
    TRANSPORT_FLUX,
    TRANSPORT_PN_F1_DRIFT_F0_DIFFUSION,
    TRANSPORT_MC_FIXED_POPULATION,
    TRANSPORT_MC_WEIGHTED_GROWTH,
    TRANSPORT_PROPAGATOR_FLUX,
)

if TYPE_CHECKING:  # pragma: no cover
    from electron_swarm.core.config import SwarmConfig
    from electron_swarm.core.results import SwarmCaseResult
    from electron_swarm.orchestration.plan import SolverPlanItem


class ProductMetadataContractError(ValueError):
    """Raised when a solver result lacks or contradicts required evidence."""


class TransportDefinition(str, Enum):
    """Canonical transport definitions accepted by product contracts."""

    FLUX = TRANSPORT_FLUX
    PN_F1_DRIFT_F0_DIFFUSION = TRANSPORT_PN_F1_DRIFT_F0_DIFFUSION
    MC_FIXED_POPULATION = TRANSPORT_MC_FIXED_POPULATION
    MC_WEIGHTED_GROWTH = TRANSPORT_MC_WEIGHTED_GROWTH
    PROPAGATOR_FLUX = TRANSPORT_PROPAGATOR_FLUX


class TransportEvidence(str, Enum):
    """Where a contract requires its transport definition to be recorded."""

    RESULT = "result"
    RESULT_AND_SOLVER_METADATA = "result_and_solver_metadata"


class ProductMetadataBuilder(Protocol):
    def __call__(
        self,
        case: SwarmCaseResult,
        config: SwarmConfig,
        item: SolverPlanItem,
        selected_angular_metadata: Mapping[str, object],
    ) -> dict[str, object]: ...


@dataclass(frozen=True, slots=True)
class SolverProductContract:
    """Product-facing evidence requirements for one canonical solver."""

    solver: str
    transport_definitions: frozenset[TransportDefinition]
    transport_evidence: TransportEvidence
    metadata_builder: ProductMetadataBuilder

    def product_metadata(
        self,
        case: SwarmCaseResult,
        *,
        config: SwarmConfig,
        item: SolverPlanItem,
        selected_angular_metadata: Mapping[str, object],
    ) -> dict[str, object]:
        if case.solver != self.solver:
            raise ProductMetadataContractError(
                f"{self.solver} product contract received result from "
                f"{case.solver!r}"
            )

        transport = case.transport
        if transport is None:
            raise ProductMetadataContractError(
                f"{self.solver} result lacks required transport evidence"
            )
        definition_text = transport.definition
        if not isinstance(definition_text, str):
            raise ProductMetadataContractError(
                f"{self.solver} transport definition must be a string; "
                f"found {type(definition_text).__name__}"
            )
        try:
            definition = TransportDefinition(definition_text)
        except ValueError as exc:
            raise ProductMetadataContractError(
                f"{self.solver} returned unknown transport definition "
                f"{definition_text!r}"
            ) from exc
        if definition not in self.transport_definitions:
            allowed = ", ".join(
                sorted(value.value for value in self.transport_definitions)
            )
            raise ProductMetadataContractError(
                f"{self.solver} returned transport definition "
                f"{definition_text!r}; expected one of: {allowed}"
            )

        if self.transport_evidence is TransportEvidence.RESULT_AND_SOLVER_METADATA:
            metadata_definition = _required_metadata(case, "transport_definition")
            if type(metadata_definition) is not str or metadata_definition != definition_text:
                raise ProductMetadataContractError(
                    f"{self.solver} metadata transport_definition "
                    f"{metadata_definition!r} does not match transport result "
                    f"{definition_text!r}"
                )

        metadata = self.metadata_builder(
            case,
            config,
            item,
            selected_angular_metadata,
        )
        metadata["transport_definition"] = definition_text
        return metadata


def _required_metadata(case: SwarmCaseResult, key: str) -> object:
    if key not in case.metadata:
        raise ProductMetadataContractError(
            f"{case.solver} result lacks required product evidence {key!r}"
        )
    return case.metadata[key]


def _checked_metadata(
    case: SwarmCaseResult,
    expected: Mapping[str, object],
) -> dict[str, object]:
    checked: dict[str, object] = {}
    for key, expected_value in expected.items():
        actual = _required_metadata(case, key)
        if type(actual) is not type(expected_value) or actual != expected_value:
            raise ProductMetadataContractError(
                f"{case.solver} product evidence {key!r} is {actual!r}; "
                f"expected {expected_value!r}"
            )
        checked[key] = actual
    return checked


def _optional_checked_metadata(
    case: SwarmCaseResult,
    expected: Mapping[str, object],
) -> dict[str, object]:
    present = {key: value for key, value in expected.items() if key in case.metadata}
    return _checked_metadata(case, present)


def _required_metadata_values(
    case: SwarmCaseResult,
    keys: tuple[str, ...],
) -> dict[str, object]:
    return {key: _required_metadata(case, key) for key in keys}


def _two_term_metadata(
    case: SwarmCaseResult,
    config: SwarmConfig,
    item: SolverPlanItem,
    selected_angular_metadata: Mapping[str, object],
) -> dict[str, object]:
    del case, config, item, selected_angular_metadata
    return {}


def _multi_term_metadata(
    case: SwarmCaseResult,
    config: SwarmConfig,
    item: SolverPlanItem,
    selected_angular_metadata: Mapping[str, object],
) -> dict[str, object]:
    del item
    method = str(config.solvers.multi_term.method)
    is_pn_dcs = method == "pn_dcs"
    metadata = _checked_metadata(
        case,
        {
            "solver_method": method,
            "lmax": int(config.solvers.multi_term.lmax),
            "direct_pn_operator": True,
            "exact_dcs_based": bool(
                is_pn_dcs and selected_angular_metadata["exact_dcs_based"]
            ),
            "ordinary_integral_xs_closure": not is_pn_dcs,
        },
    )
    metadata.update(
        _optional_checked_metadata(
            case,
            {
                "angular_model": selected_angular_metadata["angular_model"],
                "angular_moment_source": selected_angular_metadata[
                    "angular_moment_source"
                ],
                "moment_table_provenance": selected_angular_metadata[
                    "moment_table_provenance"
                ],
            },
        )
    )
    return metadata


def _monte_carlo_metadata(
    case: SwarmCaseResult,
    config: SwarmConfig,
    item: SolverPlanItem,
    selected_angular_metadata: Mapping[str, object],
) -> dict[str, object]:
    expected = {
        "magnetic_field_treatment": item.treatment("magnetic_field"),
        "ionization_source_model": str(config.physics.ionization.energy_sharing),
        "ionization_source_treatment": item.treatment("ionization_source"),
        "monte_carlo_base_seed": int(config.solvers.monte_carlo.seed),
    }
    expected.update(selected_angular_metadata)
    metadata = _checked_metadata(case, expected)
    case_seed = _required_metadata(case, "monte_carlo_case_seed")
    if type(case_seed) is not int or not 0 <= case_seed <= 0xFFFFFFFF:
        raise ProductMetadataContractError(
            "monte_carlo product evidence 'monte_carlo_case_seed' must be a "
            "uint32 integer"
        )
    tail_treatment = _required_metadata(case, "tail_refinement_treatment")
    configured = config.solvers.monte_carlo.tail_max_collisions is not None
    allowed_tail_treatments = (
        {"configured_not_triggered", "executed"} if configured else {"disabled"}
    )
    if tail_treatment not in allowed_tail_treatments:
        raise ProductMetadataContractError(
            "monte_carlo product evidence 'tail_refinement_treatment' is "
            f"{tail_treatment!r}; expected one of {sorted(allowed_tail_treatments)}"
        )
    metadata["monte_carlo_case_seed"] = case_seed
    metadata["tail_refinement_treatment"] = tail_treatment
    return metadata


def _propagator_metadata(
    case: SwarmCaseResult,
    config: SwarmConfig,
    item: SolverPlanItem,
    selected_angular_metadata: Mapping[str, object],
) -> dict[str, object]:
    del config, item, selected_angular_metadata
    metadata = _checked_metadata(
        case,
        {
            "velocity_space_representation": "axisymmetric_energy_theta_cells",
            "transport_components": "drift_velocity;mobility",
            "inelastic_angular_model": "isotropic_integral_xs_closure",
            "inelastic_radial_transfer": (
                "xs_knot_threshold_partitioned_positive_weak_projection"
            ),
            "elastic_recoil_model": (
                "momentum_transfer_driven_finite_temperature_"
                "first_mass_ratio_reversible_sg"
            ),
        },
    )
    metadata.update(
        _required_metadata_values(
            case,
            (
                "elastic_collision_xs_role",
                "elastic_total_xs_process_ids",
                "elastic_momentum_xs_process_ids",
            ),
        )
    )
    return metadata


SOLVER_PRODUCT_CONTRACTS: Mapping[str, SolverProductContract] = MappingProxyType(
    {
        "two_term": SolverProductContract(
            solver="two_term",
            transport_definitions=frozenset({TransportDefinition.FLUX}),
            transport_evidence=TransportEvidence.RESULT,
            metadata_builder=_two_term_metadata,
        ),
        "multi_term": SolverProductContract(
            solver="multi_term",
            transport_definitions=frozenset(
                {TransportDefinition.PN_F1_DRIFT_F0_DIFFUSION}
            ),
            transport_evidence=TransportEvidence.RESULT_AND_SOLVER_METADATA,
            metadata_builder=_multi_term_metadata,
        ),
        "monte_carlo": SolverProductContract(
            solver="monte_carlo",
            transport_definitions=frozenset(
                {
                    TransportDefinition.MC_FIXED_POPULATION,
                    TransportDefinition.MC_WEIGHTED_GROWTH,
                }
            ),
            transport_evidence=TransportEvidence.RESULT_AND_SOLVER_METADATA,
            metadata_builder=_monte_carlo_metadata,
        ),
        "propagator": SolverProductContract(
            solver="propagator",
            transport_definitions=frozenset({TransportDefinition.PROPAGATOR_FLUX}),
            transport_evidence=TransportEvidence.RESULT_AND_SOLVER_METADATA,
            metadata_builder=_propagator_metadata,
        ),
    }
)


def solver_product_contract(solver: str) -> SolverProductContract:
    """Return the product contract for a canonical solver id."""

    try:
        return SOLVER_PRODUCT_CONTRACTS[solver]
    except KeyError as exc:
        raise ValueError(f"Unknown canonical solver id: {solver!r}") from exc


__all__ = [
    "ProductMetadataContractError",
    "SOLVER_PRODUCT_CONTRACTS",
    "SolverProductContract",
    "TransportDefinition",
    "TransportEvidence",
    "solver_product_contract",
]
