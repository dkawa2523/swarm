from __future__ import annotations

from swarm_workflow.quality.policy import QualityThresholds
from swarm_workflow.tables import contracts
from swarm_workflow.tables import monte_carlo


def test_rate_gate_uses_all_basic_mc_anchors_before_profile_composition(
    monkeypatch,
) -> None:
    rows = [
        {
            "E_over_N_Td": 1.0,
            "aggregate_quality_passed": 1,
            "solver_transport_qualified": 0,
        },
        {
            "E_over_N_Td": 10.0,
            "aggregate_quality_passed": 1,
            "solver_transport_qualified": 1,
        },
    ]
    observed: dict[str, set[float]] = {}

    monkeypatch.setattr(
        monte_carlo._repository,
        "_load_quality",
        lambda *_args, **_kwargs: rows,
    )
    monkeypatch.setattr(
        monte_carlo,
        "read_metadata",
        lambda _connection: {"mc_transport_estimator_schema_version": "test"},
    )
    monkeypatch.setattr(
        monte_carlo,
        "mc_function_eedf_restricted_lmea_failure_reasons",
        lambda row, **_kwargs: ["mobility_unresolved"]
        if float(row["E_over_N_Td"]) == 1.0
        else [],
    )

    def assess_rates(*_args, allowed_e: set[float], **_kwargs):
        observed["allowed_e"] = set(allowed_e)
        return {}

    monkeypatch.setattr(
        monte_carlo,
        "_mc_eedf_rate_consistency_failures",
        assess_rates,
    )
    monkeypatch.setattr(
        monte_carlo,
        "_mc_censored_rate_relevance_failures",
        lambda *_args, **_kwargs: {},
    )

    evidence = monte_carlo._assess_monte_carlo(
        object(),
        0,
        quality_thresholds=QualityThresholds(),
        qualification_profile=(
            contracts.MC_QUALIFICATION_FUNCTION_EEDF_RESTRICTED_LMEA
        ),
    )

    assert observed["allowed_e"] == {1.0, 10.0}
    assert evidence.transport_eligible == {10.0}
