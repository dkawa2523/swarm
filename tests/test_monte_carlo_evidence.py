from __future__ import annotations

import ast
from pathlib import Path

import electron_swarm.solvers.monte_carlo.direct_transport as direct_transport
import electron_swarm.solvers.monte_carlo.weighted_transport as weighted_transport
import electron_swarm.solvers.monte_carlo.weighted_ensemble as weighted_ensemble
import electron_swarm.solvers.monte_carlo.evidence as evidence
from electron_swarm.solvers.monte_carlo.source_identity import (
    monte_carlo_source_sha256,
)


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE_NAMES = {
    "DIRECT_MC_TRANSPORT_CADENCE_TRIAL_PERIODS",
    "DIRECT_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION",
    "DIRECT_MC_TRANSPORT_LAG_PLANES",
    "DIRECT_MC_TRANSPORT_OBSERVATION_END_FRACTION",
    "DIRECT_MC_TRANSPORT_OBSERVATION_PLANES",
    "FIELD_PARITY_RESPONSE_ESTIMATOR_SCHEMA_VERSION",
    "MC_CASE_SEED_DERIVATION",
    "MC_EEDF_ESTIMATOR_SCHEMA_VERSION",
    "MC_MAGNETIC_EEDF_ESTIMATOR_SCHEMA_VERSION",
    "MC_SEED_DERIVATION_SCHEMA_VERSION",
    "MC_WORKFLOW_REPLICA_SEED_DERIVATION",
    "WEIGHTED_MC_TAIL_ESTIMATOR_SCHEMA_VERSION",
    "WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION",
    "WEIGHTED_MC_TRANSPORT_MIN_BLOCKS",
    "WEIGHTED_MC_TRANSPORT_MIN_EFFECTIVE_LINEAGE_COUNT",
    "WEIGHTED_MC_TRANSPORT_MIN_EFFECTIVE_LINEAGE_FRACTION",
}
PUBLIC_CONSUMERS = (
    "electron_swarm/solvers/monte_carlo/direct_transport.py",
    "electron_swarm/solvers/monte_carlo/weighted_transport.py",
    "electron_swarm/solvers/monte_carlo/case_result.py",
    "swarm_workflow/comsol/models/gec_ccp/mapping_sections.py",
    "swarm_workflow/comsol/models/gec_ccp/validation/mc_quality.py",
    "swarm_workflow/comsol/models/gec_ccp/validation/provenance.py",
    "swarm_workflow/quality/monte_carlo/contracts.py",
    "swarm_workflow/quality/monte_carlo/direct_transport.py",
    "swarm_workflow/quality/monte_carlo/weighted_transport.py",
    "swarm_workflow/tables/monte_carlo.py",
    "swarm_workflow/tables/repository.py",
    "swarm_workflow/campaign/sweep.py",
)
LEGACY_OWNERS = {
    "electron_swarm.solvers.monte_carlo.transport_moments",
    "electron_swarm.solvers.monte_carlo.direct_transport",
    "electron_swarm.solvers.monte_carlo.weighted_transport",
    "electron_swarm.solvers.monte_carlo.weighted_ensemble",
}


def _legacy_evidence_imports(tree: ast.Module) -> set[str]:
    return {
        alias.name
        for node in tree.body
        if isinstance(node, ast.ImportFrom) and node.module in LEGACY_OWNERS
        for alias in node.names
        if alias.name in EVIDENCE_NAMES
    }


def _top_level_evidence_bindings(tree: ast.Module) -> set[str]:
    bound: set[str] = set()
    for node in tree.body:
        if isinstance(node, ast.Assign):
            bound.update(
                target.id
                for target in node.targets
                if isinstance(target, ast.Name) and target.id in EVIDENCE_NAMES
            )
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            if node.target.id in EVIDENCE_NAMES:
                bound.add(node.target.id)
        elif isinstance(node, ast.ImportFrom):
            bound.update(
                alias.asname or alias.name
                for alias in node.names
                if (alias.asname or alias.name) in EVIDENCE_NAMES
            )
    return bound


def test_public_monte_carlo_evidence_contract_is_complete() -> None:
    assert set(evidence.__all__) == EVIDENCE_NAMES
    assert evidence.DIRECT_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION == (
        "direct_mc_transport.v4"
    )
    assert evidence.DIRECT_MC_TRANSPORT_OBSERVATION_PLANES == 128
    assert evidence.DIRECT_MC_TRANSPORT_LAG_PLANES == (8, 16, 32, 64)
    assert evidence.DIRECT_MC_TRANSPORT_CADENCE_TRIAL_PERIODS == 16
    assert evidence.DIRECT_MC_TRANSPORT_OBSERVATION_END_FRACTION == 0.8
    assert evidence.FIELD_PARITY_RESPONSE_ESTIMATOR_SCHEMA_VERSION == (
        "common_random_number_field_parity_response.v1"
    )
    assert evidence.MC_CASE_SEED_DERIVATION == "sha256_base_seed_e_over_n.v1"
    assert evidence.MC_WORKFLOW_REPLICA_SEED_DERIVATION == (
        "sha256_base_seed_mixture_e_over_n_replicate_solver.v1"
    )
    assert evidence.MC_SEED_DERIVATION_SCHEMA_VERSION == (
        "sha256_base_seed_mixture_e_over_n_replicate_solver.v1+"
        "sha256_base_seed_e_over_n.v1"
    )
    assert evidence.MC_EEDF_ESTIMATOR_SCHEMA_VERSION == (
        "exact_bin_residence_piecewise_linear_rate_integrals.v2"
    )
    assert evidence.MC_MAGNETIC_EEDF_ESTIMATOR_SCHEMA_VERSION == (
        "boris_endpoint_energy_residence.v1"
    )
    assert evidence.WEIGHTED_MC_TRANSPORT_ESTIMATOR_SCHEMA_VERSION == (
        "direct_mc_weighted_growth_configured_lag.v6"
    )
    assert evidence.WEIGHTED_MC_TRANSPORT_MIN_BLOCKS == 128
    assert evidence.WEIGHTED_MC_TRANSPORT_MIN_EFFECTIVE_LINEAGE_COUNT == 128.0
    assert evidence.WEIGHTED_MC_TRANSPORT_MIN_EFFECTIVE_LINEAGE_FRACTION == 0.05
    assert evidence.WEIGHTED_MC_TAIL_ESTIMATOR_SCHEMA_VERSION == (
        "reaction_kernel_weighted_ensemble_pooled_estimators.v4"
    )


def test_monte_carlo_source_identity_is_deterministic_sha256() -> None:
    first = monte_carlo_source_sha256()

    assert first == monte_carlo_source_sha256()
    assert len(first) == 64
    assert set(first) <= set("0123456789abcdef")


def test_private_modules_do_not_reexport_evidence_contract() -> None:
    for module in (direct_transport, weighted_transport, weighted_ensemble):
        assert EVIDENCE_NAMES.isdisjoint(vars(module))


def test_public_module_is_the_only_product_evidence_owner() -> None:
    owners: dict[str, set[str]] = {}
    for package in ("electron_swarm", "swarm_workflow"):
        for path in (ROOT / package).rglob("*.py"):
            bound = _top_level_evidence_bindings(
                ast.parse(path.read_text(encoding="utf-8"))
            )
            if bound:
                owners[path.relative_to(ROOT).as_posix()] = bound

    assert owners == {"electron_swarm/solvers/monte_carlo/evidence.py": EVIDENCE_NAMES}


def test_evidence_consumers_depend_on_public_owner() -> None:
    for relative_path in PUBLIC_CONSUMERS:
        tree = ast.parse((ROOT / relative_path).read_text(encoding="utf-8"))
        public_imported = any(
            alias.name == "electron_swarm.solvers.monte_carlo.evidence"
            for node in tree.body
            if isinstance(node, ast.Import)
            for alias in node.names
        )
        legacy_imports = _legacy_evidence_imports(tree)
        assert public_imported, relative_path
        assert not legacy_imports, (relative_path, legacy_imports)


def test_product_code_does_not_import_evidence_from_private_modules() -> None:
    offenders: dict[str, set[str]] = {}
    for package in ("electron_swarm", "swarm_workflow"):
        for path in (ROOT / package).rglob("*.py"):
            imported = _legacy_evidence_imports(
                ast.parse(path.read_text(encoding="utf-8"))
            )
            if imported:
                offenders[str(path.relative_to(ROOT))] = imported

    assert offenders == {}
