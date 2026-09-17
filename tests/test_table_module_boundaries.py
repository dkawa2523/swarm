from __future__ import annotations

import ast
from pathlib import Path
from types import ModuleType

import swarm_workflow.tables as public_tables
import swarm_workflow.tables.builder as composition
import swarm_workflow.tables.contracts as contracts
import swarm_workflow.tables.energy_loss as energy_loss
import swarm_workflow.tables.monte_carlo as mc_evidence
import swarm_workflow.tables.repository as repository


MOVED_TABLE_NAMES = (
    "TableBuildError",
    "TableBuildSummary",
    "_aggregate_rows_are_current",
    "_assess_monte_carlo",
    "_build_qualified_monte_carlo",
    "_field_source_policy",
    "_load_cases",
    "_load_eedf",
    "_load_elastic_energy_loss",
    "_load_quality",
    "_load_rate_evidence",
    "_load_rates",
    "_mc_eedf_rate_consistency_failures",
    "_mixture_ids",
    "_mixture_rows",
    "_validate_mc_censored_rate_relevance",
    "_validate_mc_sampling_plan_against_cases",
)


def _relative_imports(module: ModuleType) -> set[str]:
    tree = ast.parse(Path(module.__file__).read_text(encoding="utf-8"))
    return {
        str(node.module)
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.level
    }


def test_tables_does_not_reexport_moved_implementation() -> None:
    assert public_tables.build_tables is composition.build_tables
    for name in MOVED_TABLE_NAMES:
        assert name not in public_tables.__dict__


def test_table_dependencies_point_toward_owners() -> None:
    assert not _relative_imports(contracts)
    assert _relative_imports(repository).isdisjoint(
        {"builder", "energy_loss", "monte_carlo"}
    )
    assert _relative_imports(energy_loss).isdisjoint({"builder", "monte_carlo"})
    assert _relative_imports(mc_evidence).isdisjoint({"builder", "energy_loss"})


def test_table_modules_remain_bounded_cohesive_units() -> None:
    for module in (contracts, energy_loss, mc_evidence, repository, composition):
        physical_lines = len(
            Path(module.__file__).read_text(encoding="utf-8").splitlines()
        )
        assert physical_lines <= 1_200, module.__name__
