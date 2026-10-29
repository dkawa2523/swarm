"""Canonical product publication and verified product readers."""

from copy import deepcopy
from pathlib import Path

from electron_swarm import __version__
from electron_swarm.core.config import OutputConfig
from electron_swarm.core.provenance import publication_view
from electron_swarm.core.results import SwarmRunResult, validate_case_result
from .artifacts import PRODUCT_ARTIFACT_SCHEMA, ProductRunArtifact
from .publication import publish_tables, read_generation
from .reaction_rates import RATE_COLUMNS, build_reaction_rates, case_rate_sample
from .tables import result_tables, table_kinds, case_descriptions


def _write_generation(result: SwarmRunResult, output: OutputConfig):
    tables = result_tables(result)
    kinds = {key: table_kinds(key, frame.columns) for key, frame in tables.items()}
    reaction_rates = None
    if output.export is not None:
        reaction_rates = build_reaction_rates(
            [case_rate_sample(case) for case in result.cases],
            solver=output.export.solver,
            source={
                "kind": "product_run",
                "outcome": result.metadata["outcome"],
                "implementations": result.metadata.get("implementations", {}),
                "input_files": result.metadata.get("input_files", {}),
                "input_evidence": result.metadata.get("input_evidence", {}),
            },
        )
        tables["reaction_rates_csv"] = reaction_rates.table
        kinds["reaction_rates_csv"] = RATE_COLUMNS
    return tables, kinds, reaction_rates.metadata if reaction_rates else None


def write_outputs(result: SwarmRunResult, output: OutputConfig) -> dict[str, Path]:
    """Publish one numerical snapshot; display changes never create solver evidence."""
    output = deepcopy(output)
    output.validate()
    result = deepcopy(result)
    result, evidence_scope, execution = publication_view(result)
    for case in result.cases:
        validate_case_result(case)
    tables, kinds, reaction_rates = _write_generation(result, output)
    metadata = {
        "schema_version": 2,
        "artifact_role": "swarm_product_run",
        "artifact_schema": PRODUCT_ARTIFACT_SCHEMA,
        "package_version": result.metadata.get("package_version", __version__),
        "float_encoding": "binary64_roundtrip",
        "evidence_scope": evidence_scope,
        "execution": execution,
        "source_execution": result.metadata.get("source_execution"),
        "outcome": result.metadata["outcome"],
        "solver_plan": result.metadata.get("solver_plan", []),
        "comparison": result.metadata.get("comparison", {"status": "disabled"}),
        "resolved_solvers": result.metadata.get("resolved_solvers", {}),
        "input_evidence": result.metadata.get("input_evidence", {}),
        "cases": case_descriptions(result.cases),
        "quantity_definitions": {
            "eedf": "cell probability density per eV; integrate with energy_width_eV",
            "net_ionization_frequency_s": "net creation rate inferred from reaction rates, not a measured population-growth estimator",
            "effective_townsend_m2": "rate-derived proxy: net creation frequency divided by absolute reported flux drift and gas number density; uses retained transport when transport_stale",
            "weight_effective_count": "weight concentration, not an independent sample count or standard error",
        },
        "configuration": result.metadata.get("configuration"),
        "implementations": result.metadata.get("implementations", {}),
        "input_files": result.metadata.get("input_files", {}),
    }
    if reaction_rates:
        metadata["reaction_rates"] = reaction_rates
    return publish_tables(metadata, tables, kinds, output)


def read_product_run(
    manifest_path: str | Path, *, verify_figures: bool = True
) -> ProductRunArtifact:
    """Read verified numerical tables and, by default, verify image hashes too."""
    path = Path(manifest_path).resolve()
    suffix = "_manifest.json"
    if not path.name.endswith(suffix):
        raise ValueError("expected a product _manifest.json")
    return read_generation(
        path,
        path.name[: -len(suffix)],
        {PRODUCT_ARTIFACT_SCHEMA: "swarm_product_run"},
        verify_figures=verify_figures,
    )
