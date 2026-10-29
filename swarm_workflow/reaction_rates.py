"""Materialize neutral reaction coefficients from saved product or campaign data."""

from dataclasses import asdict
import json
from pathlib import Path
import sqlite3

from electron_swarm.core.identity import content_digest
from electron_swarm.core.output import OutputConfig
from electron_swarm.io import read_product_run
from electron_swarm.io.publication import OUTPUT_SUFFIXES
from electron_swarm.io.reaction_rates import (
    RateSample,
    ReactionRatesArtifact,
    build_reaction_rates,
    rates_from_product,
    write_reaction_rates,
)
from .campaign.provenance import require_current_implementations
from .campaign.repository import read_metadata
from .campaign.store import (
    validate_workflow_schema,
    mc_sampling_plan_provenance,
    provenance_hash_manifest,
)


def workflow_rates(samples, *, solver, metadata, path, case_keys):
    """Apply the same saved-campaign coefficient policy to an already read snapshot."""
    require_current_implementations(metadata, solver)
    samples = [sample for sample in samples if sample.solver == solver]
    keys = [key for key in case_keys if key[1] == solver]
    sampling = mc_sampling_plan_provenance(metadata, required=solver == "monte_carlo")
    if solver == "monte_carlo":
        for mixture in {key[0] for key in keys}:
            for point in sampling["entries"]:
                count = sum(
                    key[0] == mixture and key[2] == point["e_over_n_Td"] for key in keys
                )
                if count != point["replicas"]:
                    raise ValueError(
                        "reaction_rates export requires the planned MC replicas at every selected point"
                    )
    source = {
        "kind": "workflow_database",
        "path": str(path),
        "scope": "saved selected-solver cases",
        "records_sha256": content_digest([asdict(sample) for sample in samples]),
        "provenance": provenance_hash_manifest(metadata, source=solver),
        "implementations": json.loads(metadata["solver_implementations_json"])[solver],
        "input_evidence": json.loads(metadata.get("input_evidence_json", "{}")),
        "execution_coverage": json.loads(metadata.get("execution_coverage_json", "{}")),
        "mc_sampling_plan": sampling if solver == "monte_carlo" else None,
    }
    return build_reaction_rates(samples, solver=solver, source=source)


def _database_rates(path: Path, solver: str) -> ReactionRatesArtifact:
    connection = sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    try:
        connection.execute("BEGIN")
        validate_workflow_schema(connection)
        metadata = read_metadata(connection)
        require_current_implementations(metadata, solver)
        cases = connection.execute(
            "SELECT mixture_id, case_id, e_over_n_Td, replicate, metadata_json FROM cases WHERE solver=? ORDER BY mixture_id, e_over_n_Td, replicate",
            (solver,),
        ).fetchall()
        samples = []
        for case in cases:
            case_metadata = json.loads(case["metadata_json"])
            context = case_metadata.pop("physical_context", {})
            case_metadata.pop("transport_observables", None)
            rates = [
                dict(row)
                for row in connection.execute(
                    "SELECT species, process, process_type, threshold_eV, target_species_fraction, rate_coefficient_m3_s FROM rates WHERE mixture_id=? AND solver=? AND e_over_n_Td=? AND replicate=? ORDER BY rate_index",
                    (
                        case["mixture_id"],
                        solver,
                        case["e_over_n_Td"],
                        case["replicate"],
                    ),
                )
            ]
            samples.append(
                RateSample(
                    solver,
                    case["case_id"],
                    case["e_over_n_Td"],
                    context,
                    rates,
                    case_metadata,
                    case["replicate"],
                    case_metadata.get("monte_carlo_case_seed"),
                )
            )
        return workflow_rates(
            samples,
            solver=solver,
            metadata=metadata,
            path=path,
            case_keys=[
                (case["mixture_id"], solver, case["e_over_n_Td"], case["replicate"])
                for case in cases
            ],
        )
    finally:
        connection.close()


def export_reaction_rates(
    source_path: str | Path, output: OutputConfig, *, solver: str
) -> dict[str, Path]:
    """Export k(E/N) points without transport, Townsend, or mean-energy-axis gates."""
    path = Path(source_path).resolve()
    output.validate()
    targets = {
        Path(output.directory).resolve() / f"{output.base_name}_{suffix}"
        for suffix in OUTPUT_SUFFIXES
    }
    if path in targets:
        raise ValueError("reaction_rates output cannot overwrite its input")
    with path.open("rb") as handle:
        database = handle.read(16) == b"SQLite format 3\0"
    if database:
        artifact = _database_rates(path, solver)
    else:
        product = read_product_run(path, verify_figures=False)
        source_files = {
            path,
            *(
                path.parent / item["name"]
                for item in product.metadata["files"].values()
            ),
        }
        if source_files & targets:
            raise ValueError(
                "reaction_rates output cannot overwrite its input generation"
            )
        artifact = rates_from_product(
            product,
            solver=solver,
            source={
                "kind": "product_manifest",
                "path": str(path),
                "manifest_content_sha256": content_digest(product.metadata),
                "evidence_scope": product.metadata["evidence_scope"],
                "execution": product.metadata.get("execution"),
                "source_execution": product.metadata.get("source_execution"),
                "implementations": product.metadata.get("implementations", {}),
                "input_evidence": product.metadata.get("input_evidence", {}),
                "input_files": product.metadata.get("input_files", {}),
            },
        )
    return write_reaction_rates(artifact, output)
