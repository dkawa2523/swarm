"""Plot saved Swarm data without running, migrating or requalifying a solver."""

from copy import deepcopy
import json
from pathlib import Path

from electron_swarm.core.identity import content_digest
from electron_swarm.core.output import OutputConfig
from electron_swarm.core.results import SwarmRunResult
from electron_swarm.io.artifacts import (
    FIGURE_ARTIFACT_SCHEMA,
    PRODUCT_ARTIFACT_SCHEMA,
    ProductRunArtifact,
)
from electron_swarm.io.publication import (
    OUTPUT_SUFFIXES,
    publish_tables,
    read_generation,
)
from electron_swarm.io.reaction_rates import (
    RATE_ARTIFACT_SCHEMA,
    RATE_COLUMNS,
    case_rate_sample,
)
from electron_swarm.io.tables import result_tables, table_kinds, case_descriptions
from .campaign.identity import WORKFLOW_IMPLEMENTATION_KEY
from .campaign.store import WorkflowStore, WorkflowSchemaError
from .reaction_rates import workflow_rates


def _database_snapshot(path, output):
    with WorkflowStore(path, readonly=True) as store:
        original = store.metadata()
        keys = sorted(store.existing_case_keys())
        if (
            len(keys)
            != store.connection.execute("SELECT COUNT(*) FROM cases").fetchone()[0]
        ):
            raise WorkflowSchemaError("figure source contains incomplete cases")
        if not keys:
            raise ValueError("figure source contains no saved cases")
        cases = [
            store.read_case(
                mixture_id=mix, solver=solver, e_over_n_Td=field, replicate=rep
            )
            for mix, solver, field, rep in keys
        ]
        tables = result_tables(SwarmRunResult(cases=cases))
        kinds = {
            name: table_kinds(name, frame.columns) for name, frame in tables.items()
        }
        descriptions = case_descriptions(cases)
        for description, key in zip(descriptions, keys, strict=True):
            description.update(mixture_id=key[0], replicate=key[3])
        metadata = {
            "artifact_schema": FIGURE_ARTIFACT_SCHEMA,
            "artifact_role": "swarm_figures",
            "evidence_scope": "saved_results_analysis",
            "source": {
                "kind": "workflow_database",
                "path": str(path),
                "scope": "all saved cases; no new solver execution or input qualification",
                "records_sha256": content_digest(
                    {"keys": keys, "cases": cases, "provenance": original}
                ),
                "workflow_implementation": json.loads(
                    original.get(WORKFLOW_IMPLEMENTATION_KEY, "null")
                ),
                "implementations": json.loads(
                    original.get("solver_implementations_json", "{}")
                ),
                "calculation_sha256": original.get("calculation_sha256"),
                "workflow_config_sha256": original.get("workflow_config_sha256"),
                "execution_coverage": json.loads(
                    original.get("execution_coverage_json", "{}")
                ),
            },
            "cases": descriptions,
            "input_evidence": json.loads(original.get("input_evidence_json", "{}")),
            "solver_plan": json.loads(original.get("solver_plan_json", "{}")),
        }
        if output.export is not None:
            artifact = workflow_rates(
                [
                    case_rate_sample(
                        case,
                        replicate=key[3],
                        seed=case.metadata.get("monte_carlo_case_seed"),
                    )
                    for case, key in zip(cases, keys, strict=True)
                ],
                solver=output.export.solver,
                metadata=original,
                path=path,
                case_keys=keys,
            )
            metadata["reaction_rates"] = artifact.metadata
            tables["reaction_rates_csv"] = artifact.table
            kinds["reaction_rates_csv"] = RATE_COLUMNS
    return ProductRunArtifact(metadata, tables), kinds


def plot_results(source_path: str | Path, output: OutputConfig) -> dict[str, Path]:
    """Publish PNG/SVG and their native numerical snapshot from a manifest or DB.

    Replotting the same manifest in place is supported, including recovery of
    deleted figures. Numerical CSV hashes are always checked. Original solver
    execution identities remain source evidence and are never refreshed here.
    """
    output = deepcopy(output)
    output.validate()
    if not output.figures:
        raise ValueError("plot_results requires an explicit output.figures selection")
    path = Path(source_path).resolve()
    with path.open("rb") as handle:
        database = handle.read(16) == b"SQLite format 3\0"
    if database:
        targets = {
            Path(output.directory).resolve() / f"{output.base_name}_{suffix}"
            for suffix in OUTPUT_SUFFIXES
        }
        if path in targets:
            raise ValueError("figures cannot overwrite the source database")
        artifact, kinds = _database_snapshot(path, output)
    else:
        if output.export is not None:
            raise ValueError(
                "saved-manifest plotting preserves its existing export; configure export when running or use export-reaction-rates"
            )
        if not path.name.endswith("_manifest.json"):
            raise ValueError(
                "figure source must be a saved _manifest.json or workflow database"
            )
        artifact = read_generation(
            path,
            path.name[: -len("_manifest.json")],
            {
                PRODUCT_ARTIFACT_SCHEMA: "swarm_product_run",
                RATE_ARTIFACT_SCHEMA: "reaction_rates",
                FIGURE_ARTIFACT_SCHEMA: "swarm_figures",
            },
            verify_figures=False,
        )
        kinds = {
            name: {
                column["name"]: column["kind"] for column in item["schema"]["columns"]
            }
            for name, item in artifact.metadata["files"].items()
        }
        artifact.metadata.setdefault("publication", {})["source"] = {
            "kind": "saved_manifest",
            "path": str(path),
            "manifest_content_sha256": content_digest(artifact.metadata),
        }
    return publish_tables(artifact.metadata, artifact.tables, kinds, output)
