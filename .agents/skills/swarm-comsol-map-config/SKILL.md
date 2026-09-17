---
name: swarm-comsol-map-config
description: "Create or revise a mapping for an already implemented repository COMSOL plasma-model adapter. Use for MPH tags, Swarm-bundle bindings, closure ownership, reaction mappings, run conditions, output paths, and map validation; do not build bundles, implement a new adapter, or execute COMSOL."
---

# Configure a COMSOL Model Mapping

Apply this skill only in this repository. Read the root `AGENTS.md`, `docs/comsol_model_adapters.md`, the selected adapter package, and its maintained target document.

## Route to an implemented adapter

Current adapters are:

| Model adapter | Package | Command | Maintained contract |
| --- | --- | --- | --- |
| Positive column | `swarm_workflow/comsol/models/positive_column/` | `run-positive-column` | `docs/comsol_workflow.md` |
| Argon GEC-CCP | `swarm_workflow/comsol/models/gec_ccp/` | `run-gec-ccp` | `docs/comsol_gec_ccp.md` |

If the requested plasma model has no package and CLI command, this is not a map-only task. Use `swarm-comsol-model-adapter`; do not force its MPH tags into the positive-column or GEC schema.

## Establish the target contract

1. Inspect the actual source MPH and repository evidence. Record the component, physics, study, solution, dataset, feature, reaction, parameter, and result tags; do not invent tags.
2. State ownership before editing: which EEDF, transport, rates, and elastic-loss quantities come from Swarm, and which spatial PDEs, chemistry, fields, walls, circuit, and closures remain in COMSOL.
3. Select one qualified bundle source and bind the map to its declared table inventory, physical context, support range, estimator or deterministic qualification, and closure. Do not copy source-specific expectations across solvers.
4. Keep the source MPH immutable and ensure no output or generated artifact can overwrite it. Keep result and log locations explicit; the adapter may intentionally colocate generated Java with its run artifacts.
5. Keep optional reference/baseline solves off unless the user explicitly requests a comparison.

Both current adapters require `schema_version: 2`; their remaining schemas are model-specific. The current GEC Propagator map uses `swarm_mobility_einstein`, spreadsheet Function EEDF, and `external_solver_native` elastic loss and declares `results.role: physical_target`; it also requires the maintained core and GEC target qualification artifacts. Do not restore the obsolete `diagnostic_control` label.

For MC, require a solver-selection artifact only when the adapter schema calls for it. Map each reaction to the exact bundle process/channel and MPH feature; an individual excitation or ionization channel is not a total source.

## Validate without solving

Positive-column mapping:

~~~powershell
py -3 -c "from swarm_workflow.comsol.models.positive_column.config import load_comsol_mapping; print(load_comsol_mapping(r'<mapping.yaml>', validate_files=True))"
~~~

GEC-CCP mapping:

~~~powershell
py -3 -c "from swarm_workflow.comsol.models.gec_ccp.mapping import load_gec_ccp_mapping; print(load_gec_ccp_mapping(r'<mapping.yaml>', bundle_path=r'<bundle>'))"
~~~

Run the matching command with `--dry-run` under `swarm-comsol-run` when a bundle is available. GEC preparation includes its offline MPH inspection. Positive-column preparation validates mapping/file/table structure, quality and ranges, then generates Java; it does not replace independent bundle-hash verification or licensed property/formulation readback. A dry-run is never a physical solve.

Finish with the adapter, map path, source solver, ownership split, reaction bindings, operating conditions, output/log paths, reference setting, and validation result.
