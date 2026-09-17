---
name: swarm-comsol-model-adapter
description: "Implement or refactor repository support for a distinct COMSOL plasma model. Use when a target has no existing model adapter, or when separating model-specific mapping, Java generation, execution, and acceptance from shared COMSOL input/runtime code; do not use for ordinary map edits or runs of an existing adapter."
---

# Implement a COMSOL Plasma-Model Adapter

Apply this skill only in this repository. Read the root `AGENTS.md`, `docs/comsol_model_adapters.md`, `docs/dev/architecture.md`, and the closest existing adapter without copying its model-specific assumptions.

## Define a complete, narrow adapter

Before editing, write the target contract: source MPH, geometry/dimension, plasma formulation, component/physics/study/solution/dataset tags, operating parameters, reaction inventory, imported Swarm quantities, COMSOL-owned quantities, expected exports, and numerical/physics acceptance checks. Inspect the actual MPH or existing authoritative model description; never invent tags.

Create one snake-case package under `swarm_workflow/comsol/models/<model_id>/`; expose it as a kebab-case `run-<model-id>` CLI command. Keep model-specific code there:

- typed mapping and plan/result contracts;
- schema-v2 map loading and file validation;
- MPH inspection and contract readback;
- bundle-to-model closure validation;
- generated Java for apply, solve, readback, and export;
- thin prepare and execute orchestration;
- post-solve audits and target-specific acceptance;
- optional comparison/plot code only when it consumes accepted artifacts.

Expose only the stable prepare/execute/format/error surface from the package `__init__.py`. Do not add compatibility wrappers or re-export private helpers.

## Reuse only model-independent layers

- `swarm_workflow/comsol/input/` owns model-independent, provenance-bound bundles and Function-EEDF primitives. It must not import a model package or runtime.
- `swarm_workflow/comsol/runtime/` owns executable discovery, Java compile/batch commands, process execution, logs, and artifact provenance. It must not import a model package.
- The adapter may depend on both shared layers. Dependency direction never reverses.
- Keep COMSOL tags, closure ownership, continuation strategy, solver sequence, probes, and acceptance thresholds out of shared runtime.

Do not create a universal map schema or a central solver-style `if model == ...` dispatcher. Register a thin `run-<model-id>` command in `swarm_workflow/commands/parser.py` and its handler in `handlers.py`; the handler calls the adapter's public prepare/execute API.

If the canonical bundle lacks a required observable or source—currently including `multi_term` export—extend and qualify the table/input product separately. Never fill a missing quantity from another solver or reconstruct it in the adapter merely to satisfy COMSOL.

## Deliver the adapter end to end

An in-scope adapter is complete only when it has:

1. a schema-v2 example map with an immutable source MPH, a different output MPH, and explicit non-colliding artifact/log paths;
2. a dry-run that validates offline-checkable inputs, records bindings that need licensed readback, and materializes deterministic Java without compiling or launching COMSOL;
3. a licensed execution path through the shared runtime;
4. a canonical terminal status or summary with mapping, bundle, MPH, Java/class, COMSOL-build, stage, and acceptance provenance;
5. target documentation that states the physics ownership and exact CLI sequence;
6. tests for map validation, wrong bundle/source/closure rejection, generated Java, dry-run, mocked stage execution, status/acceptance, and import boundaries.

Use the positive-column adapter as the smaller apply/verify/run/export example and GEC-CCP as the richer Function-EEDF, source-qualified, multi-stage example. Reuse concepts, not hard-coded Ar, GEC, 1D, RF, LEA, reaction, dataset, or study assumptions.

Run focused adapter and COMSOL boundary tests, Ruff on changed Python, `git diff --check`, then the repository test suite required by `AGENTS.md`. Finish with the new package, CLI command, map/docs, reused shared services, model-specific ownership, test evidence, and any explicitly unsupported closure.
