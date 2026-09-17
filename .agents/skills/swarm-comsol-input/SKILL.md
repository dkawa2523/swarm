---
name: swarm-comsol-input
description: "Build and verify model-independent, provenance-bound COMSOL input bundles from qualified repository Swarm workflow results. Use for table materialization, source selection, bundle export, and manifest checks; do not edit a COMSOL model mapping, implement a model adapter, or run COMSOL."
---

# Build COMSOL Input Bundles

Apply this skill only in this repository. Read the root `AGENTS.md`, the target model's maintained contract, and `docs/comsol_model_adapters.md`. Inputs are a completed workflow SQLite database, one supported source solver, the qualification evidence required by that source, and a new output directory.

## Keep the bundle model-independent

`swarm_workflow/comsol/input/` owns canonical table and bundle construction. It must not import a COMSOL model adapter or execution runtime. Model tags, studies, geometry, boundary conditions, and result probes belong to the adapter map, not the bundle.

The current table-export sources are exactly `two_term`, `monte_carlo`, and `propagator`. `multi_term` is a canonical solver mode but is not yet a COMSOL table source; fail explicitly rather than relabeling it or substituting another solver.

## Build solver-pure tables

Two term:

~~~powershell
py -3 -m swarm_workflow.cli build-tables <database.sqlite> --output <tables> --source two_term
~~~

Monte Carlo:

~~~powershell
py -3 -m swarm_workflow.cli build-tables <database.sqlite> --output <tables> --source monte_carlo --mc-qualification-profile <full_transport-or-function_eedf_restricted_lmea>
~~~

Choose the MC profile from the quantities the target map consumes. `function_eedf_restricted_lmea` cannot support a full-transport claim. `--allow-unqualified-mc` produces evidence-only artifacts and must not feed a production bundle.

Propagator:

~~~powershell
py -3 -m swarm_workflow.cli build-tables <database.sqlite> --output <tables> --source propagator --solver-qualification <core-qualification.json> [--target-qualification <target-qualification.json>]
~~~

The core qualification is mandatory. Add a target qualification when the model contract requires refinement evidence over its operating range. P1 leaves diffusion and electron-energy transport unavailable; a compatible map must consume only available Propagator quantities and state what COMSOL closes. The maintained GEC map is a `physical_target` only when both qualification artifacts and all downstream audits pass.

Inspect each mixture `manifest.json`. Require `status: ok`, the requested source, current aggregate provenance, physical-context and cross-section hashes, valid ranges, exact table inventory, and all source-specific evidence.

## Consume MC selection only when the target requires it

Whole-closure MC decisions are produced by `swarm-mc-run`. Require its immutable `selection.json` when the target map calls for one. Proceed only for a terminal `accept_monte_carlo` or `select_two_term` action, and export the mixture-table directory named by `selected_solver`. A nonterminal or blocked decision returns to `swarm-mc-run` and produces no physical bundle. Keep the selection artifact outside the export directory.

## Export and verify

~~~powershell
py -3 -m swarm_workflow.cli export-comsol <table-root-or-mixture> [--selection <selection.json>] --output <bundle-output>
~~~

Verify the exported `manifest.json`, every listed SHA-256, units, independent variables, ranges, mixture and physical context, source-manifest copies, and any selection record. Never blend solver sources, fill unavailable transport, add rate floors, or hand-edit exported tables; rebuild from the database.

Finish with the source database, source solver, qualification profile/artifacts, mixture, table and bundle paths, valid ranges, and manifest/hash status. Stop before mapping or COMSOL execution.
