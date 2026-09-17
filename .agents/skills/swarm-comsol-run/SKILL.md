---
name: swarm-comsol-run
description: "Preflight, execute, and audit an implemented repository COMSOL plasma-model adapter using a validated Swarm bundle and mapping. Use for dry-runs, licensed solves, timing, and target-specific acceptance; do not create solver data, redesign maps, or relax failed evidence gates."
---

# Run a COMSOL Plasma-Model Adapter

Apply this skill only in this repository. Read the root `AGENTS.md`, `docs/comsol_model_adapters.md`, the selected map and bundle manifest, and the adapter's maintained contract.

## Select the adapter command

| Adapter | Dry-run and execution command | Terminal evidence |
| --- | --- | --- |
| Positive column | `run-positive-column` | `positive_column_run_summary.json` |
| Argon GEC-CCP | `run-gec-ccp` | `gec_ccp_run_status.json` |
| Argon GEC-ICP | `run-gec-icp` | `run_status.json` |

For another model, use only its implemented package, documented `run-<model>` command, and target-specific acceptance artifact. If no adapter exists, stop and use `swarm-comsol-model-adapter`; no current adapter is a generic fallback.

## Preflight

1. Confirm the output MPH differs from the source MPH and all result, generated-code, and log paths are intentional, writable, and non-colliding. Never modify the source MPH in place.
2. Confirm bundle source and hashes, physical context, support range, table inventory, closure profile, reaction inventory, solver selection, and required deterministic/statistical qualifications match the map.
3. Check for an active COMSOL process using the same output or result paths; do not launch a duplicate.
4. Resolve COMSOL in the maintained order: `--comsol`, `COMSOL_BATCH`, `COMSOL_EXECUTABLE`, then `comsol`/`comsolbatch` on `PATH`. Use `COMSOL_LICENSE_FILE` only when explicitly configured; do not commit machine-specific executable or license paths to a map.
5. Run the adapter dry-run first and inspect its materialized plan, source hashes, generated Java, offline file/model checks, and exact declared study sequence. Record which MPH bindings still require licensed COMSOL readback. Class compilation and live readback belong to licensed execution, not dry-run.

~~~powershell
py -3 -m swarm_workflow.cli run-positive-column <mapping.yaml> --bundle <bundle> --dry-run
py -3 -m swarm_workflow.cli run-gec-ccp <mapping.yaml> --bundle <bundle> --dry-run
py -3 -m swarm_workflow.cli run-gec-icp <mapping.yaml> --bundle <bundle> --dry-run
~~~

A dry-run proves preparation only, not license availability or physical convergence.

## Execute

~~~powershell
py -3 -m swarm_workflow.cli run-positive-column <mapping.yaml> --bundle <bundle> [--comsol <comsolbatch.exe>]
py -3 -m swarm_workflow.cli run-gec-ccp <mapping.yaml> --bundle <bundle> [--comsol <comsolbatch.exe>]
py -3 -m swarm_workflow.cli run-gec-icp <mapping.yaml> --bundle <bundle> [--comsol <comsolbatch.exe>]
~~~

Do not add an optional built-in/reference solve unless the user requested it and the map already declares it. On failure, inspect the recorded stage; do not blindly retry license, compile, stale-class, schema, provenance, or physics-gate failures.

## Apply target-specific acceptance

Common acceptance requires a completed licensed solve, immutable mapping/bundle/MPH/Java provenance, expected exports, finite physically allowed outputs, no unapproved support escape, and every audit required by the target contract. A zero process exit alone is insufficient.

For GEC `physical_target`, require both `quality_accepted_for_declared_closure = true` and `physical_target_accepted = true`, plus every applicable EEDF, transport, rate-censoring, conservation, support, and saved-model audit. A declared diagnostic or ablation run must retain its role and cannot be promoted by postprocessing.

For positive column, require completed apply, verification, run, and export stages; formulation and model-condition readback; all declared property and spatial consistency checks; and the expected result CSV. Keep apply, verify, solve, and export timings separate.

Report adapter, mapping and bundle hashes, COMSOL version/build, studies executed, per-stage and total runtime, output MPH, results/logs, terminal acceptance, and the first failing gate if rejected.
