---
name: swarm-two-term-run
description: "Execute and verify an already configured two-term Boltzmann electron-swarm direct run or sweep in this repository. Use for two-term execution and result checks; do not design MC campaigns or build/run COMSOL inputs."
---

# Run Two-Term Calculations

Apply this skill only in the electron-swarm repository. Read the root AGENTS.md and validate the supplied config or workflow with swarm-two-term-config guidance before spending compute.

## Preflight

- Confirm the only enabled solver is two_term. If other solvers are enabled but the request is two-term only, make a focused schema v2 config instead of running extra solvers.
- Inspect the solve plan and fail on unsupported requested physics according to feature_policy.
- For a workflow, resolve the database path, E/N anchors, mixtures, and expected case count before launch.
- If the same workflow and database are already running, do not start a duplicate process.

Use the installed electron-swarm and swarm-workflow commands when available. The workflow CLI can always be invoked from the repository as py -3 -m swarm_workflow.cli.

## Execute

For a direct product config:

~~~powershell
electron-swarm <config.yaml>
~~~

If the console entry point is unavailable, call the repository API:

~~~powershell
py -3 -c "from electron_swarm.runner import run_from_config; r=run_from_config(r'<config.yaml>'); print(r.metadata.get('output_paths', {}))"
~~~

For an E/N or mixture sweep:

~~~powershell
py -3 -m swarm_workflow.cli sweep <workflow.yaml>
~~~

A matching workflow database is resumable. Re-run the same immutable workflow to fill missing cases; do not delete the database or change its provenance to force a resume. Use a new database for changed physics, cross sections, solver configuration, mixtures, or anchor plans.

When elapsed time matters, measure one complete command and report its scope. Do not extrapolate a partial run without naming the assumptions.

## Verify

For direct runs, require the canonical outputs declared by output.base_name:

- *_summary.csv
- *_rates.csv
- *_eedf.csv
- *_solver_plan.csv
- *_comparison_summary.csv only when comparison is enabled

Check that every expected E/N/mixture case exists, the solver id is two_term, the solve-plan treatment metadata is present, numeric results are finite, and each EEDF integrates to one using energy_width_eV. Do not accept obsolete alias files as evidence.

For sweeps, inspect the SQLite/aggregate outputs and verify expected row counts rather than relying only on the process exit code. Report command, wall time, completed/skipped cases, output paths, and any physics or numerical limitation. Stop at qualified Swarm results; COMSOL table construction belongs to swarm-comsol-input.
