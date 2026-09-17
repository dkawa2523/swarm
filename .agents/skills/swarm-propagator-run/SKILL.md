---
name: swarm-propagator-run
description: "Execute, assess, and benchmark bounded P1 propagator calculations in this repository. Use for direct runs, deterministic sweeps, convergence evidence, grid refinement, or qualification; leave COMSOL bundle and solve work to COMSOL skills."
---

# Run Propagator Calculations

Apply this skill only in this repository. Read the root `AGENTS.md` and validate the config with `swarm-propagator-config` guidance before spending compute.

## Execute the encoded work

For a direct config:

~~~powershell
electron-swarm <config.yaml>
~~~

If the console entry point is unavailable:

~~~powershell
py -3 -c "from electron_swarm.runner import run_from_config; r=run_from_config(r'<config.yaml>'); print(r.metadata.get('output_paths', {}))"
~~~

For a deterministic E/N or mixture sweep:

~~~powershell
py -3 -m swarm_workflow.cli sweep <workflow.yaml>
~~~

A matching SQLite workflow resumes missing cases without rerunning completed points. Do not delete it, alter provenance, or start a duplicate process. Do not automatically increase the grid, cell-sweep budget, energy ceiling, or memory after a convergence failure. The exception diagnostics are evaluation evidence, not a result to persist as a passed case.

## Verify each completed case

Require:

- `converged=true`, normalized nonnegative cell population, operator residual at or below its recorded tolerance, number-balance residual at or below `1e-12`, passed tail and outer-flux gates, and memory below the configured bound;
- finite EEDF, mean energy, rates, drift, and mobility;
- energy-angle density integrating to one with `energy_width_eV * mu_width`;
- particle diffusion and electron-energy transport remaining empty/NULL/unavailable;
- product metadata naming the velocity-space representation, angular kernel, collision-XS role, recoil approximation, and transport components.

Canonical direct output includes `*_energy_angle_distribution.csv` in addition to summary, rates, EEDF, and solver-plan files. For workflows, verify `replicate=0`, nullable transport columns, `energy_angle_bins`, aggregate quality, and exact expected row counts.

## Run numerical qualification

Use the maintained bounded runner when method qualification or performance evidence is requested:

~~~powershell
py -3 tools\qualify_propagator_p1_deterministic.py
~~~

It records environment, code state, cross-section hashes, grid profiles, wall time, cell sweeps, residuals, memory, medium/fine differences, and release decisions. It is Propagator-only; cross-solver or stochastic comparisons are outside this qualification.

The deterministic P0/P1 core qualification passes. A COMSOL target may additionally require its own grid/refinement evidence. The maintained GEC map is a `physical_target` only when its core qualification, GEC target qualification, and all declared downstream audits pass; do not change a role to bypass those gates.

Finish with commands, wall time per scope, passed/failed anchors, first failing gate, output/database paths, and the qualification artifact. COMSOL input construction belongs to `swarm-comsol-input`; model-specific target qualification belongs to that adapter contract.
