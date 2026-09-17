---
name: swarm-multi-term-run
description: "Execute and verify an already configured direct PN multi-term electron-swarm calculation or deterministic sweep in this repository. Use for solver execution, lmax convergence evidence, and result checks; do not fabricate DCS provenance or build COMSOL inputs."
---

# Run Multi-Term Calculations

Apply this skill only in this repository. Read the root `AGENTS.md` and validate the config with `swarm-multi-term-config` before spending compute.

## Execute the encoded work

For a direct config:

~~~powershell
electron-swarm <config.yaml>
~~~

If the console entry point is unavailable:

~~~powershell
py -3 -c "from electron_swarm.runner import run_from_config; r=run_from_config(r'<config.yaml>'); print(r.metadata.get('output_paths', {}))"
~~~

For a deterministic workflow:

~~~powershell
py -3 -m swarm_workflow.cli sweep <workflow.yaml>
~~~

Reuse a compatible workflow database to resume missing cases. Do not delete it, change provenance, or start a duplicate process. A changed method, `lmax`, angular table, cross section, physics request, mixture, or E/N grid needs a distinct run/database.

## Verify the numerical and product result

Require finite normalized EEDF, finite F1 drift and reported rates, successful stationary-solve diagnostics, canonical solver metadata, and the expected case count. Diffusion is the declared F0-gradient approximation; do not relabel it as a full anisotropic gradient-response result.

For `pn_closure_direct`, require `ordinary_integral_xs_closure=true` and `exact_dcs_based=false`. For `pn_dcs`, verify the exact moment-table hash/provenance and set `exact_dcs_based=true` only for `dcs_derived` input.

When convergence with angular order is requested, run separate bounded configs at intentional `lmax` values and compare the same observables without changing grid, cross sections, physics, or operating conditions. Report nonconvergence as evidence; do not smooth or post-correct results to force agreement.

Canonical direct outputs are summary, rates, EEDF, and solver-plan CSVs, plus comparison output only when explicitly enabled. Workflow output remains SQLite/aggregate evidence; `multi_term` cannot currently be passed to `build-tables --source` for COMSOL.

Finish with commands, wall time, completed/failed cases, method, `lmax`, angular provenance/fidelity, key convergence differences, output/database paths, and the first failing numerical or policy gate.
