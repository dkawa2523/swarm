---
name: swarm-propagator-config
description: "Configure bounded schema v2 homogeneous-DC propagator calculations in this repository. Use for propagator grid, iteration, memory, physics, direct-run, or deterministic-sweep setup; do not launch the solver or claim P2 diffusion."
---

# Configure Propagator Calculations

Apply this skill only in this repository. Confirm the root contains `pyproject.toml`, `electron_swarm/`, `swarm_workflow/`, and `AGENTS.md`, then read the root `AGENTS.md`, `docs/product_schema_v2.md`, `docs/physics_limitations.md`, and the current Propagator qualification artifact linked there.

## Define the P1 calculation

1. Create a schema v2 product config with `run.solvers: [{id: propagator}]`. Keep field, angular scattering, ionization, e-e, and grid policy under `physics`.
2. P1 accepts homogeneous positive-E/N DC cases with `B=0`, no e-e, and no finite-k. It supports `isotropic` or `maxent_p1` angular scattering. Let `feature_policy` fail or skip unsupported requests; never remove a requested feature silently.
3. For particle-like collision kernels, distinguish `elastic_total`, `elastic_momentum_transfer`, and `effective_momentum_transfer`. Isotropic closure may use momentum transfer as total. `maxent_p1` requires explicit total and momentum-transfer inputs for every active species. An effective-only cross section cannot define collision events.
4. Use only the public `solvers.propagator` fields: `method: stationary_response`, `energy_cells`, even `polar_cells`, `max_iterations`, `convergence_tolerance`, and `max_memory_mb`.
5. Treat `max_iterations`, `physics.energy_grid_policy.max_eV_limit`, adaptive cycles, and `max_memory_mb` as hard bounds. Do not raise them simply to turn a failed point into a result.
6. P1 provides EEDF, energy-angle density, rates, growth, flux drift, mobility, and the solver-native elastic energy-loss moment. Bulk drift, particle diffusion, and electron-energy transport remain unavailable.

Use `examples/argon_propagator.yaml` as a bounded direct-run example. The GEC workflow is a target-specific pattern, not a default for another plasma model. Propagator cases use `replicate=0`; omit the workflow `mc` block when no Monte Carlo solver is enabled.

## Validate without solving

~~~powershell
py -3 -c "from electron_swarm import load_config; from electron_swarm.orchestration.plan import build_solve_plan, solver_plan_metadata; c=load_config(r'<config.yaml>'); print(solver_plan_metadata(build_solve_plan(c)))"
~~~

For a workflow:

~~~powershell
py -3 -c "from swarm_workflow.campaign import load_workflow; w=load_workflow(r'<workflow.yaml>'); print({'database': str(w.database_path), 'anchors': list(w.e_over_n_Td), 'mixtures': len(w.mixtures)})"
~~~

Estimate the configured core cell count as `energy_cells * polar_cells`; threshold insertion and append-only tail cells can increase it. Confirm the preflight memory limit is large enough, while remaining below the user's resource budget.

Finish with the config paths, E/N anchors, requested physics treatments, core cell count, cell-sweep and memory bounds, and current qualification scope. Target-specific COMSOL refinement is a separate decision made by that model adapter. Do not execute under this configuration skill.
