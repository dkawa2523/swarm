---
name: swarm-mc-config
description: "Configure reproducible, compute-bounded internal Monte Carlo electron-swarm jobs in this repository. Use for MC physics, sampling plans, replicas, seeds, quality gates, and budget review; do not launch the campaign."
---

# Configure Monte Carlo Calculations

Apply this skill only in the electron-swarm repository. Confirm the repository markers and read the root AGENTS.md, docs/product_schema_v2.md, and docs/dev/internal_monte_carlo.md.

## Define the physical calculation

1. Create a schema v2 base config with run.solvers: [{id: monte_carlo}]. Keep physics requests under physics and use an explicit feature_policy.
2. Choose the population model from the observable:
   - fixed_particle_single_daughter for fixed-population particle tracking.
   - weighted_branching when both ionization daughters and population growth must be represented.
3. Use only public MC fields accepted by `solvers.monte_carlo`: `population_model`, `seed`, `particles`, `warmup_collisions`, `max_collisions`, `tail_max_collisions`, `tail_rate_rse_trigger`, `transport_correlation_lag_barriers`, and `numeric_kernel`. The recorded global collision clock is solver provenance, not a public YAML option.
4. Prefer numeric_kernel: auto. Use python for reference or audit work. An explicit unsupported numba combination must fail.
5. Check the current solve plan for angular scattering, magnetic field, e-e, RF, ionization, and tail handling. Never convert an unsupported sampler or field model to isotropic/DC silently.
6. Keep cross-section high-energy behavior explicit. Ordinary integral cross sections do not supply differential scattering.

## Build a bounded workflow

For qualification or a sweep, create a workflow with explicit E/N anchors, mixtures, SQLite path, quality thresholds, and an integer mc.base_seed. Independent replicas must use derived distinct seeds.

Use mc.sampling_plan for field-specific work. The loader resolves a complete plan: omitted anchors inherit the base solver values, while listed rows override that anchor. Prefer full rows when reproducibility matters. Never list an anchor absent from mc.e_over_n_Td or the workflow E/N grid.

For each effective row, calculate the nominal upper bound:

~~~text
per replica = particles * (warmup_collisions + max_collisions + tail_max_collisions)
total       = sum(per replica * replicas over anchors)
~~~

The tail term is a conditional upper bound. Configure mc.convergence with finite per-replica and total particle-barrier ceilings, finite maximum_attempts, and replica limits. A bounded convergence campaign supports one mixture. Do not add reuse_database unless the user intentionally reuses a compatible prior campaign.

For weighted-growth function_eedf_restricted_lmea qualification, require positive warmup, a power-of-two correlation lag L >= 4, and production length that is a multiple of 2L with at least 128 complete blocks. Fixed-population configs cannot set the weighted-growth lag or tail controls.

Choose mc.workers no larger than the number of independent MC jobs or the useful local CPU capacity. More workers do not change the statistical sample count.

## Validate without running

~~~powershell
py -3 -c "from swarm_workflow.campaign import load_workflow; from swarm_workflow.quality.monte_carlo.policy import sampling_budget_provenance; w=load_workflow(r'<workflow.yaml>'); print({'jobs': sum(x.replicas for x in w.mc_sampling_plan), 'workers': w.mc_workers, 'budget': sampling_budget_provenance(w.mc_sampling_plan)})"
~~~

Validation must resolve the canonical sampling plan and reject any ceiling violation before execution. Finish with config/workflow paths, physics treatment, anchors, replicas, workers, nominal budget, conditional-tail assumption, and the intended qualification profile. Do not start MC under this configuration skill.
