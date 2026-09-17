---
name: swarm-two-term-config
description: "Configure schema v2 two-term Boltzmann electron-swarm calculations in this repository, including direct cases and E/N or mixture sweep YAML. Use for two-term setup or changes; leave solver execution and COMSOL work to their dedicated skills."
---

# Configure Two-Term Calculations

Apply this skill only in the electron-swarm repository. Confirm the root contains pyproject.toml, electron_swarm/, swarm_workflow/, and AGENTS.md, then read the root AGENTS.md.

## Produce the configuration

1. Identify the requested gas species and fractions, masses, temperature, pressure or number density, cross-section files, E/N points, requested physics, and output location. Do not assume argon or copy GEC-CCP values into an unrelated case.
2. Create or update a product YAML with:
   - schema_version: 2
   - run.solvers: [{id: two_term}]
   - only the canonical two_term id
   - physics requests under physics, separate from solver selection
   - an explicit feature_policy
3. Keep solvers.two_term limited to the public fields implemented by the parser: backend: native_sg, nonconservative_model, and min_momentum_cross_section_m2. Do not expose internal grid or convergence knobs in product YAML.
4. For a direct calculation, put the requested E/N values in run.e_over_n_Td.
5. For an external sweep, keep one base product config and create a workflow YAML with base_config, database, e_over_n_Td, mixtures, and intentional quality thresholds. A pure two-term workflow does not need MC sampling fields.
6. Resolve paths relative to the file that owns them. Use a new descriptive output/database path for a materially different campaign.

Use examples/argon_gec_ccp_base.yaml and examples/workflow_argon_gec_ccp_two_term.yaml only as target-specific examples. Read docs/product_schema_v2.md and docs/physics_limitations.md for the maintained public contract.

## Validate without solving

Load the product config and inspect its solve plan:

~~~powershell
py -3 -c "from electron_swarm import load_config; from electron_swarm.orchestration.plan import build_solve_plan, solver_plan_metadata; c=load_config(r'<config.yaml>'); print(solver_plan_metadata(build_solve_plan(c)))"
~~~

For a workflow, resolve it without starting any cases:

~~~powershell
py -3 -c "from swarm_workflow.campaign import load_workflow; w=load_workflow(r'<workflow.yaml>'); print({'database': str(w.database_path), 'anchors': len(w.e_over_n_Td), 'mixtures': len(w.mixtures)})"
~~~

Review every requested physics feature in the solve-plan metadata. Do not silence unsupported or degraded handling by removing the request. Fix the model, choose an explicit policy, or report the limitation.

Finish with the created paths, enabled solver, case count, physics treatments, and validation result. Do not execute the calculation under this configuration skill.
