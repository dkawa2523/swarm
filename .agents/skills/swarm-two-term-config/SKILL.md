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
3. Keep solvers.two_term limited to backend: native_sg and min_momentum_cross_section_m2. Ionization sharing belongs to physics.ionization.energy_sharing; loss_only selects a single-daughter source explicitly. Do not expose internal grid or convergence knobs in product YAML.
4. For a direct calculation, put the requested E/N values in run.e_over_n_Td.
5. For persistence or a mixture sweep, add campaign.database and intentional campaign.quality thresholds to the same schema-v2 calculation YAML. E/N remains under run.e_over_n_Td. Use aligned fraction arrays under conditions.gas_mixture for multiple mixtures. Do not create a separate base-config file. A pure two-term campaign does not need campaign.mc fields.
6. Resolve paths relative to the file that owns them. Use a new descriptive output/database path for a materially different campaign.

Use examples/workflow_argon_gec_ccp_two_term.yaml only as target-specific examples. Read docs/product_schema_v2.md and docs/physics_limitations.md for the maintained public contract.

## Validate without solving

Validate the complete document, including any campaign controls, without solving:

~~~powershell
py -3 -m swarm_workflow.cli run <config.yaml> --plan --json
~~~

For a workflow, resolve it without starting any cases:

~~~powershell
py -3 -c "from swarm_workflow.campaign import load_workflow; w=load_workflow(r'<workflow.yaml>'); print({'database': str(w.database_path), 'anchors': len(w.e_over_n_Td), 'mixtures': len(w.mixtures)})"
~~~

Review every requested physics feature in the solve-plan metadata. Do not silence unsupported or degraded handling by removing the request. Fix the model, choose an explicit policy, or report the limitation.

Finish with the created paths, enabled solver, case count, physics treatments, and validation result. Do not execute the calculation under this configuration skill.
