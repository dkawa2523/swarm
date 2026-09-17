---
name: swarm-multi-term-config
description: "Configure schema-v2 direct PN multi-term electron-swarm calculations in this repository. Use for method, angular-moment source, lmax, direct cases, or deterministic sweep setup; do not run the solver, claim ordinary integral cross sections are DCS, or configure COMSOL export."
---

# Configure Multi-Term Calculations

Apply this skill only in this repository. Read the root `AGENTS.md`, `docs/product_schema_v2.md`, `docs/physics_limitations.md`, and `docs/theory/multiterm_closure.md`.

## Choose the physical angular model

Create schema-v2 product YAML with `run.solvers: [{id: multi_term}]`, physics requests under `physics`, and an explicit `feature_policy`. Public `solvers.multi_term` fields are only `method` and `lmax`, with `lmax >= 1`.

The current direct PN scope is homogeneous axisymmetric DC with `B=0`; RF and magnetic requests are unsupported, and e-e postprocessing leaves transport stale. Preserve those requests in the solve plan and apply `feature_policy` rather than deleting them to make a run proceed.

Choose one consistent pair:

- `method: pn_closure_direct` uses the selected ordinary-integral-cross-section angular closure (`isotropic`, `momentum_power`, or `maxent_p1`). It is a direct full-PN solve, but not an exact DCS solver.
- `method: pn_dcs` requires `physics.angular_scattering.model: moment_table`. The table contains normalized Legendre moments and is DCS-derived only when its recorded provenance is `dcs_derived`; `model_derived` and `unknown` must remain labeled accordingly.

For a moment table, require strictly increasing nonnegative energy, `m0 = 1`, finite `m_l` in `[-1, 1]`, enough moments for the requested `lmax`, explicit provenance, and `extrapolation: error`. Do not infer angular moments from a momentum-transfer cross section and call them DCS.

Use `examples/multi_term_direct_lmax4.yaml` for ordinary-XS closure and `examples/pn_dcs_moment_table.yaml` for the moment-table path. A workflow may sweep E/N or mixtures, but the current COMSOL `build-tables` source contract does not include `multi_term`.

## Validate without solving

~~~powershell
py -3 -c "from electron_swarm import load_config; from electron_swarm.orchestration.plan import build_solve_plan, solver_plan_metadata; c=load_config(r'<config.yaml>'); print(solver_plan_metadata(build_solve_plan(c)))"
~~~

For a workflow:

~~~powershell
py -3 -c "from swarm_workflow.campaign import load_workflow; w=load_workflow(r'<workflow.yaml>'); print({'database': str(w.database_path), 'anchors': list(w.e_over_n_Td), 'mixtures': len(w.mixtures)})"
~~~

Confirm the plan records method, angular source/fidelity, closure assumption, `lmax`, ionization treatment, and every unsupported/degraded feature. Finish with config paths, case count, method, angular provenance, `lmax`, physics treatments, and validation status. Do not execute under this skill.
