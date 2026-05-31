# AGENTS.md

This repository implements electron swarm simulation tools for low-pressure plasma modeling. The current task is a destructive product refactor toward a solver comparison platform.

## Non-negotiable product direction

- Treat `two_term`, `multi_term`, and `monte_carlo` as canonical solver modes.
- Treat e-e collisions, magnetic fields, angular scattering, state-resolved processes, superelastic processes, ionization source models, and tail refinement as physics features.
- Separate solver-mode selection from physics-feature requests.
- Do not preserve obsolete public YAML behavior if it blocks product architecture.
- Require `schema_version: 2` for the new product schema.
- Use `run.solvers`, not `run.mode`.
- Use canonical solver ids `two_term`, `multi_term`, `monte_carlo`.
- Do not keep public ids `boltzmann_two_term` or `multiterm_boltzmann` in product schema.
- Do not keep `both`, `all`, `output.compatibility`, or legacy alias output files in product mode.
- Do not silently ignore requested physics. Fail, skip, or explicitly approximate according to `feature_policy`.
- Always emit metadata describing how each solver handled each requested physics feature.

## Phase implementation policy

- Treat each future phase as a limited-scope complete implementation, not a skeleton or placeholder.
- Keep each phase scope narrow, but make in-scope behavior work end to end: schema v2 config, solver plan, execution, result metadata, canonical outputs, docs, and behavior tests.
- If a requested physical model is not implemented in the current phase, mark it clearly as unsupported or raise `NotImplementedError`; do not add ambiguous placeholder behavior.
- Delete obsolete compatibility paths, thin one-use wrappers, unused helpers, duplicate tests, excessive diagnostics, and metadata that is not useful to product users.
- Keep README short and product-facing. Move longer theory, research notes, or development history to docs.
- After implementation, run `py -3 -m pytest -q`.

## Phase completion checklist

- Schema v2 config can express the feature.
- Solver plan records the feature or unsupported status.
- Runnable in-scope behavior executes end to end.
- Results contain only the necessary metadata.
- Canonical outputs include the product-visible result.
- README or docs explain the behavior briefly.
- Behavior tests cover config, plan, execution/output, and unsupported handling.
- Old public API and legacy outputs remain removed.

## Physics honesty rules

- Ordinary integral cross sections do not determine full differential scattering.
- `multi_term` with ordinary integral cross sections is an angular-closure direct PN solver unless DCS moments are actually provided.
- Never label integral-cross-section `multi_term` as an exact DCS-based multi-term solver.
- Monte Carlo and multi-term comparisons are meaningful only when angular model assumptions are explicit.
- Existing e-e relaxation postprocess can be reused initially, but it updates EEDF/rates and does not recompute transport. Mark transport stale.
- Axisymmetric `m=0` PN cannot handle arbitrary crossed E-B fields. Full magnetic support requires full spherical harmonics `Y_lm` or MC trajectory integration.

## Architecture rules

- Prefer typed dataclasses and explicit validation.
- Create a solve plan before executing solvers.
- Keep runner orchestration thin.
- Avoid growing ad hoc `if solver == ...` chains.
- Put capability and feature-policy decisions in dedicated modules.
- Keep output schema canonical and predictable.
- Write comparison outputs explicitly rather than relying on legacy aliases.
- Add tests for every schema, metadata, output, and policy change.
- Remove or rewrite tests whose only purpose is to preserve obsolete compatibility behavior.

## Expected test categories

- v2 config accepted and v1 config rejected.
- Old ids and old mode schema rejected with a migration error.
- Canonical solver list validation.
- Capability matrix correctness.
- Feature policy behavior: fail, skip, approximate/fallback.
- Solve plan generation.
- Metadata propagation to every `SwarmCaseResult`.
- Writer output names and schemas.
- Comparison outputs.
- Basic solver smoke tests.
- e-e postprocess metadata and transport-stale flag.
- Magnetic-field request metadata and strict failure behavior.

## Suggested first command after changes

```bash
py -3 -m pytest -q
```

If existing tests fail because they preserve removed compatibility behavior, update or delete them according to the new product schema rather than restoring obsolete behavior.
