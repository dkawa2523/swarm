# Public Product API

Stable public entry points:

- Python package: `electron_swarm`
- CLI module: `python -m electron_swarm <config.yaml>`
- Console script: `electron-swarm <config.yaml>`
- Python helpers: `load_config`, `run`, `run_from_config`, `SwarmConfig`

Stable schema surface:

- required `schema_version: 2`
- solver ids: `two_term`, `multi_term`, `monte_carlo`
- `multi_term` methods: `pn_closure_direct`, `pn_dcs`
- solver selection under `run.solvers`
- physics requests under `physics.*`
- policy under `feature_policy`

Stable output files:

- `<base>_summary.csv`
- `<base>_rates.csv`
- `<base>_eedf.csv`
- `<base>_solver_plan.csv`
- `<base>_comparison_summary.csv` when comparison is enabled

Comparison summary rows include same-angular validation fields for PN/MC:
`same_angular_model`, `angular_model_reference`, `angular_model_candidate`,
`angular_sampler_treatment`, and `angular_model_mismatch_reason`.

Stable summary metadata:

- solver method and physics level
- angular model, moment source, and moment-table provenance
- ordinary-XS closure / exact-DCS flags
- `multi_term` `lmax` and direct-PN flag
- ionization source treatment
- e-e treatment and transport-stale flag
- transport definition
- magnetic treatment
- compact tail metrics

Direct-PN residuals, negative-mass diagnostics, benchmark failure categories,
and solver-plan capability details are internal or benchmark outputs unless
documented in the product schema.
