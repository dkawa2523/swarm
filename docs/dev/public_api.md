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

Comparison summary rows include `angular_model_status`, scalar relative
differences, and optional `eedf_l1_error`.

Stable summary metadata:

- solver method and angular model/source
- ordinary-XS closure / exact-DCS flags
- `multi_term` `lmax`
- ionization source treatment
- e-e treatment and transport-stale flag
- transport definition
- magnetic treatment
- tail refinement treatment

Direct-PN residuals, negative-mass diagnostics, MC audit counters, benchmark
failure categories, and solver-plan capability details are internal or
benchmark outputs unless documented in the product schema.
