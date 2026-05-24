# Public Product API

Stable public entry points:

- Python package: `electron_swarm`
- CLI module: `python -m electron_swarm <config.yaml>`
- Console script: `electron-swarm <config.yaml>`
- Python helpers: `load_config`, `run`, `run_from_config`, `SwarmConfig`

Stable schema surface:

- required `schema_version: 2`
- solver ids: `two_term`, `multi_term`, `monte_carlo`
- solver selection under `run.solvers`
- physics requests under `physics.*`
- policy under `feature_policy`

Stable output files:

- `<base>_summary.csv`
- `<base>_rates.csv`
- `<base>_eedf.csv`
- `<base>_solver_plan.csv`
- `<base>_comparison_summary.csv` when comparison is enabled

Stable summary metadata:

- solver and physics level
- angular model and moment source
- ordinary-XS closure / exact-DCS flags
- ionization source treatment
- e-e treatment and transport-stale flag
- transport definition
- magnetic treatment
- compact tail metrics

Everything else is internal unless documented in the product schema.
