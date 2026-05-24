# electron-swarm

Electron swarm solver comparison tools for low-pressure plasma modeling.

The product interface is a schema v2 YAML file that selects canonical solvers,
requests physics features, writes canonical CSV outputs, and records how each
solver handled the requested physics.

## Install

```powershell
py -3 -m pip install -e ".[dev]"
```

For plotting helpers:

```powershell
py -3 -m pip install -e ".[plot]"
```

## Quickstart

```powershell
py -3 -m electron_swarm examples\two_term.yaml --no-write
electron-swarm examples\two_term.yaml --no-write
```

Runnable examples:

- `examples/two_term.yaml`
- `examples/multi_term_surrogate.yaml`
- `examples/compare_three_solvers.yaml`
- `examples/magnetic_mc.yaml`
- `examples/pn_dcs_moment_table.yaml`

## Solver Modes

- `two_term`: native two-term energy-space Boltzmann solver for local swarm
  transport and reaction rates.
- `multi_term`: product multi-term mode. Ordinary integral cross sections run
  through `pn_closure_surrogate`; moment-table input can run `pn_dcs`.
  `pn_closure_direct` is accepted by schema v2 but fails fast until a true
  coupled block PN solver exists.
- `monte_carlo`: external adapter or limited internal particle backend. The
  internal backend supports DC magnetic motion with a Boris pusher.

## Feature Policy

Physics features are requested separately from solver selection. Unsupported
requests are handled by `feature_policy.unsupported`:

- `fail`
- `skip_solver`
- `approximate`

Approximate fallback must be explicit in the config. Requested physics is not
silently ignored; solver-plan rows and result metadata record the treatment.

## Outputs

Canonical files:

- `<base>_summary.csv`
- `<base>_rates.csv`
- `<base>_eedf.csv`
- `<base>_solver_plan.csv`
- `<base>_comparison_summary.csv` when comparison is enabled

Summary metadata is intentionally compact: solver method, angular model,
closure or moment source, e-e treatment, transport definition, magnetic
treatment, and tail metrics.

## Limitations

- Ordinary-cross-section `multi_term` is an angular-closure solver, not an
  exact DCS solver. `pn_closure_direct` is not runnable yet because the
  independent coupled f0/f1 block operator is not implemented.
- Raw DCS angle-table parsing, RF fields, finite-k transport, full Landau e-e,
  Coulomb MC, and PN arbitrary crossed-field dynamics are not implemented.
- Validated bulk/source-gradient transport is not emitted.

## Tests

```powershell
py -3 -m pytest -q
py -3 -m pytest -q -m "not slow and not mc"
py -3 -m pytest -q -m regression
```

Benchmark configs live under `configs/benchmarks`:

```powershell
py -3 benchmarks\run_product_benchmarks.py --config configs\benchmarks\ar_simple.yaml
```

See `docs/product_schema_v2.md` and `docs/testing_and_benchmarks.md` for the
schema and test marker policy.
