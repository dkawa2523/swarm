# electron-swarm

Electron swarm solver comparison tools for low-pressure plasma modeling.

The public interface is a schema v2 YAML file. It selects canonical solver
modes, requests physics features separately, writes canonical CSV outputs, and
records the solver treatment needed to interpret each result.

## Quickstart

```powershell
py -3 -m pip install -e ".[dev]"
py -3 -m electron_swarm examples\two_term.yaml --no-write
electron-swarm examples\two_term.yaml --no-write
```

Runnable examples:

- `examples/two_term.yaml`
- `examples/multi_term_direct_lmax1.yaml`
- `examples/multi_term_direct_lmax4.yaml`
- `examples/compare_three_solvers.yaml`
- `examples/magnetic_mc.yaml`
- `examples/pn_dcs_moment_table.yaml`

## Solver Modes

- `two_term`: native two-term energy-space Boltzmann solver for local flux
  transport, EEDF, and reaction-rate convolution.
- `multi_term`: direct PN comparison mode. Ordinary integral cross sections
  use angular-closure assumptions; `pn_dcs` uses normalized Legendre moment
  tables and is DCS-based only when table provenance says `dcs_derived`.
- `monte_carlo`: lightweight internal particle backend. It supports isotropic
  scattering, `maxent_p1` when both total/effective elastic and
  momentum-transfer cross sections are supplied, DC Boris magnetic motion,
  seeded sampling, and optional weighted branching. External BOLSIG+/MCIG/MC
  references live in benchmark tools, not product YAML.

## Outputs

Canonical files:

- `<base>_summary.csv`
- `<base>_rates.csv`
- `<base>_eedf.csv`
- `<base>_solver_plan.csv`
- `<base>_comparison_summary.csv` when comparison is enabled

Summary metadata is intentionally small: solver method, angular model/source,
e-e treatment, transport definition, magnetic treatment, tail treatment, and
the minimum DCS/closure flags needed for solver comparison. Detailed
internal-MC audit tables live in development tools under `tools/`.

## Limitations

- Ordinary-XS `multi_term` is an angular-closure PN solver, not an exact DCS
  solver.
- Internal MC transport is flux-like particle tracking, not validated bulk
  transport.
- Internal MC EEDF uses midpoint time-residence histogramming.
- e-e models are EEDF/rate postprocess hooks; transport is marked stale.
- RF/time-dependent fields, finite-k transport, full Landau e-e, Coulomb MC,
  raw DCS parsing, and arbitrary crossed-field PN dynamics are not implemented.

## Tests

```powershell
py -3 -m pytest -q
py -3 -m pytest -q -m "not slow and not mc"
py -3 -m pytest -q -m regression
```

See `docs/product_schema_v2.md` and `docs/physics_limitations.md` for the
schema and interpretation rules.
