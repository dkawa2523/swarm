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
- `examples/multi_term_direct_lmax1.yaml`
- `examples/multi_term_direct_lmax4.yaml`
- `examples/compare_three_solvers.yaml`
- `examples/magnetic_mc.yaml`

- `examples/pn_dcs_moment_table.yaml`

## Solver Modes

- `two_term`: native two-term energy-space Boltzmann solver for local swarm
  transport and reaction rates.
- `multi_term`: product multi-term mode. Ordinary integral cross sections run
  through `pn_closure_direct`. The direct method is validated by the `lmax: 1`
  SG-reduction gate and supports a limited ordinary-XS angular-closure
  higher-l path. `pn_dcs` requires normalized Legendre moment-table input and
  runs the same direct PN block with the table moments as the angular source.
- `monte_carlo`: external adapter or limited internal particle backend. The
  internal backend supports DC magnetic motion with a Boris pusher. Product
  internal MC uses a fixed-particle single-daughter population model; optional
  `warmup_collisions` discards startup transients before production sampling.

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

PN vs MC comparison should use the same angular model. Set
`solvers.monte_carlo.angular_scattering: same_as_physics` to require MC output
metadata to match `physics.angular_scattering`; comparison rows include
`same_angular_model`, reference/candidate angular model fields, sampler
treatment, and mismatch reason.

## Limitations

- Ordinary-cross-section `multi_term` is an angular-closure solver, not an
  exact DCS solver. `pn_closure_direct` higher orders are limited to B=0, DC,
  axisymmetric m=0, ordinary-XS angular closure with sink-only l>0 inelastic
  treatment.
- Raw DCS angle-table parsing, RF fields, finite-k transport, full Landau e-e,
  Coulomb MC, and PN arbitrary crossed-field dynamics are not implemented.
- Validated bulk/source-gradient transport is not emitted.

## Tests

```powershell
py -3 -m pytest -q
py -3 -m pytest -q -m "not slow and not mc"
py -3 -m pytest -q -m regression
```

Product benchmark configs live under `configs/benchmarks`. Heavy Ar
BOLSIG+/MCIG comparison and MC audit commands are development-only checks; they
ingest external reference files when available and are not product solver modes.

```powershell
py -3 benchmarks\run_product_benchmarks.py --config configs\benchmarks\ar_simple.yaml
py -3 tools\benchmark_ar_eedf_consistency.py --config configs\benchmarks\ar_bolsig_eedf_consistency.yaml
py -3 tools\benchmark_ar_external_references.py --config configs\benchmarks\ar_bolsig_plus_equivalence.yaml --reference bolsig
py -3 tools\benchmark_ar_external_references.py --config configs\benchmarks\ar_mcig_reference.yaml --reference mcig
py -3 tools\benchmark_ar_bolsig_mcig_triage.py --config configs\benchmarks\ar_bolsig_mcig_triage.yaml
```

External BOLSIG+ / MCIG binaries can be run through explicit command templates,
for example `--run-bolsig "... {output} ..." --bolsig-output ...`; otherwise
the benchmarks ingest pre-generated reference CSV/text output.

See `docs/product_schema_v2.md` and `docs/testing_and_benchmarks.md` for the
schema and test marker policy. The focused Ar EEDF direct-PN gate is described
in `docs/benchmarks/ar_eedf_consistency.md`. BOLSIG+ / MCIG output ingest for
external Ar references is described in
`docs/benchmarks/ar_external_references.md`; the Ar BOLSIG+ equivalence
benchmark is described in `docs/benchmarks/ar_bolsig_plus_equivalence.md`.
The MCIG reference benchmark is described in
`docs/benchmarks/ar_mcig_reference.md`. The combined BOLSIG+ / MCIG triage is
described in `docs/benchmarks/ar_bolsig_mcig_triage.md`.
