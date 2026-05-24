# Product Schema v2

Product YAML files must set `schema_version: 2` and select solvers with
`run.solvers`.

Canonical solver ids:

- `two_term`
- `multi_term`
- `monte_carlo`

## Multi-Term Methods

- `pn_closure_surrogate`: runnable ordinary-integral-XS angular-closure path.
- `pn_closure_direct`: accepted by schema v2 but runtime fail-fast until an
  independent coupled f0/f1 block PN solver is implemented.
- `pn_dcs`: runnable only with `physics.angular_scattering.model: moment_table`.

Ordinary-integral-XS `multi_term` results record
`ordinary_integral_xs_closure=true` and `exact_dcs_based=false`.

## Angular Scattering

Valid model/closure pairs:

- `isotropic` / `zero`
- `momentum_power` / `power`
- `maxent_p1` / `maxent`
- `moment_table` / `table`

Moment-table input is a wide CSV with `energy_eV,m0,m1,...`. Energy must be
nonnegative and strictly increasing, `m0` must be 1, and higher moments must be
finite values in `[-1, 1]`. No extrapolation is performed.

```yaml
physics:
  angular_scattering:
    model: moment_table
    higher_moment_closure: table
    moment_table:
      path: moment_tables/argon_demo_moments.csv
      format: csv
      provenance: precomputed_moments
      extrapolation: error
```

## Monte Carlo

`monte_carlo.backend: external` delegates execution to a configured command or
Python API. `angular_scattering: same_as_physics` is a strict metadata
validation contract.

`monte_carlo.backend: internal` runs the limited product particle backend and
requires `angular_scattering: same_as_physics`. It supports `isotropic` and
`maxent_p1` angular samplers plus DC magnetic dynamics through a Boris pusher.

## Physics Features

- e-e: `none`, `relaxation_postprocess`, or `fp_energy`.
- ionization source: `equal`, `primary_secondary`, or `loss_only`.
- tail metrics: compact high-energy probability/rate diagnostics.
- finite-k: request is validated and policy-handled, but no solver emits
  validated finite-k coefficients yet.

Unsupported physics is handled by:

```yaml
feature_policy:
  unsupported: fail        # fail | skip_solver | approximate
  degraded: warn           # fail | warn | record_only
  allow_unsupported_fallback: false
```

## Outputs

Canonical product outputs:

- `<base>_summary.csv`
- `<base>_rates.csv`
- `<base>_eedf.csv`
- `<base>_solver_plan.csv`
- `<base>_comparison_summary.csv` when comparison is enabled

The `output` block configures only `directory`, `base_name`, and
`float_format`.
