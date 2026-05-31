# Product Schema v2

Product YAML files must set `schema_version: 2` and select solvers with
`run.solvers`.

Canonical solver ids:

- `two_term`
- `multi_term`
- `monte_carlo`

## Multi-Term Methods

- `pn_closure_direct`: runnable for `lmax: 1` as the direct SG-reduction
  regression path against the native two-term reference, and runnable for
  higher `lmax` in the limited B=0/DC/axisymmetric ordinary-XS angular-closure
  scope.
- `pn_dcs`: runnable with `physics.angular_scattering.model: moment_table`.
  It uses normalized Legendre moments from the table as the angular source for
  the direct PN block. Ordinary XS closure alone cannot run `pn_dcs`.

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
`exact_dcs_based=true` applies only when `provenance: dcs_derived`.

```yaml
physics:
  angular_scattering:
    model: moment_table
    higher_moment_closure: table
    moment_table:
      path: moment_tables/argon_demo_moments.csv
      format: normalized_legendre_moments
      provenance: model_derived
      extrapolation: error
```

## Monte Carlo

`monte_carlo.backend: external` delegates execution to a configured command or
Python API. `angular_scattering: same_as_physics` is a strict metadata
validation contract.

`monte_carlo.backend: internal` runs the limited product particle backend and
requires `angular_scattering: same_as_physics`. It supports `isotropic` and
`maxent_p1` angular samplers plus DC magnetic dynamics through a Boris pusher.
The default internal population model is
`mc_population_model=fixed_particle_single_daughter`: ionization follows one
tracked daughter, while the untracked secondary-energy gap is exposed through MC
audit metadata. `max_collisions` controls the production sampling length.
`warmup_collisions` is optional and discards initial transient flights from the
EEDF/transport estimator without changing the public solver mode.
Full ionization branching and weighted growth-population MC remain roadmap and
are not accepted in product schema.

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

The EEDF CSV always includes nullable `sample_count`, `effective_sample_count`,
and `relative_standard_error` columns. Boltzmann solvers leave them empty;
internal Monte Carlo fills them from EEDF histogram counts or weighted bin ESS
so weak tail bins are visible instead of being mistaken for converged physics.
Internal MC summary metadata also includes `mc_tail_comparison_status`:
`ok`, `weak_tail_statistics`, `energy_balance_warning`,
or `energy_balance_fail`.
Tail bins with low ESS do not automatically fail the whole comparison when
their combined probability mass is below the product threshold; the reported
`mc_tail_weak_probability_fraction` records that residual weak-tail mass.

PN-vs-MC comparison rows include `same_angular_model`,
`angular_model_reference`, `angular_model_candidate`,
`angular_sampler_treatment`, and `angular_model_mismatch_reason`. A required
comparison fails when a PN/MC row has unknown or mismatched angular metadata.

The `output` block configures only `directory`, `base_name`, and
`float_format`.

External BOLSIG+ / MCIG benchmark files can be listed under
`references.external`. They are reference sources, not solver ids, and are used
only by benchmark tooling. The supported ingest format is
`electron_swarm_reference_csv` with either `eedf_eV_inv` or `eepf_eV_m32`; both
are converted to normalized EEDF `F(E)` in `1/eV` before comparison.
