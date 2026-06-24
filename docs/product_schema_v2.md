# Product Schema v2

Product YAML files must set `schema_version: 2` and select solvers with
`run.solvers`. Each solver entry is a mapping with `id` and optional `enabled`;
string entries and labels are not part of the product schema.

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

Product `monte_carlo` always means the limited internal particle backend.
External BOLSIG+, MCIG, and external MC ingest are benchmark/reference tooling
concerns and are not accepted by product YAML.

The internal MC sampler uses and reports the angular model requested under
`physics.angular_scattering`. It supports `isotropic` and `maxent_p1` angular
samplers plus DC magnetic dynamics through a Boris pusher. Unsupported angular
samplers are handled by `feature_policy` and are never silently approximated.
`maxent_p1` requires both total/effective elastic and momentum-transfer cross
sections for each active scattering species. Magnetic trajectories are
substepped only inside a sampled null-collision trial interval; collision
acceptance still occurs once at the trial event.

Public `solvers.monte_carlo` fields are limited to `population_model`, `seed`,
`particles`, `max_collisions`, and `warmup_collisions`. The default population
model is `fixed_particle_single_daughter`; `population_model:
weighted_branching` enables a bounded weighted two-daughter ionization model.
The two daughter energies follow `physics.ionization.energy_sharing`; the
ionization threshold is always recorded as an energy loss. When the particle
ensemble grows past twice `particles`, the backend systematic-resamples back to
`particles` while preserving total weight.

Internal MC EEDF output uses midpoint time-residence histogramming. Product
metadata reports only the observable definition with `transport_definition`
plus solver-comparison interpretation fields. Population-model details such as
`swarm_population_treatment`, branching counters, and estimator diagnostics are
available only when benchmark tools request diagnostic collection; they are not
part of normal product summary metadata. Weighted branching is still a
flux-like particle-tracking estimator; validated bulk transport remains out of
scope.

## Physics Features

- e-e: `none`, `relaxation_postprocess`, or `fp_energy`.
- ionization source: `equal`, `primary_secondary`, or `loss_only`.
- tail refinement: product metadata records treatment only; compact
  high-energy probability/rate diagnostics live in `case.diagnostics`.
- finite-k: request is validated and policy-handled, but no solver emits
  validated finite-k coefficients yet.

Unsupported physics is handled by:

```yaml
feature_policy:
  unsupported: fail        # fail | skip_solver
  degraded: record         # fail | record
```

Unsupported physics is never approximated or ignored as a fallback. Degraded
support is either recorded in the solver plan or rejected with `degraded: fail`.

## Outputs

Canonical product outputs:

- `<base>_summary.csv`
- `<base>_rates.csv`
- `<base>_eedf.csv`
- `<base>_solver_plan.csv`
- `<base>_comparison_summary.csv` when comparison is enabled

The EEDF CSV always includes `energy_width_eV` plus nullable `sample_count`,
`effective_sample_count`, and `relative_standard_error` columns. Boltzmann
solvers leave the sample-quality columns empty. Internal Monte Carlo fills them
from the EEDF histogram raw counts and weighted bin ESS. The EEDF is a
normalized density in `1/eV`, so `sum(eedf * energy_width_eV) == 1` up to
floating-point tolerance.

Summary metadata includes only product-facing interpretation fields: solver
method, angular model/source, ionization/e-e/magnetic/tail/transport treatment,
and DCS/closure flags. Detailed grid, convolution, tail, e-e operator,
direct-PN, and internal-MC audit values are diagnostics or development-tool
outputs, not summary columns.

PN-vs-MC comparison rows include `angular_model_status` plus scalar relative
differences and optional `eedf_l1_error`. A required comparison can still
require known, matching angular metadata, but external reference agreement is
benchmark tooling rather than a product schema contract.

The `output` block configures only `directory`, `base_name`, and
`float_format`. Product YAML does not accept external reference inputs.
BOLSIG+ / MCIG files are benchmark-tool inputs, not product schema fields.
