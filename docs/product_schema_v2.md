# Product Schema v2

Product YAML files must set `schema_version: 2` and select solvers with
`run.solvers`.

Canonical solver ids:

- `two_term`
- `multi_term`
- `monte_carlo`

`multi_term` method values:

- `pn_closure_surrogate`: runnable product path for ordinary integral cross sections
- `pn_closure_direct`: future direct PN block operator; parsed but not runnable
- `pn_dcs`: future DCS/moment-input PN path; parsed but not runnable

Old public method names `moment_closure`, `operator`, and `hybrid` are not
schema v2 API. `allow_experimental_operator` is not a product field.
Direct PN implementation conditions are tracked in
`docs/roadmap/direct_pn_closure_operator.md`.

Angular scattering is implemented end to end for ordinary integral cross
sections with these closure models:

- `physics.angular_scattering.model`: `isotropic`, `momentum_power`, or `maxent_p1`
- `physics.angular_scattering.higher_moment_closure`: `zero`, `power`, or `maxent`

Valid pairs are fixed: `isotropic/zero`, `momentum_power/power`, and
`maxent_p1/maxent`.

Angular assumptions live only under `physics.angular_scattering`; solver-specific
angular settings are not part of the product schema.

With ordinary integral cross sections, `multi_term` records
`physics_level=surrogate`, `exact_dcs_based=false`,
`ordinary_integral_xs_closure=true`, and `direct_pn_operator=false`.

Tail-aware product metrics are controlled under `physics.energy_grid_policy`:

- `tail_metrics`: enable compact product tail metrics, default `true`
- `tail_threshold_eV`: explicit tail threshold or `null` for reaction-threshold auto selection
- `tail_rate_warning_fraction`: warning threshold for rate tail fractions

The summary CSV records tail probability, maximum reaction-rate tail fraction,
dominant tail process, high-energy cutoff rate fraction, and
`energy_grid_tail_status`. The rates CSV records per-reaction `tail_fraction`.

Electron-electron handling is schema v2 typed:

- disabled: `physics.electron_electron.enabled: false`, `model: none`
- enabled: `model: relaxation_postprocess`

`relaxation_postprocess` is a Boltzmann-family EEDF/rate postprocess. It does
not recompute transport.

Unsupported physics is controlled by `feature_policy.unsupported`:

- `fail`
- `skip_solver`
- `approximate`

`approximate` cannot ignore unsupported magnetic dynamics unless
`feature_policy.allow_unsupported_fallback: true` is explicitly set. In that
case results record `magnetic_field_treatment=ignored_fallback`.

Removed public fields:

- `run.mode`
- `both`
- `all`
- `boltzmann_two_term`
- `multiterm_boltzmann`
- `output.compatibility`
- `physics.state_resolved`
- `physics.superelastic`

State-resolved and generated-superelastic chemistry are not product schema v2
runtime features yet. Cross-section files may still contain ordinary
`superelastic` processes; the removed fields only refer to YAML-side generation
helpers.

Canonical product outputs:

- `<base>_summary.csv`
- `<base>_rates.csv`
- `<base>_eedf.csv`
- `<base>_solver_plan.csv`
- `<base>_comparison_summary.csv` when comparison is enabled

The `output` block only configures `directory`, `base_name`, and
`float_format`. Per-file write toggles and legacy alias outputs are not part of
product mode.
