# Physics Limitations

This product is a solver comparison platform, not a complete plasma-chemistry
or differential-scattering package.

## Solver Meaning

- `two_term` is a local-field two-term Boltzmann solver. It reports flux
  transport from the solved EEDF.
- Ordinary-cross-section `multi_term` is a direct PN angular-closure solver.
  Integral cross sections do not determine full differential scattering, so it
  must not be read as an exact DCS multi-term solver.
- `pn_dcs` consumes normalized Legendre moment tables. It is marked
  DCS-based only when the table provenance is `dcs_derived`.
- Internal `monte_carlo` is a lightweight particle backend. Its transport is
  flux-like particle tracking and not validated bulk/source-gradient transport.

## Angular Scattering

Supported angular models are `isotropic`, `momentum_power`, `maxent_p1`, and
`moment_table`.

`momentum_power` and `maxent_p1` are closure assumptions based on ordinary
cross sections. Internal MC can sample `isotropic` directly. It can sample
`maxent_p1` only when each active scattering species has both total/effective
elastic and momentum-transfer cross sections. `momentum_power` and
`moment_table` do not define unique internal-MC samplers in this release.

## Monte Carlo Estimators

Internal MC uses a null-collision event clock and midpoint time-residence
histogramming for EEDF output. Magnetic motion is integrated with a Boris push
between trial events. Collision acceptance is evaluated at the trial event, not
at every orbit substep.

Fixed-particle ionization tracks one daughter and records the observable as
fixed-population flux-like transport. `weighted_branching` tracks both
daughters with weights and resampling, but it is still not a validated bulk
transport estimator.

## Unsupported Physics

- RF/time-dependent fields
- finite-k hydrodynamic transport
- full Landau electron-electron collisions
- Coulomb particle-particle MC
- raw angle-resolved DCS parsing
- arbitrary crossed-field PN dynamics
- YAML-side state-resolved or generated-superelastic chemistry

Requested unsupported physics is handled by `feature_policy`; it is not silently
ignored.

## Electron-Electron Treatment

`relaxation_postprocess` and `fp_energy` update Boltzmann-family EEDF and rates
after the solver run. Transport is not recomputed and is marked stale.
