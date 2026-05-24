# Physics Limitations

The current product release is a solver comparison platform, not a complete
plasma-chemistry or full differential-scattering package.

`multi_term` with ordinary integral cross sections is an angular-closure or
surrogate PN solver. It must not be interpreted as an exact DCS-based multi-term
solver unless future metadata explicitly states that DCS moments are consumed by
the solver.

The ordinary-XS `multi_term` product method is `pn_closure_surrogate`.
`pn_dcs` is runnable only with `physics.angular_scattering.model: moment_table`,
which supplies validated precomputed Legendre moments to the angular transport
closure. DCS angle-table parsing is not implemented. `pn_closure_direct` still
fails fast because the independent coupled f0/f1 block PN operator is not
implemented. Direct PN implementation conditions live in
`docs/dev/direct_pn_closure_operator.md`.

Angular scattering is implemented for `isotropic`, `momentum_power`,
`maxent_p1`, and table-based Legendre moments. `momentum_power` derives the
first Legendre moment from total and momentum-transfer cross sections, then
uses a power closure for higher moments. `maxent_p1` uses the same first moment
and computes higher Legendre moments from a maximum-entropy P1 distribution.
These ordinary-XS models are closure assumptions, not differential-cross-section
providers. Screened Coulomb, ML-prior angular models, and raw DCS angle-table
parsing are future unsupported work.

Tail metrics are product output fields for high-energy reaction-rate quality.
They report tail probability and reaction-rate tail fractions after any EEDF
postprocess has been applied. They do not prove full distribution convergence:
mean energy can be stable while high-threshold chemistry is still tail
sensitive.

Transport columns are local-swarm transport fields. `two_term` and `multi_term`
report `transport_definition=flux`; internal Monte Carlo reports
`transport_definition=mc_particle_tracking`. Validated bulk/source-gradient
transport is not emitted. Finite-k transport requests are policy-handled and do
not produce coefficients in this phase.

Magnetic-field dynamics are supported only by the internal Monte Carlo backend
in this phase, using a Boris Lorentz pusher for DC fields. `two_term` and the
current axisymmetric `m=0` `multi_term` path remain unsupported for arbitrary
E-B geometry; PN magnetic support requires full spherical-harmonic dynamics.

Current limitations:

- no independent direct PN block operator
- no raw DCS angle-table parser
- no PN arbitrary crossed E-B dynamics
- no RF/time-dependent solver
- no finite-k hydrodynamic transport
- no full Landau electron-electron operator
- no particle-particle Coulomb Monte Carlo
- no YAML-side state-resolved or generated-superelastic chemistry

Ionization source energy-sharing is implemented for `two_term` only:
`equal`, `primary_secondary`, and `loss_only` are explicit product models.
Ordinary-XS `multi_term` computes ionization rates in its surrogate closure path,
but it does not solve an energy-sharing source operator; non-default source
models are rejected or skipped through `feature_policy`.

Supported e-e models are limited. `relaxation_postprocess` updates EEDF and
rates after Boltzmann-family solves, but transport is not recomputed and is
marked stale. `fp_energy` applies a simplified f0 energy-space Fokker-Planck
relaxation step for Boltzmann-family solvers; it is not a full Landau operator,
does not implement `l>0` angular damping, and does not support Coulomb MC.
Unsupported solver/physics combinations are still handled through
`feature_policy`.

Cross-section files may contain explicit `superelastic` processes. The removed
schema fields are only the old generation helpers, which were not a complete
product feature.

Requested physics is checked by the solver plan. Unsupported features fail or
skip according to `feature_policy`; they are not silently ignored.
