# Physics Limitations

The current product release is a solver comparison platform, not a complete
plasma-chemistry or full differential-scattering package.

`multi_term` with ordinary integral cross sections is an angular-closure direct
PN solver. It must not be interpreted as an exact DCS-based multi-term
solver unless future metadata explicitly states that DCS moments are consumed by
the solver.

The ordinary-XS `multi_term` product method is `pn_closure_direct`.
`pn_closure_direct lmax: 1` is the strict SG-reduction
gate against the native two-term reference. Higher `lmax` values are available
only in the B=0, DC, axisymmetric m=0, ordinary-XS angular-closure scope. The
current higher-l inelastic model is sink-only and does not generate
anisotropic inelastic or ionization-secondary sources. `pn_dcs` requires
`physics.angular_scattering.model: moment_table` and runs the direct PN block
with table-provided normalized Legendre moments. Raw DCS angle-table
parsing is not implemented. Direct PN implementation conditions live in
`docs/dev/direct_pn_closure_operator.md`. The explicit higher-order direct PN
scope lives in `docs/dev/direct_pn_lmax_gt1_gate.md`.

Angular scattering is implemented for `isotropic`, `momentum_power`,
`maxent_p1`, and table-based Legendre moments. `momentum_power` derives the
first Legendre moment from total and momentum-transfer cross sections, then
uses a power closure for higher moments. `maxent_p1` uses the same first moment
and computes higher Legendre moments from a maximum-entropy P1 distribution.
These ordinary-XS models are closure assumptions, not differential-cross-section
providers. Screened Coulomb, ML-prior angular models, and raw DCS angle-table
parsing are future unsupported work.

Monte Carlo can validate `same_as_physics` angular metadata and the internal MC
backend has product samplers for `isotropic` and `maxent_p1`. `momentum_power`
and `moment_table` do not define unique MC samplers in this release; strict
same-angular MC requests for those models fail or skip according to
`feature_policy` rather than using fake sampling.

The internal MC backend defaults to a fixed-particle single-daughter swarm
model. In ionization it tracks one daughter according to the configured
energy-sharing model and records the untracked secondary energy as a branching
gap. Full branching or growth-population MC is roadmap work, not a product
runnable mode. MCIG and other Monte Carlo tools should therefore be treated as
benchmark targets, not as automatic truth or fitting targets. Internal MC output
includes energy-balance and tail uncertainty metadata so model limitations and
statistically weak EEDF tail bins remain visible.

Tail metrics are product output fields for high-energy reaction-rate quality.
They report tail probability and reaction-rate tail fractions after any EEDF
postprocess has been applied. They do not prove full distribution convergence:
mean energy can be stable while high-threshold chemistry is still tail
sensitive. Internal Monte Carlo additionally writes EEDF bin counts, effective
bin counts for weighted runs, and relative standard errors; bins with low ESS
are diagnostic rather than product-grade tail evidence.
For routine BOLSIG+/MCIG-style EEDF comparison, use the fixed-particle
single-daughter Internal MC and inspect `mc_tail_comparison_status` before
interpreting the high-energy tail.

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

- no raw DCS angle-table parser
- no PN arbitrary crossed E-B dynamics
- no RF/time-dependent solver
- no finite-k hydrodynamic transport
- no full Landau electron-electron operator
- no particle-particle Coulomb Monte Carlo
- no YAML-side state-resolved or generated-superelastic chemistry

Ionization source energy-sharing is implemented for `two_term` only:
`equal`, `primary_secondary`, and `loss_only` are explicit product models.
Ordinary-XS `multi_term` computes ionization rates from the solved direct-PN
`F0`; non-default energy-sharing source models are rejected or skipped through
`feature_policy` until their higher-l source terms are defined.

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
