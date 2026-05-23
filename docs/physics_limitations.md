# Physics Limitations

The current product release is a solver comparison platform, not a complete
plasma-chemistry or full differential-scattering package.

`multi_term` with ordinary integral cross sections is an angular-closure or
surrogate PN solver. It must not be interpreted as an exact DCS-based multi-term
solver unless future metadata explicitly states that DCS moments are consumed by
the solver.

The runnable `multi_term` product method is `pn_closure_surrogate`.
`pn_closure_direct` and `pn_dcs` are roadmap methods that parse as schema v2
method values but fail fast at execution because the direct PN block operator
and DCS moment provider are not implemented. Direct PN implementation
conditions live in `docs/roadmap/direct_pn_closure_operator.md`.

Phase 4 angular scattering is a limited-scope complete implementation for
`isotropic`, `momentum_power`, and `maxent_p1`. `momentum_power` derives the
first Legendre moment from total and momentum-transfer cross sections, then
uses a power closure for higher moments. `maxent_p1` uses the same first moment
and computes higher Legendre moments from a maximum-entropy P1 distribution.
These are closure assumptions for ordinary integral cross sections, not
differential-cross-section providers. DCS, screened Coulomb, and ML-prior
angular models are future unsupported work.

Tail metrics are product output fields for high-energy reaction-rate quality.
They report tail probability and reaction-rate tail fractions after any EEDF
postprocess has been applied. They do not prove full distribution convergence:
mean energy can be stable while high-threshold chemistry is still tail
sensitive.

Current limitations:

- no direct PN block operator or DCS-based PN solve
- no arbitrary crossed E-B dynamics
- no RF/time-dependent solver
- no finite-k hydrodynamic transport
- no full electron-electron Fokker-Planck operator
- no particle-particle Coulomb Monte Carlo
- no YAML-side state-resolved or generated-superelastic chemistry

The only supported e-e model is `relaxation_postprocess`. It updates EEDF and
rates for Boltzmann-family solvers, but transport is not recomputed and is
marked stale. Magnetic-field dynamics are unsupported for all canonical solvers
in this phase unless the user explicitly requests an ignored fallback through
`feature_policy.allow_unsupported_fallback`.

Cross-section files may contain explicit `superelastic` processes. The removed
schema fields are only the old generation helpers, which were not a complete
product feature.

Requested physics is checked by the solver plan. Unsupported features fail or
skip according to `feature_policy`; they are not silently ignored.
