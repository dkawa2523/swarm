# Direct PN Scope and Numerical Gates

The supported scope is a spatially homogeneous, axisymmetric `m=0`, DC-field
PN system without magnetic field. `lmax >= 1` uses the same operator and
stationary solution path; `lmax: 1` is not a special two-term reduction.

Every accepted run must satisfy:

- finite normalized `F0` on the configured multi-term grid;
- negative `F0` mass below `1e-8`, without clipping;
- converged full-vector shape and temporal eigenvalue;
- complete PN residual below the configured residual tolerance;
- positive finite F1 flux drift and F0-derived diffusion;
- angular provenance consistent with `pn_closure_direct` or `pn_dcs`.

Tests exercise `lmax` sensitivity, every residual block, the staggered grid,
ordinary-XS closure metadata, DCS-moment metadata, and fail-closed behavior at
the iteration bound. Agreement with `two_term` remains a useful control-case
comparison, but it is not used to redefine or repair the multi-term solution.
