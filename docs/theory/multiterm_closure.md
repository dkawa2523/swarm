# Multi-Term PN Model

The `multi_term` solver expands the homogeneous, axisymmetric electron
distribution in Legendre coefficients `F_l(epsilon)`, `l=0..lmax`, and solves
all coefficients in one stationary sparse eigenproblem. Electric acceleration
couples every neighboring pair in both directions (`l-1 -> l` and
`l+1 -> l`). It is not a two-term energy-diffusion solve with higher moments
added afterward.

Even moments are stored at energy-cell centers and odd moments on internal
cell faces. This staggered finite-volume layout gives a conservative `l=0`
field divergence and avoids a collocated odd-even derivative mode. Odd moments
vanish at the lower and upper energy boundaries. The temporal growth
eigenvalue acts on every PN equation. A run is rejected unless shape,
eigenvalue, and the residual of the complete PN system converge.

The reported drift is the velocity moment of the solved `F1`. Particle
diffusion remains an `F0` gradient reconstruction, so the transport definition
is `pn_f1_flux_drift_f0_gradient_diffusion`.

## Angular input

`pn_closure_direct` uses ordinary integral cross sections and therefore needs
an explicit angular closure:

- `isotropic` sets higher elastic Legendre moments to zero.
- `momentum_power` derives `m1` from total and momentum-transfer integrals and
  closes higher moments with powers of `m1`.
- `maxent_p1` derives the same `m1` and obtains higher moments from its
  maximum-entropy P1 distribution.

These are angular-closure PN calculations, not exact differential-cross-
section calculations. Their results record
`ordinary_integral_xs_closure=true` and `exact_dcs_based=false`.

`pn_dcs` accepts only `physics.angular_scattering.model: moment_table` and uses
the supplied normalized Legendre moments in every angular damping block. It is
marked `exact_dcs_based=true` only when the table provenance is
`dcs_derived`; a model-derived table retains `exact_dcs_based=false`. Raw DCS
parsing is not implemented.
