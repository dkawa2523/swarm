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

For the energy-density convention, electric acceleration contains
`-E c epsilon^((l+1)/2) d/d epsilon [epsilon^(-l/2) F_(l-1)]`
and `-E c epsilon^(-l/2) d/d epsilon [epsilon^((l+1)/2) F_(l+1)]`,
where `c=sqrt(2 e/m_e)`. Their recurrence factors are `l/(2l-1)` and
`(l+1)/(2l+3)`. The same bidirectional operator is used for every `lmax >= 1`.

## Angular input

`pn_closure_direct` uses ordinary integral cross sections and therefore needs
an explicit angular closure:

- `isotropic` sets higher elastic Legendre moments to zero.
- `momentum_power` derives `m1` from total and momentum-transfer integrals and
  closes higher moments with powers of `m1`.
- `maxent_p1` derives the same `m1` and obtains higher moments from its
  maximum-entropy P1 distribution.

Each species is closed before mixing. Angular damping is the sum of
`N_s v sigma_0,s (1-m_l,s)` plus inelastic loss; applying a nonlinear closure
to mixture-averaged integrals would describe a different angular model.
The maximum-entropy moments use scaled modified Bessel functions with a
forward/backward limit, avoiding loss of angular mass in fixed quadrature.
The configured momentum floor regularizes transport and angular damping
consistently; elastic recoil uses the actual physical cross section.

Momentum-only isotropic data use `sigma_0=sigma_m`. Effective-only data cannot
identify elastic recoil and are rejected through the feature policy.

These are angular-closure PN calculations, not exact differential-cross-
section calculations. Their results record
`ordinary_integral_xs_closure=true` and `exact_dcs_based=false`.

`pn_dcs` accepts only `physics.angular_scattering.model: moment_table` and uses
the supplied normalized Legendre moments in every angular damping block. It is
marked `exact_dcs_based=true` only when the table provenance is
`dcs_derived`; a model-derived table retains `exact_dcs_based=false`. Raw DCS
parsing is not implemented.
Moment tables undergo a truncated Hausdorff moment test on `[-1,1]`;
positive semidefinite moment/localizing matrices ensure that an angular
measure exists. A table satisfying only individual `[-1,1]` bounds may fail.

Accepted solutions have normalized finite `F0`, negative mass below `1e-8`,
converged full-vector shape and growth eigenvalue, and a converged residual
on every original PN equation. The normalization equation is excluded only
while solving, not while evaluating that residual. Failed or singular solves
raise errors without coefficient clipping or substitution by another solver.
