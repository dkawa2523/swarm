# Direct PN Closure Operator Roadmap

This roadmap is the implementation gate for `multi_term.method:
pn_closure_direct`. The current product runtime must keep that method
`NotImplemented` until the equations, discretization, and lmax=1 regression
below are fixed and tested. The runnable ordinary-cross-section path remains
`pn_closure_surrogate`.

## 1. Candidate Operator Equation

The future direct PN solve should assemble one stationary homogeneous block
system for Legendre coefficients of the electron energy distribution:

```text
A(E, collisions, field) u = 0
normalization(f0) = 1
```

The field-coupling, collision, and inelastic source blocks must be derived from
the PN-expanded Boltzmann equation. They must not be inferred from the current
surrogate closure coefficients. The lmax=1 reduction must reproduce the
existing two-term native Scharfetter-Gummel energy-space result before the
operator is exposed as runnable product code.

## 2. Unknown Vector Layout

Use ell-major, energy-minor ordering:

```text
u = [f_0(E_0..E_N-1), f_1(E_0..E_N-1), ..., f_lmax(E_0..E_N-1)]
```

The `f0` block is the scalar energy distribution that is normalized with the
same cell widths used by the two-term solver. Higher `f_l` blocks are angular
Legendre coefficients on the same energy grid.

## 3. Field Coupling Block

The electric-field operator couples neighboring angular orders
`ell <-> ell-1` and `ell <-> ell+1`. Its energy-space discretization must be
fixed from the PN equation and validated against the existing two-term native
Scharfetter-Gummel operator at `lmax=1`.

Do not substitute the current `field_coupling_scale` or surrogate anisotropy
heuristics for this block.

## 4. Collision Damping Block

Elastic angular damping requires collision frequencies for each Legendre order.
With ordinary integral cross sections, only total and momentum-transfer
information may be available. That is not enough to claim exact DCS-based PN
damping for all `ell`.

The first implementation must either consume validated angular moments/DCS data
or explicitly document and test any closure used for higher-order damping.

## 5. Source And Sink Treatment

Inelastic, ionization, attachment, and superelastic processes need conservative
energy-loss source and sink terms. The existing two-term scalar operator uses
energy-shift deposits for `f0`.

The `l>0` source treatment is not fixed. Isotropic re-emission, momentum-loss
damping, secondary-electron sharing, and growth/source terms must be derived and
regression-tested before direct PN execution is enabled.

## 6. Boundary Conditions

For `f0`, the lmax=1 regression must preserve the two-term zero-flux
Scharfetter-Gummel boundary behavior at the lower and upper energy boundaries.

Boundary conditions for `l>0` are still uncertain. They must be chosen from the
PN flux form, not by copying the scalar `f0` boundary row without validation.
The upper boundary must also satisfy tail-quality checks from product tail
metrics.

## 7. Normalization Constraint

Replace one `f0` equation row with:

```text
sum_i f0_i * delta_E_i = 1
```

Do not normalize higher angular blocks directly. The solve must report the
normalization integral, linear residual, and bounded negative mass diagnostics.

## 8. lmax=1 Regression Condition

Before `pn_closure_direct` can run in product mode, `lmax=1` must match the
existing two-term native Scharfetter-Gummel path on the same grid, cross
sections, gas conditions, and E/N sweep.

Required tolerances:

- mean energy within 1-3%
- drift velocity and mobility within 1-3%
- dominant reaction rates within 1-3% when the same rate integrals are used
- EEDF normalization within numerical tolerance
- nonnegative `f0`, or bounded and reported negative mass
- finite residual metadata

No bulk/source-gradient output may be emitted until separately validated.

## 9. Required Tests

- schema v2 accepts `pn_closure_direct`, and execution fails fast while the
  operator is not implemented
- `pn_dcs` fails fast without a DCS angular-moment provider
- `pn_closure_surrogate` continues to run and records
  `direct_pn_operator=false`
- future direct PN matrix shape matches `(lmax + 1) * n_energy`
- lmax=1 direct PN regression against the two-term native SG solver
- normalization, residual, finite values, and negative-mass diagnostics
- boundary-condition regression on small and refined energy grids
- no direct PN transport or bulk output before validation

## 10. Physically Uncertain Parts

- exact PN electric-field coupling in energy space
- angular collision damping beyond momentum transfer when only ordinary
  integral cross sections are available
- `l>0` source/sink treatment for inelastic and ionizing collisions
- secondary-electron energy sharing in a block PN operator
- nonconservative growth terms in the direct block system
- boundary conditions for higher angular orders
- validated bulk/source-gradient transport from direct PN coefficients
