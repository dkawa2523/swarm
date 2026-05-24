# Direct PN Closure Operator Gate

`multi_term.method: pn_closure_direct` is a planned independent coupled block PN solver.
It is accepted by schema v2, but product execution fails fast until the
operator is implemented from source-of-truth PN equations.

The previous shared-SG reduction was removed from the product path because the
current shared-SG reduction is not an independent PN block solve. It reused the
two-term scalar Scharfetter-Gummel operator and constructed `f1` after the
solve. That can be useful as an internal reduction check, but it is not an
independent multi-term solver validation.

## Candidate Operator Equation

The future direct PN solve must assemble one sparse stationary block system for
Legendre coefficients:

```text
A(E, collisions, field) u = 0
normalization(f0) = 1
u = [f0(E0..EN-1), f1(E0..EN-1), ..., f_lmax(E0..EN-1)]
```

## Unknown Vector Layout

Use ell-major, energy-minor ordering:

```text
index(ell, i) = ell * n_energy + i
```

## Field Coupling Block

The electric-field operator must couple `ell <-> ell-1` and `ell <-> ell+1`
from the PN-expanded Boltzmann equation. It must not wrap the two-term scalar
operator as the direct PN implementation.

## Collision Damping Block

`ell = 1` damping must reduce to the same momentum relaxation frequency used by
the two-term solver. Higher `ell` damping must come from validated angular
moments or a clearly tested ordinary-XS closure.

## Source And Sink Treatment

The future implementation must define inelastic, ionization, attachment, and
superelastic `l>0` source/sink treatment. Fake anisotropic inelastic sources are
not allowed.

## Boundary Conditions

`f0` must satisfy the two-term zero-flux energy boundary in the lmax=1
regression. `l>0` boundary conditions must be derived from the PN flux form.

## Normalization Constraint

Replace one `f0` equation row with:

```text
sum_i f0_i * delta_E_i = 1
```

Do not normalize higher angular blocks directly.

## lmax=1 Regression Condition

A runnable `pn_closure_direct` implementation must pass lmax=1 regression
against the two-term native Scharfetter-Gummel solver without copying the
two-term EEDF:

- coupled f0/f1 sparse block solve
- no two-term scalar operator reuse as the solver
- no post-hoc f1 construction to match drift
- rates recomputed from solved `f0`
- EEDF relative L1 below 0.01
- mean energy within 0.5%
- drift velocity within 1%
- major rates within 2%
- normalization error below `1e-8`
- bounded negative mass diagnostics

Until those conditions are met, `direct_pn_operator=true` is not emitted by
product results.
In plain terms: direct_pn_operator=true is not emitted.

## Required Tests

- schema v2 accepts `pn_closure_direct`
- execution raises `NotImplementedError` while the coupled block operator is missing
- no product result reports `solver_method=pn_closure_direct`
- no product result reports `direct_pn_operator=true`
- `pn_closure_surrogate` remains runnable and records `direct_pn_operator=false`
- future lmax=1 direct PN regression against the two-term native SG solver
- lmax=1 two-term Scharfetter-Gummel regression harness

## Physically Uncertain Parts

- PN-equation-derived energy-space field coupling
- `ell >= 1` collision damping for ordinary integral XS closure
- `l>0` source/sink treatment
- `l>0` boundary conditions
- secondary-electron energy sharing in a block PN operator
- nonconservative growth terms in the direct block system
- validated bulk/source-gradient transport from PN coefficients
