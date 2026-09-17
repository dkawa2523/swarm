# Direct PN Operator

`multi_term.method: pn_closure_direct` solves the axisymmetric homogeneous-DC
PN equations directly. The implementation is divided by responsibility:

- `case.py`: product-to-numerical case construction
- `grid.py`: multi-term energy grid
- `operator.py`: collision damping and bidirectional adjacent-moment field
  blocks
- `steady.py`: normalized stationary eigenmode and full-system convergence
- `observables.py`: rates, F1 drift, and F0 diffusion approximation
- `result.py`: canonical product result
- `solver.py`: thin composition entry point

## Unknown layout

Even `F_l` values use the `n` energy-cell centers. Odd values use the `n-1`
internal faces; their two boundary values are zero. Blocks are concatenated in
increasing `l`, and `F0` obeys

```text
sum_i F0_i * DeltaE_i = 1.
```

For an energy-density coefficient, electric acceleration contains the radial
operators

```text
-E c eps^((l+1)/2) d/deps [eps^(-l/2) F_(l-1)]
-E c eps^(-l/2) d/deps [eps^((l+1)/2) F_(l+1)],
c = sqrt(2 e / m_e),
```

multiplied by the Legendre recurrence factors `l/(2l-1)` and
`(l+1)/(2l+3)`. Every available neighbor is assembled; in particular, higher
moments feed back through the chain to `F1` and `F0`.

`F0` contains the conservative inelastic/nonconservative energy redistribution
and finite-temperature elastic energy relaxation. `F_l`, `l>0`, uses angular
damping

```text
nu_l = N v sigma_total (1 - m_l) + nu_inelastic_loss.
```

An effective-momentum-only input is usable only with the explicitly isotropic
ordinary-XS closure, where its effective damping is applied consistently. It
cannot normalize DCS moments.

## Stationary solve and evidence

The solver finds `A F = gamma F`, replaces one `F0` equation by normalization,
and updates the temporal growth eigenvalue from the conservative `F0` balance.
The accepted residual is evaluated afterward on every original PN row, not on
the normalization-replaced matrix. Iteration-limit, non-finite, singular,
excess-negative-mass, or residual failures raise errors; coefficients are not
clipped or repaired.

The result reports:

- the actual multi-term grid and cell widths;
- `pn_full_relative_residual` and iteration/eigenvalue changes;
- `field_coupling=bidirectional_adjacent_legendre_moments`;
- `drift_observable=F1_velocity_moment`;
- `diffusion_observable=F0_gradient_reconstruction`.
