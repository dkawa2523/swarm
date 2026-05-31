# Direct PN Closure Operator

`multi_term.method: pn_closure_direct` is the product path for the direct PN
closure. The protected `lmax: 1` path is the SG-reduction gate. Higher `lmax`
values run in a limited ordinary-XS angular-closure scope documented below.

## Unknown Vector Layout

The product layout is:

```text
index(ell, i) = ell * n_energy + i
sum_i F0_i * delta_E_i = 1
```

The `lmax: 1` implementation solves a sparse block system for `[F0, G1]`, where
`G1` is the validated SG energy-flux auxiliary. `G1` is not treated as a
physical Legendre coefficient. The direct path uses shared kinetic projection
data from `electron_swarm.solvers.kinetic` and does not instantiate
`TwoTermSolver`, copy the two-term EEDF, or reuse two-term rates.

## Current lmax=1 Equation

For `lmax: 1`, the block is the conservative reduction of the native
Scharfetter-Gummel balance:

```text
collision(F0) + G1 = gamma F0
G1 = energy_flux(F0)
normalization(F0) = 1
```

Eliminating `G1` gives the native two-term SG energy-space balance. Keeping
`G1` inside the sparse solve fixes the direct reduction gate without post-hoc
construction.

## lmax>1 Coefficient-Space Scope

Higher `lmax` solves coefficient blocks `F_l(E_i), l=0..L` in one sparse
system. The accepted coefficient-space equation form is:

```text
0 = C_l[F_l] + S_l
    + E_{l,l-1}[F_{l-1}] + E_{l,l+1}[F_{l+1}]
    - gamma F_l
```

where the `E` blocks use documented Legendre recurrence coefficients. Fitted
field-coupling scale factors are not allowed. A runnable `lmax > 1` product
path uses:

- `ell=1` damping exactly consistent with the shared momentum frequency.
- `ell>=2` damping from angular-closure moments and shared `sigma_total_like`,
  `nu_l = N * v * sigma_total_like * (1 - m_l) + nu_inelastic_loss`, with no
  silent fallback when `sigma_total_like` is unavailable.
- `l=0` source/sink from shared kinetic projection.
- `l>0` inelastic source/sink policy fixed to sink-only for the first accepted
  ordinary-XS scope.
- low/high energy boundary conditions for every `l>0` block.
- residuals computed separately from the normalization row.

Unsupported conditions still fail fast: magnetic PN, non-DC fields,
superelastic lmax>1 treatment, invalid moments/damping, non-finite solves, or
excessive negative `F0` mass. `pn_dcs` is the moment-table entry point for the
same block; `pn_closure_direct` remains the ordinary-XS angular-closure entry
point.

## lmax=1 Regression Condition

The direct implementation is accepted only if it passes the Ar/BOLSIG
regression against the two-term native SG solver:

- coupled sparse block solve
- no two-term EEDF copy
- no post-hoc `G1` construction to match drift
- rates recomputed from solved `F0`
- EEDF relative L1 below 0.01
- mean energy within 0.5%
- drift velocity within 1%
- major rates within 2%
- normalization error below `1e-8`
- negative mass fraction below `1e-8`

## Required Tests

- schema v2 accepts `pn_closure_direct`
- `lmax: 1` reports `solver_method=pn_closure_direct`
- `lmax: 1` reports `direct_pn_operator=true`
- `lmax > 1` smoke runs for the supported ordinary-XS scope
- `pn_dcs` runs only with validated moment-table input
- Ar/BOLSIG lmax=1 two-term SG regression harness
