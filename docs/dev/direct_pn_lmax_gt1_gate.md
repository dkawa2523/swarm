# Direct PN lmax>1 Scope

`multi_term.method: pn_closure_direct` is runnable for the limited
coefficient-space scope `B=0`, DC field, axisymmetric `m=0`, ordinary integral
cross sections, and angular-closure moments.

## Implemented Model

- Layout: `index(ell, i) = ell * n_energy + i`.
- `lmax: 1` keeps the protected SG-reduction gate `[F0, G1]`; `G1` is the
  energy-flux auxiliary, not a physical Legendre coefficient.
- `lmax > 1` solves a sparse block system for `F0..FL` using Legendre
  recurrence coefficients between neighboring blocks.
- `ell=1` damping is the shared momentum frequency.
- `ell>=2` damping is
  `nu_l = N * v * sigma_total_like * (1 - m_l) + nu_inelastic_loss`.
- `l=0` source/sink uses the shared kinetic collision projection.
- `l>0` inelastic treatment is sink-only; anisotropic inelastic and
  ionization-secondary sources are not generated.
- `F0` is normalized with `sum(F0_i * delta_E_i) = 1`.

## Fail-Fast Conditions

- non-DC field or magnetic PN request
- `moment_table` is routed through `pn_dcs`, not ordinary-XS `pn_closure_direct`
- superelastic lmax>1 treatment
- missing or nonpositive `sigma_total_like`
- invalid angular moments or damping
- non-finite solve or excessive negative `F0` mass

The Ar benchmark still treats `lmax: 1` as the strict regression gate. Higher
orders must remain finite, normalized, and close to the gate in the simple Ar
control case; deviations are benchmark diagnostics, not relaxed tolerances.
